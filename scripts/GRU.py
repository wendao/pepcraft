#!/usr/bin/env python
"""Standalone GRU generator for fixed-length 2D lattice SAWs.

Only NumPy and PyTorch are required.  The loader accepts the compact format
used by RBM-tor.py (a metadata header followed by digit strings) as well as
FLR strings and whitespace-separated category ids.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F
from torch.utils.data import DataLoader, TensorDataset


@dataclass(frozen=True)
class GRUConfig:
    seq_len: int
    n_categories: int = 3
    hidden_size: int = 128
    embedding_size: int = 64
    num_layers: int = 2
    dropout: float = 0.1


class AutoregressiveGRU(nn.Module):
    def __init__(self, config: GRUConfig) -> None:
        super().__init__()
        self.config = config
        self.bos = config.n_categories
        self.embedding = nn.Embedding(config.n_categories + 1, config.embedding_size)
        self.gru = nn.GRU(
            config.embedding_size,
            config.hidden_size,
            num_layers=config.num_layers,
            dropout=config.dropout if config.num_layers > 1 else 0.0,
            batch_first=True,
        )
        self.output = nn.Linear(config.hidden_size, config.n_categories)

    def loss(self, targets: torch.Tensor) -> torch.Tensor:
        bos = torch.full(
            (targets.shape[0], 1), self.bos, dtype=torch.long, device=targets.device
        )
        inputs = torch.cat((bos, targets[:, :-1]), dim=1)
        features, _ = self.gru(self.embedding(inputs))
        logits = self.output(features)
        return F.cross_entropy(logits.flatten(0, 1), targets.flatten())

    @torch.inference_mode()
    def sample(self, count: int, temperature: float, device: torch.device) -> torch.Tensor:
        current = torch.full((count,), self.bos, dtype=torch.long, device=device)
        hidden = None
        generated: list[torch.Tensor] = []
        for _ in range(self.config.seq_len):
            features, hidden = self.gru(self.embedding(current[:, None]), hidden)
            logits = self.output(features[:, -1])
            if temperature == 0:
                current = logits.argmax(dim=-1)
            else:
                probabilities = torch.softmax(logits / temperature, dim=-1)
                current = torch.multinomial(probabilities, 1).squeeze(1)
            generated.append(current)
        return torch.stack(generated, dim=1)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def choose_device(requested: str, gpu: int | None) -> torch.device:
    if requested != "auto":
        device = torch.device(requested)
    elif gpu is not None and gpu < 0:
        device = torch.device("cpu")
    elif torch.cuda.is_available():
        index = 0 if gpu is None else gpu
        device = torch.device(f"cuda:{index}")
    elif torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable")
    if device.type == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("MPS was requested but is unavailable")
    return device


def _metadata_header(line: str) -> tuple[int, int] | None:
    fields = line.split()
    if len(fields) < 2 or not all(field.isdigit() for field in fields):
        return None
    return int(fields[0]), int(fields[1])


def _infer_seq_len(lines: Sequence[str], skip_lines: int) -> int:
    candidates = [line.strip() for line in lines[skip_lines:] if line.strip()]
    if not candidates:
        raise ValueError("Input contains no conformation records")
    header = _metadata_header(candidates[0]) if skip_lines == 0 else None
    if header is not None:
        return header[1]
    compact = "".join(candidates[0].split())
    return len(compact) if compact.isalnum() else len(candidates[0].split())


def _parse_record(
    line: str, seq_len: int, n_categories: int, relative_map: str
) -> tuple[int, ...]:
    compact = "".join(line.split())
    if len(compact) == seq_len and compact.isdigit():
        values = tuple(int(char) for char in compact)
    elif len(compact) == seq_len and set(compact.upper()) <= {"F", "L", "R"}:
        raw_id = {name: index for index, name in enumerate(relative_map)}
        values = tuple(raw_id[char] for char in compact.upper())
    else:
        fields = line.replace(",", " ").split()
        if len(fields) != seq_len or not all(field.isdigit() for field in fields):
            raise ValueError(f"expected {seq_len} category ids")
        values = tuple(int(field) for field in fields)
    if any(value < 0 or value >= n_categories for value in values):
        raise ValueError(f"category ids must be in [0, {n_categories - 1}]")
    return values


def load_data(
    path: Path,
    seq_len: int | None,
    n_categories: int,
    skip_lines: int,
    allow_invalid: bool,
    max_samples: int,
    seed: int,
    relative_map: str,
) -> tuple[list[tuple[int, ...]], dict[str, object]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    inferred = _infer_seq_len(lines, skip_lines)
    if seq_len is None:
        seq_len = inferred
    if inferred != seq_len:
        raise ValueError(f"data records have {inferred} tokens, expected {seq_len}")

    start = skip_lines
    declared_count = None
    if skip_lines == 0:
        nonblank_index = next((i for i, line in enumerate(lines) if line.strip()), None)
        if nonblank_index is not None:
            header = _metadata_header(lines[nonblank_index])
            if header is not None:
                declared_count, declared_width = header
                if declared_width != seq_len:
                    raise ValueError(
                        f"header declares {declared_width} tokens, expected {seq_len}"
                    )
                start = nonblank_index + 1

    records: list[tuple[int, ...]] = []
    invalid: list[str] = []
    for line_number, raw in enumerate(lines[start:], start=start + 1):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        try:
            record = _parse_record(stripped, seq_len, n_categories, relative_map)
            path_coords = turns_to_coords(record, relative_map)
            if len(path_coords) != len(set(path_coords)):
                raise ValueError("conformation is not self-avoiding")
            records.append(record)
        except ValueError as exc:
            message = f"line {line_number}: {exc}; raw={stripped[:80]!r}"
            invalid.append(message)
            if not allow_invalid:
                raise ValueError(message) from exc
    if not records:
        raise ValueError(f"No conformations loaded from {path}")
    if declared_count is not None and declared_count != len(records) and not allow_invalid:
        raise ValueError(
            f"header declares {declared_count} records, but {len(records)} were loaded"
        )
    original_count = len(records)
    if max_samples and len(records) > max_samples:
        rng = np.random.default_rng(seed)
        indices = rng.choice(len(records), max_samples, replace=False)
        records = [records[int(index)] for index in indices]
    report = {
        "path": str(path.resolve()),
        "n_loaded": original_count,
        "n_used": len(records),
        "seq_len": seq_len,
        "n_invalid": len(invalid),
        "invalid_examples": invalid[:10],
        "header_detected": declared_count is not None,
    }
    return records, report


def turns_to_coords(raw: Sequence[int], relative_map: str) -> tuple[tuple[int, int], ...]:
    semantic = {raw_id: name for raw_id, name in enumerate(relative_map)}
    coords: list[tuple[int, int]] = [(0, 0), (1, 0)]
    direction = 0
    vectors = ((1, 0), (0, 1), (-1, 0), (0, -1))
    for value in raw:
        turn = semantic[int(value)]
        if turn == "L":
            direction = (direction + 1) % 4
        elif turn == "R":
            direction = (direction - 1) % 4
        dx, dy = vectors[direction]
        x, y = coords[-1]
        coords.append((x + dx, y + dy))
    return tuple(coords)


def canonical_code(raw: Sequence[int], relative_map: str) -> str:
    coords = turns_to_coords(raw, relative_map)

    def variants(path: Sequence[tuple[int, int]]) -> list[tuple[tuple[int, int], ...]]:
        result = []
        for reflect in (False, True):
            for rotations in range(4):
                transformed = []
                for x0, y0 in path:
                    x, y = (-x0, y0) if reflect else (x0, y0)
                    for _ in range(rotations):
                        x, y = -y, x
                    transformed.append((x, y))
                ox, oy = transformed[0]
                result.append(tuple((x - ox, y - oy) for x, y in transformed))
        return result

    candidates = variants(coords) + variants(tuple(reversed(coords)))
    return ";".join(f"{x},{y}" for x, y in min(candidates))


def evaluate_samples(
    samples: Sequence[Sequence[int]],
    training: Sequence[Sequence[int]],
    relative_map: str,
    seed: int,
) -> dict[str, int | float]:
    valid = []
    paths = []
    for sample in samples:
        path = turns_to_coords(sample, relative_map)
        if len(path) == len(set(path)):
            valid.append(tuple(map(int, sample)))
            paths.append(path)
    valid_codes = {canonical_code(sample, relative_map) for sample in valid}
    train_codes = {canonical_code(sample, relative_map) for sample in training}
    novel_codes = valid_codes - train_codes

    rng = np.random.default_rng(seed)
    distances = []
    max_pairs = len(valid) * (len(valid) - 1) // 2
    for _ in range(min(20_000, max_pairs)):
        i, j = rng.choice(len(valid), 2, replace=False)
        a, b = valid[int(i)], valid[int(j)]
        distances.append(sum(x != y for x, y in zip(a, b)) / len(a))

    rg2, end2, contacts = [], [], []
    for path in paths:
        array = np.asarray(path, dtype=np.float64)
        rg2.append(float(np.mean(np.sum((array - array.mean(axis=0)) ** 2, axis=1))))
        end2.append(float(np.sum((array[-1] - array[0]) ** 2)))
        contact_count = 0
        for i, (x1, y1) in enumerate(path):
            for x2, y2 in path[i + 2 :]:
                contact_count += abs(x1 - x2) + abs(y1 - y2) == 1
        contacts.append(float(contact_count))

    safe = lambda numerator, denominator: float(numerator / denominator) if denominator else 0.0
    return {
        "n_generated": len(samples),
        "n_valid": len(valid),
        "n_unique_valid": len(valid_codes),
        "n_unique_novel": len(novel_codes),
        "validity": safe(len(valid), len(samples)),
        "unique_among_valid": safe(len(valid_codes), len(valid)),
        "novelty_among_unique": safe(len(novel_codes), len(valid_codes)),
        "effective_innovation_rate": safe(len(novel_codes), len(samples)),
        "mean_pairwise_hamming": float(np.mean(distances)) if distances else 0.0,
        "mean_rg2": float(np.mean(rg2)) if rg2 else 0.0,
        "mean_end_to_end2": float(np.mean(end2)) if end2 else 0.0,
        "mean_contacts": float(np.mean(contacts)) if contacts else 0.0,
    }


def split_loaders(
    records: Sequence[Sequence[int]], batch_size: int, validation: float, seed: int
) -> tuple[DataLoader, DataLoader]:
    tensor = torch.tensor(np.asarray(records), dtype=torch.long)
    generator = torch.Generator().manual_seed(seed)
    order = torch.randperm(len(tensor), generator=generator)
    n_validation = max(1, int(round(validation * len(tensor)))) if len(tensor) > 1 else 0
    validation_tensor = tensor[order[:n_validation]] if n_validation else tensor
    training_tensor = tensor[order[n_validation:]] if n_validation else tensor
    if len(training_tensor) == 0:
        raise ValueError("validation split leaves no training records")
    train_loader = DataLoader(
        TensorDataset(training_tensor),
        batch_size=min(batch_size, len(training_tensor)),
        shuffle=True,
        generator=generator,
    )
    validation_loader = DataLoader(
        TensorDataset(validation_tensor),
        batch_size=min(batch_size, len(validation_tensor)),
        shuffle=False,
    )
    return train_loader, validation_loader


def average_loss(model: AutoregressiveGRU, loader: DataLoader, device: torch.device) -> float:
    model.eval()
    total, count = 0.0, 0
    with torch.inference_mode():
        for (batch,) in loader:
            batch = batch.to(device)
            total += float(model.loss(batch)) * len(batch)
            count += len(batch)
    return total / count


def atomic_save(payload: dict[str, object], path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    torch.save(payload, temporary)
    os.replace(temporary, path)


def write_history(rows: Sequence[dict[str, int | float]], path: Path) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-d", "--data", required=True, type=Path)
    parser.add_argument("--length", type=int, help="Number of monomers; inferred when omitted")
    parser.add_argument("-v", "--visible", type=int, help="Number of relative-turn tokens (L-2)")
    parser.add_argument("-l", "--level", "--categories", dest="n_categories", type=int, default=3)
    parser.add_argument("--skip-lines", type=int, default=0)
    parser.add_argument("--relative-map", default="LFR", help="Raw ids 0,1,2 mapped to these turns")
    parser.add_argument("--allow-invalid", action="store_true")
    parser.add_argument("--max-samples", type=int, default=100_000)

    parser.add_argument("-n", "--hidden", "--hidden-size", dest="hidden_size", type=int, default=128)
    parser.add_argument("--embedding-size", type=int, default=64)
    parser.add_argument("--num-layers", type=int, default=2)
    parser.add_argument("--dropout", type=float, default=0.1)

    parser.add_argument("-t", "--train", "--epochs", dest="epochs", type=int, default=15)
    parser.add_argument("-b", "--batch_size", "--batch-size", dest="batch_size", type=int, default=256)
    parser.add_argument("-r", "--learning_rate", "--lr", dest="lr", type=float, default=3e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--grad-clip", type=float, default=1.0)
    parser.add_argument("--validation", type=float, default=0.1)
    parser.add_argument("-i", "--interval", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", default="auto")
    parser.add_argument("-x", "--gpu", type=int)

    parser.add_argument("-g", "--generate", "--eval-samples", dest="eval_samples", type=int, default=100_000)
    parser.add_argument("--eval-batch-size", type=int, default=1_000)
    parser.add_argument("--temperature", type=float, default=1.0)
    parser.add_argument("-w", "--weight", "--resume", dest="checkpoint", type=Path)
    parser.add_argument("--output", type=Path, default=Path("runs/gru"))
    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> int | None:
    if args.length is not None and args.length < 3:
        parser.error("--length must be at least 3")
    seq_len = args.length - 2 if args.length is not None else args.visible
    if args.visible is not None and seq_len is not None and args.visible != seq_len:
        parser.error("--visible must equal --length - 2")
    if args.n_categories != 3:
        parser.error("this SAW implementation requires --level=3")
    if len(args.relative_map) != 3 or set(args.relative_map.upper()) != {"F", "L", "R"}:
        parser.error("--relative-map must be a permutation of FLR")
    args.relative_map = args.relative_map.upper()
    if not 0 < args.validation < 1:
        parser.error("--validation must lie between 0 and 1")
    for name in ("hidden_size", "embedding_size", "num_layers", "batch_size", "interval", "eval_batch_size"):
        if getattr(args, name) <= 0:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if args.epochs < 0 or args.eval_samples < 0 or args.max_samples < 0:
        parser.error("epochs, generated samples, and max samples must be non-negative")
    if not 0 <= args.dropout < 1 or args.lr <= 0 or args.grad_clip <= 0:
        parser.error("require 0 <= dropout < 1, lr > 0, and grad-clip > 0")
    if args.temperature < 0:
        parser.error("--temperature must be non-negative")
    return seq_len


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    requested_seq_len = validate_args(parser, args)
    set_seed(args.seed)
    device = choose_device(args.device, args.gpu)

    records, load_report = load_data(
        args.data,
        requested_seq_len,
        args.n_categories,
        args.skip_lines,
        args.allow_invalid,
        args.max_samples,
        args.seed,
        args.relative_map,
    )
    seq_len = int(load_report["seq_len"])
    if args.length is None:
        args.length = seq_len + 2
    config = GRUConfig(
        seq_len=seq_len,
        n_categories=args.n_categories,
        hidden_size=args.hidden_size,
        embedding_size=args.embedding_size,
        num_layers=args.num_layers,
        dropout=args.dropout,
    )
    model = AutoregressiveGRU(config).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    train_loader, validation_loader = split_loaders(records, args.batch_size, args.validation, args.seed)

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    best_path, latest_path = output / "best.pt", output / "latest.pt"
    start_epoch, best_validation = 0, math.inf
    history: list[dict[str, int | float]] = []
    if args.checkpoint is not None:
        payload = torch.load(args.checkpoint, map_location=device, weights_only=False)
        saved_config = GRUConfig(**payload["config"])
        if saved_config != config:
            raise ValueError(f"checkpoint config {saved_config} does not match requested {config}")
        model.load_state_dict(payload["model_state"])
        if "optimizer_state" in payload:
            optimizer.load_state_dict(payload["optimizer_state"])
        start_epoch = int(payload.get("epoch", 0))
        best_validation = float(payload.get("best_validation", math.inf))
        history = list(payload.get("history", []))

    run_config = {
        "args": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "model_config": asdict(config),
        "load_report": load_report,
        "n_parameters": sum(parameter.numel() for parameter in model.parameters()),
        "device": str(device),
        "torch_version": torch.__version__,
    }
    (output / "run_config.json").write_text(
        json.dumps(run_config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    started = time.perf_counter()
    for epoch in range(start_epoch + 1, args.epochs + 1):
        epoch_started = time.perf_counter()
        model.train()
        total, count = 0.0, 0
        for (batch,) in train_loader:
            batch = batch.to(device)
            optimizer.zero_grad(set_to_none=True)
            loss = model.loss(batch)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), args.grad_clip)
            optimizer.step()
            total += float(loss.detach()) * len(batch)
            count += len(batch)
        train_loss = total / count
        validation_loss = average_loss(model, validation_loader, device)
        row = {
            "epoch": epoch,
            "train_loss": train_loss,
            "validation_loss": validation_loss,
            "seconds": time.perf_counter() - epoch_started,
        }
        history.append(row)
        payload = {
            "config": asdict(config),
            "model_state": model.state_dict(),
            "optimizer_state": optimizer.state_dict(),
            "epoch": epoch,
            "best_validation": min(best_validation, validation_loss),
            "history": history,
            "run_config": run_config,
        }
        atomic_save(payload, latest_path)
        if validation_loss < best_validation:
            best_validation = validation_loss
            payload["best_validation"] = best_validation
            atomic_save(payload, best_path)
        if epoch % args.interval == 0 or epoch == args.epochs:
            print(
                f"[GRU] epoch {epoch}/{args.epochs} train={train_loss:.6f} "
                f"validation={validation_loss:.6f} seconds={row['seconds']:.2f}",
                flush=True,
            )
    write_history(history, output / "history.csv")

    checkpoint_for_sampling = best_path if best_path.exists() else args.checkpoint
    if checkpoint_for_sampling is not None:
        payload = torch.load(checkpoint_for_sampling, map_location=device, weights_only=False)
        model.load_state_dict(payload["model_state"])
        best_validation = float(payload.get("best_validation", best_validation))
    model.eval()

    generated: list[list[int]] = []
    remaining = args.eval_samples
    while remaining:
        count = min(remaining, args.eval_batch_size)
        generated.extend(model.sample(count, args.temperature, device).cpu().tolist())
        remaining -= count
    with (output / "samples.txt").open("w", encoding="utf-8") as handle:
        for row in generated:
            handle.write("".join(map(str, row)) + "\n")

    sample_metrics = evaluate_samples(generated, records, args.relative_map, args.seed)
    result = {
        **run_config,
        "training_seconds": time.perf_counter() - started,
        "best_validation_score": best_validation,
        "sample_metrics": sample_metrics,
        "output_dir": str(output),
    }
    (output / "metrics.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    if args.eval_samples == 100_000:
        (output / "eval_100k.json").write_text(
            json.dumps(sample_metrics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
    print(json.dumps(sample_metrics, ensure_ascii=False, indent=2))
    print(f"GRU run complete: {output}")


if __name__ == "__main__":
    main()
