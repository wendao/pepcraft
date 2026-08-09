#!/usr/bin/env python3
"""Train and evaluate a GRU generator for fixed-length lattice SAWs.

The script is the GRU counterpart of the RBM training workflow.  It uses the
same conformation loader, geometric validity checker, novelty definition, and
output metrics as the rest of this repository.

Example
-------
python3 GRU.py \
    --data /path/to/conf1.txt \
    --length 20 \
    --skip-lines 1 \
    --relative-map LFR \
    --output runs/l20/gru
"""

from __future__ import annotations

import argparse
from pathlib import Path

from sawgen.train import run_training


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Train an autoregressive GRU on relative-turn SAW conformations, "
            "then generate and evaluate new conformations."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    data = parser.add_argument_group("data")
    data.add_argument("--data", required=True, help="Input conformation file")
    data.add_argument("--length", type=int, required=True, help="Number of monomers L")
    data.add_argument(
        "--input-format",
        default="auto",
        choices=[
            "auto",
            "coords",
            "absolute",
            "relative",
            "onehot-absolute",
            "onehot-relative",
        ],
    )
    data.add_argument("--skip-cols", type=int, default=0)
    data.add_argument("--skip-lines", type=int, default=0)
    data.add_argument(
        "--relative-map",
        default="FLR",
        help="Meaning of raw relative-turn ids 0,1,2 (for example LFR)",
    )
    data.add_argument(
        "--allow-invalid",
        action="store_true",
        help="Skip malformed/non-self-avoiding rows instead of failing",
    )
    data.add_argument(
        "--max-samples",
        type=int,
        default=100_000,
        help="Maximum number of input conformations; 0 uses all rows",
    )

    model = parser.add_argument_group("GRU model")
    model.add_argument("--hidden-size", type=int, default=128)
    model.add_argument("--embedding-size", type=int, default=64)
    model.add_argument("--num-layers", type=int, default=2)
    model.add_argument("--dropout", type=float, default=0.1)

    training = parser.add_argument_group("training")
    training.add_argument("--epochs", type=int, default=15)
    training.add_argument("--batch-size", type=int, default=256)
    training.add_argument("--lr", type=float, default=3e-4)
    training.add_argument("--weight-decay", type=float, default=1e-4)
    training.add_argument("--grad-clip", type=float, default=1.0)
    training.add_argument("--seed", type=int, default=42)
    training.add_argument(
        "--device",
        default="auto",
        help="auto, cpu, mps, cuda, or a concrete device such as cuda:0",
    )
    training.add_argument("--workers", type=int, default=0)

    sampling = parser.add_argument_group("generation and evaluation")
    sampling.add_argument("--temperature", type=float, default=1.0)
    sampling.add_argument(
        "--eval-samples",
        type=int,
        default=100_000,
        help="Number of conformations generated after training",
    )
    sampling.add_argument("--eval-batch-size", type=int, default=1_000)
    sampling.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Directory for checkpoints, history, samples, and metrics",
    )
    return parser


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    positive = {
        "--length": args.length,
        "--hidden-size": args.hidden_size,
        "--embedding-size": args.embedding_size,
        "--num-layers": args.num_layers,
        "--epochs": args.epochs,
        "--batch-size": args.batch_size,
        "--eval-samples": args.eval_samples,
        "--eval-batch-size": args.eval_batch_size,
    }
    for option, value in positive.items():
        if value <= 0:
            parser.error(f"{option} must be positive")
    if args.length < 3:
        parser.error("--length must be at least 3")
    if args.max_samples < 0:
        parser.error("--max-samples must be non-negative")
    if args.dropout < 0 or args.dropout >= 1:
        parser.error("--dropout must satisfy 0 <= dropout < 1")
    if args.lr <= 0:
        parser.error("--lr must be positive")
    if args.weight_decay < 0:
        parser.error("--weight-decay must be non-negative")
    if args.grad_clip <= 0:
        parser.error("--grad-clip must be positive")
    if args.temperature < 0:
        parser.error("--temperature must be non-negative")
    normalized_map = args.relative_map.upper()
    if len(normalized_map) != 3 or set(normalized_map) != {"F", "L", "R"}:
        parser.error("--relative-map must be a permutation of FLR")
    args.relative_map = normalized_map


def add_shared_training_defaults(args: argparse.Namespace) -> argparse.Namespace:
    """Add fields used by the shared trainer but irrelevant to a GRU."""
    args.model = "gru"
    args.output = str(args.output)

    # The shared trainer constructs one ModelConfig for every model family.
    # These values are serialized for compatibility but are unused by a GRU.
    args.latent_size = 32
    args.nhead = 4
    args.dim_feedforward = 256
    args.noise_size = 64
    args.critic_steps = 5
    args.gp_weight = 10.0
    args.gan_softmax_temperature = 1.0
    return args


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    result = run_training(add_shared_training_defaults(args))
    print(f"GRU run complete: {result['output_dir']}")


if __name__ == "__main__":
    main()
