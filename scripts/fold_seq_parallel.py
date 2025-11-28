import os
import tempfile
import subprocess
import shutil
import re
import asyncio
from concurrent.futures import ThreadPoolExecutor
import argparse
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(description='蛋白质序列折叠能级计算')
    parser.add_argument('input_file', help='输入序列文件路径')
    parser.add_argument('output_file', help='输出结果文件路径')
    parser.add_argument('--workers', type=int, default=4, help='并行工作线程数')
    return parser.parse_args()

def fold_seq(seq):
    """
    计算给定序列的最低能级(E)和简并度(N)
    
    参数:
        seq (str): 蛋白质序列字符串
    
    返回:
        tuple: (seq, E, N) 序列、最低能级和简并度
    """
    # 创建临时文件夹
    temp_dir = tempfile.mkdtemp(dir='/tmp')
    try:
        # 创建输入文件
        input_file = os.path.join(temp_dir, 'seq.inp')
        with open(input_file, 'w') as f:
            f.write(f"""RngType  2
RngSeed  5

HPModelDimension  2
OccupancyFieldType  0
MoveFractions  0.75  0.98

HistogramDims  1  1
HistogramSpecs  0  80  1

MCAlgorithm  0
PhysicalModel 0
ModFactorInit  1.0
ModFactorDivider  2.0
ModFactorThreshold  1.0e-8
HistogramCheckInterval  1000000
FlatnessMeasure  0.80
HPSequence {seq}
""")
        
        # 调用外部程序
        process = subprocess.Popen(['/home/wendao/work/Lattice/pepcraft/Sim/MCSim', 'seq.inp'],
                                  cwd=temp_dir,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        # 解析stderr中的最终结果
        stderr_output = stderr.decode('utf-8')
        final_match = re.search(r'FINAL:\s*(-?\d+)\s+(\d+)', stderr_output)
        
        if final_match:
            E = int(final_match.group(1))
            N = int(final_match.group(2))
            return (seq, E, N)
        else:
            raise ValueError(f"无法从stderr中解析FINAL结果，序列: {seq}\nstderr输出:\n{stderr_output}")
            
    finally:
        # 清理临时文件夹
        shutil.rmtree(temp_dir)

async def process_sequences(input_file, output_file, workers=4):
    """
    处理序列文件的主函数
    
    参数:
        input_file: 输入序列文件路径
        output_file: 输出结果文件路径
        workers: 并行工作线程数
    """
    # 读取所有序列
    with open(input_file, 'r') as f:
        all_seqs = [line.strip() for line in f if line.strip()]
    
    # 读取已处理的结果
    processed_seqs = set()
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.split()
                    if len(parts) >= 1:
                        processed_seqs.add(parts[0])
    
    # 过滤掉已处理的序列
    seqs_to_process = [seq for seq in all_seqs if seq not in processed_seqs]
    
    print(f"总序列数: {len(all_seqs)}, 已处理: {len(processed_seqs)}, 待处理: {len(seqs_to_process)}")
    
    # 准备输出文件
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 使用线程池并行处理
    with ThreadPoolExecutor(max_workers=workers) as executor:
        #loop = asyncio.get_event_loop()
        loop = asyncio.get_running_loop()
        futures = []
        
        # 创建进度条
        pbar = tqdm(total=len(seqs_to_process), desc="处理进度")
        
        # 提交任务
        for seq in seqs_to_process:
            future = loop.run_in_executor(executor, fold_seq, seq)
            futures.append(future)
        
        # 处理结果
        with open(output_file, 'a') as out_f:
            for future in asyncio.as_completed(futures):
                try:
                    seq, E, N = await future
                    out_f.write(f"{seq} {E} {N}\n")
                    out_f.flush()  # 确保及时写入文件
                    pbar.update(1)
                except Exception as e:
                    print(f"\n处理序列时出错: {str(e)}")
                    continue
        
        pbar.close()

def main():
    args = parse_args()
    
    asyncio.run(process_sequences(
        args.input_file,
        args.output_file,
        workers=args.workers
    ))
    

if __name__ == "__main__":
    main()
