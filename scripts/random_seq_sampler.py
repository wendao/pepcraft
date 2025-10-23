import sys
import random

def generate_sequence(L, alphabet_size):
    if alphabet_size == 2:
        alphabet = ['H', 'P']
    elif alphabet_size == 3:
        alphabet = ['H', 'P', 'I']
    else:
        raise ValueError("字母表长度必须是2或3")
    
    return ''.join(random.choice(alphabet) for _ in range(L))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("用法: python script.py <序列长度L> <字母表长度(2或3)> <条数>")
        sys.exit(1)
    
    try:
        L = int(sys.argv[1])
        alphabet_size = int(sys.argv[2])
        N = int(sys.argv[3])
    except ValueError:
        print("错误: 参数必须是整数")
        sys.exit(1)
    
    if alphabet_size not in [2, 3]:
        print("错误: 字母表长度必须是2或3")
        sys.exit(1)
    
    for i in range(N):
        sequence = generate_sequence(L, alphabet_size)
        print(sequence)
