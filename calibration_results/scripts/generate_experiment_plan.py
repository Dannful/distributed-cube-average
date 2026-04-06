#!/usr/bin/env python3
"""Generate experiment plan for SMPI calibration."""

import itertools
import random
import sys

def generate_plan(output_file, seed=42, num_sizes=100, repetitions=20):
    """
    Generate experiment plan with random log-uniform message sizes.

    Args:
        output_file: Path to output CSV file
        seed: Random seed for reproducibility
        num_sizes: Number of unique message sizes to test
        repetitions: Number of times to repeat each (op, size) pair
    """
    random.seed(seed)

    ops = ['Recv', 'Isend', 'PingPong', 'Wtime', 'Iprobe', 'Test']
    sizes = [int(10**random.uniform(0, 9)) for _ in range(num_sizes)]

    experiments = list(itertools.product(ops, sizes))
    experiments *= repetitions
    random.shuffle(experiments)

    with open(output_file, 'w') as f:
        for op, size in experiments:
            f.write(f'{op},{size}\n')

    print(f'Generated {len(experiments)} experiments:')
    print(f'  Operations: {len(ops)}')
    print(f'  Unique sizes: {len(sizes)}')
    print(f'  Repetitions: {repetitions}')
    print(f'  Output: {output_file}')

if __name__ == '__main__':
    output = sys.argv[1] if len(sys.argv) > 1 else 'exp_plan.csv'
    generate_plan(output)
