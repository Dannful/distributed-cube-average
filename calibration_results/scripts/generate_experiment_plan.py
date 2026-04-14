#!/usr/bin/env python3
"""Generate experiment plan for SMPI calibration."""

import argparse
import itertools
import math
import random


def generate_plan(output_file, seed=42, num_sizes=150, repetitions=15, max_size=134217728):
    """
    Generate experiment plan with random log-uniform message sizes.

    Args:
        output_file: Path to output CSV file
        seed: Random seed for reproducibility
        num_sizes: Number of unique message sizes to test
        repetitions: Number of times to repeat each (op, size) pair
        max_size: Maximum message size in bytes
    """
    random.seed(seed)

    ops = ['Recv', 'Isend', 'PingPong', 'Wtime', 'Iprobe', 'Test']
    log_max = math.log10(max_size)
    sizes = sorted(set(int(10**random.uniform(0, log_max)) for _ in range(num_sizes * 2)))[:num_sizes]

    experiments = list(itertools.product(ops, sizes))
    experiments *= repetitions
    random.shuffle(experiments)

    with open(output_file, 'w') as f:
        for op, size in experiments:
            f.write(f'{op},{size}\n')

    print(f'Generated {len(experiments)} experiments:')
    print(f'  Operations: {len(ops)}')
    print(f'  Unique sizes: {len(sizes)} (range: 1–{max_size})')
    print(f'  Repetitions: {repetitions}')
    print(f'  Output: {output_file}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate SMPI calibration experiment plan')
    parser.add_argument('output', nargs='?', default='exp_plan.csv')
    parser.add_argument('--num-sizes', type=int, default=150)
    parser.add_argument('--repetitions', type=int, default=15)
    parser.add_argument('--max-size', type=int, default=134217728,
                        help='Maximum message size in bytes (default: 256^3*8 = 134217728)')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    generate_plan(args.output, seed=args.seed, num_sizes=args.num_sizes,
                  repetitions=args.repetitions, max_size=args.max_size)
