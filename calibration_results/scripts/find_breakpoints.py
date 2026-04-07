#!/usr/bin/env python3
"""Binary search for SMPI semantic breakpoints."""

import subprocess
import sys
import shlex


def deadlocked(bp_search1_cmd, mpi_opt, size):
    args = ['mpirun', *mpi_opt, *shlex.split(bp_search1_cmd), str(size)]
    print(' '.join(args))
    try:
        result = subprocess.run(args, timeout=3, check=True,
                              stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        return False
    except subprocess.TimeoutExpired:
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: bp_search1 failed with exit code {e.returncode}", file=sys.stderr)
        if e.stderr:
            print(f"STDERR: {e.stderr.decode('ascii', errors='replace')}", file=sys.stderr)
        sys.exit(1)


def corrupted(bp_search2_cmd, mpi_opt, size):
    args = ['mpirun', *mpi_opt, *shlex.split(bp_search2_cmd), str(size)]
    print(' '.join(args))
    try:
        result = subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = result.stdout.decode('ascii').strip()
        return 'corrupted' in output
    except subprocess.CalledProcessError as e:
        print(f"ERROR: bp_search2 failed with exit code {e.returncode}", file=sys.stderr)
        if e.stderr:
            print(f"STDERR: {e.stderr.decode('ascii', errors='replace')}", file=sys.stderr)
        sys.exit(1)


class SearchError(Exception):
    pass


def run_search(bp_cmd, mpi_opt, test_func, min_size=1, max_size=1000000000):
    if test_func(bp_cmd, mpi_opt, min_size) or not test_func(bp_cmd, mpi_opt, max_size):
        raise SearchError('no breakpoint found in the interval [%d, %d]' % (min_size, max_size))

    def find_max():  # used to speed-up the search
        size = min_size*2
        while size < max_size and not test_func(bp_cmd, mpi_opt, size):
            size *= 2
        return size

    max_size = min(find_max(), max_size)
    min_size = max_size // 2
    # we search bp such that test_func returns False for bp and True for bp+1
    # invariant: bp is in [min_size, max_size]
    while min_size != max_size:
        size = (min_size + max_size)//2
        if test_func(bp_cmd, mpi_opt, size):
            max_size = size - 1
        else:
            if min_size == size:
                assert max_size == min_size + 1
                if test_func(bp_cmd, mpi_opt, max_size):
                    return min_size
                else:
                    return max_size
            else:
                min_size = size
    return size


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit('Syntax: %s <bp_search1_cmd> <bp_search2_cmd> <mpi_options>' % sys.argv[0])

    bp_search1_cmd = sys.argv[1]
    bp_search2_cmd = sys.argv[2]
    mpi_opts = sys.argv[3:]

    try:
        bp1 = run_search(bp_search2_cmd, mpi_opts, corrupted)
    except SearchError as e:
        bp1 = str(e)
    try:
        bp2 = run_search(bp_search1_cmd, mpi_opts, deadlocked)
    except SearchError as e:
        bp2 = str(e)
    print()
    print('First breakpoint:  %s' % bp1)
    print('Second breakpoint: %s' % bp2)
