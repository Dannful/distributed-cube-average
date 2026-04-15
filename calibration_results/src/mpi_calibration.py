#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pycewise
import plotnine
from plotnine import *
from pathlib import Path
import argparse
import yaml

print(f'pycewise version: {pycewise.__version__}')

# ==========================================
# Configuration & Constants
# ==========================================
expected_bandwidth = 937e6 / 8
expected_latency   = 22.7e-6

expected_loopback_bandwidth = 2e11 / 8
expected_loopback_latency   = 1e-7

plotnine.options.figure_size = (12, 8)

# ==========================================
# Data Loading
# ==========================================
def extract_local_dir(dir_name='.'):
    base_path = Path(dir_name)
    result = {}
    for file_path in base_path.rglob('*'):
        if not file_path.is_file():
            continue
        name = file_path.relative_to(base_path).as_posix()
        if name.endswith('.csv'):
            dataframe = pd.read_csv(file_path, names=['op', 'msg_size', 'start', 'duration'])
            dataframe['type'] = name
            dataframe['index'] = range(len(dataframe))
            result[name] = dataframe
        elif name.endswith('.yaml'):
            with open(file_path, 'r') as f:
                result[name] = yaml.load(f, Loader=yaml.BaseLoader)
        else:
            with open(file_path, 'rb') as f:
                result[name] = f.read()
    return result

def aggregate_dataframe(dataframe):
    df = dataframe.groupby('msg_size').median(numeric_only=True).reset_index()
    for col in dataframe.columns:
        if col not in df.columns:
            value = list(dataframe[col].unique())
            if len(value) == 1:
                df[col] = value[0]
    return df

def load_experiment(dir_name='.'):
    result = extract_local_dir(dir_name)

    breakpoints_key = next((k for k in result.keys() if k.endswith('breakpoints')), None)
    if breakpoints_key:
        raw_lines = result[breakpoints_key].decode().strip().split('\n')
        semantic_breaks = [int(line.strip()) for line in raw_lines if line.strip()]
        semantic_breaks.sort()
        if len(semantic_breaks) == 1:
            semantic_breaks = [semantic_breaks[0], semantic_breaks[0]]
        elif len(semantic_breaks) >= 2:
            semantic_breaks = semantic_breaks[-2:]
        else:
            semantic_breaks = [65536, 65536]
    else:
        print("WARNING: 'breakpoints' file missing. Defaulting to [65536, 65536].")
        semantic_breaks = [65536, 65536]

    pingpong_key = next((k for k in result.keys() if k.endswith('exp_PingPong.csv')), None)
    if not pingpong_key:
        raise FileNotFoundError("Could not find '*exp_PingPong.csv'.")

    df_pingpong   = result[pingpong_key]
    pingpong_send = df_pingpong[df_pingpong.op == 'MPI_Send'].reset_index(drop=True)
    pingpong_recv = df_pingpong[df_pingpong.op == 'MPI_Recv'].reset_index(drop=True)
    df_pingpong_combined = pd.DataFrame(dict(
        op       = 'PingPong',
        msg_size = pingpong_send.msg_size,
        start    = pingpong_send.start,
        duration = pingpong_send.duration + pingpong_recv.duration,
    ))

    def get_csv(name):
        return next((result[k] for k in result.keys() if k.endswith(name)), None)

    return {
        'pingpong': df_pingpong_combined,
        'send'    : pingpong_send,
        'isend'   : get_csv('exp_Isend.csv'),
        'recv'    : get_csv('exp_Recv.csv'),
        'wtime'   : get_csv('exp_Wtime.csv'),
        'test'    : get_csv('exp_Test.csv'),
        'iprobe'  : get_csv('exp_Iprobe.csv'),
    }, semantic_breaks

# ==========================================
# Regression & Formatting Functions
# ==========================================
def draw_plot(reg, filename, alpha=1):
    reg.plot_dataset(log=True, alpha=alpha)
    plt.xlabel('Message size (bytes)')
    plt.ylabel('Duration (seconds)')
    plt.savefig(filename)
    plt.clf()

def my_join(*values, inner_sep=':', outer_sep=';'):
    return outer_sep.join([inner_sep.join([str(it) for it in items]) for items in zip(*values)])

def regression_to_str(reg_df, drop_threshold=None, keys=['intercept', 'coefficient']):
    if drop_threshold is not None:
        reg_df.loc[reg_df.min_x >= drop_threshold, 'intercept']   = 0
        reg_df.loc[reg_df.min_x >= drop_threshold, 'coefficient'] = 0
    reg_df.loc[reg_df.min_x < 0, 'min_x'] = 0
    reg_df.min_x = reg_df.min_x.astype(int)
    reg_df.loc[reg_df.coefficient < 0, 'coefficient'] = 0
    reg_df.loc[reg_df.intercept < 0, 'intercept']     = 0
    reg_df = reg_df.sort_values(by='min_x')
    values = [reg_df.min_x] + [reg_df[key] for key in keys]
    return my_join(*values)

def predict(reg_df, x):
    for _, row in reg_df.iterrows():
        if row.min_x <= x < row.max_x:
            return row.coefficient * x + row.intercept
    return None

def pretty_plot(df, reg, title, n=10000):
    df = df.sample(n=min(n, len(df))).copy()
    df['pred'] = df.apply(lambda row: reg.predict(row['msg_size']), axis=1)
    df['group'] = 0
    for i, bp in enumerate(reg.breakpoints):
        df.loc[df['msg_size'] > bp, 'group'] = i + 1
    plot = ggplot(df) + geom_point(aes(x='msg_size', y='duration', color='factor(group)'), alpha=0.5, show_legend=False)
    for bp in reg.breakpoints:
        plot += geom_vline(xintercept=bp, linetype='dashed', color='gray')
    bps = [float('-inf')] + list(reg.breakpoints) + [float('inf')]
    for bp1, bp2 in zip(bps[:-1], bps[1:]):
        segment_df = df[(df.msg_size > bp1) & (df.msg_size < bp2)]
        if not segment_df.empty:
            plot += geom_line(segment_df, aes(x='msg_size', y='pred'), size=1)
    return plot + theme_bw() + scale_x_log10() + scale_y_log10() + \
           xlab('Message size (bytes)') + ylab('Duration (seconds)') + ggtitle(title)

# ==========================================
# Execution Pipeline
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SMPI calibration analysis and flag generation')
    parser.add_argument('--data-dir',   type=str,   default='./data')
    parser.add_argument('--output-dir', type=str,   default='.')
    parser.add_argument('--bandwidth',  type=float, default=937e6 / 8,
                        help='Expected network bandwidth in bytes/s (default: 937Mbps/8)')
    parser.add_argument('--latency',    type=float, default=22.7e-6,
                        help='Expected network latency in seconds (default: 22.7µs)')
    args = parser.parse_args()

    expected_bandwidth = args.bandwidth
    expected_latency   = args.latency

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading experiment data...")
    experiment, semantic_breaks = load_experiment(args.data_dir)
    drop_thresh = max(semantic_breaks)
    print(f"Semantic breakpoints: {semantic_breaks}")

    # --- Pass 1: Auto-discover Breakpoints ---
    print("\nDiscovering natural breakpoints...")
    discovered_bps = set()
    for op_key in ['send', 'isend', 'recv', 'pingpong']:
        tmp = aggregate_dataframe(experiment[op_key])
        reg = pycewise.compute_regression(x=tmp.msg_size.values, y=tmp.duration.values, mode='log').auto_simplify()
        print(f"  {op_key}: {set(reg.breakpoints)}")
        discovered_bps.update(reg.breakpoints)

    # Replace manual visual-inspection removal with a statistical check:
    # drop breakpoints where fewer than 3 data points exist on either side.
    all_x = np.sort(np.unique(np.concatenate([
        aggregate_dataframe(experiment[k]).msg_size.values
        for k in ['send', 'isend', 'recv', 'pingpong']
    ])))
    filtered_bps = set()
    for bp in sorted(discovered_bps):
        left  = int(np.sum(all_x <= bp))
        right = int(np.sum(all_x > bp))
        if left >= 3 and right >= 3:
            filtered_bps.add(bp)
        else:
            print(f" -> Dropped breakpoint {bp} (left={left}, right={right} data points)")

    breakpoints = sorted(list(filtered_bps | set(semantic_breaks)))
    print(f"\nCandidate breakpoints: {breakpoints}")

    # Retain only breakpoints that leave >= 2 data points in every sub-range for
    # every operation. This ensures segment-aligned regressions (needed for the
    # segment-level transfer model) without requiring manual visual inspection.
    def valid_bps_for(x_vals, bps):
        valid = []
        last = x_vals.min()
        for bp in sorted(bps):
            if ((x_vals > last) & (x_vals <= bp)).sum() >= 2 and (x_vals > bp).sum() >= 2:
                valid.append(bp)
                last = bp
        return set(valid)

    per_op_valid = {
        op: valid_bps_for(aggregate_dataframe(experiment[op]).msg_size.values, breakpoints)
        for op in ['send', 'isend', 'recv', 'pingpong']
    }
    breakpoints = sorted(set.intersection(*per_op_valid.values()))
    print(f"Final breakpoints (valid for all ops): {breakpoints}")

    # --- Pass 2: Final Regressions with Unified Breakpoints ---
    print("\nComputing final regressions...")
    def compute_reg(op_key, label):
        tmp = aggregate_dataframe(experiment[op_key])
        reg = pycewise.compute_regression(x=tmp.msg_size.values, y=tmp.duration.values,
                                          mode='log', breakpoints=breakpoints)
        draw_plot(reg, output_dir / f"{label}_regression.png")
        return reg

    reg_send     = compute_reg('send',     'MPI_Send')
    reg_isend    = compute_reg('isend',    'MPI_Isend')
    reg_recv     = compute_reg('recv',     'MPI_Recv')
    reg_pingpong = compute_reg('pingpong', 'PingPong')

    reg_send_str  = regression_to_str(reg_send.to_pandas(),     drop_threshold=drop_thresh)
    reg_isend_str = regression_to_str(reg_isend.to_pandas())
    reg_recv_str  = regression_to_str(reg_recv.to_pandas(),     drop_threshold=drop_thresh)

    # --- Transfer Model (segment-level subtraction, matching reference notebook) ---
    print("\nComputing Transfer Model...")
    df_reg_recv     = reg_recv.to_pandas()
    df_reg_pingpong = reg_pingpong.to_pandas()

    df_reg_transfer = df_reg_pingpong[['coefficient', 'intercept', 'max_x', 'min_x']].copy()
    df_reg_transfer['coefficient'] /= 2
    df_reg_transfer['intercept']   /= 2

    # Both regressions use the same breakpoints, so segment counts are aligned.
    pp_below   = df_reg_transfer['max_x'] <= drop_thresh
    recv_below = df_reg_recv['max_x']     <= drop_thresh
    df_reg_transfer.loc[pp_below, 'coefficient'] -= df_reg_recv.loc[recv_below, 'coefficient'].values
    df_reg_transfer.loc[pp_below, 'intercept']   -= df_reg_recv.loc[recv_below, 'intercept'].values

    # Clamp non-physical values
    df_reg_transfer.loc[df_reg_transfer['coefficient'] <= 0, 'coefficient'] = 1 / expected_bandwidth
    df_reg_transfer.loc[df_reg_transfer['intercept']   <  0, 'intercept']   = 0

    df_reg_transfer['bandwidth_factor'] = 1 / (df_reg_transfer.coefficient * expected_bandwidth)
    df_reg_transfer['latency_factor']   = df_reg_transfer.intercept / expected_latency

    bad_bw = df_reg_transfer['bandwidth_factor'].isna() | np.isinf(df_reg_transfer['bandwidth_factor'])
    if bad_bw.any():
        print(f" -> Clamped {bad_bw.sum()} unphysical bandwidth_factor value(s) to 1.0")
        df_reg_transfer.loc[bad_bw, 'bandwidth_factor'] = 1.0

    # Draw synthetic transfer model (graphical only, as in reference notebook)
    tmp_transfer = pd.DataFrame({'x': [1.1**i for i in range(220)]})
    tmp_transfer['y'] = tmp_transfer['x'].apply(lambda x: predict(df_reg_transfer, x))
    tmp_transfer = tmp_transfer.dropna()
    if len(tmp_transfer) >= 2:
        reg_transfer_plot = pycewise.compute_regression(
            x=tmp_transfer.x.values, y=tmp_transfer.y.values, breakpoints=breakpoints)
        draw_plot(reg_transfer_plot, output_dir / "Synthetic_Transfer_regression.png")

    bandwidth_str = regression_to_str(df_reg_transfer, drop_threshold=drop_thresh, keys=['bandwidth_factor'])
    latency_str   = regression_to_str(df_reg_transfer, drop_threshold=drop_thresh, keys=['latency_factor'])

    # --- Constants ---
    print("\nCalculating constants...")
    wtime_mean  = experiment['wtime'].duration.mean()
    test_mean   = experiment['test'].duration.mean()
    iprobe_mean = experiment['iprobe'].duration.mean()

    for op, data in [('wtime', experiment['wtime']), ('test', experiment['test']), ('iprobe', experiment['iprobe'])]:
        data[['duration']].boxplot()
        plt.savefig(output_dir / f"{op}_boxplot.png")
        plt.clf()

    # --- Output ---
    cli_flags = (
        f'--cfg=smpi/os:"{reg_send_str}" \\\n'
        f'--cfg=smpi/or:"{reg_recv_str}" \\\n'
        f'--cfg=smpi/ois:"{reg_isend_str}" \\\n'
        f'--cfg=smpi/bw-factor:"{bandwidth_str}" \\\n'
        f'--cfg=smpi/lat-factor:"{latency_str}" \\\n'
        f'--cfg=smpi/async-small-thresh:{semantic_breaks[0]} \\\n'
        f'--cfg=smpi/send-is-detached-thresh:{semantic_breaks[1]} \\\n'
        f'--cfg=smpi/iprobe:{iprobe_mean} \\\n'
        f'--cfg=smpi/test:{test_mean} \\\n'
        f'--cfg=smpi/host-speed:1'
    )

    flags_file = output_dir / 'smpi_flags.txt'
    with open(flags_file, 'w') as f:
        f.write(cli_flags)

    print("\n=========================================")
    print(f" READY TO USE FLAGS (Saved to {flags_file})")
    print("=========================================\n")
    print(cli_flags)

    print("Generating PDF plots...")
    pretty_plot(experiment['recv'], reg_recv, "Durations of MPI_Recv").save(output_dir / 'recv.pdf')
    pretty_plot(experiment['send'], reg_send, "Durations of MPI_Send").save(output_dir / 'send.pdf')

    print("\nCalibration Complete!")
