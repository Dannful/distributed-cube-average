#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import io
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pycewise
import plotnine
from plotnine import *
from pathlib import Path
import argparse

print(f'pycewise version: {pycewise.__version__}')

# ==========================================
# Configuration & Constants
# ==========================================
expected_bandwidth = 937e6 / 8
expected_latency   = 22.7e-6

expected_loopback_bandwidth = 2e11 / 8
expected_loopback_latency = 1e-7

# Plotnine configuration
plotnine.options.figure_size = (12, 8)

# ==========================================
# Data Extraction & Aggregation Functions
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
    # Sort by start time and drop first measurement per (op, msg_size) group (warmup)
    df_sorted = dataframe.sort_values('start')
    df_filtered = df_sorted.groupby(['op', 'msg_size'], sort=False).apply(
        lambda g: g.iloc[1:] if len(g) > 2 else g, include_groups=False
    ).reset_index(level=['op', 'msg_size'])

    # Use median for robust aggregation (resilient to outliers)
    numeric_cols = df_filtered.select_dtypes(include='number').columns.drop('msg_size', errors='ignore')
    df = df_filtered.groupby('msg_size')[numeric_cols].median().reset_index()
    for col in dataframe.columns:
        if col not in df.columns:
            value = list(dataframe[col].unique())
            if len(value) == 1:
                df[col] = value[0]
    return df

def filter_pathological_sizes(raw_df, expected_lat, expected_bw, multiplier=20):
    """Remove message sizes whose median duration is >multiplier× the physical expectation.

    Catches systematic TCP/MPI pathologies at specific sizes (e.g., a size that
    always takes 1.4ms when the expected time is 23µs). The physics-based reference
    is expected_latency + size/expected_bandwidth, which is the minimum possible time.
    """
    medians = raw_df.groupby('msg_size')['duration'].median()
    expected = medians.index.to_series().map(lambda s: expected_lat + s / expected_bw)
    ratio = medians / expected
    bad_sizes = ratio[ratio > multiplier].index.tolist()
    if bad_sizes:
        print(f" -> Removed {len(bad_sizes)} pathological size(s): {bad_sizes} "
              f"(>{multiplier}× expected physical time)")
    return raw_df[~raw_df['msg_size'].isin(bad_sizes)]

def load_experiment(dir_name='.'):
    result = extract_local_dir(dir_name)

    breakpoints_key = next((k for k in result.keys() if k.endswith('breakpoints')), None)
    if breakpoints_key:
        print(f"Found breakpoints file at: {breakpoints_key}")
        raw_lines = result[breakpoints_key].decode().strip().split('\n')
        semantic_breaks = [int(line.strip()) for line in raw_lines if line.strip()]

        semantic_breaks.sort()
        
        if len(semantic_breaks) == 1:
            semantic_breaks = [semantic_breaks[0], semantic_breaks[0]]
        elif len(semantic_breaks) >= 2:
            semantic_breaks = semantic_breaks[-2:]
        else:
            default_threshold = 65536
            print(f"WARNING: 'breakpoints' file was empty. Defaulting semantic breaks to [{default_threshold}].")
            semantic_breaks = [default_threshold, default_threshold]
    else:
        default_threshold = 65536 
        print(f"WARNING: 'breakpoints' missing. Defaulting semantic breaks to [{default_threshold}].")
        semantic_breaks = [default_threshold, default_threshold]

    pingpong_key = next((k for k in result.keys() if k.endswith('exp_PingPong.csv')), None)
    if not pingpong_key:
        raise FileNotFoundError("Could not find '*exp_PingPong.csv'.")

    df_pingpong = result[pingpong_key]
    pingpong_send = df_pingpong[df_pingpong.op=='MPI_Send'].reset_index(drop=True)
    pingpong_recv = df_pingpong[df_pingpong.op=='MPI_Recv'].reset_index(drop=True)

    # Detect and correct the polling artifact: for small messages, the calibrate
    # tool's PingPong Recv duration is a constant ~10ms (iteration interval), not
    # actual network time. For these sizes, 2*Send_duration is a better estimate
    # of the round-trip, since Send blocks until the reply arrives.
    recv_medians = pingpong_recv.groupby('msg_size')['duration'].median()
    baseline_recv = recv_medians.iloc[:3].median()  # median of smallest sizes
    artifact_threshold = baseline_recv * 2.0

    corrected_duration = []
    for _, (send_dur, recv_dur, size) in enumerate(
        zip(pingpong_send.duration, pingpong_recv.duration, pingpong_send.msg_size)
    ):
        median_recv_for_size = recv_medians.loc[size]
        if median_recv_for_size < artifact_threshold and recv_dur > send_dur * 100:
            corrected_duration.append(2 * send_dur)
        else:
            corrected_duration.append(send_dur + recv_dur)

    df_pingpong_combined = pd.DataFrame(dict(
        op       = 'PingPong',
        msg_size = pingpong_send.msg_size,
        start    = pingpong_send.start,
        duration = corrected_duration
    ))

    n_corrected = sum(1 for d, s in zip(corrected_duration, pingpong_send.duration + pingpong_recv.duration)
                      if abs(d - s) > 1e-12)
    print(f"PingPong artifact correction: {n_corrected}/{len(corrected_duration)} rows corrected "
          f"(baseline Recv={baseline_recv:.4e}s, threshold={artifact_threshold:.4e}s)")

    def get_csv(name):
        return next((result[k] for k in result.keys() if k.endswith(name)), None)

    final_result = {
        'pingpong': df_pingpong_combined,
        'send'    : pingpong_send,
        'isend'   : get_csv('exp_Isend.csv'),
        'recv'    : get_csv('exp_Recv.csv'),
        'wtime'   : get_csv('exp_Wtime.csv'),
        'test'    : get_csv('exp_Test.csv'),
        'iprobe'  : get_csv('exp_Iprobe.csv'),
    }
    return final_result, semantic_breaks

# ==========================================
# Regression & Formatting Functions
# ==========================================
def draw_plot(reg, filename, alpha=1):
    reg.plot_dataset(log=True, alpha=alpha)
    plt.xlabel('Message size (bytes)')
    plt.ylabel('Duration (seconds)')
    plt.savefig(filename)
    plt.clf()

def safe_compute_regression(raw_df, op_name, original_breakpoints, output_dir=Path('.')):
    print(f"\n=== Processing: {op_name} ===")
    tmp = aggregate_dataframe(raw_df)
    x_vals = tmp.msg_size.values
    y_vals = tmp.duration.values

    valid_bps = []
    last_bp = x_vals.min()
    for bp in sorted(original_breakpoints):
        if ((x_vals > last_bp) & (x_vals <= bp)).sum() >= 2 and (x_vals > bp).sum() >= 2:
            valid_bps.append(bp)
            last_bp = bp

    print(f" -> Validated breakpoints: {valid_bps}")
    reg_model = pycewise.compute_regression(x=x_vals, y=y_vals, mode='log', breakpoints=valid_bps)
    draw_plot(reg_model, output_dir / f"{op_name}_regression.png")
    return reg_model

def my_join(*values, inner_sep=':', outer_sep=';'):
    return outer_sep.join([inner_sep.join([str(it) for it in items]) for items in zip(*values)])

def regression_to_str(reg_df, drop_threshold=None, keys=['intercept', 'coefficient']):
    if drop_threshold is not None:
        reg_df.loc[reg_df.min_x >= drop_threshold, 'intercept'] = 0
        reg_df.loc[reg_df.min_x >= drop_threshold, 'coefficient'] = 0
    reg_df.loc[reg_df.min_x < 0, 'min_x'] = 0
    reg_df.min_x = reg_df.min_x.astype(int)
    reg_df.loc[reg_df.coefficient < 0, 'coefficient'] = 0
    reg_df.loc[reg_df.intercept < 0, 'intercept'] = 0
    reg_df = reg_df.sort_values(by='min_x')
    
    values = [reg_df.min_x] + [reg_df[key] for key in keys]
    return my_join(*values)

def predict(reg, x):
    for _, row in reg.iterrows():
        if row.min_x <= x < row.max_x:
            return row.coefficient * x + row.intercept
    return 0

def pretty_plot(df, reg, title, n=10000):
    df = df.sample(n=min(n, len(df))).copy()
    df['pred'] = df.apply(lambda row: reg.predict(row['msg_size']), axis=1)
    df['group'] = 0
    for i, bp in enumerate(reg.breakpoints):
        df.loc[df['msg_size'] > bp, 'group'] = i+1
        
    plot = ggplot(df) + geom_point(aes(x='msg_size', y='duration', color='factor(group)'), alpha=0.5, show_legend=False)
    
    for bp in reg.breakpoints:
        plot += geom_vline(xintercept=bp, linetype='dashed', color='gray')
        
    bps = [float('-inf')] + list(reg.breakpoints) + [float('inf')]
    for bp1, bp2 in zip(bps[:-1], bps[1:]):
        segment_df = df[(df.msg_size > bp1) & (df.msg_size < bp2)]
        if not segment_df.empty:
            plot += geom_line(segment_df, aes(x='msg_size', y='pred'), size=1)
            
    plot = plot + theme_bw() + scale_x_log10() + scale_y_log10() + \
           xlab('Message size (bytes)') + ylab('Duration (seconds)') + ggtitle(title)
    return plot

def validate_data_quality(experiment):
    """Warn about potential data quality issues."""
    warnings = []

    # Check PingPong duration grows with message size
    pp = aggregate_dataframe(experiment['pingpong'])
    if len(pp) > 5:
        corr = pp[['msg_size', 'duration']].corr().iloc[0, 1]
        if corr < 0.5:
            warnings.append(f"PingPong duration poorly correlated with msg_size (r={corr:.3f})")

    # Check standalone Recv monotonicity
    recv = aggregate_dataframe(experiment['recv']).sort_values('msg_size')
    if len(recv) > 5:
        violations = (recv.duration.diff().dropna() < 0).sum()
        if violations > len(recv) * 0.3:
            warnings.append(f"Recv duration non-monotonic: {violations}/{len(recv)-1} decreasing steps")

    # Check per-size variance
    for op_key in ['send', 'recv', 'isend']:
        df = experiment[op_key]
        if df is None:
            continue
        grouped = df.groupby('msg_size')['duration']
        for size, group in grouped:
            if len(group) < 3:
                continue
            cv = group.std() / group.mean() if group.mean() > 0 else 0
            if cv > 1.0:
                warnings.append(f"{op_key} msg_size={size}: CV={cv:.1f} (high variance)")

    # Check minimum data points
    for op_key in ['send', 'recv', 'isend', 'pingpong']:
        df = experiment[op_key]
        if df is None:
            continue
        sizes = df.groupby('msg_size').size()
        sparse = sizes[sizes < 5]
        if len(sparse) > 0:
            warnings.append(f"{op_key}: {len(sparse)} sizes have fewer than 5 data points")

    if warnings:
        print(f"\nDATA QUALITY WARNINGS ({len(warnings)}):")
        for w in warnings:
            print(f"  - {w}")
    else:
        print("\nData quality checks passed.")

# ==========================================
# Execution Pipeline
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='SMPI calibration analysis and flag generation')
    parser.add_argument('--data-dir', type=str, default='./data',
                        help='Directory containing experimental CSV files (default: ./data)')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='Directory for output files (default: current directory)')
    parser.add_argument('--bandwidth', type=float, default=937e6 / 8,
                        help='Expected network bandwidth in bytes/s (default: 937Mbps / 8)')
    parser.add_argument('--latency', type=float, default=22.7e-6,
                        help='Expected network latency in seconds (default: 22.7us)')
    args = parser.parse_args()

    expected_bandwidth = args.bandwidth
    expected_latency = args.latency

    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading experiment data...")
    experiment, semantic_breaks = load_experiment(args.data_dir)
    drop_thresh = max(semantic_breaks)

    validate_data_quality(experiment)

    # Pre-filter pathological message sizes from all operations
    print("\nFiltering pathological message sizes...")
    for op_key in ['send', 'recv', 'isend', 'pingpong']:
        if experiment[op_key] is not None:
            experiment[op_key] = filter_pathological_sizes(
                experiment[op_key], expected_latency, expected_bandwidth
            )

    # --- Pass 1: Auto-discover Breakpoints ---
    print("\nDiscovering natural breakpoints...")
    discovered_bps = set()
    for op_key in ['send', 'isend', 'recv', 'pingpong']:
        tmp = aggregate_dataframe(experiment[op_key])
        reg = pycewise.compute_regression(x=tmp.msg_size.values, y=tmp.duration.values, mode='log').auto_simplify()
        discovered_bps.update(reg.breakpoints)

    # Statistically filter spurious breakpoints: keep only those where both adjacent
    # segments have enough data and the slope changes meaningfully.
    all_x = np.sort(np.unique(np.concatenate([
        aggregate_dataframe(experiment[k]).msg_size.values for k in ['send', 'isend', 'recv', 'pingpong']
    ])))
    filtered_bps = set()
    for bp in sorted(discovered_bps):
        left_count = np.sum(all_x <= bp)
        right_count = np.sum(all_x > bp)
        if left_count >= 3 and right_count >= 3:
            filtered_bps.add(bp)
        else:
            print(f" -> Dropped breakpoint {bp} (left={left_count}, right={right_count} points)")

    breakpoints = sorted(list(filtered_bps | set(semantic_breaks)))

    # --- Pass 2: Final Regressions ---
    reg_send     = safe_compute_regression(experiment['send'], "MPI_Send", breakpoints, output_dir)
    reg_isend    = safe_compute_regression(experiment['isend'], "MPI_Isend", breakpoints, output_dir)
    reg_recv     = safe_compute_regression(experiment['recv'], "MPI_Recv", breakpoints, output_dir)
    reg_pingpong = safe_compute_regression(experiment['pingpong'], "PingPong", breakpoints, output_dir)

    reg_send_str  = regression_to_str(reg_send.to_pandas(), drop_threshold=drop_thresh)
    reg_isend_str = regression_to_str(reg_isend.to_pandas())
    reg_recv_str  = regression_to_str(reg_recv.to_pandas(), drop_threshold=drop_thresh)

    # --- Transfer Model Construction ---
    print("\nComputing Transfer Models...")
    df_reg_recv = reg_recv.to_pandas()
    df_reg_transfer = reg_pingpong.to_pandas()[['coefficient', 'intercept', 'max_x', 'min_x']].copy()
    
    df_reg_transfer['coefficient'] /= 2
    df_reg_transfer['intercept'] /= 2
    
    mask = df_reg_transfer['max_x'] <= drop_thresh
    df_reg_transfer.loc[mask, 'coefficient'] -= df_reg_recv[df_reg_recv['max_x'] <= drop_thresh]['coefficient'].values
    df_reg_transfer.loc[mask, 'intercept'] -= df_reg_recv[df_reg_recv['max_x'] <= drop_thresh]['intercept'].values
    
    # Clamp non-positive coefficients: coefficient=0 causes inf bw-factor, negative is unphysical
    df_reg_transfer.loc[df_reg_transfer['coefficient'] <= 0, 'coefficient'] = 1 / expected_bandwidth
    df_reg_transfer.loc[df_reg_transfer['intercept'] < 0, 'intercept'] = 0

    # Draw synthetic transfer plot with strict validation
    tmp_synth = pd.DataFrame({'x': [1.1**i for i in range(220)]})
    tmp_synth['y'] = tmp_synth.apply(lambda row: predict(df_reg_transfer, row.x), axis=1)
    
    synth_x_vals = tmp_synth.x.values
    synth_y_vals = tmp_synth.y.values
    
    valid_synth_bps = []
    last_bp = synth_x_vals.min()
    
    for bp in sorted(breakpoints):
        points_in_segment = ((synth_x_vals > last_bp) & (synth_x_vals <= bp)).sum()
        points_after = (synth_x_vals > bp).sum()
        
        if points_in_segment >= 2 and points_after >= 2:
            valid_synth_bps.append(bp)
            last_bp = bp
        else:
            print(f" -> Dropped synthetic breakpoint {bp} to preserve covariance math.")

    reg_synthetic = pycewise.compute_regression(
        x=synth_x_vals,
        y=synth_y_vals,
        breakpoints=valid_synth_bps
    )
    draw_plot(reg_synthetic, output_dir / "Synthetic_Transfer_regression.png")

    df_reg_transfer['bandwidth_factor'] = 1 / (df_reg_transfer.coefficient * expected_bandwidth)
    df_reg_transfer['latency_factor']   = df_reg_transfer.intercept / expected_latency

    # Clamp unphysical bw-factor values: inf/nan → 1.0 (no scaling = use platform bandwidth as-is)
    bad_bw = df_reg_transfer['bandwidth_factor'].isna() | np.isinf(df_reg_transfer['bandwidth_factor'])
    if bad_bw.any():
        print(f" -> Clamped {bad_bw.sum()} unphysical bandwidth_factor value(s) to 1.0")
        df_reg_transfer.loc[bad_bw, 'bandwidth_factor'] = 1.0

    bandwidth_str = regression_to_str(df_reg_transfer, drop_threshold=drop_thresh, keys=['bandwidth_factor'])
    latency_str   = regression_to_str(df_reg_transfer, drop_threshold=drop_thresh, keys=['latency_factor'])

    # --- Constant Means & Boxplots ---
    print("\nCalculating constants and generating boxplots...")
    wtime_mean  = experiment['wtime'].duration.mean()
    test_mean   = experiment['test'].duration.mean()
    iprobe_mean = experiment['iprobe'].duration.mean()

    for op, data in [('wtime', experiment['wtime']), ('test', experiment['test']), ('iprobe', experiment['iprobe'])]:
        data[['duration']].boxplot()
        plt.savefig(output_dir / f"{op}_boxplot.png")
        plt.clf()

    # --- Generate SimGrid Command-Line Arguments ---
    print("\nGenerating SMPI Command-Line Flags...")
    
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

    # --- Save High-Quality Plotnine Visuals ---
    print("Generating High-Quality PDF plots...")
    plot_recv = pretty_plot(experiment['recv'], reg_recv, "Durations of MPI_Recv")
    plot_recv.save(output_dir / 'recv.pdf')

    plot_send = pretty_plot(experiment['send'], reg_send, "Durations of MPI_Send")
    plot_send.save(output_dir / 'send.pdf')

    print("\nCalibration Complete!")
    
    # Raw Bandwidth (bps) and Latency (s) printed at the very end
    print(expected_bandwidth * 8)
    print(expected_latency)
