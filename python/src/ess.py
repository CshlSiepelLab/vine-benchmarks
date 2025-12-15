#!/usr/bin/env python3

import argparse
import numpy as np


def read_log_file(logfile, param):
    """
    Read a tab-delimited log file and return a 
    list of floats for the given parameter column.
    """
    param_index = None
    trace_values = []
    with open(logfile, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("[ID:"):
                continue
            fields = line.split()
            # Detect header (first non-comment line with at least 2 columns)
            if param_index is None:
                if param in fields:
                    param_index = fields.index(param)
                    continue
                else:
                    # This line must be header, but param is not in it
                    raise ValueError(f"Parameter '{param}' not found in header.")
            # Data line
            if param_index is not None:
                val = float(fields[param_index])
                trace_values.append(val)
    return trace_values


def calculate_running_ess(trace, burnin=0.1, max_lag=2000):
    """
    Compute running ESS of the trace, similar to BEAST/Tracer style.
    """
    trace = np.asarray(trace, dtype=float)
    
    ess_values = []
    
    for i in range(len(trace)):
        start_idx = int(burnin * (i + 1))
        subset = trace[start_idx : i + 1]  # Apply burn-in and keep subsetting for more samples each iteration
        n = len(subset) # Current samples subset size
        
        if n < 2:
            ess_values.append(0.0)
            continue
        
        centered = subset - np.mean(subset) # Autocorrelation and autocovariance must be computed on mean-zero data
        lag_limit = min(max_lag, n - 1) # Cannot compute lag > number of available samples
        
        acov = np.correlate(centered, centered, mode="full")    # Cross-correlation, yielding autocovariance for center data
        acov = acov[n - 1 : n - 1 + lag_limit + 1]  # Extract only lags 0,1,...,lag_limit
        if acov[0] == 0:
            ess_values.append(float(n))  # Each sample is effectively independent, so ess = n
            continue
        acf = acov / acov[0]  # Normalize autocovariance to get autocorrelation
        
        # Determine cutoff lag by Geyer’s initial positive sequence rule, removing noisy high-lag correlations
        cutoff = lag_limit
        for k in range(2, lag_limit + 1):
            if acf[k - 1] + acf[k] <= 0:
                cutoff = k
                break
        
        # Calculate autocorrelation time (ACT) using even lags up to cutoff
        act = acf[0]    # Start with lag 0 term
        for k in range(2, cutoff + 1, 2):   # Even lags only
            act += 2.0 * (acf[k - 1] + acf[k])
        
        if act < 1.0:   # ACT < 1 makes no sense (would imply negative correlation strong enough to increase ESS above n)
            act = 1.0
            
        ess = n / act
        ess_values.append(ess)
        
    return ess_values


def write_pruned_log_file(input_logfile, output_logfile, sample_cutoff_idx):
    """
    Write a pruned log file after burn-in.
    """
    with open(input_logfile, "r") as infile, open(output_logfile, "w") as outfile:
        line_count = 0
        for line in infile:
            line = line.strip()
            if line.startswith("#"):
                outfile.write(line + "\n")
                continue
            elif line_count <= sample_cutoff_idx:
                outfile.write(line + "\n")
            line_count += 1


def plot_ess_values(ess_values, ess_threshold, output_pdf, parameters=None, individual_parameter_ess=None):
    """
    Plot ESS values over samples and save to PDF.
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    fs=18
    plt.plot(range(1, len(ess_values) + 1), ess_values, label="Mean ESS", color="black")
    plt.axhline(y=ess_threshold, color="red", linestyle="-", label=f"ESS Threshold ({ess_threshold})")
    if individual_parameter_ess:
        colors = mpl.colormaps["tab20"].resampled(20).colors.tolist() + mpl.colormaps["tab20b"].resampled(20).colors.tolist() + mpl.colormaps["tab20c"].resampled(20).colors.tolist()
        for i, (param, ess_vals) in enumerate(individual_parameter_ess.items()):
            plt.plot(range(1, len(ess_vals) + 1), ess_vals, linestyle="--", label=f"ESS for {param}", color=colors[i])
    plt.xlabel("Number of Samples", fontsize=fs)
    plt.ylabel("Effective Sample Size (ESS)", fontsize=fs)
    plt.xticks(fontsize=fs-4)
    plt.yticks(fontsize=fs-4)
    if parameters:
        plt.title(f"Mean running ESS over samples for\n{', '.join(parameters)}", fontsize=fs-4)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=fs-4)
    plt.grid()
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close()


def parse_csv_string(s):
    if "," not in s:
        return [s.strip()]
    return [item.strip() for item in s.split(",")]

def main():
    parser = argparse.ArgumentParser(description="Compute running ESS for a BEAST2 log file (pure Python).")
    parser.add_argument("--logfile", help="Input TSV style MCMC .log file")
    parser.add_argument("--parameter", type=parse_csv_string, help="Parameter(s) (column name(s)) to compute ESS for. Format as csv string. If multiple parameters are provided, the threshold will be based on the average ESS across them.",)
    parser.add_argument("--output_logfile", help="Output pruned log file after burn-in")
    parser.add_argument("--ess_threshold", type=float, help="ESS threshold for cutoff (default: 200)", default=200.0)
    parser.add_argument("--burnin", type=float, help="Burn-in proportion (default: 0.1)", default=0.1)
    parser.add_argument("--plot_ess", action="store_true", help="Whether to plot ESS values over samples")
    args = parser.parse_args()

    logfile = args.logfile
    parameter = args.parameter
    output_logfile = args.output_logfile
    ess_threshold = args.ess_threshold
    burnin = args.burnin
    plot_ess = args.plot_ess

    all_ess_values = {}
    for param in parameter:
        trace = read_log_file(logfile, param)
        all_ess_values[param] = calculate_running_ess(trace, burnin=burnin)

    # Compute average ESS across parameters for each sample index
    num_samples = len(all_ess_values[parameter[0]])
    ess_values = []
    for i in range(num_samples):
        avg_ess = np.mean([all_ess_values[param][i] for param in parameter])
        ess_values.append(avg_ess)

    sample_end = next((i for i, ess in enumerate(ess_values) if ess >= ess_threshold), len(ess_values) - 1) + 1
    write_pruned_log_file(logfile, output_logfile, sample_end)
    if plot_ess:
        plot_ess_values(ess_values, ess_threshold, output_logfile.replace(".log", ".pdf"), parameters=parameter, individual_parameter_ess=all_ess_values)
        
        
if __name__ == "__main__":
    main()

