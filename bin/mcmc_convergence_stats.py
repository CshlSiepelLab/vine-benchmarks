import argparse
import os
import numpy as np
from scipy.stats import norm
import matplotlib as mpl
mpl.use("Agg")   # Use non-interactive backend
import matplotlib.backends.backend_pdf  # For pyinstaller inclusion
import matplotlib.pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


def read_log_file(logfile, param):
    """
    Read a tab-delimited log file and return an np.array of floats for the given parameter column.
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
    trace_values_array = np.asarray(trace_values, dtype=float)
    return trace_values_array


def calculate_one_chain_ess_and_mcse(subset, lag_limit=2000):
    """ 
    Calculate Effective Sample Size (ESS) and Monte Carlo Standard Error (MCSE).
    """
    n = len(subset) # Current samples subset size
    if n < 2:
        return np.nan, np.nan  # Cannot compute ESS with one sample
    
    max_lag = min(lag_limit, n - 1) # Cannot compute lag > number of available samples
    
    centered = subset - np.mean(subset) # Act must be computed on mean-zero data
    
    # My version
    # Not exact due to rounding and the order of operations, but fast
    full = np.correlate(centered, centered, mode="full")
    gamma = full[n - 1 : n - 1 + max_lag + 1] / np.arange(n, n - max_lag - 1, -1)
    cutoff = max_lag
    for k in range(2, max_lag + 1, 2):  # even lags only
        if gamma[k-1] + gamma[k] <= 0:  # Determine cutoff lag by Geyer’s initial positive sequence rule
            cutoff = k
            break
    varStat = gamma[0]
    for k in range(2, cutoff + 1, 2):
        varStat += 2.0 * (gamma[k-1] + gamma[k])
    
    # # Direct translation of Tracer java version from https://github.com/beast-dev/beast-mcmc/blob/main/src/dr/inference/trace/TraceCorrelation.java
    # # Exact, but slow
    # mean = np.mean(subset)
    # gamma = np.zeros(max_lag)
    # for lag in range(max_lag):
    #     s = 0.0
    #     for j in range(n - lag):
    #         s += (centered[j] * centered[j + lag])
    #     gamma[lag] = s / float(n - lag)
    #     if lag == 0:
    #         varStat = gamma[0]
    #     elif lag % 2 == 0:  # even lags only
    #         if gamma[lag - 1] + gamma[lag] > 0:
    #             varStat += 2.0 * (gamma[lag - 1] + gamma[lag])
    #         else:
    #             current_max_lag = lag
    #             break
        
    if gamma[0] == 0:
        act = 0
        ess = 1.0
    else:
        act = varStat / gamma[0]
        if act <= 0:
            ess = float(n)
        else:
            ess = n / act
    
    mcse = np.std(subset, ddof=1) / np.sqrt(ess)
    
    return ess, mcse


def calculate_one_chain_running_stats(trace, burnin=0.1, calc_freq=1000):
    """
    Compute running stats of one chain.
    """
    ess_values = []
    mcse_values = []
    chain_len = len(trace)
    
    for i in range(0, chain_len, calc_freq):
        print(i)
        start_idx = int(burnin * (i + 1))
        subset = trace[start_idx : i + 1]  # Apply burn-in and keep subsetting for more samples each iteration
        ess, mcse = calculate_one_chain_ess_and_mcse(subset)
        ess_values.append(ess)
        mcse_values.append(mcse)
        
    return ess_values, mcse_values


def rank_normalize_chains(chains):
    """
    Apply Vehtari et al. 2020 and Blom 1958 rank-normalization.
    Preserves original chain structure.
    """
    # Flatten all chains into one vector while remembering boundaries
    flat = np.concatenate(chains)
    S = len(flat)
    
    # Compute pooled ranks
    order = np.argsort(flat)
    ranks = np.empty(S, dtype=float)
    ranks[order] = np.arange(1, S+1)
    
    # Transform ranks to normal scores using inverse transform and fractional offset
    z = norm.ppf((ranks - 3/8) / (S + 1/4)) # Blom 1958 formula; Vehtari 2020 uses (- 1/4) in denominator, but neither influences results significantly
    
    # Restore original chains
    out = []
    idx = 0
    for c in chains:
        n = len(c)
        out.append(z[idx:idx+n])
        idx += n
        
    return out


def calculate_rhat(chains, split=True, rank_normalize=True):
    """
    Calculate R-hat (optionally split and/or rank normalized) for multiple chains.
    Calculations taken from Vehtari et. al. 2020 paper at https://sites.stat.columbia.edu/gelman/research/published/Vehtari_etal_2020_rhat_ess.pdf
    Assumes burnin has already been applied to the input.
    """
    # Optionally split chains in half
    if split:
        split_chains = []
        for chain in chains:
            m = len(chain) // 2
            if m < 2:
                return np.nan
            split_chains.append(chain[:m])
            split_chains.append(chain[m:])
    else:
        split_chains = chains
    
    # Number of chains
    M = len(split_chains)
    
    # Ensure a common length for all chains
    N = min(chain.shape[0] for chain in split_chains)
    split_chains = [c[:N] for c in split_chains]
    
    # Optionally rank-normalize chains
    if rank_normalize:
        split_chains = rank_normalize_chains(split_chains)
    
    # Calculate between chains variance
    chain_means = np.array([np.mean(chain) for chain in split_chains])
    overall_mean = np.mean(chain_means)
    B = ((N) / (M - 1)) * np.sum((chain_means - overall_mean) ** 2)
    
    # Calculate within chain variance
    W = np.mean([np.var(chain, ddof=1) for chain in split_chains])
    
    var_hat = (((N-1)/N)*W) + (B/N)
    r_hat = np.sqrt(var_hat / W)
    return r_hat


def calculate_multiple_chain_running_stats(chains, burnin=0.1, calc_freq=1000):
    """
    Compute running stats of multiple matched chains.
    """
    rhat_values = []
    chain_len = len(chains[0])
    for i in range(0, chain_len, calc_freq):
        start_idx = int(burnin * (i + 1))
        subset = [chain[start_idx : i + 1] for chain in chains]  # Apply burn-in and keep subsetting for more samples each iteration
        if len(subset[0]) < 2:
            rhat = np.nan
        else:
            rhat = calculate_rhat(subset)
        rhat_values.append(rhat)
        
    return rhat_values


def plot_values(mean_values, median_values, threshold, crossing_point, calc_freq, raw_values, output_pdf, ylabel):
    """
    Plot ESS or MCSE values over samples and save to PDF.
    """
    plt.figure(figsize=(10, 6))
    fs = 18
    # Get x-axis points
    n = len(raw_values[list(raw_values.keys())[0]]) # Raw values must be provided; Assumes all data have same mcmc sample length
    if calc_freq is None:
        x = lambda n: range(1, n + 1)
    else:
        x = lambda n: np.arange(1, n + 1) * calc_freq
    # Plot mean line
    if mean_values is not None:
        plt.plot(x, mean_values, label="Mean", color="black", linestyle="-", linewidth=4)
    # Plot median line
    if median_values is not None:
        plt.plot(x, median_values, label="Median", color="grey", linestyle="-", linewidth=4)
    # Plot threshold line
    if threshold is not None:
        plt.axhline(y=threshold, color="red", linestyle="-", label=f"Threshold ({threshold})", linewidth=2)
    # Plot specific ESS values
    colors = mpl.colormaps["tab20"].resampled(20).colors.tolist() + \
             mpl.colormaps["tab20b"].resampled(20).colors.tolist() + \
             mpl.colormaps["tab20c"].resampled(20).colors.tolist()
    for i, (id, values) in enumerate(raw_values.items()):
        plt.plot(x, values, linestyle="--", label=f"{id}", color=colors[i], linewidth=1)
    # Configure plot appearance
    plt.xlabel("MCMC sample number", fontsize=fs)
    plt.ylabel(ylabel, fontsize=fs)
    plt.xticks(fontsize=fs - 4)
    plt.yticks(fontsize=fs - 4)
    if calc_freq is not None:
        crossing_point_adjusted = crossing_point * calc_freq
    else:
        crossing_point_adjusted = crossing_point
    plt.title(f"All lines cross threshold at sample {crossing_point_adjusted}", fontsize=fs - 4)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=fs - 4)
    plt.grid()
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close()


def parse_csv_string(s):
    if "," not in s:
        return [s.strip()]
    return [item.strip() for item in s.split(",")]

parser = argparse.ArgumentParser(description="Compute running ESS for a BEAST2 log file (pure Python).")
parser.add_argument("--logfiles", type=parse_csv_string, help="Input TSV style MCMC .log files, which must all match in chain length.")
parser.add_argument("--parameters", type=parse_csv_string, help="Parameter(s) (column name(s)) to compute ESS for. Format as csv string. If multiple parameters are provided, the threshold will be based on the average ESS across them.",)
parser.add_argument("--outputprefix", help="Output prefix for stat files.")
parser.add_argument("--burnin", type=float, help="Burn-in proportion for ESS and MCSE calculations (default: 0.1)", default=0.1)
parser.add_argument("--rhat_burnin", type=float, help="Burn-in proportion for R-hat calculation (default: 0.0)", default=0.0)
parser.add_argument("--calc_freq", type=int, help="Frequency of ESS/MCSE/R-hat calculations (default: 1000)", default=1000)
args = parser.parse_args()

logfiles = args.logfiles
parameters = args.parameters
outputprefix = args.outputprefix
burnin = args.burnin
rhat_burnin = args.rhat_burnin
calc_freq = args.calc_freq

all_ess_datasets = {}
all_mcse_datasets = {}
mean_ess_datasets = {}
mean_mcse_datasets = {}

for logfile in logfiles:
    id = os.path.basename(logfile)
    id = f"{os.path.basename(os.path.dirname(logfile))}_{os.path.basename(logfile)}"    # If basename is not unique enough
    all_ess_datasets[id] = {}
    all_mcse_datasets[id] = {}
    for param in parameters:
        trace = read_log_file(logfile, param)
        trace = trace[0:150000]
        all_ess_datasets[id][param], all_mcse_datasets[id][param] = calculate_one_chain_running_stats(trace, burnin=burnin, calc_freq=calc_freq)
        
    # Summarize across parameters for this dataset
    ess_arr = np.array([all_ess_datasets[id][param] for param in parameters])
    mcse_arr = np.array([all_mcse_datasets[id][param] for param in parameters])
    mean_ess_datasets[id]  = ess_arr.mean(axis=0)
    mean_mcse_datasets[id] = mcse_arr.mean(axis=0)

# Summarize across datasets
mean_ess_arr = np.array(list(mean_ess_datasets.values()))
mean_mcse_arr = np.array(list(mean_mcse_datasets.values()))
mean_ess_values = mean_ess_arr.mean(axis=0)
mean_mcse_values = mean_mcse_arr.mean(axis=0)
median_ess_values = np.median(mean_ess_arr, axis=0)
median_mcse_values = np.median(mean_mcse_arr, axis=0)

# Output raw data values
with open(f"{outputprefix}_ess_mcse_raw.tsv", "w") as f:
    param_header = "\t".join([f"ess_{param}\tmcse_{param}" for param in parameters]) + "\tess_mean\tmcse_mean"
    f.write(f"logfile_id\tsample_index\t{param_header}\n")
    # Write per-dataset values
    for dataset_id in all_ess_datasets.keys():
        num_samples = len(mean_ess_datasets[dataset_id])
        for i in range(num_samples):
            ess_mcse_str = ""
            for param in parameters:
                ess_val = all_ess_datasets[dataset_id][param][i]
                mcse_val = all_mcse_datasets[dataset_id][param][i]
                ess_mcse_str += f"{ess_val}\t{mcse_val}\t"
            ess_mean_val = mean_ess_datasets[dataset_id][i]
            mcse_mean_val = mean_mcse_datasets[dataset_id][i]
            ess_mcse_str += f"{ess_mean_val}\t{mcse_mean_val}"
            f.write(f"{dataset_id}\t{i + 1}\t{ess_mcse_str}\n")

# Output mean and median ESS/MCSE values across datasets
with open(f"{outputprefix}_ess_mcse_summary.tsv", "w") as f:
    f.write("sample_index\tmean_ess\tmedian_ess\tmean_mcse\tmedian_mcse\n")
    num_samples = len(mean_ess_values)
    for i in range(num_samples):
        f.write(f"{i + 1}\t{mean_ess_values[i]}\t{median_ess_values[i]}\t{mean_mcse_values[i]}\t{median_mcse_values[i]}\n")

# Calculate between chain stats
if len(logfiles) > 1:
    all_r_hats = {}
    for param in parameters:
        param_traces = []
        for logfile in logfiles:
            trace = read_log_file(logfile, param)
            param_traces.append(trace)
        all_r_hats[param] = calculate_multiple_chain_running_stats(param_traces, burnin=rhat_burnin, calc_freq=calc_freq)

# Output R-hat values
with open(f"{outputprefix}_rhat_raw.tsv", "w") as f:
    param_header = "\t".join([f"rhat_{param}" for param in parameters])
    f.write(f"sample_index\t{param_header}\n")
    num_samples = len(all_r_hats[parameters[0]])
    for i in range(num_samples):
        rhat_str = ""
        for param in parameters:
            rhat_str += f"{all_r_hats[param][i]}\t"
        f.write(f"{i + 1}\t{rhat_str.strip()}\n")
        

# Find the convergence point
ess_threshold = 400.0
mcse_threshold = 0.05
rhat_threshold = 1.01

all_ess_datasets_flat = {f"{dataset}_{param}": all_ess_datasets[dataset][param] for dataset in all_ess_datasets.keys() for param in parameters}
crossing_point_ess = next((i for i in range(len(mean_ess_values)) if all(ess_values[i] >= ess_threshold for ess_values in all_ess_datasets_flat.values())), len(mean_ess_values) - 1) + 1

all_mcse_datasets_flat = {f"{dataset}_{param}": all_mcse_datasets[dataset][param] for dataset in all_mcse_datasets.keys() for param in parameters}
crossing_point_mcse = next((i for i in range(len(mean_mcse_values)) if all(mcse_values[i] <= mcse_threshold for mcse_values in all_mcse_datasets_flat.values())), len(mean_mcse_values) - 1) + 1
# crossing_point_mcse = 1 # Ignore MCSE threshold for now, as it is often not met

num_samples = len(mean_ess_values)

if len(logfiles) > 1:
    crossing_point_rhat = next((i for i in range(num_samples) if all(all_r_hats[param][i] <= rhat_threshold for param in parameters)), len(mean_ess_values) - 1) + 1
    crossing_point = max(crossing_point_ess, crossing_point_mcse, crossing_point_rhat)   # Max of ESS, MCSE, and R-hat crossing points to ensure all criteria are met
else:
    crossing_point = max(crossing_point_ess, crossing_point_mcse)  # Max of ESS and MCSE crossing points

scale_factor = crossing_point / num_samples
with open(f"{outputprefix}_convergence_scale_factor.txt", "w") as f:
    f.write(f"{scale_factor}\n")


# Plot ESS values
# plot_values(mean_ess_values, median_ess_values, ess_threshold, crossing_point_ess, all_ess_datasets_flat, output_pdf=f"{outputprefix}_ess.pdf", ylabel="Effective Sample Size (ESS)")
plot_values(None, None, ess_threshold, crossing_point_ess, calc_freq, all_ess_datasets_flat, output_pdf=f"{outputprefix}_ess.pdf", ylabel="Effective Sample Size (ESS)")

# Plot MCSE values (ignoring threshold and crossing point for now)
# plot_values(mean_mcse_values, median_mcse_values, mcse_threshold, crossing_point_mcse, all_mcse_datasets_flat, output_pdf=f"{outputprefix}_mcse.pdf", ylabel="Monte Carlo Standard Error (MCSE)")
plot_values(None, None, mcse_threshold, crossing_point_mcse, calc_freq, all_mcse_datasets_flat, output_pdf=f"{outputprefix}_mcse.pdf", ylabel="Monte Carlo Standard Error (MCSE)")

# Plot R-hat values if multiple chains
if len(logfiles) > 1:
    all_rhat_datasets_flat = {f"{param}": all_r_hats[param] for param in parameters}
    # mean_rhat_values = np.array([np.mean([all_r_hats[param][i] for param in parameters]) for i in range(num_samples)])
    # median_rhat_values = np.array([np.median([all_r_hats[param][i] for param in parameters]) for i in range(num_samples)])
    # plot_values(mean_rhat_values, median_rhat_values, rhat_threshold, crossing_point_rhat, all_rhat_datasets_flat, output_pdf=f"{outputprefix}_rhat.pdf", ylabel="R-hat")
    plot_values(None, None, rhat_threshold, crossing_point_rhat, calc_freq, all_rhat_datasets_flat, output_pdf=f"{outputprefix}_rhat.pdf", ylabel="R-hat")
