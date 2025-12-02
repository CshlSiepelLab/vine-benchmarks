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


def calculate_one_chain_running_trace_stats(trace, burnin=0.1):
    """
    Compute running stats of one trace.
    """
    ess_values = []
    mcse_values = []
    
    for i in range(len(trace)):
        start_idx = int(burnin * (i + 1))
        subset = trace[start_idx : i + 1]  # Apply burn-in and keep subsetting for more samples each iteration
        ess, mcse = calculate_one_chain_ess_and_mcse(subset)
        ess_values.append(ess)
        mcse_values.append(mcse)
        
    return ess_values, mcse_values


def rank_normalize_chains(chains):
    """
    Apply Vehtari et al. (2020) + Blom (1958) rank-normalization.
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
    z = norm.ppf((ranks - 3.0/8.0) / (S - 1.0/4.0))
    
    # Restore original chains
    out = []
    idx = 0
    for c in chains:
        n = len(c)
        out.append(z[idx:idx+n])
        idx += n
        
    return out


def calculate_rhat(chains, split=True, rank_normalize=False):
    """
    Calculate split R-hat (optionally rank normalized) for multiple chains.
    Calculations taken from Vehtari et. al. 2020 paper at https://sites.stat.columbia.edu/gelman/research/published/Vehtari_etal_2020_rhat_ess.pdf
    Assumes burnin has already been applied to the input.
    """
    # Optionally split chains in half
    if split:
        split_chains = []
        for chain in chains:
            m = len(chain) // 2
            if m < 2:
                print("Not enough samples to split chains for R-hat calculation.")
                return np.nan
            split_chains.append(chain[:m])
            split_chains.append(chain[m:])
    else:
        split_chains = chains
    
    # Number of chains
    M = len(split_chains)
    
    # Ensure a common length for all chains
    N = min(chain.shape[0] for chain in split_chains)
    split_chains = [c[:N] for c in chains]
    
    # Optionally rank-normalize chains
    if rank_normalize:
        chains = rank_normalize_chains(split_chains)
    
    # Calculate between chains variance
    chain_means = np.array([np.mean(chain) for chain in split_chains])
    overall_mean = np.mean(chain_means)
    B = ((N) / (M - 1)) * np.sum((chain_means - overall_mean) ** 2)
    
    # Calculate within chain variance
    W = np.mean([np.var(chain, ddof=1) for chain in split_chains])
    
    var_hat = (((N-1)/N)*W) + (B/N)
    r_hat = np.sqrt(var_hat / W)
    return r_hat


# def plot_values(mean_values, median_values, threshold, crossing_point, raw_dataset_values, output_pdf, ylabel="Effective sample size (ESS)"):
#     """
#     Plot ESS or MCSE values over samples and save to PDF.
#     """
#     plt.figure(figsize=(10, 6))
#     fs = 18
#     # Plot mean ESS values
#     plt.plot(range(1, len(mean_values) + 1), mean_values, label="Mean", color="black", linestyle="-", linewidth=4)
#     # Plot median ESS values
#     plt.plot(range(1, len(median_values) + 1), median_values, label="Median", color="grey", linestyle="-", linewidth=4)
#     # Plot ESS threshold
#     plt.axhline(y=threshold, color="red", linestyle="-", label=f"Threshold ({threshold})", linewidth=2)
#     # Plot dataset-specific ESS values
#     colors = mpl.colormaps["tab20"].resampled(20).colors.tolist() + \
#              mpl.colormaps["tab20b"].resampled(20).colors.tolist() + \
#              mpl.colormaps["tab20c"].resampled(20).colors.tolist()
#     for i, (dataset, values) in enumerate(raw_dataset_values.items()):
#         plt.plot(range(1, len(values) + 1), values, linestyle="--", label=f"{dataset}", color=colors[i], linewidth=1)
#     # Configure plot appearance
#     plt.xlabel("Number of samples", fontsize=fs)
#     plt.ylabel(ylabel, fontsize=fs)
#     plt.xticks(fontsize=fs - 4)
#     plt.yticks(fontsize=fs - 4)
#     plt.title(f"Mean crosses threshold at sample {crossing_point}", fontsize=fs - 4)
#     plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=fs - 4)
#     plt.grid()
#     plt.tight_layout()
#     plt.savefig(output_pdf)
#     plt.close()


def parse_csv_string(s):
    if "," not in s:
        return [s.strip()]
    return [item.strip() for item in s.split(",")]

# parser = argparse.ArgumentParser(description="Compute running ESS for a BEAST2 log file (pure Python).")
# parser.add_argument("--logfiles", type=parse_csv_string, help="Input TSV style MCMC .log files, which must all match in chain length.")
# parser.add_argument("--parameters", type=parse_csv_string, help="Parameter(s) (column name(s)) to compute ESS for. Format as csv string. If multiple parameters are provided, the threshold will be based on the average ESS across them.",)
# parser.add_argument("--outputfile", help="Output file for the scale factor for the chain when the threshold is crossed")
# parser.add_argument("--ess_threshold", type=float, help="ESS threshold for cutoff (default: 200)", default=200.0)
# parser.add_argument("--burnin", type=float, help="Burn-in proportion (default: 0.1)", default=0.1)
# args = parser.parse_args()

# logfiles = args.logfiles
# parameters = args.parameters
# outputfile = args.outputfile
# ess_threshold = args.ess_threshold
# burnin = args.burnin

# For testing purposes, hardcode arguments here
logfiles = [
    "/local/storage/no-backup/vine-benchmarks/dna2/testdata.10/tree.1.beast.log",
    "/local/storage/no-backup/vine-benchmarks/dna2/testdata.10/tree.2.beast.log"
]
parameters = ["posterior", "likelihood", "Tree.height", "Tree.Length"]
outputfile="test_mcmc_convergence/test_raw.txt"
ess_threshold = 200.0
burnin = 0.1


all_ess_datasets = {}
all_mcse_datasets = {}
mean_ess_datasets = {}
mean_mcse_datasets = {}

for logfile in logfiles:
    id = os.path.basename(logfile)
    all_ess_datasets[id] = {}
    all_mcse_datasets[id] = {}
    for param in parameters:
        trace = read_log_file(logfile, param)
        all_ess_datasets[id][param], all_mcse_datasets[id][param] = calculate_one_chain_running_trace_stats(trace, burnin=burnin)
        
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
with open(outputfile, "w") as f:
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
with open(outputfile.replace(".txt", "_summary_across_logfiles.txt"), "w") as f:
    f.write("sample_index\tmean_ess\tmedian_ess\tmean_mcse\tmedian_mcse\n")
    num_samples = len(mean_ess_values)
    for i in range(num_samples):
        f.write(f"{i + 1}\t{mean_ess_values[i]}\t{median_ess_values[i]}\t{mean_mcse_values[i]}\t{median_mcse_values[i]}\n")


#TODO: Update this old code to reinstate plotting and scale factor output
# # Determine crossing point for mean ESS
# crossing_point = next((i for i, ess in enumerate(mean_ess_values) if ess >= ess_threshold), len(mean_ess_values) - 1) + 1
# scale_factor = crossing_point / num_samples
# with open(outputfile, "w") as f:
#     f.write(f"{scale_factor}\n")

# # Plot ESS values
# plot_values(mean_ess_values, median_ess_values, ess_threshold, crossing_point, avg_ess_datasets, output_pdf=outputfile.replace(".txt", ".pdf"))

# # Plot MCSE values (ignoring threshold and crossing point for now)
# plot_values(mean_mcse_values, median_mcse_values, 0, 0, avg_mcse_datasets, output_pdf=outputfile.replace(".txt", "_mcse.pdf"), ylabel="Monte Carlo Standard Error (MCSE)")


# Calculate between chain stats
if len(logfiles) > 1:
    all_r_hat = {}
    for param in parameters:
        param_traces = []
        for logfile in logfiles:
            trace = read_log_file(logfile, param)
            param_traces.append(trace)
            all_r_hat[param] = calculate_rhat(param_traces)
        
