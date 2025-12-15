import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")   # Use non-interactive backend
import matplotlib.backends.backend_pdf  # For pyinstaller inclusion
import matplotlib as mpl
import matplotlib.pyplot as plt


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
    mcse_values = []
    
    for i in range(len(trace)):
        start_idx = int(burnin * (i + 1))
        subset = trace[start_idx : i + 1]  # Apply burn-in and keep subsetting for more samples each iteration
        n = len(subset) # Current samples subset size
        
        if n < 2:
            ess_values.append(0.0)
            mcse_values.append(np.nan)
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
        
        # Running mcse add on
        sd = np.std(subset, ddof=1) # ddof=1 for sample standard deviation
        mcse = sd / np.sqrt(ess)
        mcse_values.append(mcse)
        
    return ess_values, mcse_values


def plot_values(mean_values, median_values, threshold, crossing_point, raw_dataset_values, output_pdf, ylabel="Effective sample size (ESS)"):
    """
    Plot ESS or MCSE values over samples and save to PDF.
    """
    plt.figure(figsize=(10, 6))
    fs = 18
    # Plot mean ESS values
    plt.plot(range(1, len(mean_values) + 1), mean_values, label="Mean", color="black", linestyle="-", linewidth=4)
    # Plot median ESS values
    plt.plot(range(1, len(median_values) + 1), median_values, label="Median", color="grey", linestyle="-", linewidth=4)
    # Plot ESS threshold
    plt.axhline(y=threshold, color="red", linestyle="-", label=f"Threshold ({threshold})", linewidth=2)
    # Plot dataset-specific ESS values
    colors = mpl.colormaps["tab20"].resampled(20).colors.tolist() + \
             mpl.colormaps["tab20b"].resampled(20).colors.tolist() + \
             mpl.colormaps["tab20c"].resampled(20).colors.tolist()
    for i, (dataset, values) in enumerate(raw_dataset_values.items()):
        plt.plot(range(1, len(values) + 1), values, linestyle="--", label=f"{dataset}", color=colors[i], linewidth=1)
    # Configure plot appearance
    plt.xlabel("Number of samples", fontsize=fs)
    plt.ylabel(ylabel, fontsize=fs)
    plt.xticks(fontsize=fs - 4)
    plt.yticks(fontsize=fs - 4)
    plt.title(f"Mean crosses threshold at sample {crossing_point}", fontsize=fs - 4)
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
parser.add_argument("--outputfile", help="Output file for the scale factor for the chain when the threshold is crossed")
parser.add_argument("--ess_threshold", type=float, help="ESS threshold for cutoff (default: 200)", default=200.0)
parser.add_argument("--burnin", type=float, help="Burn-in proportion (default: 0.1)", default=0.1)
args = parser.parse_args()

logfiles = args.logfiles
parameters = args.parameters
outputfile = args.outputfile
ess_threshold = args.ess_threshold
burnin = args.burnin

avg_ess_datasets = {}
avg_mcse_datasets = {}
for logfile in logfiles:
    logfile_basename = os.path.basename(logfile)
    all_ess_values = {}
    all_mcse_values = {}
    for param in parameters:
        trace = read_log_file(logfile, param)
        all_ess_values[param], all_mcse_values[param] = calculate_running_ess(trace, burnin=burnin)
        
    # Compute average ESS across parameters for each sample index
    num_samples = len(all_ess_values[parameters[0]])
    # Vectorize across parameters to avoid Python loops
    ess_matrix = np.array(
        [all_ess_values[param] for param in parameters],
        dtype=float,
    )
    mcse_matrix = np.array(
        [all_mcse_values[param] for param in parameters],
        dtype=float,
    )
    ess_values = ess_matrix.mean(axis=0).tolist()
    mcse_values = mcse_matrix.mean(axis=0).tolist()
    
    avg_ess_datasets[logfile_basename] = ess_values
    avg_mcse_datasets[logfile_basename] = mcse_values

# Get ess_values for mean and median across datasets
num_samples = len(next(iter(avg_ess_datasets.values())))
# Vectorize across datasets: rows=datasets, cols=samples
ess_matrix = np.array(
    list(avg_ess_datasets.values()),
    dtype=float,
)
mcse_matrix = np.array(
    list(avg_mcse_datasets.values()),
    dtype=float,
)
mean_ess_values = ess_matrix.mean(axis=0).tolist()
median_ess_values = np.median(ess_matrix, axis=0).tolist()
mean_mcse_values = mcse_matrix.mean(axis=0).tolist()
median_mcse_values = np.median(mcse_matrix, axis=0).tolist()

# Determine crossing point for mean ESS
crossing_point = next((i for i, ess in enumerate(mean_ess_values) if ess >= ess_threshold), len(ess_values) - 1) + 1
scale_factor = crossing_point / num_samples
with open(outputfile, "w") as f:
    f.write(f"{scale_factor}\n")

# Plot ESS values
plot_values(mean_ess_values, median_ess_values, ess_threshold, crossing_point, avg_ess_datasets, output_pdf=outputfile.replace(".txt", ".pdf"))

# Plot MCSE values (ignoring threshold and crossing point for now)
plot_values(mean_mcse_values, median_mcse_values, 0, 0, avg_mcse_datasets, output_pdf=outputfile.replace(".txt", "_mcse.pdf"), ylabel="Monte Carlo Standard Error (MCSE)")

