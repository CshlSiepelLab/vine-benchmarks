import numpy as np
import sys

# ---- USER SETTING ----
m = 100    # number of replicates per experiment (set this appropriately)
# -----------------------

fname = sys.argv[1]

# Read file, ignoring the dashed line
raw = []
with open(fname) as f:
    lines = f.readlines()
    
for line in lines[1:]:
        if line.strip().startswith("-----"):
            continue
        if line.strip() == "":
            continue
        raw.append(line.strip().split())

# Convert to float
data = np.array(raw, dtype=float)

# Last row = grand means (ignore final SDs)
grand_row = data[-1, :]

# All experiment rows
exp = data[:-1, :]

nrows, ncols = exp.shape
n = nrows  # number of experiments
assert ncols % 2 == 0, "Expected pairs of (mean, sd) columns."

npairs = ncols // 2

grand_means = []
global_sds = []

for k in range(npairs):
    mean_col = exp[:, 2*k]      # mu_i
    sd_col   = exp[:, 2*k + 1]  # sigma_i

    mu = grand_row[2*k]         # grand mean already provided
    mu_i = mean_col
    sigma_i = sd_col

    # sample variance of the mu_i
    S2 = np.var(mu_i, ddof=1)

    # compute global variance
    numerator = np.sum((m-1)*(sigma_i**2)) + m*(n-1)*S2
    denom = n*m - 1
    s_global = np.sqrt(numerator / denom)

    grand_means.append(mu)
    global_sds.append(s_global)

# Output final row: pairs of (grand mean, global sd)
out = []
for gm, gs in zip(grand_means, global_sds):
    out.append(f"{gm:.6f}")
    out.append(f"{gs:.6f}")

print("\t".join(out))
