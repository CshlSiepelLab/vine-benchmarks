import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["DejaVu Sans"]
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


tsv_input = sys.argv[1]
pdf_output = sys.argv[2]

# Read in thresholded precison-recall data
df = pd.read_csv(tsv_input, sep="\t", index_col=False)

# Drop the f1 score column
df = df.drop(columns=["f1"])

# Get the mean across sim for each method and threshold
df_mean = df.groupby(["method", "threshold"]).mean().reset_index()
df_mean = df_mean.drop(columns=["sim"])

# Use formatted method names
df_mean["method"] = df_mean["method"].replace({"METIENT": "Metient"})

# Order method in order of release date
df_mean["method"] = pd.Categorical(
    df_mean["method"],
    categories=["Metient", "MACH2", "BEAM", "VINE"],
    ordered=True,
)

# Plot the precision-recall curve
plt.figure(figsize=(5, 4))
sns.lineplot(data=df_mean, x="recall", y="precision", hue="method", marker="", linewidth=2)

fs=18
plt.xlabel("Recall", fontsize=fs)
plt.ylabel("Precision", fontsize=fs)
plt.title("", fontsize=fs)
plt.xticks(fontsize=fs-4)
plt.yticks(fontsize=fs-4)

plt.legend(
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    fontsize=fs-4
)

plt.tight_layout()
plt.savefig(pdf_output)
plt.close()
