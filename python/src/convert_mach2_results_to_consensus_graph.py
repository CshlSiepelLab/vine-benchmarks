
import sys

combined_graphs_file = sys.argv[1]
outfile = sys.argv[2]

migration_counts = {}
results_checked = set()
with open(combined_graphs_file, "r") as f:
    for line in f.readlines():
        # skip the header line
        if "result_num" in line:
            continue
        result_num, source, target, count = line.strip().split(",")
        results_checked.add(result_num)
        for i in range(1, int(count) + 1):
            migration = f"{source}_{target}_{i}"
            if migration in migration_counts:
                migration_counts[migration] += 1
            else:
                migration_counts[migration] = 1

total_num_results = len(results_checked)

# normalize the counts to probabilities by dividing by the total number of results
for migration in migration_counts:
    migration_counts[migration] /= total_num_results

# write the consensus graph to a file
with open(outfile, "w") as f:
    for migration in migration_counts:
        f.write(f"{migration},{migration_counts[migration]}\n")
