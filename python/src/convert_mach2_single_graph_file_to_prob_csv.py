

inputfile = "/local/storage/no-backup/vine-benchmarks/tissue_migration_real_data/quinn.4.cass.mach2.LL-G-0.graph"
outputfile = inputfile + ".csv"

with open(inputfile, "r") as f, open(outputfile, "w") as out:
    for line in f:
        source, target, edge_num = line.strip().split()
        edge_num = int(edge_num)
        for i in range(1, edge_num + 1):
            out.write(f"{source}_{target}_{i},1.0\n")
            