
import argparse

parser = argparse.ArgumentParser(description="Add MrBayes block to a Nexus file.")
parser.add_argument("--in_nexus", required=True, help="Path to the input Nexus file")
parser.add_argument("--out_nexus", required=True, help="Path to the desired output Nexus file")
parser.add_argument("--mcmc_length", type=int, required=False, default=1000000, help="Number of generations for MCMC")
parser.add_argument("--sample_freq", type=int, default=10000, required=False, help="Sampling frequency for MCMC")
parser.add_argument("--print_freq", type=int, default=10000, required=False, help="Print frequency for MCMC")
parser.add_argument("--diagn_freq", type=int, default=10000, required=False, help="Diagnostic frequency for MCMC")
parser.add_argument("--burnin", type=float, default=0.25, required=False, help="Burn-in fraction for summarizing trees")
parser.add_argument("--nruns", type=int, default=1, required=False, help="Number of independent runs for MCMC. MrBayes runs 2 runs by default, but for simplicity we set to 1 here.")
parser.add_argument("--sump", action="store_true", help="If set, run mrbayes sump command for summarizing parameters after MCMC")
parser.add_argument("--sumt", action="store_true", help="If set, run mrbayes sumt command for summarizing trees after MCMC")
args = parser.parse_args()

in_nexus = args.in_nexus
out_nexus = args.out_nexus
mcmc_length = args.mcmc_length
sample_freq = args.sample_freq
print_freq = args.print_freq
diagn_freq = args.diagn_freq
burnin = args.burnin
nruns = args.nruns
sump = args.sump
sumt = args.sumt


with open(in_nexus, "r") as in_f:
    with open(out_nexus, "w") as f:
        # Copy over all lines from the input Nexus file to the output Nexus file
        for line in in_f:
            f.write(line)
            
        # Add model specs for MrBayes
        f.write("\nBegin mrbayes;\n")
        
        # Generic run parameters
        f.write("\tset autoclose=yes nowarn=yes;\n")
        
        # JC69 model (nst=1; rates=equal, no among-site rate variation)
        f.write("\tlset nst=1 rates=equal;\n")
        
        # JC69 uses equal base frequencies
        f.write("\tprset statefreqpr=fixed(equal);\n")
        
        
        
        # MCMC parameters
        f.write(f"\tmcmc nruns={nruns} ngen={mcmc_length} samplefreq={sample_freq} printfreq={print_freq} diagnfreq={diagn_freq};\n")
        
        # Post MCMC summary commands
        if sump:
            f.write(f"\tsump relburnin=yes burninfrac={burnin};\n")
        if sumt:
            f.write(f"\tsumt relburnin=yes burninfrac={burnin};\n")
            
        f.write("End;\n")