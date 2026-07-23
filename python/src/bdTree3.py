#!/usr/bin/env python3
"""
Generate a Newick-formatted tree using a birth–death process, with options to
(1) rescale to a target crown age (time units),
(2) ρ-sample from an oversampled tree,
(3) enforce a minimum edge length,
(4) scale to substitutions/site via either a fixed strict clock rate OR a target expected height,
(5) optionally output a rooted tree with no leading root edge (--no-stem),
(6) introduce relaxed-clock rate variation across branches (UCLN, UCG, or autocorrelated) with optional terminal jitter,
(7) normalize post-rate root-to-tip height to a target statistic for comparability.
"""

import argparse
import math
import random
import sys
from statistics import median

import dendropy as dp
from dendropy.simulate import treesim

# ---------- helpers ----------
def tree_height(tree):
    # max root-to-tip distance
    return max(leaf.distance_from_root() for leaf in tree.leaf_node_iter())

def rescale_to_height(tree, target_height):
    cur = tree_height(tree)
    if cur > 0:
        tree.scale_edges(target_height / cur)

def enforce_min_edge(tree, eps):
    """Apply a floor to all non-root edges."""
    if eps is None or eps <= 0:
        return
    for e in tree.postorder_edge_iter():
        # skip the root edge (head_node is the seed/root)
        if e.head_node is tree.seed_node:
            continue
        if e.length is not None and e.length < eps:
            e.length = eps

def rho_sample_extant(tree, target_n, suppress=True):
    """Randomly downsample extant tips to exactly target_n."""
    leaves = list(tree.leaf_node_iter())
    n = len(leaves)
    if target_n >= n:
        return  # nothing to do
    keep = set(random.sample(leaves, target_n))
    for lf in leaves:
        if lf not in keep:
            tree.prune_subtree(lf, suppress_unifurcations=suppress)

def ensure_rooted(tree, suppress=True):
    if not tree.is_rooted:
        tree.reroot_at_midpoint(update_splits=False, suppress_unifurcations=suppress)
        tree.is_rooted = True

def maybe_root_or_unroot(tree, rooted, suppress=True):
    if rooted:
        ensure_rooted(tree, suppress=suppress)
    else:
        tree.is_rooted = False

def remove_root_stem(tree, suppress=True):
    """Make the crown node the root and omit any root/stem edge length."""
    ensure_rooted(tree, suppress=suppress)
    # If root has a single child, lift to first bifurcation
    while len(tree.seed_node.child_nodes()) == 1:
        child = tree.seed_node.child_nodes()[0]
        tree.seed_node = child
    if suppress:
        tree.suppress_unifurcations()
    # Ensure no printed length on the root edge
    if tree.seed_node.edge is not None:
        tree.seed_node.edge.length = None

def root_to_tip_lengths(tree):
    return [lf.distance_from_root() for lf in tree.leaf_node_iter()]

def rescale_to_target_stat(tree, target_value, stat="median"):
    if target_value is None:
        return
    rtt = root_to_tip_lengths(tree)
    if not rtt:
        return
    cur = (median(rtt) if stat == "median" else sum(rtt) / len(rtt))
    if cur > 0:
        tree.scale_edges(target_value / cur)

# ---------- relaxed clock models ----------
def apply_ucln(tree, log_sd):
    """
    Uncorrelated LogNormal: independent lognormal rates per edge with mean 1.
    If X ~ LogNormal(mu, sigma), mean = exp(mu + 0.5*sigma^2).
    Set mu = -0.5*sigma^2 to make mean 1.
    """
    mu = -0.5 * (log_sd ** 2)
    for e in tree.postorder_edge_iter():
        if e.length is None or e.head_node is tree.seed_node:
            continue
        r = random.lognormvariate(mu, log_sd)
        e.length *= r

def apply_ucg(tree, k):
    """
    Uncorrelated Gamma: Gamma(k, theta=1/k) so mean=1, var=1/k.
    """
    theta = 1.0 / k
    for e in tree.postorder_edge_iter():
        if e.length is None or e.head_node is tree.seed_node:
            continue
        r = random.gammavariate(k, theta)
        e.length *= r

def apply_autocorrelated(tree, log_sd):
    """
    Autocorrelated (Thorne–Kishino): rate on child ~ LogNormal(log(parent_rate), sigma*sqrt(time)).
    Root rate fixed at 1.0.
    """
    # Do a preorder traversal so parent rates are set before children
    # DendroPy does not have explicit preorder_edge_iter; use node iteration.
    # We'll compute child rates then scale each child edge by its rate.
    # Store temporary __rate on nodes.
    for n in tree.preorder_node_iter():
        if n is tree.seed_node:
            n.__rate = 1.0
            continue
        p = n.parent_node
        parent_rate = getattr(p, "__rate", 1.0)
        t = n.edge.length if (n.edge and n.edge.length is not None) else 0.0
        t = max(t, 1e-12)
        mu = math.log(parent_rate)
        sigma = log_sd * math.sqrt(t)
        r = random.lognormvariate(mu, sigma)
        n.__rate = r
        if n.edge and n.edge.length is not None:
            n.edge.length *= r
    # cleanup
    for n in tree.postorder_node_iter():
        if hasattr(n, "__rate"):
            delattr(n, "__rate")

def add_terminal_jitter(tree, lam):
    """
    Add Exp(lam) noise to pendant (terminal) edges; mean increment = 1/lam.
    """
    if lam and lam > 0:
        for lf in tree.leaf_node_iter():
            e = lf.edge
            if e.length is not None:
                e.length += random.expovariate(lam)

# ---------- CLI ----------
brate_def = 1.0
drate_def = 0.5
ntips_def = 100

desc = "Generate a Newick tree from a birth–death process with realism knobs."
p = argparse.ArgumentParser(description=desc)
p.add_argument("-b","--brate", type=float, default=brate_def,
               help=f"birth rate λ (default {brate_def})")
p.add_argument("-d","--drate", type=float, default=drate_def,
               help=f"death rate μ (default {drate_def})")
p.add_argument("-n","--ntips", type=int, default=ntips_def,
               help=f"target number of tips after ρ-sampling (default {ntips_def})")

# realism knobs
p.add_argument("--height", type=float, default=None,
               help="target crown age in TIME units (rescale edges so root→tip max distance = HEIGHT)")
p.add_argument("--rho", type=float, default=None,
               help="ρ-sampling fraction of extant tips to keep; simulate ~ntips/ρ tips, then downsample to exactly ntips")
p.add_argument("--oversample-k", type=float, default=None,
               help="Alternative to --rho: simulate k×ntips tips then downsample to ntips")
p.add_argument("--min-edge", type=float, default=None,
               help="minimum edge length ε to enforce post hoc (applied after rate changes, before output)")

# strict clock / expected-height (mutually exclusive)
clock = p.add_mutually_exclusive_group()
clock.add_argument("--strict-clock-rate", type=float, default=None,
                   help="multiply all branch lengths by this rate to convert time→subs/site")
clock.add_argument("--expected-height", type=float, default=None,
                   help="TARGET final crown age in substitutions/site; used for post-rate normalization")

# relaxed clock options
rc = p.add_argument_group("relaxed_clock", "Branch-specific rate variation")
rcm = rc.add_mutually_exclusive_group()
rcm.add_argument("--ucln-sd", type=float, default=None,
                 help="Uncorrelated lognormal relaxed clock; log-sd (sigma) of rates (mean 1).")
rcm.add_argument("--ucg-k", type=float, default=None,
                 help="Uncorrelated gamma relaxed clock; shape k (mean 1, var=1/k).")
rcm.add_argument("--ac-sd", type=float, default=None,
                 help="Autocorrelated Brownian relaxed clock; log-sd per sqrt(time).")
rc.add_argument("--target-stat", choices=["mean","median"], default="median",
                help="When normalizing to expected height, match this root-to-tip statistic (default: median).")
rc.add_argument("--terminal-jitter", type=float, default=0.0,
                help="Add Exp(lambda) noise to pendant edges with mean=1/lambda (0 = off).")

p.add_argument("--seed", type=int, default=None, help="random seed for reproducibility")

# output/cleanup
p.add_argument("--suppress-unifurcations", dest="suppress_uni", action="store_true",
               help="suppress unifurcations after pruning/rerooting (default)")
p.add_argument("--no-suppress-unifurcations", dest="suppress_uni", action="store_false",
               help="do NOT suppress unifurcations after pruning/rerooting")
p.set_defaults(suppress_uni=True)

p.add_argument("--no-stem", action="store_true",
               help="Rooted output with no leading root/stem edge (root edge omitted).")

grp = p.add_mutually_exclusive_group()
grp.add_argument("--rooted", action="store_true", help="output rooted tree")
grp.add_argument("--unrooted", action="store_true", help="output unrooted tree (default)")
p.set_defaults(rooted=False)

args = p.parse_args()

# ---------- RNG & taxon namespace ----------
if args.seed is not None:
    random.seed(args.seed)
    try:
        dp.random.seed(args.seed)  # some versions expose this
    except Exception:
        pass

tns = dp.TaxonNamespace()

# ---------- determine oversampling size ----------
if args.rho is not None and args.oversample_k is not None:
    print("ERROR: use either --rho or --oversample-k, not both.", file=sys.stderr)
    sys.exit(2)

if args.rho is not None:
    if not (0 < args.rho <= 1.0):
        print("ERROR: --rho must be in (0,1].", file=sys.stderr)
        sys.exit(2)
    sim_ntips = max(args.ntips, int(round(args.ntips / args.rho)))
elif args.oversample_k is not None:
    if args.oversample_k < 1.0:
        print("ERROR: --oversample-k must be >= 1.", file=sys.stderr)
        sys.exit(2)
    sim_ntips = max(args.ntips, int(round(args.ntips * args.oversample_k)))
else:
    sim_ntips = args.ntips

# ---------- simulate ----------
tree = treesim.birth_death_tree(
    birth_rate=float(args.brate),
    death_rate=float(args.drate),
    num_extant_tips=int(sim_ntips),
    taxon_namespace=tns
)

# ---------- ρ-sample down to exactly ntips (if oversampled) ----------
if sim_ntips > args.ntips:
    rho_sample_extant(tree, args.ntips, suppress=args.suppress_uni)

# ---------- rescale to target crown age in TIME units (if requested) ----------
if args.height is not None:
    rescale_to_height(tree, float(args.height))

# ---------- rooting / no-stem preparation ----------
want_rooted = args.rooted or args.no_stem
if want_rooted:
    ensure_rooted(tree, suppress=args.suppress_uni)
# If removing the stem, do it BEFORE any expected-height normalization,
# so normalization targets the crown height.
if args.no_stem:
    remove_root_stem(tree, suppress=args.suppress_uni)

# ---------- validate clock choices ----------
relaxed = (args.ucln_sd is not None) or (args.ucg_k is not None) or (args.ac_sd is not None)
if relaxed and args.strict_clock_rate is not None:
    print("ERROR: Relaxed clock flags are incompatible with --strict-clock-rate.", file=sys.stderr)
    sys.exit(2)

# ---------- apply rate model & normalization ----------
# At this point, lengths are still in TIME units (after --height and --no-stem).
if relaxed:
    # Apply relaxed rates (convert time to subs/site via edge-specific multipliers)
    if args.ucln_sd is not None:
        if args.ucln_sd <= 0:
            print("ERROR: --ucln-sd must be > 0.", file=sys.stderr)
            sys.exit(2)
        apply_ucln(tree, float(args.ucln_sd))
    elif args.ucg_k is not None:
        if args.ucg_k <= 0:
            print("ERROR: --ucg-k must be > 0.", file=sys.stderr)
            sys.exit(2)
        apply_ucg(tree, float(args.ucg_k))
    elif args.ac_sd is not None:
        if args.ac_sd <= 0:
            print("ERROR: --ac-sd must be > 0.", file=sys.stderr)
            sys.exit(2)
        apply_autocorrelated(tree, float(args.ac_sd))

    # Optional pendant-edge jitter
    if args.terminal_jitter and args.terminal_jitter > 0:
        add_terminal_jitter(tree, float(args.terminal_jitter))

    # Final normalization to target expected height in subs/site (if requested)
    if args.expected_height is not None:
        rescale_to_target_stat(tree, float(args.expected_height), stat=args.target_stat)

elif args.strict_clock_rate is not None:
    # strict clock conversion time -> subs/site
    r = float(args.strict_clock_rate)
    if r != 1.0:
        tree.scale_edges(r)
elif args.expected_height is not None:
    # strict normalization to expected height (subs/site) assuming mean rate ~ 1
    rescale_to_target_stat(tree, float(args.expected_height), stat="median")

# ---------- clean up small edges (floor does not affect the root edge) ----------
enforce_min_edge(tree, args.min_edge)

# ---------- set final rooted/unrooted choice & bipartitions ----------
if not want_rooted and not args.unrooted:
    tree.is_rooted = False
# If user explicitly asked for unrooted and also --no-stem, honor no-stem (which implies rooted)
if args.unrooted and args.no_stem:
    tree.is_rooted = True

tree.update_bipartitions()

# ---------- emit ----------
sys.stdout.write(tree.as_string(schema="newick"))
