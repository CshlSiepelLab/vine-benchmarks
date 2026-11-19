#!/usr/bin/env python3
"""Split (clade) posterior probabilities and supporting substitutions for two
posterior tree samples

  splits.py --reference beast.nwk --estimate vine.nwk --sites 29903 > splits.tsv

Each input has one Newick tree per line, branch lengths in substitutions/site.
Splits are unrooted bipartitions of the taxa; trivial splits (one side with
fewer than two taxa) are ignored, and the two edges below a bifurcating root,
which induce the same bipartition, count as one edge whose length is their sum.

For every split seen in either sample, writes (tab-separated, with header):
  split_id       index, ordered by decreasing reference posterior probability
  size           taxa on the smaller side
  ref_pp         fraction of reference trees containing the split
  ref_mean_len   mean length of its edge over the reference trees containing it
  ref_support    expected supporting substitutions under the reference =
                 ref_pp * ref_mean_len * sites
  est_pp, est_mean_len, est_support   the same for the estimate sample
  taxa           taxa on the smaller side, comma-separated
"""
import argparse
import re
import sys

TOKEN = re.compile(r"\(|\)|,|;|:[^,();]+|[^(),:;]+")


def parse(newick):
    """Yield (leaf_names, edges): edges is a list of (child_node_id, parent_id,
    length) over a tree stored as parent pointers; leaves carry names."""
    tokens = TOKEN.findall(newick.strip())
    parent, length, name = [-1], [0.0], [None]
    stack, cur = [], 0
    for tok in tokens:
        if tok == "(":
            parent.append(cur); length.append(0.0); name.append(None)
            stack.append(cur); cur = len(parent) - 1
        elif tok == ",":
            p = stack[-1]
            parent.append(p); length.append(0.0); name.append(None)
            cur = len(parent) - 1
        elif tok == ")":
            cur = stack.pop()
        elif tok == ";":
            break
        elif tok.startswith(":"):
            length[cur] = float(tok[1:])
        else:
            name[cur] = tok.strip().strip("'\"")
    # node 0 is the root (the outermost parenthesis); each "(" creates the
    # first child of the current node and descends into it
    return parent, length, name


def tree_splits(line, index, full):
    """Return {canonical split mask: edge length} for one tree."""
    parent, length, name = parse(line)
    n = len(parent)
    children = [[] for _ in range(n)]
    for v in range(1, n):
        children[parent[v]].append(v)
    root = 0   # the outermost "(...)" is node 0
    mask = [0] * n
    order = []
    stack = [root]
    while stack:
        v = stack.pop(); order.append(v); stack.extend(children[v])
    for v in reversed(order):
        if not children[v]:
            if name[v] not in index:
                sys.exit(f"unknown taxon {name[v]!r}")
            mask[v] = 1 << index[name[v]]
        else:
            m = 0
            for c in children[v]:
                m |= mask[c]
            mask[v] = m
    if mask[root] != full:
        sys.exit("a tree does not contain every taxon")
    ntax = full.bit_length()
    out = {}
    for v in order:
        if v == root:
            continue
        m = mask[v]
        canon = m if not (m & 1) else full ^ m      # side without taxon 0
        k = bin(canon).count("1")
        if k < 2 or k > ntax - 2:
            continue
        out[canon] = out.get(canon, 0.0) + length[v]  # root edges merge here
    return out


def tally(path, index, full):
    count, total, ntrees = {}, {}, 0
    with open(path) as fh:
        for line in fh:
            if "(" not in line:
                continue
            ntrees += 1
            for s, bl in tree_splits(line, index, full).items():
                count[s] = count.get(s, 0) + 1
                total[s] = total.get(s, 0.0) + bl
    return count, total, ntrees


def leaf_names(path):
    with open(path) as fh:
        for line in fh:
            if "(" in line:
                parent, _, name = parse(line)
                kids = set(parent[1:])
                return sorted(name[v] for v in range(1, len(parent)) if v not in kids)
    sys.exit(f"no trees in {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reference", required=True, help="reference posterior (e.g. BEAST 2)")
    ap.add_argument("--estimate", required=True, help="posterior to compare (e.g. VINE)")
    ap.add_argument("--sites", type=int, required=True, help="alignment length")
    a = ap.parse_args()

    taxa = leaf_names(a.reference)
    index = {t: i for i, t in enumerate(taxa)}
    full = (1 << len(taxa)) - 1
    rc, rt, rn = tally(a.reference, index, full)
    ec, et, en = tally(a.estimate, index, full)
    print(f"{a.reference}: {rn} trees; {a.estimate}: {en} trees; {len(taxa)} taxa", file=sys.stderr)

    rows = []
    for s in set(rc) | set(ec):
        r_pp = rc.get(s, 0) / rn
        r_len = rt[s] / rc[s] if s in rc else 0.0
        e_pp = ec.get(s, 0) / en
        e_len = et[s] / ec[s] if s in ec else 0.0
        side = s if bin(s).count("1") <= len(taxa) / 2 else full ^ s
        rows.append((r_pp, r_len, e_pp, e_len, side))
    rows.sort(key=lambda r: (-r[0], -r[2]))
    print("split_id\tsize\tref_pp\tref_mean_len\tref_support\test_pp\test_mean_len\test_support\ttaxa")
    for i, (r_pp, r_len, e_pp, e_len, side) in enumerate(rows, 1):
        members = [taxa[j] for j in range(len(taxa)) if side >> j & 1]
        print(f"{i}\t{len(members)}\t{r_pp:.6g}\t{r_len:.6g}\t{r_pp * r_len * a.sites:.6g}\t"
              f"{e_pp:.6g}\t{e_len:.6g}\t{e_pp * e_len * a.sites:.6g}\t{','.join(members)}")


if __name__ == "__main__":
    main()
