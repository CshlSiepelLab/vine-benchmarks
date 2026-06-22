#!/usr/bin/env python3

"""
Convert posterior trees from a Newick text file into a NEXUS .trees file.

The input can contain non-tree log lines. Any line that does not look like a
complete Newick tree is ignored.
"""

import argparse
from pathlib import Path


def looks_like_newick(line: str) -> bool:
    """Return True when a line appears to be a complete Newick tree."""
    stripped = line.strip()
    if not stripped:
        return False
    if not stripped.startswith("("):
        return False
    if not stripped.endswith(";"):
        return False
    if "(" not in stripped or ")" not in stripped:
        return False

    # Reject lines with unbalanced parentheses (common in truncated trees).
    balance = 0
    for char in stripped:
        if char == "(":
            balance += 1
        elif char == ")":
            balance -= 1
            if balance < 0:
                return False
    return balance == 0


def extract_newick_trees(input_path: Path) -> list[str]:
    """Extract Newick trees from input file, one tree per matching line."""
    trees: list[str] = []
    with input_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            if looks_like_newick(raw_line):
                trees.append(raw_line.strip())
    return trees


def write_trees_file(
    output_path: Path,
    trees: list[str],
    state_start: int,
    state_step: int,
) -> None:
    """Write trees into a NEXUS .trees file."""
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("#NEXUS\n\n")
        handle.write("Begin trees;\n")
        state = state_start
        for tree in trees:
            handle.write(f"tree STATE_{state} = [&R] {tree}\n")
            state += state_step
        handle.write("End;\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert posterior trees from a .nwk file into a .trees file."
        )
    )
    parser.add_argument("input_nwk", help="Input Newick file path")
    parser.add_argument(
        "-o",
        "--output",
        help="Output .trees file path (default: input basename + .trees)",
    )
    parser.add_argument(
        "--state-start",
        type=int,
        default=0,
        help="Starting STATE value (default: 0)",
    )
    parser.add_argument(
        "--state-step",
        type=int,
        default=1,
        help="Increment between STATE values (default: 1)",
    )
    args = parser.parse_args()

    input_path = Path(args.input_nwk)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    output_path = (
        Path(args.output)
        if args.output
        else input_path.with_suffix(".trees")
    )

    trees = extract_newick_trees(input_path)
    if not trees:
        raise ValueError(
            f"No Newick trees found in input file: {input_path}"
        )

    write_trees_file(
        output_path=output_path,
        trees=trees,
        state_start=args.state_start,
        state_step=args.state_step,
    )
    print(
        f"Wrote {len(trees)} trees to {output_path}"
    )


if __name__ == "__main__":
    main()
