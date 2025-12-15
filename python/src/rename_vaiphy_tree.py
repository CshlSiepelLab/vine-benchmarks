#!/usr/bin/env python3
"""
Rename VaiPhy trees from internal numeric IDs to original T-names.

Inputs:
- --log: VaiPhy run log file that contains lines like:
    "Renaming: taxon0 T10"
- --in:  Input Newick file with numeric leaf IDs (from VaiPhy)
- --out: Output Newick file with original T-names

The script:
- Parses "Renaming:" lines to build a mapping.
- Rewrites any leaf token appearing before a ':' that is:
  * a zero-based index (e.g., 0)
  * a one-based index (e.g., 10)
  * a "taxonN" token
  * already a T-name (left unchanged)
- Avoids touching internal node labels that occur after ')'.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


RENAME_RE = re.compile(
    r"^Renaming:\s+taxon(\d+)\s+T(\d+)\s*$"
)

# Match a leaf name (not an internal node label), i.e. a token that:
# - starts at BOF or right after '(' or ','
# - is followed by ':' (branch length)
# - does not include parentheses, colons, or commas
LEAF_NAME_RE = re.compile(
    r"(?P<prefix>(?:^|[(,])\s*)(?P<name>[^():,]+)(?=\s*:)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rename VaiPhy tree leaves to original T-names."
    )
    parser.add_argument(
        "--log", required=True, help="Path to VaiPhy log file."
    )
    parser.add_argument(
        "--in", dest="inp", required=True,
        help="Input Newick file from VaiPhy."
    )
    parser.add_argument(
        "--out", required=True, help="Output Newick file."
    )
    return parser.parse_args()


def build_mapping(log_path: Path) -> dict[str, str]:
    """
    Build a map from various forms to T-names:
    - "n" (0-based) -> "Tt"
    - "n1" (1-based) -> "Tt"
    - "taxon{n}" -> "Tt"
    - "Tt" -> "Tt"
    """
    mapping: dict[str, str] = {}
    with log_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            m = RENAME_RE.match(line.strip())
            if not m:
                continue
            n = int(m.group(1))
            t = int(m.group(2))
            tname = f"T{t}"
            # 0-based, 1-based, taxonN, and identity
            mapping[str(n)] = tname
            mapping[str(n + 1)] = tname
            mapping[f"taxon{n}"] = tname
            mapping[tname] = tname
    return mapping


def rename_newick(nwk_text: str, mapping: dict[str, str]) -> str:
    """
    Replace leaf names using mapping. Leaves are tokens matched by
    LEAF_NAME_RE. Internal node labels (after ')') are not matched.
    """
    def repl(match: re.Match[str]) -> str:
        prefix = match.group("prefix")
        name = match.group("name").strip()
        new_name = mapping.get(name, name)
        return f"{prefix}{new_name}"

    # Apply substitution across the entire string
    renamed = LEAF_NAME_RE.sub(repl, nwk_text)
    return renamed


def main() -> None:
    args = parse_args()
    log_path = Path(args.log)
    in_path = Path(args.inp)
    out_path = Path(args.out)

    mapping = build_mapping(log_path)
    if not mapping:
        raise SystemExit(
            "No 'Renaming:' lines found in log; cannot build mapping."
        )

    nwk_text = in_path.read_text(encoding="utf-8", errors="ignore")
    # Keep first non-empty line as tree if file contains multiple lines
    if "\n" in nwk_text:
        for line in nwk_text.splitlines():
            line = line.strip()
            if line:
                nwk_text = line
                break

    renamed = rename_newick(nwk_text, mapping)
    out_path.write_text(renamed + ("\n" if not renamed.endswith("\n") else ""),
                        encoding="utf-8")


if __name__ == "__main__":
    main()


