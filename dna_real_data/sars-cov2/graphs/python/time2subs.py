# time2subs.py
import re

# If your global clock.rate != 1.0, set it here
CLOCK_RATE = 1.0

# Match:  ...[& meta stuff including rate=... ] : <length>
P_META_BEFORE_LEN = re.compile(
    r"\[\&(?P<meta>[^\]]*?)\]\s*:\s*(?P<t>[0-9eE.+-]+)"
)

# Fallback (rare):  ... ) : <length> [& meta ...]
P_LEN_BEFORE_META = re.compile(
    r":\s*(?P<t>[0-9eE.+-]+)\s*\[\&(?P<meta>[^\]]*?)\]"
)

def apply_rate_to_length(meta: str, t_str: str) -> str:
    # find rate=... inside meta (allow commas or start)
    m = re.search(r"(?:^|,)\s*rate\s*=\s*([0-9eE.+-]+)", meta)
    if not m:
        # no rate -> leave length unchanged
        return t_str
    rate = float(m.group(1)) * CLOCK_RATE
    t = float(t_str)
    return f"{t*rate:.10g}"

def convert_line(line: str) -> str:
    # Case 1: metadata before length
    def repl1(m):
        meta = m.group("meta")
        t_new = apply_rate_to_length(meta, m.group("t"))
        return f"[&{meta}]:{t_new}"

    out = P_META_BEFORE_LEN.sub(repl1, line)

    # Case 2 (fallback): length before metadata
    def repl2(m):
        meta = m.group("meta")
        t_new = apply_rate_to_length(meta, m.group("t"))
        return f":{t_new}[&{meta}]"

    out = P_LEN_BEFORE_META.sub(repl2, out)
    return out

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Convert BEAST time-tree to substitution-length tree")
    parser.add_argument("input", help="input .trees file")
    parser.add_argument("output", help="output .trees/.nex file")
    parser.add_argument("--skip", type=int, default=0, metavar="N",
                        help="skip first N tree states (burnin removal)")
    args = parser.parse_args()

    trees_seen = 0
    with open(args.input, "r") as fin, open(args.output, "w") as fout:
        for line in fin:
            s = line.lstrip()
            if s.startswith("tree ") or s.startswith("("):
                trees_seen += 1
                if trees_seen <= args.skip:
                    continue  # skip burnin state
                fout.write(convert_line(line))
            else:
                # pass through translate blocks, headers, etc.
                fout.write(line)

if __name__ == "__main__":
    main()
