#!/usr/bin/env python3
"""
hotspot_scan.py

Detect plastome divergence hotspots by sliding-window nucleotide
diversity (pi) and number of segregating sites (S), following the
DnaSP sliding-window approach.

Steps
  1. Keep only gap-free sites (columns where every taxon has A/C/G/T;
     equivalent to DnaSP "gaps not considered" / complete deletion).
  2. Slide a window of --window gap-free sites with a step of --step
     sites. For each window compute
        pi = mean pairwise differences / number of sites (Nei & Li 1979)
        S  = number of segregating sites
  3. Flag windows where pi or S exceeds mean + --sd * SD.
  4. Merge flagged windows that share at least one site into hotspots.
  5. Report each hotspot in alignment columns and in ungapped
     coordinates of the reference taxon, with peak pi, peak S, which
     statistic flagged it, and (optionally) the overlapping or flanking
     genes from the reference GFF3.

Outputs
  <out>_windows.tsv   every window with pi, S and flags (for plotting and
                      for cross-checking a few windows against DnaSP)
  <out>_hotspots.tsv  merged hotspots with details
  <out>_regions.tsv   name<TAB>start<TAB>end in reference coordinates,
                      ready for hotspot_haplotypes.py --regions

Example
  python3 hotspot_scan.py --aln aln17.fasta --ref Dendrocalamus_strictus \
      --gff Dstrictus_curated.gff3 --out scan
"""

import argparse
import statistics
import sys
from itertools import combinations

BASES = set("ACGT")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf).upper()
    if len({len(s) for s in seqs.values()}) != 1:
        sys.exit("ERROR: sequences in the alignment differ in length")
    return seqs


def read_gff_genes(path):
    """Return sorted list of (start, end, name) for gene features."""
    genes = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            name = attrs.get("Name", attrs.get("gene", attrs.get("ID", "?")))
            genes.append((int(f[3]), int(f[4]), name))
    return sorted(genes)


def annotate(genes, s, e):
    """Describe a reference interval by overlapping or flanking genes."""
    over = [g[2] for g in genes if g[0] <= e and g[1] >= s]
    if over:
        return "within/overlapping: " + ", ".join(dict.fromkeys(over))
    left = [g for g in genes if g[1] < s]
    right = [g for g in genes if g[0] > e]
    l = max(left, key=lambda g: g[1])[2] if left else "start"
    r = min(right, key=lambda g: g[0])[2] if right else "end"
    return "IGS: %s - %s" % (l, r)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--gff", help="GFF3 of the reference (optional)")
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--step", type=int, default=25)
    ap.add_argument("--sd", type=float, default=3.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    seqs = read_fasta(a.aln)
    if a.ref not in seqs:
        sys.exit("ERROR: reference '%s' not found. Headers: %s"
                 % (a.ref, ", ".join(seqs)))
    taxa = list(seqs)
    n = len(taxa)
    L = len(seqs[taxa[0]])
    pairs = list(combinations(range(n), 2))
    rows = [seqs[t] for t in taxa]

    # ungapped reference position for every alignment column
    refpos, p = [], 0
    for ch in seqs[a.ref]:
        if ch != "-":
            p += 1
        refpos.append(p)

    # gap-free sites and per-site statistics
    sites = []          # alignment column (0-based)
    diffs = []          # number of differing pairs at the site
    seg = []            # 1 if segregating
    for c in range(L):
        col = [r[c] for r in rows]
        if any(x not in BASES for x in col):
            continue
        sites.append(c)
        diffs.append(sum(1 for i, j in pairs if col[i] != col[j]))
        seg.append(1 if len(set(col)) > 1 else 0)
    print("Taxa: %d  Alignment columns: %d  Gap-free sites: %d"
          % (n, L, len(sites)))

    npairs = len(pairs)
    windows = []
    for i in range(0, len(sites) - a.window + 1, a.step):
        j = i + a.window
        pi = sum(diffs[i:j]) / npairs / a.window
        S = sum(seg[i:j])
        windows.append([i, j - 1, pi, S])

    pis = [w[2] for w in windows]
    Ss = [w[3] for w in windows]
    tpi = statistics.mean(pis) + a.sd * statistics.pstdev(pis)
    tS = statistics.mean(Ss) + a.sd * statistics.pstdev(Ss)
    print("Windows: %d  mean pi %.5f  threshold pi %.5f  mean S %.2f  "
          "threshold S %.2f" % (len(windows), statistics.mean(pis), tpi,
                                statistics.mean(Ss), tS))

    with open(a.out + "_windows.tsv", "w") as fw:
        fw.write("win_first_site\twin_last_site\taln_start\taln_end\t"
                 "ref_start\tref_end\tref_mid\tpi\tS\tflag_pi\tflag_S\n")
        for i, j, pi, S in windows:
            cs, ce = sites[i], sites[j]
            fw.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.5f\t%d\t%d\t%d\n"
                     % (i + 1, j + 1, cs + 1, ce + 1, refpos[cs], refpos[ce],
                        (refpos[cs] + refpos[ce]) // 2, pi, S,
                        pi > tpi, S > tS))

    # merge flagged windows that share sites
    flagged = [w + [w[2] > tpi, w[3] > tS] for w in windows
               if w[2] > tpi or w[3] > tS]
    hot = []
    for w in flagged:
        if hot and w[0] <= hot[-1]["last"]:
            h = hot[-1]
            h["last"] = max(h["last"], w[1])
            h["pi"] = max(h["pi"], w[2])
            h["S"] = max(h["S"], w[3])
            h["by_pi"] |= w[4]
            h["by_S"] |= w[5]
        else:
            hot.append({"first": w[0], "last": w[1], "pi": w[2], "S": w[3],
                        "by_pi": w[4], "by_S": w[5]})
    hot.sort(key=lambda h: (-h["pi"], -h["S"]))

    genes = read_gff_genes(a.gff) if a.gff else []
    with open(a.out + "_hotspots.tsv", "w") as fh, \
         open(a.out + "_regions.tsv", "w") as fr:
        fh.write("rank\tref_start\tref_end\taln_start\taln_end\taln_cols\t"
                 "gapfree_sites\tpeak_pi\tpeak_S\tdetected_by\tannotation\n")
        fr.write("# name\tstart\tend (ungapped %s coordinates)\n" % a.ref)
        for k, h in enumerate(hot, 1):
            cs, ce = sites[h["first"]], sites[h["last"]]
            rs, re_ = refpos[cs], refpos[ce]
            by = ", ".join(x for x, f in (("pi", h["by_pi"]), ("S", h["by_S"])) if f)
            ann = annotate(genes, rs, re_) if genes else ""
            fh.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.5f\t%d\t%s\t%s\n"
                     % (k, rs, re_, cs + 1, ce + 1, ce - cs + 1,
                        h["last"] - h["first"] + 1, h["pi"], h["S"], by, ann))
            fr.write("HS%02d\t%d\t%d\n" % (k, rs, re_))
    print("Hotspots: %d. Wrote %s_windows.tsv, %s_hotspots.tsv, %s_regions.tsv"
          % (len(hot), a.out, a.out, a.out))


if __name__ == "__main__":
    main()
