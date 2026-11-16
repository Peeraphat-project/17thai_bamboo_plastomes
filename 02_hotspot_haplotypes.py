#!/usr/bin/env python3
"""
hotspot_haplotypes.py

Count haplotypes in plastome divergence hotspots and test how well the
hotspots discriminate the sampled taxa.

For each hotspot (given in reference coordinates, e.g. D. strictus) the
script extracts the matching alignment columns and reports:
  - variable sites and parsimony-informative sites (gap-free columns only)
  - haplotypes with gaps treated as characters
  - haplotypes from substitutions only (columns with a gap or ambiguous
    base in any taxon are removed, i.e. complete deletion as in DnaSP)
  - which taxa share each haplotype
It then finds a small set of hotspots that gives every taxon a unique
haplotype (greedy forward selection), and lists taxon pairs that no
hotspot can separate.

Input:
  --aln      aligned FASTA (all taxa, same length)
  --ref      name of the reference taxon in the alignment (exact header)
  --regions  tab-separated file: name<TAB>start<TAB>end  (1-based,
             inclusive, ungapped reference coordinates)
  --out      output prefix

Outputs:
  <out>_summary.tsv    one row per hotspot
  <out>_groups.tsv     taxa sharing a haplotype, per hotspot and mode
  <out>_report.txt     greedy combinations and unresolved pairs

Example:
  python3 hotspot_haplotypes.py --aln aln17.fasta \
      --ref "Dendrocalamus_strictus" --regions hotspots.tsv --out hs
"""

import argparse
import itertools
import sys

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
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        sys.exit("ERROR: sequences in the alignment differ in length")
    return seqs


def ref_to_columns(ref_seq):
    """Map 1-based ungapped reference position -> 0-based alignment column."""
    mapping, pos = {}, 0
    for col, ch in enumerate(ref_seq):
        if ch != "-":
            pos += 1
            mapping[pos] = col
    return mapping


def read_regions(path):
    regions = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            name, start, end = line.rstrip("\n").split("\t")[:3]
            regions.append((name, int(start), int(end)))
    return regions


def column_stats(cols):
    """Variable and parsimony-informative sites among gap-free columns."""
    var = pi = 0
    for col in cols:
        if any(c not in BASES for c in col):
            continue
        counts = {}
        for c in col:
            counts[c] = counts.get(c, 0) + 1
        if len(counts) > 1:
            var += 1
            if sum(1 for n in counts.values() if n >= 2) >= 2:
                pi += 1
    return var, pi


def haplotype_groups(taxa, seqs, start_col, end_col, subs_only):
    """Return list of taxon groups sharing identical haplotypes."""
    cols = range(start_col, end_col + 1)
    if subs_only:
        keep = [c for c in cols
                if all(seqs[t][c] in BASES for t in taxa)]
    else:
        # drop columns that are gaps in every taxon (empty columns)
        keep = [c for c in cols
                if any(seqs[t][c] != "-" for t in taxa)]
    hap = {}
    for t in taxa:
        key = "".join(seqs[t][c] for c in keep)
        hap.setdefault(key, []).append(t)
    return list(hap.values())


def partition_labels(groups):
    lab = {}
    for i, g in enumerate(groups):
        for t in g:
            lab[t] = i
    return lab


def greedy_min_set(taxa, partitions):
    """Greedy forward selection of hotspots until all taxa are unique.

    partitions: dict region -> {taxon: haplotype label}
    """
    chosen = []
    current = {t: () for t in taxa}

    def n_unique(labels):
        seen = {}
        for v in labels.values():
            seen[v] = seen.get(v, 0) + 1
        return sum(1 for t in taxa if seen[labels[t]] == 1)

    while True:
        best, best_score, best_labels = None, n_unique(current), None
        for reg, lab in partitions.items():
            if reg in chosen:
                continue
            trial = {t: current[t] + (lab[t],) for t in taxa}
            score = n_unique(trial)
            if score > best_score:
                best, best_score, best_labels = reg, score, trial
        if best is None:
            break
        chosen.append(best)
        current = best_labels
        if best_score == len(taxa):
            break
    return chosen, n_unique(current)


def unresolved_pairs(taxa, partitions):
    pairs = []
    for a, b in itertools.combinations(taxa, 2):
        if all(lab[a] == lab[b] for lab in partitions.values()):
            pairs.append((a, b))
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = read_fasta(args.aln)
    if args.ref not in seqs:
        sys.exit("ERROR: reference '%s' not found. Headers: %s"
                 % (args.ref, ", ".join(seqs)))
    taxa = list(seqs)
    pos2col = ref_to_columns(seqs[args.ref])
    regions = read_regions(args.regions)
    print("Taxa: %d  Alignment length: %d  Hotspots: %d"
          % (len(taxa), len(seqs[taxa[0]]), len(regions)))

    part = {"gaps": {}, "subs": {}}
    with open(args.out + "_summary.tsv", "w") as fs, \
         open(args.out + "_groups.tsv", "w") as fg:
        fs.write("region\tref_start\tref_end\taln_start\taln_end\taln_cols\t"
                 "variable\tparsimony_inf\thap_gaps\tunique_gaps\t"
                 "hap_subs\tunique_subs\n")
        fg.write("region\tmode\thaplotype\tn_taxa\ttaxa\n")
        for name, s, e in regions:
            if s not in pos2col or e not in pos2col:
                print("WARNING: %s outside reference range, skipped" % name)
                continue
            cs, ce = pos2col[s], pos2col[e]
            cols = [[seqs[t][c] for t in taxa] for c in range(cs, ce + 1)]
            var, pi = column_stats(cols)
            row = [name, s, e, cs + 1, ce + 1, ce - cs + 1, var, pi]
            for mode, subs in (("gaps", False), ("subs", True)):
                groups = haplotype_groups(taxa, seqs, cs, ce, subs)
                groups.sort(key=len, reverse=True)
                part[mode][name] = partition_labels(groups)
                row += [len(groups), sum(1 for g in groups if len(g) == 1)]
                for i, g in enumerate(groups, 1):
                    fg.write("%s\t%s\tH%d\t%d\t%s\n"
                             % (name, mode, i, len(g), ", ".join(g)))
            fs.write("\t".join(map(str, row)) + "\n")

    with open(args.out + "_report.txt", "w") as fr:
        for mode, label in (("gaps", "gaps as characters"),
                            ("subs", "substitutions only")):
            chosen, n = greedy_min_set(taxa, part[mode])
            pairs = unresolved_pairs(taxa, part[mode])
            fr.write("== %s ==\n" % label)
            fr.write("Greedy combination (%d hotspots): %s\n"
                     % (len(chosen), " + ".join(chosen)))
            fr.write("Taxa with unique haplotype: %d of %d\n" % (n, len(taxa)))
            fr.write("Pairs not separated by any hotspot: %d\n" % len(pairs))
            for a, b in pairs:
                fr.write("  %s  |  %s\n" % (a, b))
            fr.write("\n")
    print("Done. Wrote %s_summary.tsv, %s_groups.tsv, %s_report.txt"
          % (args.out, args.out, args.out))


if __name__ == "__main__":
    main()
