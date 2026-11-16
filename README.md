# Thai bamboo plastomes

Scripts for the divergence-hotspot analyses of 17 complete plastomes
of Thai bamboos, described in:

> [Authors]. [Year]. [Title]. *Plant Diversity* [volume, pages]. https://doi.org/[DOI]

Python 3 (>= 3.8), standard library only. Run a script without arguments
(or with `-h`) to see its options.

| Script | Purpose | Paper |
|---|---|---|
| `hotspot_scan.py` | Nucleotide diversity (pi, Nei and Li 1979) and number of segregating sites (S) in sliding windows of 100 gap-free sites (step 25); windows above mean + 3 SD are merged into hotspots | Section 2.4, Table 3, Fig. 4 |
| `hotspot_haplotypes.py` | Haplotypes at each hotspot, with gaps as characters and from substitutions only; forward selection of the smallest set of hotspots that gives every accession a unique haplotype | Section 2.4, Table 3, Fig. 5 |

Input: a FASTA alignment of the 17 plastomes (MAFFT v7).

## License

MIT (see `LICENSE`).
