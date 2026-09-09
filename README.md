# cross_annotate

A lightweight Bash pipeline linking genomic BED regions in non-model species to model-organism protein homologues. It expands regions, identifies overlapping genes from a GFF3 annotation, extracts gene-body DNA and searches a model proteome using DIAMOND BLASTX.

## Install

With Conda installed, replace the repository URL below with the URL hosting the script:

```bash
git clone <repository-url> cross_annotate
cd cross_annotate

conda create -n cross_annotate -c conda-forge -c bioconda \
    bedtools diamond samtools curl gzip -y
conda activate cross_annotate
```

Commands below assume the script is named `cross_annotate.sh` in the repository root. Adjust its path if needed. The script should use the activated environment; remove any hardcoded activation of another environment.

## Run

```bash
bash cross_annotate.sh peaks.bed reference.fa annotation.gff \
    zebrafish 5000 1e-10 results
```

Arguments, in order:

1. Peak BED file (0-based coordinates).
2. Reference genome FASTA (uncompressed).
3. Matching GFF3 annotation; only `gene` features are retained.
4. Model: `human`, `zebrafish`, `fruitfly`, `celegans`, `mouse` or `frog` (Xenopus tropicalis).
5. Bases added to each side of each peak.
6. E-value cutoff.
7. Optional output directory; default: `ca<model>`.

Reference FASTA, GFF3 and BED must use the same assembly and sequence identifiers.

Proteomes are downloaded from UniProt when absent and cached with DIAMOND databases in `./protein_databases/<model>/`. Internet access is required for the first download. Existing downloads are reused.

Optional settings:

```bash
THREADS=16 DATABASE_ROOT=/path/to/protein_databases \
    bash cross_annotate.sh peaks.bed reference.f(n)a annotation.gff \
    zebrafish 5000 1e-10 results
```

## Outputs

- `regions_slopped.bed`: expanded peaks.
- `regions_intersected.bed`: peak–gene associations (peak BED4 + gene BED6).
- `filter_intersected.bed`: unique intersecting gene intervals.
- `gene_names.txt`: associated gene names.
- `gen_regions.fa`: extracted gene-body sequences.
- `<model>_hits`: DIAMOND alignments in default 12-column tabular format.
- `<model>_proteins`: best-hit UniProt accessions, one per matched query; accessions may repeat.


Matches are candidate homologues, not confirmed orthologues.
