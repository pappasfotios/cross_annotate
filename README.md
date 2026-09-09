# cross_annotate

A lightweight Bash pipeline linking genomic regions from non-model species (in BED files) to model-organism protein homologues. It expands regions (flanking genomic intervals), identifies overlapping genes from a GFF3 annotation, extracts gene-body DNA and searches a model proteome using DIAMOND BLASTX.

## Install

```bash
git clone https://github.com/pappasfotios/cross_annotate.git cross_annotate
cd cross_annotate

conda create -n cross_annotate -c conda-forge -c bioconda bedtools diamond samtools curl gzip -y
conda activate cross_annotate
```

The script should use the activated environment.

## Run

```bash
bash cross_annotate.sh peaks.bed reference.fa annotation.gff zebrafish 5000 1e-10 results
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

Proteomes are downloaded from UniProt when absent and cached with DIAMOND databases in `./protein_databases/<model>/`. Existing downloads are reused.

## Outputs

- `regions_slopped.bed`: expanded peaks.
- `regions_intersected.bed`: peak–gene associations (peak BED4 + gene BED6).
- `filter_intersected.bed`: unique intersecting gene intervals.
- `gene_names.txt`: associated gene names.
- `gen_regions.fa`: extracted gene-body sequences.
- `<model>_hits`: DIAMOND alignments in default 12-column tabular format.
- `<model>_proteins`: best-hit UniProt accessions, one per matched query; accessions may repeat.


Matches are candidate homologues, not confirmed orthologues.

## Inspiration

Earlier versions of this workflow were used in the following studies:

- Pappas, F., Kurta, K., Vanhala, T., Jeuthe, H., Hagen, Ø., Beirão, J., & Palaiokostas, C. (2023). Whole-genome re-sequencing provides key genomic insights in farmed Arctic charr (Salvelinus alpinus) populations of anadromous and landlocked origin from Scandinavia. Evolutionary Applications, 16, 797–813. https://doi.org/10.1111/eva.13537
  
- Pappas, F., Johnsson, M., Andersson, G. et al. Sperm DNA methylation landscape and its links to male fertility in a non-model teleost using EM-seq. Heredity 134, 293–305 (2025). https://doi.org/10.1038/s41437-025-00756-y
