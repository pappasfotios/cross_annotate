#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  cross_annotate.sh PEAKS.bed REFERENCE.f(n)a ANNOTATION.gff MODEL_SPECIES SLOP_VALUE EVALUE [OUTPUT_DIR]

MODEL:
  human | zebrafish | fruitfly | celegans | mouse | frog

Details:
  - The GFF may contain any feature types; only rows with type "gene" are used.
  - "frog" corresponds to Xenopus tropicalis.
  - Model proteomes and DIAMOND databases are cached under:
      ${DATABASE_ROOT:-./protein_databases}/MODEL/
  - THREADS controls DIAMOND threads (default: 8).
EOF
}

if [[ $# -lt 6 || $# -gt 7 ]]; then
    usage >&2
    exit 1
fi

peaks_bed=$1
reference_fasta=$2
genes_gff=$3
model_input=${4,,}
slop_value=$5
evalue_cutoff=$6
output_dir=${7:-}

for program in awk bedtools curl diamond gzip samtools sort; do
    command -v "$program" >/dev/null 2>&1 || {
        echo "Error: required program not found: $program" >&2
        exit 1
    }
done

for input_file in "$peaks_bed" "$reference_fasta" "$genes_gff"; do
    [[ -s "$input_file" ]] || {
        echo "Error: input file is missing or empty: $input_file" >&2
        exit 1
    }
done

[[ "$slop_value" =~ ^[0-9]+$ ]] || {
    echo "Error: SLOP must be a non-negative integer." >&2
    exit 1
}

[[ "$evalue_cutoff" =~ ^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]] || {
    echo "Error: EVALUE must be numeric, for example 1e-10." >&2
    exit 1
}

case "$model_input" in
    human|hs)
        model=human
        proteome_id=UP000005640
        ;;
    zebrafish|zf)
        model=zebrafish
        proteome_id=UP000000437
        ;;
    fruitfly|fruit-fly|fly|dm)
        model=fruitfly
        proteome_id=UP000000803
        ;;
    celegans|c.elegans|worm|ce)
        model=celegans
        proteome_id=UP000001940
        ;;
    mouse|mm)
        model=mouse
        proteome_id=UP000000589
        ;;
    frog|xenopus|x.tropicalis|xt)
        model=frog
        proteome_id=UP000008143
        ;;
    *)
        echo "Error: unsupported MODEL: $model_input" >&2
        usage >&2
        exit 1
        ;;
esac

if [[ -z "$output_dir" ]]; then
    output_dir="ca${model}"
fi

database_root=${DATABASE_ROOT:-./protein_databases}
database_dir="${database_root}/${model}"
proteome_fasta="${database_dir}/${model}.proteome.faa.gz"
diamond_db="${database_dir}/${model}"

mkdir -p "$output_dir" "$database_dir"

peaks_clean="${output_dir}/peaks.bed"
genes_bed="${output_dir}/genes.bed"
gene_metadata="${output_dir}/gene_metadata.tsv"
regions_slopped="${output_dir}/regions_slopped.bed"
regions_intersected="${output_dir}/regions_intersected.bed"
filter_intersected="${output_dir}/filter_intersected.bed"
gene_names="${output_dir}/gene_names.txt"
gene_fasta="${output_dir}/gen_regions.fa"
hits="${output_dir}/${model}_hits"
proteins="${output_dir}/${model}_proteins"

# Retain BED coordinates and add a stable peak identifier.
awk 'BEGIN {OFS="\t"}
     !/^#/ && NF >= 3 {
         if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 > $3) {
             print "Error: invalid BED coordinates on input line " NR > "/dev/stderr"
             exit 1
         }
         print $1, $2, $3, "peak_" ++n
     }
     END {
         if (n == 0) {
             print "Error: no BED records found" > "/dev/stderr"
             exit 1
         }
     }' "$peaks_bed" > "$peaks_clean"

# Convert only GFF gene features from 1-based closed coordinates to BED6.
# Use the GFF ID as the unique sequence identifier and retain its display name.
awk -v metadata="$gene_metadata" '
    BEGIN {
        FS=OFS="\t"
        print "gene_id", "gene_name" > metadata
    }
    function attribute(key, value, n, fields, i, prefix) {
        value=$9
        n=split(value, fields, ";")
        prefix=key "="
        for (i=1; i<=n; i++) {
            sub(/^[[:space:]]+/, "", fields[i])
            if (index(fields[i], prefix) == 1) {
                return substr(fields[i], length(prefix) + 1)
            }
        }
        return ""
    }
    !/^#/ && $3 == "gene" {
        id=attribute("ID")
        name=attribute("Name")
        if (name == "") name=attribute("gene")
        if (id == "") id="gene_" ++fallback_id
        if (name == "") name=id
        print $1, $4-1, $5, id, ".", $7
        print id, name >> metadata
        n_genes++
    }
    END {
        if (n_genes == 0) {
            print "Error: no features with GFF type gene were found" > "/dev/stderr"
            exit 1
        }
    }' "$genes_gff" > "$genes_bed"

if [[ ! -s "${reference_fasta}.fai" ]]; then
    samtools faidx "$reference_fasta"
fi

bedtools slop \
    -i "$peaks_clean" \
    -b "$slop_value" \
    -g "${reference_fasta}.fai" \
    > "$regions_slopped"

# Output columns: peak BED4 followed by gene BED6.
bedtools intersect \
    -a "$regions_slopped" \
    -b "$genes_bed" \
    -wa -wb \
    > "$regions_intersected"

# Keep each intersecting gene once.
cut -f5-10 "$regions_intersected" \
    | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n -k4,4 -u \
    > "$filter_intersected"

if [[ ! -s "$filter_intersected" ]]; then
    : > "$gene_names"
    : > "$gene_fasta"
    : > "$hits"
    : > "$proteins"
    echo "No genes intersected the slopped peaks. Empty result files were created in: $output_dir"
    exit 0
fi

# Write gene display names in same ordr as filter_intersected.bed.
awk 'BEGIN {FS=OFS="\t"}
     NR==FNR {
         if (FNR > 1) name[$1]=$2
         next
     }
     {print (($4 in name) ? name[$4] : $4)}' \
    "$gene_metadata" "$filter_intersected" \
    > "$gene_names"

# DIAMOND blastx searches both strands.
bedtools getfasta \
    -fi "$reference_fasta" \
    -bed "$filter_intersected" \
    -nameOnly \
    -fo "$gene_fasta"

if [[ ! -s "$proteome_fasta" ]]; then
    echo "Downloading $model reference proteome ($proteome_id)..."
    temporary_fasta="${proteome_fasta}.part"
    trap 'rm -f "$temporary_fasta"' EXIT

    curl --fail --location --retry 3 --get \
        --data-urlencode "compressed=true" \
        --data-urlencode "format=fasta" \
        --data-urlencode "query=(proteome:${proteome_id})" \
        --output "$temporary_fasta" \
        "https://rest.uniprot.org/uniprotkb/stream"

    gzip -t "$temporary_fasta"
    mv "$temporary_fasta" "$proteome_fasta"
    trap - EXIT
fi

if [[ ! -s "${diamond_db}.dmnd" ]]; then
    echo "Building DIAMOND database for $model..."
    diamond makedb --in "$proteome_fasta" --db "$diamond_db"
fi

diamond blastx \
    --db "$diamond_db" \
    --query "$gene_fasta" \
    --evalue "$evalue_cutoff" \
    --sensitive \
    --threads "${THREADS:-8}" \
    --out "$hits"

# Keep the lowest-E-value hit per query. Break ties using the highest bit score.
LC_ALL=C sort -t $'\t' -k1,1 -k11,11g -k12,12gr "$hits" \
    | awk -F'\t' '!seen[$1]++ {
          id=$2
          if (id ~ /^[^|]+\|[^|]+\|/) {
              split(id, part, "|")
              id=part[2]
          }
          sub(/\.[0-9]+$/, "", id)
          print id
      }' \
    > "$proteins"

echo "Annotation complete: $output_dir"
echo "Best-hit protein accessions: $proteins"
