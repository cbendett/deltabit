#!/bin/bash


# Check arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <GFF_FILE> <WINDOW_SIZE>"
    exit 1
fi

GFF_FILE="$1"
WINDOW_SIZE="$2"

# ---- Step 1: Extract gene coordinates from GFF ----
echo "Extracting gene coordinates from $GFF_FILE..."

awk -F '\t' '
    $3 == "gene" {
        contig = $1
        start = $4
        end = $5
        alias = "-"
        n = split($9, attrs, ";")
        for (i = 1; i <= n; i++) {
            if (attrs[i] ~ /Alias=/) {
                split(attrs[i], kv, "=")
                alias = kv[2]
                break
            }
        }
        print contig "\t" start "\t" end "\t" alias
    }
' "$GFF_FILE" | sort -k1,1 -k2,2n > sorted_genes.tsv

# ---- Step 2: Cluster genes into windows ----
echo "Clustering genes with window size $WINDOW_SIZE..."

awk -v window="$WINDOW_SIZE" '
BEGIN {
    OFS="\t"
    print "window_id", "scaffold", "start", "end", "gene_count", "gene_names"
    cluster_size = 0
    window_num = 1
}
{
    scaffold = $1
    start = $2
    end = $3
    gene_id = $4
    if (scaffold == last_scaffold && start - last_end <= window) {
        cluster_end = (end > cluster_end ? end : cluster_end)
        gene_ids[cluster_size++] = gene_id
    } else {
        if (cluster_size >= 2) {
            print "window_" window_num++, last_scaffold, cluster_start, cluster_end, cluster_size, join(gene_ids, cluster_size)
        }
        delete gene_ids
        cluster_size = 1
        cluster_start = start
        cluster_end = end
        gene_ids[0] = gene_id
    }
    last_scaffold = scaffold
    last_end = end
}
END {
    if (cluster_size >= 2) {
        print "window_" window_num++, last_scaffold, cluster_start, cluster_end, cluster_size, join(gene_ids, cluster_size)
    }
}
function join(arr, len,   out, i) {
    out = arr[0]
    for (i = 1; i < len; i++) {
        out = out "," arr[i]
    }
    return out
}
' sorted_genes.tsv > gene_windows.tsv

# ---- Step 3: Extract first and last gene of each window ----
echo "Extracting first and last genes of each window..."

awk -F'\t' 'NR > 1 {
    n = split($6, parts, ",");
    print $1 "\t" parts[1] "\t" parts[n];
}' gene_windows.tsv > processed_genes.txt

# ---- Step 4: Process each row and retrieve genes ----
echo "Processing each row..."

declare -A ome_counts

tail -n +2 processed_genes.txt | while IFS=$'\t' read -r col1 col2 col3 rest; do
  [[ -z "$col1" || -z "$col2" || -z "$col3" ]] && continue

  first_gene=$(acc2locus -a "$col2" -p 9 | head -n 1)
  last_gene=$(acc2locus -a "$col3" -p 9 | tail -n 2)

  ome="${col2%%_*}"
  count=${ome_counts[$ome]:-0}
  count=$((count + 1))
  ome_counts[$ome]=$count

  locus_file="${ome}_window${count}_locus.txt"

  echo "Processing $col1"
  echo "  First gene: $first_gene"
  echo "  Last gene: $last_gene"
  echo "  Output file: $locus_file"

  if [[ -n "$first_gene" && -n "$last_gene" ]]; then
    acc2locus -a "${first_gene} ${last_gene}" -b > "$locus_file"
  else
    echo "  Skipped $col1: One or both genes not found"
    continue
  fi
done

# ---- Step 5: Generate GBK files ----
echo "Generating GBK files..."

for txt_file in *_window*_locus.txt; do
  gbk_file="${txt_file%.txt}.gbk"
  acc2gbk -i "$txt_file" > "$gbk_file"
done

echo "Done."

