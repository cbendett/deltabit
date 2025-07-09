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

# ---- Step 2: Cluster genes into non-overlapping maximal windows ----
echo "Clustering genes with window size $WINDOW_SIZE..."

awk -v window="$WINDOW_SIZE" '
BEGIN {
    OFS = "\t"
    print "window_id", "scaffold", "start", "end", "gene_count", "gene_names"
}

{
    scaffold[NR] = $1
    start[NR] = $2
    end[NR] = $3
    gene[NR] = $4
    total = NR
}

END {
    window_num = 1
    i = 1
    while (i <= total) {
        cluster_start = start[i]
        cluster_end = end[i]
        cluster_genes = gene[i]
        gene_count = 1
        j = i + 1
        while (j <= total && scaffold[j] == scaffold[i] && start[j] - end[j - 1] <= window) {
            cluster_end = (end[j] > cluster_end ? end[j] : cluster_end)
            cluster_genes = cluster_genes "," gene[j]
            gene_count++
            j++
        }
        if (gene_count >= 2) {
            print "window_" window_num++, scaffold[i], cluster_start, cluster_end, gene_count, cluster_genes
        }
        i = j
    }
}
' sorted_genes.tsv > gene_windows.tsv

# ---- Step 2.5: Merge overlapping windows if any genes are shared ----
echo "Merging overlapping windows if they share genes..."

awk -F'\t' 'NR > 1 {
    window_id = $1
    scaffold = $2
    start = $3
    end = $4
    count = $5
    split($6, genes, ",")
    for (i in genes) {
        g = genes[i]
        window_for_gene[g] = (window_for_gene[g] ? window_for_gene[g] "," window_id : window_id)
    }
    window_info[window_id] = $0
}

END {
    # Group windows that share genes
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (g in window_for_gene) {
        split(window_for_gene[g], ws, ",")
        # Create union of all windows sharing at least one gene
        n = asort(ws)
        for (i = 1; i <= n; i++) {
            group[ws[i]]
        }
        key = ""
        for (i = 1; i <= n; i++) {
            if (ws[i] in grouped) continue
            key = key ? key "," ws[i] : ws[i]
            grouped[ws[i]]
        }
        if (key) merge_groups[key]
    }

    if (length(merge_groups) == 0) {
        # No overlaps, just print original windows
        for (w in window_info) print window_info[w]
        next
    }

    id = 1
    seen = ""
    for (gk in merge_groups) {
        split(gk, members, ",")
        merged_genes = ""
        merged_start = ""
        merged_end = ""
        merged_scaffold = ""
        gene_set[""] = ""
        delete gene_set

        for (i in members) {
            win = members[i]
            if (index(seen, win)) continue
            seen = seen "," win
            split(window_info[win], fields, "\t")
            if (!merged_scaffold) merged_scaffold = fields[2]
            s = fields[3]
            e = fields[4]
            split(fields[6], gs, ",")
            for (j in gs) gene_set[gs[j]] = 1
            if (!merged_start || s < merged_start) merged_start = s
            if (!merged_end || e > merged_end) merged_end = e
        }

        n = asorti(gene_set, sorted_genes)
        for (i = 1; i <= n; i++) merged_genes = (merged_genes ? merged_genes "," sorted_genes[i] : sorted_genes[i])

        print "merged_window_" id++, merged_scaffold, merged_start, merged_end, n, merged_genes
    }
}
' gene_windows.tsv > merged_gene_windows.tsv


# ---- Step 3: Extract first and last gene of each window ----
echo "Extracting first and last genes of each window..."

awk -F'\t' 'NR > 1 {
    n = split($6, parts, ",");
    print $1 "\t" parts[1] "\t" parts[n];
}' merged_gene_windows.tsv > processed_genes.txt

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
  # Extract the OME (prefix before _window) and window number
  base_name=$(basename "$txt_file" .txt)
  ome=$(echo "$base_name" | cut -d'_' -f1)
  num=$(echo "$base_name" | grep -oP '(?<=_window)\d+')
  gbk_file="${ome}_${num}.gbk"

  acc2gbk -i "$txt_file" > "$gbk_file"
done

echo "Done."

