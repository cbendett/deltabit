#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(argparse)

# Set up argument parser
parser = ArgumentParser()
parser$add_argument("-p", "--prefix", help = "Prefix for output files", required = TRUE)
parser$add_argument("-f", "--file", help = "Path to input blast TSV file", required = TRUE)
parser$add_argument("-t", "--top", help = "How many of the genomes with the most blast hits to export", default = 10)
parser$add_argument("-m", "--minimum", help = "Minimum number of blast hits required for a genome to exported", default = 3)
parser$add_argument("-b", "--bitscore", help = "Minimum bitscore to count as a significant hit", default = 40)
parser$add_argument("-s", "--sorted", help = "How many of the sorted top blast hits to report per query gene", default = 20)

# Parse the command-line arguments
args = parser$parse_args()

# Assign command_line arguments to variables

prefix = args$prefix
df_path = args$file
top_n = as.integer(args$top)
minimum_n = as.integer(args$minimum)
bitscore_n = as.numeric(args$bitscore)
sorted_n = as.integer(args$sorted)

# Read blast results
df = read.delim(df_path, sep = "\t", header = FALSE)
colnames(df) = c("qseqid", "sseqid", "pident", "ppos", "sstart", "send", "evalue", "bitscore")
df$ome = sub("_.*", "", df[[2]])

df_summary <- df %>%
  filter(bitscore >= bitscore_n) %>%
  group_by(qseqid) %>%
  slice_max(bitscore, n = sorted_n, with_ties = TRUE) %>%
  ungroup() %>%
  select(qseqid, ome) %>%
  distinct() %>%
  arrange(ome, qseqid) %>%
  group_by(ome) %>%
  summarise(
    qseqid_group_count = n(),
    qseqids = toString(qseqid)  # Creates a comma-separated string of qseqid names
  ) %>%
  left_join(df %>% select(ome) %>% distinct(), by = "ome") %>%
  arrange(desc(qseqid_group_count))


# Write list of best omes to text file
top_omes = df_summary %>%
  filter(qseqid_group_count >= minimum_n) %>%
  slice_head(n=top_n) %>%
  pull(ome)

write_lines(top_omes, "ome_list.txt")

for (ome_val in top_omes) {
  df_filtered <- df %>%
    filter(ome == ome_val) %>%
    group_by(qseqid) %>%
    slice_max(bitscore, with_ties = FALSE) %>%
    ungroup()
  
  # Assign the result to a new variable named after the ome value
  assign(ome_val, df_filtered)
}

for (ome_val in top_omes) {
  df_current = get(ome_val)
  sseqids = df_current$sseqid
  
  # Write to a text file named after the ome value
  writeLines(as.character(sseqids), paste0(prefix, "_", ome_val, "_homologs.txt"))
}
