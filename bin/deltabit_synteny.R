#!/usr/bin/env Rscript

library(geneviewer)
library(Biostrings)
library(pwalign)
library(parallel)
library(argparse)
library(dplyr)
library(webshot2)

parser = ArgumentParser()
parser$add_argument("-f", "--folder", help = "Path to folder containing gbk files", required = TRUE)
parser$add_argument("-m", "--minimum", help = "Percent identity similarity minimum for homology", default = 30)
parser$add_argument("-q", "--query", help = "Name of focal gbk to align to", required = TRUE)
parser$add_argument("-w", "--width", help = "Width in pixels to export plot in", default = "800px")
parser$add_argument("-H", "--height", help = "Width in pixels to export plot in", default = "1000px")
parser$add_argument("-p", "--prefix", help = "Prefix for output files", required = TRUE)

args = parser$parse_args()

folder_path = args$folder
minimum_n = as.numeric(args$minimum)

BlastP_results = geneviewer::protein_blast(
  folder_path,
  query = args$query,
  id = "protein_id", 
  cluster = "cluster",
  identity = minimum_n,
  parallel = FALSE
)

BlastP_results = BlastP_results %>%
  mutate(sense = ifelse(grepl("complement", region), "complement", "forward"))

plot = GC_chart(
  data = BlastP_results,
  cluster = "cluster",
  strand = "sense",
  group = "BlastP",
  width = args$width,
  height = args$height
) %>%
  #GC_labels(label = "protein_id", cluster = 1) %>%
  GC_scale(axis_type = "range") %>%
  GC_links(group = "BlastP", measure = "identity", labelStyle = list(fontSize = "9px")) %>%
  GC_clusterLabel() %>%
  GC_legend(FALSE)

# Set up file name
filename_output = paste0(args$prefix, "_deltabit_synteny_plot.html")

# Save as html
htmlwidgets::saveWidget(plot, filename_output, selfcontained = TRUE)
