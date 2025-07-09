#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(slider)
library(ggExtra)
library(tidyr)
library(ggpubr)
library(argparse)

# Set up argument parser
parser = ArgumentParser()
parser$add_argument("-p", "--prefix", help = "Prefix for output files", required = TRUE)
parser$add_argument("-d", "--deltabit", help = "Path to whole contig deltabit output tsv", required = TRUE)
parser$add_argument("-b", "--busco", help = "Path to the busco deltabit output tsv", required = TRUE)
parser$add_argument("-S", "--smooth", help = "How many surrounding genes on each side to average when smoothing deltabit.", default = 10)
parser$add_argument("-s", "--bit_smooth", help = "How many surrounding genes on each side to average when smoothing raw bitcores.", default = 3)
parser$add_argument("-t", "--threshold", help = "Significance threshold at which to call a significant peak, derived from normal distribution of buscos", default = .01)

# Parse the command-line arguments
args = parser$parse_args()

# Assign command_line arguments to variables

prefix = args$prefix
df_path = args$deltabit
smooth_n = as.integer(args$smooth)
bit_smooth_n = as.integer(args$bit_smooth)
threshold_n = as.numeric(args$threshold)

df = read.delim(df_path, sep = "\t", header = TRUE)

busco_df = read.delim(args$busco, sep = "\t", header = TRUE)

mean_bit = mean(busco_df$delta_bit, na.rm = TRUE)
sd_bit = sd(busco_df$delta_bit, na.rm = TRUE)

busco_p_01 = qnorm(.01, mean = mean_bit, sd = sd_bit)
busco_p_05 = qnorm(.05, mean = mean_bit, sd = sd_bit)

busco_df$log_bit <- sign(busco_df$delta_bit) * log1p(abs(busco_df$delta_bit)) / log(2)

mean_log_bit = mean(busco_df$log_bit, na.rm = TRUE)
sd_log_bit = sd(busco_df$log_bit, na.rm = TRUE)

busco_log_p_01 = qnorm(threshold_n, mean = mean_log_bit, sd = sd_log_bit)



# Sort data by distance_from_captain
df <- df %>% arrange(distance_from_captain)

df$log_ingroup = sign(df$ingroup_bitscore) * log1p(abs(df$ingroup_bitscore)) / log(2)

df$log_outgroup = sign(df$outgroup_bitscore) * log1p(abs(df$outgroup_bitscore)) / log(2)

df_upper_log_smoothed_in <- df %>%
  mutate(log_ingroup_smooth = slide_dbl(
    .x = log_ingroup,
    .before = bit_smooth_n,
    .after = bit_smooth_n,
    .f = mean,
    .complete = FALSE  # Ensures no excessive NA values at the edges
  ))

df_upper_log_smoothed_both = df_upper_log_smoothed_in %>%
  mutate(log_outgroup_smooth = slide_dbl(
    .x = log_outgroup,
    .before = bit_smooth_n,
    .after = bit_smooth_n,
    .f = mean, 
    complete = FALSE
  ))

df_long = df_upper_log_smoothed_both %>%
  pivot_longer(cols = c(log_ingroup_smooth, log_outgroup_smooth), 
               names_to = "group",
               values_to = "log_value")

rawplot = ggplot(df_long, aes(x = distance_from_captain / 1e3, y = log_value, color = group)) +
  geom_line() + 
  guides(color = "none") +
  scale_color_manual(values = c(
    "log_ingroup_smooth" = "#CC79A7",
    "log_outgroup_smooth" = "#56B4E9"
  )) +
  theme_pubr() + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_line(color = "black", size = .1), axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5))

# Distribution of busco deltabit to make sure it's normal
svg(paste0(prefix, "_busco_log_deltabit_distribution.svg"))
hist(busco_df$log_bit, breaks = 50)
dev.off()

df$log_bit <- sign(df$delta_bit) * log1p(abs(df$delta_bit)) / log(2)

# Sort data by distance_from_captain
df <- df %>% arrange(distance_from_captain)

# Compute rolling mean with dynamic window based on distance
df_log_rolling_filtered <- df %>%
  mutate(avg_log_delta_bit = slide_dbl(
    .x = log_bit,
    .before = smooth_n,
    .after = smooth_n,
    .f = mean,
    .complete = FALSE  # Ensures no excessive NA values at the edges
  ))

scatter_plot = ggplot(df_log_rolling_filtered, aes(x = distance_from_captain / 1e3, y = avg_log_delta_bit, color = log_bit < busco_log_p_01)) +
  geom_point(alpha = .5, size = 2.5) + 
  geom_hline(yintercept = c(busco_log_p_01), color = "red", linetype = "dashed", linewidth = .7) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
  guides(color = "none") +
  labs(
    x = "Distance from Captain (kb)",
    y = expression("Average log"[2]*"Deltabit")
  ) +
  theme_minimal()


filename_no_peak = paste0(args$prefix, "_deltabit_no_peakcalling.svg")

svg(filename_no_peak)
ggarrange(rawplot, scatter_plot, nrow = 2, align = "v", widths = c(5, 1), heights = c(.2, 1))
dev.off()

# PEAKCALLING

df <- df_log_rolling_filtered  # just to keep things short

# Step 1: Find the gene closest to distance_from_captain == 0
start_index <- which.min(abs(df$distance_from_captain))

# Step 2: Check if that gene is part of a cluster
if (df$avg_log_delta_bit[start_index] < busco_log_p_01) {
  
  # Expand backward from start_index
  cluster_start <- start_index
  for (i in (start_index - 1):1) {
    if (df$avg_log_delta_bit[i] < busco_log_p_01) {
      cluster_start <- i
    } else {
      break
    }
  }
  
  # Expand forward from start_index
  cluster_end <- start_index
  for (i in (start_index + 1):nrow(df)) {
    if (df$avg_log_delta_bit[i] < busco_log_p_01) {
      cluster_end <- i
    } else {
      break
    }
  }
  
} else {
  # Otherwise, find cluster in either direction like before
  
  find_cluster_start <- function(df, start_index, direction, threshold) {
    if (direction == "forward") {
      for (i in start_index:nrow(df)) {
        if (df$avg_log_delta_bit[i] < threshold) return(i)
      }
    } else {
      for (i in start_index:1) {
        if (df$avg_log_delta_bit[i] < threshold) return(i)
      }
    }
    return(NA)
  }
  
  forward_start <- find_cluster_start(df, start_index, "forward", busco_log_p_01)
  backward_start <- find_cluster_start(df, start_index, "backward", busco_log_p_01)
  
  dist_forward <- if (!is.na(forward_start)) abs(forward_start - start_index) else Inf
  dist_backward <- if (!is.na(backward_start)) abs(backward_start - start_index) else Inf
  
  if (dist_forward < dist_backward) {
    cluster_start <- forward_start
    cluster_end <- forward_start
    for (i in (cluster_start + 1):nrow(df)) {
      if (df$avg_log_delta_bit[i] < busco_log_p_01) {
        cluster_end <- i
      } else {
        break
      }
    }
  } else if (dist_backward < Inf) {
    cluster_end <- backward_start
    cluster_start <- backward_start
    for (i in (cluster_start - 1):1) {
      if (df$avg_log_delta_bit[i] < busco_log_p_01) {
        cluster_start <- i
      } else {
        break
      }
    }
  } else {
    stop("No cluster found in either direction.")
  }
}

# Step 3: Extend 5 genes beyond both sides
final_start <- max(cluster_start - 5, 1)
final_end <- min(cluster_end + 5, nrow(df))

# Step 4: Extract all rows from final_start to final_end
extracted_genes <- df[final_start:final_end, ]


# Export gene_name column to a text file
writeLines(extracted_genes$gene_name, paste0(prefix, "_htr.txt"))


extracted_max = max(extracted_genes$distance_from_captain)
extracted_min = min(extracted_genes$distance_from_captain)

#Add extrcted genes to scatter plot
scatter_extracted = ggplot(df_log_rolling_filtered, aes(x = distance_from_captain / 1e3, y = avg_log_delta_bit, color = log_bit < busco_log_p_01)) +
  geom_rect(xmin = extracted_min / 1e3, xmax = extracted_max / 1e3, ymin = -Inf, ymax = Inf, alpha = .3, fill = NA, color = "#AFDEC7", linewidth = 1.3) +
  geom_point(alpha = .5, size = 2.5) + 
  geom_hline(yintercept = c(busco_log_p_01), color = "red", linetype = "dashed", linewidth = .7) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"))+
  guides(color = "none") +
  labs(
    x = "Distance from Captain (kb)",
    y = expression("Average log"[2]*"Deltabit")
  ) +
  theme_minimal()

filename_peak = paste0(args$prefix, "_deltabit_peakcalling.svg")
svg(filename_peak)
ggarrange(rawplot, scatter_extracted, nrow = 2, align = "v", widths = c(5, 1), heights = c(.2, 1))
dev.off()

#filter genes from extracted_genes with significant log_bits for blasting
sig_genes = extracted_genes[extracted_genes$log_bit < busco_log_p_01, ]

#write sig_genes to a text file
writeLines(sig_genes$gene_name, paste0(prefix, "_sig_genes.txt"))
