#!/usr/bin/env Rscript
# spectratype.R
# V-gene usage, CDR3 length distribution, and germline divergence plots.
# Usage: Rscript spectratype.R <sample_id> <chain> <input_dir> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript spectratype.R <sample_id> <chain> <input_dir> <output_dir>")
}

sample_id <- args[1]
chain <- args[2]
input_dir <- args[3]
output_dir <- args[4]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Read clone table
input_file <- file.path(input_dir, paste0(sample_id, ".", chain, ".downsampled.tsv"))
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

clones <- read_tsv(input_file, show_col_types = FALSE)

# --- V-gene usage plot ---
if ("allVHitsWithScore" %in% colnames(clones)) {
  v_col <- "allVHitsWithScore"
} else if ("allVHits" %in% colnames(clones)) {
  v_col <- "allVHits"
} else {
  v_col <- NULL
}

if (!is.null(v_col)) {
  v_genes <- clones %>%
    mutate(v_gene = gsub("\\*.*", "", !!sym(v_col))) %>%
    mutate(v_family = gsub("[-/].*", "", v_gene)) %>%
    count(v_family, wt = cloneFraction, name = "frequency") %>%
    arrange(desc(frequency))

  p_vgene <- ggplot(v_genes, aes(x = reorder(v_family, frequency), y = frequency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = paste0(sample_id, " - ", chain, " V-gene Usage"),
      x = "V-gene Family",
      y = "Clone Fraction"
    ) +
    theme_minimal()

  ggsave(
    file.path(output_dir, paste0(sample_id, "_", chain, "_Vgene_usage.pdf")),
    p_vgene, width = 10, height = 8
  )
}

# --- CDR3 length distribution ---
if ("aaSeqCDR3" %in% colnames(clones)) {
  cdr3_len <- clones %>%
    filter(!is.na(aaSeqCDR3), aaSeqCDR3 != "") %>%
    mutate(cdr3_length = nchar(aaSeqCDR3)) %>%
    count(cdr3_length, wt = cloneFraction, name = "frequency")

  p_cdr3 <- ggplot(cdr3_len, aes(x = cdr3_length, y = frequency)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    labs(
      title = paste0(sample_id, " - ", chain, " CDR3 Length Distribution"),
      x = "CDR3 Length (aa)",
      y = "Clone Fraction"
    ) +
    theme_minimal()

  ggsave(
    file.path(output_dir, paste0(sample_id, "_", chain, "_CDR3_length.pdf")),
    p_cdr3, width = 10, height = 6
  )
}

# --- Germline divergence (SHM for BCR) ---
if ("allVHitsWithScore" %in% colnames(clones) && "vGeneMutations" %in% colnames(clones)) {
  divergence <- clones %>%
    filter(!is.na(vGeneMutations)) %>%
    count(vGeneMutations, wt = cloneFraction, name = "frequency")

  p_div <- ggplot(divergence, aes(x = vGeneMutations, y = frequency)) +
    geom_bar(stat = "identity", fill = "firebrick") +
    labs(
      title = paste0(sample_id, " - ", chain, " Germline Divergence"),
      x = "V-gene Mutations",
      y = "Clone Fraction"
    ) +
    theme_minimal()

  ggsave(
    file.path(output_dir, paste0(sample_id, "_", chain, "_divergence.pdf")),
    p_div, width = 10, height = 6
  )
}

cat(paste0("Plots generated for ", sample_id, " ", chain, "\n"))
