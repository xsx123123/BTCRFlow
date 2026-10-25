#!/usr/bin/env Rscript
# vj_pairing_plot.R
# V-J pairing chord diagram using circlize.
# Usage: Rscript vj_pairing_plot.R <sample_id> <chain> <input_dir> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript vj_pairing_plot.R <sample_id> <chain> <input_dir> <output_dir>")
}

sample_id <- args[1]
chain <- args[2]
input_dir <- args[3]
output_dir <- args[4]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(circlize)
library(dplyr)
library(readr)

input_file <- file.path(input_dir, paste0(sample_id, ".", chain, ".downsampled.tsv"))
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

clones <- read_tsv(input_file, show_col_types = FALSE)

# Extract V and J gene names
if ("allVHitsWithScore" %in% colnames(clones)) {
  v_col <- "allVHitsWithScore"
} else if ("allVHits" %in% colnames(clones)) {
  v_col <- "allVHits"
} else {
  stop("No V-gene column found")
}

if ("allJHitsWithScore" %in% colnames(clones)) {
  j_col <- "allJHitsWithScore"
} else if ("allJHits" %in% colnames(clones)) {
  j_col <- "allJHits"
} else {
  stop("No J-gene column found")
}

vj_matrix <- clones %>%
  mutate(
    v_gene = gsub("\\*.*", "", !!sym(v_col)),
    j_gene = gsub("\\*.*", "", !!sym(j_col))
  ) %>%
  filter(v_gene != "", j_gene != "") %>%
  count(v_gene, j_gene, wt = cloneFraction, name = "frequency") %>%
  pivot_wider(names_from = j_gene, values_from = frequency, values_fill = 0)

mat <- as.matrix(vj_matrix[, -1])
rownames(mat) <- vj_matrix$v_gene

# Filter low-frequency pairs
mat[mat < 0.001] <- 0

pdf(
  file.path(output_dir, paste0(sample_id, "_", chain, "_VJ_pairing.pdf")),
  width = 10, height = 10
)

circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0))
chordDiagram(
  mat,
  annotationTrack = c("grid"),
  preAllocateTracks = list(
    list(track.height = 0.05)
  ),
  transparency = 0.3
)
title(paste0(sample_id, " - ", chain, " V-J Pairing"))

dev.off()

cat(paste0("V-J pairing plot generated for ", sample_id, " ", chain, "\n"))
