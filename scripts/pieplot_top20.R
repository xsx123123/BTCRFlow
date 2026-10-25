#!/usr/bin/env Rscript
# pieplot_top20.R
# Top-20 clonotype pie charts from downsampled clone tables.
# Usage: Rscript pieplot_top20.R <sample_id> <chain> <input_file> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript pieplot_top20.R <sample_id> <chain> <input_file> <output_dir>")
}

sample_id <- args[1]
chain <- args[2]
input_file <- args[3]
output_dir <- args[4]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(ggplot2)
library(dplyr)
library(readr)

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

clones <- read_tsv(input_file, show_col_types = FALSE)

if (!"aaSeqCDR3" %in% colnames(clones) || !"cloneFraction" %in% colnames(clones)) {
  stop("Required columns aaSeqCDR3 and cloneFraction not found")
}

top20 <- clones %>%
  filter(!is.na(aaSeqCDR3), aaSeqCDR3 != "") %>%
  arrange(desc(cloneFraction)) %>%
  head(20) %>%
  mutate(label = ifelse(nchar(aaSeqCDR3) > 15, paste0(substr(aaSeqCDR3, 1, 12), "..."), aaSeqCDR3))

others_fraction <- 1 - sum(top20$cloneFraction)
if (others_fraction > 0) {
  top20 <- bind_rows(top20, data.frame(
    aaSeqCDR3 = "Others",
    cloneFraction = others_fraction,
    label = "Others",
    stringsAsFactors = FALSE
  ))
}

p <- ggplot(top20, aes(x = "", y = cloneFraction, fill = label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(
    title = paste0(sample_id, " - ", chain, " Top 20 Clonotypes"),
    fill = "CDR3"
  ) +
  theme_void() +
  theme(legend.position = "right", legend.text = element_text(size = 7))

ggsave(
  file.path(output_dir, paste0(sample_id, "_", chain, "_top20_pie.pdf")),
  p, width = 12, height = 8
)

cat(paste0("Top-20 pie chart generated for ", sample_id, " ", chain, "\n"))
