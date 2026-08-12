# Regenerate sessionInfo.txt — a record of the exact package versions used to
# produce the results in the manuscript.
#
# Run from the repository root:  Rscript scripts/sessionInfo.R

pkgs <- c("tidyverse", "phyloseq", "microViz", "vegan", "decontam", "ape",
          "nlme", "phytools", "corHMM", "ggpubr", "cowplot", "RColorBrewer",
          "gridExtra", "reshape2", "vroom", "funrar", "smplot2")

missing <- character(0)
for (p in pkgs) {
  ok <- suppressWarnings(suppressPackageStartupMessages(
    require(p, character.only = TRUE, quietly = TRUE)))
  if (!ok) missing <- c(missing, p)
}

out <- file("sessionInfo.txt", "wt")
writeLines(c(
  "# Session info for SB_Pipeline",
  "# Package versions used to produce the results in the manuscript.",
  "# Generated with sessionInfo() after loading every package required by scripts 01-03.",
  "# Regenerate with: Rscript scripts/sessionInfo.R",
  ""), out)

if (length(missing)) {
  writeLines(c("## Declared in the scripts but NOT installed in this environment:",
               paste("#  -", missing),
               ""), out)
}

sink(out); print(sessionInfo()); sink(); close(out)

cat("Wrote sessionInfo.txt\n")
if (length(missing)) cat("Missing packages:", paste(missing, collapse = ", "), "\n")
