suppressPackageStartupMessages({
  library(phyloseq); library(microViz); library(tidyverse); library(cowplot)
})

setwd("/Users/nickolevillabona/Desktop/Ch1/SB_Pipeline")
phyloseqOBJ <- readRDS("data/processed/phyloseqOBJ.rds")

meta <- data.frame(sample_data(phyloseqOBJ), stringsAsFactors = FALSE) %>%
  rownames_to_column("SampleID")

# ── Split genera into two rows ──────────────────────────────────────────────
all_genera <- sort(unique(meta$Bee_Genus))
# Row 1: Austroplebeia → Scaptotrigona (52 samples)
# Row 2: Scaura → Trigonisca        (55 samples)
cutpoint   <- "Scaura"
row1_genera <- all_genera[all_genera < cutpoint]
row2_genera <- all_genera[all_genera >= cutpoint]

cat("Row 1 genera:", paste(row1_genera, collapse=", "), "\n")
cat("Row 2 genera:", paste(row2_genera, collapse=", "), "\n")

# ── Shared settings ─────────────────────────────────────────────────────────
genus_abbrev <- c(
  "Austroplebeia"  = "Aust", "Duckeola"       = "Duck",
  "Eomelipona"     = "Eome", "Frieseomelitta" = "Frme",
  "Geotrigona"     = "Geot", "Heterotrigona"  = "Hete",
  "Lestrimelitta"  = "Lest", "Melikerria"     = "Mlkr",
  "Meliponula"     = "Mlpl", "Michmelia"      = "Miml",
  "Nannotrigona"   = "Nann", "Oxytrigona"     = "Oxyt",
  "Sundatrigona"   = "Sund", "Tetragona"      = "Ttgo",
  "Tetragonisca"   = "Trgs", "Tetragonula"    = "Trgu",
  "Trigona"        = "Tgna", "Trigonisca"     = "Tgca",
  "Scaptotrigona"  = "Scap", "Scaura"         = "Scau"
)

region_cols <- c(
  "Neotropics"             = "grey80",
  "AfroTropics"            = "grey45",
  "Indo-Malay-Australasia" = "grey15"
)

# ── Helper: build one row ────────────────────────────────────────────────────
make_row <- function(ps_sub, genera_in_row, show_legend = FALSE) {
  meta_sub <- data.frame(sample_data(ps_sub), stringsAsFactors = FALSE) %>%
    rownames_to_column("SampleID")

  ordered_samples <- meta_sub %>%
    arrange(Bee_Genus, SampleID) %>%
    pull(SampleID)

  p_base <- ps_sub %>%
    comp_barplot(
      tax_level          = "Genus",
      n_taxa             = 10,
      other_name         = "Other",
      sample_order       = ordered_samples,
      taxon_renamer      = function(x) str_remove(x, " [ae]t rel."),
      palette            = distinct_palette(n = 10, add = "grey90"),
      merge_other        = FALSE,
      bar_outline_colour = "white"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    labs(x = NULL, y = "Relative Abundance") +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.ticks.x       = element_line(linewidth = 0.3),
      panel.grid         = element_blank(),
      legend.text        = element_text(face = "italic", size = 11),
      plot.margin        = margin(28, 5, 2, 5)
    )

  actual_order <- ggplot_build(p_base)$layout$panel_params[[1]]$x$limits

  sample_genus <- tibble(SampleID = actual_order) %>%
    left_join(meta_sub %>% select(SampleID, Bee_Genus, Region), by = "SampleID")

  genus_seq <- unique(sample_genus$Bee_Genus)
  order_with_gaps <- c()
  for (i in seq_along(genus_seq)) {
    samps <- sample_genus$SampleID[sample_genus$Bee_Genus == genus_seq[i]]
    order_with_gaps <- c(order_with_gaps, samps)
    if (i < length(genus_seq)) order_with_gaps <- c(order_with_gaps, paste0("GAP_", i, "a"), paste0("GAP_", i, "b"), paste0("GAP_", i, "c"))
  }

  xpos_df <- tibble(
    SampleID = order_with_gaps,
    x_pos    = seq_along(order_with_gaps),
    is_gap   = startsWith(order_with_gaps, "GAP_")
  ) %>% left_join(meta_sub %>% select(SampleID, Bee_Genus, Region), by = "SampleID")

  breaks <- xpos_df %>%
    filter(!is_gap) %>%
    group_by(Bee_Genus) %>%
    summarise(x_mid = mean(x_pos), .groups = "drop") %>%
    arrange(x_mid) %>%
    mutate(label = genus_abbrev[Bee_Genus])

  p_annot <- p_base +
    scale_x_discrete(limits = order_with_gaps, breaks = actual_order) +
    annotate("text",
             x = breaks$x_mid, y = 1.02,
             label = breaks$label,
             size = 4.0, angle = 0, hjust = 0.5, fontface = "italic") +
    coord_cartesian(clip = "off")

  annot_reg <- xpos_df %>%
    filter(!is_gap) %>%
    mutate(SampleID = factor(SampleID, levels = order_with_gaps))

  p_region <- ggplot(annot_reg, aes(x = SampleID, y = 1, fill = Region)) +
    geom_col(width = 1) +
    scale_x_discrete(limits = order_with_gaps) +
    scale_fill_manual(values = region_cols, name = "Region") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      legend.text     = element_text(size = 11),
      legend.title    = element_text(size = 11, face = "bold"),
      plot.margin     = margin(0, 5, 2, 5)
    )

  left_col <- plot_grid(
    p_annot  + theme(legend.position = "none"),
    p_region + theme(legend.position = "none"),
    ncol = 1, align = "v", axis = "lr",
    rel_heights = c(9, 0.5)
  )

  list(left = left_col, taxa_legend = get_legend(p_base), region_legend = get_legend(p_region))
}

# ── Build subsetted phyloseq objects ─────────────────────────────────────────
ps_row1 <- subset_samples(phyloseqOBJ, Bee_Genus %in% row1_genera)
ps_row2 <- subset_samples(phyloseqOBJ, Bee_Genus %in% row2_genera)

row1 <- make_row(ps_row1, row1_genera)
row2 <- make_row(ps_row2, row2_genera)

# ── Combine legends ───────────────────────────────────────────────────────────
legends <- plot_grid(
  row1$taxa_legend,
  row1$region_legend,
  ncol = 1, rel_heights = c(3, 1)
)

# ── Stack rows ────────────────────────────────────────────────────────────────
both_rows <- plot_grid(
  row1$left,
  row2$left,
  ncol = 1, align = "v"
)

fig <- plot_grid(both_rows, legends, ncol = 2, rel_widths = c(10, 1.5))

out <- "figures/raw/composition_barplot_2rows.pdf"
ggsave(out, plot = fig, device = "pdf", width = 16, height = 10)
cat("Guardado:", out, "\n")
