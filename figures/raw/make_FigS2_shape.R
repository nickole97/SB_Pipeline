suppressPackageStartupMessages({
  library(phyloseq); library(microViz); library(vegan)
  library(tidyverse); library(wesanderson); library(cowplot)
})

jb_dir <- "/Users/nickolevillabona/Library/Mobile Documents/com~apple~CloudDocs/Desktop/StinglesBees/MeliponiniSeqs/QiimeAnalisys"

# ── 1. Load data ──────────────────────────────────────────────────────────────
obj_joyce <- readRDS(file.path(jb_dir, "physeq_picosELN.rds"))

taxonomy_jb <- read_tsv(
  file.path(jb_dir, "taxa/data/taxonomy.tsv"),
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)

# ── 2. Taxonomy cleaning ──────────────────────────────────────────────────────
tax_df_jb <- taxonomy_jb %>%
  separate_wider_delim(
    Taxon, delim = ";",
    names = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    too_few = "align_start"
  ) %>%
  mutate(across(Domain:Species, ~ sub("^[a-zA-Z]__", "", trimws(.x)))) %>%
  mutate(FeatureID = trimws(FeatureID)) %>%
  distinct(FeatureID, .keep_all = TRUE) %>%
  as.data.frame()

rownames(tax_df_jb) <- tax_df_jb$FeatureID
desired <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_df_jb$Kingdom <- tax_df_jb$Domain
tax_df_jb$Domain  <- NULL
for (col in setdiff(desired, colnames(tax_df_jb))) tax_df_jb[[col]] <- NA_character_
tax_df_jb <- tax_df_jb[, desired, drop = FALSE]

ids <- taxa_names(obj_joyce)
missing_asvs <- setdiff(ids, rownames(tax_df_jb))
if (length(missing_asvs)) {
  old_tax <- as.data.frame(tax_table(obj_joyce), stringsAsFactors = FALSE)
  for (col in setdiff(colnames(tax_df_jb), colnames(old_tax))) old_tax[[col]] <- NA_character_
  for (col in setdiff(colnames(old_tax), colnames(tax_df_jb))) tax_df_jb[[col]] <- NA_character_
  comunes <- intersect(rownames(old_tax), rownames(tax_df_jb))
  old_tax[comunes, colnames(tax_df_jb)] <- tax_df_jb[comunes, colnames(tax_df_jb)]
  new_mat <- as.matrix(old_tax[ids, , drop = FALSE])
} else {
  new_mat <- as.matrix(tax_df_jb[ids, , drop = FALSE])
}
storage.mode(new_mat) <- "character"
tax_table(obj_joyce)  <- tax_table(new_mat)

# LCA corrections + Giliamella typo fix
tax <- tax_table(obj_joyce)
lca_corrections <- c(
  "2f5bd1cf50b6d70e150bc1829af4ee7e" = "Gilliamella",
  "9c489be807c8d254d8af389af7b2ad60" = "Snodgrassella",
  "c1d35f5da6018f72638c83306839f516" = "Gilliamella",
  "07dc3714ead232d81471f5fa6f65f4a8" = "Gilliamella"
)
for (asv in names(lca_corrections)) {
  if (asv %in% rownames(tax)) tax[asv, "Genus"] <- lca_corrections[asv]
}
tax[tax[, "Genus"] %in% "Giliamella", "Genus"] <- "Gilliamella"
tax_table(obj_joyce) <- tax

# ── 3. Filter & rarefy ───────────────────────────────────────────────────────
set.seed(123)
physeq_clean <- rarefy_even_depth(
  obj_joyce, sample.size = 2500,
  rngseed = 123, replace = FALSE, trimOTUs = TRUE, verbose = FALSE
) %>%
  subset_samples(!is.na(Organism)) %>%
  tax_fix() %>%
  phyloseq_validate() %>%
  tax_fix(unknowns = c("uncultured"))

# ── 4. CLR transform ─────────────────────────────────────────────────────────
physeq_clr <- physeq_clean %>%
  tax_transform(trans = "clr", rank = "Genus")

# ── 5. Compute Gilliamella P/A and add to sample_data ────────────────────────
physeq_ang_raw <- subset_samples(physeq_clean, Organism == "Tetragonisca_angustula")

ps_rel_ang  <- transform_sample_counts(physeq_ang_raw, function(x) x / sum(x))
ps_gill_ang <- ps_rel_ang %>%
  subset_taxa(Genus == "Gilliamella") %>%
  tax_glom(taxrank = "Genus")

gill_df <- psmelt(ps_gill_ang) %>%
  select(Sample, Abundance) %>%
  rename(SampleID = Sample, gill_abun = Abundance)

meta_ang <- data.frame(sample_data(physeq_ang_raw), stringsAsFactors = FALSE)
meta_ang$SampleID <- rownames(meta_ang)
meta_ang <- meta_ang %>%
  left_join(gill_df, by = "SampleID") %>%
  mutate(
    gill_abun      = replace_na(gill_abun, 0),
    Gilliamella_PA = ifelse(gill_abun > 0.01, "Present", "Absent")
  )
rownames(meta_ang) <- meta_ang$SampleID

# Add P/A to CLR-transformed object too (for PCA)
physeq_angustula <- subset_samples(physeq_clr, Organism == "Tetragonisca_angustula")
meta_ang_clr <- data.frame(sample_data(physeq_angustula), stringsAsFactors = FALSE)
meta_ang_clr$SampleID <- rownames(meta_ang_clr)
meta_ang_clr <- meta_ang_clr %>%
  left_join(meta_ang %>% select(SampleID, gill_abun, Gilliamella_PA), by = "SampleID")
rownames(meta_ang_clr) <- meta_ang_clr$SampleID
sample_data(physeq_angustula) <- sample_data(meta_ang_clr[sample_names(physeq_angustula), ])

physeq_orizabaensis <- subset_samples(physeq_clr, Organism == "Partamona_orizabaensis")

cat("T. angustula samples:", nsamples(physeq_angustula), "\n")
cat("Gilliamella P/A table:\n")
print(table(meta_ang_clr$Gilliamella_PA))

# ── 6. PCA ordinations ────────────────────────────────────────────────────────
pca_all         <- ord_calc(physeq_clr,         method = "PCA")
pca_angustula   <- ord_calc(physeq_angustula,   method = "PCA")
pca_orizabaensis <- ord_calc(physeq_orizabaensis, method = "PCA")

sp_colors_org <- c(
  "Partamona_orizabaensis"         = "#FC4E07",
  "Scaptotrigona_subobscuripennis" = "#1B7B3D",
  "Tetragonisca_angustula"         = "#00AFBB",
  "Tetragona_zieglieri"            = "#9d4edd",
  "Trigona_fulviventris"           = "#E7B800"
)

# Panel A: all species
p_pca_all <- ord_plot(pca_all, color = "Organism", size = 3, alpha = 0.7) +
  scale_color_manual(values = sp_colors_org,
                     labels = function(x) gsub("_", " ", x)) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_text(face = "italic")) +
  labs(color = "Species")

# Panel B: T. angustula — color = colony, shape = Gilliamella P/A (●=Present ○=Absent)
colony_cols <- wes_palette(n = 4, name = "AsteroidCity3")

p_angustula <- ord_plot(pca_angustula,
                        colour = "ColonyName",
                        shape  = "Gilliamella_PA",
                        size = 3.5, alpha = 0.8) +
  scale_color_manual(values = colony_cols, name = "Colony") +
  scale_shape_manual(values = c("Present" = 16, "Absent" = 1),
                     name = "Gilliamella") +
  theme_classic(base_size = 10) +
  ggtitle(expression(italic("Tetragonisca angustula"))) +
  guides(
    color = guide_legend(override.aes = list(shape = 16)),
    shape = guide_legend(override.aes = list(color = "grey30", size = 3))
  )

# Panel C: P. orizabaensis — color = colony only
p_orizabaensis <- ord_plot(pca_orizabaensis, colour = "ColonyName",
                           size = 3.5, alpha = 0.7) +
  scale_color_manual(values = c("#FC4E07", "#ffa987")) +
  theme_classic(base_size = 10) +
  ggtitle(expression(italic("Partamona orizabaensis"))) +
  labs(color = "Colony")

# ── 7. Assemble figure ────────────────────────────────────────────────────────
right_col <- plot_grid(p_angustula, p_orizabaensis,
                       ncol = 1, align = "hv", axis = "tblr",
                       labels = c("B", "C"), label_size = 14)

fig_s2 <- plot_grid(p_pca_all + theme(legend.position = "right"),
                    right_col,
                    labels = c("A", ""), rel_widths = c(1.6, 1), align = "h")

out_raw  <- "/Users/nickolevillabona/Desktop/Ch1/SB_Pipeline/figures/raw/FigS2_ordination_shape.pdf"
out_figs <- "/Users/nickolevillabona/Desktop/Ch1/resubmission/Figs/FigS2_ordination.pdf"
ggsave(out_raw,  plot = fig_s2, device = "pdf", width = 32, height = 18, units = "cm")
ggsave(out_figs, plot = fig_s2, device = "pdf", width = 32, height = 18, units = "cm")
cat("Guardado:", out_raw, "\n")
cat("Guardado:", out_figs, "\n")
