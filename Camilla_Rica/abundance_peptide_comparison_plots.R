# Five candidate proteins: total abundance, phosphopeptide detection, and FDR
# Place DE_total.xlsx and combined_site_STYH.tsv in your project folder.
# If your TSV is named combined_site_STYH(4).tsv, change phospho_file below.
# Install once: install.packages(c('readxl','readr','dplyr','tidyr','stringr','ggplot2','patchwork'))
#all done with chat gpt

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)


# EDIT THESE PATHS if needed
excel_file <- 'DE_total.xlsx'
phospho_file <- 'combined_site_STYH.tsv'
output_dir <- 'five_peptide_plots'
dir.create(output_dir, showWarnings = FALSE)

candidates <- tibble(
  id = c('ARG45065','ARG45665','ARG45823','ARG46060','ARG46318'),
  name = c('L7-L12','RsmH','Homoserine dehydrogenase','ScpA','Putative peptidase'),
  peptide = c('ANhEQIIEAIK','EVLKDHLhK','DLLAVhLK','LDAFNGPLDLLLhLIQK','HLVEQALNQAIDQINhADK')
)

# Total protein: normalized abundance matrix (already log2 scale)
total <- read_excel(excel_file, sheet = 'normalized_abundances_matrix') %>%
  filter(str_detect(protein_Id, paste(candidates$id, collapse = '|'))) %>%
  mutate(id = str_extract(protein_Id, 'ARG[0-9]+')) %>%
  select(id, matches('^(JE2|2anc_|2stk_|2stp_)')) %>%
  pivot_longer(-id, names_to = 'sample', values_to = 'log2_abundance') %>%
  mutate(group = case_when(
    str_detect(sample, 'pulseSN') ~ 'SNpulse',
    str_detect(sample, '_TSB_') ~ 'TSB',
    str_detect(sample, '_SN_') ~ 'SN',
    TRUE ~ NA_character_),
    strain = case_when(
      str_detect(sample, 'anc') ~ 'anc',
      str_detect(sample, 'stk') ~ 'stk',
      str_detect(sample, 'stp') ~ 'stp',
      TRUE ~ NA_character_),
    condition = factor(paste(strain, group, sep = ' '),
      levels = c('anc TSB','stk TSB','stp TSB',
                 'anc SN','stk SN','stp SN',
                 'anc SNpulse','stk SNpulse','stp SNpulse')))

# Phosphopeptide table: candidate rows and enriched intensity columns
phos <- read_tsv(phospho_file, show_col_types = FALSE)
phos_focus <- phos %>%
  mutate(id = str_extract(`Protein ID`, 'ARG[0-9]+')) %>%
  inner_join(candidates, by = 'id') %>%
  filter(str_detect(Peptide, fixed(peptide)))

if (nrow(phos_focus) != 5) {
  warning('Expected 5 phosphopeptide rows; found ', nrow(phos_focus),
          '. Check candidate peptide sequences and duplicates.')
}

# Intensity > 0 indicates detected phosphopeptide signal.
# This does NOT prove confident site localization in that replicate.
phos_long <- phos_focus %>%
  select(id, matches('_Enriched Intensity$')) %>%
  pivot_longer(-id, names_to = 'column', values_to = 'intensity') %>%
  mutate(sample = str_extract(column, 'JE2(?:anc|stk|stp)_(?:TSB|SN|pulseSN)_[123]'),
         group = case_when(
           str_detect(sample, 'pulseSN') ~ 'SNpulse',
           str_detect(sample, '_TSB_') ~ 'TSB',
           str_detect(sample, '_SN_') ~ 'SN',
           TRUE ~ NA_character_),
         strain = str_extract(sample, '(?<=JE2)(anc|stk|stp)'),
         condition = factor(paste(strain, group), levels = levels(total$condition)),
         detected = !is.na(intensity) & intensity > 0) %>%
  filter(!is.na(condition))

# Group detection, including 0/3 groups
all_conditions <- factor(levels(total$condition), levels = levels(total$condition))
detection <- phos_long %>%
  group_by(id, condition) %>%
  summarise(n_detected = sum(detected), n_measured = n(), .groups = 'drop') %>%
  complete(id = candidates$id, condition = all_conditions,
           fill = list(n_detected = 0, n_measured = 0))

# Differential expression results: FDR and signed log2 contrast
# 'diff' follows the order in 'contrast', e.g. TSB_vs_SN = TSB minus SN.
de <- read_excel(excel_file, sheet = 'diff_exp_analysis') %>%
  filter(str_detect(protein_Id, paste(candidates$id, collapse = '|'))) %>%
  mutate(id = str_extract(protein_Id, 'ARG[0-9]+')) %>%
  select(id, contrast, diff, FDR, conf.low, conf.high)

# Selected biologically relevant comparisons; change this list if desired.
chosen <- c('stp_vs_anc_at_TSB', 'stk_vs_anc_at_SN',
            'TSB_vs_SN_at_anc', 'TSB_vs_SN_at_stp')
de_focus <- de %>% filter(contrast %in% chosen)

make_plot <- function(target_id) {
  label <- candidates$name[match(target_id, candidates$id)]
  a <- total %>% filter(id == target_id)
  b <- detection %>% filter(id == target_id)
  c <- de_focus %>% filter(id == target_id)

  p_a <- ggplot(a, aes(condition, log2_abundance)) +
    geom_point(position = position_jitter(width = 0.12, height = 0),
               size = 2.2, alpha = 0.75, na.rm = TRUE) +
    stat_summary(fun = mean, geom = 'point', shape = 18, size = 3.5,
                 colour = 'steelblue', na.rm = TRUE) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = 'errorbar', width = 0.25, colour = 'steelblue',
                 na.rm = TRUE) +
    labs(title = 'A. Total protein', x = NULL, y = 'Normalized log2 abundance') +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_b <- ggplot(b, aes(condition, n_detected)) +
    geom_col(fill = 'steelblue', width = 0.7) +
    geom_text(aes(label = paste0(n_detected, '/', n_measured)),
              vjust = -0.4, size = 3.4) +
    scale_y_continuous(limits = c(0, 3.7), breaks = 0:3) +
    labs(title = 'B. Candidate phosphopeptide detected',
         x = NULL, y = 'Replicates with intensity > 0') +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_c <- ggplot(c, aes(x = diff, y = reorder(contrast, diff),
                       colour = !is.na(FDR) & FDR < 0.05)) +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey60') +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_point(size = 3) +
    geom_text(aes(label = paste0('FDR=', signif(FDR, 2))),
              hjust = -0.1, colour = 'black', size = 3, check_overlap = TRUE) +
    scale_colour_manual(values = c('TRUE' = 'forestgreen', 'FALSE' = 'grey55'),
                        guide = 'none') +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.40))) +
    labs(title = 'C. Total protein: differential abundance',
         x = 'Signed log2 difference (first group minus second)', y = NULL) +
    theme_minimal(base_size = 12)

  combined <- (p_a / p_b / p_c) +
    plot_annotation(title = paste0(label, ' (', target_id, ')'),
                    subtitle = 'Total-protein FDR does not test phosphorylation changes',
                    theme = theme(plot.title = element_text(hjust = 0.5,
                                                           face = 'bold', size = 17)))
  combined
}

plots <- lapply(candidates$id, make_plot)
names(plots) <- candidates$name

for (i in seq_along(plots)) {
  file_stub <- paste0(candidates$id[i], '_', gsub('[^A-Za-z0-9]+', '_', candidates$name[i]))
  ggsave(file.path(output_dir, paste0(file_stub, '.png')),
         plots[[i]], width = 11, height = 12, dpi = 300)
}

pdf(file.path(output_dir, 'five_protein_comparisons.pdf'), width = 11, height = 12)
for (p in plots) print(p)
dev.off()

write_csv(de_focus %>% left_join(candidates %>% select(id, name), by = 'id'),
          file.path(output_dir, 'selected_total_protein_FDR.csv'))

# Display ScpA in the RStudio Plots pane:
print(plots[['ScpA']])
