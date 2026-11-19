# first idea
#FINAL plot 
# library(dplyr)
# library(readr)
# library(ggplot2)
# 
# 
# library(here)
# source(here("R", "setup.R"))
# data <- read_tsv(
#   here("Camilla_Rica", "combined_site_STYH.tsv"))
# 
# 
# library(ggplot2)
# 
# # Create table 
# filter_data <- data.frame(
#   Step = factor(
#     c("Histidine phosphorylation",
#       "LP >= 0.75",
#       "No S/T/Y on peptide",
#       "Detected in >=2/3 replicates"),
#     levels = c(
#       "Histidine phosphorylation",
#       "LP >= 0.75",
#       "No S/T/Y on peptide",
#       "Detected in >=2/3 replicates"
#     )
#   ),
#   Count = c(
#     nrow(his_all),
#     nrow(his_075),
#     nrow(his_only),
#     5
#   )
# )
# 
# # 
# ggplot(filter_data, aes(x = Step, y = Count)) +
#   geom_col(fill = "steelblue", width = 0.6) +
#   geom_text(aes(label = Count),
#             vjust = -0.5, size = 5) +
#   labs(
#     title = "Filtering potential histidine phosphorylation sites",
#     x = NULL,
#     y = "Number of peptide candidates"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 25, hjust = 1)
#   )
# 
# 
# 
# #size good for the poster 
# 
# ggplot(filter_data, aes(x = Step, y = Count)) +
#   geom_col(fill = "steelblue", width = 0.6) +
#   geom_text(aes(label = Count),
#             vjust = -0.5, size = 5) +
#   labs(
#     title = "Filtering potential histidine phosphorylation sites",
#     x = NULL,
#     y = "Number of peptide candidates"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, size = 18, face = "bold"
#     ),
#     axis.text.x = element_text(
#       angle = 25, hjust = 1, size = 14
#     ),
#     axis.text.y = element_text(size = 14),
#     axis.title.y = element_text(size = 16)
#   )
# 







# FINAL DIAGNOSTIC PLOT: HISTIDINE PHOSPHORYLATION

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(here)
library(tidyr)

# 1. Load data
data <- read_tsv(
  here("Camilla_Rica", "combined_site_STYH.tsv"))

# 2. Identify all histidine phosphorylation sites

his_all_check <- data %>%
  filter(str_detect(Index, "_H\\d+$")) %>%
  distinct(Index, .keep_all = TRUE)

# 3. Check and remove contaminants

contaminant_check <- his_all_check %>%
  filter(
    str_detect(Protein, regex("contam_", ignore_case = TRUE)))

his_no_contam <- his_all_check %>%
  filter(
    !str_detect(Protein, regex("contam_", ignore_case = TRUE)))

# 4. Prepare evidence table from the earlier analysis

his_evidence <- his %>%
  group_by(site_id) %>%
  summarise(
    best_loc_prob = if (all(is.na(best_loc_prob))) {
      NA_real_
    } else {
      max(best_loc_prob, na.rm = TRUE)
    },
    no_alt_acceptor = all(no_alt_acceptor),
    seen_in_STY_search = all(seen_in_STY_search),
    max_reps_in_group = max(max_reps_in_group, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(Index = site_id)

# 5. Start with all His sites without contaminants

filter_start <- his_no_contam %>%
  left_join(his_evidence, by = "Index")

# 6. Apply the four confidence filters sequentially

# Filter 1: localisation probability >= 0.75

his_075_check <- filter_start %>%
  filter(best_loc_prob >= 0.75)

# Filter 2: no alternative S/T/Y acceptor

his_only_check <- his_075_check %>%
  filter(no_alt_acceptor == TRUE)

# Filter 3: peptide not explained by the STY search

his_not_sty_check <- his_only_check %>%
  filter(seen_in_STY_search == FALSE)

# Filter 4: reproducible detection in at least 2 out of 3 replicates in one group

his_reproducible_check <- his_not_sty_check %>%
  filter(max_reps_in_group >= 2)

# 7. Check duplicates

duplicate_check <- his_only_check %>%
  count(Index) %>%
  filter(n > 1)

# 8. Print diagnostic results

cat("\n========== HISTIDINE SITE ANALYSIS ==========\n")

cat("\nAll His sites:", n_distinct(his_all_check$Index))

cat("\nContaminants:", n_distinct(contaminant_check$Index))

cat("\nHis sites after contaminant removal:",
    n_distinct(his_no_contam$Index))

cat("\nAfter localisation probability >= 0.75:",
    n_distinct(his_075_check$Index))

cat("\nAfter removing alternative S/T/Y acceptors:",
    n_distinct(his_only_check$Index))

cat("\nAfter excluding peptides explained by STY:",
    n_distinct(his_not_sty_check$Index))

cat("\nAfter reproducibility filter:",
    n_distinct(his_reproducible_check$Index))

cat("\nDuplicate sites:", nrow(duplicate_check))

# Show contaminant entries

cat("\n========== CONTAMINANTS ==========\n")

print(
  contaminant_check %>%
    select(Index, Protein, Peptide)
)

# Show duplicate entries

cat("\n========== DUPLICATES ==========\n")

print(duplicate_check)

# Show final candidates

cat("\n========== FINAL CANDIDATES ==========\n")

print(
  his_reproducible_check %>%
    select(
      Index,
      Protein,
      Peptide,
      best_loc_prob,
      no_alt_acceptor,
      seen_in_STY_search,
      max_reps_in_group))

# 9. Create filtering table

filter_data <- tibble(
  Step = factor(
    c(
      "All His sites\n(no contaminants)",
      "Localisation\nprobability >= 0.75",
      "No alternative\nS/T/Y acceptor",
      "Not explained\nby STY search",
      "Reproducible\nin >= 2/3 replicates"),
    levels = c(
      "All His sites\n(no contaminants)",
      "Localisation\nprobability >= 0.75",
      "No alternative\nS/T/Y acceptor",
      "Not explained\nby STY search",
      "Reproducible\nin >= 2/3 replicates")),
  Count = c(
    n_distinct(his_no_contam$Index),
    n_distinct(his_075_check$Index),
    n_distinct(his_only_check$Index),
    n_distinct(his_not_sty_check$Index),
    n_distinct(his_reproducible_check$Index)))

print(filter_data)

# 10. Final plot

ggplot(filter_data, aes(x = Step, y = Count)) +
  
  geom_col(
    fill = "steelblue",
    width = 0.65) +
  
  geom_text(
    aes(label = Count),
    vjust = -0.5,
    size = 7,
    fontface = "bold") +
  
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 6),
    expand = expansion(mult = c(0, 0.15))
  ) +
  
  labs(
    title = "Filtering potential histidine phosphorylation sites",
    subtitle = "Number of sites remaining after each confidence criterion",
    x = NULL,
    y = "Number of His-site candidates") +
  
  theme_minimal(base_size = 20) +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      colour = "grey85",
      linewidth = 0.5),
    
    plot.title = element_text(
      hjust = 0.5,
      size = 32,
      face = "bold"),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 21),
    
    axis.text.x = element_text(
      angle = 20,
      hjust = 1,
      size = 19,
      colour = "black"),
    
    axis.text.y = element_text(
      size = 19,
      colour = "black"),
    
    axis.title.y = element_text(
      size = 23,
      face = "bold"),
    
    plot.margin = margin(15, 20, 20, 15))
