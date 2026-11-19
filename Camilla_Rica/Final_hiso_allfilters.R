#FINAL plot 
library(dplyr)
library(readr)
library(ggplot2)


library(here)
source(here("R", "setup.R"))
data <- read_tsv(
  here("Camilla_Rica", "combined_site_STYH.tsv"))


library(ggplot2)

# Create table 
filter_data <- data.frame(
  Step = factor(
    c("Histidine phosphorylation",
      "LP >= 0.75",
      "No S/T/Y on peptide",
      "Detected in >=2/3 replicates"),
    levels = c(
      "Histidine phosphorylation",
      "LP >= 0.75",
      "No S/T/Y on peptide",
      "Detected in >=2/3 replicates"
    )
  ),
  Count = c(
    nrow(his_all),
    nrow(his_075),
    nrow(his_only),
    5
  )
)

# 
ggplot(filter_data, aes(x = Step, y = Count)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_text(aes(label = Count),
            vjust = -0.5, size = 5) +
  labs(
    title = "Filtering potential histidine phosphorylation sites",
    x = NULL,
    y = "Number of peptide candidates"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )



#size good for the poster 

ggplot(filter_data, aes(x = Step, y = Count)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_text(aes(label = Count),
            vjust = -0.5, size = 5) +
  labs(
    title = "Filtering potential histidine phosphorylation sites",
    x = NULL,
    y = "Number of peptide candidates"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5, size = 18, face = "bold"
    ),
    axis.text.x = element_text(
      angle = 25, hjust = 1, size = 14
    ),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16)
  )








# FINAL diagnostic plot: Histidine phosphorylation filtering

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(here)

# Load data
data <- read_tsv(
  here("Camilla_Rica", "combined_site_STYH.tsv")
)

# --------------------------------------------------
# 1. Identify all histidine phosphorylation sites
# --------------------------------------------------

his_all_check <- data %>%
  filter(str_detect(Index, "_H\\d+$")) %>%
  distinct(Index, .keep_all = TRUE)

# --------------------------------------------------
# 2. Check contaminants
# --------------------------------------------------

contaminant_check <- his_all_check %>%
  filter(
    str_detect(Protein, regex("contam_", ignore_case = TRUE))
  )

his_no_contam <- his_all_check %>%
  filter(
    !str_detect(Protein, regex("contam_", ignore_case = TRUE))
  )

# --------------------------------------------------
# 3. Localization probability >= 0.75
# --------------------------------------------------

his_075_check <- his_no_contam %>%
  filter(`Best Localization Probability` >= 0.75)

# --------------------------------------------------
# 4. Remove peptides containing phosphorylated S/T/Y
# Lowercase s, t, y indicate modified residues
# --------------------------------------------------

his_only_check <- his_075_check %>%
  filter(!str_detect(Peptide, "[sty]"))

# --------------------------------------------------
# 5. Check duplicates
# --------------------------------------------------

duplicate_check <- his_only_check %>%
  count(Index) %>%
  filter(n > 1)

# --------------------------------------------------
# 6. Compare with previous objects
# --------------------------------------------------

comparison <- tibble(
  Dataset = c(
    "Original his_all",
    "Original his_075",
    "Original his_only"
  ),
  Rows = c(
    nrow(his_all),
    nrow(his_075),
    nrow(his_only)
  ),
  Unique_sites = c(
    n_distinct(his_all$Index),
    n_distinct(his_075$Index),
    n_distinct(his_only$Index)
  )
)

# --------------------------------------------------
# 7. Print diagnostic results
# --------------------------------------------------

cat("\n========== HISTIDINE SITE ANALYSIS ==========\n")

cat("\nAll histidine sites:", nrow(his_all_check))
cat("\nUnique histidine sites:", n_distinct(his_all_check$Index))

cat("\n\nContaminants:", nrow(contaminant_check))
cat("\nSites after contaminant removal:", nrow(his_no_contam))

cat("\n\nLP >= 0.75:", nrow(his_075_check))

cat("\n\nNo phosphorylated S/T/Y on peptide:",
    nrow(his_only_check))

cat("\n\nDuplicate sites:", nrow(duplicate_check))

cat("\n\n========== COMPARISON WITH ORIGINAL OBJECTS ==========\n")

print(comparison)

# Show contaminant entries
cat("\n========== CONTAMINANTS ==========\n")
print(
  contaminant_check %>%
    select(Index, Protein, Peptide)
)

# Show duplicate entries
cat("\n========== DUPLICATES ==========\n")
print(duplicate_check)

# --------------------------------------------------
# 8. Create filtering table
# --------------------------------------------------
# --------------------------------------------------
# 8. Create filtering table
# --------------------------------------------------

filter_data <- tibble(
  Step = factor(
    c(
      "Histidine phosphorylation",
      "LP >= 0.75",
      "No S/T/Y on peptide",
      "Detected in >=2/3 replicates"
    ),
    levels = c(
      "Histidine phosphorylation",
      "LP >= 0.75",
      "No S/T/Y on peptide",
      "Detected in >=2/3 replicates"
    )
  ),
  Count = c(
    n_distinct(his_no_contam$Index),
    n_distinct(his_075_check$Index),
    n_distinct(his_only_check$Index),
    5
  )
)

# --------------------------------------------------
# 9. Final plot
# --------------------------------------------------

ggplot(filter_data, aes(x = Step, y = Count)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_text(
    aes(label = Count),
    vjust = -0.5,
    size = 5
  ) +
  labs(
    title = "Filtering potential histidine phosphorylation sites",
    x = NULL,
    y = "Number of peptide candidates"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      size = 14,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )
