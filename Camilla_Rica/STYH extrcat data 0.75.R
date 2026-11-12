#Exclude proteins containing contam_sp or FGCZ.
library(dplyr)
library(readr)

# 1. Read original STYH dataset

                
library(here)
source(here("R", "setup.R"))
data <- read_tsv(
  here("Camilla_Rica", "combined_site_STYH.tsv"))

# 2. Select peptides containing lowercase h
# and exclude contaminant proteins
histidine <- data %>%
  filter(
    grepl("h", Peptide, fixed = TRUE),
    !grepl("contam_sp|FGCZ", Protein,
           ignore.case = TRUE)
  ) %>%
  select(
    Protein,
    `Protein ID`,
    Peptide,
    `Best Localization Probability`
  )

# Number of initial candidates
nrow(histidine)  # 624

# 3. Keep localization probability >= 0.75
histidine_075 <- histidine %>%
  filter(`Best Localization Probability` >= 0.75) %>%
  arrange(desc(`Best Localization Probability`))

# Number of remaining candidates
nrow(histidine_075)  # 308

# 4. Export as CSV
write_csv(histidine_075, "histidine_075.csv")


library(ggplot2)

ggplot(histidine, aes(x = `Best Localization Probability`)) +
  geom_histogram(binwidth = 0.05,
                 fill = "steelblue",
                 color = "white") +
  geom_vline(xintercept = 0.75,
             color = "red",
             linetype = "dashed",
             linewidth = 1) +
  labs(
    title = "Potential histidine phosphorylation sites",
    x = "Best Localization Probability",
    y = "Number of peptide entries"
  ) +
  theme_minimal() # IMPORTANT; here, still STY can be phosphorilated.
















