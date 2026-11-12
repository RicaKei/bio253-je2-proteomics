#compare two data sets 
#Peptide contains lowercase h.
#Peptide contains no S, T, or Y, uppercase or lowercase.
#Best Localization Probability ≥ 0.75.
#Exclude proteins containing contam_sp or FGCZ.
#Check whether the same unmodified peptide sequence appears in the STY-only dataset.

library(dplyr)
library(readr)

# 1. Read both datasets
STYH <- read_tsv("C:\\Users\\Rica\\Documents\\UZH\\Bio253\\bio253-je2-proteomics\\Camilla_Rica\\combined_site_STYH.tsv")
STY <- read_tsv("C:\\Users\\Rica\\Documents\\UZH\\Bio253\\bio253-je2-proteomics\\Camilla_Rica\\combined_site_STY.tsv")

# 2. Filter histidine phosphorylation candidates
histidine <- STYH %>%
  filter(
    grepl("h", Peptide, fixed = TRUE),
    !grepl("[STYsty]", Peptide),
    !grepl("contam_sp|FGCZ", Protein,
           ignore.case = TRUE),
    `Best Localization Probability` >= 0.75
  ) %>%
  select(
    Protein,
    `Protein ID`,
    Peptide,
    `Best Localization Probability`
  ) %>%
  arrange(desc(`Best Localization Probability`))

# 3. Count remaining entries
nrow(histidine)

# 4. Remove phosphorylation notation
# (convert peptide sequences to uppercase)
histidine <- histidine %>%
  mutate(sequence = toupper(Peptide))

STY <- STY %>%
  mutate(sequence = toupper(Peptide))

# 5. Compare sequences with STY dataset
histidine <- histidine %>%
  mutate(
    found_in_STY = sequence %in% STY$sequence
  )

# 6. Keep only peptides NOT found in STY
final <- histidine %>%
  filter(found_in_STY == FALSE)

# 7. View results
View(final)
nrow(final)

# 8. Save results
write_csv(final, "pHis_no_STY_075.csv")
