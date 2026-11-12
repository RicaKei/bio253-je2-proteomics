
library(dplyr)
library(readr)
library(ggplot2)

# 1. Read your original STYH dataset
data <- read_tsv("C:\\Users\\Rica\\Documents\\UZH\\Bio253\\bio253-je2-proteomics\\Camilla_Rica\\combined_site_STYH.tsv")

# 2. Find all lowercase h and exclude contaminants
his_all <- data %>%
  filter(
    grepl("h", Peptide, fixed = TRUE),
    !grepl("contam_sp|FGCZ", Protein,
           ignore.case = TRUE)
  )

# 3. Select probability >= 0.75
his_075 <- his_all %>%
  filter(`Best Localization Probability` >= 0.75)

# 4. Exclude peptides containing S/T/Y
his_only <- his_075 %>%
  filter(!grepl("[STYsty]", Peptide))

# Check how many entries remain
nrow(his_only)

# 5. Create dot plot
ggplot(his_only,
       aes(x = `Best Localization Probability`,
           y = reorder(Peptide,
                       `Best Localization Probability`))) +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(
    xintercept = 0.75,
    linetype = "dashed",
    color = "red"
  ) +
  scale_x_continuous(
    limits = c(0.7, 1.02),
    breaks = seq(0.75, 1, 0.05)
  ) +
  labs(
    title = "Potential histidine phosphorylation sites",
    subtitle = "Peptides without S/T/Y residues",
    x = "Best Localization Probability",
    y = "Peptide sequence"
  ) +
  theme_minimal()

