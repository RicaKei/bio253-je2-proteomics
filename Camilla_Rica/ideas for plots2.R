
library(dplyr)
library(readr)
library(ggplot2)

# Read original dataset
data <- read_tsv("C:\\Users\\Rica\\Documents\\UZH\\Bio253\\bio253-je2-proteomics\\Camilla_Rica\\combined_site_STYH.tsv")

# All lowercase h, excluding contaminants
his_all <- data %>%
  filter(
    grepl("h", Peptide, fixed = TRUE),
    !grepl("contam_sp|FGCZ", Protein,
           ignore.case = TRUE)
  )

# Probability >= 0.75
his_075 <- his_all %>%
  filter(`Best Localization Probability` >= 0.75)

# Histidine without any S/T/Y
his_only <- his_075 %>%
  filter(!grepl("[STYsty]", Peptide))

# Plot 1: Below vs above threshold
plot1_data <- his_all %>%
  mutate(Group = if_else(
    `Best Localization Probability` < 0.75,
    "Below 0.75", "0.75 or higher"
  )) %>%
  count(Group)

ggplot(plot1_data, aes(x = Group, y = n, fill = Group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5) +
  labs(title = "Localization probability of pHis candidates",
       x = NULL, y = "Peptide entries") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 2: h with S/T/Y vs h without S/T/Y
plot2_data <- his_075 %>%
  mutate(Group = if_else(
    grepl("[STYsty]", Peptide),
    "h + S/T/Y", "h without S/T/Y"
  )) %>%
  count(Group)

ggplot(plot2_data, aes(x = Group, y = n, fill = Group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5) +
  labs(title = "High-localization pHis candidates",
       x = NULL, y = "Peptide entries") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 3: Filtering funnel
plot3_data <- data.frame(
  Step = factor(
    c("All lowercase h", "Probability >= 0.75",
      "No S/T/Y"),
    levels = c("All lowercase h",
               "Probability >= 0.75", "No S/T/Y")
  ),
  Count = c(nrow(his_all),
            nrow(his_075),
            nrow(his_only))
)

ggplot(plot3_data, aes(x = Step, y = Count)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5) +
  labs(title = "Filtering potential pHis sites",
       x = NULL, y = "Peptide entries") +
  theme_minimal()

#
# Continue after creating his_only

library(ggplot2)

# Create table with filtering results
filter_data <- data.frame(
  Step = factor(
    c("Lowercase h",
      "Probability >= 0.75",
      "No S/T/Y",
      "Detected in >=2/3 replicates"),
    levels = c(
      "Lowercase h",
      "Probability >= 0.75",
      "No S/T/Y",
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

# Plot filtering progress
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

#
# Final 5 peptides identified in your analysis
final_peptides <- c(
  "ANhEQIIEAIK;MANhEQIIEAIK",
  "EVLKDHLhK",
  "DLLAVhLK",
  "LDAFNGPLDLLLhLIQK",
  "HLVEQALNQAIDQINhADK"
)

# Filter the 18 candidates down to the final 5
his_final <- his_only %>%
  filter(Peptide %in% final_peptides)

# Plot the final 5
ggplot(his_final,
       aes(x = `Best Localization Probability`,
           y = reorder(Peptide,
                       `Best Localization Probability`))) +
  geom_point(size = 3, color = "steelblue") +
  geom_vline(
    xintercept = 0.75,
    color = "red",
    linetype = "dashed"
  ) +
  scale_x_continuous(
    limits = c(0.70, 1.03),
    breaks = seq(0.75, 1, 0.05)
  ) +
  labs(
    title = "Final potential pHis candidates",
    subtitle = "Detected in at least 2 of 3 replicates",
    x = "Best Localization Probability",
    y = "Peptide sequence"
  ) +
  theme_minimal(base_size = 14)
