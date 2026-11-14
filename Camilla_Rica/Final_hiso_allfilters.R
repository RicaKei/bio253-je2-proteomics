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
