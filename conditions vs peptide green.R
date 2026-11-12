#plot to show contitions and peptide name

library(ggplot2)
library(dplyr)
library(tidyr)

peptides <- c(
  "ANhEQIIEAIK;MANhEQIIEAIK",
  "EVLKDHLhK",
  "DLLAVhLK",
  "LDAFNGPLDLLLhLIQK",
  "HLVEQALNQAIDQINhADK"
)

conditions <- c(
  "Anc TSB", "Anc SN", "Anc pulse",
  "stk1 TSB", "stk1 SN", "stk1 pulse",
  "stp TSB", "stp SN", "stp pulse"
)

# Only previously confirmed counts
counts <- matrix(NA_integer_, nrow = 5, ncol = 9)

counts[1, c(1, 2, 5, 8)] <- c(2, 3, 3, 2)
counts[2, 5] <- 2
counts[3, 7] <- 2
counts[4, 7] <- 2
counts[5, 2] <- 2

rownames(counts) <- peptides
colnames(counts) <- conditions

plot_data <- as.data.frame(as.table(counts))
names(plot_data) <- c("Peptide", "Condition", "Replicates")

plot_data <- plot_data %>%
  mutate(
    Condition = factor(Condition, levels = conditions),
    Peptide = factor(Peptide, levels = rev(peptides)),
    Detection = case_when(
      Replicates == 3 ~ "3/3",
      Replicates == 2 ~ "2/3",
      TRUE ~ "Not confirmed"
    )
  )

ggplot(plot_data,
       aes(x = Condition, y = Peptide, fill = Detection)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = ifelse(is.na(Replicates), "",
                       paste0(Replicates, "/3"))),
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      "3/3" = "#167b65",
      "2/3" = "#78b7a4",
      "Not confirmed" = "grey90"
    )
  ) +
  labs(
    title = "Detection of potential pHis peptides",
    subtitle = "S. aureus JE2 ancestor, stk1 and stp mutants",
    x = "Experimental condition",
    y = "Peptide sequence",
    fill = "Detection"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
