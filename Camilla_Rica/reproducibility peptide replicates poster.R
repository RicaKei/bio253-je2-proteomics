
# Poster heatmap with larger font sizes

library(ggplot2)
library(dplyr)
library(tidyr)

peptides <- c(
  "Peptide 1",
  "Peptide 2",
  "Peptide 3",
  "Peptide 4",
  "Peptide 5"
)

conditions <- c(
  "Anc TSB", "Anc SN", "Anc pulse",
  "stk1 TSB", "stk1 SN", "stk1 pulse",
  "stp TSB", "stp SN", "stp pulse"
)

# Confirmed replicate counts
counts <- matrix(NA_integer_, nrow = 5, ncol = 9)

counts[1, c(1, 2, 5, 8)] <- c(2, 3, 3, 2)
counts[2, 5] <- 2
counts[3, 7] <- 2
counts[4, 7] <- 2
counts[5, 2] <- 2

rownames(counts) <- peptides
colnames(counts) <- conditions

# Prepare data
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

# Create heatmap
p <- ggplot(
  plot_data,
  aes(x = Condition, y = Peptide, fill = Detection)
) +
  geom_tile(
    color = "white",
    linewidth = 1.5
  ) +
  
  # Larger replicate numbers inside tiles
  geom_text(
    aes(label = ifelse(
      is.na(Replicates), "",
      paste0(Replicates, "/3")
    )),
    size = 9,
    fontface = "bold",
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "3/3" = "#167b65",
      "2/3" = "#78b7a4",
      "Not confirmed" = "grey90"
    ),
    breaks = c("3/3", "2/3", "Not confirmed")
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    fill = "Replicates"
  ) +
  
  theme_minimal(base_size = 22) +
  
  theme(
    # X-axis: experimental conditions
    axis.text.x = element_text(
      size = 24,
      face = "bold",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    
    # Y-axis: peptide names
    axis.text.y = element_text(
      size = 26,
      face = "bold"
    ),
    
    # Legend title
    legend.title = element_text(
      size = 25,
      face = "bold"
    ),
    
    # Legend labels
    legend.text = element_text(
      size = 23
    ),
    
    # Larger legend colour boxes
    legend.key.size = grid::unit(1.2, "cm"),
    
    # Remove background grid
    panel.grid = element_blank(),
    
    # Margins
    plot.margin = margin(25, 35, 25, 25)
  )

# Display heatmap
print(p)

# Save as high-resolution PNG
ggsave(
  "pHis_heatmap_poster.png",
  plot = p,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300
)

# Save as vector PDF
ggsave(
  "pHis_heatmap_poster.pdf",
  plot = p,
  width = 16,
  height = 9,
  units = "in"
)
