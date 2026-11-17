
library(readr)
library(dplyr)
library(ggplot2)

# 1. Read both datasets
sty <- read_tsv("combined_site_STY.tsv",
                show_col_types = FALSE)

styh <- read_tsv("combined_site_STYH.tsv",
                 show_col_types = FALSE)

# 2. Exclude contaminants, FGCZ, reverse and decoys
remove_unwanted <- function(df) {
  annotation <- df %>%
    select(any_of(c("Gene", "Protein", "Protein ID")))
  
  unwanted <- apply(annotation, 1, function(x) {
    any(grepl("contaminant|fgcz|reverse|decoy",
              ifelse(is.na(x), "", as.character(x)),
              ignore.case = TRUE))
  })
  
  df[!unwanted, , drop = FALSE]
}

sty <- remove_unwanted(sty)
styh <- remove_unwanted(styh)

# 3. Match peptides between STY and STYH
sty_keys <- sty %>%
  transmute(
    protein_id = `Protein ID`,
    peptide_key = toupper(Peptide)
  ) %>%
  distinct()

h_candidates <- styh %>%
  filter(grepl("h", Peptide, fixed = TRUE)) %>%
  mutate(
    protein_id = `Protein ID`,
    peptide_key = toupper(Peptide)
  ) %>%
  left_join(
    sty_keys %>% mutate(in_STY = TRUE),
    by = c("protein_id", "peptide_key")
  ) %>%
  mutate(
    group = ifelse(
      !is.na(in_STY),
      "H candidate: also in STY",
      "H candidate: no STY match"
    )
  )

# 4. Prepare density plot data
plot_data <- bind_rows(
  sty %>%
    transmute(
      prob = as.numeric(`Best Localization Probability`),
      group = "STY reference"
    ),
  h_candidates %>%
    transmute(
      prob = as.numeric(`Best Localization Probability`),
      group
    )
) %>%
  filter(!is.na(prob), prob >= 0, prob <= 1) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "STY reference",
        "H candidate: also in STY",
        "H candidate: no STY match"
      )
    )
  )

# 5. Show counts in console only
print(
  plot_data %>%
    group_by(group) %>%
    summarise(
      n = n(),
      n_above_075 = sum(prob >= 0.75),
      percent_above_075 = round(
        100 * mean(prob >= 0.75), 1
      ),
      .groups = "drop"
    )
)

# 6. Create density plot
p <- ggplot(
  plot_data,
  aes(x = prob, fill = group, color = group)
) +
  geom_density(
    alpha = 0.34,
    linewidth = 1.05,
    adjust = 1.25,
    bounds = c(0, 1)
  ) +
  geom_vline(
    xintercept = 0.75,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 0.735,
    y = Inf,
    label = "LP = 0.75",
    hjust = 1.05,
    vjust = 1.5,
    size = 5
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25)
  ) +
  scale_fill_manual(
    values = c(
      "STY reference" = "#A5AAA8",
      "H candidate: also in STY" = "#C77D70",
      "H candidate: no STY match" = "#38866D"
    )
  ) +
  scale_color_manual(
    values = c(
      "STY reference" = "#747B78",
      "H candidate: also in STY" = "#AD6056",
      "H candidate: no STY match" = "#246F57"
    )
  ) +
  labs(
    title = "How confident are the phosphosite assignments?",
    subtitle = "JE2 phosphoproteomics | STY vs STYH searches",
    x = "Best localization probability",
    y = "Density",
    fill = NULL,
    color = NULL,
    caption = paste(
      "STY matching is based on protein ID and peptide sequence,",
      "not verified site localization.\n",
      "A high localization probability does not validate pHis."
    )
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(
      size = 22, face = "bold"
    ),
    plot.subtitle = element_text(
      size = 15, color = "grey35"
    ),
    axis.title = element_text(
      size = 19, face = "bold"
    ),
    axis.text = element_text(
      size = 16
    ),
    legend.position = "bottom",
    legend.text = element_text(
      size = 15
    ),
    plot.caption = element_text(
      size = 12, hjust = 0, color = "grey40"
    ),
    panel.grid.major = element_line(
      color = "grey92", linewidth = 0.3
    )
  )

# 7. Display and save
print(p)

ggsave(
  "JE2_STY_STYH_density_poster.png",
  p, width = 13, height = 8, dpi = 300
)

ggsave(
  "JE2_STY_STYH_density_poster.pdf",
  p, width = 13, height = 8
)
