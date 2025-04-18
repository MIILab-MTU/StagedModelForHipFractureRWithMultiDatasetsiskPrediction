library(ggplot2)
library(dplyr)
library(tidyr)

# Combined variable list from both groups
plot_data <- data.frame(
  Variable = c(
    "HA_FND", "HA_GRIPAVG", "HA_THD", "HA_WLKSPED", "HA_LSD", "HA_BMI", "HA_MMSE",
    "HA_AGE", "HA_CALCIUM", "HA_TRLBTS", "HA_HEIGHT", "HA_IADL51", "HA_SMOKE", "HA_KIDNYST",
    "HA_SLDFX", "HA_WRSTFX"
  ),
  MrOS_Male = c(1.00, 2.00, 3.00, 4.00, 5.26, 6.30, 6.62, 7.82, 9.38, 10.34, 10.64, 11.64, 13.02, 13.98, NA, NA),
  SOF_Female = c(1.00, 8.21, 2.00, 5.13, 3.45, 11.99, 6.37, 9.87, NA, 10.65, 9.28, 13.00, 14.00, NA, 3.56, 6.49)
)

# Convert to long format
plot_data_long <- plot_data %>%
  pivot_longer(cols = c("MrOS_Male", "SOF_Female"),
               names_to = "Group",
               values_to = "RankProportion")

# Reorder by decreasing importance (i.e., increasing RankProportion)
# Using the maximum RankProportion per variable for consistent sorting
max_ranks <- plot_data_long %>%
  group_by(Variable) %>%
  summarize(MaxRank = max(RankProportion, na.rm = TRUE))

plot_data_long <- left_join(plot_data_long, max_ranks, by = "Variable")

# Plot
ggplot(plot_data_long, aes(x = RankProportion, y = reorder(Variable, MaxRank), fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("MrOS_Male" = "steelblue", "SOF_Female" = "orange")) +
  labs(
    title = "Top 15 Variables for MrOS Male and SOF Female (Ensemble 2)",
    x = "Rank Proportion (Lower = More Important)",
    y = "Variable"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(), legend.position = "right")
