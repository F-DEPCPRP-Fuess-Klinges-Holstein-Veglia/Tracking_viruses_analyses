# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the data
data <- read.csv("rnavsrnaviralspades_comp.csv")

# Reshape the data to include the Total column
data_long <- data %>%
  pivot_longer(cols = c("Total", "Duplodnaviria", "Monodnaviria", "Riboviria", "Varidnaviria"),
               names_to = "Viral_Group",
               values_to = "Sequence_Count")

# Reorder the Viral_Group factor so Total appears first
data_long$Viral_Group <- factor(data_long$Viral_Group, 
                                levels = c("Total", "Duplodnaviria", "Monodnaviria", "Riboviria", "Varidnaviria"))

# Perform pairwise t-tests for each viral group
stat_results <- data_long %>%
  group_by(Viral_Group) %>%
  summarise(
    t_test_p_value = t.test(Sequence_Count ~ algorithm, data = .)$p.value
  )

# Adjust p-values for multiple testing (e.g., Bonferroni correction)
stat_results <- stat_results %>%
  mutate(adjusted_p_value = p.adjust(t_test_p_value, method = "bonferroni"))

# Print the statistical results
print(stat_results)

# Create the facet-wrapped box plot
ggplot(data_long, aes(x = algorithm, y = Sequence_Count, fill = algorithm)) +
  geom_boxplot() +
  facet_wrap(~ Viral_Group, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Comparison of Virus Sequences by Algorithm and Viral Group",
    x = "Assembly Algorithm",
    y = "Sequence Count",
    fill = "Algorithm"
  ) +
  theme(
    strip.text = element_text(size = 14),  # Adjust facet label font size
    axis.text.x = element_text(angle = 45, hjust = 1)  # Adjust x-axis text
  )
