# Load required libraries
library(ggplot2)

# Read data from the provided CSV file
data <- read.csv("violin_input_file.csv")

# Reorder factor levels of 'Group' variable
data$Group <- factor(data$Group, levels = c("Ex-miRs", "All_miRs"))

# Create a violin plot
p <- ggplot(data, aes(x = Group, y = Values)) +
    geom_violin(fill = "orange", color = "darkred") +
    geom_boxplot(width = 0.1, fill = "white", color = "black") +  # Add boxplot for better visualization
    ylim(0,100) +       ## Set y-axis limit
    theme_minimal() +
    labs(x = "Group",
         y = "GC content (%)",
         fill = "Group") +
    theme(axis.title.x = element_text(size = 14, face = "bold"),  # Increase font size and make it bold for X-axis label
          axis.title.y = element_text(size = 14, face = "bold"),  # Increase font size and make it bold for Y-axis label
          axis.text.x = element_text(size = 12, face = "bold"),    # Increase font size and make it bold for X-axis text
          axis.text.y = element_text(size = 12, face = "bold"))    # Increase font size and make it bold for Y-axis text

# Perform one-sided Wilcoxon rank sum test with V1 as Ex-miRs and V2 as All_miRs, alternative = "less"
wilcox_test_result <- wilcox.test(data[data$Group == "Ex-miRs", "Values"],
                                  data[data$Group == "All_miRs", "Values"],
                                  alternative = "greater")

# Extract p-value
p_value <- wilcox_test_result$p.value

# Plot p-value on the graph
p <- p + annotate("text", x = 1.5, y = max(data$Values),
                  label = paste("p-value:", round(p_value, 4)),
                  size = 6, fontface = "bold")  # Increase font size and make it bold

# Save the plot as JPG with 600 dpi resolution
ggsave("plot.jpg", plot = p, width = 8, height = 6, units = "in", dpi = 600, type = "jpg")
