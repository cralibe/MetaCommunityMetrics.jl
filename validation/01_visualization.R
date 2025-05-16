library(tidyverse)
library(ggplot2)
library(cowplot)
#Load in the dnci values
dnci_df<-read.csv("validation_output/dnci_combined_df.csv")

dnci_plot<-dnci_df %>%
  unite(group_pair, c("group1", "group2"), sep = "-") %>%
  ggplot(aes(x = group_pair, y = DNCI, color = programming_language, fill = programming_language)) +
  geom_violin(position = position_dodge(width = 0.9), trim = FALSE, alpha = 0.4) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.7) +
  labs(x = "Group Pair", y = "DNCI", color = "Programming Language", fill = "Programming Language") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")+
  ylim(-3.5, 0)+
  theme_cowplot()

ggsave("validation_output/dnci_plot.png", dnci_plot, dpi=300, width = 8, height = 5, bg ="white")
  