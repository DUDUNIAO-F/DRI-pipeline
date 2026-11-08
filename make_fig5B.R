# make_fig5B.R
# Redraw Figure 5B from Table_S14_spatial_null.csv (single source of truth)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

d <- read.csv("output/Table_S14_spatial_null.csv")
d$group <- sub("-.*", "", sub(".*_", "", d$Sample))
d$group <- factor(d$group, levels = c("N", "PTC", "LPTC", "ATC"))
d$sig <- ifelse(d$p_empirical == 0, "p < 0.001", sprintf("p = %.3f", d$p_empirical))

grp <- d %>% group_by(group) %>%
  summarise(n = n(), mean_rho = round(mean(rho_partial), 3),
            median_rho = round(median(rho_partial), 3), .groups = "drop")
cat("Group summary:\n"); print(as.data.frame(grp))
cat("\nOverall: mean =", round(mean(d$rho_partial), 3),
    "| median =", round(median(d$rho_partial), 3), "\n")

p <- ggplot(d, aes(x = group, y = rho_partial, fill = group)) +
  geom_hline(yintercept = -0.581, linetype = "dashed", color = "firebrick") +
  annotate("text", x = -Inf, y = -0.581, label = "bulk TCGA (rho = -0.581)",
           hjust = -0.1, vjust = -0.5, color = "firebrick", size = 3.2) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(aes(shape = sig), width = 0.08, size = 2.6, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 4, size = 3.2, stroke = 1.3) +
  scale_fill_manual(values = c(N = "grey70", PTC = "steelblue",
                               LPTC = "orange", ATC = "tomato")) +
  scale_shape_manual(values = c("p < 0.001" = 16, "p = 0.128" = 17)) +
  labs(x = NULL, y = "Partial Spearman rho (DRI vs TDS)",
       title = "Spatial DRI-TDS correlation across 16 Visium sections",
       subtitle = "All sections shown; crosses = group means",
       shape = "Permutation test") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("output/Fig5B_spatial_rho.pdf", p, width = 6.5, height = 4.8)
cat("\nSaved: output/Fig5B_spatial_rho.pdf\n")
