# Count cCREs containing >=1 base of overlap with HAR-fixed, HAQER-fixed,
# or Human-FLSS annotations and render a barplot.


suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# Input: per-cCRE table produced by ccre_har_haqer_fix.py
ccre <- read.csv("ccre_har_haqer_fix.csv")
total <- nrow(ccre)

# Count cCREs with >=1 bp of overlap (or >=1 fHSS variant) per annotation class
summary_df <- data.frame(
  Category = c("HAR-fixed", "HAQER-fixed", "Human-FLSS"),
  Count    = c(sum(ccre$har_count   > 0),
               sum(ccre$haqer_count > 0),
               sum(ccre$fhss_count  > 0)),
  stringsAsFactors = FALSE
)
summary_df$Percent <- round(summary_df$Count / total * 100, 2)

# Persist the summary CSV (same schema as the prior Python output)
write.csv(summary_df, "ccre_har_haqer_fix_overlap.csv", row.names = FALSE)

cat(sprintf("Total cCREs: %s\n", comma(total)))
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  %s: %s (%s%%)\n",
              summary_df$Category[i],
              comma(summary_df$Count[i]),
              summary_df$Percent[i]))
}

# ── Plot ─────────────────────────────────────────────────────────────────
summary_df$Category <- factor(summary_df$Category,
                              levels = c("HAR-fixed", "HAQER-fixed", "Human-FLSS"))

colors <- c("HAR-fixed"   = "#b2182b",
            "HAQER-fixed" = "#e08214",
            "Human-FLSS"  = "#2166ac")

p <- ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = colors) +
  scale_y_log10(labels = comma) +
  geom_text(aes(label = paste0(comma(Count), "\n(", Percent, "%)")),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  labs(x = NULL, y = "Number of cCREs", title = "cCREs containing overlaps") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12)
  ) +
  coord_cartesian(ylim = c(500, 5e6))

ggsave("ccre_har_haqer_fix_overlap.pdf", p, width = 5, height = 5)
ggsave("ccre_har_haqer_fix_overlap.png", p, width = 5, height = 5, dpi = 300)

cat("\nWrote ccre_har_haqer_fix_overlap.csv, .pdf, .png\n")
