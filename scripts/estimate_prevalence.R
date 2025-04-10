#! /usr/bin/env Rscript
library(PooledInfRate)
library(ggplot2)

virus_data <- read.delim("pathogen_presence.tsv", header=TRUE, sep="\t")
virus_columns <- colnames(virus_data)[3:ncol(virus_data)]

#===============================================
# Estimate prevalence for each virus
#===============================================

results <- lapply(virus_columns, function(virus) {
  x <- as.numeric(virus_data[[virus]])
  m <- as.numeric(virus_data$pool_size)

  if (sum(x) == 0) return(data.frame(Virus = virus, Prevalence = 0, Lower_CI = 0, Upper_CI = 0))

  result <- pooledBin(x, m)
  data.frame(Virus = virus, Prevalence = result$P, Lower_CI = result$Lower, Upper_CI = result$Upper)
})

prevalence_df <- do.call(rbind, results)
write.table(prevalence_df, "virus_prevalence_estimates.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#===============================================
# Plot prevalence estimates
#===============================================

p <- ggplot(prevalence_df) +
 geom_point(aes(x = Prevalence, y = Virus), size = 3) +
  geom_point(aes(x = Lower_CI, y = Virus), size = 1, color = "darkgray") +
  geom_point(aes(x = Upper_CI, y = Virus), size = 1, color = "darkgray") +
  geom_segment(aes(x = Lower_CI, xend = Upper_CI, y = Virus, yend = Virus),
              color = "darkgray", alpha = 0.5) +
  scale_x_log10(labels = scales::percent_format(accuracy = 0.1)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Virus Prevalence Estimates with Confidence Intervals", x = "Prevalence (%)", y = NULL)

ggsave("prevalence_dotplot.pdf", p, width = 10, height = 6)
