# Statistical analysis to review data, mostly between control
# and experiments, same time points.

# Packages and libraries --------------------------------------------------
# Install packages
install.packages("dplyr")
install.packages("ggplot2")

# Load libraries
{
  library(dplyr) # organize data
  library(ggplot2) # plotting
  library(scales)
}

# Read data ---------------------------------------------------------------
{
  obs.data <- read.csv("Data/RTqPCR.csv")
  # Extract relevant columns
  cDNA <- obs.data[, c("Experiment", "Percent_biochar",
                       "Group", "time", "Replicate", "cDNA")]
  DNA <- obs.data[, c("Experiment", "Percent_biochar",
                      "Group", "time", "Replicate", "DNA")]
  tgr <- obs.data[, c("Experiment", "Percent_biochar",
               "Group", "time", "Replicate", "transcript_gene_ratio")]
}

# Select
cDNA$Group_biochar <- paste(cDNA$Group, cDNA$Percent_biochar, sep = "-")
DNA$Group_biochar <- paste(DNA$Group, cDNA$Percent_biochar, sep = "-")
tgr$Group_biochar <- paste(tgr$Group, cDNA$Percent_biochar, sep = "-")

# Plot
plot.cdna <- ggplot(cDNA, aes(x = time, y = cDNA, color = Group_biochar)) +
  geom_point(size = 2.5, alpha = 0.8,
             position = position_jitter(width = 0.5)) +
  theme_bw() +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(10^2, 10^8)) +
  labs(x = expression(bold("Time (day)")),
       y = bquote(bold("cDNA (copies/g)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

# See plot
plot.cdna

plot.dna <- ggplot(DNA, aes(x = time, y = DNA, color = Group_biochar)) +
  geom_point(size = 2.5, alpha = 0.8,
             position = position_jitter(width = 0.5)) +
  theme_bw() +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(1, 10^8)) +
  labs(x = expression(bold("Time (day)")),
       y = bquote(bold("DNA (copies/g)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

# See plot
plot.dna

plot.tgr <- ggplot(tgr, aes(x = time, y = transcript_gene_ratio, color = Group_biochar)) +
  geom_point(size = 2.5, alpha = 0.8,
             position = position_jitter(width = 0.5)) +
  theme_bw() +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(1, 10^6)) +
  labs(x = expression(bold("Time (day)")),
       y = bquote(bold("Transcript Gene Ratio (copies/g)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

# See plot
plot.tgr

