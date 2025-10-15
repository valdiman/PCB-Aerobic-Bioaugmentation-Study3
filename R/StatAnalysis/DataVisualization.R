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
}

# Read data ---------------------------------------------------------------
{
  obs.data <- read.csv("Data/uncoated_biochar_V2.csv")
  # Select individual congener from datasets
  pcb.ind <- "LCPCB"
  yvar <- sym(pcb.ind) 
  # Extract relevant columns
  pcbi <- obs.data[, c("Sample_medium", "Experiment", "percent_biochar",
                       "Group", "time", "Replicate", pcb.ind)]
}

# Select SPME & time series data
pcbi.spme <- pcbi %>%
  filter(Sample_medium =="SPME", Experiment =="biochar_timeseries")

pcbi.spme$Group_biochar <- paste(pcbi.spme$Group,
                                 pcbi.spme$percent_biochar, sep = "_")

# Plot
plot.spme <- ggplot(pcbi.spme, aes(x = time, y = !!yvar, color = Group_biochar)) +
  geom_point(size = 2.5, alpha = 0.8,
             position = position_jitter(width = 0.5)) +
  theme_bw() +
  labs(x = expression(bold("Time (day)")),
       y = bquote(bold(.(pcb.ind) ~ "(ng/cm)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

# See plot
print(plot.spme)

# Save plot in folder
ggsave(paste0("Output/Plots/TimeSeries/SPME/SPME_", pcb.ind, ".png"),
       plot = plot.spme, width = 10, height = 5, dpi = 500)

# Select PUF & time series data
pcbi.puf <- pcbi %>%
  filter(Sample_medium =="PUF", Experiment =="biochar_timeseries")

pcbi.puf$Group_biochar <- paste(pcbi.puf$Group,
                                 pcbi.puf$percent_biochar, sep = "_")

# Plot
plot.puf <- ggplot(pcbi.puf, aes(x = time, y = !!yvar, color = Group_biochar)) +
  geom_point(size = 2.5, alpha = 0.8,
             position = position_jitter(width = 0.5)) +
  theme_bw() +
  labs(x = expression(bold("Time (day)")),
       y = bquote(bold(.(pcb.ind) ~ "(ng/PUF)"))) +
  theme(axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.x = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10))

# See plot
print(plot.puf)

# Save plot in folder
ggsave(paste0("Output/Plots/TimeSeries/PUF/PUF_", pcb.ind, ".png"),
       plot = plot.puf, width = 10, height = 5, dpi = 500)
