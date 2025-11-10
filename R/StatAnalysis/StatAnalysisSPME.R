# Script to analysis SMPE data using ANOVA and Tukey's Honest test
# for lab data from bioaugmentation experiments

# Install Packages
{
  install.packages("dplyr")
  install.packages("tidyr")
}

# Load Libraries
{
  library(dplyr)
  library(tidyr)
}

# Read data ---------------------------------------------------------------
obs.data <- read.csv("Data/uncoated_biochar_V2.csv")

# Format data -------------------------------------------------------------
# Select SPME & time series data
pcbi.spme <- obs.data %>%
  filter(Sample_medium =="SPME", Experiment =="biochar_timeseries")

pcbi.spme$Group_biochar <- paste(pcbi.spme$Group,
                                 pcbi.spme$percent_biochar, sep = "_")
pcbi.spme$Group_biochar <- factor(pcbi.spme$Group_biochar)

# ANOVA -------------------------------------------------------------------
pcb.cols <- grep("^PCB_", names(pcbi.spme), value = TRUE)
time_points <- unique(pcbi.spme$time)

all_pcb <- c(pcb.cols, "tPCB", "LCPCB")
anova_pvalues <- data.frame(
  Variable = character(),
  Time = numeric(),
  p_value = numeric()
)

for (time_pt in time_points) {
  df_time <- pcbi.spme %>% filter(time == time_pt)
  
  for (col in all_pcb) {
    fit <- aov(as.formula(paste(col, "~ Group_biochar")), data = df_time)
    pval <- summary(fit)[[1]][["Pr(>F)"]][1]
    
    anova_pvalues <- rbind(anova_pvalues, data.frame(
      Variable = col,
      Time = time_pt,
      p_value = pval
    ))
  }
}

anova_sig <- anova_pvalues %>%
  filter(p_value < 0.05)

# Check
print(anova_sig)

# Export data
write.csv(anova_sig, file = "Output/Data/Stat/ANOVA_SPME.csv")

# Tukey's test ------------------------------------------------------------
# Initialize empty data frame
tukey_sig_df <- data.frame()

# Loop through each significant ANOVA result
for (i in 1:nrow(anova_sig)) {
  var <- anova_sig$Variable[i]
  t <- anova_sig$Time[i]
  
  df_time <- pcbi.spme %>% filter(time == t)
  fit <- aov(as.formula(paste(var, "~ Group_biochar")), data = df_time)
  tuk <- TukeyHSD(fit)
  
  df <- as.data.frame(tuk$Group_biochar) %>%
    mutate(
      Variable = var,
      Time = t,
      Comparison = rownames(.)
    ) %>%
    rename(p.adj = `p adj`) %>%
    select(Variable, Time, Comparison, diff, lwr, upr, p.adj)
  
  # Split the comparison into two groups (e.g., "Control_0-Treatment_5")
  df <- df %>%
    tidyr::separate(Comparison, into = c("Group1", "Group2"), sep = "-") %>%
    mutate(
      Higher_Group = ifelse(diff > 0, Group1, Group2)  # which has higher mean
    )
  
  # Keep only significant comparisons
  df_sig <- df %>% filter(p.adj < 0.05)
  
  tukey_sig_df <- bind_rows(tukey_sig_df, df_sig)
}

# Arrange and inspect
tukey_sig_df <- tukey_sig_df %>%
  arrange(Variable, Time, p.adj)

print(tukey_sig_df)

# Export results
write.csv(tukey_sig_df, file = "Output/Data/Stat/Tukey_SPME.csv",
          row.names = FALSE)
