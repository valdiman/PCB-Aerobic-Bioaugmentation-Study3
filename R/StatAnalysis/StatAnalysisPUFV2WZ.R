# Script to analysis PUF data using ANOVA and Tukey's Honest test
# for lab data from bioaugmentation experiments

# Install Packages
{
  install.packages("dplyr")
  install.packages("tidyr")
  install.packages("tibble")
}

# Load Libraries
{
  library(dplyr)
  library(tidyr)
  library(tibble)
}

# Read data ---------------------------------------------------------------
obs.data <- read.csv("Data/uncoated_biochar_V2.csv")

# Format data -------------------------------------------------------------
# Select PUF & time series data
pcbi.puf <- obs.data %>%
  filter(Sample_medium =="PUF", Experiment =="biochar_timeseries")

pcbi.puf$Group_biochar <- paste(pcbi.puf$Group,
                                 pcbi.puf$percent_biochar, sep = "_")
pcbi.puf$Group_biochar <- factor(pcbi.puf$Group_biochar)

# ANOVA -------------------------------------------------------------------
pcb.cols <- grep("^PCB_", names(pcbi.puf), value = TRUE)
time_points <- unique(pcbi.puf$time)

all_pcb <- c(pcb.cols, "tPCB", "LCPCB")
anova_pvalues <- data.frame(
  Variable = character(),
  Time = numeric(),
  p_value = numeric()
)

for (time_pt in time_points) {
  df_time <- pcbi.puf %>% filter(time == time_pt)
  
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
write.csv(anova_sig, file = "Output/Data/Stat/ANOVA_PUF.csv")

# Tukey's test ------------------------------------------------------------
# Initialize empty data frame
tukey_sig_df <- data.frame()

# Loop through each significant ANOVA result
for (i in seq_len(nrow(anova_sig))) {
  var <- anova_sig$Variable[i]
  t <- anova_sig$Time[i]
  
  df_time <- pcbi.puf %>% filter(time == t)
  fit <- aov(as.formula(paste(var, "~ Group_biochar")), data = df_time)
  tuk <- TukeyHSD(fit)
  
  # Tukey results to df; 'Comparison' from rownames
  df <- as.data.frame(tuk$Group_biochar) %>%
    rownames_to_column(var = "Comparison") %>%
    mutate(
      Variable = var,
      Time = t
    ) %>%
    rename(p.adj = `p adj`) %>%
    select(Variable, Time, Comparison, diff, lwr, upr, p.adj)
  
  # Split the comparison into two groups (e.g., "Control_0-Treatment_5")
  df <- df %>%
    tidyr::separate(Comparison, into = c("Group1", "Group2"), sep = "-") %>%
    mutate(Higher_Group = ifelse(diff > 0, Group1, Group2))
  
  # Compute group means for this time and variable
  # Use .data[[var]] to programmatically access column by name
  means_df <- df_time %>%
    group_by(Group_biochar) %>%
    summarize(mean_val = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  
  # Join means onto df by matching Group1 and Group2
  df <- df %>%
    left_join(means_df %>% rename(Group1 = Group_biochar, Mean_Group1 = mean_val),
              by = "Group1") %>%
    left_join(means_df %>% rename(Group2 = Group_biochar, Mean_Group2 = mean_val),
              by = "Group2")
  
  # Keep only significant comparisons
  df_sig <- df %>% filter(p.adj < 0.05)
  
  # Append to output
  tukey_sig_df <- bind_rows(tukey_sig_df, df_sig)
}

# Arrange and inspect
tukey_sig_df <- tukey_sig_df %>%
  arrange(Variable, Time, p.adj)

print(tukey_sig_df)

# Export results
write.csv(tukey_sig_df, file = "Output/Data/Stat/Tukey_PUF.csv",
          row.names = FALSE)
