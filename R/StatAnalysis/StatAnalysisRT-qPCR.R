# Script to analysis PUF data using ANOVA and Tukey's Honest test
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
  library(stringr)
}

# Read data ---------------------------------------------------------------
obs.data <- read.csv("Data/RTqPCR.csv")

# cDNA --------------------------------------------------------------------
cDNA <- obs.data[, c("Experiment", "Percent_biochar",
                     "Group", "time", "Replicate", "cDNA")]
cDNA$Group_biochar <- paste(cDNA$Group, cDNA$Percent_biochar, sep = "_")
# Remove time 0 Treatment_5
cDNA <- cDNA %>%
  filter(!(time == 0 & Group_biochar == "Treatment_5"))

# cDNA ANOVA --------------------------------------------------------------
time_points <- unique(cDNA$time)

anova_pvalues <- data.frame(
  Time = numeric(),
  p_value = numeric()
)

for (time_pt in time_points) {
  df_time <- cDNA %>% filter(time == time_pt)
  
  # Check number of unique groups
  if (length(unique(df_time$Group_biochar)) > 1) {
    fit <- aov(cDNA ~ Group_biochar, data = df_time)
    pval <- summary(fit)[[1]][["Pr(>F)"]][1]
    
    anova_pvalues <- rbind(anova_pvalues, data.frame(
      Time = time_pt,
      p_value = pval
    ))
  } else {
    # Skip this time point or record NA
    anova_pvalues <- rbind(anova_pvalues, data.frame(
      Time = time_pt,
      p_value = NA
    ))
  }
}

anova_sig <- anova_pvalues %>% filter(!is.na(p_value) & p_value < 0.05)
anova_sig

# Export data
write.csv(anova_sig, file = "Output/Data/Stat/ANOVA_cDNA.csv")

# cDNA Tukey's test -------------------------------------------------------
# Initialize empty data frame
tukey_sig_df <- data.frame()

for (i in 1:nrow(anova_sig)) {
  t <- anova_sig$Time[i]
  
  # Subset data for this time point
  df_time <- cDNA %>% filter(time == t)
  
  # Fit ANOVA
  fit <- aov(cDNA ~ Group_biochar, data = df_time)
  
  # Tukey HSD
  tuk <- TukeyHSD(fit)
  
  # Convert to data frame
  df <- as.data.frame(tuk$Group_biochar) %>%
    mutate(Time = t,
           Comparison = rownames(.)) %>%
    rename(p.adj = `p adj`) %>%
    select(Time, Comparison, diff, lwr, upr, p.adj)
  
  # Split the comparison into two groups (replace '-' with '_' if needed)
  df <- df %>%
    tidyr::separate(Comparison, into = c("Group1", "Group2"), sep = "-") %>%
    mutate(Higher_Group = ifelse(diff > 0, Group1, Group2))  # which has higher mean
  
  # Keep only significant comparisons
  df_sig <- df %>% filter(p.adj < 0.05)
  
  # Combine
  tukey_sig_df <- bind_rows(tukey_sig_df, df_sig)
}

# Arrange and inspect
tukey_sig_df <- tukey_sig_df %>%
  arrange(Time, p.adj)

tukey_sig_df

# Export results
write.csv(tukey_sig_df, file = "Output/Data/Stat/Tukey_cDNA.csv",
          row.names = FALSE)

# DNA ---------------------------------------------------------------------
DNA <- obs.data[, c("Experiment", "Percent_biochar",
                     "Group", "time", "Replicate", "DNA")]
DNA$Group_biochar <- paste(DNA$Group, DNA$Percent_biochar, sep = "-")

# DNA ANOVA ---------------------------------------------------------------
anova_pvalues <- data.frame(
  Time = numeric(),
  p_value = numeric()
)

for (time_pt in time_points) {
  df_time <- DNA %>% filter(time == time_pt)
  
  # Check number of unique groups
  if (length(unique(df_time$Group_biochar)) > 1) {
    fit <- aov(DNA ~ Group_biochar, data = df_time)
    pval <- summary(fit)[[1]][["Pr(>F)"]][1]
    
    anova_pvalues <- rbind(anova_pvalues, data.frame(
      Time = time_pt,
      p_value = pval
    ))
  } else {
    # Skip this time point or record NA
    anova_pvalues <- rbind(anova_pvalues, data.frame(
      Time = time_pt,
      p_value = NA
    ))
  }
}

anova_sig <- anova_pvalues %>% filter(!is.na(p_value) & p_value < 0.05)
anova_sig

# Export data
write.csv(anova_sig, file = "Output/Data/Stat/ANOVA_DNA.csv")

# Tukey's test ------------------------------------------------------------
# Initialize empty data frame
tukey_sig_df <- data.frame()

# Loop through each significant ANOVA result
for (i in 1:nrow(anova_sig)) {
  t <- anova_sig$Time[i]  # time point
  
  # Subset data for this time point
  df_time <- cDNA %>% filter(time == t)
  
  # Fit ANOVA
  fit <- aov(cDNA ~ Group_biochar, data = df_time)
  
  # Tukey HSD
  tuk <- TukeyHSD(fit)
  
  # Convert to data frame
  df <- as.data.frame(tuk$Group_biochar) %>%
    mutate(
      Time = t,
      Comparison = rownames(.)
    ) %>%
    rename(p.adj = `p adj`) %>%
    select(Time, Comparison, diff, lwr, upr, p.adj)
  
  # Split the comparison into two groups
  df <- df %>%
    tidyr::separate(Comparison, into = c("Group1", "Group2"), sep = "-") %>%
    mutate(
      Higher_Group = ifelse(diff > 0, Group1, Group2)  # which has higher mean
    )
  
  # Keep only significant comparisons
  df_sig <- df %>% filter(p.adj < 0.05)
  
  # Combine
  tukey_sig_df <- bind_rows(tukey_sig_df, df_sig)
}

# Arrange and inspect
tukey_sig_df <- tukey_sig_df %>%
  arrange(Time, p.adj)

print(tukey_sig_df)

# Export results
write.csv(tukey_sig_df, file = "Output/Data/Stat/Tukey_DNA.csv",
          row.names = FALSE)

# Transcript gene ratio ---------------------------------------------------
tgr <- obs.data[, c("Experiment", "Percent_biochar",
                     "Group", "time", "Replicate", "transcript_gene_ratio")]
tgr$Group_biochar <- paste(trg$Group, tgr$Percent_biochar, sep = "-")

# ANOVA -------------------------------------------------------------------
anova_pvalues <- data.frame(
  Time = numeric(),
  p_value = numeric()
)

for (time_pt in time_points) {
  df_time <- tgr %>% filter(time == time_pt)
  
  # Check number of unique groups
  if (length(unique(df_time$Group_biochar)) > 1) {
    fit <- aov(tgr ~ Group_biochar, data = df_time)
    pval <- summary(fit)[[1]][["Pr(>F)"]][1]
    
    anova_pvalues <- rbind(anova_pvalues, data.frame(
      Time = time_pt,
      p_value = pval
    ))
  } else {
    # Skip this time point or record NA
    anova_pvalues <- rbind(anova_pvalues, data.frame(
      Time = time_pt,
      p_value = NA
    ))
  }
}

anova_sig <- anova_pvalues %>% filter(!is.na(p_value) & p_value < 0.05)
anova_sig

# Export data
write.csv(anova_sig, file = "Output/Data/Stat/ANOVA_cDNA.csv")

# Tukey's test ------------------------------------------------------------
# Initialize empty data frame
tukey_sig_df <- data.frame()

# Loop through each significant ANOVA result
for (i in 1:nrow(anova_sig)) {
  t <- anova_sig$Time[i]  # time point
  
  # Subset data for this time point
  df_time <- cDNA %>% filter(time == t)
  
  # Fit ANOVA
  fit <- aov(cDNA ~ Group_biochar, data = df_time)
  
  # Tukey HSD
  tuk <- TukeyHSD(fit)
  
  # Convert to data frame
  df <- as.data.frame(tuk$Group_biochar) %>%
    mutate(
      Time = t,
      Comparison = rownames(.)
    ) %>%
    rename(p.adj = `p adj`) %>%
    select(Time, Comparison, diff, lwr, upr, p.adj)
  
  # Split the comparison into two groups
  df <- df %>%
    tidyr::separate(Comparison, into = c("Group1", "Group2"), sep = "-") %>%
    mutate(
      Higher_Group = ifelse(diff > 0, Group1, Group2)  # which has higher mean
    )
  
  # Keep only significant comparisons
  df_sig <- df %>% filter(p.adj < 0.05)
  
  # Combine
  tukey_sig_df <- bind_rows(tukey_sig_df, df_sig)
}

# Arrange and inspect
tukey_sig_df <- tukey_sig_df %>%
  arrange(Time, p.adj)

print(tukey_sig_df)

# Export results
write.csv(tukey_sig_df, file = "Output/Data/Stat/Tukey_PUF.csv",
          row.names = FALSE)







