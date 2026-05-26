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
# Latest file is version 5. It includes the biochar_percentage 5% for 90 days.
obs.data <- read.csv("Data/uncoated_biochar_V5.csv")

# Biochar time series analysis --------------------------------------------
# Format data
# Select SPME & time series data
pcbi.spme <- obs.data %>%
  filter(Sample_medium == "SPME",
         (Experiment == "biochar_timeseries") |
           (Experiment == "biochar_percentage" &
              percent_biochar %in% c(0, 5) &
              time == 90))

pcbi.spme$Group_biochar <- paste(pcbi.spme$Group, pcbi.spme$percent_biochar,
                                sep = "_")

pcbi.spme$Group_biochar <- factor(pcbi.spme$Group_biochar)

# ANOVA 1------------------------------------------------------------------
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
write.csv(anova_sig, file = "Output/Data/Stat/ANOVA_SPMEV3.csv")

# Tukey's test 1 ------------------------------------------------------------
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
write.csv(tukey_sig_df, file = "Output/Data/Stat/Tukey_SPMEV3.csv",
          row.names = FALSE)

# Biochar percentage analysis: Control -----------------------------------
# Select PUF & time series data
pcbi.spme.bc <- obs.data %>%
  filter(Sample_medium =="SPME", Experiment =="biochar_percentage",
         Group == "Control")

pcbi.spme.bc$Group_biochar <- paste(pcbi.spme.bc$Group,
                                   pcbi.spme.bc$percent_biochar, sep = "_")
pcbi.spme.bc$Group_biochar <- factor(pcbi.spme.bc$Group_biochar)

# ANOVA 2 ------------------------------------------------------------------
pcb.cols <- grep("^PCB_", names(pcbi.spme.bc), value = TRUE)

all_pcb <- c(pcb.cols, "tPCB", "LCPCB")

anova_pvalues <- data.frame(
  Variable = character(),
  p_value = numeric()
)

for (col in all_pcb) {
  fit <- aov(
    as.formula(paste(col, "~ factor(percent_biochar)")),
    data = pcbi.spme.bc %>%
      filter(!is.na(.data[[col]]))
  )
  pval <- summary(fit)[[1]][["Pr(>F)"]][1]
  anova_pvalues <- rbind(
    anova_pvalues,
    data.frame(
      Variable = col,
      p_value = pval
    )
  )
}

anova_sig <- anova_pvalues %>%
  filter(p_value < 0.05)

print(anova_sig)

# Export data
write.csv(anova_sig, file = "Output/Data/Stat/ANOVA_SPMEBCV1.csv")

# Tukey's test 2 -----------------------------------------------------------
# Initialize empty data frame
tukey_sig_df <- data.frame()

# Loop through each significant ANOVA result
for (i in 1:nrow(anova_sig)) {
  var <- anova_sig$Variable[i]
  # Remove NA values
  df_tmp <- pcbi.spme.bc %>%
    filter(!is.na(.data[[var]]))
  # ANOVA model
  fit <- aov(
    as.formula(paste(var, "~ factor(percent_biochar)")),
    data = df_tmp
  )
  # Tukey test
  tuk <- TukeyHSD(fit)
  # Convert to dataframe
  df <- as.data.frame(tuk$`factor(percent_biochar)`) %>%
    mutate(
      Variable = var,
      Comparison = rownames(.)
    ) %>%
    rename(p.adj = `p adj`) %>%
    select(Variable, Comparison, diff, lwr, upr, p.adj)
  # Separate comparison groups
  df <- df %>%
    tidyr::separate(
      Comparison,
      into = c("Biochar1", "Biochar2"),
      sep = "-"
    ) %>%
    mutate(
      Higher_Group = ifelse(diff > 0, Biochar1, Biochar2)
    )
  # Keep significant comparisons only
  df_sig <- df %>%
    filter(p.adj < 0.05)
  tukey_sig_df <- bind_rows(tukey_sig_df, df_sig)
}

# Arrange results
tukey_sig_df <- tukey_sig_df %>%
  arrange(Variable, p.adj)

print(tukey_sig_df)

# Export results
write.csv(tukey_sig_df, file = "Output/Data/Stat/Tukey_SPMEBCV1.csv",
          row.names = FALSE)

# Biochar percentage analysis: LB400 -------------------------------------
pcbi.spme.bc.lb400 <- obs.data %>%
  filter(Sample_medium =="PUF", Experiment =="biochar_percentage",
         Group == "Treatment")

pcbi.spme.bc.lb400$Group_biochar <- paste(pcbi.spme.bc.lb400$Group,
                                         pcbi.spme.bc.lb400$percent_biochar,
                                         sep = "_")
pcbi.spme.bc.lb400$Group_biochar <- factor(pcbi.spme.bc.lb400$Group_biochar)

# ANOVA 2 ------------------------------------------------------------------
pcb.cols <- grep("^PCB_", names(pcbi.spme.bc.lb400), value = TRUE)

all_pcb <- c(pcb.cols, "tPCB", "LCPCB")

anova_pvalues <- data.frame(
  Variable = character(),
  p_value = numeric()
)

for (col in all_pcb) {
  fit <- aov(
    as.formula(paste(col, "~ factor(percent_biochar)")),
    data = pcbi.spme.bc.lb400 %>%
      filter(!is.na(.data[[col]]))
  )
  pval <- summary(fit)[[1]][["Pr(>F)"]][1]
  anova_pvalues <- rbind(
    anova_pvalues,
    data.frame(
      Variable = col,
      p_value = pval
    )
  )
}

anova_sig <- anova_pvalues %>%
  filter(p_value < 0.05)

print(anova_sig)

# Export data
write.csv(anova_sig, file = "Output/Data/Stat/ANOVA_SPMEBCLB400V1.csv")

# Tukey's test 3 -----------------------------------------------------------
# Initialize empty data frame
tukey_sig_df <- data.frame()

# Loop through each significant ANOVA result
for (i in 1:nrow(anova_sig)) {
  var <- anova_sig$Variable[i]
  # Remove NA values
  df_tmp <- pcbi.spme.bc.lb400 %>%
    filter(!is.na(.data[[var]]))
  # ANOVA model
  fit <- aov(
    as.formula(paste(var, "~ factor(percent_biochar)")),
    data = df_tmp
  )
  # Tukey test
  tuk <- TukeyHSD(fit)
  # Convert to dataframe
  df <- as.data.frame(tuk$`factor(percent_biochar)`) %>%
    mutate(
      Variable = var,
      Comparison = rownames(.)
    ) %>%
    rename(p.adj = `p adj`) %>%
    select(Variable, Comparison, diff, lwr, upr, p.adj)
  # Separate comparison groups
  df <- df %>%
    tidyr::separate(
      Comparison,
      into = c("Biochar1", "Biochar2"),
      sep = "-"
    ) %>%
    mutate(
      Higher_Group = ifelse(diff > 0, Biochar1, Biochar2)
    )
  # Keep significant comparisons only
  df_sig <- df %>%
    filter(p.adj < 0.05)
  tukey_sig_df <- bind_rows(tukey_sig_df, df_sig)
}

# Arrange results
tukey_sig_df <- tukey_sig_df %>%
  arrange(Variable, p.adj)

print(tukey_sig_df)

# Export results
write.csv(tukey_sig_df, file = "Output/Data/Stat/Tukey_SPMEBCLB400V1.csv",
          row.names = FALSE)

# Biochar percentage analysis: %BC vs LB400 -------------------------------
# Select PCB columns
pcb.cols <- grep("^PCB_", names(obs.data), value = TRUE)

all_pcb <- c(pcb.cols, "tPCB", "LCPCB")

# Filter dataset:
# - SPME only
# - biochar_percentage experiment only
# - remove 0% biochar
pcbi.spme.bc.lb400.2 <- obs.data %>%
  filter(
    Sample_medium == "SPME",
    Experiment == "biochar_percentage",
    percent_biochar != 0
  )

# Ensure factors
pcbi.spme.bc.lb400.2$Group <- factor(pcbi.spme.bc.lb400.2$Group)
pcbi.spme.bc.lb400.2$percent_biochar <- factor(pcbi.spme.bc.lb400.2$percent_biochar)

# Initialize results dataframe
anova4_results <- data.frame()

# Loop through all PCB variables
for (col in all_pcb) {
  # Remove NA values
  df_tmp <- pcbi.spme.bc.lb400.2 %>%
    filter(!is.na(.data[[col]]))
  # Two-way ANOVA
  fit <- aov(
    as.formula(
      paste(col, "~ Group * percent_biochar")
    ),
    data = df_tmp
  )
  # Extract ANOVA table
  anova_tab <- summary(fit)[[1]]
  # Store results
  tmp <- data.frame(
    Variable = col,
    Effect = rownames(anova_tab),
    p_value = anova_tab$`Pr(>F)`
  )
  anova4_results <- bind_rows(anova4_results, tmp)
}

# Keep significant results
anova4_sig <- anova4_results %>%
  filter(
    !is.na(p_value),
    p_value < 0.05
  )

# View
print(anova4_sig)

# Notes. If p-value < 0.05
# Group: Control ≠ Treatment overall
# percent_biochar: Different biochar percentages produce different PCB levels
# Group:percent_biochar: The effect of biochar depends on the Group

# Export data
write.csv(anova4_sig, file = "Output/Data/Stat/2_way_ANOVA_SPMEBCLB400V1.csv")

# Tukey's test 4 ----------------------------------------------------------
# Create combined factor
pcbi.spme.bc.lb400.2$Group_Biochar <- interaction(
  pcbi.spme.bc.lb400.2$Group,
  pcbi.spme.bc.lb400.2$percent_biochar, sep = "_")

# Initialize dataframe
tukey_sig_df <- data.frame()

# Loop through all PCB variables 
for (col in all_pcb) {
  # Remove NA values
  df_tmp <- pcbi.spme.bc.lb400.2 %>%
    filter(!is.na(.data[[col]]))
  # One-way ANOVA using combined factor
  fit <- aov(
    as.formula(
      paste(col, "~ Group_Biochar")
    ),
    data = df_tmp
  )
  # Tukey test
  tuk <- TukeyHSD(fit)
  # Convert Tukey output to dataframe
  df <- as.data.frame(tuk$Group_Biochar) %>%
    mutate(
      Variable = col,
      Comparison = rownames(.)
    ) %>%
    rename(p.adj = `p adj`) %>%
    select(
      Variable,
      Comparison,
      diff,
      lwr,
      upr,
      p.adj
    )
  # Split comparison names
  df <- df %>%
    separate(
      Comparison,
      into = c("Group1", "Group2"),
      sep = "-"
    )
  # Identify higher mean group
  df <- df %>%
    mutate(
      Higher_Group = ifelse(
        diff > 0,
        Group1,
        Group2
      )
    )
  # Keep significant comparisons only
  df_sig <- df %>%
    filter(p.adj < 0.05)
  # Store results
  tukey_sig_df <- bind_rows(
    tukey_sig_df,
    df_sig
  )
}

# Arrange results
tukey_sig_df <- tukey_sig_df %>%
  arrange(Variable, p.adj)

#View results
print(tukey_sig_df)

# Export Tukey results
write.csv(tukey_sig_df, file = "Output/Data/Stat/Tukey_2_way_ANOVA_SPMEBCLB400V1.csv",
          row.names = FALSE)




