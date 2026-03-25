# Code to calculate at equilibrium, the fraction of individual PCB at:
# sediment, biochar, water, spme, air and puf

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


# Read chemical data
chem.data <- read.csv("Data/PhysicalChemicalProperties.csv")

chem.data <- chem.data %>%
  mutate(
    Kaw = 10^logKaw, # [Lw/La]
    Kow = 10^logKow, # [Lw/Loc]
    Koa = 10^logKoa, # [La/Loc]
    
    Koc = 10^(0.94*logKow + 0.42), # [Lw/Kgoc]
    Kspme = 10^(1.06*logKow - 1.16), # [Lw/Lspme]
    Kpuf = 10^(0.6366*logKoa - 3.1774), # [m3a/gpuf]
    Kpuf = Kpuf * 21300, # [La/Lpuf], density of puf in gpuf/m3a
    Kbc = 10^(0.8737*log10(Koc) + 0.263)) # [Lw/Kgbc], regression from Qin, only 3 congeners (4, 18, 52)

ms    <- 0.01 # [Kgs]
mbc   <- ms * 0.05 # [Kgbc] 5% of sediment
Vw    <- 0.1 # [Lw]
Vspme_Lspme <- 6.9e-8 # [Lspme/cmspme]
Lspme <- 1 # [cmspme] considering 1 cm length
Vspme <- Vspme_Lspme * Lspme # [Lspme]
Va    <- 0.125 # [La]
Vpuf  <- 0.029 # [Lpuf]
foc   <- 0.064 # [Kgbc/Kgs]

# No biochar
fract <- chem.data %>%
  mutate(
    # capacity in each phase in relation to water
    As    = foc * Koc * ms / Vw,
    Aspme = Kspme * Vspme / Vw,
    Aair  = Kaw * Va / Vw,
    Apuf  = Kpuf * Kaw * Vpuf / Vw,
    D     = 1 + As + Aspme + Aair + Apuf,
    
    fractw    = 1 / D,
    fracts    = As / D,
    fractspme = Aspme / D,
    fracta    = Aair / D,
    fractpuf  = Apuf / D,
    sumfrac   = fractw + fracts + fractspme + fracta + fractpuf) %>%
  select(congener, fracts, fractw, fractspme, fracta, fractpuf, sumfrac)

# Export data
write.csv(fract, file = "Output/Data/EqPModel/fractionNobc.csv")

# With biochar
fract_bc <- chem.data %>%
  mutate(
    # capacity in each phase in relation to water
    As    = foc * Koc * ms / Vw,
    Abc   = Kbc * mbc / Vw,
    Aspme = Kspme * Vspme / Vw,
    Aair  = Kaw * Va / Vw,
    Apuf  = Kpuf * Kaw * Vpuf / Vw,
    D     = 1 + As + Abc + Aspme + Aair + Apuf,
    
    fractw    = 1 / D,
    fracts    = As / D,
    fractbc   = Abc / D,
    fractspme = Aspme / D,
    fracta    = Aair / D,
    fractpuf  = Apuf / D,
    sumfrac   = fractw + fracts + fractbc + fractspme + fracta + fractpuf) %>%
  select(congener, fracts, fractbc, fractw, fractspme, fracta, fractpuf, sumfrac)

# Export data
write.csv(fract_bc, file = "Output/Data/EqPModel/fractionbc.csv")

# With biochar 1.5%
mbc1.5   <- ms * 0.015 # [Kgbc] 1.5% of sediment

fract_bc1.5 <- chem.data %>%
  mutate(
    # capacity in each phase in relation to water
    As    = foc * Koc * ms / Vw,
    Abc   = Kbc * mbc1.5 / Vw,
    Aspme = Kspme * Vspme / Vw,
    Aair  = Kaw * Va / Vw,
    Apuf  = Kpuf * Kaw * Vpuf / Vw,
    D     = 1 + As + Abc + Aspme + Aair + Apuf,
    
    fractw    = 1 / D,
    fracts    = As / D,
    fractbc   = Abc / D,
    fractspme = Aspme / D,
    fracta    = Aair / D,
    fractpuf  = Apuf / D,
    sumfrac   = fractw + fracts + fractbc + fractspme + fracta + fractpuf) %>%
  select(congener, fracts, fractbc, fractw, fractspme, fracta, fractpuf, sumfrac)

# Export data
write.csv(fract_bc1.5, file = "Output/Data/EqPModel/fractionbc1_5.csv")

# With biochar 10%
mbc10   <- ms * 0.1 # [Kgbc] 10% of sediment

fract_bc10 <- chem.data %>%
  mutate(
    # capacity in each phase in relation to water
    As    = foc * Koc * ms / Vw,
    Abc   = Kbc * mbc10 / Vw,
    Aspme = Kspme * Vspme / Vw,
    Aair  = Kaw * Va / Vw,
    Apuf  = Kpuf * Kaw * Vpuf / Vw,
    D     = 1 + As + Abc + Aspme + Aair + Apuf,
    
    fractw    = 1 / D,
    fracts    = As / D,
    fractbc   = Abc / D,
    fractspme = Aspme / D,
    fracta    = Aair / D,
    fractpuf  = Apuf / D,
    sumfrac   = fractw + fracts + fractbc + fractspme + fracta + fractpuf) %>%
  select(congener, fracts, fractbc, fractw, fractspme, fracta, fractpuf, sumfrac)

# Export data
write.csv(fract_bc10, file = "Output/Data/EqPModel/fractionbc10.csv")


