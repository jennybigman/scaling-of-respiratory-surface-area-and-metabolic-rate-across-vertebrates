# Bigman et al. metabolic rate and respiratory surface area scaling Science Advances

# 00-data wrangling

library(lme4)
library(lmerTest)
library(MuMIn)
library(ggplot2)
library(scales)
library(calibrate)
library(magicaxis)
library(grid)
library(dplyr)
library(scales)
library(AICcmodavg)
library(arm)
library(ape)
library(geiger)
library(phytools)
library(nlme)
library(brms)
library(MCMCglmm)
library(tidybayes)
library(tidyr)
library(modelr)
library(bayesplot)
library(rstan) # observe startup messages
library(gdata)
library(beepr)
library(loo)
library(here)

dataMR <- read.csv(here("PairedMeanMRData.csv"))
dataRSA <- read.csv(here("PairedRSAData.csv"))

## add a column with both genus & species names ##
dataMR$Binomial <- paste(dataMR$Genus, dataMR$Species,sep = " ")
str(dataMR)
names(dataMR)

dataRSA$Binomial <- paste(dataRSA$Genus, dataRSA$Species,sep = " ")
str(dataRSA)

## standardize units of metabolic rate ##

dataMR <- dataMR %>% mutate(MeanMassKG = MeanMassG/1000)

dataMR <- dataMR %>% mutate(MeanWholeOrganismMRTime = case_when(
  MeanMRUnits %in% "mg/h" ~ MeanMR,
  MeanMRUnits %in% "mg/g/h" ~ MeanMR * MeanMassG,
  MeanMRUnits %in% "mg/kg/h" ~ MeanMR * MeanMassKG,
  MeanMRUnits %in% "ml/h" ~ MeanMR,
  MeanMRUnits %in% "ml/min" ~ MeanMR,
  MeanMRUnits %in% "mg/min" ~ MeanMR,
  MeanMRUnits %in% "ml/g/h" ~ MeanMR * MeanMassG,
  MeanMRUnits %in% "ml/kg/h" ~ MeanMR * MeanMassKG,
  MeanMRUnits %in% "ml/kg/min" ~ MeanMR * MeanMassKG,
  MeanMRUnits %in% "joules/h" ~ MeanMR,
  MeanMRUnits %in% "watts" ~ MeanMR))

dataMR <- dataMR %>% mutate(MeanWholeOrganismMRTimeUnits = case_when(
  MeanMRUnits %in% "mg/h" ~ "mg/h",
  MeanMRUnits %in% "mg/g/h" ~ "mg/h",
  MeanMRUnits %in% "mg/kg/h" ~ "mg/h",
  MeanMRUnits %in% "ml/h" ~ "ml/h",
  MeanMRUnits %in% "ml/min" ~ "ml/min",
  MeanMRUnits %in% "mg/min" ~ "mg/min",
  MeanMRUnits %in% "ml/g/h" ~ "ml/h",
  MeanMRUnits %in% "ml/kg/h" ~ "ml/h",
  MeanMRUnits %in% "ml/kg/min" ~ "ml/min",
  MeanMRUnits %in% "joules/h" ~ "joules/h",
  MeanMRUnits %in% "watts" ~ "watts"))

dataMR <- dataMR %>% mutate(
  MeanWholeOrganismMRWatts = case_when(
    MeanWholeOrganismMRTimeUnits %in% "ml/min" ~ (MeanWholeOrganismMRTime / 60) * 20.1,
    MeanWholeOrganismMRTimeUnits %in% "ml/h" ~ (MeanWholeOrganismMRTime / 3600) * 20.1,
    MeanWholeOrganismMRTimeUnits %in% "mg/min" ~ (MeanWholeOrganismMRTime / 60) * 14.1,
    MeanWholeOrganismMRTimeUnits %in% "mg/h" ~ (MeanWholeOrganismMRTime / 3600) * 14.1,
    MeanWholeOrganismMRTimeUnits %in% "joules/h" ~ (MeanWholeOrganismMRTime / 3600),
    MeanWholeOrganismMRTimeUnits %in% "watts" ~ MeanWholeOrganismMRTime))

names(dataMR)

## convert all respiratory surface area data to cm2 ##

dataRSA <- dataRSA %>% mutate(
  RSAcm2 = case_when(
    MeanRSAUnits %in% "cm2" ~ MeanRSA,
    MeanRSAUnits %in% "mm2" ~ (MeanRSA / 100)
    ))

## data transformations ##

# 1.log respiratory surface area
dataRSA <- dataRSA %>% mutate(LogRSAcm2 = log(RSAcm2))

# 2.log body mass associated with respiratory surface area 
dataRSA <- dataRSA %>% mutate(RSA.LogMeanMassG = log(MeanMassG))

# 3.log metabolic rate (in watts)
dataMR <- dataMR %>% mutate(LogMeanWholeOrganismMRWatts = log(MeanWholeOrganismMRWatts))

# 4. log body mass associated with metabolic rate 
dataMR <- dataMR %>% mutate(MR.LogMeanMassG = log(MeanMassG))

names(dataMR)
names(dataRSA)

## center log mass for both metabolic rate & respiratory surface area ##

dataMR <- dataMR %>% mutate(LogCenteredMeanMassMR = MR.LogMeanMassG - (mean(dataMR$MR.LogMeanMassG)))

dataRSA <- dataRSA %>% mutate(LogCenteredMeanMassRSA = RSA.LogMeanMassG - (mean(dataRSA$RSA.LogMeanMassG)))

## Temperature transformations

#Boltzmann constant
Bol <- 8.617343e-05

## convert metabolic rate temp into kelvin and then into inverse temperature ##
dataMR <- dataMR %>% mutate(MR.TempK = Temperature + 273.15)

dataMR <- dataMR %>% mutate(MR.InverseTemp = (1/(Bol*MR.TempK)))


## standardize temperature associated with metabolic rate ##

CS.MR.Temp <- arm::rescale(dataMR$Temperature, "full")
dataMR <- dataMR %>% mutate(Scaled.Centered.MR.Temperature = CS.MR.Temp)

CS.MR.InverseTemp <- arm::rescale(dataMR$MR.InverseTemp, "full")
dataMR <- dataMR %>% mutate(Scaled.Centered.MR.InverseTemp = CS.MR.InverseTemp)

STD_mass_MR <- arm::rescale(dataMR$MR.LogMeanMassG, "full")
dataMR <- dataMR %>% mutate(STD_mass_MR = STD_mass_MR)

STD_mass_RSA <- arm::rescale(dataRSA$RSA.LogMeanMassG, "full")
dataRSA <- dataRSA %>% mutate(STD_mass_RSA = STD_mass_RSA)

## standardize data using one sd (z score) ##

dataMR <- dataMR %>% mutate(std_temp = (MR.InverseTemp - mean(MR.InverseTemp))/ sd(MR.InverseTemp))

dataMR <- dataMR %>% mutate(std_MMR = (MR.LogMeanMassG - mean(MR.LogMeanMassG))/ sd(MR.LogMeanMassG))

dataRSA <- dataRSA %>% mutate(std_MRSA = (RSA.LogMeanMassG - mean(RSA.LogMeanMassG))/ sd(RSA.LogMeanMassG))


## rename some metabolic rate columns for when databases are merged ##
names(dataMR)
colnames(dataMR)[colnames(dataMR) == "Taxon"] <- "MR.Taxon"
colnames(dataMR)[colnames(dataMR) == "MeanMassG"] <- "MR.MeanMassG"
colnames(dataMR)[colnames(dataMR) == "Genus"] <- "MR.Genus"

## rename some respiratory surface area columns for when databases are merged ##
names(dataRSA)
colnames(dataRSA)[colnames(dataRSA) == "Taxon"] <- "RSA.Taxon"
colnames(dataRSA)[colnames(dataRSA) == "MeanMassG"] <- "RSA.MeanMassG"

## merge datasets ##
OverlapMR <- dataMR %>% dplyr::select(MR.Taxon, MR.Genus, Binomial, CommonName,
                                      ThermoStrat, MR.MeanMassG, MR.LogMeanMassG,
                                      MeanWholeOrganismMRWatts, STD_mass_MR, std_temp,
                                      LogMeanWholeOrganismMRWatts, std_MMR,
                                      MR.InverseTemp, LogCenteredMeanMassMR,
                                      Scaled.Centered.MR.InverseTemp)

OverlapRSA <- dataRSA %>% dplyr::select(RSA.Taxon, Binomial, LungsGills,
                                        RSA.MeanMassG, RSA.LogMeanMassG,
                                        RSAcm2, LogRSAcm2, STD_mass_RSA,
                                        std_MRSA, LogCenteredMeanMassRSA)

Overlap.MR.RSA <- merge(OverlapMR, OverlapRSA, by = "Binomial",
                        stringsAsFactors = FALSE)

names(Overlap.MR.RSA)

## dummy code thermoregulatory strategy & respiratory organ ##

Overlap.MR.RSA$ThermoStrat[Overlap.MR.RSA$ThermoStrat == "ectotherm"] <- 0
Overlap.MR.RSA$ThermoStrat[Overlap.MR.RSA$ThermoStrat == "endotherm"] <- 1

Overlap.MR.RSA$LungsGills[Overlap.MR.RSA$LungsGills == "gills"] <- 0
Overlap.MR.RSA$LungsGills[Overlap.MR.RSA$LungsGills == "lungs"] <- 1

str(Overlap.MR.RSA)

Overlap.MR.RSA$ThermoStrat <- as.numeric(Overlap.MR.RSA$ThermoStrat)
Overlap.MR.RSA$LungsGills <- as.numeric(Overlap.MR.RSA$LungsGills)


summary <- Overlap.MR.RSA %>%
           group_by(MR.Taxon) %>%
           summarise(n())
