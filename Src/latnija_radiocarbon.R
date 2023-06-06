# load libraries
library(tidyverse)
library(readxl)

# data wrangling

# pull in the data from Excel sheets and CSVs as needed

# get sheet names
data_path <- "./Data/radiocarbon_latnija_062023.xlsx"
sheets <- excel_sheets(data_path)
sheets

# pull the variables sheet ("Measures")
latnija <- read_excel(data_path, 
                    sheet = sheets[1])

# extract necessary c14 data for later modelling in OxCal

charcoal_idx <- grep("char|Char", latnija$user_desc1)

cols <- c("sample_nr", "C14_age", "C14_age_sig")

latnija_charcoal <- latnija[charcoal_idx, cols]

latnija_charcoal[, c(2, 3)] <- round(latnija_charcoal[, c(2, 3)])

write.csv(latnija_charcoal, 
        file = "Data/lantnija_charcoal.csv", 
        row.names = FALSE)

# for comparing contexts
cols <- c("sample_nr", "C14_age", "C14_age_sig", "Context")

latnija[charcoal_idx, cols]

# phasing

# context numbers grouped by phases

phase_v_e <- c(54, 51, 50, 49, 47, 46, 43, 42, 40, 39, 33)
phase_v_d <- c(44, 45, 48)
phase_v_c <- 
