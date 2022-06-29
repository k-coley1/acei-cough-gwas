# Author: KC
# Description: R script to define ACEI-induced cough phenotype using ACEI to ARB switch.

library(tidyverse)
library(data.table)
library(lubridate)
`%!in%` <- negate(`%in%`)

setwd("~/acei-cough-gwas")

# Load data ----
## Dates data
d <- read_tsv("data/UKB_issue_dob_dod_doe_record_length.txt") %>%
  rename(ACEI_1st_issue_date = issue_date)

## Data frame containing switch pool (individuals with both ACEI and ARB Rx)
ACEI_ARB.df <- read_tsv("data/UKB_all_ACEI_ARB_prescriptions.txt") %>%
  left_join(d, by = "eid") %>%
  arrange(eid, issue_date) %>%
  filter(length_of_records_years >= 2)

## Data frame of all ACEI and ARB Rx 
ACEI <- read_tsv("data/UKB_all_ACEI_prescriptions.txt") %>%
  mutate(drug_type = "ACEI") %>%
  select(-ACEI)
ARB <- read_tsv("data/UKB_all_ARB_prescriptions.txt") %>%
  filter(eid %in% ACEI$eid) %>%           # Filter for only individuals with an ACEI Rx
  mutate(drug_type = "ARB") %>%
  select(-ARB)

## Combined data frame of all ACEI and ARB data
all_ACEI_ARB.df <- bind_rows(ACEI, ARB) %>%
  arrange(eid, issue_date) %>%
  left_join(d, by = "eid") %>%
  filter(length_of_records_years >= 2)


# Data preparation ----
## Identify individuals prescribed their first ACEI and ARB on the same date.
exclude_timing <- ACEI_ARB.df %>%
  group_by(eid) %>%
  filter(issue_date == min(issue_date)) %>%      
  mutate(count.total = n()) %>%           #Events of 2 prescriptions on first date
  filter(count.total > 1) %>%              # More than one prescription on date of first ARB/ACEI
  select(eid, issue_date, drug_type, count.total) %>%
  mutate(count.ACEI = sum(drug_type == "ACEI"),
         count.ARB = sum(drug_type == "ARB")) %>%
  filter(count.ACEI >= 1 & count.ARB >= 1) %>%
  distinct(eid)

## Identify individuals who received at least 1 ARB Rx before their first ACEI Rx.
exclude_ARB_first <- ACEI_ARB.df %>%
  distinct(eid, .keep_all = TRUE) %>%
  filter(drug_type == "ARB") %>%
  distinct(eid)

## Date of first ARB Rx
ARB_1st_issue_date <- ACEI_ARB.df %>%
  filter(drug_type == "ARB") %>%
  group_by(eid) %>%
  mutate(ARB_1st_issue_date = min(issue_date)) %>%
  select(eid, ARB_1st_issue_date) %>%
  distinct(.keep_all = TRUE)


# Perform exclusions  ----
ACEI_ARB_filt.df <- ACEI_ARB.df %>%
  filter(eid %!in% exclude_timing$eid) %>%
  filter(eid %!in% exclude_ARB_first$eid) %>%
  left_join(ARB_1st_issue_date, by = "eid")


# Define phenotype ----
## Cases: Switchers from ACEI to ARB within 12 months of index date 
twelve_months_switch_date <- ACEI_ARB_filt.df %>%
  filter(ARB_1st_issue_date <= (ACEI_1st_issue_date %m+% months(12)),
         drug_type == "ARB") %>%
  distinct(eid, .keep_all = TRUE) %>%
  select(eid, switch_date = issue_date)

### Identify individuals who received another ACEI in the following year after switch date ----
exclude_twelve_months_switch <- twelve_months_switch_date %>%
  left_join(ACEI_ARB_filt.df, by = "eid") %>%
  filter(drug_type == "ACEI",
         issue_date >= switch_date) %>%
  filter(as.numeric(difftime(issue_date, switch_date), unit = "weeks") <= 52) %>% 
  distinct(eid)

### Perform exclusions ----
twelve_months <- twelve_months_switch_date %>%
  filter(eid %!in% exclude_twelve_months_switch$eid) %>%
  distinct(eid) %>%
  mutate(ACEI.cough = 1)

## Controls: continuous users of ACEIs
### ACEI Rx within 12 months of index date AND no ARB Rx within 12 months after index date.
controls <- all_ACEI_ARB.df %>%
  select(eid, issue_date, drug_type, ACEI_1st_issue_date) %>%
  filter(eid %!in% exclude_timing$eid,
         eid %!in% exclude_ARB_first$eid,
         eid %!in% twelve_months$eid,           # Remove individuals already identified as a case
         eid %!in% exclude_twelve_months_switch$eid) %>%     # Remove individuals who switched but didn't fit full case criteria
  mutate(ACEI_cut_off = ACEI_1st_issue_date %m+% months(12)) %>%
  filter(issue_date <= ACEI_cut_off,
         drug_type == "ACEI") %>%
  distinct(eid, issue_date, drug_type, .keep_all = TRUE) %>%         # Remove events on same date
  group_by(eid) %>%
  summarise(count.ACEI = n()) %>%
  filter(count.ACEI >= 2) %>%         #Individuals with at least 2 prescriptions in 12 months.
  mutate(ACEI.cough = 0) %>%
  select(-count.ACEI)


# Create phenotype file ----
bind_rows(twelve_months, controls) %>%
  arrange(eid) %>%
  fwrite("data/UKB_ACEI_cough_pheno.txt", sep = "\t", col.names = TRUE)
