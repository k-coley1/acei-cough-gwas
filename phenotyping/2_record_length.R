# Author: KC
# Description: R script to generate a text file with length of records data.

library(tidyverse)
library(data.table)
library(lubridate)

setwd("~/acei-cough-gwas")

# Load data ----
## First ACEI prescription
entry1_ACEI <- read_tsv("data/UKB_first_ACEI_prescription.txt") %>%
  select(eid, issue_date)

## Baseline age data generated from questionnaire data: contains eid and date of birth.
dob <- read_tsv("data/ukb_dob.txt")

## Death registry data
deaths <- as_tibble(fread("data/ukb_death_data.tab.gz")) %>%
  select(f.eid, f.40000.0.0, f.40000.1.0) %>%
  mutate(date_deceased = coalesce(f.40000.0.0, f.40000.1.0)) %>%
  rename(eid = f.eid) %>%
  select(eid, date_deceased) %>%
  drop_na()


# Calculate age at first ACEI prescription ----
entry1_ACEI_age <- entry1_ACEI %>%
  left_join(dob, by = "eid") %>%
  mutate(int1 = interval(ymd(dob), ymd(issue_date)),
         age_presc = trunc(time_length(int1, "year"))) %>%
  select(eid, issue_date, age_presc)


# Determine length of records ----
## Firstly, assume that extraction date was 6 months before notable decrease in records before extraction date i.e. 2015-12-01.
## Determine the event which happened earliest in time (i.e. first) out of date of death and record extraction.
## Add column "record_end" which is either DOD or DOE, whichever comes first in time. Also, determine the length of records until this point, in days and years.
d <- entry1_ACEI_age %>%
  mutate(date_extracted = as.Date("2015-12-01")) %>%
  left_join(deaths, by = "eid") %>%
  mutate(record_end = pmin(date_extracted, date_deceased, na.rm = TRUE), #Earliest date out of date of death (if present) and date of extraction.
         length_of_records_days = as.numeric(difftime(record_end, issue_date, units = "days")),
         length_of_records_years = as.numeric(difftime(as.Date(record_end), as.Date(issue_date), unit = "weeks"))/52.25)

## Write file
write_tsv(d, "data/UKB_issue_dob_dod_doe_record_length.txt")