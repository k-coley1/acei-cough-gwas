# Author: KC
# Description: R script to extract ACEI and ARB prescription events from UK Biobank primary care data.

library(tidyverse)
library(data.table)
`%!in%` <- negate(`%in%`)

setwd("~/acei-cough-gwas")

# Read data ----
## Withdrawn individuals
withdrawn <- read_csv("data/withdrawn.csv", col_names = FALSE)

## Prescription data
Rx <- read_tsv("data/gp_scripts.txt.gz",
               col_types = cols(eid = col_integer(),
                                data_provider = col_factor(levels = c("1", "2", "3", "4")),
                                issue_date = col_date(format = "%d/%m/%Y"),
                                read_2 = col_character(),
                                bnf_code = col_character(),
                                dmd_code = col_character(),
                                drug_name = col_character(),
                                quantity = col_character())) %>%
  filter(eid %!in% withdrawn$X1)

## ACEI code lists
ACEI_read2_codes <- read_csv("code_lists/ACEI_codes_read_v2.csv")
ACEI_BNF_codes <- read_csv("code_lists/ACEI_codes_bnf.csv")
ACEI_dmd_codes <- read_csv("code_lists/ACEI_codes_dmd.csv")
ACEI_names <- read_csv("code_lists/ACEI_drug_names.csv")

## ARB code lists
ARB_read2_codes <- read_csv("code_lists/ARB_codes_read_v2.csv")
ARB_BNF_codes <- read_csv("code_lists/ARB_codes_bnf.csv")
ARB_dmd_codes <- read_csv("code_lists/ARB_codes_dmd.csv")
ARB_names <- read_csv("code_lists/ARB_drug_names.csv")


#------------------------------
# ACEI prescription events ----
#------------------------------

# Extract ACEI prescription events ----
## Generate vectors for each coding system and terms
ACEI_read2_search <- paste(as_vector(ACEI_read2_codes$Read_v2), collapse = "|")
ACEI_terms_search <- paste0('(?:\\s|^)', as_vector(ACEI_names$name), '(?:\\s|$)', collapse = "|")

## Extract ACEI Rx events with Read v2 codes
ACEI_read2 <- Rx %>%
  filter(grepl(ACEI_read2_search, read_2, perl = TRUE))

## Extract ACEI Rx events with drug names
ACEI_terms <- Rx %>%
  filter(grepl(ACEI_terms_search, drug_name, perl = TRUE, ignore.case = TRUE))

## Bind rows then remove duplicate entries obtained by more than 1 search
ACEI <- rbind(ACEI_read2, ACEI_terms) %>%
  arrange(eid, issue_date) %>%
  distinct()


# Cleaning ----
## Add drug name (either using regular expressions or using read version 2 [data provider 4 only])
ACEI <- ACEI %>%
  mutate(drug_name = str_to_sentence(drug_name, locale = "en"),          # Title case to allow mapping
         ACEI = str_extract(drug_name, '[A-Za-z]+'),          # Extract first word of drug name
         ACEI = case_when(ACEI %in% (ACEI_names %>% filter(ACEI == "Benazepril"))$name ~ "Benazepril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Captopril"))$name ~ "Captopril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Cilazapril"))$name ~ "Cilazapril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Enalapril"))$name ~ "Enalapril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Fosinopril"))$name ~ "Fosinopril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Imidapril"))$name ~ "Imidapril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Lisinopril"))$name ~ "Lisinopril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Moexipril"))$name ~ "Moexipril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Perindopril"))$name ~ "Perindopril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Quinapril"))$name ~ "Quinapril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Ramipril"))$name ~ "Ramipril",
                          ACEI %in% (ACEI_names %>% filter(ACEI == "Trandolapril"))$name ~ "Trandolapril"))

## Note NA present for those without drug descriptions (data provider 4), add ACEI name using code lists.
ACEI_read2_codes <- ACEI_read2_codes %>%
  rename("read_2" = "Read_v2",
         "drug_name" = "Description") %>%
  select(read_2, drug_name, ACEI)

## Entries defined using Read 2 codes
ACEI_read2_defined <- ACEI %>%
  filter(is.na(drug_name)) %>%
  select(-drug_name, -ACEI) %>%
  left_join(ACEI_read2_codes, by = "read_2") %>%
  mutate(drug_name = str_to_sentence(drug_name, locale = "en")) %>%
  relocate(drug_name, .before = quantity)

## Entries defined using terms
ACEI_term_defined <- ACEI %>%
  filter(!is.na(drug_name))


# Final datasets ----
## All prescriptions
ACEI <- rbind(ACEI_term_defined, ACEI_read2_defined) %>%
  arrange(eid, issue_date) %>%
  mutate(ACEI = case_when(drug_name %in% "Felodipine 5mg modified-release / ramipril 5mg tablets" | drug_name %in% "Felodipine 2.5mg modified-release / ramipril 2.5mg tablets" | drug_name %in% "Felodipine + ramipril mr tab 5mg + 5mg" ~ "Ramipril",
                          drug_name %in% "Verapamil 180mg modified-release / trandolapril 2mg capsules" | drug_name %in% "Verapamil hcl + trandolapril mr cap 180mg + 2mg" ~ "Trandolapril",
                          drug_name %in% "Co-zidocapt (hydrochlorothiazide & captopril) tabs 25mg+50mg" | drug_name %in% "Co-zidocapt 12.5mg/25mg tablets" | drug_name %in% "Co-zidocapt 25/50 tablets" | drug_name %in% "Co-zidocapt 25mg/50mg tablets" | drug_name %in% "Co-zidocapt tabs 25mg + 50mg" | drug_name %in% "Co-zidocapt tabs 50mg + 25mg" ~ "Captopril",
                          T ~ ACEI))            # Sort out residual drug names without correct ACEI description

## First prescription
entry1_ACEI <- ACEI %>%
  distinct(eid, .keep_all = TRUE)


# Write files ----
fwrite(ACEI, "data/UKB_all_ACEI_prescriptions.txt", sep = "\t", col.names = TRUE)
fwrite(entry1_ACEI, "data/UKB_first_ACEI_prescription.txt", sep = "\t", col.names = TRUE)


#-----------------------------
# ARB prescription events ----
#-----------------------------

# Extract ARB prescription events ----
## Generate vectors for each coding system and terms
ARB_read2_search <- paste(as_vector(ARB_read2_codes$Read_v2), collapse = "|")
ARB_terms_search <- paste0('(?:\\s|^)', as_vector(ARB_names$name), '(?:\\s|$)', collapse = "|")

## Extract ARB Rx events with Read v2 codes
ARB_read2 <- Rx %>%
  filter(grepl(ARB_read2_search, read_2, perl = TRUE))

## Extract ARB Rx events with drug names
ARB_terms <- Rx %>%
  filter(grepl(ARB_terms_search, drug_name, perl = TRUE, ignore.case = TRUE))

## Bind rows then remove duplicate entries obtained by more than 1 search
ARB <- rbind(ARB_read2, ARB_terms) %>%
  arrange(eid, issue_date) %>%
  distinct()


# Cleaning ----
## Add drug name (either using regular expressions or using read version 2 [data provider 4 only])
ARB <- ARB %>%
  mutate(drug_name = str_to_sentence(drug_name, locale = "en"),
         ARB = str_extract(drug_name, '[A-Za-z]+'),
         ARB = case_when(ARB %in% (ARB_names %>% filter(ARB == "Azilsartan"))$name ~ "Azilsartan",
                         ARB %in% (ARB_names %>% filter(ARB == "Candesartan"))$name ~ "Candesartan",
                         ARB %in% (ARB_names %>% filter(ARB == "Eprosartan"))$name ~ "Eprosartan",
                         ARB %in% (ARB_names %>% filter(ARB == "Irbesartan"))$name ~ "Irbesartan",
                         ARB %in% (ARB_names %>% filter(ARB == "Losartan"))$name ~ "Losartan",
                         ARB %in% (ARB_names %>% filter(ARB == "Olmesartan"))$name ~ "Olmesartan",
                         ARB %in% (ARB_names %>% filter(ARB == "Telmisartan"))$name ~ "Telmisartan",
                         ARB %in% (ARB_names %>% filter(ARB == "Valsartan"))$name ~ "Valsartan"))

## Note NA present for those without drug descriptions (data provider 4), add ARB name using code lists.
ARB_read2_codes <- ARB_read2_codes %>%
  rename("read_2" = "Read_v2",
         "drug_name" = "Description") %>%
  select(read_2, drug_name, ARB)

## Entries defined using Read 2 codes
ARB_read2_defined <- ARB %>%
  filter(is.na(drug_name)) %>%
  select(-drug_name, -ARB) %>%
  left_join(ARB_read2_codes, by = "read_2") %>%
  mutate(drug_name = str_to_sentence(drug_name, locale = "en")) %>%
  relocate(drug_name, .before = quantity)

## Entries defined using terms
ARB_term_defined <- ARB %>%
  filter(!is.na(drug_name))


# Final datasets ----
## All prescriptions
ARB <- rbind(ARB_term_defined, ARB_read2_defined) %>%
  arrange(eid, issue_date)

## First prescription
entry1_ARB <- ARB %>%
  distinct(eid, .keep_all = TRUE)           


## Write files ----
fwrite(ARB, "data/UKB_all_ARB_prescriptions.txt", sep = "\t", col.names = TRUE)
fwrite(entry1_ARB, "data/UKB_first_ARB_prescription.txt", sep = "\t", col.names = TRUE)


#------------------------
# ACEI and ARB users ----
#------------------------

# Utilise ARB data ----
## Add ARB_user column which denotes whether an ACEI user also have a ARB prescription (yes, 1; no, 0)
ACEI_ARB_user <- entry1_ACEI %>%
  mutate(ARB_user = case_when(eid %in% entry1_ARB$eid ~1,
                              T ~ 0)) %>%
  filter(ARB_user == 1) %>%          # Select for ACEI and ARB users
  select(eid)

## Combine ACEI and ARB data and filter for individuals with both, only
ACEI <- ACEI %>%
  mutate(drug_type = "ACEI") %>%
  select(-ACEI)

ARB <- ARB %>%
  mutate(drug_type = "ARB") %>%
  select(-ARB)

ACEI_ARB <- rbind(ACEI, ARB) %>%
  arrange(eid, issue_date) %>%
  filter(eid %in% ACEI_ARB_user$eid)

## Write file
fwrite(ACEI_ARB, "data/UKB_all_ACEI_ARB_prescriptions.txt", sep = "\t", col.names = TRUE)
