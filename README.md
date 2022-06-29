# Genome-wide association study of ACE inhibitor-induced cough

This repository contains code to define ACE inhibitor (ACEI)-induced cough in primary care data linked to UK Biobank participants.

## Repository structure

**./code_lists** 

Includes .csv files containing:

* BNF, dm+d and Read v2 codes for ACEI and angiotensin-II receptor blocker (ARB) prescription.  
* Derived search terms (drug names) for ACEI and ARB prescription.
* Read v2 and Read v3 codes for ACEI-related ADR/allergy.  
* Read v2 and Read v3 codes ACEI-related side effects, including cough and angioedema.

**./phenotyping**   

* **/1_event_extraction.R** - R script to extract ACEI and ARB prescription events from primary care data, using a combination of Read v2 codes and search terms, depending on data provider.  

* **/2_record_length.R** - R script to determine duration of complete data i.e. time between date of first ACEI prescription (*index date*) and end of records (defined as date of record extraction, or date of death, whichever occurred first).  

* **/3_switch_phenotyping.R** - R script to define cases and controls for ACEI-induced cough phenotype. Cases are individuals who switched from an ACEI to an ARB within 12 months of *index date*, and did not receive a further ACEI prescription in the following year after switch date. Controls are individuals who received at least one further ACEI prescription and no ARB prescriptions within 12 months of index date. Preliminary exclusions included individuals who; (1) had less than 2 years of complete data, or (2) received an ARB prescription before, or on the same date, as *index date*. Cases are coded as 1, controls are coded as 0.
