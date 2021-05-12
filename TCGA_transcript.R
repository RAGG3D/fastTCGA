install.packages("data.table")

options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(tidyverse)
library(tidybulk)
library(data.table)

#Aggregate samples and do TMM normalization
TCGA_transcript <- function(cancer){
  setwd(paste0("D:/Isoforms/TCGA-", cancer, "/Transcript/"))
  as_tibble(
    data.table::as.data.table(
      map_dfr(as.list(list.files(paste0("D:/Isoforms/TCGA-", cancer, "/Transcript/"))), 
              function(i){
                data.table::fread(i) %>%
                  mutate(sample = i)
              }) %>%
        tidybulk::rename(raw_count = `V2`) %>%
        mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
        inner_join(toTable(org.Hs.egENSEMBL)) %>%
        inner_join(toTable(org.Hs.egSYMBOL)) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0("D:/Isoforms/TCGA-", cancer, "/gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%  #Sum the raw_counts of duplicated rows
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}

#Combine clinical data
clinical_combine <- function(x, cancer) {
  x %>% 
    inner_join(read.csv(paste0("D:/Isoforms/TCGA-",cancer, "/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))}

