install.packages("data.table")

options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(tidyverse)
library(tidybulk)
library(survminer)
library(survival)
library(gapminder)
library(foreach)
library(cowplot)
library(ggsci)
library(GGally)
library(gridExtra)
library(grid)
library(reshape)
library(Hmisc)
library(data.table)


g2s=toTable(org.Hs.egSYMBOL)
g2t=toTable(org.Hs.egENSEMBLTRANS)
g2e=toTable(org.Hs.egENSEMBL)

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
        inner_join(g2e) %>%
        inner_join(g2s) %>%
        dplyr::select(sample, symbol, raw_count) %>%
        mutate(sample = gsub("counts.*", "counts.gz", sample)) %>%
        inner_join(
          read_csv(paste0("D:/Isoforms/TCGA-", cancer, "/gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
        mutate(sample = `Case ID`) %>%
        dplyr::select(sample, symbol, raw_count)
    )[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%
    tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol)
}



clinical_combine <- function(x, cancer) {
  x %>% 
    inner_join(read.csv(paste0("D:/Isoforms/TCGA-",cancer, "/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))}
