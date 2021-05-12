source("function.R")

#比较foreach循环和map循环
t1 <- proc.time()
setwd("D:/Isoforms/TCGA-LUSC/Transcript/")
foreach(i = list.files("D:/Isoforms/TCGA-LUSC/Transcript/"), .combine = bind_rows) %do% {
  read_table2(i, col_names = F) %>%
    mutate(sample = i)
} 
t_foreach <- proc.time() - t1

t1 <- proc.time()
setwd(paste0("D:/Isoforms/TCGA-LUSC/Transcript/"))
map_dfr(as.list(list.files("D:/Isoforms/TCGA-LUSC/Transcript/")), 
                function(i){
                  read_table2(i, col_names = F) %>%
                    mutate(sample = i)
                })
t_map <- proc.time() - t1

#比较fread和read_table2
t1 <- proc.time()
setwd(paste0("D:/Isoforms/TCGA-LUSC/Transcript/"))
as_tibble(
  data.table::as.data.table(map_dfr(as.list(list.files("D:/Isoforms/TCGA-LUSC/Transcript/")), 
        function(i){
          data.table::fread(i) %>%
            mutate(sample = i)
        }) %>%
  tidybulk::rename(raw_count = `V2`) %>%
  mutate(ensembl_id = gsub("\\..*", "", V1)) %>%
  inner_join(g2e) %>%
  inner_join(g2s) %>%
  dplyr::select(sample, symbol, raw_count) %>%
  mutate(sample = gsub("counts", "counts.gz", sample)) %>%
  inner_join(
    read_csv(paste0("D:/Isoforms/TCGA-LUSC/gdc_sample_sheet.csv")) %>% mutate(sample = `File Name`)) %>%
  mutate(sample = `Case ID`) %>%
  dplyr::select(sample, symbol, raw_count))[, list(raw_count = sum(raw_count)), by = c("sample", "symbol")]) %>%
  tidybulk::scale_abundance(.sample = sample, .abundance = raw_count, .transcript = symbol) 

t_run <- proc.time() - t1

