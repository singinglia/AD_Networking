library(biomaRt)
library(dplyr)
# library(edgeR) #scoped on cpm

DO_CPM <- TRUE

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

todenCounts = read.delim("D:/Classes/AD Project/toden_counts.txt", row.names = 1, sep = "\t")

if (DO_CPM){
  todenCounts = edgeR::cpm(todenCounts)
  
  todenCounts = round(todenCounts*10000, 0)
  
  todenCounts = data.frame(todenCounts)
}

geneRef = getBM(attributes = c("ensembl_gene_id","external_gene_name", "hgnc_symbol"), 
                filters = 'ensembl_gene_id', 
                values = row.names(todenCounts), 
                mart = ensembl)
# head(geneRef)

names = geneRef %>% select(ensembl_gene_id, hgnc_symbol) %>% filter(hgnc_symbol != "")

todenCounts = todenCounts %>% tibble::rownames_to_column("ens_id")

gather = names %>% inner_join(todenCounts, by = c("ensembl_gene_id" = "ens_id"))

gather = gather %>% select(-ensembl_gene_id) %>%  
  group_by(hgnc_symbol) %>% 
  summarise_if(is.numeric, sum, na.rm = TRUE)

genes = gather$hgnc_symbol
row.names(gather) <- gather$hgnc_symbol

gather2 = gather %>% select(-hgnc_symbol)

gather2 = data.frame(t(gather2))

names(gather2) <- genes

meta = readr::read_csv("D:/Classes/AD Project/toden_metadata.csv")

final_full = gather2 %>% tibble::rownames_to_column("Run") %>% left_join(meta, by="Run")
editMeta = meta %>% mutate(hasAD = ifelse(Disease == "NCI", FALSE, TRUE)) %>%
  select(Run, hasAD)
final_state_only = gather2 %>% tibble::rownames_to_column("Run") %>% left_join(editMeta, by="Run")

if (DO_CPM){
  write_name = "D:/Github/AD_Networking/data/cleaned_tolden_data_cpm_times10000.csv"
} else {
  write_name = "D:/Github/AD_Networking/data/cleaned_tolden_data_counts.csv"
}

readr::write_csv(final_state_only, write_name)

