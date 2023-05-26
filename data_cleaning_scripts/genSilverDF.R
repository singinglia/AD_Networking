library(biomaRt)
library(dplyr)
# library(edgeR) #scoped on cpm

DO_CPM <- TRUE

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

rnaCounts = read.delim("D:/Classes/AD Project/silver_seq_counts.txt", row.names = 1, sep = "\t")

if (DO_CPM){
  rnaCounts = edgeR::cpm(rnaCounts)
  
  rnaCounts = round(rnaCounts*10000, 0)
  
  rnaCounts = data.frame(rnaCounts)
}

geneRef = getBM(attributes = c("ensembl_gene_id","external_gene_name", "hgnc_symbol"), 
                filters = 'ensembl_gene_id', 
                values = row.names(rnaCounts), 
                mart = ensembl)
# head(geneRef)

names = geneRef %>% select(ensembl_gene_id, hgnc_symbol) %>% filter(hgnc_symbol != "")

rnaCounts = rnaCounts %>% tibble::rownames_to_column("ens_id")

gather = names %>% inner_join(rnaCounts, by = c("ensembl_gene_id" = "ens_id"))

gather = gather %>% select(-ensembl_gene_id) %>%  
  group_by(hgnc_symbol) %>% 
  summarise_if(is.numeric, sum, na.rm = TRUE)

genes = gather$hgnc_symbol
row.names(gather) <- gather$hgnc_symbol

gather2 = gather %>% select(-hgnc_symbol)

gather2 = data.frame(t(gather2))

names(gather2) <- genes

meta = readr::read_csv("D:/Classes/AD Project/silver_seq_metadata.csv")

final_full = gather2 %>% tibble::rownames_to_column("Run") %>% 
  left_join(meta, by=c("Run" = "sample_id_alias"))
editMeta = meta %>% mutate(hasAD = ifelse(donor_group == "N", FALSE, TRUE)) %>%
  select(sample_id_alias, hasAD)
final_state_only = gather2 %>% tibble::rownames_to_column("Run") %>% 
  left_join(editMeta, by="Run")

if (DO_CPM){
  write_name = "D:/Github/AD_Networking/data/cleaned_silver_data_cpm_times10000.csv"
} else {
  write_name = "D:/Github/AD_Networking/data/cleaned_silver_data_counts.csv"
}

readr::write_csv(final_state_only, write_name)

