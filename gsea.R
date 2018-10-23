library(tidyverse)
  
load_edb_data = function(path) {
  

## Load .edb data
path.edb <- list.files(path = file.path(path, "edb"),
                       pattern = ".edb$", full.names = TRUE)
gsea.edb <- read.delim(file = path.edb,
                       header = FALSE, stringsAsFactors = FALSE)
gsea.edb <- unlist(gsea.edb)
gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
gsea.metric <- unlist(strsplit(gsea.metric, " "))
gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
gsea.metric <- gsub("METRIC=", "", gsea.metric)
gsea.edb <- gsea.edb[grep("<DTG", gsea.edb)]

# # Select the right gene set
# if (length(gsea.edb) == 0) {
#   stop(paste("The gene set name was not found, please provide",
#              "a correct name"))
# }
# 
# if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)) > 1) {
#   warning(paste("More than 1 gene set matched the gene.set",
#                 "argument; the first match is plotted"))
# }
# gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)[1]]

# Get template name
#gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)

gsea.form = map(gsea.edb, function(gsea.edb1) {
  #gsea.edb1 <- unlist(strsplit(gsea.edb1, "\t"))
  #gsea.template <- gsea.edb1[1]
  #gsea.template = gsub("(?<=TEMPLATE\\=)(.*)(?= GENE)", "\\2",  gsea.edb1, perl=T)
  gsea.rnk = str_extract(gsea.edb1, pattern = "(?<=RANKED_LIST=)(.*)(?= TEMPLATE)")
  gsea.template = str_extract(gsea.edb1, pattern = "(?<=TEMPLATE\\=)(.*)(?= GENE)")
  # Get gene set name
  #gsea.gene.set <- gsea.edb1[2]
  #gsea.gene.set <- gsub("GENESET=gene_sets.gmt#(.*) ES", "\\1", gsea.gene.set)
  gsea.gene.set = str_extract(gsea.edb1, "(?<=GENESET\\=gene_sets.gmt#)(.*)(?= ES=)")
  # Get enrichment score
  #gsea.enrichment.score <- gsea.edb1[3]
  #gsea.enrichment.score <- gsub("ES=", "", gsea.enrichment.score)
  gsea.enrichment.score = as.numeric(str_extract(gsea.edb1, "(?<=ES\\=)(.*)(?= NES=)"))
  # Get gene set name
  #gsea.normalized.enrichment.score <- gsea.edb1[4]
  #gsea.normalized.enrichment.score <- gsub("NES=", "",
  #                                         gsea.normalized.enrichment.score)
  gsea.normalized.enrichment.score = as.numeric(str_extract(gsea.edb1, "(?<=NES\\=)(.*)(?= NP=)"))
  # Get nominal p-value
  #gsea.p.value <- gsea.edb1[5]
  #gsea.p.value <- gsub("NP=", "", gsea.p.value)
  #gsea.p.value <- as.numeric(gsea.p.value)
  gsea.p.value = as.numeric(str_extract(gsea.edb1, "(?<=NP=)(.*)(?= FDR=)"))
  # Get FDR
  #gsea.fdr <- gsea.edb1[6]
  #gsea.fdr <- gsub("FDR=", "", gsea.fdr)
  #gsea.fdr <- as.numeric(gsea.fdr)
  gsea.fdr = as.numeric(str_extract(gsea.edb1, "(?<=FDR=)(.*)(?= FWER=)"))
  gsea_df = data.frame(rnk = gsea.rnk,
                       template = gsea.template,
                       gene.set = gsea.gene.set,
                       ES = gsea.enrichment.score,
                       NES = gsea.normalized.enrichment.score,
                       pvalue = gsea.p.value,
                       fdr = gsea.fdr,
                       stringsAsFactors = F)
  return(gsea_df)
}) %>% bind_rows()
return(gsea.form)
}