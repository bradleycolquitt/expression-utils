load_umi_tools_counts = function(umi_dir,
                                 bc_file = "~/lib/barcodes/bc48.txt",
                                 genes_file="/media/data2/assembly/lonStrDom1/ncbi/lonStrDom1.ncbi.apollo6/merge_inter_apollo_processed_mt.gtf") {
  count_files = list.files(umi_dir, pattern="counts.tsv.gz", recursive = T, full.names = T)
  count_prefix = map_chr(count_files, function(x) {
    x1 = unlist(str_split(x, "/"))
    x2 = x1[(length(x1)-1)]
    x3 = str_extract( x2, "(?<=pool)(.*)")
    x3
  }) 
  names(count_files) = count_prefix
  data = map(count_files, ~read_delim(.x, delim="\t", col_names=c("transcript_id", "cell", "unique_umi"), skip=1, )) %>% bind_rows(.id="pool") %>% mutate(pool = as.numeric(pool))
  
  bc = read_delim(bc_file, delim="\t", col_names=c("bc", "cell"))
  bc = bc %>% mutate(barcode = sub("bc", "", bc))
  genes = import(genes_file)
  genes_df = as.data.frame(genes[,c("gene_id", "gene_name", "transcript_id")]) %>% distinct(gene_id, transcript_id, gene_name)
  data = data %>% left_join(bc) %>% left_join(genes_df) %>%
    select(-gene_id) %>% 
    rename(gene_id = gene_name) %>%
    filter(!is.na(barcode))
  
  return(data)
}