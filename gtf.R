#library(stringr)
#library(dplyr)
library(parallel)
library(tidyverse)
select = dplyr::select
map = purrr::map

intron_ex = function(d) {
  out = vector("numeric", nrow(d)-1)
  for (i in 1:length(out)) {
    out[i] = d[i+1, 3] - d[i, 2] 
  }
  return(out)
}

read_gtf = function(fname) {
  colClasses = c("character", 
                 "character", 
                 "character", 
                 "numeric", 
                 "numeric", 
                 "character", 
                 "character", 
                 "character", 
                 "character")
  read_delim(fname, delim="\t", col_names=F)
  
}

id_split = function(x) {
  fr = str_split(gtf[1,9,with=F], ";")[[1]]
  fr = sapply(fr, function(x) str_split(str_trim(x), " ")[[1]][1])
  fr = fr[1:(length(fr)-1)] 
}


parse_gtf = function(gtf, ncores=6, adjust=TRUE) {
  std_cols = c("chrom", "build" ,"feature", "start", "end", "flag", "strand", "flag2", "attrs")
  colnames(gtf) = std_cols
  gtf1 = gtf %>% mutate(a = map(attrs, function(d) strsplit(d, ";"))) %>%
    mutate(b = map(a, function(d) map(d, trimws))) %>% 
    mutate(c = map(b, function(d) map(d, function(d1) strsplit(d1, " ")))) %>%
    mutate(d = map(c, function(x)  {
      ns = map_chr(x[[1]], function(x1) x1[[1]][[1]])
      x1 = set_names(x[[1]], ns)
      x1 = map(x1, function(x2) x2[2])
      x1 = map(x1, function(x2) gsub("\"", "", x2))
      x1
    })) %>%
    mutate(d1 = map(d, bind_rows)) %>% 
    unnest(d1)
  gtf1 %>% select(-a, -b, -c, -d)
}

write_parsed_gtf = function(gtf, fname, character_cols=NULL, cols_order=NULL) {
  gtf = as.data.frame(gtf)
  cnames = colnames(gtf)[9:ncol(gtf)]
  cnames_l = length(cnames)
  if (is.null(character_cols)) character_cols = 1:cnames_l
  
  if (is.null(cols_order)) cols_order = 1:(ncol(gtf)-8)
  gtf[,9:ncol(gtf)] = gtf[, 8+cols_order]
  
  id = apply(gtf, 1, function(d) {
    tmp = sapply(1:cnames_l, function(x) {
      if (x %in% character_cols) {
        paste(cnames[x], " \"", d[x+8], "\"", sep="")
      } else {
        paste(cnames[x], d[x+8], sep=" ") 
      }
    })
    paste(tmp, collapse="; ")
    })
  write.table(cbind(gtf[,1:8], id), file=fname, quote=F, sep="\t", row.names=F, col.names=F)
}

parse_gff = function(gff) {
  std_cols = c("chrom", "build" ,"feature", "start", "end", "flag", "strand", "flag2")
  colnames(gff)[1:8] = std_cols
  id = str_split(gff[,9], ";")
  id2 = lapply(id, function(x) str_split(x, "="))
  cols = unlist(lapply(id2[[1]], function(x) x[1]))
  vals = lapply(id2, function(x) unlist(lapply(x, function(y) y[2])))
  vals1 = do.call("rbind", vals)
  vals1[,ncol(vals1)] = str_replace(vals1[,ncol(vals1)], ";", "")
  colnames(vals1)  = cols
  return( data.frame(gff[,1:8], vals1) )
}

gtf_stats = function(gtf_p) {
  #gtf_p = parse_gtf(gtf)
  # transcript number
  tn = nrow(gtf_p[gtf_p[,3]=="transcript",])
  
  # exon distribution
  en = gtf_p %>% group_by(transcript_id) %>% summarize(exon_num = n()) 
  
  # number of homologs
  hom = apply(gtf, 1, function(x) grep("Similar", x[9]))
  
  # transcript length distribution
  tl = 0
  if (tn > 0) {
    tl = gtf_p %>% filter(feature=="transcript") %>% mutate(tlen = end-start)  
  } else {
    tl = gtf_p %>% group_by(transcript_id) %>% summarize(tlen = max(end) - min(start))
  }
  
  df = merge(en, tl)
  return( df)
}

gtf_stats_maker = function(gtf_p) {
  print("Parsing gtf...")
  #gtf_p = parse_gtf(gtf)
  
  # transcript number
  #tn = length(unique(gtf_p[,9]))
  
  # exon distribution
  gtf_p %>% group_by(transcript_id, gene_name) %>% summarize(exon_num = sum(feature=="exon"), 
                                                             cds_num = sum(feature=="CDS"),
                                                             length = max(end) - min(start)) 
  
  # number of homologs
  #hom = apply(gtf, 1, function(x) grep("Similar", x[9]))
  
  # transcript length distribution
  #tl = 0
  #if (tn > 0) {
  #  tl = gtf_p %>% filter(feature=="transcript") %>% mutate(tlen = end-start)  
  #} else {
  #  tl = gtf_p %>% group_by(transcript_id) %>% summarize(tlen = max(end) - min(start))
  #}
  
  #df = merge(en)
  #return( df)
}

aggregate_to_gene = function(gtf_p, gene_name_present=T) {
  gtf_p = gtf_p %>% group_by(gene_id) %>% do({
    min_exon = which.min(.$start)
    max_exon = which.max(.$end)
    if (gene_name_present) {
    return(data.frame(chrom=.$chrom[1], 
               build=.$build[1], 
               feature="gene", 
               start=.$start[min_exon], 
               end=.$end[max_exon], 
               flag=.$flag[1], 
               strand=.$strand[1], 
               flag2=.$flag2[1],
               transcript_id=.$transcript_id[1], 
               gene_id=.$gene_id[1], 
               gene_name=.$gene_name[1]))
    } else {
      return(data.frame(chrom=.$chrom[1], 
                        build=.$build[1], 
                        feature="gene", 
                        start=.$start[min_exon], 
                        end=.$end[max_exon], 
                        flag=.$flag[1], 
                        strand=.$strand[1], 
                        flag2=.$flag2[1],
                        transcript_id=.$transcript_id[1], 
                        gene_id=.$gene_id[1]))
    }
  })
}

parse_ncbi_headers = function(fname, species) {
  a = read.delim(fname, header=F)
  a1 = str_split(a[,1], " ")
  a1a = unlist(lapply(a1, function(x) str_replace(x[1], ">", "")))
  a1a1 = str_split(a1a, "\\|")
  gi = unlist(lapply(a1a1, function(x) x[2]))
  ref = unlist(lapply(a1a1, function(x) x[4]))
  
  a2 = str_split(a[,1], paste(" \\[", species, sep=""))
  a2a = unlist(lapply(a2, function(x) x[1]))
  a2a1 = str_split(a2a, "\\| ")
  a2a1a = unlist(lapply(a2a1, function(x) x[2]))
  a2a1a1 = str_split(a2a1a, " ")
  a2a1a1 = unlist(lapply(a2a1a1, paste, collapse="_"))
  
  return(data.frame(id=a1a, name=a2a1a1, db="ncbi", species=species))
}

#' species is of form CHICK
parse_uniprot_headers = function(fname, species=NULL) {
  
  a = scan(fname, what="character", sep = "\n")
  a1 = str_split(a, " ")
  a1a = unlist(lapply(a1, function(x) str_replace(x[1], ">", "")))
  
  gn_ind = unlist(lapply(a1, function(x) grep("\\=", x)[1]))
  product = unlist(lapply(1:length(a1), function(i) {
    paste(a1[[i]][2:((gn_ind[i])-1)], collapse=" ")
  }))
  
  a1a1 = str_split(a1a, "\\|")
  id = unlist(lapply(a1a1, function(x) x[2]))
  if (is.null(species)) {
    a1a2 = str_split(a1a, "_")
    species = unlist(lapply(a1a2, function(x) x[2]))
  }
  
  # dummy GN for those without
  ind = grep("GN\\=", a)
  indn = 1:length(a)
  indn = indn[-ind]
  b = id[indn]
  dummy_gn = paste("GN=", b, sep="")
  a2 = lapply(1:length(indn), function(i) c(a1[[indn[i]]], dummy_gn[i]))

  a3 = vector(mode = "character", length = length(a1))
  a3[ind] = a1[ind]
  a3[-ind] = a2
  gn = unlist(lapply(a3, function(x) x[grep("GN\\=", x)]))
  gn1 = str_replace(gn, "GN\\=", "")
  df = data.frame(id=a1a, name=gn1, db="uniprot", product=product, species=species)
  #df = data.frame(apply(df, 2, as.character))
  return(df)
}

form_uniprot_headers = function(df) {
  formed = df %>% 
    mutate(form = sprintf("%s %s OS=%s GN=%s", 
                          id, 
                          product_std, 
                          species, 
                          name))
  return(formed$form)
  
}

