library(stringr)
library(dplyr)

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
  read.delim(fname, colClasses=colClasses, header=F)
  
}

id_split = function(x) {
  fr = str_split(gtf[1,9,with=F], ";")[[1]]
  fr = sapply(fr, function(x) str_split(str_trim(x), " ")[[1]][1])
  fr = fr[1:(length(fr)-1)] 
}


parse_gtf = function(gtf, ncores=6) {
  std_cols = c("chrom", "build" ,"feature", "start", "end", "flag", "strand", "flag2")
  colnames(gtf)[1:8] = std_cols
  fr = str_split(gtf[1,9,with=F], ";")[[1]]
  fr = sapply(fr, function(x) str_split(str_trim(x), " ")[[1]][1])
  #fr = fr[1:(length(fr)-1)]
  gtf2 = gtf[,c("A","B", "C", "D") := tstrsplit(V9, ";", fixed=T)]
  id = str_split(gtf[,9, with=F], ";")
  id2 = mclapply(id, function(x) str_split(str_trim(x[1:(length(x)-1)]), " "), mc.cores=ncores)
  cols = unlist(lapply(id2[[1]], function(x) x[1]))
  vals = mclapply(id2, function(x) unlist(lapply(x, function(y) y[2])), mc.cores=ncores)
  vals1 = do.call("rbind", vals)
  #vals1[,ncol(vals1)] = str_replace(vals1[,ncol(vals1)], ";", "")
  
  colnames(vals1)  = cols
  
  return( data.frame(gtf[,1:8], vals1) )
}

write_parsed_gtf = function(gtf, fname) {
  cnames = colnames(gtf)[9:ncol(gtf)]
  cnames_l = length(cnames)
  id = apply(gtf, 1, function(d) {
    tmp = sapply(1:cnames_l, function(x) {
      paste(cnames[x], " \"", d[x+8], "\"", sep="")
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

aggregate_to_gene = function(gtf_p) {
  gtf_p = gtf_p %>% group_by(gene_id) %>% do({
    min_exon = which.min(.$start)
    max_exon = which.max(.$end)
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
  })
}



