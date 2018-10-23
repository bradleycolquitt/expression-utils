library(data.table)
library(reshape2)
library(tidyr)
library(dplyr)
source("~/src/seqAnalysis/R/boot.R")

round_to_arbitrary = function(numbers, arbitrary.numbers, range=1) {
  arbitrary.numbers<-sort(arbitrary.numbers)          # need them sorted
  range <- range*1.000001                             # avoid rounding issues
  nearest <- findInterval(numbers, arbitrary.numbers - range) # index of nearest
  nearest <- c(-Inf, arbitrary.numbers)[nearest + 1]  # value of nearest
  diff <- numbers - nearest                           # compute errors
  snap <- diff <= range                               # only snap near numbers
  numbers[snap] <- nearest[snap]                      # snap values to nearest
  numbers
}


intersect_data_anno = function(data, anno, wsize) {
  data1 = lapply(data, function(d) {
    data.frame(pos=1:length(d) * wsize, d)
  })
 # data1 = bind_rows(data, .id="chrom")
  anno = as.data.frame(fread(anno, header=F))
  colnames(anno) = c("chrom", "start", "end", "name", "score", "strand")
  anno = anno %>% 
    group_by(chrom) %>% 
    mutate(start_break = round_to_arbitrary(start, data1[[chrom[1]]]$pos, range=wsize/2))
  anno
}


load_annotation_profile = function(fname) {
  data = as.data.frame(fread(fname, header=F))
  colnames(data) = c("chr", "start", "end", "name", "position", "strand", "value")
  data
}

group_by_position = function(data, func="mean") {
 # data = read.delim(fname, header=F)
  params = NULL
  if (func=="mean") {
    params = list(trim=.001, na.rm=T)
  }
  #datac = acast(data, position~name, value.var="value", fun.aggregate = mean)
  #datac = data %>% spread(position, name, value)
  data1 = data %>% group_by(position) %>% summarize(value=do.call(func, c(list(value), params)))
}