library(tidyverse)
library(rjson)

load_hvc_json = function(fname) {
  d = fromJSON(file=fname)
  feats = d$features
  feats_df = data.frame(do.call(rbind, feats))
  colnames(feats_df) = make.names(d$feature_list)
  feats_df = feats_df %>% mutate(labels = d$pred_labels,
                                 prob_ent = d$pred_probs_ent,
                                 onsets = d$onsets_s,
                                 offsets = d$offsets_s,
                                 wav = d$songfiles[d$songfile_IDs+1], # +1 for python 0-index
                                 mat = paste0(wav, ".not.mat"),
                                 bird = d$bird_ID)
  
  return(feats_df)
}

load_hvc_db = function(db_fname, collect=T) {
  print(db_fname)
  #json_fname = gsub("\\.db", "\\.json", db_fname)
  #d = fromJSON(file=json_fname)
  if (collect) {
  feats_df = collect(tbl(src_sqlite(db_fname), "data"), n=Inf)
  } else{
    feats_df = tbl(src_sqlite(db_fname), "data")
  }
  #feats_df = feats_df %>% mutate(#labels = d$pred_labels,
                                 #onsets = d$onsets_s,
                                 #offsets = d$offsets_s,
                                 #wav = d$songfiles[d$songfile_IDs+1], # +1 for python 0-index
                                 #mat = paste0(wav, ".not.mat"),
  #                               bird = d$bird_ID)
  #
  return(feats_df)
}