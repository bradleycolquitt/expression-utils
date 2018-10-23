library(tidyverse)
library(googlesheets)
library(dplyr)
library(RMySQL)
library(lubridate)
library(stringr)
library(lazyeval)
tmp_gs_dir = "/media/data2/rstudio/birds/old_googlesheets"

#' http://stackoverflow.com/questions/5090082/handling-field-types-in-database-interaction-with-r
dbGetTypes = function(con,table){
  statement <- paste("DESCRIBE ",table,sep="")
  desc <- dbGetQuery(con=con,statement)[,1:2]
  return(desc)
}

dbReadMap <- function(con, table) {
  # strip row_names if exists because it's an attribute and not real column
  # otherweise it causes problems with the row count if the table has a row_names col
  if(length(grep(pattern="row_names",x=desc)) != 0){
    x <- grep(pattern="row_names",x=desc)
    desc <- desc[-x,]
  }
  
  
  
  # replace length output in brackets that is returned by describe
  desc[,2] <- gsub("[^a-z]","",desc[,2])
  
  # building a dictionary 
  fieldtypes <- c("int","tinyint","bigint","float","double","date","character","varchar","text")
  rclasses <- c("as.numeric","as.numeric","as.numeric","as.numeric","as.numeric","as.Date","as.character","as.character","as.character") 
  fieldtype_to_rclass = cbind(fieldtypes,rclasses)
  
  map <- merge(fieldtype_to_rclass,desc,by.x="fieldtypes",by.y="Type")
  map$rclasses <- as.character(map$rclasses)
  #get data
  res <- dbReadTable(con=con,table)
  
  
  
  i=1
  for(i in 1:length(map$rclasses)) {
    cvn <- call(map$rclasses[i],res[,map$Field[i]])
    res[map$Field[i]] <- eval(cvn)
  }
  
  
  return(res)
}

## GOOGLE / MYSQL -------------------------------------------------------------------------------------------------------------------
insert_and_update = function(con, name, value) {
  value = as.data.frame(value)
  types = dbGetTypes(con, name)
  
  # Ensure that dates are in ISO format
  date_ind = which(types$Type=="date")
  for (i in date_ind) {
    inds = grep("/", value[,i])
    if (length(inds)>0)
      value[inds,i] = as.character(as.Date(value[inds,i], format="%m/%d/%Y"))
  }
  
  # Add quotes around character, date, and set fields
  character_cols = grep("char|date|set|enum", types$Type)
  for (i in character_cols) {
      value[,i] = paste("'", value[,i], "'", sep="")
  }
  
  #res = vector("logical", nrow(value))
  cnames = colnames(value)
  for (i in 1:nrow(value)) {
    print(value[i,])
    print(i)
    #print(paste(value[i,], sep=","))
    nna_ind = !is.na(value[i,])
    sql = paste("INSERT INTO", name, "(", paste(colnames(value)[nna_ind], collapse=","), 
                ") VALUES (", paste(value[i,nna_ind], collapse=","), ") ON DUPLICATE KEY UPDATE",
                paste(paste(cnames[nna_ind], "= VALUES (", cnames[nna_ind], ")", sep=" "), collapse=","), sep=" ")
    #print(sql)
    dbSendQuery(con, sql)
  }
  
}

get_google_data = function(gs_name, table_name=1) {
  ss = gs_title(gs_name)
  data = ss %>% gs_read(table_name)
  return(data)
}

#' Reads data from Google Sheets and inserts into mysql database
#' Note: will update records with identical primary keys
update_mysql_from_google = function(db_name="sample_db", gs_name, table_name) {
  data = get_google_data(gs_name, table_name)
  print(dim(data))
  print(colnames(data))
  mdb = src_mysql(db_name, username = 'brad', password = 'Eu23ler1')
  types = dbGetTypes(mdb$con, table_name)
  types_list = as.list(types$Type)
  names(types_list) = types$Field
  

  insert_and_update(mdb$con, table_name, value=data)
}

update_mysql_from_google2 = function(db_name="sample_db", expt, table_name) {
  gs_name = paste(expt, table_name, "current", sep="_")
  data = get_google_data(gs_name)
  print(dim(data))
  print(colnames(data))
  mdb = src_mysql(db_name, username = 'brad', password = 'Eu23ler1')
  types = dbGetTypes(mdb$con, table_name)
  types_list = as.list(types$Type)
  names(types_list) = types$Field
  data = data.frame(data)
  
  # not_all_na = function(x) {
  #   res = is.na(x)
  #   sum(res) < length(x)
  # }
  # not_na = function(x) {
  #   res = !is.na(x)
  #   res
  # }
  # data1 = data.frame(data) %>% filter_all(all_vars(not_na(.)))
  insert_and_update(mdb$con, table_name, value=data)
}


update_google_from_mysql_by_expt = function(db_name="sample_db", expt, table_name=NULL) {
  mdb = src_mysql(db_name, username = 'brad', password = 'Eu23ler1')
  table_names = c("sectioning", "lcm", "purification", "library_prep", "library_pooled", "library_sequenced", "library", "sequencing")
  if (!is.null(table_name))
    table_names = table_name
  
  comb = load_sample_info(db_name, expt)
  for (tab in table_names) {
    type = dbGetTypes(mdb$con, tab)
    id_name = paste("id", tab, sep="")
    d = comb %>% dplyr::select_(.dots=type$Field) %>% 
      distinct() %>% 
      arrange_(.dots=id_name) %>% 
      dplyr::filter_(.dots=interp(paste("!is.na(", id_name, ")", sep="")))

    export_data_to_google(expt, tab, d, expt=expt)
    #export_data_to_google2(paste(expt, tab, sep="_"), d)
  }
}

update_google_from_mysql_by_expt2 = function(db_name="sample_db", expt, table_name=NULL) {
  mdb = src_mysql(db_name, username = 'brad', password = 'Eu23ler1')
  table_names = c("sectioning", "lcm", "purification", "library_prep", "library_pooled", "library_sequenced", "library", "sequencing")

  comb = load_sample_info(db_name, expt)
  if (!is.null(table_name)) 
    table_names = table_name
  
  for (tab in table_names) {
    type = dbGetTypes(mdb$con, tab)
    id_name = paste("id", tab, sep="")
    d = comb %>% dplyr::select_(.dots=type$Field) %>% 
      distinct() %>% 
      arrange_(.dots=id_name) %>% 
      dplyr::filter_(.dots=lazyeval::interp(paste("!is.na(", id_name, ")", sep="")))
    
    export_data_to_google2(paste(expt, tab, sep="_"), d)
  }
  
}

update_google_from_mysql = function(db_name="sample_db", gs_name, table_name) {
  mdb = src_mysql(db_name)
  data = tbl(mdb, "view1")
  export_data_to_google(gs_name, table_name, data)
}

delete_mysql_by_expt = function(db_name="sample_db", expt, table_name) {
  mdb = src_mysql(db_name, username = 'brad', password = 'Eu23ler1')
  info = load_sample_info(expt=expt)
  filter_column = paste("id", table_name, sep="")
  filter_tags = as.data.frame(unique(na.omit(info[,filter_column])[,1]))
  filter_tags = str_replace(filter_tags[,1], "\\\n", "")
  sql = paste("DELETE FROM", table_name, "WHERE", filter_column, "IN (", paste(filter_tags, collapse=","), ")", sep=" ")
  print(sql)
  dbSendQuery(mdb$con, sql)
}

delete_mysql_by_expt_and_bird = function(db_name="sample_db", expt, table_name, birds) {
  mdb = src_mysql(db_name,username = 'brad', password = 'Eu23ler1')
  info = load_sample_info(expt=expt)
  #expt_tags = unique(info$tags)
  filter_tags = sapply(birds, function(x) paste("'", x, "'", sep=""))
  filter_column = "tags"
  #filter_tags = birds
  sql = paste("DELETE FROM", table_name, "WHERE", filter_column, "IN ", "(", paste(filter_tags, collapse=","), ")", sep=" ")
  print(sql)
  dbSendQuery(mdb$con, sql)
}

delete_mysql_by_bird = function(db_name="sample_db",table_name, birds) {
  mdb = src_mysql(db_name,username = 'brad', password = 'Eu23ler1')
  #info = load_sample_info(expt=expt)
  #expt_tags = unique(info$tags)
  filter_tags = sapply(birds, function(x) paste("'", x, "'", sep=""))
  filter_column = "tags"
  #filter_tags = birds
  sql = paste("DELETE FROM", table_name, "WHERE", filter_column, "IN ", "(", paste(filter_tags, collapse=","), ")", sep=" ")
  print(sql)
  dbSendQuery(mdb$con, sql)
}



export_data_to_google = function(gs_name, table_name, data, expt="generic") {
  gs = gs_title(gs_name)
  gs_ws_info = gs$ws
  orig = NULL

  if (table_name %in% gs_ws_info$ws_title) {
    orig = gs %>% gs_read(table_name)
    saveRDS(orig, file=paste(tmp_gs_dir, paste(gs_name, table_name, Sys.Date(), sep="-"), sep="/"))
    
    tmp_name = paste(table_name, "tmp", sep="_")
    tmp_title = gs_ws_info$ws_title[grep(tmp_name, gs_ws_info$ws_title)]
    if (length(tmp_title)>0) {
      max_tmp_num = str_split(tmp_title, "_")
      max_tmp_num = max(as.numeric(unlist(lapply(max_tmp_num, function(x) x[length(x)]))))
      tmp_name = paste(tmp_name, max_tmp_num+1, sep="_")
    } else {
      tmp_name = paste(tmp_name, 1, sep="_")
    }
    gs = gs_ws_rename(gs, table_name, tmp_name)
  }
  gs = gs_ws_new(gs, ws_title=table_name, input=data)
  return(orig)
}

# individual sheets instead of worksheets within a common sheet
export_data_to_google2 = function(gs_name, data) {
  gs_name_current = sprintf("%s_current", gs_name)
  all_gs = gs_ls()
  if (gs_name_current %in% all_gs$sheet_title) {
    gs = gs_title(gs_name_current)
    gs_ws_info = gs$ws
    orig = NULL
    orig = gs %>% gs_read()
    saveRDS(orig, file=paste(tmp_gs_dir, paste(gs_name, Sys.Date(), sep="-"), sep="/"))
    
    tmp_name = paste(gs_name, "tmp", sep="_")
    
    tmp_title = all_gs$sheet_title[grep(tmp_name,all_gs$sheet_title)]
    if (length(tmp_title)>0) {
      max_tmp_num = str_split(tmp_title, "_")
      max_tmp_num = max(as.numeric(unlist(lapply(max_tmp_num, function(x) x[length(x)]))))
      tmp_name = paste(tmp_name, max_tmp_num+1, sep="_")
    } else {
      tmp_name = paste(tmp_name, 1, sep="_")
    }
    gs_copy(gs, tmp_name)
    
  }
  tmp_file = sprintf("%s.csv", tempfile())
  write_csv(data, tmp_file)
  gs_upload(tmp_file, gs_name_current, overwrite=TRUE)
}

## MYSQL ----------------------------------------------------------------------------------------------------------------------------
load_sample_info = function(db_name="sample_db", expt, include_seq_info=TRUE) {
  select = dplyr::select
  
  
  nname = Sys.info()["nodename"]
  if (nname == "lyre") {
    mdb = src_mysql(db_name, host="localhost", user="brad", pass = "Eu23ler1")
  } else if (nname == "osprey") {
    mdb = src_mysql(db_name, host="lyre.cin.ucsf.edu", user="osprey", pass = "Bigbird1?") 
  }
  experiments = collect(tbl(mdb, "experiments"))
  birds = NULL
  if (db_name=="sample_db") {
    birds = collect(tbl(mdb, "bird"))
  } else if (db_name == "rubenstein") {
    birds = collect(tbl(mdb, "animal"))
  }

  embedded = collect(tbl(mdb, "embedded"))
  
  if (include_seq_info) {
    sectioning = collect(tbl(mdb, "sectioning"))
    lcm = collect(tbl(mdb, "lcm"))
    purification = collect(tbl(mdb, "purification"))
    prep = collect(tbl(mdb, "library_prep"))
    pooled = collect(tbl(mdb, "library_pooled"))
    library = collect(tbl(mdb, "library"))
    libseq = collect(tbl(mdb, "library_sequenced"))
    sequencing = collect(tbl(mdb, "sequencing"))
    
    comb = experiments %>% 
      
      inner_join(embedded, by="experiment") %>%
      inner_join(birds, by="idbirds") %>%
      full_join(sectioning, by="tags") %>% #select(-idsectioning) %>%
      full_join(lcm, by=c("section_date", "tags", "lcm_slide", "slide_section")) %>% #select(-idlcm) %>%
      full_join(purification, by=c("lcm_date", "lcm_plate")) %>% #select(-idpurification) %>%
      full_join(prep, by=c("pure_date", "pure_plate", "lcm_row", "lcm_column")) %>%
      full_join(pooled, by=c("lib_date", "lib_plate", "pool")) %>%
      full_join(library, by=c("idlibrary")) %>%
      full_join(libseq, by=c("idlibrary")) %>%
      full_join(sequencing %>% select(-idlibrary), by=c("idsequencing")) %>%
      filter(experiment==expt)
  } else {
    comb = experiments %>%
      inner_join(embedded, by="experiment") %>%
      inner_join(birds, by="idbirds") %>% 
      filter(experiment==expt)
  }
  colnames(comb) = make.names(colnames(comb))
  return(comb)
}

umi_dirname = "~/data/umi_db"

get_birds_by_experiment = function(db_name="sample_db", expt) {
  embed = tbl(mdb, "embedded")
  birds = tbl(mdb, "bird")
  a = birds %>% left_join(embed, by="idbirds")
  ac = collect(a, n=Inf)
  ac1 = ac %>% filter(experiment==expt)
  ac1t = ac1$tags
  return(ac1t)
}

## SQLITE ---------------------------------------------------------------------------------------------------------------------------
load_umi_db = function(db_names, db_seq, table_name = "collapse_nf", grouping=NULL, collect=TRUE) {
  
  if (length(db_names)> 1) {
    dbs = lapply(db_names, function(d) {
      src_sqlite(paste(umi_dirname, d, sep="/"))
    })
    
    data = lapply(1:length(dbs), function(i) {
      data_tbl = tbl(dbs[[i]], table_name)
      #data_tbl = data_tbl %>% filter(barcode>0)
      
      if (!is.null(grouping)) {
        data_tbl = data_tbl %>% group_by_(.dots=grouping) %>% 
          summarize(unique_umi=sum(unique_umi),
                    total_umi=sum(total_umi))
      }
      d = collect(data_tbl, n=Inf)
      d$idsequencing = db_seq[i]
      d
    })
    data = bind_rows(data)
  } else {
    db =src_sqlite(paste(umi_dirname, db_names, sep="/"))
    
    data_tbl = tbl(db, table_name)
    #data_tbl = data_tbl %>% filter(barcode>0)
    if (!is.null(grouping)) {
      data_tbl = data_tbl %>% group_by_(.dots=grouping) %>% 
        summarize(unique_umi=sum(unique_umi),
                  total_umi=sum(total_umi))
    }
    if (collect) {
    data = collect(data_tbl, n=Inf)
    } else {
      return(data_tbl)
    }
    data$idsequencing = db_seq
  }
  data$index1 = factor(data$index1)
  data$index2 = factor(data$index2)
  data$barcode = factor(data$barcode)
  data$idsequencing = factor(data$idsequencing)
  return(data)
}

load_umi_info = function(info_name, table_name="db_info") {
  db_info = src_sqlite(paste(umi_dirname, info_name, sep="/"))
  info = collect(tbl(db_info, table_name))
  info = info %>% filter(position != "NA")
  info$barcode = factor(info$barcode)
  info$index1 = factor(info$index1)
  info$index2 = factor(info$index2)
  info$idsequencing = factor(info$idsequencing)
  info$id = info$idlibrary_prep
  return(info)
}

sqlite_insert = function(data, db, table_name, replace=F) {
  if (!db_has_table(db$con, table_name)) {
    copy_to(db, data, name=table_name, temporary=F)
  } else {
    if (replace) {
      db_drop_table(db$con, table_name)
      copy_to(db, data, name=table_name, temporary=F)
    } else {
      db_insert_into(db$con, table_name, data) 
    }
  }
}

get_unprocessed_mat_files = function(db_name, table_name, info, replace=F) {
  if (!file.exists(db_name)) {
    db = src_sqlite(db_name, create=T)
    info1 = info
  } else {
    if (!replace) {
      db = src_sqlite(db_name, create=F)
      d01 = tbl(db, table_name)
      mats = d01 %>% dplyr::select(mat) %>% collect() %>% .$mat
      info1 = info %>% filter(!(mat %in% mats))
    } else {
      info1 = info
    }
  }
  return(info1)
}

delete_by_key = function(db, table_name, key, value) {
  if (db_has_table(db$con, table_name)) {
    sql = sprintf("DELETE FROM %s WHERE %s='%s'", table_name, key, value)
    dbSendQuery(db$con, sql)
  } else{
    stop("No table.")
  }
}
replace_by_key = function(db, table_name, data, key) {
  data = ungroup(data)
  value = unique(data[,key])
  if (db_has_table(db$con, table_name)) {
    sql = sprintf("DELETE FROM %s WHERE %s='%s'", table_name, key, value)
    dbSendQuery(db$con, sql)
    db_insert_into(db$con, table_name, data, temporary=F)
  } else {
    copy_to(db, data, table_name, temporary=F)    
  }
}

collapse_umis = function(data) {
  data %>% group_by(umi) %>% summarize(count = sum(count))
}

get_umi_groups = function(umis, thresh, counts) {
  names(counts) = umis
  n_umis = length(umis)
  umi_dist = as.matrix(stringdistmatrix(umis, method="hamming", useNames = T))
  umi_dist[lower.tri(umi_dist)] = 0
  umi_dist_ind = which(umi_dist <= thresh & umi_dist > 0, arr.ind=T)

  
  if (nrow(umi_dist_ind) > 0) {
    #umi_dist_filt = umi_dist_filt[1:(nrow(umi_dist_filt)/2),]
    #umi_dist_filt = umi_dist_filt[seq(1, nrow(umi_dist_filt), 2),]
    for (i in 1:nrow(umi_dist_ind)) {
      if(counts[umi_dist_ind[i,2]] > (counts[umi_dist_ind[i,1]] * 2 - 1)) {
        tmp = umi_dist_ind[i,2]
        umi_dist_ind[i,2] = umi_dist_ind[i,1]
        umi_dist_ind[i,1] = tmp
      }
    } 
    #umi_dist_filt = bind_rows(umi_dist1, umi_dist_filt)
    ig = graph_from_data_frame(umi_dist_ind, directed = T, vertices = 1:n_umis)
    out = data.frame(count = n_umis, count_filt=count_components(ig))
  } else {
    out = data.frame(count = n_umis, count_filt=n_umis)
  }
  return(out)
}

get_umi_groups_batch = function(data, db) {
 
  table_name = "counts_filt"
  
  if (db_has_table(db$con, table_name))
    db_drop_table(db$con, table_name)
  
  
  bcs = data %>% select(bc) %>% distinct(bc) %>% filter(bc>0) %>% collect() %>% unlist() %>% set_names(NULL)
  index1s = data %>% select(i7) %>% distinct(i7) %>% collect() %>% unlist() %>% set_names(NULL)
  
  foreach(barcode=bcs) %do% {
    data1 = data %>% filter(bc==barcode) %>% collect(n=Inf)
    system.time({
      a = foreach(index1 = index1s) %do% {
        res = data1 %>% 
          filter(i7==index1) %>%
          group_by(i7, gene_id) %>% 
          filter(!grepl("NNNN", umi)) %>% 
          do(collapse_umis(.)) %>% 
          do({
            get_umi_groups(.$umi, 1, .$count)
          })
        res$bc = barcode
        db_insert_into(db$con, "counts_filt", res)
        rm(res)
      } %>% set_names(index1s) %>% bind_rows(.id="i7")
    })
  } %>% set_names(bcs) %>% bind_rows(.id="bc")
}