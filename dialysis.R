source("~/src/songanalysis/song_util.R")

library(googlesheets)
filter = dplyr::filter

load_song_info = function(bird_name,  dates, data_dir="/mnt/bengal_home/song", subdir="songs", batch_file=NULL,  tz="PST", dial_log=NULL, evtaf_log=NULL, surgery_log=NULL, file_ex="wav.not.mat", adjust_time=0) {
  wdirs = sapply(dates, function(x) {
    paste(data_dir, bird_name, x, subdir, sep="/")
  })
  x  = 1
  infos = lapply(wdirs, function(w) {
  #  x = 1
    print(w)
   #batch_file = file.path(w, batch_file)
    infos= load_mat_info(w, file_ex = file_ex, batch_file=batch_file)
#   infos= list(load_mat_info(wdirs[[1]], file_ex = file_ex, batch_file=batch_file))
  })

  info = bind_rows(infos)
  info$bird = bird_name
  
  ## Force timezone
  info = info %>% mutate(mtime = as.POSIXct(mtime, tz=tz),
                         date = as.Date(date))
  ## Parse time info --------------------------------------------------------------------------------
  if (file_ex == "wav.not.mat") {
    info$parsed_time_epoch = parse_fname_for_timestamp(info$wav)
    info$parsed_time = as.POSIXct(parse_timestamp(info$wav), tz = tz)
  } else if (file_ex == "cbin\\.not\\.mat") {
    info$parsed_time = as.POSIXct(parse_timestamp_cbin(info$wav, prefix = sprintf("%s_", bird_name)))
  }

  info = info %>% mutate(time_h = as.numeric(time_h))
  info$date_noon = paste(info$date, "12:00:00", sep=" ")
  info$date_noon = parse_date_time(info$date_noon, orders="ymd hms", tz = tz)
  info$mtime_hour = floor_date(info$mtime, unit="hour")

  ## Combine with drug info ------------------------------------------------------------------------
  if (is.null(dial_log)) {
    dial_log = load_dialysis_gs(bird_name)
  } else {
    dial_log = dial_log %>% filter(log.bird == bird_name)
  }
  if (is.null(surgery_log)) {
    surgery_log = load_event_gs(bird_name)
  } else {
    surgery_log = surgery_log %>% filter(bird == bird_name )
  }
  # dial_log = gs_title("RA cannulation: log") %>% gs_read("dialysis log") %>% filter(bird==bird_name)
  # colnames(dial_log) = paste("log", colnames(dial_log), sep=".")
  # dial_log$drug_time = with(dial_log,paste(log.date, log.time))
  # dial_log$drug_time = as.POSIXct(dial_log$drug_time, format="%Y-%m-%d %I:%M:%S %p")

  if (nrow(dial_log)==0) {
    info$drug_time = NA
  } else {


    info$drug_start = as.POSIXct(cut(info$mtime, breaks=c(dial_log$drug_start, now()), right = T, include.lowest=F))
    #info$drug_time_stop = as.POSIXct(cut(info$mtime, breaks=c(dial_log$drug_time_stop, now()), right = T, include.lowest=F))
    info = left_join(info, dial_log, by=c("drug_start"))
    info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
    

    # info1 = info %>% mutate(dummy=TRUE) %>%
    #   left_join(dial_log %>% mutate(dummy=TRUE)) %>%
    #   filter(parsed_time >= datetime_start, parsed_time <= datetime_stop) %>% select(-dummy)
    # 
    # info2 = info %>% anti_join(info1, by="mat_base")
    # info = bind_rows(info1, info2)
    # info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))

  }
  
  info = info %>% left_join(surgery_log)
  
  # To make compatible across machines
  info = info %>% mutate(mat_base = basename(mat))
  
  if (!is.null(evtaf_log)) {
    info1 = info %>% mutate(dummy=TRUE) %>%
      left_join(evtaf_log %>% mutate(dummy=TRUE) %>% select(-labels)) %>%
    filter(parsed_time >= datetime_start, parsed_time <= datetime_stop) %>% select(-dummy)

    info2 = info %>% anti_join(info1, by="mat_base")
    info = bind_rows(info1, info2)

    #info$evtaf_time = as.POSIXct(cut(info$mtime, breaks=c(evtaf_log$_time, now()), right = T, include.lowest=F))
    #info = left_join(info, dial_log, by="drug_time")
    #info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
  }


  return(info)
}

load_song_info_evtaf = function(bird_name,  dates, data_dir="/mnt/bengal_home/song", subdir="songs", batch_file=NULL,  tz="PST", evtaf_log=NULL, file_ex="cbin.not.mat") {
  wdirs = sapply(dates, function(x) {
    paste(data_dir, bird_name, x, subdir, sep="/")
  })
  x  = 1
  infos = lapply(wdirs, function(w) {
    #  x = 1
    print(w)
    #batch_file = file.path(w, batch_file)
    infos= load_mat_info(w, file_ex = file_ex, batch_file=batch_file)
    #   infos= list(load_mat_info(wdirs[[1]], file_ex = file_ex, batch_file=batch_file))
  })
  
  info = bind_rows(infos)
  info$bird = bird_name
  
  ## Force timezone
  info = info %>% mutate(mtime = as.POSIXct(mtime, tz=tz),
                         date = as.Date(date))
  ## Parse time info --------------------------------------------------------------------------------
  if (file_ex == "wav.not.mat") {
    info$parsed_time_epoch = parse_fname_for_timestamp(info$wav)
    info$parsed_time = as.POSIXct(parse_timestamp(info$wav), tz = tz)
  } else if (file_ex == "cbin\\.not\\.mat") {
    info$parsed_time = as.POSIXct(parse_timestamp_cbin(info$wav, prefix = sprintf("%s_", bird_name)))
  }
  
  info = info %>% mutate(time_h = as.numeric(time_h))
  info$date_noon = paste(info$date, "12:00:00", sep=" ")
  info$date_noon = parse_date_time(info$date_noon, orders="ymd hms", tz = tz)
  info$mtime_hour = floor_date(info$mtime, unit="hour")
  
  
  
  # if (nrow(dial_log)==0) {
  #   info$drug_time = NA
  # } else {
  #   
  #   
  #   info$drug_start = as.POSIXct(cut(info$mtime, breaks=c(dial_log$drug_start, now()), right = T, include.lowest=F))
  #   #info$drug_time_stop = as.POSIXct(cut(info$mtime, breaks=c(dial_log$drug_time_stop, now()), right = T, include.lowest=F))
  #   info = left_join(info, dial_log, by=c("drug_start"))
  #   info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
  #   
  #   
  #   # info1 = info %>% mutate(dummy=TRUE) %>%
  #   #   left_join(dial_log %>% mutate(dummy=TRUE)) %>%
  #   #   filter(parsed_time >= datetime_start, parsed_time <= datetime_stop) %>% select(-dummy)
  #   # 
  #   # info2 = info %>% anti_join(info1, by="mat_base")
  #   # info = bind_rows(info1, info2)
  #   # info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
  #   
  # }
  # 
  # info = info %>% left_join(surgery_log)
  
  # To make compatible across machines
  info = info %>% mutate(mat_base = basename(mat))
  
  if (!is.null(evtaf_log)) {

    info1 = info %>% mutate(dummy=TRUE) %>%
      left_join(evtaf_log %>% mutate(dummy=TRUE) %>% select(-labels)) %>%
      filter(parsed_time >= datetime_start, parsed_time <= datetime_stop) %>% select(-dummy)
    
    info2 = info %>% anti_join(info1, by="mat_base")
    info = bind_rows(info1, info2)
    
    #info$evtaf_time = as.POSIXct(cut(info$mtime, breaks=c(evtaf_log$_time, now()), right = T, include.lowest=F))
    #info = left_join(info, dial_log, by="drug_time")
    #info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
  }
  
  
  return(info)
}

load_song_info_full = function(bird_name, data_dir="/mnt/bengal_home/song", subdir="songs", batch_file=NULL,  tz="PST", dial_log=NULL, evtaf_log=NULL, file_ex="wav.not.mat", adjust_time=0) {
  wdir = paste(data_dir, bird_name, subdir, sep="/")
  info = load_mat_info(wdir, file_ex = file_ex, batch_file=batch_file)
  info$bird = bird_name
  
  ## Force timezone
  info = info %>% mutate(mtime = as.POSIXct(mtime, tz=tz, origin="1970-01-01"),
                         date = as.Date(date))
  ## Parse time info --------------------------------------------------------------------------------
  if (file_ex == "wav.not.mat") {
    info$parsed_time = as.POSIXct(parse_timestamp(info$wav), tz = tz)
  } else if (file_ex == "cbin\\.not\\.mat") {
    info$parsed_time = as.POSIXct(parse_timestamp_cbin(info$wav, prefix = sprintf("%s_", bird_name)))
  }
  
  info = info %>% mutate(time_h = as.numeric(time_h))
  info$date_noon = paste(info$date, "12:00:00", sep=" ")
  info$date_noon = parse_date_time(info$date_noon, orders="ymd hms", tz = tz)
  info$mtime_hour = floor_date(info$mtime, unit="hour")
  
  ## Combine with drug info ------------------------------------------------------------------------
  if (is.null(dial_log)) {
    dial_log = load_dialysis_gs(bird_name)
  }
  # dial_log = gs_title("RA cannulation: log") %>% gs_read("dialysis log") %>% filter(bird==bird_name)
  # colnames(dial_log) = paste("log", colnames(dial_log), sep=".")
  # dial_log$drug_time = with(dial_log,paste(log.date, log.time))
  # dial_log$drug_time = as.POSIXct(dial_log$drug_time, format="%Y-%m-%d %I:%M:%S %p")
  
  if (nrow(dial_log)==0) {
    info$drug_time = NA
  } else {
    
    
    info$drug_start = as.POSIXct(cut(info$mtime, breaks=c(dial_log$drug_start, now()), right = T, include.lowest=F))
    #info$drug_time_stop = as.POSIXct(cut(info$mtime, breaks=c(dial_log$drug_time_stop, now()), right = T, include.lowest=F))
    info = left_join(info, dial_log, by="drug_start")
    info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
    
    # info1 = info %>% mutate(dummy=TRUE) %>%
    #   left_join(dial_log %>% mutate(dummy=TRUE)) %>%
    #   filter(parsed_time >= datetime_start, parsed_time <= datetime_stop) %>% select(-dummy)
    # 
    # info2 = info %>% anti_join(info1, by="mat_base")
    # info = bind_rows(info1, info2)
    # info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
    
  }
  
  # To make compatible across machines
  info = info %>% mutate(mat_base = basename(mat))
  
  if (!is.null(evtaf_log)) {
    info1 = info %>% mutate(dummy=TRUE) %>%
      left_join(evtaf_log %>% mutate(dummy=TRUE) %>% select(-labels)) %>%
      filter(parsed_time >= datetime_start, parsed_time <= datetime_stop) %>% select(-dummy)
    
    info2 = info %>% anti_join(info1, by="mat_base")
    info = bind_rows(info1, info2)
    
    #info$evtaf_time = as.POSIXct(cut(info$mtime, breaks=c(evtaf_log$_time, now()), right = T, include.lowest=F))
    #info = left_join(info, dial_log, by="drug_time")
    #info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
  }
  
  
  return(info)
}



load_song_info_event = function(bird_name,dates, data_dir="/mnt/bengal_home/song", subdir="songs", batch_file=NULL, tz="PST", file_ex="wav.not.mat", log_data=NULL, evtaf_log=NULL) {
  wdirs = sapply(dates, function(x) {
    paste(data_dir, bird_name, x, subdir, sep="/")
  })
  infos = lapply(wdirs, function(w) {
    infos= load_mat_info(wdirs, file_ex = file_ex, batch_file=batch_file)
    #   infos= list(load_mat_info(wdirs[[1]], file_ex = file_ex, batch_file=batch_file))
  })
  
  info = bind_rows(infos)
  info$bird = bird_name
  
  ## Force timezone
  info = info %>% mutate(mtime = as.POSIXct(mtime, tz=tz))
  
  ## Parse time info --------------------------------------------------------------------------------
  if (file_ex == "wav.not.mat" | file_ex == "wav$") {
    info$parsed_time_epoch = parse_fname_for_timestamp(info$wav)
    info$parsed_time = as.POSIXct(parse_timestamp(info$wav), tz = tz)
  } else if (file_ex == "cbin.not.mat") {
    info$parsed_time = as.POSIXct(parse_timestamp_cbin(info$wav, prefix = sprintf("%s_", bird_name)))
  }

  #info = info %>% mutate(date = as.Date(parsed_time))
  info = info %>% mutate(date = floor_date(parsed_time, unit = "day"))
  info = info %>% mutate(time_h = as.numeric(time_h))
  info$date_noon = paste(info$date, "12:00:00", sep=" ")
  info$date_noon = parse_date_time(info$date_noon, orders="ymd hms")
  
  info$mtime_hour = floor_date(info$mtime, unit="hour")
  
  #info$manipulation = manipulation
  
  #info$date_noon = paste(info$date, "12:00:00", sep=" ")
  #info$date_noon = parse_date_time(info$date_noon, orders="ymd HMS", tz = "PST")
  
  #if (manipulation=="deaf") {
  #  info$deafened = info$date>event_date
  #} else {
  #  info$deafened = FALSE
  #}
  
  if (is.logical(log_data) && !log_data) {
    
    return(info)
  } else {
    if (is.null(log_data)) {
      log_data = load_knockdown_gs(bird_name)
    } 
    info = info %>% left_join(log_data)
    if ("log.date" %in% colnames(log_data))
      info = info %>% rename(event.date = log.date)
    
    info = info %>% mutate(rel_mtime = as.numeric(difftime(parsed_time, as.POSIXct(paste(event.date, "0:00:00 AM", sep=" ")), units="days")),
                           rel_date = as.numeric(difftime(as.Date(date), as.Date(event.date, tz=tz )), units="days"),
                      
                           rel_date_noon = as.numeric(difftime(as.Date(date_noon), as.Date(paste(event.date, "12:00:00 PM", sep=" ")), units="days")))
    
    info = info %>%
      select(-size, -isdir, -mode, -ctime, -atime, -uid, -gid, -uname, -grname)
    return(info)
  }
  
}

load_song_info_day = function(bird_name,  dates, data_dir="/mnt/bengal_home/song", subdir="songs", batch_file=NULL,  tz="PST", log_data=NULL, evtaf_log=NULL, file_ex="wav.not.mat", adjust_time=0) {
  wdirs = sapply(dates, function(x) {
    paste(data_dir, bird_name, x, subdir, sep="/")
  })
  x  = 1
  infos = lapply(wdirs, function(w) {
    #  x = 1
    print(w)
    infos= load_mat_info(wdirs, file_ex = file_ex, batch_file=batch_file)
    #   infos= list(load_mat_info(wdirs[[1]], file_ex = file_ex, batch_file=batch_file))
  })
  
  info = bind_rows(infos)
  info$bird = bird_name
  
  ## Force timezone
  info = info %>% mutate(mtime = as.POSIXct(mtime, tz=tz),
                         date = as.Date(date))
  ## Parse time info --------------------------------------------------------------------------------
  if (file_ex == "wav.not.mat") {
    info$parsed_time_epoch = parse_fname_for_timestamp(info$wav)
    info$parsed_time = as.POSIXct(parse_timestamp(info$wav), tz = tz)
  } else if (file_ex == "cbin\\.not\\.mat") {
    info$parsed_time = as.POSIXct(parse_timestamp_cbin(info$wav, prefix = sprintf("%s_", bird_name)))
  }
  
  info = info %>% mutate(time_h = as.numeric(time_h))
  info$date_noon = paste(info$date, "12:00:00", sep=" ")
  info$date_noon = parse_date_time(info$date_noon, orders="ymd hms", tz = tz)
  info$mtime_hour = floor_date(info$mtime, unit="hour")
  
  ## Combine with drug info ------------------------------------------------------------------------
  if (is.null(log_data)) {
    log_data = load_knockdown_gs(bird_name)
  }
  # dial_log = gs_title("RA cannulation: log") %>% gs_read("dialysis log") %>% filter(bird==bird_name)
  # colnames(dial_log) = paste("log", colnames(dial_log), sep=".")
  # dial_log$drug_time = with(dial_log,paste(log.date, log.time))
  # dial_log$drug_time = as.POSIXct(dial_log$drug_time, format="%Y-%m-%d %I:%M:%S %p")
  
  # if (nrow(dial_log)==0) {
  #   info$drug_time = NA
  # } else {
  #   
  #   
  #   #info$drug_start = as.POSIXct(cut(info$mtime, breaks=c(dial_log$drug_start, now()), right = T, include.lowest=F))
  #   #info$drug_time_stop = as.POSIXct(cut(info$mtime, breaks=c(dial_log$drug_time_stop, now()), right = T, include.lowest=F))
  #   info = left_join(info, dial_log, by="date")
  #   #info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
  #   
  #   # info1 = info %>% mutate(dummy=TRUE) %>%
  #   #   left_join(dial_log %>% mutate(dummy=TRUE)) %>%
  #   #   filter(parsed_time >= datetime_start, parsed_time <= datetime_stop) %>% select(-dummy)
  #   # 
  #   # info2 = info %>% anti_join(info1, by="mat_base")
  #   # info = bind_rows(info1, info2)
  #   # info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
  #   
  # }
  
  info = info %>% left_join(log_data)
  # To make compatible across machines
  info = info %>% mutate(mat_base = basename(mat))
  
  if (!is.null(evtaf_log)) {
    info1 = info %>% mutate(dummy=TRUE) %>%
      left_join(evtaf_log %>% mutate(dummy=TRUE) %>% select(-labels)) %>%
      filter(parsed_time >= datetime_start, parsed_time <= datetime_stop) %>% select(-dummy)
    
    info2 = info %>% anti_join(info1, by="mat_base")
    info = bind_rows(info1, info2)
    
    #info$evtaf_time = as.POSIXct(cut(info$mtime, breaks=c(evtaf_log$_time, now()), right = T, include.lowest=F))
    #info = left_join(info, dial_log, by="drug_time")
    #info = info %>% mutate(log.drug = ifelse(is.na(log.drug), "none", log.drug))
  }
  
  
  return(info)
}

load_dialysis_gs = function(bird_name) {
  dial_log = gs_title("RA cannulation: log") %>% gs_read("dialysis log") %>% filter(bird==bird_name)
  colnames(dial_log) = paste("log", colnames(dial_log), sep=".")
  dial_log$drug_start = with(dial_log,paste(log.date, log.time))
  dial_log$drug_start = as.POSIXct(dial_log$drug_start, format="%Y-%m-%d %H:%M:%S")
  
  dial_log$drug_stop = with(dial_log,paste(log.date_stop, log.time_stop))
  dial_log$drug_stop = as.POSIXct(dial_log$drug_stop, format="%Y-%m-%d %H:%M:%S")
  
  ## Adjust time with delivery time (tubing length * rate)
  dial_log = dial_log %>% mutate(drug_start = drug_start + seconds(`log.delivery time`)) %>% do({
    
    d = . 
    for (i in 1:(nrow(d)-1)) {
      d$drug_stop[i] = d$drug_stop[i] + seconds(d$`log.delivery time`[i+1])
    }
    d
  })

  return(dial_log)
}

load_evtaf_gs = function(bird_name) {
  lg = gs_title("evtaf") %>% gs_read("log") %>% filter(bird==bird_name)
  colnames(lg) = make.names(tolower(colnames(lg)))
  #colnames(lg) = paste("evtaf", colnames(lg), sep=".")
  lg = lg %>% mutate(datetime_start = as.POSIXct(paste(date.start, time.start, sep=" "), format="%Y-%m-%d %H:%M:%S"),
                     datetime_stop = ifelse(!is.na(date.stop),
                                            paste(date.stop, time.stop, sep =" "),
                                            as.character(now()))
                                            #Inf
                                            ) %>% 
    mutate(datetime_stop = as.POSIXct(datetime_stop, format="%Y-%m-%d %H:%M:%S"))
  
  #lg$datetime_stop = as.POSIXct("2100-01-01 01:00:00", format="%Y-%m-%d %H:%M:%S")
  #for (i in 1:(nrow(lg)-1)) {
  #  lg$datetime_stop[i]  = as.POSIXct(as.numeric(lg$datetime[i+1]), origin="1970-01-01")

  #}
  
  lg = lg %>% mutate(ff.threshold = ifelse(ff.threshold > 0, ff.threshold / 1000, ff.threshold))
  return(lg)
}

load_knockdown_gs = function(bird_name) {
  dial_log = gs_title("transfection/gene manipulation") %>% gs_read("siRNA injection log")
  colnames(dial_log) = tolower(make.names(colnames(dial_log)))
  if (length(bird_name) > 1) {
    dial_log = dial_log %>% filter(bird %in% bird_name)
  } else {
    dial_log = dial_log %>% filter(bird==bird_name)
  }
  colnames(dial_log) = paste("log", tolower(make.names(colnames(dial_log))), sep=".")
  dial_log = dial_log %>% rename(bird = log.bird)
  return(dial_log)
}

load_event_gs = function(bird_name, title="RA cannulation: log", sheet="log") {
  dial_log = gs_title(title) %>% gs_read(sheet)
  colnames(dial_log) = tolower(make.names(colnames(dial_log)))
  if (length(bird_name) > 1) {
    dial_log = dial_log %>% filter(bird %in% bird_name)
  } else {
    dial_log = dial_log %>% filter(bird==bird_name)
  }
  colnames(dial_log) = paste("event", tolower(make.names(colnames(dial_log))), sep=".")
  dial_log = dial_log %>% rename(bird = event.bird)
  return(dial_log)
}
  