teensy_load_log = function(bird, subdir="teensyTAF/TAFlogs", log_file=NULL, log_ind=1, fft_len=128) {
  data_dir = sprintf("/mnt/bengal_home/song/%s/%s", bird, subdir)
  data_files = sort(list.files(data_dir, full.names = T ), decreasing = T)
  
  if (is.null(log_file)) {
    log_file = file.path(data_files[log_ind])
  }
  print(log_file)
  
  cnames = c("time", "event_num",  "teensy_usec", "teensy_event_num", "teensy_event", paste("freqbin", 1:fft_len, sep=""))
  ctypes = c(list(col_character(), col_character(), col_double(), col_integer(), col_character()), map(1:fft_len, ~col_double()))
  logs = map(log_file, ~read_csv(.x, col_names=cnames, col_types = ctypes , skip=0)) %>% bind_rows()
  
  params = logs %>% filter(is.na(teensy_usec)) %>%
    select(time, event_num) %>%
    rename(variable = time,
           value = event_num)
  logs = logs %>% filter(!is.na(teensy_usec)) %>%
    mutate(time = as.numeric(time),
           event_num = as.integer(event_num))
  
  ## Reregister times
  logs = logs %>% group_by(teensy_event_num) %>%
    filter(any(teensy_event %in% "amp")) %>%
    mutate(time = time[teensy_event=="amp"] + (teensy_usec - teensy_usec[teensy_event=="amp"]) / 1E6)
  
  logs = logs %>% arrange(time) %>%
    distinct(teensy_event_num, teensy_event, .keep_all=T)
  
  list(params,logs)
}

teensy_extract_point_long = function(logs) {
  logs_point =  logs_filt %>% filter(teensy_event != "PSD") %>% 
    filter(!is.na(teensy_event_num)) %>%
    filter(teensy_event_num>0) %>%
    group_by(teensy_event_num) %>%
    select(time, event_num, teensy_usec, teensy_event_num, teensy_event, freqbin1) %>%
    rename(event_value = freqbin1) 
}
teensy_extract_point_stats = function(logs) {
  ## Extract point statistics
  logs_point = teensy_extract_point_long(logs)
  # logs_point =  logs %>% filter(teensy_event != "PSD") %>% 
  #   filter(!is.na(teensy_event_num)) %>%
  #   filter(teensy_event_num>0) %>%
  #   group_by(teensy_event_num) %>%
  #   select(time, event_num, teensy_usec, teensy_event_num, teensy_event, freqbin1) %>%
  #   rename(event_value = freqbin1) 
  
  logs_dp = logs_point %>%
    filter(teensy_event %in% c("amp", "DP", "action", "FF", "FREQTHRESH")) %>%
    select(teensy_event, teensy_event_num, event_value) %>%
    spread(teensy_event, event_value) %>%
    mutate(action = if_else(is.na(action), -1, action),
           FF = if_else(is.na(FF), -1, FF),
           FREQTHRESH = if_else(is.na(FREQTHRESH), -1, FREQTHRESH))
  
  logs_dp_time = logs_point %>%
    filter(teensy_event %in% c("DP")) %>%
    select(teensy_event_num, time)
  
  logs_dp = logs_dp %>% 
    left_join(logs_dp_time) %>%
    ungroup() %>%
    mutate(action = factor(action))
  logs_dp
}


# 
teensy_extract_dp = function(logs_point) {
  logs_dp = logs_point %>%
    filter(teensy_event %in% c("DP", "action")) %>%
    select(teensy_event, teensy_event_num, event_value) %>%
    spread(teensy_event, event_value) %>%
    mutate(action = if_else(is.na(action), -1, action))
  
  logs_dp_time = logs_point %>%
    filter(teensy_event %in% c("DP")) %>%
    select(teensy_event_num, time)
  
  logs_dp = logs_dp %>% 
    left_join(logs_dp_time) %>%
    ungroup() %>%
    mutate(action = factor(action))
  
  logs_dp
  
}