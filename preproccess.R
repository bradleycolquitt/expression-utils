library(lubridate)

calc_simple_song_stats = function(tags, refdate=NULL, reftime_end=NULL) {
  dir = paste("/mnt/bengal_home/song", tags, "songs", sep="/")
  files = list.files(dir, pattern = "*wav$", full.names = T)
  info = file.info(files)
  
  last_song_datetime = max(info$mtime)
  info$date = floor_date(info$mtime, unit = "day")
  num_songs_on_euth_date = info %>% dplyr::filter(date==max(date)) %>% summarize(num_songs_on_euth_date=n(),
                                                                                 first_song_datetime = min(mtime))
  

  
  info$hour = floor_date(info$mtime, unit = "hour")
  rate_per_hour = info %>% group_by(hour) %>% summarize(num_songs=n()) %>% dplyr::filter(num_songs>0) %>% summarize(song_rate_med = median(num_songs))

  if (!is.null(refdate)) {
    refdate = as.POSIXct(as.character(refdate))
    nsongs_since_refdate = nrow(info[info$date>refdate,])
    if (is.null(reftime_end)) {
      return(data.frame(last_song_datetime=last_song_datetime,
                        rate_per_hour,
                        num_songs_on_euth_date,
                        nsongs_since_refdate))
      
    } else if (!is.null(reftime_end)) {
      nsongs_in_last_hour = info %>% dplyr::filter(mtime >= (reftime_end - hours(1))) %>% summarize(nsongs_in_last_hour = n())
      nsongs_in_last_two_days = info %>% dplyr::filter(mtime >= (reftime_end - days(2))) %>% summarize(nsongs_in_last_two_days = n())
      nsongs_in_last_week = info %>% dplyr::filter(mtime >= (reftime_end - days(7))) %>% summarize(nsongs_in_last_two_days = n())
      nsongs_per_day_pre = info %>% dplyr::filter(mtime <= refdate) %>% group_by(date) %>% summarize(n = n()) %>% ungroup() %>% summarize(nsongs_per_day_pre = mean(n))
      nsongs_per_day_post = info %>% dplyr::filter(mtime > refdate) %>% group_by(date) %>% summarize(n = n()) %>% ungroup() %>% summarize(nsongs_per_day_post = mean(n))
      
      minutes_from_last_song = difftime(reftime_end, last_song_datetime, units = "mins")
      minutes_from_last_song = ifelse(minutes_from_last_song<0, 0, minutes_from_last_song)
      # minutes_from_last_song = info %>%  dplyr::mutate(minutes_from_last_song = difftime(reftime_end, last_song_datetime, units = "mins")) %>% 
      #   dplyr::mutate(minutes_from_last_song = ifelse(minutes_from_last_song<0, 0, minutes_from_last_song))
      return(data.frame(
        last_song_datetime=last_song_datetime,
        minutes_from_last_song,
        rate_per_hour,
        num_songs_on_euth_date,
        nsongs_in_last_hour,
        nsongs_in_last_two_days,
        nsongs_in_last_week,
        nsongs_per_day_pre,
        nsongs_per_day_post,
        nsongs_since_refdate))
    }
  } else {
    return(data.frame(last_song_datetime=last_song_datetime,
                      rate_per_hour,
                      num_songs_on_euth_date))

  }
}