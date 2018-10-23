library(lubridate)

calc_simple_song_stats = function(tags, refdate=NULL) {
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
    return(data.frame(last_song_datetime=last_song_datetime,
                      rate_per_hour,
                      num_songs_on_euth_date,
                      nsongs_since_refdate))
  } else {
    return(data.frame(last_song_datetime=last_song_datetime,
                      rate_per_hour,
                      num_songs_on_euth_date))

  }
}