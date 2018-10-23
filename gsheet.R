library(XML)
library(googlesheets)

# sleep_time is kludge to get around 429 too many requests error
load_lcm_sample_info = function(ssname, sleep_time=1) {
  gdb = gs_title(ssname)
  lcm = gdb %>% gs_read("lcm")
  Sys.sleep(sleep_time)
  positions = gdb %>% gs_read("positions")
  Sys.sleep(sleep_time)
  orient = gdb %>% gs_read("orientation")
  Sys.sleep(sleep_time)
  i7 = gdb %>% gs_read("i7")
  Sys.sleep(sleep_time)
  ind = grep("division.id", colnames(i7))
  val = unique(i7[,grep("division$", colnames(i7))])
  colnames(i7)[ind] = paste("lib", val, sep=".")
  bc = gdb %>% gs_read("bc")
  Sys.sleep(sleep_time)
  prep = gdb %>% gs_read("prep")
  Sys.sleep(sleep_time)

  info = inner_join(lcm, positions)
  info = inner_join(info, orient)
  info = inner_join(info, prep)
  info = inner_join(info, bc)
  info = inner_join(info, i7)
  info$i7 = factor(info$i7)
  info$bc = factor(info$bc)
  info$lib.column = factor(info$lib.column)
  info$id = 1:nrow(info)
  rm(gdb)
  return(info)
}


cleanGoogleTable <- function(dat, table=1, skip=0, ncols=NA, nrows=-1, header=TRUE, dropFirstCol=NA){
  if(!is.data.frame(dat)){
    dat <- dat[[table]]
  }
  if(is.na(dropFirstCol)) {
    firstCol <- na.omit(dat[[1]])
    if(all(firstCol == ".") || all(firstCol== as.character(seq_along(firstCol)))) {
      dat <- dat[, -1]
    }
  } else if(dropFirstCol) {
    dat <- dat[, -1]
  }
  if(skip > 0){
    dat <- dat[-seq_len(skip), ]
  }
  if(nrow(dat) == 1) return(dat)
  if(nrow(dat) >= 2){
    if(all(is.na(dat[2, ]))) dat <- dat[-2, ]
  }
  if(header && nrow(dat) > 1){
    header <- as.character(dat[1, ])
    names(dat) <- header
    dat <- dat[-1, ]
  }
  # Keep only desired columns
  if(!is.na(ncols)){
    ncols <- min(ncols, ncol(dat))
    dat <- dat[, seq_len(ncols)]
  }
  # Keep only desired rows
  if(nrows > 0){
    nrows <- min(nrows, nrow(dat))
    dat <- dat[seq_len(nrows), ]
  }
  # Rename rows
  rownames(dat) <- seq_len(nrow(dat))
  dat
}

readGoogleSheet <- function(url, na.string="", header=TRUE){
  stopifnot(require(XML))
  # Suppress warnings because Google docs seems to have incomplete final line
  suppressWarnings({
    doc <- paste(readLines(url), collapse=" ")
  })
  if(nchar(doc) == 0) stop("No content found")
  htmlTable <- gsub("^.*?(<table.*</table).*$", "\\1>", doc)
  ret <- readHTMLTable(htmlTable, header=header, stringsAsFactors=FALSE, as.data.frame=TRUE)
  lapply(ret, function(x){ x[ x == na.string] <- NA; x})
}
