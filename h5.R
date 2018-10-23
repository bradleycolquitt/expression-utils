
test_fn = function(data) {
  groups(data)
}

modify_track2 = function(data, h5, target_name, operation, ...) {
  new_group_name = paste("", target_name, sep="/")
  print(new_group_name)
  
  
  if (existsGroup(h5, new_group_name))
    h5unlink(h5, new_group_name)
  createGroup(h5, new_group_name)
  
  # pull out attributes -----------------------------------
  scaf = list.datasets(data$h5, full.names=F)[1:5]

  # loop through scaffolds ---------------------------------
  res = map(scaf, function(s) {
    print(s)
    ptr = data$h5[s]
    #array = readDataSet(openDataSet(data$h5, s))
    print("array loaded")
    # apply function -----------------------------------------------------
    if (operation=="shrink")
      fnc = shrink_track1
    if (operation=="smooth")
      fnc = smooth_track
    
    res1 = fnc(ptr, ...)
    print("function applied")
    new_dataset = paste(new_group_name, s, sep="/")
    h5[new_dataset] = res1
    h5n = h5[new_dataset]
    
    # copy attributes ------------------------------------
    #h5attr(h5n, "Resolution") = common_attrs[["Resolution"]]
  })
  return(1)
}

reduce_tracks2 = function(data, h5, target_name, fnc) {
  new_group_name = paste("", target_name, sep="/")
  print(new_group_name)
  if (existsGroup(h5, new_group_name))
    h5unlink(h5, new_group_name)
  createGroup(h5, new_group_name)
  
  # pull out attributes -----------------------------------
  scaf = list.datasets(data$h5[[1]], full.names=F)
  #src_names = data[,src_col]
  #example_scaf1 = h5[paste(src_names[1], scaf[1], sep="/")]
  
  #cnames = c("Resolution", "AStart", "Start")
  #common_attrs = map(cnames, ~ h5attr(example_scaf1, .x)) %>% set_names(cnames)
  #unames = c("Stop", "AStop")
  
  # loop through scaffolds ---------------------------------
  res = map(scaf, function(s) {
    print(s)
    print(d)
    
    arrays = map(data$h5, function(d) {
      #to_get = paste(d, s, sep="/")
      readDataSet(openDataSet(d, s))
    })
    print("array loaded")
    # apply function -----------------------------------------------------
    res1 = do.call(fnc, arrays)
    print("function applied")
    new_dataset = paste(new_group_name, s, sep="/")
    h5[new_dataset] = res1
    h5n = h5[new_dataset]
    
    # copy attributes ------------------------------------
    h5attr(h5n, "Resolution") = common_attrs[["Resolution"]]
  })
  return(1)
}


reduce_tracks = function(data, h5, src_names, target_name, fnc) {
  new_group_name = paste("", target_name, sep="/")
  print(new_group_name)
  if (existsGroup(h5, new_group_name))
    h5unlink(h5, new_group_name)
  createGroup(h5, new_group_name)
  
  # pull out attributes -----------------------------------
  scaf = list.datasets(h5[src_names[1]], full.names=F)
  example_scaf1 = h5[paste(src_names[1], scaf[1], sep="/")]
  
  cnames = c("Resolution", "AStart", "Start")
  common_attrs = map(cnames, ~ h5attr(example_scaf1, .x)) %>% set_names(cnames)
  unames = c("Stop", "AStop")
  
  # loop through scaffolds ---------------------------------
  res = map(scaf, function(s) {
    arrays = map(src_names, function(d) {
      to_get = paste(d, s, sep="/")
      readDataSet(openDataSet(h5, to_get))
    })
    
    # apply function -----------------------------------------------------
    res1 = do.call(fnc, arrays)
    
    new_dataset = paste(new_group_name, s, sep="/")
    h5[new_dataset] = res1
    h5n = h5[new_dataset]
    
    # copy attributes ------------------------------------
    h5attr(h5n, "Resolution") = common_attrs[["Resolution"]]
  })
  return(1)
}


average_tracks = function(d_list) {
  out = vector(mode = "numeric", length = length(d_list[[1]]))
  for (i in 1:length(d_list)) {
    out = out + d_list[[i]]
  }
  out / length(d_list)
}

add_tracks = function(d_list) {
  print("here")
  print("here")
  out = vector(mode = "numeric", length = length(d_list[[1]]))
  for (i in 1:length(d_list)) {
    out = out + d_list[[i]]
  }
  out
}

shrink_track1 = function(vals, factor, fnc=mean) {
  vals = readDataSet(ptr)
  groups = rep(1:floor(length(vals) / factor), each=factor)
  groups = c(groups, rep(max(groups)+1, each=length(vals) - length(groups)))
  tapply(vals, groups, fnc)
  
  # inds = seq(1, length(vals)-1, factor)
  # sapply(inds, function(i) {
  #   fnc(vals[i:(i+factor-1)])
  # })
}

shrink_track2 = function(vals, factor, fnc=mean) {
  #groups = rep(1:floor(length(vals) / factor), each=factor)
  #groups = c(groups, rep(max(groups)+1, each=length(vals) - length(groups)))
  #tapply(vals, groups, fnc)
  
  inds = seq(1, length(vals)-1, factor)
  sapply(inds, function(i) {
    fnc(vals[i:(i+factor-1)])
  })
}
