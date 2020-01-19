library(xml2)

flip_geometry_y = function(x, max_val) {
  (x - c(0, max_val)) * c(1,-1)
}

flip_geometry_xy = function(x, max_val) {
  (x - c(max_val, max_val)) * c(-1,-1)
}

rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
rotate_image = function(x, max_val, alpha) {
  centroid = c(max_val/2, max_val/2)
  ((x-centroid) * rot(alpha)) + centroid
}

is_touching = function(t1_radius, t2_radius) {
  t1_radius + t2_radius # touching
}

is_within = function(t1_radius, t2_radius) {
  ifelse(t1_radius >= t2_radius, t1_radius, t2_radius) # centroid smaller within larger blob
}

distance_funcs = list(touching = is_touching, within = is_within)

intersect_by_centroid_min = function(t1, t2, method, to_return="ids") {
  if (nrow(t1) == 0 || nrow(t2) == 0) {
    return(data.frame())
  }
  
  res_dist = st_distance(t1, t2)
  distance_tol_mat = matrix(0,  nrow=nrow(t1), ncol=nrow(t2))
  
  #distance_func = case_when(method == "touching" ~ is_touching,
  #                          method == "within" ~ is_within)
  distance_func = distance_funcs[[method]]
  for (i in 1:nrow(t1)) {
    for (j in 1:nrow(t2)) {
      distance_tol_mat[i,j] = distance_func(t1$radius[i], t2$radius[j])
      #distance_tol_mat[i,j] = ifelse(t1$radius[i] >= t2$radius[j], t1$radius[i], t2$radius[j]) # centroid smaller within larger blob
      #distance_tol_mat[i,j] = t1$radius[i] + t2$radius[j] # touching
    }
  }
  

  res_dist_filt = res_dist < distance_tol_mat# * res_dist_filt_min
  if (to_return=="matrix") {
    res_dist_ind = which(res_dist_filt, arr.ind = T)
    out = matrix(0, nrow = nrow(res_dist), ncol=ncol(res_dist))
    out[res_dist_ind] = res_dist[res_dist_ind]
    return(out)
  }
  
  res_dist_filt_min = apply(res_dist, 2, function(x) x == min(x))
  res_dist_filt = res_dist_filt & res_dist_filt_min

  res_dist_ind = which(res_dist_filt, arr.ind = T)
 
  
  #if (nrow(res_dist_ind == 0))
  res = data.frame(dapi_id = t1$id[res_dist_ind[,1]], probe_id = t2$id[res_dist_ind[,2]], probe=t2$probe[1])
  #print(res)
  return(res)
}

intersect_by_centroid = function(data, probe1, probe2, area_filter=NULL, method = c("touching", "within")) {
  if (!is.null(area_filter))
    data = data %>% filter(area_sf >= area_filter)
  
  data %>%
    filter(probe %in% c(probe1, probe2)) %>% 
    group_by(nucleus2, image) %>% do({
      d = .
      nuc = as.character(d$nucleus2[1])
      im = d$image[1]
      print(im)
      print(nuc)
      t1 = d %>% filter(probe==probe1)
      t2 = d %>% filter(probe==probe2)
      intersect_by_centroid_min(t1, t2, method=method)
      
    })
}


intersect_by_centroid_across_z = function(data, probe_to_use, method=c("touching", "within")) {
  data1 = data %>% filter(probe==probe_to_use)
  zs = unique(data$Z)
  comb_zs = map(1:(length(zs)-1), function(x) c(x, x+1))
  
  res = map(comb_zs, function(zs1) {
    t1 = data1 %>% filter(Z==zs1[1]) %>% st_zm(drop=T, what="ZM")
    t2 = data1 %>% filter(Z==zs1[2]) %>% st_zm(drop=T, what="ZM")
    intersect_by_centroid_min(t1, t2, method=method)
  })
}


calc_classes = function( data, ref_probe, group_by_class = T) {
  probes = as.character(unique(data$probe))
  compare_probes = probes[!(probes %in% c(ref_probe, "dapi"))]
  data_ref_probe = data %>% filter(probe==ref_probe) %>% ungroup() %>% select(dapi_id, image) %>%
    mutate(#dapi_id = as.character(dapi_id),
           image = as.character(image))
  data1 = data %>% inner_join(data_ref_probe)
  data2 = data1 %>%
    #distinct(dapi_id, image, nucleus2, probe_id, probe) %>% 
    #group_by(dapi_id, image, nucleus2) %>% 
    distinct(dapi_id, bird, brain_side, nucleus2, probe_id, probe) 
  
  if (group_by_class) {
   # data2 = data2 %>%
  #    filter(dapi_id==938) 
    data2 = data2 %>% 
      group_by(dapi_id, image, bird, brain_side, nucleus2)  %>% 
      summarize(class = case_when(all(unique(probe) %in% ref_probe) ~ "none",
                                  all(c(ref_probe, compare_probes) %in% unique(probe)) ~ "both",
                                  all(c(ref_probe, compare_probes[1]) %in% unique(probe)) ~ compare_probes[1],
                                  all(c(ref_probe, compare_probes[2]) %in% unique(probe)) ~ compare_probes[2],
        
                                  TRUE ~ "blah")) 
    
  } else {
    x = 1
    data2 = data2 %>%
      #filter(probe != ref_probe) %>%
      mutate(class = probe)
    #group_by(dapi_id, bird, brain_side, nucleus2)  %>% 
    #summarize(class = if_else(n()==1, "none", probe[probe != ref_probe]))
  }
  data2 = data2 %>% 
    mutate(probe_ref = ref_probe)
  data2
}

group_by_probe_ref = function(data, ref_probe) {
  probes = unique(data$probe)
  compare_probes = probes[!(probes %in% c(ref_probe, "dapi"))]
  data_ref_probe = data %>% filter(probe==ref_probe) %>% ungroup() %>% select(dapi_id, image)
  data1 = data %>% inner_join(data_ref_probe)
  data1 = data %>% mutate(probe_ref = ref_probe)
  data1
}

calc_intersection_statistics = function(data_inter, class_levels=NULL, group_by_class = T, remove_both=F) {

  if (group_by_class) {
    data_inter1 = data_inter %>% 
      group_by(brain_side, bird, nucleus2, class, probe_ref)
  } else{
    data_inter1 = data_inter %>% 
      group_by(brain_side, bird, nucleus2, probe, probe_ref)
  }

  if (remove_both) {
    data_inter1 = data_inter1 %>% filter(class!="both")
  }
  data_inter1 = data_inter1 %>% 
    summarize(n_objs=n()) 
  #data_inter1_probe_ref = data_inter1 %>% filter(class == probe_ref) %>%
  #  rename(n_probe_ref = n_objs) %>% 
  #  ungroup() %>% 
  #  select(brain_side, bird, nucleus2, probe_ref, n_probe_ref)
  
  #data_inter1 = data_inter1 %>% left_join(data_inter1_probe_ref)
  data_inter1 = data_inter1 %>%
    ungroup() %>%
    group_by(brain_side, bird, nucleus2, probe_ref) %>% 
    mutate(n_objs = if_else(class == probe_ref, n_objs - sum(n_objs[class != probe_ref]), n_objs))
  #%>% 
  data_inter1 = data_inter1 %>%
    ungroup() %>% 
    #group_by(image, nucleus2, probe_ref) %>%
    group_by(brain_side, bird, nucleus2, probe_ref)
  if (group_by_class) {
    data_inter1 = data_inter1 %>% 
      mutate(class = factor(class, levels=class_levels)) %>% 
      complete(class, fill = list(n_objs = 0))
  }
  #data_inter1 = data_inter1 
  data_inter1 = data_inter1 %>%
    #group_by(image, nucleus2, probe_ref) %>% 
    group_by(brain_side, bird, nucleus2, probe_ref) %>% 
    mutate(frac = n_objs / sum(n_objs),
           perc = frac * 100)
  
  #if (group_by_class) {
  #  data_inter1 = data_inter1 %>% filter(class != probe_ref)
  #}
  data_inter1
  
}

xml_to_df = function(fname) {
  xm = read_xml(fname)
  xm1 = xml_children(xm)
  xm2l = as.list(xml_children(xm1[[2]]))#[[2]]
  #xm2l = as.list(xm2)
  d = map_dfr(xm2l, function(xm2) {
    xm_type = as.numeric(xml_text(xml_find_all(xm2, "Type")))
    if (length(xm_type)==0) {
      return(data.frame())
    }
    xm_mark = xml_find_all(xm2, "Marker")
    
    if (length(xm_mark)==0) {
      return(data.frame())
    }
    
    xm_markx = as.numeric(xml_text(xml_find_all(xm_mark, "MarkerX")))
    xm_marky = as.numeric(xml_text(xml_find_all(xm_mark, "MarkerY")))
    
    xm_markz_ind = xml_find_all(xm_mark, "MarkerZ")
    
    if (length(xm_markz_ind)>0) {
      xm_markz = as.numeric(xml_text(xm_markz_ind))
      d1 = data.frame(X = xm_markx, Y = xm_marky, Z = xm_markz, measure_channel = xm_type)
    } else {
      d1 = data.frame(X = xm_markx, Y = xm_marky, measure_channel = xm_type)
    }
    d1
  })
  
  
  image_fname_base = str_split(basename(fname), "\\.")[[1]][[1]]
  image_fname_base = str_split(image_fname_base, "_counts")[[1]][[1]]
  d = d %>% mutate(image = image_fname_base)
  d
}

# <?xml version="1.0" encoding="UTF-8"?>
#   <CellCounter_Marker_File>
#   <Image_Properties>
#   <Image_Filename>Counter Window - pk74rd100_sl5_sec5_dcx_nectin3_neurod6_hvc_left_20x_tile_stitch.lsm</Image_Filename>
#   </Image_Properties>
#   <Marker_Data>
#   <Current_Type>2</Current_Type>
#   <Marker_Type>
#   <Type>1</Type>
#   </Marker_Type>
#   <Marker_Type>
#   <Type>2</Type>
#   <Marker>
#   <MarkerX>1225</MarkerX>
#   <MarkerY>650</MarkerY>
#   <MarkerZ>1</MarkerZ>
sf_to_xml = function(ds, image_file_ext=".lsm", out_fname) {
  #test1 = test %>% filter(probe=="cdh9")
  #coords = st_coordinates(test1)
  #ids = test1$id
  image_fname = paste0(ds$image[1], image_file_ext)

  header_list =  list(Image_Filename = list(image_fname))

  
  probes = levels(ds$probe)
  if(length(probes)==0) {
    stop("Probes must have levels")
  }
  
  marker_types = map(1:length(probes), function(i) {
    print(i)
    ds_p = ds %>% filter(probe==probes[i])
    coords = st_coordinates(ds_p)
    coords[,2] = (coords[,2] - ds$image_y_dim[1]) * -1
    print(nrow(coords))
    coords_list = lapply(1:nrow(coords), function(j) {
      list(MarkerX=list(coords[j,1]),
           MarkerY=list(coords[j,2]),
           MarkerZ=list(i))
    })
    names(coords_list) = rep("Marker", times=length(coords_list))
    type_list = c(list(Type=list(i)), coords_list)
    type_list
  })
  
  marker_types_empty = map(((length(probes)+1):8), function(i) {
    #names(coords_list) = rep("Marker", times=length(coords_list))
    type_list = list(Type=list(i))
    type_list
  })
  marker_types = c(marker_types, marker_types_empty)
  names(marker_types) = rep("Marker_Type", times=8)
  marker_data = c(list(Current_Type=list(2)), marker_types)
  #names(marker_data) = "Marker_Data"
  data_list = list(Image_Properties = header_list, Marker_Data = marker_data)
  out_list = list(CellCounter_Marker_File = data_list)
  #out_list = list(Marker_Data = list(Current_Type=2, Marker_Type=list(Type=2, Marker=list(MarkerX=1225, MarkerY=650, MarkerZ=1))))
  #out_list = list(Marker_data = list(Current_Type = list(`2` = list())))
  #out_list = list(Marker_data = list(Current_Type = list(2),
  #                                  Marker_Type=list(Type=list(1)),
  #                                  Marker_Type=list(Type = list(2), Marker = list(MarkerX = list(1225),
  #                                                                                 MarkerY = list(650),
  #                                                                                 MarkerZ = list(1)))))
  
  xml = as_xml_document(out_list)
  write_xml(xml, file = out_fname)
  # xm = read_xml(fname)
  # xm1 = xml_children(xm)
  # xm2l = as.list(xml_children(xm1[[2]]))#[[2]]
  # #xm2l = as.list(xm2)
  # d = map_dfr(xm2l, function(xm2) {
  #   xm_type = as.numeric(xml_text(xml_find_all(xm2, "Type")))
  #   if (length(xm_type)==0) {
  #     return(data.frame())
  #   }
  #   xm_mark = xml_find_all(xm2, "Marker")
  #   
  #   if (length(xm_mark)==0) {
  #     return(data.frame())
  #   }
  #   
  #   xm_markx = as.numeric(xml_text(xml_find_all(xm_mark, "MarkerX")))
  #   xm_marky = as.numeric(xml_text(xml_find_all(xm_mark, "MarkerY")))
  #   
  #   xm_markz_ind = xml_find_all(xm_mark, "MarkerZ")
  #   
  #   if (length(xm_markz_ind)>0) {
  #     xm_markz = as.numeric(xml_text(xm_markz_ind))
  #     d1 = data.frame(X = xm_markx, Y = xm_marky, Z = xm_markz, measure_channel = xm_type)
  #   } else {
  #     d1 = data.frame(X = xm_markx, Y = xm_marky, measure_channel = xm_type)
  #   }
  #   d1
  # })
  # 
  # 
  # image_fname_base = str_split(basename(fname), "\\.")[[1]][[1]]
  # image_fname_base = str_split(image_fname_base, "_counts")[[1]][[1]]
  # d = d %>% mutate(image = image_fname_base)
  # d
}

sf_to_xml2 = function(ds_cur, image_file_ext=".lsm", out_fname) {
  require(XML)
  
  probes = levels(ds_cur$probe)
  if(length(probes)==0) {
    stop("Probes must have levels")
  }
  
  image_fname = paste0(ds_cur$image[1], image_file_ext)
  
  xml <- xmlTree()
  # names(xml)
  xml$addTag("CellCounter_Marker_File", close=FALSE)
  xml$addTag("Image_Properties", close=FALSE)
  xml$addTag("Image_Filename", image_fname)
  xml$closeTag(name = "Image_Properties")
  
  xml$addTag("Marker_Data", close=F)
  xml$addTag("Current_Type", 2)
  
  for(i in 1:length(probes)) {
    print(i)
    xml$addTag("Marker_Type", close=F)
    xml$addTag("Type", i)
    
    ds_p = ds_cur %>% filter(probe==probes[i])
    coords = st_coordinates(ds_p)
    coords[,2] = (coords[,2] - ds$image_y_dim[1]) * -1
    print(nrow(coords))
    for (j in 1:nrow(coords)) {
      print(j)
      xml$addTag("Marker", close=F)
      xml$addTag("MarkerX", coords[j,1])
      xml$addTag("MarkerY", coords[j,2])
      xml$addTag("MarkerZ", i)
      xml$closeTag(name="Marker") # Marker
    }
    xml$closeTag(name="Marker_Type") # Marker_Type
  }
 
  
  for(i in ((length(probes)+1):8)) function(i) {
    xml$addTag("Marker_Type", close=F)
    xml$addTag("Type", i)
    xml$closeTag(name="Marker_Type") # Marker_Type
  }
  xml$closeTag(name="Marker_Data") # Marker_Data
  xml$closeTag(name="CellCounter_Marker_File") # CellCounter_Marker_File
  saveXML(xml, out_fname)
  cat(as(xml, "character"), file=out_fname, sep="\n")
  # 
  #     names(coords_list) = rep("Marker", times=length(coords_list))
#     type_list = c(list(Type=list(i)), coords_list)
#     type_list
  }
  # for (i in 1:nrow(data)) {
  #   xml$addTag("page", close=FALSE)
  #   for (j in names(data)) {
  #     xml$addTag(j, data[i, j])
  #   }
  #   xml$closeTag()
  # }
  # xml$closeTag()
  # xml$closeTag()
  # cat(saveXML(xml))

  
  #header_list =  list(Image_Filename = list(image_fname))
  
  

  
  
  
  
  # marker_types = c(marker_types, marker_types_empty)
  # names(marker_types) = rep("Marker_Type", times=8)
  # marker_data = c(list(Current_Type=list(2)), marker_types)
  # #names(marker_data) = "Marker_Data"
  # data_list = list(Image_Properties = header_list, Marker_Data = marker_data)
  # out_list = list(CellCounter_Marker_File = data_list)
  # #out_list = list(Marker_Data = list(Current_Type=2, Marker_Type=list(Type=2, Marker=list(MarkerX=1225, MarkerY=650, MarkerZ=1))))
  # #out_list = list(Marker_data = list(Current_Type = list(`2` = list())))
  # #out_list = list(Marker_data = list(Current_Type = list(2),
  # #                                  Marker_Type=list(Type=list(1)),
  # #                                  Marker_Type=list(Type = list(2), Marker = list(MarkerX = list(1225),
  # #                                                                                 MarkerY = list(650),
  # #                                                                                 MarkerZ = list(1)))))
  # 
  # xml = as_xml_document(out_list)
  # write_xml(xml, file = out_fname)
  # xm = read_xml(fname)
  # xm1 = xml_children(xm)
  # xm2l = as.list(xml_children(xm1[[2]]))#[[2]]
  # #xm2l = as.list(xm2)
  # d = map_dfr(xm2l, function(xm2) {
  #   xm_type = as.numeric(xml_text(xml_find_all(xm2, "Type")))
  #   if (length(xm_type)==0) {
  #     return(data.frame())
  #   }
  #   xm_mark = xml_find_all(xm2, "Marker")
  #   
  #   if (length(xm_mark)==0) {
  #     return(data.frame())
  #   }
  #   
  #   xm_markx = as.numeric(xml_text(xml_find_all(xm_mark, "MarkerX")))
  #   xm_marky = as.numeric(xml_text(xml_find_all(xm_mark, "MarkerY")))
  #   
  #   xm_markz_ind = xml_find_all(xm_mark, "MarkerZ")
  #   
  #   if (length(xm_markz_ind)>0) {
  #     xm_markz = as.numeric(xml_text(xm_markz_ind))
  #     d1 = data.frame(X = xm_markx, Y = xm_marky, Z = xm_markz, measure_channel = xm_type)
  #   } else {
  #     d1 = data.frame(X = xm_markx, Y = xm_marky, measure_channel = xm_type)
  #   }
  #   d1
  # })
  # 
  # 
  # image_fname_base = str_split(basename(fname), "\\.")[[1]][[1]]
  # image_fname_base = str_split(image_fname_base, "_counts")[[1]][[1]]
  # d = d %>% mutate(image = image_fname_base)
  # d
}