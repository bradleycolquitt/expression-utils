library(dtw)
library(tidyr)

align_by_areas = function(info, groups, plot=T, ref_id=NULL) {
  dots = lapply(groups, as.symbol)
  areas = info.ra %>% group_by_(.dots=dots) %>% arrange(um_from_dorsal_local)
  areas$group_id = group_indices(areas)
  #info$group_id = group_indices(areas)
  group_id = unique(areas$group_id)
  
  sp = asymmetricP0
  sp = rabinerJuangStepPattern(7, "c", FALSE )
  
  if (is.null(ref_id)) {
    areas.align = lapply(group_id, function(ref_id) {
      tmp = lapply(group_id[group_id!=ref_id], function(x) {
        d1 = areas[areas$group_id==ref_id,]; d2 = areas[areas$group_id==x,]
        res = simple_dtw(d1, d2, "area_norm", stepPattern = sp, plot=F)
        #res = with(areas, dtw(d1$area_norm, d2$area_norm, open.end = T, open.begin = T, step.pattern = sp))
        return(res)
      })
      names(tmp) = group_id[group_id!=ref_id]
      return(tmp)
    })
    names(areas.align) = group_id
  
  
  ### Plot alignments
  if (plot) {
    par(mfrow=c(2,length(areas.align)/2))
    for (ref_id in names(areas.align)) {
      d1 = areas[areas$group_id==ref_id,]
      plot(d1$area_norm[areas.align[[ref_id]][[1]]$index1], type="l", main=ref_id, ylab="area_norm")
      for (x in group_id[group_id!=ref_id]) {
        d2 = areas[areas$group_id==x,]
        lines(d2$area_norm[areas.align[[ref_id]][[as.character(x)]]$index2], col=2)
      }
    }
  }

  names(areas.align) = group_id
  return(areas.align)
  }
  
  ### ref_id specified
  #areas.align = lapply(group_id, function(ref_id) {
  #return(areas)
  areas.align = lapply(group_id[group_id!=ref_id], function(x) {
      d1 = areas[areas$group_id==ref_id,]; d2 = areas[areas$group_id==x,]
      res1 = simple_dtw(d1, d2, "area_norm", stepPattern = sp, plot=F)
      res2 = simple_dtw(d2, d1, "area_norm", stepPattern = sp, plot=F)
      df1 = na.omit(data.frame(ref=d1$id[res1$index1], query=d2$id[res1$index2], group_id=x, set=1))
      df1$pos1 = (1:nrow(df1)) / nrow(df1)
      df2 = na.omit(data.frame(ref=d2$id[res2$index1], query=d1$id[res2$index2], group_id=x, set=2, um_from_dorsal_local=d2$um_from_dorsal_local[res2$index1]))
      df2$ref_pos = df1$pos1[match(df2$ref, df1$query)]
      splinef = splinefun(x = df2$um_from_dorsal_local, y = df2$ref_pos)
      df2$ref_pos_inter = splinef(df2$um_from_dorsal_local)
      df2$ref_pos_inter[df2$ref_pos_inter<0] = 0

      m = df2 %>% gather(compare, id, ref:query) %>% group_by(group_id, compare, id) %>% summarize(dtw_pos = mean(ref_pos_inter))
      return(m)
  })
    df = do.call(rbind, areas.align)
    df = df %>% group_by(id) %>% summarize(dtw_pos = mean(dtw_pos))
    df$ref_id = ref_id
    #return(df)
    return(list(info=areas, alignments=df))
  #})
  #df = do.call(rbind, areas.align)
  #return(list(info=info, alignments=df))
}

simple_dtw = function(ref_data, query_data, variable, stepPattern = asymmetricP0, plot = F) {
  res = dtw(ref_data$area_norm, query_data$area_norm, 
            open.end = T, open.begin = T, 
            step.pattern = stepPattern)
  if (plot) plot_dtw_result(ref_data, query_data, variable, res)
  return(res)
  
}

plot_dtw_result = function(ref_data, query_data, variable, dtw_result) {
  plot(as.data.frame(ref_data[dtw_result$index1, variable])[,1], type="l", main="reference", ylab=variable)
  lines(query_data[dtw_result$index2, variable], col=2)
}


