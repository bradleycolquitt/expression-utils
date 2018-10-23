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

simple_dtw_mat = function(cost_mat, stepPattern=asymmetricP0, plot=F) {
  res = dtw(cost_mat, 
            open.end = T, open.begin = T, 
            step.pattern = stepPattern)
  if (plot) plot_dtw_result(ref_data, query_data, variable, res)
  return(res)
}

plot_dtw_result = function(ref_data, query_data, variable, dtw_result) {
  plot(as.data.frame(ref_data[dtw_result$index1, variable])[,1], type="l", main="reference", ylab=variable)
  lines(query_data[dtw_result$index2, variable], col=2)
}


dtw_expr_sub = function(ref_data, query_data, stepPattern, return_val="ind") {
  require(dtw)
  res = dtw(ref_data, query_data,
            open.end = F, open.begin = F, 
            step.pattern = stepPattern)
  
  if (return_val == "ind") {
    #df1 = res
    df1 = data.frame(ref=ref_data[res$index1], query=query_data[res$index2], ind=1:length(res$index2), index1=res$index1, index2=res$index2)
  } else if (return_val == "dist") {
    df1 = res$distance
  }
  return(df1)
}

test_dtw_expr_ind = function() {
  td = readRDS("~/data2/rstudio/birds/deafen/song_analysis/dialysis/wh27pk57/pitch_exprs.rds")
  td1 = td %>% ungroup() %>% filter(labels=="f")  %>% group_by(id) %>% 
    mutate(ff_norm = ff - mean(ff)) %>%
    ungroup()
  ids = unique(td1$id)
  td1 = td1 %>% filter(id %in% ids[1:5])
  out = dtw_expr(td1, stepPattern = asymmetric, ref_ind = 1, return_val = "ind")
  return(out)
}

test_dtw_expr_dist = function() {
  td = readRDS("~/data2/rstudio/birds/deafen/song_analysis/dialysis/wh27pk57/pitch_exprs.rds")
  td1 = td %>% ungroup() %>% filter(labels=="f") %>% top_n(20, wt = "id") %>% group_by(id) %>% 
    mutate(ff_norm = ff - mean(ff)) %>%
    ungroup()
  ids = unique(td1$id)
  out = dtw_expr(td1, 1, return_val="dist")
  return(out)
  par(mfrow=c(1,2))
  plot(ref_data, type="l")
  lines(query_data, col=2)
  plot(out[,c("ind", "ref")], type="l")
  lines(out[,c("ind", "query")], col=2)
}

# against one reference
dtw_expr = function(data, stepPattern=symmetric1, ref_ind=1, ref_id=NULL, return_val="ind") {
  if (is.null(ref_id)) {
    groups = unique(data$group)
    ref_id = groups[ref_ind]
  }
  ref_data = data  %>% filter(group == ref_id) %>% select(value) %>% unlist()
  out = data %>% group_by(group) %>% do({
    query_data = .$value
    out1 = dtw_expr_sub(ref_data, query_data, stepPattern=stepPattern, return_val=return_val)
    if (return_val=="dist") {
      data.frame(dtw_dist=out1)
    } else if (return_val=="ind") {
      data.frame(value_dtw=.$value[out1$index2], ind=1:length(out1$index2), pos = .$pos[out1$index2])
    }
  })
  return(out)
}

test_dtw_expr2_dist = function() {
  td = readRDS("~/data2/rstudio/birds/deafen/song_analysis/dialysis/wh27pk57/pitch_exprs.rds")
  td1 = td %>% ungroup() %>% filter(labels=="f")  %>% group_by(id) %>% 
    mutate(ff_norm = ff - mean(ff)) %>%
    ungroup()
  ids = unique(td1$id)
  td1 = td1 %>% filter(id %in% ids[1:5])
  print(length(unique(td1$id)))
  #return()
  out = dtw_expr2(td1, return_val="dist", parallel=T)
  return(out)
}

# against all references
dtw_expr2_sub = function(data, i, stepPattern=symmetric1, return_val="ind") {
  ref_data = data  %>% filter(group == i) %>% select(value) %>% unlist()

  out1 = data %>% group_by(group) %>% do({
    query_data = .$value
    out2 = dtw_expr_sub(ref_data, query_data, stepPattern=stepPattern, return_val=return_val)
    if (return_val=="dist") {
      data.frame(dtw_dist=out2)
    } else {
      out2
    }
  })
  return(out1)
}

dtw_expr2 = function(data, stepPattern=symmetric1, return_val="ind", parallel=F) {
  groups = unique(data$group)
  if (parallel) {
    out = foreach(i=groups) %dopar% dtw_expr2_sub(data, i, stepPattern=stepPattern, return_val=return_val)
  } else {
    out = foreach(i=groups) %do% dtw_expr2_sub(data, i, stepPattern=stepPattern, return_val=return_val)
  }
  names(out) = groups
  out_df = bind_rows(out, .id="ref_id")
  out_df
}

#' taken from dtw package
dtw_custom = function (x, y = NULL, dist.method = "Euclidean", step.pattern = symmetric2, 
                     window.type = "none", keep.internals = FALSE, distance.only = FALSE, 
                     open.end = FALSE, open.begin = FALSE, ...) 
{
  lm <- NULL
  print("here")
  if (is.null(y)) {
    if (!is.matrix(x)) 
      stop("Single argument requires a global cost matrix")
    lm <- x
  }
  else if (is.character(dist.method)) {
    x <- as.matrix(x)
    y <- as.matrix(y)
    lm <- proxy::dist(x, y, method = dist.method)
  }
  else if (is.function(dist.method)) {
    stop("Unimplemented")
  }
  else {
    stop("dist.method should be a character method supported by proxy::dist()")
  }
  wfun <- .canonicalizeWindowFunction(window.type)
  dir <- step.pattern
  norm <- attr(dir, "norm")
  if (!is.null(list(...)$partial)) {
    warning("Argument `partial' is obsolete. Use `open.end' instead")
    open.end <- TRUE
  }
  n <- nrow(lm)
  m <- ncol(lm)
  if (open.begin) {
    if (is.na(norm) || norm != "N") {
      stop("Open-begin requires step patterns with 'N' normalization (e.g. asymmetric, or R-J types (c)). See papers in citation().")
    }
    lm <- rbind(0, lm)
    np <- n + 1
    precm <- matrix(NA, nrow = np, ncol = m)
    precm[1, ] <- 0
  }
  else {
    precm <- NULL
    np <- n
  }
  gcm <- globalCostMatrix(lm, step.matrix = dir, window.function = wfun, 
                          seed = precm, ...)
  gcm$N <- n
  gcm$M <- m
  gcm$call <- match.call()
  gcm$openEnd <- open.end
  gcm$openBegin <- open.begin
  gcm$windowFunction <- wfun
  lastcol <- gcm$costMatrix[np, ]
  if (is.na(norm)) {
  }
  else if (norm == "N+M") {
    lastcol <- lastcol/(n + (1:m))
  }
  else if (norm == "N") {
    lastcol <- lastcol/n
  }
  else if (norm == "M") {
    lastcol <- lastcol/(1:m)
  }
  gcm$jmin <- m
  if (open.end) {
    if (is.na(norm)) {
      stop("Open-end alignments require normalizable step patterns")
    }
    gcm$jmin <- which.min(lastcol)
  }
  gcm$distance <- gcm$costMatrix[np, gcm$jmin]
  if (is.na(gcm$distance)) {
    stop("No warping path exists that is allowed by costraints")
  }
  if (!is.na(norm)) {u
    gcm$normalizedDistance <- lastcol[gcm$jmin]
  }
  else {
    gcm$normalizedDistance <- NA
  }
  if (!distance.only) {
    mapping <- backtrack(gcm)
    gcm <- c(gcm, mapping)
  }
  if (open.begin) {
    gcm$index1 <- gcm$index1[-1] - 1
    gcm$index1s <- gcm$index1s[-1] - 1
    gcm$index2 <- gcm$index2[-1]
    gcm$index2s <- gcm$index2s[-1]
    lm <- lm[-1, ]
    gcm$costMatrix <- gcm$costMatrix[-1, ]
    gcm$directionMatrix <- gcm$directionMatrix[-1, ]
  }
  if (!keep.internals) {
    gcm$costMatrix <- NULL
    gcm$directionMatrix <- NULL
  }
  else {
    gcm$localCostMatrix <- lm
    if (!is.null(y)) {
      gcm$query <- x
      gcm$reference <- y
    }
  }
  class(gcm) <- "dtw"
  return(gcm)
}

.canonicalizeWindowFunction <- function(w) {
  if(is.function(w)) {
    return(w);
  }
  
  # else 
  wt<-pmatch(w,c("none","sakoechiba","itakura","slantedband"));
  if(is.na(wt)) {
    stop("Ambiguous or unsupported char argument for window.type");
  } 
  
  wfun<-switch(wt,
               noWindow,
               sakoeChibaWindow,
               itakuraWindow,
               slantedBandWindow);
  
  return(wfun);
}

globalCostMatrix <-
  function(lm,
           step.matrix=symmetric1,
           window.function=noWindow,
           native=TRUE,
           seed=NULL,
           ...) {
    
    
    ## sanity check - be extra cautions w/ binary
    if (!is.stepPattern(step.matrix))
      stop("step.matrix is no stepMatrix object");
    
    
    
    # i = 1 .. n in query sequence, on first index, ie rows
    # j = 1 .. m on reference sequence, on second index, ie columns
    #   Note:  reference is usually drawn vertically, up-wise
    
    n <- nrow(lm);
    m <- ncol(lm);
    
    
    # number of individual steps (counting all patterns)
    nsteps<-dim(step.matrix)[1];
    
    
    # clear the cost and step matrix
    # these will be the outputs of the binary
    # for  cm use  seed if given
    if(!is.null(seed)) {
      cm <- seed;
    } else {
      cm <- matrix(NA,nrow=n,ncol=m);
      cm[1,1] <- lm[1,1];                 # Questionable.
    }
    
    sm <- matrix(NA,nrow=n,ncol=m);
    
    
    if(is.loaded("computeCM") && native){
      ## precompute windowing
      wm <- matrix(FALSE,nrow=n,ncol=m);
      wm[window.function(row(wm),col(wm),
                         query.size=n, reference.size=m,
                         ...)]<-TRUE;
      
      if(FALSE) {
        ## this call could be optimized. Copies are killing perf.
        out<-.C(computeCM,
                NAOK=TRUE,
                PACKAGE="dtw",
                ## IN
                as.integer(dim(cm)),               # int *s
                as.logical(wm),                    # int *wm
                as.double(lm),                     # double *lm
                as.integer(nsteps),                # int *nstepsp
                as.double(step.matrix),            # double *dir
                ## IN+OUT
                costMatrix=as.double(cm),                 # double *cm
                ## OUT
                directionMatrix=as.integer(sm));               # int *sm
        
        ## Hopefully avoids a copy
        dim(out$costMatrix) <- c(n,m);     
        dim(out$directionMatrix) <- c(n,m);
        warning("You should not be here");
      } else {
        storage.mode(wm) <- "logical";
        storage.mode(lm) <- "double";
        storage.mode(cm) <- "double";
        storage.mode(step.matrix) <- "double";
        out <- .Call("computeCM_Call",
                     wm,lm,cm,step.matrix);
      }
      
    } else {
      
      ####################
      ## INTERPRETED PURE-R IMPLEMENTATION
      warning("Native dtw implementation not available: using (slow) interpreted fallback");
      # now walk through the matrix, column-wise and row-wise,
      # and recursively compute the accumulated distance. Unreachable
      # elements are handled via NAs (removed)
      dir <- step.matrix;
      npats <- attr(dir,"npat");
      for (j in 1:m) {
        for (i in 1:n) {
          ## It is ok to window on the arrival point (?)
          if(!window.function(i,j, query.size=n, reference.size=m, ...)) { next; }
          
          ## Skip if already initialized
          if(!is.na(cm[i,j])) { next; }
          
          clist<-numeric(npats)+NA;
          for (s in 1:nsteps) {
            ## current pattern
            p<-dir[s,1];
            ## ii,jj is the cell from where potentially we could
            ## have come from. 
            ii<-i-dir[s,2];                 # previous step in inp
            jj<-j-dir[s,3];                 # previous step in tpl
            if(ii>=1 && jj>=1) {            # element exists?
              cc<-dir[s,4];                 # step penalty
              if(cc == -1) {                #  -1? cumulative cost:
                clist[p]<-cm[ii,jj];	#  there must be exactly 1 per pattern
              } else {			#  a cost for 
                clist[p]<-clist[p]+cc*lm[ii,jj];
              }
            }
          }
          
          
          ## no NAs in clist at this point BUT clist can be empty
          ## store in cost matrix
          minc<-which.min(clist);           # pick the least cost
          if(length(minc)>0) {          	# false if clist has all NAs
            cm[i,j]<-clist[minc];
            sm[i,j]<-minc;			# remember the pattern picked
          }
        }
      }
      out <- list(costMatrix=cm,directionMatrix=sm);
    }
    
    ## END PURE-R IMPLEMENTATION
    ####################
    
    ## At this point out$cmo and out$smo should be set
    out$stepPattern <- step.matrix;
    return(out);
  }

is.stepPattern <- function(x) {
  return(inherits(x,"stepPattern"));
}


