deseq2_sub = function(dat, info, predictors, contrasts, norm_genes=NULL, plot_prefix=NULL) {
  predictors1 = as.formula(paste("~", paste(predictors, collapse="+"), sep=""))
  dds = NULL
    tryNULL({
      dds = DESeqDataSetFromMatrix(countData = dat,
                                    colData = info,
                                    design = predictors1)
      
      if (is.null(norm_genes)) {
        dds = estimateSizeFactors(dds)        
      } else {
        print("here")
        dds = estimateSizeFactors(dds, controlGenes=norm_genes)
      }
      dds = estimateDispersions(dds)
      dds = DESeq(dds, parallel = TRUE)
    })
  
  
  if (is.null(dds))
    return(NULL)
    
    
  ### CONTRASTS -------------------------------------------
  res = lapply(contrasts, function(x) {
    tryNULL({
      res = results(dds, name = x)
      
      if (!is.null(plot_prefix)) {
        plot_fname = sprintf("%s_%s.pdf", plot_prefix, x)
        cairo_pdf(plot_fname, width=5, height=5)
        my_ma_plot(res, "baseMean", "log2FoldChange", "padj", thresh=.1)
        dev.off()
      }
      
      res = as.data.frame(res)
      res$gene_id = rownames(res)
      res$coef = x
      res
    })
  })
  res = res[unlist(lapply(res, function(x) !is.null(x)))]
  res = bind_rows(res)
  
  #res[,factors[2]] = l2
  return(res)
}

deseq2_one_level = function(dat, info, predictors, factors, contrasts, outdir, model_fname, blocking=NULL, 
                          normalize.method="none", idcol="id") {
  
  if (!all(predictors %in% colnames(info)))
    stop("ERROR: Not all predictors in info data.frame.")
  
  prefix1 = paste(outdir, model_fname, sep="/")
  dir.create(prefix1, recursive=T)
  
  info_df = as.data.frame(info)
  levels1 = levels(info_df[,factors[1]])
  #levels2 = levels(info_df[,factors[2]])
  
  fname = paste(prefix1, "results.rds", sep="/")
  fname_pred = paste(prefix1, "predictors.txt", sep="/")
  
  predictors1 = as.formula(paste("~", paste(predictors, collapse="+"), sep=""))
  
  write.table(predictors, fname_pred)
  curr = lapply((1:length(levels1)), function(i1) {
    print(levels1[i1])
    l1 = levels1[i1]
    #curr = lapply(1:length(levels2), function(i2) {
    #  l2 = as.character(levels2[i2])
    #  print(l2)
    
    ### Setup metadata  --------------------------------------
    dots = lazyeval::interp(~f1==l1, .values=list(f1=as.name(factors[1])))
    info_curr = droplevels(info %>% filter_(.dots=dots ))
    info_curr_df = as.data.frame(info_curr)
    dat1 = dat[,match(info_curr_df[,idcol], colnames(dat))]
    
    res = voom_sub(dat1, info_curr, predictors, contrasts, normalize.method, blocking)
    if (is.null(res))
      return(NULL)
    res[,factors[1]] = l1
    return(res)
  })
  curr = curr[unlist(lapply(curr, function(x) !is.null(x)))]
  curr = bind_rows(curr)
  if (is.null(curr))
    return(NULL)
  curr = bind_rows(curr)
  dots = list(lazyeval::interp(~factor(f1, levels=levels1), .values=list(f1=as.name(factors[1]))))
  names(dots) = c(factors[1])
  curr = curr %>% mutate_(.dots=dots)
  saveRDS(curr, fname)
  curr
  
}

deg_two_levels = function(dat, info, predictors, 
                          factors, 
                          contrasts, 
                          method="voom",
                          outdir, 
                          model_fname, 
                          blocking=NULL, 
                          normalize.method="none",
                          norm_genes = NULL,
                          idcol="id") {
  
  if (!all(predictors %in% colnames(info)))
    stop("ERROR: Not all predictors in info data.frame.")
  
  prefix1 = paste(outdir, model_fname, sep="/")
  dir.create(prefix1, recursive=T)
  
  info_df = as.data.frame(info)
  levels1 = levels(info_df[,factors[1]])
  levels2 = levels(info_df[,factors[2]])
  
  fname = paste(prefix1, "results.rds", sep="/")
  fname_pred = paste(prefix1, "predictors.txt", sep="/")
  
  predictors1 = as.formula(paste("~", paste(predictors, collapse="+"), sep=""))
  
  write.table(predictors, fname_pred)
  
  ma_dir = paste(prefix1, "MA", sep="/")
  dir.create(ma_dir)
  
  res = lapply((1:length(levels1)), function(i1) {
    print(levels1[i1])
    l1 = levels1[i1]
    curr = lapply(1:length(levels2), function(i2) {
      l2 = as.character(levels2[i2])
      print(l2)
      
      ### Setup metadata  --------------------------------------
      dots = lazyeval::interp(~f1==l1 & f2==l2, .values=list(f1=as.name(factors[1]),
                                                             f2=as.name(factors[2])))
      info_curr = droplevels(info %>% filter_(.dots=dots ))
      info_curr_df = as.data.frame(info_curr)
      dat1 = dat[,match(info_curr_df[,idcol], colnames(dat))]
      
      res = NULL
      if (method == "voom") {
        res = voom_sub(dat1, info_curr, predictors, contrasts, normalize.method, blocking)        
      } else if (method == "deseq2") {
        plot_prefix = paste(ma_dir, sprintf("MA_%s_%s", l1, l2), sep="/")
        res = deseq2_sub(dat1, info_curr, predictors, contrasts, norm_genes, plot_prefix)
      }

      if (is.null(res))
        return(NULL)
      res[,factors[2]] = l2
      return(res)
    })
    curr = curr[unlist(lapply(curr, function(x) !is.null(x)))]
    curr = bind_rows(curr)
    if (is.null(curr))
      return(NULL)
    curr[,factors[1]] = l1
    curr
  })
  res = bind_rows(res)
  dots = list(lazyeval::interp(~factor(f1, levels=levels1), .values=list(f1=as.name(factors[1]))), 
              lazyeval::interp(~factor(f2, levels=levels2), .values=list(f2=as.name(factors[2]))))
  names(dots) = c(factors[1], factors[2])
  res = res %>% mutate_(.dots=dots)
  saveRDS(res, fname)
  res
  
}

my_ma_plot = function(data, avg_value, fc_value, p_value, thresh=.1, ylim=NULL) {
  data = as.data.frame(data)
  dots = list(lazyeval::interp(~v<thresh, v=as.name(p_value)),
              lazyeval::interp(~log10(v+1), v=as.name(avg_value)))
  names(dots) = c("sig", "avg_log")
  data = data %>% mutate_(.dots=dots)
  alpha_value = 1/10
  gg = ggplot(data, aes_string("avg_log", fc_value, color="sig"))
  gg = gg + scale_color_manual(values=c("black", "red"))
  gg = gg + geom_point(alpha=alpha_value)
  gg = gg + stat_smooth(color="blue")
  gg = gg + labs(x="mean expression")
  gg = gg + theme(legend.position="none")
  if (!is.null(ylim))
    gg = gg + lims(y=ylim)
  print(gg)
  
  
}