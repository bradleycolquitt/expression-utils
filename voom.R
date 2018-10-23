source("~/data2/rstudio/birds/utils/common_aesthetics.R")
library(NMF)
library(limma)
library(gridExtra)
library(stringr)
library(reshape2)
rename = dplyr::rename

multilevel_voom_fit = function(data, info, design, block_var) {
  nf = estimateSizeFactorsForMatrix(data)
  lib_size = colSums(data) * nf
  v1 = voom(data, design, lib.size=lib_size)
  corfit = duplicateCorrelation(v1, design,  block=info[,block_var])
  v1 = voom(data, 
            design, 
            lib.size=lib_size, 
            block=info[,block_var], 
            correlation=corfit$consensus.correlation)
  fit = lmFit(v1, 
              design, 
              block=info[,block_var], 
              correlation=corfit$consensus.correlation)
  fit
}

voom_sub = function(dat, info, predictors, contrasts, normalize.method="none", lib.size=NULL, blocking=NULL) {
  predictors1 = as.formula(paste("~", paste(predictors, collapse="+"), sep=""))
  fit = NULL
  dat = dat[,match(info$id, colnames(dat))]
  
  if (normalize.method=="TMM") {
    dge = DGEList(dat)
    dge= calcNormFactors(dge,method =c("TMM"))
    normalize.method="none"
  }
  if (!is.null(blocking)) {
    # info_tab = table(info[, predictors[1:2]])
    # info_tab[info_tab == 1] = 0
    # info_tab_sel = apply(info_tab, 2, prod)
    # info_tab_names = names(info_tab_sel)[info_tab_sel>0]
    # 
    # dots = lazyeval::interp(~v %in% info_tab_names, v=as.name(predictors[2]))
    # info = droplevels(info %>% filter_(dots))
    info_df = as.data.frame(info)
   
    
    #tryCatch({
      ### DESIGN -----------------------------------------------
      design = model.matrix(predictors1, info)
     
      ### FIT --------------------------------------------------
      if (is.null(lib.size)) {
        v1 = voom(dge,design,plot=F, normalize.method=normalize.method)
      } else {
        v1 = voom(dge,design,plot=F, normalize.method=normalize.method, lib.size=lib.size)       
      }
 
      corfit = duplicateCorrelation(v1, design, block = info_df[,blocking])
      fit = lmFit(v1, design, block = info_df[,blocking], correlation =
                    corfit$consensus.correlation)
    #}, error = function(e) {
    #  return(NULL)
    #})
  } else {
    tryCatch({
    ### DESIGN -----------------------------------------------
    design = model.matrix(predictors1, info)
    
    ### FIT --------------------------------------------------
    if (is.null(lib.size)) {
      v1 = voom(dat,design,plot=F, normalize.method=normalize.method)
    } else {
      v1 = voom(dat,design,plot=F, normalize.method=normalize.method, lib.size=lib.size)       
    }
    fit = lmFit(v1, design)    
    }, error = function(e) {
      return(NULL)
    })
  }
  
  if (is.null(fit))
    return(NULL)
  
  ### CONTRASTS --------------------------------------------
  
  res = lapply(contrasts, function(x) {
    tryCatch({
      fit2 = contrasts.fit(fit, coefficients = x)   
      fit2 = eBayes(fit2)
      res = topTable(fit2, p.value=1, number=nrow(dat), confint=T)
      res$gene_id = rownames(res)
      res$coef = x
      res
    }, error = function(e) {
      return(NULL)
    })
  })
  res = res[unlist(lapply(res, function(x) !is.null(x)))]
  res = bind_rows(res)
  
  #res[,factors[2]] = l2
  return(res)
}

voom_one_level = function(dat, info, predictors, factors, contrasts, outdir, model_fname, blocking=NULL, 
                           normalize.method="none", idcol="id",  scde_path=NULL) {
  
  # if (!all(predictors %in% colnames(info)))
  #   stop("ERROR: Not all predictors in info data.frame.")
  
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
      
      if (!is.null(scde_path)) {
        norms = readRDS(paste(scde_path, "model.rds", sep="/"))
        norms$id = as.numeric(rownames(norms))
        norms = norms %>% left_join(info_curr_df, by="id")
        norm_cur = norms %>% filter(id %in% info_curr_df$id)
        norm_cur = norm_cur[match(info_curr_df$id, norm_cur$id),]
        norm_cur = norm_cur[!is.na(norm_cur$id),]
        norm_vals = norm_cur$corr.b - min(norm_cur$corr.b)
        dat1 = dat1[,match(norm_cur$id, colnames(dat1))]
        info_curr_df = info_curr_df[match(colnames(dat1), info_curr_df$id),]
        res = voom_sub(dat1, info_curr_df, predictors, contrasts, normalize.method, lib.size=norm_vals, blocking=blocking)
      } else {
        res = voom_sub(dat1, info_curr_df, predictors, contrasts, normalize.method, blocking)        
      }
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

voom_two_levels = function(dat, info, predictors, factors, contrasts, outdir, model_fname, blocking=NULL, 
                           normalize.method="none", idcol="id", scde_path=NULL, scde_obj=NULL, scde_grouping_var="position") {
  
  if (!all(predictors %in% colnames(info)))
    stop("ERROR: Not all predictors in info data.frame.")
  
  prefix1 = paste(outdir, model_fname, sep="/")
  dir.create(prefix1, recursive=T)
  
  info_df = as.data.frame(info)
  levels1 = levels(info_df[,factors[1]])
  levels2 = levels(info_df[,factors[2]])
  
  if (is.null(levels1) || is.null(levels2))
    stop("Factors must have levels.")
  
  fname = paste(prefix1, "results.rds", sep="/")
  fname_pred = paste(prefix1, "predictors.txt", sep="/")

  predictors1 = as.formula(paste("~", paste(predictors, collapse="+"), sep=""))

  write.table(predictors, fname_pred)
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
      dat1 = dat[,colnames(dat) %in% info_curr_df[,idcol]]
      info_curr_df = info_curr_df[info_curr_df[,idcol] %in% colnames(dat1),]
      #info_curr_df = info_curr_df[match(colnames(dat), info_curr_df[,idcol]),]
      dat1 = dat1[,match(info_curr_df[,idcol], colnames(dat1))]
    
      
      if (!is.null(scde_path) || !is.null(scde_obj)) {
        if (is.null(scde_obj)) {
          scde_subdir = paste(scde_path, unique(info_curr_df[,scde_grouping_var]), sep="/")
          norms = readRDS(paste(scde_subdir, "model.rds", sep="/"))
        }
        if (is.null(scde_path))
          norms = scde_obj
        norms$id = as.numeric(rownames(norms))
        norms = norms %>% left_join(info_curr_df, by="id")
        norm_cur = norms %>% filter(id %in% info_curr_df$id)
        norm_cur = norm_cur[match(info_curr_df$id, norm_cur$id),]
        norm_cur = norm_cur[!is.na(norm_cur$id),]
        norm_vals = norm_cur$corr.b - min(norm_cur$corr.b)
        dat1 = dat1[,match(norm_cur$id, colnames(dat1))]
        info_curr_df = info_curr_df[match(colnames(dat1), info_curr_df$id),]
        res = voom_sub(dat1, info_curr_df, predictors, contrasts, normalize.method, lib.size=norm_vals, blocking=blocking)
      } else {
        #dge = DGEList(counts=dat1, samples=info_curr_df)
        #norm_vals = getOffset(dge)
        res = voom_sub(dat1, info_curr_df, predictors, contrasts, normalize.method=normalize.method, lib.size=norm_vals, blocking=blocking)        
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

voom_two_levels_boot = function(dat, info, predictors, factors, contrasts, outdir, model_fname, blocking=NULL, 
                           normalize.method="none", reps=100, randomize=F, ncores=4) {
  
  if (!all(predictors %in% colnames(info)))
    stop("ERROR: Not all predictors in info data.frame.")
  
  prefix1 = paste(outdir, model_fname, sep="/")
  dir.create(prefix1, recursive=T)
  
  to_select = c(predictors, factors, "id", "tags")
  info = info %>% select(one_of(to_select))
  info_df = as.data.frame(info)
  levels1 = levels(info_df[,factors[1]])
  levels2 = levels(info_df[,factors[2]])
  
  fname = paste(prefix1, "results.rds", sep="/")
  fname_pred = paste(prefix1, "predictors.txt", sep="/")
  
  predictors1 = as.formula(paste("~", paste(predictors, collapse="+"), sep=""))
  
  write.table(predictors, fname_pred)
  
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
      
      curr1 = mclapply(1:reps, function(r) {
      #curr1 = lapply(1:reps, function(r) {
        
        if (!randomize) {
          
          info_curr1 = info_curr %>% 
            group_by_(lazyeval::interp(~c, c=as.name(predictors[1]))) %>% 
            sample_n(nrow(.), replace=T)
        } else {
          info_curr1 = info_curr %>% select(tags, id) %>% 
            mutate(tags = sample(tags, length(tags), replace=F)) %>% 
            inner_join(info_curr %>% select(-id), by="tags")
        }
        info_curr_df = as.data.frame(info_curr1)
        dat1 = dat[,match(info_curr1$id, colnames(dat))]
        
        res = voom_sub(dat1, info_curr1, predictors1, contrasts, normalize.method, blocking)
        res[,factors[2]] = l2
        res$rep = r
        return(res)
      }, mc.cores=ncores)
      #})
      bind_rows(curr1)
    })
    curr = curr[unlist(lapply(curr, function(x) !is.null(x)))]
    curr = bind_rows(curr)
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

voom_two_levels_pool = function(dat, info, predictors, factors, contrasts, outdir, model_fname, blocking=NULL, 
                                normalize.method="none", reps=100, randomize=F, ncores=4) {
  
  if (!all(predictors %in% colnames(info)))
    stop("ERROR: Not all predictors in info data.frame.")
  
  prefix1 = paste(outdir, model_fname, sep="/")
  dir.create(prefix1, recursive=T)
  
  to_select = c(predictors, factors, "id", "tags")
  info = info %>% select(one_of(to_select))
  info_df = as.data.frame(info)
  levels1 = levels(info_df[,factors[1]])
  levels2 = levels(info_df[,factors[2]])
  
  fname = paste(prefix1, "results.rds", sep="/")
  fname_pred = paste(prefix1, "predictors.txt", sep="/")
  
  predictors1 = as.formula(paste("~", paste(predictors, collapse="+"), sep=""))
  
  write.table(predictors, fname_pred)
  
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
      
      curr1 = mclapply(1:reps, function(r) {
        #curr1 = lapply(1:reps, function(r) {
        
        if (!randomize) {
          
          info_curr1 = info_curr %>% group_by_(lazyeval::interp(~c, c=as.name(predictors[1]))) %>% 
            sample_n(nrow(.), replace=T)
        } else {
          info_curr1 = info_curr %>% select(tags, id) %>% 
            mutate(tags = sample(tags, length(tags), replace=F)) %>% 
            inner_join(info_curr %>% select(-id), by="tags")
        }
        info_curr_df = as.data.frame(info_curr1)
        dat1 = dat[,match(info_curr1$id, colnames(dat))]
        
        res = voom_sub(dat1, info_curr1, predictors1, contrasts, normalize.method, blocking)
        res[,factors[2]] = l2
        res$rep = r
        return(res)
      }, mc.cores=ncores)
      #})
      bind_rows(curr1)
    })
    curr = curr[unlist(lapply(curr, function(x) !is.null(x)))]
    curr = bind_rows(curr)
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

plot_top_genes_one_level = function(res, data, factors, contrast, coefficient, thresh, value="voom_mean", include_sem=T, outdir) {
  
  outdir_figures = paste(outdir, "figures", sep="/")
  dir.create(outdir_figures)
  info_df = as.data.frame(res)
  levels1 = levels(info_df[,factors[1]])
  
  # Geoms
  if (include_sem) {
    dots = list(lazyeval::interp(~v, v=as.name(value)),
                lazyeval::interp(~v, v=as.name(paste0(str_replace(value, "_[a-z]+", ""), "_sem"))))
    names(dots) = c("value", "value_sem")
    pr =  geom_pointrange(aes(y=value,
                              ymin=(value - value_sem), 
                              ymax=(value + value_sem)),
                          position=position_dodge(width=.5), size=.2)
  } else {
    dots = list(lazyeval::interp(~v, v=as.name(value)))
    names(dots) = c("value")
    pr =  geom_point(aes(y=value), alpha=.5,
                     position=position_jitterdodge(dodge.width=.5, jitter.width=.2), size=.2)
  }
  data = data %>% rename_(.dots=dots)
  # pr =  geom_pointrange(aes(y=value,
  #                           ymin=(value - value_sem), 
  #                           ymax=(value + value_sem)),
  #                       position=position_dodge(width=.5), size=.2)
  col = scale_color_manual(values=c("black", "red"))
  th = theme(strip.text.y = element_text(angle=0))
  for (l1 in levels1) {
    print(l1)

      dots = lazyeval::interp(~f1==l1 & coef==contrast, 
                              .values=list(f1=as.name(factors[1])))
      
      res_cur = ungroup(res) %>% 
        filter_(.dots=dots) %>% 
        filter(adj.P.Val<thresh) %>%
        top_n(-9, adj.P.Val)
      toplot = res_cur$gene_id
      if (length(toplot) == 0)
        next()
      
      tmp = data %>% filter(gene_id %in% toplot)
      gg = ggplot(tmp , aes_string(factors[1], color=coefficient))
      gg = gg + pr
      gg = gg + col
      gg = gg + th
      gg = gg + facet_grid(gene_id~., scales="free_y")
      gg
      ggsave(paste(outdir_figures, sprintf("topgenes_%s_by_position.pdf", l1), sep="/"), width=10, height=7.5)
      
      # tmp_res = res %>% filter(gene_id %in% toplot, coef==contrast)
      # gg = ggplot(tmp_res, 
      #             aes(position, logFC, ymin=CI.L, ymax=CI.R, color=position))
      # gg = gg + geom_hline(yintercept=0, linetype=1)
      # gg = gg + geom_pointrange(size=.2)
      # gg = gg + facet_grid(gene_id~duration.of.experiment)
      # gg = gg + scale_color_manual(values=position_colors)
      # gg = gg + theme_bw()
      # 
      # gg
      # ggsave(paste(outdir_figures, sprintf("topgenes_%s_coef_by_duration.pdf", l1), sep="/"), width=6, height=7.5)      
  }
}


plot_top_genes_two_levels = function(res, data, factors, contrast, coefficient, thresh, value="voom_mean", include_sem=T, outdir) {
  
  outdir_figures = paste(outdir, "figures", sep="/")
  dir.create(outdir_figures)
  info_df = as.data.frame(res)
  levels1 = levels(info_df[,factors[1]])
  levels2 = levels(info_df[,factors[2]])
  
  # Geoms
  if (include_sem) {
    
    
    dots = list(lazyeval::interp(~v, v=as.name(value)),
                lazyeval::interp(~v, v=as.name(paste0(str_replace(value, "_[a-z]+", ""), "_sem"))))
    names(dots) = c("value", "value_sem")
    pr =  geom_pointrange(aes(y=value,
                              ymin=(value - value_sem), 
                              ymax=(value + value_sem)),
                          position=position_dodge(width=.5), size=.2)
  } else {
    dots = list(lazyeval::interp(~v, v=as.name(value)))
    names(dots) = c("value")
    pr =  geom_point(aes(y=value), alpha=.5,
                     position=position_jitterdodge(dodge.width=.5, jitter.width=.2), size=.2)
  }
  data = data %>% rename_(.dots=dots)
  
  #lin = geom_line(aes_string(y="value", group="set", color=coefficient), position=position_dodge(width=.5))
  col = scale_color_manual(values=c("black", "red"))
  th = theme(strip.text.y = element_text(angle=0))
  for (l1 in levels1) {
    print(l1)
    for (l2 in levels2) {
      print(l2)
      dots = lazyeval::interp(~f1==l1 & f2==l2 & coef==contrast, 
                    .values=list(f1=as.name(factors[1]),
                                 f2=as.name(factors[2])))
      
      res_cur = ungroup(res) %>% 
        filter_(.dots=dots) %>% 
        filter(adj.P.Val<thresh) %>%
        top_n(-9, adj.P.Val)
      toplot = res_cur$gene_id
      if (length(toplot) == 0)
        next()
      tmp = data %>% filter(gene_id %in% toplot)
      gg = ggplot(tmp , aes_string(factors[1], color=coefficient))
      gg = gg + pr
      gg = gg + col
      #gg = gg + lin
      gg = gg + th
      gg = gg + facet_grid(as.formula(sprintf("gene_id~%s", factors[2])), scales="free_y")
      gg
      ggsave(paste(outdir_figures, sprintf("topgenes_%s_%s_by_%s.pdf", l1, l2, factors[2]), sep="/"), width=10, height=7.5)
      
      gg = ggplot(tmp, aes_string(factors[2], color=coefficient))
      gg = gg + pr
      gg = gg + col
      gg = gg + th
      gg = gg + facet_grid(as.formula(sprintf("gene_id~%s", factors[1])), scales="free_y")
      gg
      ggsave(paste(outdir_figures, sprintf("topgenes_%s_%s_by_%s.pdf", l1, l2, factors[1]), sep="/"), width=10, height=7.5)
      # 
      # tmp_res = res %>% filter(gene_id %in% toplot, coef==contrast)
      # gg = ggplot(tmp_res, 
      #             aes(position, logFC, ymin=CI.L, ymax=CI.R, color=position))
      # gg = gg + geom_hline(yintercept=0, linetype=1)
      # gg = gg + geom_pointrange(size=.2)
      # gg = gg + facet_grid(gene_id~duration.of.experiment)
      # gg = gg + scale_color_manual(values=position_colors)
      # gg = gg + theme_bw()
    
      gg
      ggsave(paste(outdir_figures, sprintf("topgenes_%s_%s_coef_by_duration.pdf", l1, l2), sep="/"), width=6, height=7.5)      
    }
  }
}

plot_foldchange_two_levels = function(res, res2=NULL, factors, contrast, coefficient, thresh, value="voom_mean", outdir) {
  outdir_figures = paste(outdir, "figures", sep="/")
  dir.create(outdir_figures)
  info_df = as.data.frame(res)
  levels1 = levels(info_df[,factors[1]])
  levels2 = levels(info_df[,factors[2]])
  
  for (l1 in levels1) {
    print(l1)
    for (l2 in levels2) {
      print(l2)
      dots = lazyeval::interp(~f1==l1 & f2==l2 & coef==contrast, 
                    .values=list(f1=as.name(factors[1]),
                                 f2=as.name(factors[2])))
      
      if (!is.null(res2)) {
        res_cur = ungroup(res2) %>% 
          filter_(.dots=dots) %>%
          top_n(-9, adj.P.Val_mean)
          #arrange(desc(nsig)) %>% 
          #slice(1:9)
      } else {
        res_cur = ungroup(res) %>% 
          filter_(.dots=dots) %>% 
          filter(adj.P.Val<thresh) %>%
          top_n(-9, adj.P.Val)
      }

      toplot = res_cur$gene_id
      if (length(toplot) == 0)
        next()
      alpha_factor = 5/100 
      #toplot = ungroup(res2) %>% filter(position=="lman", duration.of.experiment=="4", coef=="procedure2deaf") %>% top_n(-12, logFC_mean)
      gg = ggplot(res %>% filter(gene_id %in% toplot, coef=="procedure2deaf"), aes(position, logFC))
      gg = gg + geom_hline(yintercept=0, linetype=3)
      gg = gg + geom_point(alpha=alpha_factor)
      gg = gg + facet_grid(gene_id~duration.of.experiment)
      gg = gg + theme(strip.text.y= element_text(angle=0))
      gg = gg + ylim(-5, 5)
      gg
      ggsave(paste(outdir_figures, sprintf("logFC_%s_%s_coef_by_duration.pdf", l1, l2), sep="/"), width=6, height=7.5)
    }
  }
}

plot_num_deg_genes_bar = function(deg_df, outdir, coefficients, factors, colors) {
  plot_dir = paste(outdir, "figures", sep="/")
  dir.create(plot_dir)
  
  for (co in coefficients) {
    print(co)
    tmp = deg_df %>% filter(coef==co)
    if (nrow(tmp) < 1)
      next()
    print(factors)
    if (length(factors)==2) {
      
      map(seq_along(factors), function(f) {
        gg = ggplot(tmp, aes_string(factors[f], "nsig", fill=factors[f])) + geom_bar(stat="identity")
        gg = gg + theme_bw()
        gg = gg + scale_fill_manual(values=colors[[f]])
        gg = gg + theme(axis.text.x = element_text(angle=45, hjust=1))
        gg = gg + facet_grid(as.formula(sprintf(".~%s", factors[setdiff(1:length(factors), f)])))
        gg = gg + labs(y="Number of genes", x="Duration")
        gg
        ggsave(paste(plot_dir, sprintf("ngenes_%s_bar_by_%s.pdf", co, factors[f]), sep="/"), width=8, height=6)
    })
    } 
  }
}

plot_heatmap_threshed_split = function(res, info, factors, thresh, contrast, coefficient, value, colors, outdir, lims, deg_un=NULL, factors2=NULL) {
  plot_dir = paste(outdir, "figures", sep="/")

  
  if (is.null(factors2)) {
    
    info_df = as.data.frame(res)
  levels1 = levels(info_df[,factors[1]])
  levels2 = levels(info_df[,factors[2]])
  } else {
    info_df = as.data.frame(deg_un)
    levels1 = levels(info_df[,factors2[1]])
    levels2 = levels(info_df[,factors2[2]])
}
  for (l1 in levels1) {
    print(l1)
    for (l2 in levels2) {
      print(l2)
      plot_dir2 = paste(plot_dir, sprintf("%s_%s", l1, l2), sep="/")
      dir.create(plot_dir2)

      if (!is.null(factors2)) {
        dots = lazyeval::interp(~f1==l1 & f2==l2 & coef==contrast, 
                                .values=list(f1=as.name(factors2[1]),
                                             f2=as.name(factors2[2])))
      } else {
        dots = lazyeval::interp(~f1==l1 & f2==l2 & coef==contrast, 
                                .values=list(f1=as.name(factors[1]),
                                             f2=as.name(factors[2])))
      }
      if (!is.null(deg_un)) {

        res1 = deg_un
      } else {

        res1 = res
      }
      res_cur = ungroup(res1) %>% 
        filter_(.dots=dots) %>% 
        filter(adj.P.Val<thresh)
      if (nrow(res_cur)<10)
        next()
      plot_heatmap_threshed(res, info, factors, thresh, contrast, value, lims=lims, colors=colors, outdir=plot_dir2, deg_un=res_cur)
    }
  }
}

plot_heatmap_threshed = function(res, info, factors, thresh, coefficient, value, lims=c(-5,5), colors, outdir, deg_un=NULL) {
  
  plot_dir = paste(outdir, "figures", sep="/")
  dir.create(plot_dir, recursive=T)
  
  dotsm = list(lazyeval::interp(~paste(f1, f2, sep="_"),
                                f1=as.name(factors[1]), 
                                f2=as.name(factors[2])))
  names(dotsm) = "condition"
  info_res = info %>% 
    distinct_(.dots=list(lazyeval::interp(~v, v=as.name(factors[1])),
                   lazyeval::interp(~v, v=as.name(factors[2])))) %>%
  arrange_(.dots=list(lazyeval::interp(~v, v=as.name(factors[1])),
                        lazyeval::interp(~v, v=as.name(factors[2])))) %>%
  mutate_(.dots=dotsm)
  
  res = res %>% 
    mutate(sig = adj.P.Val<thresh) %>%
    mutate_(.dots=dotsm)
  if (is.null(deg_un)) {
    deg_un = res %>% filter(sig, coef==coefficient) %>% distinct(gene_id)
  }
  print(nrow(deg_un))
    
  res = res %>% filter(coef==coefficient, gene_id %in% deg_un$gene_id)
  missing = setdiff(info_res$condition, unique(res$condition))
  if (length(missing)>0) {
    to_add = res[1:length(missing),]
    to_add = to_add %>% mutate(condition = missing) %>% 
      mutate(factor1 = unlist(lapply(str_split(condition, "_"), function(x) x[1])),
             factor2 = unlist(lapply(str_split(condition, "_"), function(x) x[2])),
             sig=FALSE,
             coef=coefficient)
    
    to_add = to_add %>% select(-one_of(factors))
    dots = list(lazyeval::interp(~factor1),
                lazyeval::interp(~factor2))
    names(dots) = c(factors[1], factors[2])
    to_add = to_add %>% rename_(.dots=dots)
    res = rbind(res, to_add)
  }
  tmp = res
  #tmp = res %>% filter(coef==coefficient, gene_id %in% deg_un$gene_id)
  
  resa = acast(tmp, as.formula(sprintf("gene_id~%s+%s", factors[1], factors[2])), value.var=value, fill = 0)
  aheatmap(resa, Colv=T, color="-RdBu:100", breaks = seq(lims[1], lims[2], length.out=101), hclustfun = "ward",
           annCol=info_res,
           annColors=colors,
           labCol=NA, 
           filename = paste(plot_dir, sprintf("heatmap_%s_by_%s_clust_thresh%s.pdf", value, factors[1], thresh), sep="/"),
           width=10,
           height=8,
           info=T)
  aheatmap(resa, Colv=NA, color="-RdBu:100", breaks = seq(lims[1], lims[2], length.out=101), hclustfun = "ward",
           annCol=info_res,
           annColors=colors,
           labCol=NA, 
           filename = paste(plot_dir, sprintf("heatmap_%s_by_%s_thresh%s.pdf", value, factors[1], thresh), sep="/"),
           width=10,
           height=8,
           info=T)
  
  ## Threshold for significance
  sig = acast(tmp,as.formula(sprintf("gene_id~%s+%s", factors[1], factors[2])), value.var="sig")
  resa[!sig]  = 0
  aheatmap(resa, Colv=T, color="-RdBu:100", breaks = seq(lims[1], lims[2], length.out=101), hclustfun = "ward",
           annCol=info_res,
           annColors=colors,
           labCol=NA, 
           filename = paste(plot_dir, sprintf("heatmap_%s_thresh%sfilter_by_%s.pdf", value, thresh, factors[1]), sep="/"),
           width=10,
           height=8,
           info=T)
  
  # Reverse order
  info_res = info %>% 
    distinct_(.dots=list(lazyeval::interp(~v, v=as.name(factors[1])),
                         lazyeval::interp(~v, v=as.name(factors[2])))) %>%
    arrange_(.dots=list(lazyeval::interp(~v, v=as.name(factors[2])),
                        lazyeval::interp(~v, v=as.name(factors[1])))) %>%
    select_(.dots=list(lazyeval::interp(~v, v=as.name(factors[2])),
                       lazyeval::interp(~v, v=as.name(factors[1]))))
  resa = acast(tmp,  as.formula(sprintf("gene_id~%s+%s", factors[2], factors[1])), value.var=value)
  aheatmap(resa, Colv=NA, color="-RdBu:100", breaks = seq(lims[1], lims[2], length.out=101), hclustfun = "ward",
           annCol=info_res,
           annColors=rev(colors),
          labCol=NA, 
           filename = paste(plot_dir, sprintf("heatmap_%s_by_%s.pdf", value, factors[2]), sep="/"),
           width=10,
           height=8,
           info=T)
  
  sig = acast(tmp,  as.formula(sprintf("gene_id~%s+%s", factors[2], factors[1])), value.var="sig")
  resa[!sig]  = 0
  aheatmap(resa, Colv=NA, color="-RdBu:100", breaks = seq(lims[1], lims[2], length.out=101), hclustfun = "ward",
           annCol=info_res,
           annColors=rev(colors),
           labCol=NA, 
           filename = paste(plot_dir, sprintf("heatmap_%s_thresh%sfilter_by_%s.pdf", value, thresh, factors[2]), sep="/"),
           width=10,
           height=8,
           info=T)
}

calc_number_sig_genes = function(res, factors, thresh=.1) {
  res1 = res %>% filter(adj.P.Val<thresh)
  # res1_sum = res1 %>%
  #   group_by_(.dots=list(lazyeval::interp(~coef), 
  #                        lazyeval::interp(~v, v=as.name(factors[1])),
  #                        lazyeval::interp(~v, v=as.name(factors[2])))) %>% 
  #   summarize(nsig=n())
  res1_sum = res1 %>%
    group_by(coef)
  for (i in 1:length(factors)) {
    res1_sum = res1_sum %>%
      group_by_(.dots=list(lazyeval::interp(~v, v=as.name(factors[i]))), add=T)
  }
  res1_sum = res1_sum %>% 
    summarize(nsig=n())
  return(res1_sum)
}

generate_deg_plots_one_level = function(res, data, info, thresh, factors, contrasts, coefficient, factor_colors, plot_value="voom",
                                        include_sem = F,
                                        p_value="adj.P.Val",
                                        fc_value="logFC",
                                        t_value="t", 
                                        outdir1outdir1) {
  if (!dir.exists(outdir1))
    dir.create(outdir1, recursive=T)
  
  dots = list(lazyeval::interp(~v, v=as.name(p_value)),
              lazyeval::interp(~v, v=as.name(fc_value)),
              lazyeval::interp(~v, v=as.name(t_value)))
  names(dots) = c("adj.P.Val", "logFC", "t")
  res = res %>% mutate_(.dots=dots)
  
  res1 = res %>% filter(adj.P.Val<thresh)
  res1_sum = res1 %>%
    group_by_(.dots=list(lazyeval::interp(~coef),
                         lazyeval::interp(~v, v=as.name(factors[1])))) %>%
    summarize(nsig=n())
  res1_sum = calc_number_sig_genes(res, factors, thresh)
  
  vp.page <- viewport(x = 0.5, y = 0.5,
                      width = unit(x = 8.5, units = "inches"),
                      height = unit(x = 11, units = "inches"))
  vp.marg <- viewport(x = 0.5, y = 0.5,
                      width = (7.5 / 8.5), height = (10 / 11))
  col1 <- tableGrob(res1_sum, rows = NULL)
  cairo_pdf(paste(outdir1, "nsig_table.pdf", sep="/"))
  pushViewport(vp.page); pushViewport(vp.marg)
  grid.draw(col1)
  #grid.table(res1_sum)
  dev.off()
  
  write.table(res1_sum, paste(outdir1, sprintf("nsig_at_%s.txt", thresh), sep="/"))
  
  cat("--Plotting bar graphs\n")
  plot_num_deg_genes_bar(res1_sum, outdir1, coefficients=contrasts, factors=factors, colors=factor_colors)
  cat("--Plotting gene plots\n")
  plot_top_genes_one_level(res, data,
                            factors=factors,
                            contrast=contrasts[1],
                            coefficient=coefficient,
                            thresh=.1,
                            value=plot_value,
                           include_sem = include_sem,
                            outdir1)
  # 
  cat("--Plotting heatmaps")
  plot_heatmap_threshed(res, info, 
                        factors, thresh=.1, 
                        coefficient=contrasts[1], 
                        value="logFC",
                        lims=c(-3, 3), 
                        colors=factor_colors,
                        outdir = outdir1)
  plot_heatmap_threshed(res, info, 
                        factors, thresh=.1, 
                        coefficient=contrasts[1], 
                        value="t",
                        lims=c(-5, 5), 
                        colors=factor_colors,
                        outdir = outdir1)
}


generate_deg_plots = function(res, data, info, thresh, 
                              factors, contrasts, coefficient, 
                              factor_colors, plot_value="voom", 
                              include_sem=TRUE, 
                              p_value="adj.P.Val",
                              fc_value="logFC",
                              t_value="t", 
                              outdir1
                              ) {
  if (!dir.exists(outdir1))
    dir.create(outdir1, recursive=T)
  
  dots = list(lazyeval::interp(~v, v=as.name(p_value)),
              lazyeval::interp(~v, v=as.name(fc_value)),
              lazyeval::interp(~v, v=as.name(t_value)))
  names(dots) = c("adj.P.Val", "logFC", "t")
  res = res %>% mutate_(.dots=dots)
  res1 = res %>% filter(adj.P.Val<thresh)
  res1_sum = res1 %>%
    group_by_(.dots=list(lazyeval::interp(~coef),
                         lazyeval::interp(~v, v=as.name(factors[1])),
                         lazyeval::interp(~v, v=as.name(factors[2])))) %>%
    summarize(nsig=n())
  res1_sum = calc_number_sig_genes(res, factors, thresh)

  cairo_pdf(paste(outdir1, "nsig_table.pdf", sep="/"))
  grid.table(res1_sum)
  dev.off()

  write.table(res1_sum, paste(outdir1, sprintf("nsig_at_%s.txt", thresh), sep="/"))

  cat("--Plotting bar graphs\n")
  plot_num_deg_genes_bar(res1_sum, outdir1, coefficients=contrasts, colors=factor_colors, factors=factors)
  
  cat("--Plotting gene plots\n")
  plot_top_genes_two_levels(res, data,
                            factors=factors,
                              contrast=contrasts[1],
                            coefficient=coefficient,
                            thresh=.1,
                            value=plot_value,
                            include_sem=include_sem,
                            outdir1)

  cat("--Plotting heatmaps")
  plot_heatmap_threshed(res, info, 
                        factors, thresh=.1, 
                        coefficient=contrasts[1], 
                        value="logFC",
                        lims=c(-1, 1), 
                        colors=factor_colors,
                        outdir = outdir1)
  plot_heatmap_threshed(res, info, 
                        factors, thresh=.1, 
                        coefficient=contrasts[1], 
                        value="t",
                        lims=c(-5, 5), 
                        colors=factor_colors,
                        outdir = outdir1)
}

write_out_deg_genes = function(res, full, factors, outdir, prefix ) {
  gene_dir = paste(outdir, "genes", sep="/")
  dir.create(gene_dir, recursive=T)

  conditions = unique(res$condition)
  
  for (con in conditions) {
    print(con)
    fname = paste(gene_dir, sprintf("%s_%s.txt", prefix, con), sep="/")
    res_cur = res$gene_id[res$condition==con]
    cat(sprintf("Num. genes: %s\n", length(res_cur)))
    g = process_gene_ids(res_cur)
    g = gene_symbol_to_entrez(g)
    write.table(g, fname, quote=F, sep="\t", row.names=F, col.names=F)
  }
  
  genes = gene_symbol_to_entrez(process_gene_ids(unique(full$gene_id)))
  write.table(genes, paste(gene_dir, sprintf("%s_total_genes.txt", prefix), sep="/"), quote=F, sep="\t", row.names=F, col.names=F)
}


