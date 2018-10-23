library(scde)

get_go_env = function(dat) {
  library(org.Hs.eg.db)
  # translate gene names to ids
  genes = rownames(dat)
  genes = str_replace(genes, "\\*", "")
  genes = str_replace(genes, "\\.[0-9]+", "")
  ids = unlist(lapply(mget(genes, org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
  rids = names(ids); names(rids) = ids 
  # convert GO lists from ids to gene names
  gos.interest = unique(c(ls(org.Hs.egGO2ALLEGS))) 
  go.env = lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
  go.env = clean.gos(go.env) # remove GOs with too few or too many genes
  go.env = list2env(go.env) # convert to an environment
  return(go.env)
}

knn_sub = function(dat, groups=NULL, subpops=2, outdir) {
  if (is.null(groups)) {
  knn = knn.error.models(dat, 
                         k = ncol(dat)/subpops, 
                         n.cores = 2, 
                         min.count.threshold = 2, 
                         min.nonfailed = 5, 
                         max.model.plots = 10,
                         save.model.plots = T)
  } else {
    knn = knn.error.models(dat, 
                           groups=groups, 
                           k = ncol(dat)/subpops, 
                           n.cores = 2, 
                           min.count.threshold = 2, 
                           min.nonfailed = 5, 
                           max.model.plots = 10,
                           save.model.plots = T)
  }
  saveRDS(knn, paste(outdir, "error_model.rds", sep="/"))
  return(knn)
}

varnorm_sub = function(knn, dat, batch=NULL, outdir) {
  fname = paste(outdir, "varinfo.rds", sep="/")
  varinfo = pagoda.varnorm(knn, 
                           counts = dat,
                           batch = batch,
                           trim = 3/ncol(dat), 
                           max.adj.var = 5, 
                           n.cores = 2, 
                           plot = TRUE)
  saveRDS(varinfo, fname)
  return(varinfo)
}

wpca_sub = function(varinfo, go.env, outdir) {
  fname = paste(outdir, "wpca.rds", sep="/")
  pwpca = pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 2)
  saveRDS(pwpca, fname)
  return(pwpca)
}

gene_clusters_sub = function(varinfo, outdir) {
  fname = paste(outdir, "gene_clusters.rds", sep="/")
  clpca = pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = 2, plot = TRUE)
  saveRDS(clpca, fname)
  return(clpca)
}

scde_routine = function(dat, groups, subpops, batch, go.env, outdir) {
  
  cat("KNN\n")
  knn = knn_sub(dat=dat, 
               # groups=groups, 
                subpops=subpops,
                outdir=outdir)
  
  cat("varnorm\n")
  varinfo = varnorm_sub(knn = knn, dat = dat, batch = batch, outdir = outdir)
  varinfo = pagoda.subtract.aspect(varinfo, colSums(dat[, rownames(knn)]>0))
  cat("wpca\n")
  pwpca = wpca_sub(varinfo = varinfo, go.env = go.env, outdir = outdir)
  cat("clpca\n")
  clpca = gene_clusters_sub(varinfo = varinfo, outdir = outdir)
  return(list(knn=knn, varinfo=varinfo, pwpca=pwpca, clpca=clpca))
}

pagoda_two_levels = function(dat, info, go.env=NULL, factors, group, batch, idcol, outdir) {
  
  if (is.null(go.env))
    go.env = get_go_env(dat)
  
  #prefix1 = paste(outdir, model_fn, sep="/")
  dir.create(outdir, recursive=T)
  
  info_df = as.data.frame(info)
  levels1 = levels(info_df[,factors[1]])
  levels2 = levels(info_df[,factors[2]])
  
  res = lapply((1:length(levels1)), function(i1) {
    print(levels1[i1])
    l1 = levels1[i1]
    curr = lapply(1:length(levels2), function(i2) {
      l2 = as.character(levels2[i2])
      print(l2)
      
      outdir1 = paste(outdir, sprintf("%s_%s", l1, l2), sep="/")
      dir.create(outdir1, recursive=T)
      setwd(outdir1)
      ### Setup metadata  --------------------------------------
      dots = lazyeval::interp(~f1==l1 & f2==l2, .values=list(f1=as.name(factors[1]),
                                                             f2=as.name(factors[2])))
      info_curr = droplevels(info %>% filter_(.dots=dots ))
      info_curr_df = as.data.frame(info_curr)
      dat_cur = dat[,match(info_curr_df[,idcol], colnames(dat))]
      
      groups = info_curr_df[,group]
      batch = info_curr_df[,batch]
      scde_routine(dat_cur, groups = groups, 
                   subpops = length(levels(groups)), 
                   batch = batch, 
                   go.env = go.env, 
                   outdir = outdir1)
    })
  })
}

pagoda_pathways = function(info)
pagoda_reduce = function(info, corr.power=4, color_var=NULL, outdir) {
  varinfo = readRDS(paste(outdir, "varinfo.rds", sep="/"))
  pwpca = readRDS(paste(outdir, "wpca.rds", sep="/"))
  clpca = readRDS(paste(outdir, "gene_clusters.rds", sep="/"))
  # get full info on the top aspects
  tam = pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
  # determine overall cell clustering
  hc = pagoda.cluster.cells(tam, varinfo, method="ward")
  tamr = pagoda.reduce.loading.redundancy(tam, pwpca, clpca, corr.power=corr.power)
  
  plot_dir = paste(outdir, "pathways", sep="/")
  dir.create(plot_dir)
  fname = paste(plot_dir, sprintf("overall_%s.pdf", color_var), sep="/")
  cairo_pdf(fname)
  
  if (is.null(color_var)) {
  tamr2 = pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, 
                                   #labRow = NA, 
                                   labCol = NA, 
                                   #col.cols = info[,color_var],
                                   #col.cols= info_cur[match(hc$labels, info_cur$id), color_var],
                                   #col.cols = info_cur$deafcol[hc$order],
                                   box = TRUE, margins = c(1, 20), trim = 0)
  } else {
    tamr2 = pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, 
                                     #labRow = NA, 
                                     labCol = NA, 
                                     #col.cols = as.character(info[,color_var]),
                                     #col.cols= info[match(hc$labels, info$id), color_var],
                                     col.cols = info[hc$order, color_var],
                                     box = TRUE, margins = c(1, 20), trim = 0) 
  }
  dev.off()
  return(tamr2$df)
}

pagoda_reduce_two_levels = function(info, factors, corr.power=4, color_var=NULL, outdir) {
  general_two_levels(info, 
                     factors, 
                     outdir=outdir, 
                     fnc="pagoda_reduce", 
                     args =list(corr.power=corr.power, color_var=color_var))
}

general_two_levels = function(info, factors, outdir, fnc, args) {
  
  info_df = as.data.frame(info)
  levels1 = levels(info_df[,factors[1]])
  levels2 = levels(info_df[,factors[2]])
  
  res = lapply((1:length(levels1)), function(i1) {
    print(levels1[i1])
    l1 = levels1[i1]
    res1 = lapply(1:length(levels2), function(i2) {
      l2 = as.character(levels2[i2])
      print(l2)
      
      outdir1 = paste(outdir, sprintf("%s_%s", l1, l2), sep="/")
      dir.create(outdir1, recursive=T)
      setwd(outdir1)
      ### Setup metadata  --------------------------------------
      dots = lazyeval::interp(~f1==l1 & f2==l2, .values=list(f1=as.name(factors[1]),
                                                             f2=as.name(factors[2])))
      info_curr = droplevels(info %>% filter_(.dots=dots ))
      info_curr_df = as.data.frame(info_curr)
      tryNULL({
        do.call(fnc, c(list(info=info_curr_df, outdir=outdir1), args))
      })
    })
    names(res1) = levels2
    res1 = res1[unlist(lapply(res1, function(x) !is.null(x)))]
    bind_rows(res1, .id=factors[2])
  }) 
  names(res) = levels1
  bind_rows(res, .id=factors[1])
}

scde_weight_corr = function(model, dat, mags) {
  # load boot package for the weighted correlation implementation
  require(boot)
  k = 0.95
  ids = colnames(dat)
  datl = log10(dat+1)
  reciprocal.dist = stats::as.dist(1 - do.call(rbind, lapply(ids, function(nam1) {

    unlist(lapply(ids, function(nam2) {
      # reciprocal probabilities
      # TODO figure how to flip this
      f1 = scde.failure.probability(models = model[nam1,,drop = FALSE], magnitudes = mags[, nam2])
      f2 = scde.failure.probability(models = model[nam2,,drop = FALSE], magnitudes = mags[, nam1])
      # weight factor
      pnf = sqrt((1-f1)*(1-f2))*k +(1-k); 
      #boot::corr(log10(cbind(dat[, nam1], dat[, nam2])+1), w = pnf)
      boot::corr(datl[,c(nam1, nam2)], w=pnf)
    }))
  })),
 # mc.cores = 6)), 
  upper = FALSE)
  recm = as.matrix(reciprocal.dist)
  dimnames(recm) = list(ids, ids)
  return(recm)
}

scde_precompute_reciprocal_probs = function(model, dat, mags) {
  ids = colnames(dat)
  # precompute failure prob matrix
  range1 = seq_len(length(ids)-1)
  res1 = lapply(range1, function(i1) {
    #ids2 = ids[-which(ids==nam1)]
    range2 = seq(i1+1,length(ids))
    res2 = lapply(range2, function(i2) {
      # reciprocal probabilities
      # TODO figure how to flip this
      f1 = scde.failure.probability(models = model[ids[i1],,drop = FALSE], magnitudes = mags[, ids[i2]])
      f2 = scde.failure.probability(models = model[ids[i2],,drop = FALSE], magnitudes = mags[, ids[i1]])
      res = data.frame(f1, f2)
      colnames(res) = c("p1", "p2")
      res$gene_id = rownames(f1)
      return(res)
    })
    names(res2) = ids[range2]
    res2 = bind_rows(res2, .id="id2")
    res2
  })
  names(res1) = ids[range1]
  res1 = bind_rows(res1, .id="id1")
  return(res1)
}

scde_precompute_failure_probs = function(model, dat, mags) {
  ids = colnames(dat)
  # precompute failure prob matrix
  mat = matrix(nrow=nrow(dat), ncol=ncol(dat), dimnames = dimnames(dat))
  range1 = seq_len(length(ids))
  for (i1 in range1) {
  #res1 = lapplyrange1, function(i1) {
      f1 = scde.failure.probability(models = model[ids[i1],,drop = FALSE], magnitudes = mags[, ids[i1]])
      mat[,i1] = f1[,1]
      #return(mat)
      #res = data.frame(f1)
      #colnames(res) = c("prob.f")
      #res$gene_id = rownames(f1)
      #return(res)
  #})
}
  return(mat)
  #names(res1) = ids[range1]
  #res1 = bind_rows(res1, .id="id")
  #return(res1)
}

scde_weight_corr_by_gene = function(dat, mags, probs) {
  k = 0.95
  datt = t(dat)
  dattl = log10(datt+1)
  magt = t(mags)
  magt[is.infinite(magt)] = NA
  genes = colnames(datt)
  res1 = mclapply(genes, function(nam1) {
    #genes1 = genes[-which(genes==nam1)]
    unlist(lapply(genes, function(nam2) {
      #probs_cur = probs %>% filter(gene_id %in% c(nam1, nam2))
      f1 = probs[nam1,]
      f2 = probs[nam2,]
      pnf = sqrt((1-f1)*(1-f2))*k +(1-k)
      d = na.omit(magt[,c(nam1, nam2)])
      if (nrow(d) < (nrow(magt)*.5))
        return(NA)
      pnf = pnf[names(pnf) %in% rownames(d)]
      #boot::corr(log10(cbind(dat[, nam1], dat[, nam2])+1), w = pnf)
      boot::corr(d, w=pnf)
    }))
  #}) 
  }, mc.cores=4)
  res1 = do.call(rbind, res1)
  dimnames(res1) = list(genes, genes)
  res1
}