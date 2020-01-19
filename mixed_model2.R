options(stringsAsFactors=F)
library(tidyverse)
library(reshape2)
library(foreach)
library(itertools)
library(doMC)
registerDoMC(cores = 8)
library(lme4)
library(arm)
library(lmerTest)
library(stringr)
library(gridExtra)
library(lazyeval)
library(rlang)
library(car)

filter_by_var = function(dat, thresh) {
  gene_var = apply(dat, 1, var)
  total = sum(gene_var)
  gene_var = sort(gene_var, decreasing = T) 
  gene_var_cum = cumsum(gene_var) / total
  
  gene_var_cum1 = gene_var_cum[gene_var_cum<=thresh]
  return(names(gene_var_cum1))
}

filter_lme_for_plot = function(coef_data, 
                               ps, 
                               n,
                               dir="abs", 
                               gene_ids=NULL,
                               slope_name="um.from.reference",
                               pvalue_filter=NULL) {
  
  
  if (is.null(gene_ids)) {
    if (is.numeric(pvalue_filter)) {
      ### Filter for genes with FDR less than pvalue_filter
      
      tmp = coef_data %>% filter(coef==slope_name, level=="fixed", reference2_position3==ps)
      p_thresh = tmp %>% 
        filter(coef_type=="Pr(>|t|)") %>%
        mutate(value_fdr = p.adjust(value, method="BH")) %>%
        filter(value_fdr<pvalue_filter)
      tmp = tmp %>% filter(gene_id %in% p_thresh$gene_id)
      
      ### Generate signed p-values
      tp = tmp %>%
        group_by(gene_id) %>% 
        #summarize(n=n())
        mutate(pvalue_signed = sign(value[coef_type=="Estimate"]) * -1 * log10(value[coef_type=="Pr(>|t|)"])) %>%
        mutate(pvalue_signed = if_else(is.infinite(pvalue_signed), -10, pvalue_signed)) %>%
        ungroup() %>%
        distinct(gene_id, pvalue_signed) %>%
        rename(value = pvalue_signed)
    }
    else {
      tp =  coef_data %>% 
        dplyr::filter(reference2_position3==ps, coef==slope_name, coef_type=="Estimate", level=="fixed")
    }
  } else {
    tp = coef_data %>% filter(gene_id %in% gene_ids) %>% 
      dplyr::filter(reference2_position3==ps, coef==slope_name, coef_type=="Estimate", level=="fixed")
  }
  
  tp_genes = tp %>% distinct(gene_id, .keep_all=T)
  #if (is.null(gene_ids)) {
  if (dir=="abs") {
    tp = tp %>%
      top_n(n, abs(value))
  } else if (dir=="up") {
    tp_genes_filt = tp_genes %>% filter(value>0)
    if (nrow(tp_genes_filt) >= n) { 
      tp = tp %>% 
        top_n(n, value)
    } else {
      tp = tp_genes_filt
    }
  } else if (dir=="down") {
    tp_genes_filt = tp_genes %>% filter(value<0)
    if (nrow(tp_genes_filt) >= n) { 
      tp = tp %>%
        top_n(-1 * n, value)
    } else {
      tp = tp_genes_filt
    }
  } else if (dir=="both") {
    tp_up = tp %>% top_n(n/2, value) %>% filter(value>0)
    tp_down = tp %>% top_n(n/2, desc(value)) %>% filter(value<0)
    tp = bind_rows(tp_up, tp_down)
  }
  print(tp)
  return(tp)
  
}

plot_v_pos_rand = function(coef_data, main_data, ps, n, dir="abs", 
                           randoms=1, gene_ids=NULL,
                           slope_name="um.from.reference",
                           y_name="voom1",
                           group_name="brain_side",
                           pvalue_filter=NULL) {
  coef_data = coef_data %>% 
    mutate(gene_id = str_replace(gene_id, "\\.", ",")) %>%
    distinct(coef, coef_type, level, reference2_position3, tags, brain_side, gene_id, .keep_all=T)
  
  tp = filter_lme_for_plot(coef_data = coef_data, 
                           ps = ps, 
                           n = n, 
                           dir = dir, 
                           gene_ids = gene_ids, 
                           slope_name = slope_name, 
                           pvalue_filter = pvalue_filter)
 
  if (nrow(tp)==0)
    return(NULL)
  tp2 = coef_data %>% 
    dplyr::filter(reference2_position3==ps, gene_id %in% tp$gene_id)
  tp2 = tp2 %>% dplyr::rename(brain_side=brain_side)
  coef = NULL
  if (randoms==2) {
    coef = tp2 %>% 
      dplyr::filter(coef_type=="Estimate") %>% 
      #ungroup(.) %>%
      dplyr::group_by(gene_id, coef) %>%
      dplyr::mutate(value_diff = value[level=="fixed"] + value) %>%
      dplyr::filter(level!="fixed") %>%
      dplyr::select(-value) %>%
      spread(coef, value_diff)
    colnames(coef)[colnames(coef)=="(Intercept)"] = "intercept"
    colnames(coef)[colnames(coef)==slope_name] = "slope"
  } else if (randoms==1) {
    inters = tp2 %>% 
      dplyr::filter(coef=="(Intercept)", coef_type=="Estimate") %>% 
      dplyr::group_by(gene_id) %>%
      dplyr::mutate(value_diff = value[level=="fixed"] + value) %>%
      dplyr::filter(level!="fixed")
    
    #dots = list(slope_name, "vst_combat")
    slopes = tp2 %>% 
      dplyr::filter(level=="fixed", coef==slope_name, coef_type=="Estimate")
    coef = left_join(inters[,c("gene_id", "reference2_position3", "value_diff", "brain_side", "tags")], 
                     slopes[,c("gene_id", "reference2_position3", "value")], by=c("gene_id", "reference2_position3"))
    colnames(coef)[grep("value", colnames(coef))] = c("intercept", "slope")
    coef = coef %>% 
      mutate_(.dots=setNames(list(0,0), c(slope_name, "vst_combat")))
  }
  
  
  tmp = NULL
  if (grepl("-", ps)) {
    reg = str_split(ps, "-")[[1]][2]
    ref = str_split(ps, "-")[[1]][1]
    tmp = main_data %>% 
      dplyr::filter(region==reg, reference2==ref, gene_id %in% tp$gene_id)
  } else {
    tmp = main_data %>% 
      dplyr::filter(region==ps, gene_id %in% tp$gene_id)
  }
  #dots = list(interp(~s-min(s), s=as.name(slope_name)))
  #names(dots) = "shifted_slope"
  #tmp = tmp %>% group_by(tags, brain_side) %>% mutate_(.dots=dots)
  
  gg = ggplot(tmp, aes_string(slope_name, y_name, color="reference2_position3")) 
  #gg = ggplot(tmp, aes_string("shifted_slope", y_name, color="reference2_position3")) 
  gg = gg + geom_point(size=.2) 
  #gg = gg + facet_grid(gene_id~tags, scales="free")
  gg = gg + geom_abline(data=coef, aes_string(slope="slope", intercept="intercept"), alpha=1/2)
  gg = gg + facet_grid(gene_id~tags+brain_side, scales="free_y")
  gg = gg + scale_color_manual(values=c("darkblue", "darkred"))
  #gg = gg + scale_y_continuous(expand=c(0,1))
  gg = gg + scale_x_continuous(expand=c(0,200), breaks=seq(round(min(tmp[,slope_name]), -3), 
                                                           round(max(tmp[,slope_name]), -3),
                                                           by=1000))
  #gg = gg + geom_text(aes(label=id))
  if (ps=="posterior-ra") { 
    gg = gg + labs(x="um from posterior surface", y="log2 CPM", title=toupper(ps))
  } else if (ps=="dorsal-ra") {      
    gg = gg + labs(x="um from dorsal surface", y="log2 CPM", title=toupper(ps))
  } else if (grepl("lateral", ps)) {
    gg = gg + labs(x="um from lateral surface", y="log2 CPM", title=toupper(ps))
  }
  gg = gg + theme_classic()
  gg = gg + theme(strip.text.y=element_text(angle=0),
                  axis.text.x=element_text(angle=90, hjust=1))
  #gg = gg + theme(panel.margin.x=unit(1, "lines"))
  return(gg)
}

plot_v_factor = function(coef_data, main_data, ps, n, dir="abs", randoms=1, gene_ids=NULL,
                         slope_name="um.from.reference",
                         y_name="voom1",
                         group_name="brain_side") {
  tp =  coef_data %>% 
    dplyr::filter(reference2_position3==ps, coef==slope_name, coef_type=="Estimate", level=="fixed")
  
  if (is.null(gene_ids)) {
    if (dir=="abs") {
      tp = tp %>%
        top_n(n, abs(value))
    } else if (dir=="up") {
      tp = tp %>% 
        top_n(n, value)
    } else if (dir=="down") {
      tp = tp %>%
        top_n(n, desc(value))
    }
  } else {
    tp = tp %>% filter(gene_id %in% gene_ids)
  }
  
  tp = tp %>% mutate(gene_id = str_replace(gene_id, "\\.", ""))
  print(tp)
  tp2 = coef_data %>% 
    dplyr::filter(reference2_position3==ps, gene_id %in% tp$gene_id)
  
  coef = NULL
  if (randoms==2) {
    coef = tp2 %>% 
      dplyr::filter(coef_type=="Estimate") %>% 
      group_by(gene_id, coef) %>%
      mutate(value_diff = value[level=="fixed"] + value) %>%
      dplyr::filter(level!="fixed") %>%
      dplyr::select(-value) %>%
      spread(coef, value_diff)
    
    colnames(coef)[colnames(coef)=="(Intercept)"] = "intercept"
    colnames(coef)[colnames(coef)==slope_name] = "slope"
  } else if (randoms==1) {
    inters = tp2 %>% 
      dplyr::filter(coef=="(Intercept)", coef_type=="Estimate") %>% 
      group_by(gene_id) %>%
      mutate(value_diff = value[level=="fixed"] + value) %>%
      dplyr::filter(level!="fixed")
    
    #dots = list(slope_name, "vst_combat")
    slopes = tp2 %>% 
      dplyr::filter(level=="fixed", coef==slope_name, coef_type=="Estimate")
    coef = left_join(inters[,c("gene_id", "reference2_position3", "value_diff", "brain_side", "tags")], 
                     slopes[,c("gene_id", "reference2_position3", "value")], by=c("gene_id", "reference2_position3"))
    colnames(coef)[grep("value", colnames(coef))] = c("intercept", "slope")
    coef = coef %>% 
      mutate_(.dots=setNames(list(0,0), c(slope_name, "vst_combat")))
  }
  
  tmp = NULL
  if (grepl("-", ps)) {
    reg = str_split(ps, "-")[[1]][2]
    ref = str_split(ps, "-")[[1]][1]
    tmp = main_data %>% 
      dplyr::filter(region==reg, reference2==ref, gene_id %in% tp$gene_id)
  } else {
    tmp = main_data %>% 
      dplyr::filter(region==ps, gene_id %in% tp$gene_id)
  }
  
  #tmp = tmp %>% mutate(um_from_reference = um_from_reference / 1000)
  
  gg = ggplot(tmp, aes_string(slope_name, y_name, shape=group_name, color="reference2_position3")) + geom_point() + facet_grid(gene_id~tags, scales="free")
  gg = gg + geom_abline(data=coef, aes_string(slope="slope", intercept="intercept"))
  gg = gg + scale_color_manual(values=c("darkblue", "darkred"))
  gg = gg + scale_y_continuous(expand=c(0,1))
  gg = gg + scale_x_continuous(expand=c(0,200), breaks=seq(round(min(tmp[,slope_name]), -3), 
                                                           round(max(tmp[,slope_name]), -3),
                                                           by=1000))
  #gg = gg + geom_text(aes(label=id))
  if (ps=="posterior-ra") { 
    gg = gg + labs(x="um from posterior surface", y="log2 CPM", title=toupper(ps))
  } else if (ps=="dorsal-ra") {      
    gg = gg + labs(x="um from dorsal surface", y="log2 CPM", title=toupper(ps))
  } else if (grepl("lateral", ps)) {
    gg = gg + labs(x="um from lateral surface", y="log2 CPM", title=toupper(ps))
  }
  gg = gg + theme_classic()
  gg = gg + theme(strip.text.y=element_text(angle=0))
  gg = gg + theme(panel.margin.x=unit(1, "lines"))
  return(gg)
}

plot_coefs = function(coef_data, ps, n, dir="abs", pvalue_filter=NULL,
                      slope_name="um.from.reference", gene_ids=NULL) {
  
  coef_data = coef_data %>%
    distinct(coef, coef_type, level, reference2_position3, tags, brain_side, gene_id, .keep_all=T)
  tp = filter_lme_for_plot(coef_data = coef_data, 
                           ps = ps, 
                           n = n, 
                           dir = dir, 
                           gene_ids = gene_ids, 
                           slope_name = slope_name, 
                           pvalue_filter = pvalue_filter)
  # if (is.numeric(pvalue_filter)) {
  #   p_thresh = coef_data %>% 
  #     filter(coef_type == "Pr(>|t|)", coef==slope_name, level=="fixed", reference2_position3==ps, value<pvalue_filter)
  #   coef_data = coef_data %>% filter(gene_id %in% p_thresh$gene_id)
  # }
  # tp = coef_data %>% dplyr::filter(reference2_position3==ps, coef==slope_name, coef_type=="Estimate", level=="fixed")
  # 
  # if (dir=="abs") {
  #   tp = tp %>%
  #     top_n(n, abs(value))
  # } else if (dir=="up") {
  #   tp = tp %>% 
  #     top_n(n, value)
  # } else if (dir=="down") {
  #   tp = tp %>%
  #     top_n(n, desc(value))
  # } else if (dir=="both") {
  #   tp_up = tp %>% top_n(n/2, value)
  #   tp_down = tp %>% top_n(n/2, desc(value))
  #   tp = bind_rows(tp_up, tp_down)
  # }
  #print(head(tp))
  tp2 = coef_data %>%
    dplyr::filter(reference2_position3==ps, gene_id %in% tp$gene_id, coef==slope_name, level=="fixed") 
  
  if (nrow(tp2)==0)
    return(NULL)
  tp2 = tp2 %>% 
    spread(coef_type, value) %>%
    arrange(Estimate) %>%
    mutate(estimate_rank = rank(Estimate),
           gene_id_rank = factor(estimate_rank, labels=gene_id)) %>%
    mutate(value_fdr = ifelse(is.na(value_fdr), 1, value_fdr))
  #print(head(tp2))
  colnames(tp2) = make.names(colnames(tp2))
  
  gg_coef = ggplot() +geom_pointrange(aes(x=gene_id_rank, y=Estimate*1000, ymin=(Estimate*1000-(2*Std..Error*1000)), ymax=(Estimate*1000+2*Std..Error*1000),
                                          color=(value_fdr<.1)))
  gg_coef = gg_coef + geom_hline(yintercept=0, linetype=2)
  gg_coef = gg_coef + scale_color_manual(values=c("black", "red"))
  gg_coef = gg_coef + coord_flip()
  gg_coef = gg_coef + theme_classic()
  gg_coef = gg_coef + labs(x="", y="log2 fold change per mm")
  gg_coef = gg_coef + scale_y_continuous(breaks=seq(-5, 5, 1))
  
  return(gg_coef %+% tp2 + labs(title=toupper(ps)))
}

plot_mixed_combined = function(ps, 
                               models, 
                               model_info, 
                               main_data, 
                               slope_name="um_from_reference", 
                               group_name="brain_side",
                               sum_width=9, 
                               sum_height=9,
                               sum_n = 20,
                               regress_n = 10,
                               save=TRUE, 
                               fig_path,
                               fname_addition="",
                               pvalue_filter=NULL,
                               gene_ids=NULL) {
  
  ggs = lapply(1:length(models), function(i) {
    gg = plot_coefs(models[[i]], ps, sum_n, slope_name=slope_name, dir="both", pvalue_filter=pvalue_filter, gene_ids = gene_ids)
    if (!is.null(gg))
      gg = gg + labs(title=names(models)[i])
    gg
  })
  ggs = compact(ggs)
  fname = paste(ps, fname_addition, "coef_summaries.pdf", sep="-")
  fname = paste(fig_path, fname, sep="/")
  
  cairo_pdf(fname, width = sum_width, height=sum_height)
  do.call(grid.arrange, c(ggs, ncol=2))
  dev.off()
  
  direction = "down"
  for (i in 1:length(models)) {
    fname = paste(paste(c(ps,   fname_addition, c(model_info[i,c("database_table", "value_name", "shortcall")]), direction), collapse="-"), ".pdf", sep="")
    fname = paste(fig_path, fname, sep="/")
    print(fname)
    if (save) 
      cairo_pdf(fname)
    randoms = 1
    if (model_info[i,"shortcall"]=="fixSlope_randInterSlope" )
      randoms = 2
    gg = plot_v_pos_rand(models[[i]], 
                         main_data[[i]], 
                         ps, 
                         regress_n, 
                         direction, 
                         randoms=randoms, 
                         slope_name=slope_name, 
                         group_name = group_name,
                         pvalue_filter = pvalue_filter, gene_ids = gene_ids)
    if (!is.null(gg)) {
      
      gg = gg + labs(title=sprintf("%s : %s", ps, names(models)[i]))
      print(gg)
    }
    if (save)
      dev.off()
  }
  
  direction = "up"
  for (i in 1:length(models)) {
    fname = paste(paste(c(ps,  fname_addition, c(model_info[i,c("database_table", "value_name", "shortcall")]), direction, fname_addition), collapse="-"), ".pdf", sep="")
    fname = paste(fig_path, fname, sep="/")
    print(fname)
    cairo_pdf(fname)
    randoms = 1
    if (model_info[i,"shortcall"]=="fixSlope_randInterSlope")
      randoms = 2
    gg = plot_v_pos_rand(models[[i]],
                         main_data[[i]], 
                         ps, 
                         regress_n, 
                         direction, 
                         randoms=randoms, 
                         slope_name=slope_name, 
                         group_name = group_name,
                         pvalue_filter = pvalue_filter, gene_ids=gene_ids)
    if (!is.null(gg)) {
      gg = gg + labs(title=sprintf("%s : %s", ps, names(models)[i]))
      print(gg)
    }
    dev.off()
  }
}


position_lme_load_data = function(db_name, db_table, info_db=NULL, info_name=NULL, info_inner_join=NULL, value_name, to_voom=TRUE) {
  db_prefix = "~/data/umi_db"
  db = src_sqlite(paste(db_prefix, db_name, sep="/"))
  tbls = src_tbls(db)
  
  stopifnot(db_table %in% tbls)
  full2 = collect(tbl(db, db_table), n=Inf)
  
  info = NULL
  # if (is.null(info_db))
  #   info_db = db
  
  if (is.null(info_db)) {
    info = data.frame(collect(tbl(db, "db_info"), n=Inf))
  } else {
    db_info = src_sqlite(paste(db_prefix, info_db, sep="/"))
    if (is.null(info_name)) {
      info = data.frame(collect(tbl(db_info, "db_info"), n=Inf))
    } else {
      info = data.frame(collect(tbl(db_info, info_name), n=Inf))
    }

  }
  info$group = factor(paste(info$bird, info$position, sep="_"))
  info$position = factor(info$position)
  info$region = factor(info$region)
  info$barcode = as.character(info$barcode)
  
  if (!is.null(info_inner_join))
    info = info %>% inner_join(info_inner_join)
  full2 = inner_join(full2, info, by=c("barcode", "id2", "index1", "id"))
  full2 = full2 %>% distinct(id, gene_id, .keep_all=T)
  dat = acast(full2, gene_id~id, value.var=value_name, fill=0)
  
  if (grepl("combat", value_name)) {
    dat = 2^(dat)
  }
  
  info2 = info
  info2 = info2 %>% filter(id %in% colnames(dat))
  dat = dat[,match(info2$id, colnames(dat))]
  info2$tags_side_position = paste(info2$tags, info2$brain_side, info2$position3, sep="_")
  
  if (to_voom) {
    nf =  estimateSizeFactorsForMatrix(dat)
    design = model.matrix(~1, data=info2)
    v = voom(dat, design, lib.sizes=colSums(dat) * nf, plot=F)
  } else {
    v = list(E=dat, weights=matrix(1, nrow=nrow(dat), ncol=ncol(dat)))
  }
  return(list(dat=dat, info=info2, v=v))
}


position_lme_sub = function(data, dat1, weights, gene_id, right_side, model) {
  data = data
  if (model=="linear") {
    mod = lmer(as.formula(paste(gene_id, right_side)), data=data, weights=weights[rownames(dat1)==gene_id,], control = lmerControl(calc.derivs = FALSE))
  } else if (model=="logistic") {
    gene_id_var = sym(gene_id)
    data = data %>% 
      filter((!! gene_id_var) > 0)
    dat1 = dat1[,match(data$id, colnames(dat1))]
    weights = weights[,match(data$id, colnames(weights))]

    data = data %>% mutate(value = if_else((!!gene_id_var) > median(!!gene_id_var), 1, 0))
    mod = glmer(as.formula(paste("value", right_side)),
                data=data,
                weights=weights[rownames(dat1)==gene_id,],
                family=binomial)
  }
  
  # Process fixed effects
  s = summary(mod)$coefficients
  # sm = melt(s)[1:4,]
  sm = melt(s)
  colnames(sm) = c("coef", "coef_type", "value")
  sm$level = "fixed"
  
  # Process random effects
  r = as.matrix(ranef(mod)$tags_side_position)
  rse = se.ranef(mod)$tags_side_position
  rm = melt(r)
  rm$coef_type = "Estimate"
  rsem = melt(rse)
  rsem$coef_type = "Std. Error"
  l2 = rbind(rm, rsem)
  colnames(l2)[1:2] = c("level", "coef")
  
  # Variance component
  vc1 = as.data.frame(VarCorr(mod))
  vc = data.frame(level="random", coef="(Intercept)", coef_type="variance component", value=vc1[1,4]/sum(vc1[,4]))
  df = bind_rows(sm, l2, vc)
  df$gene_id = gene_id
  df
}

position_lm_sub = function(data, dat1, weights, gene_id, right_side, model) {
  
  if (model=="linear") {
    mod = lm(as.formula(paste(gene_id, right_side)), data=data, weights=weights[rownames(dat1)==gene_id,])
  } else if (model=="logistic") {
    mod = glmer(as.formula(paste(gene_id, right_side)), 
                data=data, 
                weights=weights[rownames(dat1)==gene_id,], 
                family=binomial)
  }

  s = summary(mod)$coefficients
  sm = melt(s)
  colnames(sm) = c("coef", "coef_type", "value")
  sm$level = "fixed"
  
  df = sm
  df$gene_id = gene_id
  df
}


position_lme = function(db_name, 
                        db_table,  
                        info_db=NULL,
                        info_name=NULL, 
                        info_inner_join=NULL,
                        value_name,
                        variance_thresh, 
                        expression_thresh,
                        model="linear",
                        model_type="mixed",
                        permute=FALSE,
                        to_test=NULL,
                        right_side, 
                        to_voom=TRUE, 
                        outdir="~/data2/rstudio/birds/topography/models/mixed/",
                        songsystem_only = F) {
  
  data = position_lme_load_data(db_name, db_table, info_db, info_name, info_inner_join, value_name, to_voom=to_voom)
  info2 = data$info
  print(table(info2$songsystem))

  v = data$v

  if (is.null(to_test)) 
    to_test = unique(info2$position3)
  
    mods = lapply(to_test, function(rp) {
      print(rp)
      
      inf = info2 %>% filter(reference2_position3==rp)
      if (nrow(inf)==0) {
        print("reference2_position3 not found")
        return(NULL)
      }
      
      if (songsystem_only) {
        if (!("songsystem" %in% colnames(inf))) {
          print("Need to specify songsystem column.")
          return(NULL)
        }
        inf = inf %>% mutate(songsystem = if_else(songsystem==1, TRUE, FALSE)) %>% filter(songsystem)
      }
      
      inf = inf %>% distinct(id, .keep_all=T)
      print(table(inf$songsystem))
      dat = v$E
      dat = dat[,match(inf$id, colnames(dat))]
      
      gene_mean = apply(dat, 1, mean, na.rm=T)
      q = quantile(gene_mean, expression_thresh)
      exprs_sel = names(gene_mean)[gene_mean>=q]
      gene_var = apply(dat, 1, var, na.rm=T)
      total = sum(gene_var)
      gene_var = sort(gene_var, decreasing = T)
      gene_var_cum = cumsum(gene_var) / total
      
      gene_var_cum1 = gene_var_cum[gene_var_cum<=variance_thresh]
      sel = names(gene_var_cum1)
      sel = intersect(sel, exprs_sel)
      m = match(inf$id, colnames(v$E))
      
      dat1 = v$E[,m]
      m2 = rownames(dat1) %in% sel
      dat1 = dat1[m2, ]
      w1 = v$weights[,m]
      w1 = w1[m2,]
      dimnames(w1) = list(rownames(dat1), colnames(dat1))
      
      cat(sprintf("Number of genes: %s", nrow(dat1)))
      if (model=="linear") {
        tmp = data.frame(t(dat1), inf, check.names=T)        
      } else if (model=="logistic") {
      #   dat2 = t(apply(dat1, 1, function(x) {
      #     (x - min(x)) / (max(x) - min(x))
      #   }))
        dat2 = dat1
        tmp = data.frame(t(dat2), inf, check.names=T)  
      }

      rownames(dat1) = colnames(tmp)[1:nrow(dat1)]
     #a = lapply(c("PCP4"), function(gene_id) {
     
      if (model_type=="classical") {
        id_list = split(inf, inf$tags_side_position)
      }
      #a = foreach(gene_id = c("PCP4", "NEUROD6")) %dopar%  {
      #a = lapply(colnames(tmp)[6:10], function(gene_id) {
      a = foreach(gene_id = isplitVector(colnames(tmp)[1:nrow(dat1)], chunkSize=1)) %dopar%  {
        if (model_type=="mixed") {
          #tryCatch({
          position_lme_sub(tmp, dat1, w1, gene_id, right_side, model)
          
        } else if (model_type=="classical") {
          res = lapply(id_list, function(i) {
            ids = i$id
            m = match(ids, colnames(dat1))
            position_lm_sub(tmp %>% filter(id %in% ids),
                            dat1[,m],
                            w1[,m],
                            gene_id, right_side, model)
          })
          bind_rows(res, .id="tags_side_position")
        }
      }
    #})
      a = compact(a)
      a = bind_rows(a)
      a$reference2_position3 = rp
      #a$pos_orient = paste(pos, unique(inf$orientation.y), sep="_")
      a
    })
      
  mods1 = bind_rows(mods)
  
  s = str_split(mods1$level, "_")
  s1 = lapply(s, function(x) {
    if (length(x)==1) {
      return(c(x,NA,NA))
    } else {
      return(x)
    }
  })
  #s1 = bind_rows(s1)
  s1 = do.call(rbind, s1)
  s1 = as.data.frame(s1)
  colnames(s1) = c("tags", "brain_side", "position3")
  mods2 = data.frame(mods1, s1)
  
  #save = T
  prefix = paste(outdir, str_replace(basename(db_name), ".db", ""), sep="/")
  dir.create(prefix, recursive=T)
  
  models = list.files(prefix)
  models = models[-grep("info", models)]
  if (length(models)==0) {
    next_model = 1
  } else {
    models = unlist(lapply(str_split(models, "model"), function(x) x[2]))
    next_model = max(as.numeric(str_replace(models, ".rds", ""))) + 1
  }
  model = paste("model", next_model, ".rds", sep="")
  
  out = list(data=v, coef=mods2)
  saveRDS(out, paste(prefix, model, sep="/"))
  
  
  models_file_fname = paste(prefix, "models_info.txt", sep="/")
  new_entry = data.frame(datetime = as.character(now()),
                         database = db_name,
                         database_table = db_table,
                         value_name = value_name,
                         call = right_side,
                         regions = paste(to_test, collapse=","),
                         variance_thresh = variance_thresh,
                         expression_thresh = expression_thresh,
                         fname = model)
  print(new_entry)
  if (file.exists(models_file_fname)) {
    models_file = read.delim(models_file_fname)
    models_file = bind_rows(models_file, new_entry)
  } else {
    models_file = new_entry
  }
  
  write.table(models_file, models_file_fname, quote=F, sep="\t", row.names=F)
  rm(data)
  gc()
  return(mods)
}

# http://stats.stackexchange.com/questions/168650/how-to-set-custom-contrasts-with-lmer-in-r
lme_across_positions = function(dat,
                                weights,
                                info,
                                variance_thresh=NULL, 
                                positions=NULL,
                                right_side,
                                contrasts=NULL,
                                outdir) {
  
  
  if (is.null(positions)) 
    positions = unique(info$position)
  
  system.time({
    coefs = lapply(positions, function(pos) {
      print(pos)
      
      info2 = droplevels(info %>% filter(position==pos))
      print(dim(info2))
      if (nrow(info2)==0) {
        print("position not found")
        return(NULL)
      }
      #inf_sub = inf %>% filter(position==pos)
      #inf2 = info2 %>% filter(position)
      #dat = v$E
      col_ind = match(as.character(info2$id), colnames(dat))
      dat1 = dat[,col_ind]
      weights1 = weights[,col_ind]
      
      # Filter by variance --------------------------------
      if (!is.null(variance_thresh)) {
        gene_var = apply(dat1, 1, var)
        total = sum(gene_var)
        gene_var = sort(gene_var, decreasing = T)
        gene_var_cum = cumsum(gene_var) / total
        
        gene_var_cum1 = gene_var_cum[gene_var_cum<=variance_thresh]
        sel = names(gene_var_cum1)
        row_ind = match(sel, rownames(dat1))
        dat1 = dat1[row_ind, ]
        weights1 = weights1[row_ind,]
      }
      
      print(paste("N genes = ", nrow(dat1), sep=))
      
      tmp = data.frame(t(dat1), info2, check.names=T)
      
      rownames(dat1) = colnames(tmp)[1:nrow(dat1)]
      print(dim(tmp))
      #a = lapply(c("PCP4", "NEUROD6"), function(gene_id) {
      #a = foreach(gene_id = c("PCP4", "NEUROD6")) %do%  {
      
      #prefix = paste(outdir, pos, sep="/")
      #dir.create(prefix, recursive=T)
      
      #a = lapply(colnames(tmp)[1:nrow(dat1)], function(gene_id) {
      
      #a = foreach(gene_id = isplitVector(colnames(tmp)[1:nrow(dat1)], chunkSize=1)) %do%  {
      a = foreach(gene_id = isplitVector(colnames(tmp)[1:nrow(dat1)], chunkSize=1)) %dopar%  {
        if (is.null(contrasts)) {
          mod = lmer(as.formula(paste(gene_id, right_side)), 
                     data=tmp, 
                     weights=weights1[rownames(dat1)==gene_id,], 
                     control = lmerControl(calc.derivs = FALSE))
        } else {
          
          mod = lmer(as.formula(paste(gene_id, right_side)), 
                     data=tmp, 
                     weights=weights1[rownames(dat1)==gene_id,], 
                     control = lmerControl(calc.derivs = FALSE),
                     contrasts = contrasts )
        }
        extract_coefs(mod)
        #fname = paste(prefix, paste(gene_id, ".rds", sep=""), sep="/")
        #saveRDS(mod, fname)
        #return("a")
        #return(mod)
        # Process fixed effects
        # s = summary(mod)$coefficients
        # # sm = melt(s)[1:4,]
        # sm = melt(s)
        # colnames(sm) = c("coef", "coef_type", "value")
        # sm$level = "fixed"
        # 
        # # Process random effects
        # r = as.matrix(ranef(mod)[[1]])
        # rse = se.ranef(mod)[[1]]
        # rm = melt(r)
        # rm$coef_type = "Estimate"
        # rsem = melt(rse)
        # rsem$coef_type = "Std. Error"
        # l2 = rbind(rm, rsem)
        # colnames(l2)[1:2] = c("level", "coef")
        # 
        # # Variance component
        # vc1 = as.data.frame(VarCorr(mod))
        # vc = data.frame(level="random", coef="(Intercept)", coef_type="variance component", value=vc1[1,4]/sum(vc1[,4]))
        # df = bind_rows(sm, l2, vc)
        # df$gene_id = gene_id
        # df
        #})
      }
    #})
      names(a) = rownames(dat1)
      coefs1 = bind_rows(a, .id="gene_id")
      #a$pos_orient = paste(pos, unique(inf$orientation.y), sep="_")
      coefs1
      #names(a) = rownames(dat1)
      # names(a) = c("PCP4","NEUROD6")
      #a
      
    })
    names(coefs) = positions
    coefs_df = bind_rows(coefs, .id="position")
    coefs_df
    # names(mods) = positions
  })
  return(coefs_df)
}

extract_coefs_batch = function(mods, levels=NULL) {
  
  #if(is.null(levels))
  #  stop("Must specify list structure")
  coefs = lapply(mods, function(mod1) {
    
    coefs1 = foreach(mod2 = mod1, chunkSize=1) %dopar%  {
      extract_coefs(mod2)
    }
    #})
    coefs1 = rbind_all(coefs1, .id="gene_id")
    #a$pos_orient = paste(pos, unique(inf$orientation.y), sep="_")
    coefs1
  })
  coefs_df = rbind_all(coefs, .id="position")
  coefs_df
}

extract_coefs = function(mod) {
  
  # Process fixed effects
  s = summary(mod)$coefficients
  # sm = melt(s)[1:4,]
  sm = melt(s)
  colnames(sm) = c("coef", "coef_type", "value")
  sm$level = "fixed"
  
  # Process random effects
  r_list = ranef(mod)
  r = bind_rows(lapply(r_list, function(x) {
    data.frame(x) %>% rownames_to_column(var="level")
  }))
  
  rse_list = se.ranef(mod)
  rse = bind_rows(lapply(rse_list, function(x) {
    data.frame(x) %>% rownames_to_column(var="level")
  }))
  rm = melt(r)
  rm$coef_type = "Estimate"
  rsem = melt(rse)
  rsem$coef_type = "Std. Error"
  l2 = rbind(rm, rsem)
  colnames(l2)[1:2] = c("level", "coef")
  l2$level = as.character(l2$level)
  
  # Variance component
  vc1 = as.data.frame(VarCorr(mod))
  vc = data.frame(level="random", coef="(Intercept)", coef_type=vc1$grp, value=vc1$vcov/sum(vc1$vcov))
  df = bind_rows(sm, l2, vc)
  df
}

process_lme_across_positions = function(db_name, db_table, mods) {
  s = str_split(mods$level, "_")
  s1 = lapply(s, function(x) {
    if (length(x)==1) {
      return(c(x,NA,NA))
    } else {
      return(x)
    }
  })
  s1 = do.call(rbind, s1)
  s1 = as.data.frame(s1)
  #colnames(s1) = c("tags", "brain_side", "position3")
  mods2 = data.frame(mods1, s1)
  
  #save = T
  prefix = paste(outdir, str_replace(db_name, ".db", ""), sep="/")
  dir.create(prefix)
  
  models = list.files(prefix)
  models = models[-grep("info", models)]
  if (length(models)==0) {
    next_model = 1
  } else {
    models = unlist(lapply(str_split(models, "model"), function(x) x[2]))
    next_model = max(as.numeric(str_replace(models, ".rds", ""))) + 1
  }
  model = paste("model", next_model, ".rds", sep="")
  
  out = list(data=v, coef=mods2)
  saveRDS(out, paste(prefix, model, sep="/"))
  
  
  models_file_fname = paste(prefix, "models_info.txt", sep="/")
  new_entry = data.frame(datetime = as.character(now()),
                         database = db_name,
                         database_table = db_table,
                         value_name = value_name,
                         call = right_side,
                         regions = paste(to_test, collapse=","),
                         variance_thresh = variance_thresh,
                         fname = model)
  print(new_entry)
  if (file.exists(models_file_fname)) {
    models_file = read.delim(models_file_fname)
    models_file = rbind(models_file, new_entry)
  } else {
    models_file = new_entry
  }
  
  write.table(models_file, models_file_fname, quote=F, sep="\t", row.names=F)
}

plot_heatmaps = function(dat, grouping_factor, sample_factor, value.var="Estimate", fname=NULL, colors="") {
  
  print(sprintf("Plotting %s...", value.var))
  groups = levels(dat[,grouping_factor])
  if(is.null(groups))
    stop("Grouping_factor must have levels.")
  
  if (!(sample_factor %in% colnames(dat)))
    stop("Sample_factor not in dat")
  
  if (!is.null(fname))
    cairo_pdf(fname, width=12, height=8)
  
  nrow = 2
  ncol = ceiling(length(groups)/ nrow)
  par(mfcol=c(nrow, ncol), mar=c(0,1,0,1))
  
  dots = list(lazyeval::interp(quote(s), s=as.name(sample_factor)))
  names(dots) = "sample"
  dat = dat %>% rename_(.dots=dots)
  
  if (colors=="viridis") {
    colors = viridis_pal()(20)
    b = seq(min(dat[,value.var]), max(dat[,value.var]), length.out=51)
  } else {
    colors = "-RdYlBu2:100"
    b = seq(min(dat[,value.var]), max(dat[,value.var]), length.out=101)
  }
  
  #try({
  lapply(groups, function(p) {
    print(sprintf("...%s",p))
    dots = lazyeval::interp(quote(g==p), .values=list(g=as.name(grouping_factor),
                                            p=p))
    tmp1 = dat %>% filter_(.dots=dots)
    print(nrow(tmp1))
    if (nrow(tmp1) < 2)
      return(NULL)
    tmp = acast(tmp1,
                response~sample, 
                value.var=value.var)
    aheatmap(tmp, main=p, Colv=NA, breaks=b, color = colors)
  })
  #})
  if (!is.null(fname))
    dev.off()
}

lmer_across_list = function(dat_list,
                            info_list,
                            weights_list=NULL,
                            
                            predictors, 
                            rand_predictors,
                            
                            group_name="group", 
                            formatted_vars_df=NULL,
                            
                            heatmap=T, 
                            heatmap_fname=NULL,
                            
                            sample_factor="coef", 
                            value.var="Estimate",
                            
                            base_dir=NULL) {
  
  test = unlist(lapply(info_list, function(x) {
    predictors %in% colnames(x)
  }))
  if(!all(test))
    stop("All predictors not found in info_list.")
  
  if(is.null(names(dat_list)))
    stop("dat_list must be named.")
  
  predictors = c(predictors, rand_predictors)
  preds = paste(predictors, collapse="+")
  des = as.formula(paste("response~",preds, sep=""))
  mods = list()
  coefs = list()
  for(p in 1:length(dat_list)) {
    print(names(dat_list)[p])
    info_cur = info_list[[p]]
    
    #print(dim(info2))
    dat = dat_list[[p]]
    
    weights = NULL
    if (!is.null(weights_list)) {
      weights = weights_list[[p]] 
    } else {
      weights = NULL
    }
    res = lapply((1:ncol(dat)), function(i) {
    #res = mclapply((1:ncol(dat)), function(i) {
    #res = foreach(i=isplitVector(1:ncol(dat), chunkSize=1)) %dopar% {
    #for (i in 1:ncol(dat)) {
      tmp = data.frame(response=dat[,i], info_cur[match(rownames(dat), info_cur$id),], check.names=T)
      #tmp$tags = factor(tmp$tags)
      print(i)
      
      if (!is.null(weights_list)) {
        print("here")
        print(dim(weights))
        print(i)
        print(head(weights[,i]))
        ws1 = weights[,1]
        print(eval(i))
        print(class(ws1))
        d = list(i=i)
        mod = lmer(des, data=tmp, weights=base::eval(i, sys.frame()), control = lmerControl(calc.derivs = FALSE))
      } else {
        mod = lmerTest::lmer(des, data=tmp, control = lmerControl(calc.derivs = FALSE))        
      }
      coef = extract_coefs(mod)
      coef$response = colnames(dat)[i]
      coef[,group_name] = names(dat_list)[p]
      #mods = c(mods, list(mod))
      #coefs = c(coefs, list(coef))
      #return(list(mod=mod, coef=coef))
      return(list(coef=coef))
    #}
    })
    # }, mc.cores=6)
    #mod1 = lapply(res, function(x) x$mod)
    coef1 = bind_rows(lapply(res, function(x) x$coef))
    coefs = c(coefs, list(coef1))
  }
  
  coefs1 = bind_rows(coefs)
  dots = list(lazyeval::interp("factor(gr, levels=na)",.values=list(gr=as.name(group_name), 
                                                          na=names(dat_list))))
  names(dots) = group_name
  coefs1 = coefs1 %>% mutate_(.dots=dots)
  coefs3 = coefs1 %>% filter(coef != "(Intercept)") %>% spread(coef_type, value) %>% mutate(ci95.hi = Estimate + 2*`Std. Error`,
                                                                                            ci95.lo = Estimate - 2*`Std. Error`,
                                                                                            includes_zero = (ci95.hi * ci95.lo) <0)
  
  if (!is.null(formatted_vars_df)) {
    coefs3 = coefs3 %>% left_join(pretty_coef)    
  }
  
  if (heatmap) {
    print(levels(coefs3[,group_name]))
    plot_heatmaps(coefs3, grouping_factor=group_name, sample_factor=sample_factor, value.var=value.var, fname=heatmap_fname)
    
    if (!is.null(heatmap_fname))
      dev.off()
  }
  
  # Write out
  if (!is.null(base_dir)) {
    dir = paste(base_dir, "models", sep="/")
    dir.create(dir, recursive=T, showWarnings = T)
    fname = paste(paste(predictors, collapse="-"), ".rds", sep="")
    saveRDS(mods, paste(dir, fname, sep="/"))
    
    dir = paste(base_dir, "coefs", sep="/")
    dir.create(dir, recursive=T, showWarnings = F)
    fname = paste(paste(predictors, collapse="-"), ".rds", sep="")
    saveRDS(coefs3, paste(dir, fname, sep="/"))
  }
  
  return(coefs3)
}

flip_across_list = function(dat_list, 
                            info_list,
                            predictors, 
                            covariates,
                            rand_predictors,
                            
                            group_name="group", 
                            formatted_vars_df=NULL,
                            
                            heatmap=T, 
                            heatmap_fname=NULL,
                            
                            sample_factor="coef", 
                            value.var="Estimate",
                            
                            base_dir=NULL) {
  
  
  test = unlist(lapply(info_list, function(x) {
    predictors %in% colnames(x)
  }))
  if(!all(test))
    stop("All predictors not found in info_list.")
  
  if(is.null(names(dat_list)))
    stop("dat_list must be named.")
  
  #predictors = c(predictors, rand_predictors)
  preds = paste(predictors, collapse="+")
  des = as.formula(paste("~",preds, sep=""))
  
  covar = paste(covariates, collapse="+")
  des_covar = as.formula(paste("~", covar, sep=""))
  
  units = as.formula(paste("~", rand_predictors, sep=""))
  
  mods = list()
  coefs = list()

  for(p in 1:length(dat_list)) {
    print(names(dat_list)[p])
    info_cur = info_list[[p]]
    #print(dim(info2))
    dat = dat_list[[p]]
    try({
    for (i in 1:ncol(dat)) {
      tmp = na.omit(data.frame(response=dat[,i], info_cur[match(rownames(dat), info_cur$id),], check.names=T))
      tmp = data.frame(response=dat[,i], info_cur[match(rownames(dat), info_cur$id),], check.names=T)
      
      if ((i %% 100) == 0)
        cat("=")
      #fm = flipMix(modelWithin=response~1, X=des, Z=des_covar, units=units, data=tmp, testType = "permutation")
      jk = lmm.jack(response~um_from_reference + (um_from_reference|tags_side_position), data=tmp, method="MINQUE")
      res = fm@res
      coef_name = unlist(lapply(str_split(rownames(res), "\\|_"), function(x) x[2]))
      coef = data.frame(coef=coef_name, tvalue=res[,2], pvalue=res[,4])
      coef$response = colnames(dat)[i]
      coef[,group_name] = names(dat_list)[p]
      #mods = c(mods, list(mod))
      coefs = c(coefs, list(coef))
    }
      cat("\n")
    })
  }

  coefs1 = bind_rows(coefs)
  dots = list(lazyeval::interp("factor(gr, levels=na)",.values=list(gr=as.name(group_name), 
                                                          na=names(dat_list))))
  names(dots) = group_name
  coefs1 = coefs1 %>% mutate_(.dots=dots)
  
  
  if (!is.null(formatted_vars_df)) {
    coefs1 = coefs1 %>% left_join(pretty_coef)    
  }
  
  # Write out
  if (!is.null(base_dir)) {
    dir = paste(base_dir, "permutation", sep="/")
    dir.create(dir, recursive=T, showWarnings = T)
    fname = paste(paste(predictors, collapse="-"), ".rds", sep="")
    saveRDS(coefs1, paste(dir, fname, sep="/"))
  }
  
  return(coefs1)
}

lmer_across_positions_with_plots = function(dat_list, 
                                            info_list, 
                                            weights_list = NULL,
                                            
                                            predictors, 
                                            rand_predictors,
                                            
                                            group_name, 
                                            formatted_vars_df,
                                            
                                            base_dir) {
  a = lmer_across_list(dat_list,
                       info_list,
                       weights_list,
                       
                       predictors=predictors, 
                       rand_predictors = rand_predictors, 
                       
                       group_name = group_name, 
                       formatted_vars_df = formatted_vars_df, 
                       
                       heatmap=F,
                       
                       base_dir=base_dir)
  
  coef_var = "coef_nice"
  if (is.null(formatted_vars_df))
    coef_var = "coef"
  try({
  hfname_dir = paste(base_dir, "heatmaps", sep="/")
  dir.create(hfname_dir, recursive = T, showWarnings = F)
  hfname = paste(hfname_dir, paste(predictors, collapse="-"), sep="/")
  
  hfname1 = paste(hfname, "estimate.pdf", sep="_")
  plot_heatmaps(a, "position", coef_var, "Estimate", fname=hfname1)
  
  hfname1 = paste(hfname, "tvalue.pdf", sep="_")
  plot_heatmaps(a, "position", coef_var, "t value", fname=hfname1)
  
  hfname1 = paste(hfname, "pvalue.pdf", sep="_")
  a = a %>%
    mutate(`Pr(>|t|)` = ifelse(`Pr(>|t|)`==0, .Machine$double.neg.eps, `Pr(>|t|)`)) %>%
    mutate(pvalue_neg_log = -1 * log10(`Pr(>|t|)`)) #%>%
    #mutate(pvalue_neg_log = ifelse(is.infinite(pvalue_neg_log),, pvalue_neg_log))
  plot_heatmaps(a, "position", coef_var, "pvalue_neg_log", fname=hfname1, colors="viridis")
  
  hfname1 = paste(hfname, "fdr.pdf", sep="_")
  a = a %>% mutate(fdr = p.adjust(`Pr(>|t|)`, method="BH"),
                   fdr_neg_log = -1 * log10(fdr))
  plot_heatmaps(a, "position", coef_var, "fdr_neg_log", fname=hfname1, colors="viridis")
  
  hfname1 = paste(hfname, "qvalue.pdf", sep="_")
  a = a %>% mutate(q = qvalue(`Pr(>|t|)`)$qvalues,
                   q_neg_log= -1 * log10(q))
  plot_heatmaps(a, "position", coef_var, "q_neg_log", fname=hfname1, colors="viridis")
  
  
  hfname1 = paste(hfname, "pvalue_cut.pdf", sep="_")
  thresh = -1 * log10(.1)
  a = a %>% mutate(pvalue_neg_log_cut = ifelse(pvalue_neg_log>thresh, pvalue_neg_log, 0))
  plot_heatmaps(a, "position", coef_var, "pvalue_neg_log_cut", fname=hfname1, colors="viridis")
  })
  return(a)
}

flip_across_positions_with_plots = function(dat_list, 
                                            info_list, 
                                            
                                            predictors, 
                                            covariates,
                                            rand_predictors,
                                            
                                            group_name="group", 
                                            formatted_vars_df=NULL,
                                            
                                            heatmap=T, 
                                            heatmap_fname=NULL,
                                            
                                            sample_factor="coef", 
                                            value.var="Estimate",
                                            
                                            base_dir=NULL) {
  a = flip_across_list(dat_list,
                       info_list,
                       
                       predictors=predictors, 
                       covariates=covariates,
                       rand_predictors = rand_predictors, 
                       
                       group_name = group_name, 
                       formatted_vars_df = formatted_vars_df, 
                       
                       heatmap=F,
                       
                       base_dir=base_dir)
  
  hfname_dir = paste(base_dir, "permutation/heatmaps", sep="/")
  dir.create(hfname_dir, recursive = T, showWarnings = F)
  hfname = paste(hfname_dir, paste(c(predictors, covariates), collapse="-"), sep="/")
  
  hfname1 = paste(hfname, "tvalue.pdf", sep="_")
  plot_heatmaps(a, "position", "coef_nice", "tvalue", fname=hfname1)
  
  hfname1 = paste(hfname, "pvalue.pdf", sep="_")
  a = a %>% mutate(pvalue_neg_log = -1 * log10(pvalue))
  plot_heatmaps(a, "position", "coef_nice", "pvalue_neg_log", fname=hfname1, colors="viridis")
  
  hfname1 = paste(hfname, "pvalue_cut.pdf", sep="_")
  thresh = -1 * log10(.05)
  a = a %>% mutate(pvalue_neg_log_cut = ifelse(pvalue_neg_log>thresh, pvalue_neg_log, 0))
  plot_heatmaps(a, "position", "coef_nice", "pvalue_neg_log_cut", fname=hfname1, colors="viridis")
  
  return(a)
}
