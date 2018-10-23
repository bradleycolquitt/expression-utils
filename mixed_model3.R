options(stringsAsFactors=F)
library(dplyr)
library(tidyr)
filter_by_var = function(dat, thresh) {
  gene_var = apply(dat, 1, var)
  total = sum(gene_var)
  gene_var = sort(gene_var, decreasing = T)
  gene_var_cum = cumsum(gene_var) / total
  
  gene_var_cum1 = gene_var_cum[gene_var_cum<=thresh]
  return(names(gene_var_cum1))
}

plot_v_pos_rand = function(coef_data, main_data, ps, n, dir="abs", randoms=1, gene_ids=NULL,
                           slope_name="um.from.reference",
                           y_name="voom1",
                           group_name="brain.side") {
  tp =  coef_data %>% 
    dplyr::filter(position3==ps, coef==slope_name, coef_type=="Estimate", level=="fixed")
  
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
  print(tp)
  tp2 = coef_data %>% 
    dplyr::filter(position3==ps, gene_id %in% tp$gene_id)
  
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
    coef = left_join(inters[,c("gene_id", "position3", "value_diff", "brain.side", "tags")], 
                     slopes[,c("gene_id", "position3", "value")], by=c("gene_id", "position3"))
    colnames(coef)[grep("value", colnames(coef))] = c("intercept", "slope")
    coef = coef %>% 
      mutate_(.dots=setNames(list(0,0), c(slope_name, "vst_combat")))
  }
  
  # colnames(coef)[grep("value", colnames(coef))] = c("intercept", "slope")
  # coef = coef %>% 
  #   mutate(um.from.reference=0, vst_combat=0)
  
  tmp = main_data %>% 
    dplyr::filter(region==ps, gene_id %in% tp$gene_id)
  
  gg = ggplot(tmp, aes_string(slope_name, y_name, color=group_name))
  gg = gg + geom_point()
  gg = gg + geom_point(aes(alpha=factor(songsystem)), fill="white", shape=21) + scale_alpha_manual(values=c(1, 0))
 
  gg = gg + facet_grid(gene_id~tags, scales="free")
  
  #gg = gg + geom_abline(data=coef, aes_string(slope="slope", intercept="intercept", color=group_name))
  gg = gg + geom_abline(data=coef, aes_string(slope="slope", intercept="intercept"))
  gg = gg + scale_color_manual(values=c("darkblue", "darkred"))
  #`gg = gg + scale_fill_manual(values=c("white", ""))
  gg = gg + scale_y_continuous(expand=c(0,1))
  gg = gg + scale_x_continuous(expand=c(0,200), breaks=seq(round(min(tmp[,slope_name]), -3), 
                                                           round(max(tmp[,slope_name]), -3),
                                                           by=1000))
  #gg = gg + geom_text(aes(label=id))
  if (ps=="ra" ) { 
    if (grepl("coronal", tmp$orientation.y[1])) {
      gg = gg + labs(x="um from posterior surface", y="log2 CPM", title=toupper(ps))
    } else if (grepl("horizontal", tmp$orientation.y[1])) {
      gg = gg + labs(x="um from dorsal surface", y="log2 CPM", title=toupper(ps))
    }
  } else {      
    gg = gg + labs(x="um from lateral surface", y="log2 CPM", title=toupper(ps))
  }
  gg = gg + theme_classic()
  gg = gg + theme(strip.text.y=element_text(angle=0))
  return(gg)
}

plot_coefs = function(coef_data, ps, n, dir="abs",
                      slope_name="um.from.reference") {
  tp = coef_data %>% filter(position3==ps, coef==slope_name, coef_type=="Estimate", level=="fixed")
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
  #print(head(tp))
  tp2 = coef_data %>% filter(position3==ps, gene_id %in% tp$gene_id, coef==slope_name, level=="fixed") %>% 
    spread(coef_type, value) %>%
    arrange(Estimate) %>%
    mutate(estimate_rank = rank(Estimate),
           gene_id_rank = factor(estimate_rank, labels=gene_id))
  #print(head(tp2))
  colnames(tp2) = make.names(colnames(tp2))
  
  gg_coef = ggplot() + geom_pointrange(aes(x=gene_id_rank, y=Estimate*1000, ymin=(Estimate*1000-(2*Std..Error*1000)), ymax=(Estimate*1000+2*Std..Error*1000)))
  gg_coef = gg_coef + geom_hline(yintercept=0, linetype=2)
  gg_coef = gg_coef + coord_flip()
  gg_coef = gg_coef + theme_classic()
  gg_coef = gg_coef + labs(x="", y="log2 fold change per mm")
  gg_coef = gg_coef + scale_y_continuous(breaks=seq(-5, 5, 1))
  
  return(gg_coef %+% tp2 + labs(title=toupper(ps)))
}


position_lme_load_data = function(db_name, db_table, value_name) {
  db_prefix = "/media/data/umi_db"
  #db_name = "ra_dv_lonStrDom1_3utr.db" 
  #to_load = "filtered_data_700000"
  
  db = src_sqlite(paste(db_prefix, db_name, sep="/"))
  tbls = src_tbls(db)
  
  stopifnot(db_table %in% tbls)
  full2 = collect(tbl(db, db_table))
  
  info = data.frame(collect(tbl(db, "info")))
  info$group = factor(paste(info$bird, info$position, sep="_"))
  info$position = factor(info$position)
  info$region = factor(info$region)
  
  full2 = inner_join(full2, info, by=c("bc", "i7"))
  
  #value_name = "unique_umi"
  dat = acast(full2, gene_id~id, value.var=value_name)
  dat[is.na(dat)] = 0
  
  if (grepl("combat", value_name)) {
    dat = 2^(dat)
  }
  
  info2 = info
  info2 = info2 %>% filter(id %in% colnames(dat))
  dat = dat[,match(info2$id, colnames(dat))]
  info2$tags_side_position = paste(info2$tags, info2$brain.side, info2$position3, sep="_")
  
  nf =  estimateSizeFactorsForMatrix(dat)
  design = model.matrix(~1, data=info2)
  v = voom(dat, design, lib.sizes=colSums(dat) * nf, plot=F)
  return(list(dat=dat, info=info2, v=v))
}



position_lme = function(db_name, db_table, value_name, variance_thresh, to_test=NULL,
                        right_side) {
  
  data = position_lme_load_data(db_name, db_table, value_name)
  #dat = data$dat
  info2 = data$info
  v = data$v
  
  if (is.null(to_test)) 
    to_test = unique(info2$position3)
  
  system.time({
    mods = lapply(to_test, function(pos) {
      print(pos)
      
      inf = info2 %>% filter(position==pos)
      #inf_sub = inf %>% filter(position==pos)
      #inf2 = info2 %>% filter(position)
      dat = v$E
      gene_var = apply(dat, 1, var)
      total = sum(gene_var)
      gene_var = sort(gene_var, decreasing = T)
      gene_var_cum = cumsum(gene_var) / total
      
      gene_var_cum1 = gene_var_cum[gene_var_cum<=variance_thresh]
      sel = names(gene_var_cum1)
      
      m = match(inf$id, colnames(v$E))
      
      dat1 = v$E[,m]
      m2 = rownames(dat1) %in% sel
      dat1 = dat1[m2, ]
      w1 = v$weights[,m]
      w1 = w1[m2,]
      tmp = data.frame(t(dat1), inf, check.names=T)
      rownames(dat1) = colnames(tmp)[1:nrow(dat1)]
      #a = foreach(gene_id = c("PCP4", "NEUROD6")) %dopar%  {
      a = foreach(gene_id = isplitVector(colnames(tmp)[1:nrow(dat1)], chunkSize=1)) %dopar%  {
        mod = lmer(as.formula(paste(gene_id, right_side)), data=tmp, weights=w1[rownames(dat1)==gene_id,], control = lmerControl(calc.derivs = FALSE))
        
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
      a = rbind_all(a)
      a$position3 = pos
      a$pos_orient = paste(pos, unique(inf$orientation.y), sep="_")
      a
    })
  })
  mods1 = rbind_all(mods)
  
  s = str_split(mods1$level, "_")
  s1 = lapply(s, function(x) {
    if (length(x)==1) {
      return(c(x,NA,NA))
    } else {
      return(x)
    }
  })
  s1 = do.call(rbind, s1)
  s1 = as.data.frame(s1)
  colnames(s1) = c("tags", "brain.side", "position3")
  mods2 = data.frame(mods1, s1)
  
  #save = T
  prefix = "/media/data2/rstudio/birds/topography/models/mixed/"
  
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
  
  return(mods)
}

