library(wrapr)
library(gridExtra)
library(igraph)

# WGCNA -------------------------------------------------------------------
plot_scale_free = function(sft) {
  par(mfrow = c(1,2), mar = c(5, 2, 2, 2))
  cex1 = 0.9
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,sAigned R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  
  abline(h=0.90,col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}


# General -----------------------------------------------------------------
plot_sub_network = function(ig, gene_id, neighborhood_size=1, info, color_column, label_column, layout=NULL, ...) {
  info = info[match(V(ig)$name, info$gene_id),] # Ensure ordered
  
  ig_ind = neighborhood(ig, neighborhood_size, nodes=V(ig)[grep(gene_id, V(ig)$name)])
  info1 = info[ig_ind[[1]],]
  ig1 = induced_subgraph(ig, ig_ind[[1]])
  info1 = info1[match(V(ig1)$name, info1$gene_id),]
  if (is.null(layout))
    layout = layout_nicely(ig1)
  
  par(mar=c(0,0,0,0))
  #pos = "hvc"
  plot(ig1, 
       #vertex.size=as.numeric(info1$node.degree)/50,
       
       #vertex.size=info1$node.degree/10,
       #vertex.color=ifelse(abs(gene_info2[,paste("cor",pos,sep=".")])>.25, 
       #                    gene_info2[,paste("cor_color",pos,sep=".")], 
       #                     "#00000000"),
       vertex.color=info1[,color_column],
       #vertex.label=NA,
       vertex.label=info1[,label_column], 
       vertex.label.color=ifelse(info1$gene_id == gene_id, "red", "black"), 
       vertex.label.family="sans",
       #vertex.label.cex=sqrt(info1$node.degree/50),
       vertex.frame.color=NA, 
       edge.arrow.size=0, 
       edge.width=E(ig1)$weight,
       margin=c(0,0,0,0),
       layout=layout,
       ...)
  return(layout)
}



plot_sub_network_by_module = function(ig, 
                                      module, 
                                      info, 
                                      color_column, 
                                      label_column, 
                                      layout=NULL, 
                                      save=F,
                                      ...) {
  info = info[match(V(ig)$name, info$gene_id),] # Ensure ordered
  mod_ids = info$gene_id[info[,module]=="YES"]
  ig_ind = which(V(ig)$name %in% mod_ids)
  #ig_ind = neighborhood(ig, neighborhood_size, nodes=V(ig)[V(ig)$name %in% mod_ids])
  info1 = info[ig_ind,]
  ig1 = induced_subgraph(ig, ig_ind)
  info1 = info1[match(V(ig1)$name, info1$gene_id),]
  if (is.null(layout))
    layout = layout_nicely(ig1)
  
  if (save) {
    dir.create("plots", showWarnings = F)
    cairo_pdf(file=paste("plots", paste(module, ".pdf", sep=""), sep="/"), width=12, height=8)
  }
  
  par(mar=c(0,0,0,0))
  #pos = "hvc"
  plot(ig1, 
       #vertex.size=as.numeric(gene_info1$node.degree)/50,
       #vertex.size=info1$node.degree/10,
       #vertex.color=ifelse(abs(gene_info2[,paste("cor",pos,sep=".")])>.25, 
       #                    gene_info2[,paste("cor_color",pos,sep=".")], 
       #                     "#00000000"),
       vertex.color=info1[,color_column],
       #vertex.label=NA,
       vertex.label=info1[,label_column], 
       #vertex.label.color=ifelse(info1$gene_id == module, "red", "black"), 
       vertex.label.family="sans",
       #vertex.label.cex=sqrt(info1$node.degree/50),
       vertex.frame.color=NA, 
       edge.arrow.size=0, 
       edge.width=E(ig1)$weight,
       margin=c(0,0,0,0),
       layout=layout,
       ...)
  if (save)
    dev.off()
  return(layout)
}

plot_sub_network_by_select = function(ig, 
                                      to_plot,
                                      color_column, 
                                      label_column, 
                                      layout=NULL, 
                                      save=F,
                                      ...) {
  #info = info[match(V(ig)$name, info$gene_id),] # Ensure ordered
  #mod_ids = info$gene_id[info[,module]=="YES"]
  ig_ind = which(V(ig)$name %in% to_plot)
  #ig_ind = neighborhood(ig, neighborhood_size, nodes=V(ig)[V(ig)$name %in% mod_ids])
  #info1 = info[ig_ind,]
  ig1 = induced_subgraph(ig, ig_ind)
  #info1 = info1[match(V(ig1)$name, info1$gene_id),]
  if (is.null(layout))
    layout = layout_nicely(ig1)
  
  if (save) {
    dir.create("plots", showWarnings = F)
    cairo_pdf(file=paste("plots", paste(module, ".pdf", sep=""), sep="/"), width=12, height=8)
  }
  
  par(mar=c(0,0,0,0))
  #pos = "hvc"
  plot(ig1, 
       #vertex.size=as.numeric(gene_info1$node.degree)/50,
       #vertex.size=info1$node.degree/10,
       #vertex.color=ifelse(abs(gene_info2[,paste("cor",pos,sep=".")])>.25, 
       #                    gene_info2[,paste("cor_color",pos,sep=".")], 
       #                     "#00000000"),
       #vertex.color=info1[,color_column],
       #vertex.label=NA,
       #vertex.label=info1[,label_column], 
       #vertex.label.color=ifelse(info1$gene_id == module, "red", "black"), 
       vertex.label.family="sans",
       #vertex.label.cex=sqrt(info1$node.degree/50),
       vertex.frame.color=NA, 
       edge.arrow.size=0, 
       edge.width=E(ig1)$weight,
       margin=c(0,0,0,0),
       layout=layout,
       ...)
  if (save)
    dev.off()
  return(layout)
}

plot_sub_network_cor = function(ig,
                                gene_correlations,
                                gene, 
                                neighborhood_size=1, 
                                info, 
                                color_column, 
                                label_column, 
                                layout=NULL, 
                                save=F,
                                ...) {

  if (save) {
    dir.create("plots", showWarnings = F)
    cairo_pdf(file=paste("plots", paste(gene, ".pdf", sep=""), sep="/"), width=12, height=8)
  }
  par(mfrow=c(1,2))
  l = plot_sub_network(ig=ig, 
                       gene_id=gene, 
                       neighborhood_size=neighborhood_size, 
                       info=info, 
                       color_column=color_column, 
                       label_column=label_column, 
                       ...)
  
  plot.new()
  vps = baseViewports()
  pushViewport(vps$figure)
  vp1 = plotViewport(c(0,0,0,1))
  gg = ggplot(gene_correlations %>% filter(gene_id==gene), aes(position_nice, cor)) + geom_point()
  gg = gg + theme_bw()
  gg = gg + geom_hline(yintercept=0, linetype=2)
  gg = gg + coord_flip()
  gg = gg + labs(x="", y="PCC")
  print(gg, vp=vp1)
  if (save)
    dev.off()
}

plot_sub_network_by_module_cor = function(ig,
                                gene_correlations,
                                module, 
                                info, 
                                color_column, 
                                label_column, 
                                layout=NULL, 
                                save=F,
                                ...) {
  
  if (save) {
    dir.create("plots", showWarnings = F)
    cairo_pdf(file=paste("plots", paste(module, ".pdf", sep=""), sep="/"), width=12, height=8)
  }
  par(mfrow=c(1,2))
  l = plot_sub_network_by_module(ig=ig, 
                       module=module, 
                       info=info, 
                       color_column=color_column, 
                       label_column=label_column, 
                       ...)
  
  plot.new()
  vps = baseViewports()
  pushViewport(vps$figure)
  vp1 = plotViewport(c(0,0,0,1))
  mod_ids = gene_info1$gene_id[gene_info1[,module]=="YES"]
  gg = ggplot(gene_pos_cors_m %>% 
                filter(gene_id %in% mod_ids) %>% 
                group_by(position_nice) %>% 
                summarize(cor_mean=mean(cor), cor_sd=sd(cor)), aes(position_nice, cor_mean)) + geom_pointrange(aes(ymin=(cor_mean-cor_sd),ymax=(cor_mean+cor_sd)))
  gg = gg + theme_bw()
  gg = gg + geom_hline(yintercept=0, linetype=2)
  gg = gg + coord_flip()
  gg = gg + labs(x="", y="PCC")
  print(gg, vp=vp1)
  if (save)
    dev.off()
}



plot_gene_by_position = function(data, genes, value) {
  gg = ggplot(data %>% filter(gene_id %in% genes), 
              aes_string("position3", value))
  gg = gg + geom_point(position=position_jitter(width=.2))
  gg = gg +  facet_wrap(~gene_id)
  gg = gg + theme_classic()
  print(gg)
}



plot_gmm_radar = function(data) {
  require(ggradar)
  label = data[1,1]
  gg = ggradar(data, 
               grid.label.size=4,
               axis.label.size=4, 
               group.line.width=.5,
               group.point.size=3, 
               legend.title=label,
               plot.legend=F,
               grid.min = 0,
               grid.mid = .5,
               grid.max= 1)
  gg = gg + theme(legend.position="none",
                  text = element_text(size = 20, family = "arial"))
  gg = gg + labs(title=label)
  gg
}

plot_reference_gmms_radar = function(data, filter_col, filter_val, genes, ref_genes, fname=NULL) {
  #let(list(fc = filter_col, fv = filter_val),
  #    {
        tmp = data %>% 
          mutate(gmms_mod1 = map(gmms_mod1, function(d) d[rownames(d) %in% genes,])) %>%
          select(gmms_mod1)
  #    })
  tmp = data.frame(tmp[[1]][[1]])
  tmp = tmp %>% rownames_to_column(var="gene_id")
  
  tmp1 = tibble(data=split(tmp, tmp$gene_id))
  tmp1 = tmp1 %>% mutate(ggs=map(data, plot_gmm_radar))
  
  if (!is.null(fname)) 
  cairo_pdf(fname, height=10, width=10)

  do.call(grid.arrange, tmp1$ggs)
  
  if (!is.null(fname))
  dev.off()
}
