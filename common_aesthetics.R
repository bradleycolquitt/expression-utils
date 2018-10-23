#position_colors = data.frame(position=c("ra", "arco", "hvc", "ncl", "lman", "nido", "x", "stri", "meso", "ov"),
#                             color=c("#7f0000", "#ef6548",  "#662506", "#cc4c02", "#08306b", "#4292c6", "#49006a", "#ae017e", "#004529", "#41ab5d"))
#position_colors = c("#e31a1c", "#ef6548",  "#662506", "#cc4c02", "#08306b", "#4292c6", "#49006a", "#ae017e", "darkblue", "#004529", "#41ab5d")
position_levels = c("ra", "arco", "hvc", "ncl", "lman", "nido", "x", "stri", "fieldl", "ov", "meso")
position_pretty_levels = c("RA", "Arco", "HVC", "NCL", "LMAN", "Nido", "Area X", "Striatum", "Field L", "Ovoid", "Meso")
#names(position_colors) = c("ra", "arco", "hvc", "ncl", "lman", "nido", "x", "stri", "lfs", "meso", "ov")

position_colors = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
names(position_colors) = c("nido", "lman", "stri", "x", "arco", "ra", "ncl", "hvc", "fieldl", "ov", "dummy", "meso")

position_pretty_colors = position_colors
names(position_pretty_colors) = c("Nido", "LMAN", "Striatum","Area X", "Arco", "RA", "NCL", "HVC", "Field L", "Ovoid", "dummy", "Meso")

position_table = data.frame(position=names(position_colors), position_pretty=names(position_pretty_colors), colors=position_pretty_colors)
position_table$position_pretty = factor(position_table$position_pretty, levels=position_pretty_levels)
position_table$position = factor(position_table$position, levels=position_levels)
deaf_colors = c("darkred", "darkblue")
names(deaf_colors) = c("deaf", "intact")

#duration.of.experiment_colors = c("#1f78b4", "#33a02c", '#e31a1c')
duration.of.experiment_colors = c("grey80", "grey50", "grey10")
names(duration.of.experiment_colors) = c("4", "9", "14")

lesion_group_colors = c("#1f78b4", "#33a02c", "#e31a1c", "#6a3d9a")
names(lesion_group_colors) = c("intact-FALSE", "intact-TRUE", "deaf-FALSE", "deaf-TRUE")

lesion_group_colors2 = c("#1f78b4", "#a6cee3", "#e31a1c", "#fb9a99" )
names(lesion_group_colors2) = c("hearing-contra", "hearing-ipsi", "deaf-contra", "deaf-ipsi")