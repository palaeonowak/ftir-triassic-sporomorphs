# load required packages
library(corrplot)
library(viridisLite)
library(prospectr)
library(dendextend)

# load spectra and associated metadata
spectra<-read.csv("spectra.csv", fileEncoding = "UTF-8", header=TRUE)
meta<-read.csv("meta.csv", fileEncoding = "UTF-8", header=TRUE)

nwn<-ncol(spectra)
wavenumbers<-as.numeric(unlist(strsplit(names(spectra),"wn_"))[(1:nwn)*2])

# short names for plot labels
simplespecs<-unique(meta[,c(6,11)])
simplespecs$Taxon_shortform<-unname(sapply(simplespecs$Taxon, \(x) switch(x,
                                                                          "Asterotheca cf. thalensis"="As. cf. thalensis (E)",
                                                                          "Asterotheca merianii"="As. merianii (L)",
                                                                          "Asterotheca sp. 1"="As. sp. 1 (L)",
                                                                          "Asterotheca sp. 2"="As. sp. 2 (L)",
                                                                          "Merianopteris augusta"="M. augusta (K)",
                                                                          "Anomopteris mougeotii"="An. mougeotii (D)",
                                                                          "Gordonopteris lorigae"="G. lorigae (D)",
                                                                          "Scolopendrites grauvogelii"="S. grauvogelii (D)",
                                                                          "Scolopendrites scolopendrioides"="S. scolopendrioides (D)",
                                                                          "Scolopendrites sp."="S. sp. (D)",
                                                                          "Todites linnaeifolius"="T. linnaeifolius (L)",
                                                                          "Todites sp."="T. sp. (L)",
                                                                          "Dictyophyllum serratum"="D. serratum (L)",
                                                                          "Voltzia recubariensis"="V. recubariensis (D)",
                                                                          "Willsiostrobus acuminatus"="W. acuminatus (G)",
                                                                          "Willsiostrobus cf. willsii"="W. cf. willsii (E)")))



### define functions for analysis

# calculate mean spectrum from selected rows
avg_group_spectra<-function(main_column, row_select=TRUE){
  col_content<-meta[, names(meta)==main_column]
  groups<-unique(col_content[row_select])
  wn_mean<-Vectorize(function(x,y){mean(spectra[row_select & col_content==groups[y], x])})
  group_spectra<-outer(1:nwn, 1:length(groups), FUN = wn_mean)
  colnames(group_spectra)<-groups
  return(group_spectra)
}

# mean correlation coefficient per row in a correlation matrix, without diagonal identity cells (always 1)
cormean<-function(cormat){apply(cormat, 1, FUN=function(x) (sum(x)-1)/(length(x)-1))}

# variance per row in a correlation matrix, without diagonal identity cells
corvar<-function(cormat){sapply(1:nrow(cormat), FUN=function(x) var(cormat[x,-x]))}

# statistics for selected spectra
avg_comp<-function(sel_rows, sel_names, main_group, title){
  sgw<-11 # window size for Savitzky-Golay smoothing
  wn_mean<-Vectorize(function(x,y){mean(spectra[sel_rows[[y]],x])})
  mean_specs<-outer(1:nwn, 1:length(sel_rows), FUN = wn_mean)
  colnames(mean_specs)<-sel_names
  rownames(mean_specs)<-colnames(spectra)
  smoothed<-apply(mean_specs, 2, function(x) savitzkyGolay(x, m=0, p=2, w=sgw))
  smoothed_2nd_diffs<-apply(mean_specs, 2, function(x) savitzkyGolay(x, m=2, p=2, w=sgw))
  results<-list(spectra=as.data.frame(t(mean_specs)),
                smoothed_spectra=as.data.frame(t(smoothed)),
                smoothed_2nd_diffs=as.data.frame(t(smoothed_2nd_diffs)),
                selected_rows=sel_rows,
                Kendall_tau=cor(smoothed, method = "kendall"),
                Kendall_tau_2nd=cor(smoothed_2nd_diffs, method = "kendall"),
                Kendall_p=cor.mtest(smoothed, method="kendall")$p)
  results$key_values<-data.frame(batch=title,
                                 main_group=main_group,
                                 sub_group=colnames(mean_specs),
                                 count=sapply(sel_rows, length),
                                 min=apply(mean_specs, 2, min),
                                 max=apply(mean_specs, 2, max),
                                 variance=apply(mean_specs, 2, var),
                                 Kendall_mean=cormean(results$Kendall_tau),
                                 Kendall_2nd_mean=cormean(results$Kendall_tau_2nd),
                                 Kendall_var=corvar(results$Kendall_tau),
                                 Kendall_2nd_var=corvar(results$Kendall_tau_2nd))
  return(results)
}

# statistics for groups of selected spectra
avg_group_comp<-function(title, main_column, group_column, row_select=TRUE){
  main_col_content<-meta[, names(meta)==main_column]
  group_col_content<-meta[, names(meta)==group_column]
  sel_names<-unique(main_col_content[row_select])
  sel_rows<-lapply(sel_names, FUN = function(x) which(main_col_content==x & row_select))
  groups<-sapply(sel_rows, FUN = function(x) unique(group_col_content[x]))
  
  results<-lapply(unique(groups), function(x) avg_comp(sel_rows[which(groups==x)], sel_names[which(groups==x)], x, title))
  return(results)
}

# normalization (for plotting spectra)
normalize<-function(x){
  x<-x-mean(x)
  x<-x/max(x)
  return(x)
}

# draw axis ticks within plot area
internal_axis_ticks <- function(draw.axes = 1:4, number.axes = 1:2, ...){
  for(n in draw.axes){
    axis(n, tck = 0.01, labels = NA, ...)
  }
  for(n in number.axes){
    axis(n, lwd = 0, line = -0.5, las = 1, ...)
  }
  box(lwd = 1.5)
}

# plot rescaled average spectra and hierarchical cluster dendrogram
# from groups of selected spectra
avgplotspectrum_x_combo<-function(row_select_list){
  ns<-length(row_select_list)
  wn_mean<-Vectorize(function(x,y){mean(spectra[row_select_list[[y]],x])})
  mean_specs<-outer(1:nwn, 1:length(row_select_list), FUN = wn_mean)
  mean_specs<-apply(mean_specs, 2, normalize)
  colnames(mean_specs)<-names(row_select_list)
  rownames(mean_specs)<-colnames(spectra)
  Kendall_dist=as.dist(1-cor(mean_specs, method = "kendall"))
  cols<-turbo(ns, begin = 0.1, end = 0.9)
  linetypes<-c("solid", "44", "13", "1343", "73", "2262", "131343")
  dendro<-as.dendrogram(hclust(Kendall_dist, method = "average"))
  dendro<-set(set(set(dendro, "leaves_pch", 19), "leaves_col", cols[order(names(row_select_list))][rank(labels(dendro))]), "leaves_cex", 1.2)
  
  par(mfrow = c(1, 2), mgp=c(2.2,1,0))
  par(mar=c(3.5,3.5,1,1))
  plot(wavenumbers, mean_specs[,1], type="l", ylab="rescaled relative absorbance", xlab = expression(paste("wavenumber [cm"^"-1","]")), xlim = c(max(wavenumbers), min(wavenumbers)), ylim = c(min(specs), max(specs)), axes = FALSE, col=cols[1], lty=linetypes[1], lwd=2)
  for(il in 2:ns){
    lines(wavenumbers, mean_specs[,il], type="l", col=cols[il], lty=linetypes[il], lwd=2)
  }
  text(x=1800, y=1, labels=substitute(paste(bold('A'))), adj = c(0,0.5))
  legend("bottomright",
         legend = paste(names(row_select_list), " [n=", sapply(row_select_list, length), "]", sep = ""),
         col = cols,
         lty = linetypes,
         lwd = 2)
  internal_axis_ticks(1:2, 1:2)
  
  par(mar=c(4,0,0,13))
  plot(dendro, horiz = TRUE)
  text(x=max(get_branches_heights(dendro)), y=count_terminal_nodes(dendro), labels=substitute(paste(bold('B'))), adj = c(0,0))
}

# save plot as PDF
pdf_comparison<-function(row_select_list, title){
  pdf(file = paste0("Comparison_", title, ".pdf"), 13, 6)
  avgplotspectrum_x_combo(row_select_list)
  dev.off()
}


### calculate results

## grouped selections

# group genera by class and compare
results_gen_per_class<-avg_group_comp("Genera per class", "Genus", "Class", meta$Treatment=="Macerated")

# select taxa from genera with at least two taxa
results_spp_per_gen<-avg_group_comp("Species per genus", "Taxon", "Genus", meta$Treatment=="Macerated" & meta$Genus %in% c("Asterotheca", "Scolopendrites", "Todites", "Willsiostrobus"))

# select taxa from formations with at least two taxa
results_spp_per_fm<-avg_group_comp("Species per formation", "Taxon", "Formation", meta$Treatment=="Macerated" & meta$Formation %in% c("Erfurt Fm.", "Lunz Fm.", "Dont Fm."))

# only Asterotheca merianii and Gordonopteris lorigae have multiple specimens in the dataset:
#table(unique(meta[,c(6,4)])$Taxon)
# note that some specimens are not represented with macerated data and one only has one scan
results_specim_per_sp<-avg_group_comp("Specimens per species", "specimen", "Taxon", meta$Treatment=="Macerated" & meta$Taxon %in% c("Asterotheca merianii", "Gordonopteris lorigae") & meta$specimen!="1882/0013/3059")

# select specimens with more than one macerated sample
sample_table<-table(unique(meta[meta$Treatment=="Macerated", c(4,5)])$specimen)
results_samp_per_specim<-avg_group_comp("Samples per specimen", "specimen_sample", "specimen", meta$Treatment=="Macerated" & meta$specimen %in% names(sample_table)[sample_table>1])

# select macerated samples with more than one scan
meta$id2<-1:nrow(meta)
scan_table<-table(meta$specimen_sample[meta$Treatment=="Macerated"])
results_scan_per_samp<-avg_group_comp("Scan per sample", "id2", "specimen_sample", meta$Treatment=="Macerated" & meta$specimen_sample %in% names(scan_table)[scan_table>1])


## non-grouped results

#classes
sel_names<-unique(meta$Class[meta$Treatment=="Macerated"])
sel_rows<-lapply(sel_names, FUN = function(x) which(meta$Class==x & meta$Treatment=="Macerated"))
results_classes<-avg_comp(sel_rows, sel_names, "Class", "All classes")

#genera
sel_names<-unique(meta$Genus[meta$Treatment=="Macerated"])
sel_rows<-lapply(sel_names, FUN = function(x) which(meta$Genus==x & meta$Treatment=="Macerated"))
results_genera<-avg_comp(sel_rows, sel_names, "Genus", "All genera")

#species
sel_names<-unique(meta$Taxon[meta$Treatment=="Macerated"])
sel_rows<-lapply(sel_names, FUN = function(x) which(meta$Taxon==x & meta$Treatment=="Macerated"))
results_species<-avg_comp(sel_rows, sel_names, "Taxon", "All taxa")

#specimens
sel_names<-unique(meta$specimen[meta$Treatment=="Macerated"])
sel_rows<-lapply(sel_names, FUN = function(x) which(meta$specimen==x & meta$Treatment=="Macerated"))
results_specimens<-avg_comp(sel_rows, sel_names, "Specimen", "All specimens")

#samples
sel_names<-unique(meta$specimen_sample[meta$Treatment=="Macerated"])
sel_rows<-lapply(sel_names, FUN = function(x) which(meta$specimen_sample==x & meta$Treatment=="Macerated"))
results_samples<-avg_comp(sel_rows, sel_names, "Sample", "All samples")

#scans
sel_names<-unique(meta$id2[meta$Treatment=="Macerated"])
sel_rows<-lapply(sel_names, FUN = function(x) which(meta$id2==x & meta$Treatment=="Macerated"))
results_scans<-avg_comp(sel_rows, sel_names, "Scan", "All scans")


### save results

resultlist<-mget(apropos("results_"))

saveRDS(resultlist, "results.R")


keyvals<-list()
for(i in 1:length(resultlist)){
  if(is.null(names(resultlist[[i]]))){
    for(ii in 1:length(resultlist[[i]])){
      keyvals<-append(keyvals, list(resultlist[[i]][[ii]]$key_values))
    }
  }else{
    keyvals<-append(keyvals, list(resultlist[[i]]$key_values))
  }
}
keyvals<-do.call(rbind, keyvals)

write.csv(keyvals, "results_key_values.csv", row.names = FALSE, fileEncoding = "UTF-8")


### plotting

## correlation plots

pdf(file = "correlation_kendall_classes.pdf", 3, 3)
corrplot(results_classes$Kendall_tau, method = "color", type = "full", col.lim = c(0,1), tl.pos = "lt", addCoef.col="white", tl.col="black", diag = FALSE)
dev.off()

pdf(file = "correlation_kendall_genera.pdf", 10, 10)
corrplot(results_genera$Kendall_tau, method = "color", type = "full", col.lim = c(0,1), tl.pos = "lt", addCoef.col="white", tl.col="black", diag = FALSE) |>
  corrRect(name = c("Asterotheca", "Anomopteris", "Voltzia", "Willsiostrobus"), col = "grey80", lwd = 4)
dev.off()

cormat_ks<-results_species$Kendall_tau
colnames(cormat_ks)<-simplespecs$Taxon_shortform
rownames(cormat_ks)<-simplespecs$Taxon_shortform

pdf(file = "correlation_kendall_species.pdf", 10, 10)
corrplot(cormat_ks, method = "color", type = "full", col.lim = c(0,1), tl.pos = "lt", addCoef.col="white", tl.col="black", diag = FALSE) |>
  corrRect(name = c("As. cf. thalensis (E)", "An. mougeotii (D)", "V. recubariensis (D)", "W. cf. willsii (E)"), col = "grey80", lwd = 4)
dev.off()


## plot rescaled mean spectra and hierarchical cluster dendrograms

sel_rows<-list("Marattiopsida"=which(meta$Class=="Marattiopsida" & meta$Treatment=="Macerated"),
               "Polypodiopsida"=which(meta$Class=="Polypodiopsida" & meta$Treatment=="Macerated"),
               "Pinopsida"=which(meta$Class=="Pinopsida" & meta$Treatment=="Macerated"))
pdf_comparison(sel_rows, "Classes")

sel_rows<-list("Asterotheca"=which(meta$Genus=="Asterotheca" & meta$Treatment=="Macerated"),
               "Merianopteris"=which(meta$Genus=="Merianopteris" & meta$Treatment=="Macerated"),
               "Anomopteris"=which(meta$Genus=="Anomopteris" & meta$Treatment=="Macerated"),
               "Dictyophyllum"=which(meta$Genus=="Dictyophyllum" & meta$Treatment=="Macerated"),
               "Gordonopteris"=which(meta$Genus=="Gordonopteris" & meta$Treatment=="Macerated"),
               "Scolopendrites"=which(meta$Genus=="Scolopendrites" & meta$Treatment=="Macerated"),
               "Todites"=which(meta$Genus=="Todites" & meta$Treatment=="Macerated"))
pdf_comparison(sel_rows, "Fern genera")

sel_rows<-list("V. recubariensis (D)"=which(meta$Taxon=="Voltzia recubariensis" & meta$Treatment=="Macerated"),
               "W. acuminatus (G)"=which(meta$Taxon=="Willsiostrobus acuminatus" & meta$Treatment=="Macerated"),
               "W. cf. willsii (E)"=which(meta$Taxon=="Willsiostrobus cf. willsii" & meta$Treatment=="Macerated"),
               "A. merianii (L)"=which(meta$Taxon=="Asterotheca merianii" & meta$Treatment=="Macerated"))
pdf_comparison(sel_rows, "Conifers+Asterotheca")

sel_rows<-list("A. merianii macerated"=which(meta$Taxon=="Asterotheca merianii" & meta$Treatment=="Macerated"),
               "D. serratum macerated"=which(meta$Taxon=="Dictyophyllum serratum" & meta$Treatment=="Macerated"),
               "Todites sp. macerated"=which(meta$Taxon=="Todites sp."  & meta$Treatment=="Macerated"),
               "A. merianii non-macerated"=which(meta$Taxon=="Asterotheca merianii" & meta$Treatment=="Unmacerated"),
               "D. serratum non-macerated"=which(meta$Taxon=="Dictyophyllum serratum"  & meta$Treatment=="Unmacerated"),
               "Todites sp. non-macerated"=which(meta$Taxon=="Todites sp."  & meta$Treatment=="Unmacerated"))
pdf_comparison(sel_rows, "Maceration")

sel_rows<-list("A. cf. thalensis (E)"=which(meta$Taxon=="Asterotheca cf. thalensis" & meta$Treatment=="Macerated"),
               "W. cf. willsii (E)"=which(meta$Taxon=="Willsiostrobus cf. willsii" & meta$Treatment=="Macerated"),
               "Dont ferns (D)"=which(meta$Class=="Polypodiopsida" & meta$Treatment=="Macerated"),
               "V. recubariensis (D)"=which(meta$Taxon=="Voltzia recubariensis" & meta$Treatment=="Macerated"))
pdf_comparison(sel_rows, "Preservation")

sel_rows<-list("1885/0012/3939"=which(meta$specimen=="1885/0012/3939" & meta$Treatment=="Macerated"),
               "1885/0012/4069"=which(meta$specimen=="1885/0012/4069" & meta$Treatment=="Macerated"),
               "2019/0185/0003"=which(meta$specimen=="2019/0185/0003" & meta$Treatment=="Macerated"),
               "MB.Pb.2003/1099"=which(meta$specimen=="MB.Pb.2003/1099" & meta$Treatment=="Macerated"))
pdf_comparison(sel_rows, "Asterotheca specimens")
