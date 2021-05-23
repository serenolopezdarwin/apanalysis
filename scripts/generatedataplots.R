library(ggplot2)
library(gplots) 
library(LSD)
library(pheatmap)
library(plotly)
library(plyr)
# Loaded after plyr for better function masking.
library(dplyr)
library(processx)
library(RColorBrewer)
library(reshape2)
library(viridis)


## Reads in and labels our data.
args <- commandArgs(trailingOnly=TRUE)
dataset <- args[1]
working.path <- args[2]
setwd(working.path)
data.in.path <- paste0(dataset, "/cell_data.tsv")
data <- read.table(data.in.path)
colnames(data) <- c("cluster", "tsne1", "tsne2", "subcluster", "stsne1", "stsne2", 
	                "trajectory", "umap1", "umap2", "umap3", "straj", "rumap1", "rumap2", "rumap3", 
	                "age", "utr", "clusterutr", "trajectoryutr", "strajutr", "ageutr",
	                "cellutr1", "cellutr2", "geneutr3", "geneutr4", "ratio", "pseudotime", "count")
data <- data[!is.na(data$utr),]
gene.stats.in.path <- paste0(dataset, "/gene_stats.txt")
binned.umap.in.path <- paste0(dataset, "/binned_umap_data.txt")
gene.expr.in.path <- paste0(dataset, "/gene_expr_data.txt")


## Generates a plot of our read counts at each filtration level
# We get the unfiltered read count with cat cellbed/*.bed | wc -l 
gene.stats <- read.table(gene.stats.in.path)
read.counts <- log10(c(2596423721, sum(gene.stats[,1]), sum(gene.stats[,2]), sum(gene.stats[,3])))
df <- data.frame(filter=c("All Reads", "Gene-Overlapping Reads", "3'UTR-Overlapping Reads", 
						  "PAS-Overlapping Reads"), reads=read.counts)
df$filter <- factor(df$filter, levels=unique(df$filter))
a <- ggplot(data=df, aes(x=filter, y=reads, group=1)) +
    geom_line() +
    geom_point() +
    theme_bw() + 
    ggtitle("Read Counts At Each Round of Filtration") +
    xlab("Filtration") + ylab("Read Count (log10)") +
	theme(axis.line=element_blank(), axis.text.x=element_text(angle=60, hjust=1),
          panel.border=element_blank(), panel.grid.major=element_blank(),
		  panel.grid.minor=element_blank(), legend.key.size=unit(1.5, "cm"))
gene.stats.out.path <- paste0(dataset, "/figures/gene_stats_", dataset, ".pdf")
ggsave(filename=gene.stats.out.path, plot=a, width=5, height=5, limitsize=FALSE)

## Generates a plot of cell read counts vs cell utr deviation for each cell.
pdf(paste0(dataset, "/figures/expr_vs_utr.pdf"))
heatscatter(data$count, data$utr, main="", xlab="", ylab="", ylim=c(-1000,1000))
title(main="Cell UMIs vs Mean UTR Deviation", xlab="UMIs", ylab="Mean UTR Deviation")
dev.off()
# Evaluates to -0.04366728
print("UMI vs UTR Spearman:")
print(cor(data$count, data$utr, method='spearman'))
# Evaluates to -0.04362961
print("UMI vs UTR Pearson:")
print(cor(data$count, data$utr, method='pearson'))

## Makes plots comparing gene expression data at various rounds of filtration and with the baseline.
# Reads in gene expression and count data.
exp.data <- read.table(gene.expr.in.path,T, row.names=1)
# Log-transforms our expression data.
exp.data$baseline <- log10(exp.data$baseline+0.01)
exp.data$original <- log10(exp.data$original+1)
exp.data$pas <- log10(exp.data$pas+1)
# Draws a heatmap scatter of unfiltered gene counts vs baseline expression.
unfiltered.path <- paste0(dataset, "/figures/expr_baseline_vs_unfiltered.pdf")
pdf(unfiltered.path)
par(mar=c(7,7,4,1)+.1)
heatscatter(exp.data$original, exp.data$baseline, bty='n', cex.axis=2, cex.lab=2, las=1, 
	        xlim=c(1,7), ylim=c(-1.5,4), xlab="", ylab="", main="")
# We need a custom title setting to avoid typesetting issues.
title(main='Unfiltered', 
	  ylab=expression(atop("Median mRNA Level", "(log"[10]*" FPKM)")),
	  cex.main=1.5, cex.lab=1.5)
# Separate title expression for the x label to move it downwards.
title(xlab=expression(atop("Gene Body Counts", "sci-RNA-seq3 "~(log[10]))), cex.lab=1.5, line=5)
dev.off()
# Draws a heatmap scatter of filtered gene counts vs baseline expression.
filtered.path <- paste0(dataset, "/figures/expr_baseline_vs_filtered.pdf")
pdf(filtered.path)
par(mar=c(7,7,4,1)+.1)
heatscatter(exp.data$pas, exp.data$baseline, bty='n', cex.axis=2, cex.lab=2, las=1, 
	        xlim=c(1,7), ylim=c(-1.5,4), xlab="", ylab="", main="")
title(main='Filtered', 
	  ylab=expression(atop("Median mRNA Level", "(log"[10]*" FPKM)")),
	  cex.main=1.5, cex.lab=1.5)
title(xlab=expression(atop("Gene Body Counts", "sci-RNA-seq3 "~(log[10]))), cex.lab=1.5, line=5)
dev.off()
# Draws a few line plots of spearman correlations between filtered and unfiltered data and the baseline.
comparison.path <- paste0(dataset, "/figures/expr_filtered_vs_unfiltered.pdf")
pdf(comparison.path)
exp.data <- exp.data[order(exp.data$baseline, decreasing=T),]
idxs <- seq(200, nrow(exp.data), 100)
plot(idxs, round(unlist(
	 lapply(idxs, function(x){
                             	cor(exp.data$pas[1:x], exp.data$baseline[1:x], method='spearman')
                             } )), 3), 
	 xlab="Top Ranked Genes(by FPKM)", ylab='Spearman Correlation', 
	 col='red', type='l', lwd=2, ylim=c(0,0.8), bty='n')
lines(idxs, round(unlist(
	  lapply(idxs, function(x){
                              	cor(exp.data$original[1:x], exp.data$baseline[1:x], method='spearman')
                              } )), 3), 
	  col='blue', lwd=2, ylim=c(0,0.8))
dev.off()
print("Unfiltered Spearman:")
print(cor(exp.data$baseline, exp.data$original, method="spearman"))
print("Filtered Spearman:")
print(cor(exp.data$baseline, exp.data$pas, method="spearman"))
print("Unfiltered Pearson:")
print(cor(exp.data$baseline, exp.data$original, method="pearson"))
print("Filtered Pearson:")
print(cor(exp.data$baseline, exp.data$pas, method="pearson"))



## Plots dreme data
# We get these values through running dreme analysis, which our head script does for us.
# dreme <- data.frame(Dataset=c("PAS", "No PAS"), Proportion=c(5201/10000, 3478/10000))
# a <- ggplot(dreme, aes(x=Dataset, y=Proportion)) + 
#    geom_bar(stat="identity") +
#    coord_flip() + 
#    theme_bw() + 
#	 theme(panel.border=element_blank(), panel.grid.major=element_blank(),
#		   panel.grid.minor=element_blank())
# dreme.out.path.1 <- paste0(dataset, "/figures/dreme2_AAUAAA.pdf") 
# ggsave(filename=dreme.out.path.1, plot=a, width=5, height=2)
# dreme2 <- data.frame(Dataset=c("PAS", "No PAS"), Proportion=c(1586/10000, 1877/10000))
# a <- ggplot(dreme2, aes(x=Dataset, y=Proportion)) + 
#     geom_bar(stat="identity") +
#     coord_flip() + 
#     theme_bw() + 
#	 theme(panel.border=element_blank(), panel.grid.major=element_blank(),
#		   panel.grid.minor=element_blank()) +
#	 scale_y_continuous(limits=c(0,0.5))
# dreme.out.path.2 <- paste0(dataset, "/figures/dreme2_AAAAAA.pdf")
# ggsave(filename=dreme.out.path.2, plot=a, width=5, height=2)


## Generates a dataframe of trajectory means by ages for our data.
split.data <- split(data, list(data$age, data$trajectory), drop=TRUE)
split.stats <- ldply(split.data, function(sublist){
	subframe <- data.frame(sublist)
	frame.trajectory <- subframe[1,7]
	frame.age <- subframe[1,15]
	frame.median <- median(subframe[,16])
	data.frame(frame.trajectory, frame.age, frame.median)
	})
split.stats <- split.stats[,c(2:4)]
# Recasts our built dataframe to heatmap dimensions
fixed.stats <- acast(split.stats, frame.trajectory~frame.age, value.var="frame.median")
fixed.stats[is.na(fixed.stats)] <- NaN
rownames(fixed.stats) <- gsub("_", " ", rownames(fixed.stats))
# Builds a palette based on pre-defined breaks (these breaks may be adjusted as needed)
# Note that this pallete and breaks are used for later heatmap plots as well.
viridis.palette <- colorRampPalette(viridis(152))
col.breaks <- c(seq(-50, -15,length=51), seq(-14.9,14.9,length=50), seq(15,50,length=50))
## Params for UTR/Age heatmap. Use margins=c(9, 50) to generate properly cropped row labels.
traj.out.path <- paste0(dataset, "/figures/trajectory_heatmap_", dataset, ".pdf")
pdf(traj.out.path, width=25, height=25)
par(oma=c(1,1,1,15))
heatmap.2(fixed.stats, main="Deviation From Mean Utr Length", notecol="black", Colv="NA", cex.main=3,
cexRow=3, margins=c(9,30), density.info="none", trace="none", col=viridis.palette, lhei=c(1,10), 
breaks=col.breaks, dendrogram="row", cexCol=3)
dev.off()
pvals <- apply(fixed.stats, 2, function(col){(col - mean(col)) / sd(col)})
traj.out.path <- paste0(dataset, "/figures/trajectory_heatmap_pvalue_", dataset, ".pdf")
pval.breaks <- c(seq(-2.0, -0.76,length=51), seq(-0.75,0.75,length=50), seq(0.76,2.0,length=50))
pdf(traj.out.path, width=25, height=25)
par(oma=c(1,1,1,15))
heatmap.2(pvals, main="Deviation From Mean Utr Length", notecol="black", Colv="NA", cex.main=3,
cexRow=3, margins=c(9,30), density.info="none", trace="none", col=viridis.palette, lhei=c(1,10), 
breaks=pval.breaks, dendrogram="row", cexCol=3)
dev.off()


## Generates a heatmap of subtrajectories by ages, colored on the rows by parent trajectory
data1 <- data[!is.na(data$trajectory),]
split.data <- split(data1, list(data1$age, data1$straj), drop=TRUE)
n <- length(unique(data1$trajectory))
rm(data1)
split.stats <- ldply(split.data, function(sublist){
	subframe <- data.frame(sublist)
	frame.straj <- subframe[1,11]
	frame.age <- subframe[1,15]
	frame.counts <- length(subframe[,16])
	if (frame.counts < 10) {
	    frame.median <- NA
	} else {
		frame.median <- median(subframe[,16])
	}
	data.frame(frame.straj, frame.age, frame.median)
	})
split.stats <- split.stats[,c(2:4)]
fixed.stats <- acast(split.stats, frame.straj~frame.age, value.var="frame.median")
fixed.stats[is.na(fixed.stats)] <- NaN
# Writes a vector with indices of each heatmap row corresponding to which trajectory it belongs to.
cluster.vector <- as.numeric(as.factor(gsub("\\..*","",rownames(fixed.stats))))
# Writes a vector of colors to correspond to each heatmap row based on trajectory
qual.col.pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col.palette <- unlist(mapply(brewer.pal, qual.col.pals$maxcolors, rownames(qual.col.pals)))
cluster.colors <- col.palette[cluster.vector]
straj.out.path <- paste0(dataset, "/figures/subtrajectory_heatmap_", dataset, ".pdf")
pdf(straj.out.path, width=15, height=45)
par(mar=c(7,4,4,2)+0.1, oma=c(10,1,1,1))
heatmap.2(fixed.stats, colRow='white', main="Deviation From Mean UTR Length", notecol="black", Colv="NA", 
		 lhei=c(1,35), density.info="none", trace="none", col=viridis.palette, margins=c(9,20), cexCol=10, 
		 labRow=FALSE, breaks=col.breaks, dendrogram="none", Rowv=FALSE, RowSideColors=cluster.colors, cexRow=2)
legend('left', legend=unique(gsub("\\..*","",rownames(fixed.stats))), fill=unique(cluster.colors))
dev.off()
# cexCol=10, oma=c(10,1,1,1)


## Generates a heatmap of clusters by ages for our data, same as above but using clusters.
split.data <- split(data, list(data$age, data$cluster), drop=TRUE)
split.stats <- ldply(split.data, function(sublist){
	subframe <- data.frame(sublist)
	frame.cluster <- subframe[1,1]
	frame.age <- subframe[1,15]
	frame.counts <- length(subframe[,16])
	if (frame.counts < 20) {
	    frame.median <- NA
	} else {
		frame.median <- median(subframe[,16])
	}
	data.frame(frame.cluster, frame.age, frame.median)
	})
split.stats <- split.stats[,c(2:4)]
fixed.stats <- acast(split.stats, frame.cluster~frame.age, value.var="frame.median")
fixed.stats[is.na(fixed.stats)] <- NaN
rownames(fixed.stats) <- gsub("_", " ", rownames(fixed.stats))
# Params for UTR/Age heatmap. Use margins=c(9, 50) to generate properly cropped row labels.
cluster.out.path <- paste0(dataset, "/figures/cluster_heatmap_", dataset, ".pdf")
pdf(cluster.out.path, width=25, height=25)
par(oma=c(1,1,1,10))
heatmap.2(fixed.stats, main="Deviation From Mean Utr Length", notecol="black", Colv="NA", cex.main=3,
cexRow=3, margins=c(9,30), density.info="none", trace="none", col=viridis.palette, lhei=c(1,10),
breaks=col.breaks, dendrogram="row", cexCol=3)
dev.off()
pvals <- apply(fixed.stats, 2, function(col){(col - mean(col)) / sd(col)})
cluster.out.path <- paste0(dataset, "/figures/cluster_heatmap_pvalue_", dataset, ".pdf")
pval.breaks <- c(seq(-2.0, -0.76,length=51), seq(-0.75,0.75,length=50), seq(0.76,2.0,length=50))
pdf(cluster.out.path, width=25, height=25)
par(oma=c(1,1,1,15))
heatmap.2(pvals, main="Deviation From Mean Utr Length", notecol="black", Colv="NA", cex.main=3,
cexRow=3, margins=c(9,30), density.info="none", trace="none", col=viridis.palette, lhei=c(1,10), 
breaks=pval.breaks, dendrogram="row", cexCol=3)
dev.off()

## Generates a heatmap of subclusters by ages, colored on the rows by parent clusters.
data1 <- data[!is.na(data$cluster),]
split.data <- split(data1, list(data1$age, data1$subcluster), drop=TRUE)
n <- length(unique(data1$cluster))
rm(data1)
split.stats <- ldply(split.data, function(sublist){
	subframe <- data.frame(sublist)
	frame.subcluster <- subframe[1,4]
	frame.age <- subframe[1,15]
	frame.counts <- length(subframe[,16])
	if (frame.counts < 10) {
	    frame.median <- NA
	} else {
		frame.median <- median(subframe[,16])
	}
	data.frame(frame.subcluster, frame.age, frame.median)
	})
split.stats <- split.stats[,c(2:4)]
fixed.stats <- acast(split.stats, frame.subcluster~frame.age, value.var="frame.median")
fixed.stats[is.na(fixed.stats)] <- NaN
# Writes a vector with indices of each heatmap row corresponding to which cluster it belongs to.
cluster.vector <- as.numeric(as.factor(gsub("\\..*","",rownames(fixed.stats))))
# Writes a vector of colors to correspond to each heatmap row based on cluster.
qual.col.pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col.palette <- unlist(mapply(brewer.pal, qual.col.pals$maxcolors, rownames(qual.col.pals)))
cluster.colors <- col.palette[cluster.vector]
scluster.out.path <- paste0(dataset, "/figures/subcluster_heatmap_", dataset, ".pdf")
pdf(scluster.out.path, width=15, height=45)
par(mar=c(7,4,4,2)+0.1, oma=c(10,1,1,1))
heatmap.2(fixed.stats, colRow='white', main="Deviation From Mean UTR Length", notecol="black", Colv="NA", 
	     lhei=c(1,35), density.info="none", trace="none", col=viridis.palette, margins=c(9,20), cexCol=10, 
	     labRow=FALSE, breaks=col.breaks, dendrogram="none", Rowv=FALSE, RowSideColors=cluster.colors, cexRow=2)
legend('left', legend=unique(gsub("\\..*","",rownames(fixed.stats))), fill=unique(cluster.colors))
dev.off()


## Calculates and plots correlation between heatmap cell counts and values.
neuronal.clusters <- c('Cholinergic_neurons', 'Inhibitory_interneurons', 'Granule_neurons',
'Sensory_neurons', 'Excitatory_neurons', 'Inhibitory_neurons', 'Inhibitory_neuron_progenitors',
'Postmitotic_premature_neurons')
blood.clusters <- c('Primitive_erythroid_lineage', 'Definitive_erythroid_lineage', 'Neutrophils', 'White_blood_cells')
split.stats <- ldply(split.data, function(sublist){
	subframe <- data.frame(sublist)
	frame.median <- median(subframe[,16])
	frame.counts <- length(subframe[,16])
	frame.cluster <- subframe[1,1]
	frame.age <- subframe[1,15]
	if (frame.cluster %in% neuronal.clusters) {
		frame.color <- "red"
	}
	else if (frame.cluster %in% blood.clusters) {
		frame.color <- "blue"
	}
	else {
		frame.color <- "black"
	}
	data.frame(frame.median, frame.counts, frame.color, frame.age)
	})
split.medians <- split.stats$frame.median
split.counts <- log10(split.stats$frame.counts)
split.colors <- as.character(split.stats$frame.color)
split.ages <- split.stats$frame.age
# Evaluates to 0.030
print("Heatmap count pearson:")
print(cor(split.counts, split.medians, method="pearson"))
# Evaluates to .00048
print("Heatmap count spearman:")
print(cor(split.counts, split.medians, method="spearman"))
pdf(paste0(dataset, "/figures/cluster_count_vs_utr.pdf"))
plot(split.counts, split.medians, main="", xlab="", ylab="", pch=19, col=split.colors)
title(main="Cluster & Age Cell Count vs UTR Deviation", xlab="Cell Count", ylab="Deviation From Mean UTR")
dev.off()
# P-value evaluates to 7.9*10^-15
split.aov <- aov(frame.median ~ frame.age, data=split.stats)


## Performs analysis of variance test on each heatmap row for cells within age categories.
cluster.split.data <- split(data, data$cluster, drop=TRUE)
c.split.stats <- ldply(cluster.split.data, function(sublist){
	subframe <- data.frame(sublist)
	frame.counts <- length(subframe[,16])
	frame.cluster <- subframe[1,1]
	frame.aov <- aov(utr ~ age, data=subframe)
	frame.pvalue <- summary(frame.aov)[[1]][["Pr(>F)"]][1]
	data.frame(frame.cluster, frame.pvalue, frame.counts)
	})
adjusted.values = p.adjust(c.split.stats$frame.pvalue, method="bonferroni")
# Gives bonferroni-corrected values for each cluster.
c.split.stats$adjusted.pvalue <-adjusted.values
write.table(c.split.stats, file=paste0(dataset, "/cluster_heatmap_anova.tsv"), quote=FALSE, sep='\t', col.names=NA)

## Does the same as above with trajectory rows rather than cluster rows.
traj.split.data <- split(data, data$trajectory, drop=TRUE)
t.split.stats <- ldply(traj.split.data, function(sublist){
	subframe <- data.frame(sublist)
	frame.counts <- length(subframe[,16])
	frame.cluster <- subframe[1,1]
	frame.aov <- aov(utr ~ age, data=subframe)
	frame.pvalue <- summary(frame.aov)[[1]][["Pr(>F)"]][1]
	data.frame(frame.cluster, frame.pvalue, frame.counts)
	})
adjusted.values = p.adjust(t.split.stats$frame.pvalue, method="bonferroni")
t.split.stats$adjusted.pvalue <-adjusted.values
write.table(t.split.stats, file=paste0(dataset, "/traj_heatmap_anova.tsv"), quote=FALSE, sep='\t', col.names=NA)


## Generates 3D UMAP plots of our data.
umap.data <- read.table(binned.umap.in.path)
names(umap.data) <- c("umap1", "umap2", "umap3", "utr", "traj", "age")
# Adjust these cutoffs as fits your data, we put ours near our 25th and 75th percentiles due to point density.
cutoff1 <- -50
cutoff2 <- 50
umap.data$utr[umap.data$utr <= cutoff1] <- cutoff1
umap.data$utr[umap.data$utr >= cutoff2] <- cutoff2
ncmt <- umap.data[umap.data$traj=="Neural_crest_melanocytes_trajectory",]
ncmt <- arrange(ncmt, age)
ncmt <- ncmt[ncmt$umap1<1500,]
ncmt$age <- as.character(ncmt$age)
ntant <- umap.data[umap.data$traj=="Neural_tube_and_notochord_trajectory",]
ntant <- arrange(ntant, age)
ntant$age <- as.character(ntant$age)
ht <- umap.data[umap.data$traj=="Haematopoiesis_trajectory",]
ht <- arrange(ht, age)
ht$age <- as.character(ht$age)
# Plots Neural Crest and Melanocytes. The initial axis and marker settings remain the same for each plot.
axis <- list(zeroline=FALSE, showline=FALSE, showticklabels=FALSE, title="", showgrid=FALSE)
markers <- list(size=2, color=~utr, colorscale='Viridis', showscale=TRUE)
viridis.discrete <- c("#3b528bff", "#21918cff", "#5ec962ff", "#fde725ff", "#440154ff")
# The camera (and therefore scene) setting, however, changes for each plot.
# Our camera settings give the initial view we use for images in the publication, but the actual generated
# files are manually manipulated 3d plots.
camera <- list(eye=list(x=-1.670061125143662, y=-1.213551691987196, z=0.6524477980336022),
			   up=list(x=-0.10432230768447087, y=0.5785027391444107, z=0.8089817284227132))
scene <- list(xaxis=axis, yaxis=axis, zaxis=axis, camera=camera)
plot <- plotly::plot_ly(ncmt, x=~umap1, y=~umap2, z=~umap3, marker=markers, text=~utr) %>% 
        layout(scene=scene) %>%
        add_markers()
widget.path <- paste0(working.path, "/", dataset, "/figures/Neural_crest_PNS_trajectory_umap.html")
htmlwidgets::saveWidget(as_widget(plot), widget.path)
# Plots Neural Tube and Notochord.
camera <- list(eye=list(x=-1.5037249407487514, y=1.3862811297253783, z=-0.7103069279807821),
	           up=list(x=-0.15862394497834553, y=0.30850374208857667, z=0.93790398506289))
scene <- list(xaxis=axis, yaxis=axis, zaxis=axis, camera=camera)
plot <- plotly::plot_ly(ntant, x=~umap1, y=~umap2, z=~umap3, marker=markers, text=~utr) %>% 
                layout(scene=scene) %>%
                add_markers()
widget.path <- paste0(working.path, "/", dataset, "/figures/Neural_tube_and_notochord_trajectory_umap.html")
htmlwidgets::saveWidget(as_widget(plot), widget.path)
# Plots Hematopoiesis.
camera <- list(eye=list(x=-0.7952917709395577, y=-1.000999029613577, z=1.7472870232988333))
scene <- list(xaxis=axis, yaxis=axis, zaxis=axis, camera=camera)
plot <- plotly::plot_ly(ht, x=~umap1, y=~umap2, z=~umap3, marker=markers, text=~utr) %>% 
        layout(scene=scene) %>%
        add_markers()
widget.path <- paste0(working.path, "/", dataset, "/figures/Haematopoiesis_trajectory_umap.html")
htmlwidgets::saveWidget(as_widget(plot), widget.path)



## Generates Unbinned 3D Umap plots of our data's ages.
umap <- data.frame(data$rumap1, data$rumap2, data$rumap3, data$age, data$trajectory, data$utr)
names(umap) <- c("umap1", "umap2", "umap3", "age", "traj", "utr")
cutoff1 <- -50
cutoff2 <- 50
umap$utr[umap$utr <= cutoff1] <- cutoff1
umap$utr[umap$utr >= cutoff2] <- cutoff2
ncmt <- umap[umap$traj=="Neural_crest_melanocytes_trajectory",]
ncmt <- ncmt[complete.cases(ncmt),]
ncmt <- ncmt[ncmt$umap1<1.5,]
ntant <- umap[umap$traj=="Neural_tube_and_notochord_trajectory",]
ntant <- ntant[complete.cases(ntant),]
ht <- umap[umap$traj=="Haematopoiesis_trajectory",]
ht <- ht[complete.cases(ht),]
# Resets our markers to have a discrete colorscale and smaller points. Our axes remain the same.
markers <- list(size=1, color=~age, showscale=TRUE, colorscale="Viridis")
# Neural Crest Melanocytes
camera <- list(eye=list(x=-1.670061125143662, y=-1.213551691987196, z=0.6524477980336022),
			   up=list(x=-0.10432230768447087, y=0.5785027391444107, z=0.8089817284227132))
scene <- list(xaxis=axis, yaxis=axis, zaxis=axis, camera=camera)
plot <- plotly::plot_ly(ncmt, x=~umap1, y=~umap2, z=~umap3, marker=markers, text=~age) %>% 
        layout(scene=scene) %>%
        add_markers()
widget.path <- paste0(working.path, "/", dataset, "/figures/Neural_crest_PNS_trajectory_age.html")
htmlwidgets::saveWidget(as_widget(plot), widget.path)
# Neural Tube and Notochord
camera <- list(eye=list(x=-1.5037249407487514, y=1.3862811297253783, z=-0.7103069279807821),
	           up=list(x=-0.15862394497834553, y=0.30850374208857667, z=0.93790398506289))
scene <- list(xaxis=axis, yaxis=axis, zaxis=axis, camera=camera)
plot <- plotly::plot_ly(ntant, x=~umap1, y=~umap2, z=~umap3, marker=markers, text=~age) %>% 
        layout(scene=scene) %>%
        add_markers()
widget.path <- paste0(working.path, "/", dataset, "/figures/Neural_tube_and_notochord_trajectory_age.html")
htmlwidgets::saveWidget(as_widget(plot), widget.path)
# Haematopoiesis
camera <- list(eye=list(x=-0.7952917709395577, y=-1.000999029613577, z=1.7472870232988333))
scene <- list(xaxis=axis, yaxis=axis, zaxis=axis, camera=camera)
plot <- plotly::plot_ly(ht, x=~umap1, y=~umap2, z=~umap3, marker=markers, text=~age) %>% 
        layout(scene=scene) %>%
        add_markers()
widget.path <- paste0(working.path, "/", dataset, "/figures/Haematopoiesis_trajectory_age.html")
htmlwidgets::saveWidget(as_widget(plot), widget.path)

stop()
