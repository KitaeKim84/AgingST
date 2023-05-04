###After preprocessing with space ranger 
wt_old <- Load10X_Spatial("/data2/Visium/20230322_visium/preprocessing/WT_old/outs/",
                        filename="filtered_feature_bc_matrix.h5",slice="s3")
wt_old<-RenameCells(object=wt_old,new.names=paste(names(Idents(wt_old)),"-S",3,sep=""))

wt_young <- Load10X_Spatial("/data2/Visium/20230322_visium/preprocessing/WT_young/outs/",
                        filename="filtered_feature_bc_matrix.h5",slice="s4")
wt_young<-RenameCells(object=wt_young,new.names=paste(names(Idents(wt_young)),"-S",4,sep=""))

brain<-merge( wt_old, wt_young)
dim(brain)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT-")
brain <- subset(brain, subset = nFeature_Spatial > 200 & nFeature_Spatial < 2500 & percent.mt < 5)
brain <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(brain)
brain <- ScaleData(brain, features = all.genes)
brain <- RunPCA(brain, features = VariableFeatures(object = brain))
brain <- JackStraw(brain, num.replicate = 100)
brain <- ScoreJackStraw(brain, dims = 1:20)
brain <- FindNeighbors(brain, dims = 1:10)
brain <- FindClusters(brain, resolution = 0.5)
brain <- RunUMAP(brain, dims = 1:10)
brain <- RunTSNE(brain, dims = 1:10)
names(brain@images)<-c("Old","Young")

pal_d3("category10", alpha = 0.7)(10)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE, pt.size=1,label.size=7) 
p2 <- DimPlot(brain, reduction = "tsne", label = TRUE, pt.size=1,label.size=7) 

p3 <- SpatialDimPlot(brain, label = F, label.size = 7)

DimPlot(brain, reduction = "umap", label = F, pt.size=1.5,label.size=8,cols=pal_d3("category10", alpha = 1)(10)) +
  theme(legend.position="top", 
        legend.direction="horizontal",
        legend.text=element_text(size=20,face="bold"),
        legend.title = element_blank(),
        axis.text=element_text(size=20),
        axis.text.x = element_text(size = 20),
        axis.title=element_text(size=20,face="bold"),
        plot.margin = margin(2, 2, 2, 3, "cm")) 

##### cluster markers #####
total_exp<-as.matrix(GetAssayData(brain))
org.markers <- FindAllMarkers(object = brain, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

org.markers$gene<-rownames(org.markers)
summary(as.factor(org.markers$cluster))

sum.exp<-apply(total_exp,1,sum)
not_zero_gene<-names(sum.exp[sum.exp!=0])
markers <- org.markers[ org.markers$p_val_adj < 0.05 & org.markers$avg_log2FC>log2(1.5), ]

markers<-markers%>%filter(gene%in%intersect(not_zero_gene, markers$gene))
markers<-markers[!grepl("^.*Rik+|^Gm[0-9]+$|^AC[0-9]+.[0-9]$|^BC[0-9]+|ps-|mt-|-ps",markers$gene), ]
summary(as.factor(markers$cluster))

top.markers <-  markers %>%
  group_by(cluster) %>%
  arrange(avg_log2FC, .by_group = TRUE) %>%
  top_n(5) %>% data.frame()

DoHeatmap(brain, features = top.markers$gene)+ scale_fill_viridis() + 
theme(axis.line = element_line(colour = 'black', size = 0.5, linetype = "solid"),
          axis.text=element_text(size=13,face="bold",colour="black"),#gene symbol size
          axis.text.x = element_text(size = 15,face="bold",colour="black"),
          legend.title=element_text(size=15,face="bold"),
          legend.text=element_text(size=15)) 



### basic statistics for clusters

s.data<-brain
# total statistics
a1<-summary(Idents(s.data))
a2<-round(summary(Idents(s.data))/sum(summary(Idents(s.data)))*100,2)
paste(a1," (",a2,"%)",sep="")


# sample wise statistics
cell_info<-data.frame(cluster=Idents(s.data), sample=sapply(strsplit(as.character(names(Idents(s.data))), split='-'), "[", 3))
cell_info %>% 
  group_by(cluster,sample) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

library(moonBook)
out=mytable(sample~cluster,data=cell_info)
out

## line plot
##cluster
library(RColorBrewer)

cell.table<-as.data.frame(table(cell_info$cluster, cell_info$sample))
cell.table<-data.frame(cell.table%>%group_by(Var2)%>%mutate(prop = 100*Freq / sum(Freq)) )
cell.table$Var2<-factor(cell.table$Var2,
                        levels=c("S3","S4"),
                        labels=c("Old","Young"))

##barplot
ggplot() + theme_bw() +
  geom_bar(aes(y = prop, x = Var1, fill = Var2), data = cell.table, stat="identity") +
  scale_fill_brewer(palette="Pastel1") +
  theme(legend.position="top", 
        legend.direction="horizontal",
        legend.text=element_text(size=30,face="bold"),
        legend.title = element_blank(),
        axis.text=element_text(size=30,color="black"),
        axis.text.x = element_text(size = 30,color="black"),
        axis.title=element_text(size=30,face="bold"),
        plot.margin = margin(2, 2, 2, 3, "cm")) +
  labs(x="Clusters", y="Proportion (%, in each sample)") 

####### Slingshot : trajectory 
#ref: https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html
library(SingleCellExperiment)
library(slingshot)


brain_sce <- as.SingleCellExperiment(brain)

#gene filter
geneFilter <- apply(assays(brain_sce)$counts,1,function(x){
    sum(x >= 3) >= 10
})
brain_sce <- brain_sce[geneFilter, ]

#normalization
FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
}

assays(brain_sce)$norm <- FQnorm(assays(brain_sce)$counts)

#dimensional reduction
pca <- prcomp(t(log1p(assays(brain_sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

library(uwot)

rd2 <- uwot::umap(t(log1p(assays(brain_sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')
plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)


reducedDims(brain_sce) <- SimpleList(PCA = rd1, UMAP = rd2)

#clustering cells
library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(brain_sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)


cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

brain_sce <- slingshot(brain_sce, clusterLabels = 'GMM', reducedDim = 'PCA')

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(brain_sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(brain_sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(brain_sce), lwd=2, col='black')

plot(reducedDims(brain_sce)$PCA, col = brewer.pal(9,'Set1')[brain_sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(brain_sce), lwd=2, type = 'lineages', col = 'black')


######## trajectory : tradeseq #########
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)

seu<-brain

## sds = crv
crv <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters, start.clus = 4, stretch = 0)
brain.sce <- as.SingleCellExperiment(brain)
brain.counts<-GetAssayData(object = brain, slot = "counts")
brain.crv<-readRDS("brain_crv.rds")


cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}


r_col=pal_d3("category10", alpha = 0.7)(10)

library(scales)
cell_colors <- cell_pal(Idents(seu), pal_d3("category10", alpha = 0.7))
cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())

plot(reducedDim(brain.crv), col = cell_colors, pch = 20, cex = 1)
lines(brain.crv, lwd = 2, type = 'lineages', col = 'black')

pseudotime <- slingPseudotime(brain.crv, na = FALSE)
cellWeights <- slingCurveWeights(brain.crv)
sce <- fitGAM(counts = brain.counts, pseudotime = pseudotime, cellWeights = cellWeights,
                 nknots = 6, verbose = FALSE)
table(rowData(sce)$tradeSeq$converged)

assoRes <- associationTest(sce)
head(assoRes)

startRes <- startVsEndTest(sce)

assoRes<-readRDS("assoRes.rds")
startRes<-readRDS("startRes.rds")

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, brain.counts, gene = sigGeneStart)+
        theme(axis.text=element_text(size=20),
              axis.text.x = element_text(size = 20),
              axis.title=element_text(size=20,face="bold"),
              legend.title=element_text(size=20,face="bold"),
              legend.text=element_text(size=20,face="bold"),
              plot.margin = margin(2, 2, 2, 3, "cm")) 

########### cell cell communication #########
#ref: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
library(CellChat)
library(patchwork)
library(Future)
options(stringsAsFactors = FALSE)

#seurat.data <- subset(brain, subset = orig.ident == "Aged")
seurat.data <- subset(brain, subset = orig.ident == "Young")
data.input <- GetAssayData(seurat.data, , slot = "counts") # normalized data matrix
labels <- Idents(seurat.data)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

#CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
### project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- projectData(cellchat, PPI.mouse)

### long time consume
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

####signaling role identification
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- rankNetPairwise(cellchat)


#--------------------------------------------------------------------------------------------------------------
###Perform NicheNet analysis starting from a Seurat object
#--------------------------------------------------------------------------------------------------------------
library(nichenetr)
library(Seurat)
library(tidyverse)
library(readr)
library(ggh4x)
library(cowplot)

ligand_target_matrix = readRDS("/exdisk/sda1/database/nichenet/ligand_target_matrix.rds")
lr_network = readRDS("/exdisk/sda1/database/nichenet/lr_network.rds")
weighted_networks = readRDS("/exdisk/sda1/database/nichenet/weighted_networks.rds")

nichenet_output.old = nichenet_seuratobj_aggregate(seurat_obj = brain, 
  receiver = c("MB","HPF", "isoCTX", "PIR", "NF", "VEN","TH"),
  condition_colname = "orig.ident", condition_oi = "Old", condition_reference = "Young", 
  sender = "all", 
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "mouse")

nichenet_output.young = nichenet_seuratobj_aggregate(seurat_obj = brain, 
  receiver = c("MB","HPF", "isoCTX", "PIR", "NF", "VEN","TH"),
  condition_colname = "orig.ident", condition_oi = "Young", condition_reference = "Old", 
  sender = "all", 
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "mouse")

nichenet_output.old$ligand_expression_dotplot + 
theme(axis.text=element_text(size=20),
              axis.text.x = element_text(size = 20),
              axis.title=element_text(size=20,face="bold"),
              legend.title=element_text(size=20,face="bold"),
              legend.text=element_text(size=20,face="bold")) 

nichenet_output.young$ligand_expression_dotplot + 
theme(axis.text=element_text(size=20),
              axis.text.x = element_text(size = 20),
              axis.title=element_text(size=20,face="bold"),
              legend.title=element_text(size=20,face="bold"),
              legend.text=element_text(size=20,face="bold")) 

