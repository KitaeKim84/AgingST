###After preprocessing with space ranger 

wt_old <- Load10X_Spatial("/data2/Visium/20230322_visium/preprocessing/WT_old/outs/",
                        filename="filtered_feature_bc_matrix.h5",slice="s3")
wt_old<-RenameCells(object=wt_old,new.names=paste(names(Idents(wt_old)),"-S",3,sep=""))

wt_young <- Load10X_Spatial("/data2/Visium/20230322_visium/preprocessing/WT_young/outs/",
                        filename="filtered_feature_bc_matrix.h5",slice="s4")
wt_young<-RenameCells(object=wt_young,new.names=paste(names(Idents(wt_young)),"-S",4,sep=""))


brain.wt.merge <- merge( wt_old, wt_young)
brain<-brain.wt.merge
dim(brain)


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT-")
#brain <- subset(brain, subset = nFeature_Spatial > 200 & nFeature_Spatial < 2500 & percent.mt < 5)

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


names(brain@images)<-c("Aged","Young")
saveRDS(brain, file = "brain_wt.rds")


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

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,
    5, 8)), facet.highlight = TRUE, ncol = 4)

#specific image selection (combine = FALSE)
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain), 
                facet.highlight = TRUE, ncol = 5, images="Aged",cols.highlight = c("#DE2D26","#00000000"))    

options(browser="google-chrome")


##### cluster markers #####
total_exp<-as.matrix(GetAssayData(brain))
org.markers <- FindAllMarkers(object = brain, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

org.markers$gene<-rownames(org.markers)
summary(as.factor(org.markers$cluster))


b<-apply(total_exp,1,sum)
not_zero_gene<-names(b[b!=0])
markers <- org.markers[ org.markers$p_val_adj < 0.05 & org.markers$avg_log2FC>log2(1.5), ]


markers<-markers%>%filter(gene%in%intersect(not_zero_gene, markers$gene))
markers<-markers[!grepl("^.*Rik+|^Gm[0-9]+$|^AC[0-9]+.[0-9]$|^BC[0-9]+|ps-|mt-|-ps",markers$gene), ]
summary(as.factor(markers$cluster))

markers<-read.table(file="fdr0.05_fc1.5_reging_markers_all.txt", header=T, sep="\t")

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


#####sc type
# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = brain[["Spatial"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(brain@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(brain@meta.data[brain@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(brain@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"


