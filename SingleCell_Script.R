#Single Cell RNA expression evaluation
#tissue: neo-cortex (SRA667466)
#species: mouse
library(Seurat)
library(dplyr)
library(patchwork) 
setwd("D:\\OneDrive - Politecnico di Milano\\PoliMirco\\Primo anno\\Secondo semestre\\Genomics and Transcriptomics\\Trascriptomics\\TransProject\\SingleCell")
load("SRA667466_SRS3059977.sparse.RData")
head(sm@Dimnames[[1]],2000)
rownames(sm) <- gsub("[_].*", replacement = "", x = rownames(sm))
#prefiltering 
neocortex <- CreateSeuratObject(counts = sm, project = "NeoC3k",
                             min.cells=10, min.features=200)#min genes 200, min cells 3 
neocortex # only those genes which are found to be expressed count > 0 in  at least 3 cells are kept
#fish out rows starting with "mt" and compute the percentage column by column
neocortex[["percent.mt"]] <- PercentageFeatureSet(neocortex, pattern = "^mt-")
head(neocortex)#1st column: name project, 2nd:overall count of that cell 3rd: number of genes expressed by that cell

#---- VISUALIZATION and Quality Control
x11()
VlnPlot(neocortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
x11()
plot1 <- FeatureScatter(neocortex, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(neocortex, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#--- FILTERING
dim(neocortex)
neocortex <- subset(neocortex, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10 & nCount_RNA < 9000)
dim(neocortex)

#---Normailzation
neocortex <- NormalizeData(neocortex, normalization.method = "LogNormalize", scale.factor = 10000)
#choose the more meaningful genes: most variable genes wrt average across cells
neocortex <- FindVariableFeatures(neocortex, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(neocortex), 10)
top10
# plot variable features with and without labels

plot1 <- VariableFeaturePlot(neocortex)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1
#plot2
# Visualization after filtering
x11()
VlnPlot(neocortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
#neocortex[["RNA"]]@#explore

# ---SCALING
all.genes <- rownames(neocortex)
neocortex <- ScaleData(neocortex, features = all.genes)
#head(neocortex[["RNA"]]@scale.data[1,])  

#neocortex <- ScaleData(neocortex, vars.to.regress = "percent.mt")
#mt genes regressed out because interfering with DE genes analysis!

#s phase and g2m phase are problematic, address them:
cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
neocortex <- CellCycleScoring(neocortex, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(neocortex)


#PCA
neocortex <- RunPCA(neocortex, features = VariableFeatures(object = neocortex))
neocortex[["pca"]] # MOST VARIABLE FREATURES
print(neocortex[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(neocortex, dims = 1:2, reduction = "pca")
DimPlot(neocortex, reduction = "pca")
#cell cycle seems not to be a key factor in grouping of the cells,
# hence, we don't have to account for the cell cycle phase 
DimPlot(neocortex, reduction = "pca", dims = c(3,4))
DimHeatmap(neocortex, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(neocortex, dims = 1:9, cells = 500, balanced = TRUE)

# how many components do i need to choose? elbow plot is a intuitive solution
#plot standard deviation for each principal component
x11()
ElbowPlot(neocortex, ndims = 50)

neocortex10  <- FindNeighbors(neocortex, dims = 1:10)
neocortex15 <- FindNeighbors(neocortex, dims = 1:15)
#neocortex16 <- FindNeighbors(neocortex, dims = 1:16)
neocortex20  <- FindNeighbors(neocortex, dims = 1:20)
neocortex30 <- FindNeighbors(neocortex, dims = 1:30)
neocortex40 <- FindNeighbors(neocortex, dims = 1:40)
neocortex50 <- FindNeighbors(neocortex, dims = 1:50)

neocortex10.2 <- FindNeighbors(neocortex40, dims = 1:10)
neocortex10.3 <- FindNeighbors(neocortex10, dims = 1:10)

neocortex10 <- FindClusters(neocortex10, resolution = 0.5) 
summary(neocortex40@meta.data$seurat_clusters)
neocortex15 <- FindClusters(neocortex15, resolution = 0.5)
#neocortex16 <- FindClusters(neocortex16, resolution = 0.5)
neocortex20 <- FindClusters(neocortex20, resolution = 0.5) 
neocortex30 <- FindClusters(neocortex30, resolution = 0.5)
neocortex40 <- FindClusters(neocortex40, resolution = 0.5)
neocortex50 <- FindClusters(neocortex50, resolution = 0.35)
#summary(neocortex30@meta.data$seurat_clusters)
neocortex10.2 <-FindClusters(neocortex10.2, resolution = 0.5)
summary(neocortex10.3@meta.data$seurat_clusters)
# higher resolution(0.4-1.2): higher number of cluster and viceversa
neocortex10.3 <- FindClusters(neocortex10.3, resolution = 0.5) 

#x11()
#DimPlot(neocortexother, reduction = "pca")
x11()
DimPlot(neocortex10.2, reduction = "pca")
#x11()
#DimPlot(neocortex15.1, reduction = "pca")
#x11()
#DimPlot(neocortex20, reduction = "pca")


#best choice is to use UMAP
#best was neo10 + 1:10 dims
neocortex <- RunUMAP(neocortex10, dims = 1:10)
neocortex30 <- RunUMAP(neocortex30, dims = 1:30)
neocortex40 <- RunUMAP(neocortex40, dims = 1:40)
neocortex50 <- RunUMAP(neocortex50, dims = 1:50)
neocortex10.2 <- RunUMAP(neocortex10.2, dims = 1:10)
neocortex10.3 <- RunUMAP(neocortex10.3, dims = 1:10)


DimPlot(neocortex40, reduction = "umap")


x11()
DimPlot(neocortex40, reduction = "umap")
#saveRDS(neocortex, file = "./neocortex_tutorial.rds")

#Find DE genes
cluster1.markers <- FindMarkers(neocortex, ident.1 = 1, min.pct = 0.25)#use higher min.pct
head(cluster1.markers, n = 5)


#pct1 = cells in cluster 1 expressing this gene
#pct 2 = cells in other clusters expressing this gene

neocortex.mark_wilcox <- FindAllMarkers(neocortex, only.pos = TRUE, 
                                        min.pct = 0.25, min.diff.pct = 0.5, logfc.threshold = 0.25)
neocortex.mark_wilcox

neocortex.mark40<- FindAllMarkers(neocortex40, only.pos = TRUE, 
                                        min.pct = 0.25,  logfc.threshold = 0.25)

neocortex.mark10.2<- FindAllMarkers(neocortex10.2, only.pos = TRUE, 
                                  min.pct = 0.25,  logfc.threshold = 0.25)

#neocortex.mark_bimod <- FindAllMarkers(neocortex, only.pos = TRUE, 
                                    #   min.pct = 0.25, logfc.threshold = 0.25,
                                    #   test.use = "bimod", verbose = TRUE)
#neocortex.mark_bimod

#neocortex.mark_poisson <-  FindAllMarkers(neocortex, only.pos = TRUE, 
 #                                        min.pct = 0.25, logfc.threshold = 0.25,
  #                                       test.use = "poisson", verbose = TRUE)
#neocortex.mark_poisson

#neocortex.mark_negbinom <-FindAllMarkers(neocortex, only.pos = TRUE, 
                                    #     min.pct = 0.25, logfc.threshold = 0.25,
                                    #     test.use = "negbinom", verbose = TRUE)
#neocortex.mark_negbinom

neocortex.mark_wilcox %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)

neocortex.mark10.2 %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)

kl <- neocortex.mark40 %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)


#neocortex.mark_wilcox %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

#neocortex.mark_wilcox %>% group_by(cluster) %>% top_n(n = 5, wt = p_val)


#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5742025/

#marker_genes <- c("Slc17a7","Gad1", "Olig1","Aldoc", "Cldn5", "Vtn", "Cx3cr1","Aif1","Mrc1")
#11()
#VlnPlot(neocortex, features = marker_genes)
#remove Mrc1
marker_genes <- c("Slc17a7","Gad1","Pvalb","Olig1","Aldoc", "Cldn5", "Vtn", "Cx3cr1","Plp1","Synpr","Aif1")
x11()
VlnPlot(neocortex, features = marker_genes, pt.size = 0)

x11()
VlnPlot(neocortex, features = c("Mbp",          "Aldoc"  ,      "Rtn1"      ,   "Atp1a1"      , "RP23-81C12.1" ,"Mfge8",        "Sst" ,"Cldn5"    ,   "hexb"      ,   "Npy"     ,     "C1ql1"      ,  "Ctss"     ,    "Plp1" ,        "Vtn"  ,
                                "Lyz2") , pt.size = 0)

x11()
VlnPlot(neocortex, features = c(top_genes$gene) , pt.size = 0)


x11()
#top10 <- neocortex %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(neocortex, features = marker_genes) + NoLegend()

x11()
FeaturePlot(neocortex, features = marker_genesNCBI)
#FeaturePlot(neocortex, features = "Mbp")
#Aldoc was unknown -> astrocytes soma


new.cluster.ids <- c("1", "2", "3", "4", "5", "6", 
                     "7", "8", "9","10","11",
                     "12","13","14","15")

new.cluster.ids <- c("Inh.Neuron", "Inh.Neuron", "Astrocytes", "Excitatory neurons", "Astrocytes", "H.cells(Neuron)", 
                     "Excitatory neurons", "Microglia", "Endothelial Cells","Oligo1","Oligo2",
                     "Macrophages","Excitatory neurons","Endothelial Cells","Unknown")

names(new.cluster.ids) <- levels(neocortex)

neocortex <- RenameIdents(neocortex, new.cluster.ids)
x11()
DimPlot(neocortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


