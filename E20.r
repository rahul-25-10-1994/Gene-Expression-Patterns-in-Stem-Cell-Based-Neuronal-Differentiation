library(Seurat)

# Load your data (replace with your actual data path)
setwd("C:/Users/panka/OneDrive/Desktop/adipose/CYTOSCAPE FILES/der/ESC/DAY20/SC1")
cts1<-ReadMtx(mtx="matrix.mtx.gz",                             
              features = "features.tsv.gz",
              cells = "barcodes.tsv.gz")
setwd("C:/Users/panka/OneDrive/Desktop/adipose/CYTOSCAPE FILES/der/ESC/DAY20/SC2")
cts2<-ReadMtx(mtx="matrix.mtx.gz",                             
              features = "features.tsv.gz",
              cells = "barcodes.tsv.gz")

seurat_object1 <- CreateSeuratObject(counts = cts1, project = "GeneExpressionAnalysis")
seurat_object2 <- CreateSeuratObject(counts = cts2, project = "GeneExpressionAnalysis")


merged_seurat_object <- merge(seurat_object1, 
                              y = c(seurat_object2), 
                              add.cell.ids = c("SC1", "SC2"))

# Quality Control for the merged object
merged_seurat_object[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_object, pattern = "^MT-")
small_seurat <- subset(merged_seurat_object, cells = sample(colnames(merged_seurat_object), 1000))
VlnPlot(merged_seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
FeatureScatter(merged_seurat_object,feature1 = "nCount_RNA",feature2 = "nFeature_RNA") + geom_smooth(method='lm')

# Filter cells based on quality metrics
merged_seurat_object <- subset(merged_seurat_object, 
                               subset = nFeature_RNA > 200 & 
                                 nFeature_RNA < 2500 & 
                                 percent.mt < 5 & 
                                 nCount_RNA > 800 & 
                                 nCount_RNA < 20000)

# Normalization
merged_seurat_object<- NormalizeData(merged_seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable features
merged_seurat_object <- FindVariableFeatures(merged_seurat_object, selection.method = "vst", nfeatures = 2000)
mito.genes <- grep("^MT-", VariableFeatures(merged_seurat_object), value = TRUE)
VariableFeatures(merged_seurat_object) <- setdiff(VariableFeatures(merged_seurat_object), mito.genes)

top10<-head(VariableFeatures(merged_seurat_object),10)
plot1<-VariableFeaturePlot(merged_seurat_object)
label_params<-list(max.overlap=100000, xnudge=0,ynudge=0)
LabelPoints(plot=plot1, points=top10,repel = TRUE,parameters=label_params,xnudge = 0,ynudge = 0)

# Scale the data
merged_seurat_object <- ScaleData(merged_seurat_object, vars.to.regress = "percent.mt")

# PCA
merged_seurat_object <- RunPCA(merged_seurat_object, features = VariableFeatures(object = merged_seurat_object))
DimHeatmap(object=merged_seurat_object,dims = 1:5)
ElbowPlot(merged_seurat_object)

# Clustering
merged_seurat_object <- FindNeighbors(merged_seurat_object, dims = 1:10)
merged_seurat_object <- FindClusters(merged_seurat_object, resolution = 0.5)

# UMAP for visualization
merged_seurat_object <- RunUMAP(merged_seurat_object, dims = 1:10)
DimPlot(merged_seurat_object, reduction = "umap")

merged_seurat_object <- JoinLayers(merged_seurat_object)
gene_of_int<-c( "SOX2", "BDNF")
VlnPlot(merged_seurat_object, features = gene_of_int, group.by = "seurat_clusters" )+ theme_minimal()

# Differential expression analysis
cluster_markers <- FindAllMarkers(merged_seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Visualize expression of specific marker genesAD
FeaturePlot(merged_seurat_object, features = c("NGF", "SOX2", "BDNF"))

# Save the merged Seurat object
saveRDS(merged_seurat_object, file = "C:/Users/panka/OneDrive/Desktop/ESC/escmerged20_seurat_object.rds")


