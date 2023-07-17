# Raw data have been deposited to GEO under accession number GSE211381.
# The following objects and codes were used to analyze the data, and to generate graphs presented in the associated paper.
# Load all required packages
suppressMessages(library(reticulate))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(sctransform))
suppressMessages(library(harmony))
suppressMessages(library(pals))

# Create all the 10X mats ####
datasets.10X.directories <- c("PATH/filtered_feature_bc_matrix",
                              "PATH/filtered_feature_bc_matrix",
                              "PATH/filtered_feature_bc_matrix",
                              "PATH/filtered_feature_bc_matrix",
                              "PATH/filtered_feature_bc_matrix",
                              "PATH/filtered_feature_bc_matrix"
)

datasets.matrix <- list()

for(i in 1:length(datasets.10X.directories)){
  if(file.exists(datasets.10X.directories[[i]])){
    datasets.matrix[[i]] <- Read10X(data.dir = datasets.10X.directories[[i]])
  }else{
    cat("Data not found for", datasets.10X.directories[[i]], "\n\n")
  }
}

# add metadata ####
datasets.sample <- c('WT_Telo_1',
                     'WT_Telo_2',
                     'WT_Ana_1',
                     'WT_Ana_2',
                     'Alk1KO_Telo_1',
                     'Alk1KO_Telo_2')
datasets.stage<- c('WT_Telo',
                  'WT_Telo',
                  'WT_Ana',
                  'WT_Ana',
                  'Alk1KO_Telo',
                  'Alk1KO_Telo')

# Run seurat pipeline ####
n.pcs <- 20  
datasets.seurat.list <- list()
for(i in 1:length(datasets.matrix)){ #everything but D0_C, from Andrea's data
  cat(' ###################################\n',
      '### Processing dataset number ', i, '###\n',
      '###################################\n')
  
  # Initialize seurat object
  datasets.seurat.list[[i]] <- CreateSeuratObject(counts = datasets.matrix[[i]], 
                                                  min.features = 200,
                                                  min.cells = 5,
                                                  project = 'Endo_Data')
  # Add meta data
  datasets.seurat.list[[i]]$sample  <- datasets.sample[[i]]
  datasets.seurat.list[[i]]$stage <- datasets.stage[[i]]
  
  # add %MT...
  datasets.seurat.list[[i]][["percent.mt"]]  <- PercentageFeatureSet(datasets.seurat.list[[i]], pattern = "mt-") 
  
  # Filter out low quality cells according to the metrics defined above
  datasets.seurat.list[[i]] <- subset(datasets.seurat.list[[i]],
                                      subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &
                                        # nFeature_RNA < 2500 & 
                                        percent.mt < 10
  )
  
  DefaultAssay(datasets.seurat.list[[i]]) <- 'RNA'
  datasets.seurat.list[[i]] <- NormalizeData(datasets.seurat.list[[i]], 
                                             normalization.method = "LogNormalize", 
                                             scale.factor = 10000)
  datasets.seurat.list[[i]] <- FindVariableFeatures(object = datasets.seurat.list[[i]], 
                                                    assay='RNA',
                                                    selection.method = "vst", 
                                                    nfeatures = 2000, 
                                                    verbose = TRUE)
  datasets.seurat.list[[i]] <- ScaleData(datasets.seurat.list[[i]],
                                         assay='RNA')
  
  datasets.seurat.list[[i]] <- RunPCA(object = datasets.seurat.list[[i]],
                                      assay='RNA',
                                      reduction.name = 'pca_RNA',
                                      reduction.key= 'pca_RNA',
                                      verbose = TRUE,
                                      npcs = n.pcs
  )
  
  datasets.seurat.list[[i]] <- FindNeighbors(object = datasets.seurat.list[[i]], 
                                             reduction = "pca_RNA", 
                                             dims = 1:n.pcs, 
                                             k.param = 30, # <---
                                             force.recalc = TRUE,
                                             verbose = FALSE)
  
  datasets.seurat.list[[i]] <- FindClusters(object = datasets.seurat.list[[i]], 
                                            reduction = "pca_RNA", 
                                            resolution = 0.3)
  datasets.seurat.list[[i]][['RNA_snn_res.0.3']] <- as.vector(datasets.seurat.list[[i]]@active.ident)
  
  
  # Check number of cells for UMAP- make sure n.neighbors is less than num_cells
  if(dim(datasets.seurat.list[[i]])[2] < 30){
    n.neighbors <- 30
  }else{
    n.neighbors <- dim(datasets.seurat.list[[i]])[2] - 1
  }
  
  datasets.seurat.list[[i]] <- RunUMAP(object = datasets.seurat.list[[i]], 
                                       reduction = "pca_RNA", 
                                       dims = 1:n.pcs,
                                       n.neighbors = n.neighbors,
                                       reduction.name='umap_RNA')
  
  
  
}

object.seurat <- merge(
  datasets.seurat.list[[1]],
  y = c(datasets.seurat.list[[2]],
        datasets.seurat.list[[3]],
        datasets.seurat.list[[4]],
        datasets.seurat.list[[5]],
        datasets.seurat.list[[6]]
  ),
  add.cell.ids = datasets.sample
)

# QC.vln <- 
VlnPlot(
  Total_Telo_Ana,
  features = c(
    'nCount_RNA',
    'nFeature_RNA',
    'percent.mt'
  ),
  group.by = 'sample',
  # group.by = 'chemistry',
  # y.max = 250,
  # log = TRUE,
  pt.size = 0#.001
)

# Seurat non-integrated data ####
DefaultAssay(object.seurat) <- 'RNA'
object.seurat <- FindVariableFeatures(object.seurat,
                                      assay = 'RNA',
                                      selection.method = 'vst',
                                      nfeatures = 2000,
                                      verbose=TRUE)

object.seurat <- ScaleData(object.seurat,
                           assay = 'RNA',
                           verbose=TRUE)

object.seurat <- RunPCA(object = object.seurat,
                        assay='RNA',
                        #reduction.name = 'pca_RNA',
                        reduction.key= 'pca_RNA',
                        verbose = TRUE, 
                        npcs = n.pcs
)

# Harmony Integration
TotalEndo <- object.seurat %>% 
  RunHarmony(group.by.vars=c('sample'), 
             assay='RNA',               
             plot_convergence = TRUE,
             verbose=TRUE)

# Re-run most of the seurat pipeline
ElbowPlot(TotalEndo)
n.pcs = 20

TotalEndo <- FindNeighbors(object = TotalEndo, 
                               reduction = 'harmony',
                               dims = 1:n.pcs, 
                               k.param = 30, 
                               force.recalc = TRUE,
                               verbose = FALSE)

TotalEndo <- FindClusters(object = TotalEndo, 
                              reduction = 'harmony',
                              resolution = 0.3)
TotalEndo[['harmony_res.0.3']] <- as.vector(TotalEndo@active.ident) 

TotalEndo <- RunUMAP(object = TotalEndo, 
                         reduction = 'harmony',
                         dims = 1:n.pcs,
                         reduction.name='umap_harmony'
)

# Visualization

DimPlot(TotalEndo, reduction = "umap_harmony", label = TRUE, label.size = 5, pt.size = 1) + FontSize(x.title = 12, y.title = 12)

TotalIdentMarkers <- c("Cdh5", "Pecam1", "Cd34", "Emcn", "Kdr", "Cxcl12", "Cd36", "Efnb2", "Gja4", "Eln", "Nr2f2", "Klf5", "Il33", "Prox1", "Rgcc", "Ccl21a", "Pdpn", "Lyve1", "Vwf", "Ackr1", "Selp", "Mki67", "Top2a", "Acta2", "Rgs4", "Pdgfrb", "Des", "Pdgfra", "Pcolce2", "Mmp2", "Loxl1", "Dpt", "Mbp", "Plp1", "Sox10", "Ptprc", "Cd74", "Cd14", "Krt5", "Krt14", "Krt15")

DotPlot(TotalEndo, features = TotalIdentMarkers) + RotatedAxis()

# Add cluster cell types to the metadata

cluster_cell_types <- RenameIdents(TotalEndo, '0' = "Endo1", '1' = "Endo2", '2' = "Perineurial cells", '3' = "Endo3", '4' = "Immune1", '5' = "Endo4", '6' = "Endo5", '7' = "Endo6", '8' = "Fibroblasts", '9' = "vSMC/Pericytes", '10' = "Immune2", '11' = "Immune3", '12' = "Immune4", '13' = "Immune5", '14' = "Endo7", '15' = "Schwann cells") 

# Add the information above to the metadata

TotalEndo$metadata_column_name <- Idents(cluster_cell_types)

# Check to make sure the clusters were labeled correctly

head(x = TotalEndo@meta.data, 50)

# Arrange the order of clusters
Totnew.cluster.ids <- c("Endo1", "Endo2", "Perineurial cells", "Endo3", "Immune1", "Endo4", "Endo5", "Endo6", "Fibroblasts", "vSMC/Pericytes", "Immune2", "Immune3", "Immune4", "Immune5", "Endo7", "Schwann cells")
names(Totnew.cluster.ids) <- levels(TotalEndo)
TotalEndo <- RenameIdents(TotalEndo, Totnew.cluster.ids)

# Subset for Total_Telo_Ana
# For Total_Telo_Ana
Idents(TotalEndo) <- TotalEndo$stage
Total_Telo_Ana <- subset(TotalEndo, idents = c("WT_Telo", "WT_Ana"), invert = FALSE)
Idents(Total_Telo_Ana) <- "metadata_column_name"
Idents(Total_Telo_Ana) <- factor(Idents(Total_Telo_Ana), levels= my_levels3)
DimPlot(Total_Telo_Ana, reduction = "umap_harmony", label = FALSE, split.by = "sample", label.size = 4, pt.size = 0.6, cols = c("red", "yellow3", 'orange', 'purple', 'hotpink', 'brown', 'green', "blue", "dimgrey", "cornsilk3", "cornsilk3", "cornsilk3", "cornsilk3", "cornsilk3", "cyan2", "black")) + FontSize(x.title = 12, y.title = 12)

# For Total_KO_vs_WT_Telo
Idents(TotalEndo) <- TotalEndo$stage
Total_WT_KO <- subset(TotalEndo, idents = c("WT_Telo", "Alk1KO_Telo"), invert = FALSE)
Idents(Total_WT_KO) <- "metadata_column_name"
Idents(Total_WT_KO) <- factor(Idents(Total_WT_KO), levels= my_levels3)
DimPlot(Total_WT_KO, reduction = "umap_harmony", label = TRUE, label.size = 4, pt.size = 0.8, cols = c('orange', "red", 'brown', 'purple', 'hotpink', "yellow3", 'green', "blue", "dimgrey", "cornsilk3", "black", "darkslategray4", "darkslategray4", "darkslategray4", "darkslategray4", "darkslategray4")) + FontSize(x.title = 12, y.title = 12)
DimPlot(Total_WT_KO, reduction = "umap_harmony", label = TRUE, label.size = 4, split.by = "stage", pt.size = 0.8, cols = c('orange', "red", 'brown', 'purple', 'hotpink', "yellow3", 'green', "blue", "dimgrey", "cornsilk3", "black", "darkslategray4", "darkslategray4", "darkslategray4", "darkslategray4", "darkslategray4")) + FontSize(x.title = 12, y.title = 12)

# Remove non endothelial cells for further endothelial specific analysis

Idents(TotalEndo) <- "metadata_column_name"
MatureEndo <- subset(x = TotalEndo, idents = c("Immune1", "Immune2", "Immune3", "Immune4", "Immune5", "Perineurial cells", "Fibroblasts", "vSMC/Pericytes", "Schwann cells"), invert = TRUE)
ncol(MatureEndo)
MatureEndo <-subset(x = MatureEndo, Ptprc <0.5)
ncol(MatureEndo)

# Check if the cluster is removed from the UMAP

DimPlot(MatureEndo, reduction = "umap_harmony", group.by = 'metadata_column_name', label = FALSE, pt.size = 0.6) + FontSize(x.title = 12, y.title = 12)

# Check number of remaining cells

ncol(MatureEndo)

#re-run most of the seurat pipeline (don't need to NormalizeData, ScaleData, or RunPCA)

ElbowPlot(MatureEndo)
n.pcs = 18

MatureEndo <- FindVariableFeatures(MatureEndo,
                                      assay = 'RNA',
                                      selection.method = 'vst',
                                      nfeatures = 2000,
                                      verbose=TRUE)

MatureEndo <- FindNeighbors(object = MatureEndo, 
                                   reduction = 'harmony',
                                   dims = 1:n.pcs, 
                                   k.param = 30, 
                                   force.recalc = TRUE,
                                   verbose = FALSE)

MatureEndo <- FindClusters(object = MatureEndo, 
                                  reduction = 'harmony',
                                  resolution = 0.35)
MatureEndo[['harmony_res.0.35']] <- as.vector(MatureEndo@active.ident) #Saves the cluster idents to the object

MatureEndo <- RunUMAP(object = MatureEndo, 
                             reduction = 'harmony',
                             dims = 1:n.pcs,
                             reduction.name='umap_harmony' # I always explicitly name my reductions, so things don't get confusing
)

DimPlot(MatureEndo, reduction = "umap_harmony", label = TRUE, label.size = 4, pt.size = 0.8) + FontSize(x.title = 12, y.title = 12)

# Rename clusters
new.cluster_cell_types <- RenameIdents(MatureEndo, '0' = "Capillary-Kdr", '1' = "Artery", '2' = "Lymph-Col", '3' = "Capillary-Rgcc", '4' = "Lymph-Cap", '5' = "Vein", '6' = "Proli-EC") 

# Add the information above to the metadata

MatureEndo$metadata_column_name <- Idents(new.cluster_cell_types)

# Check to make sure the clusters were labeled correctly

head(x = MatureEndo@meta.data, 50)


# For Telo_Ana
Idents(MatureEndo) <- MatureEndo$stage
Telo_Ana <- subset(MatureEndo, idents = c("WT_Telo", "WT_Ana"), invert = FALSE)
Idents(Telo_Ana) <- "metadata_column_name"
Idents(Telo_Ana) <- factor(Idents(Telo_Ana), levels= my_levels2)
DimPlot(Telo_Ana, reduction = "umap_harmony", label = FALSE, label.size = 4, pt.size = 0.5, cols = c("red", "yellow3", 'orange', 'purple', 'brown', 'hotpink', "green")) + FontSize(x.title = 12, y.title = 12)
DimPlot(Telo_Ana, reduction = "umap_harmony", label = FALSE, label.size = 4, split.by = "stage", pt.size = 0.5, cols = c("red", "yellow3", 'orange', 'purple', 'brown', 'hotpink', "green")) + FontSize(x.title = 12, y.title = 12)

# For KO_vs_WT_Telo
Idents(MatureEndo) <- MatureEndo$stage
KO_vs_WT_Telo <- subset(MatureEndo, idents = c("WT_Telo", "Alk1KO_Telo"), invert = FALSE)  
Idents(KO_vs_WT_Telo) <- "metadata_column_name"
Idents(KO_vs_WT_Telo) <- factor(Idents(KO_vs_WT_Telo), levels= my_levels2)
DimPlot(KO_vs_WT_Telo, reduction = "umap_harmony", label = FALSE, label.size = 4, split.by = "stage", pt.size = 0.8, cols = c("red", "yellow3", 'orange', 'purple', 'brown', 'hotpink', "green")) + FontSize(x.title = 12, y.title = 12)


# FindAllMarker, if you run this, you will get errors for other codes. so just skip it for now. and you can do this later.
Idents(MatureEndo) <- "metadata_column_name"
Idents(MatureEndo) <- factor(Idents(MatureEndo), levels= my_levels2)
Markers_MatureEndo <- FindAllMarkers(MatureEndo, test.use = "wilcox", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
Markers_MatureEndo %>% group_by(cluster) %>% top_n(5)
write.table(Markers_MatureEndo, "PATH/FILENAME.txt", sep="\t")

Idents(Telo_Ana) <- "metadata_column_name"
Idents(Telo_Ana) <- factor(Idents(Telo_Ana), levels= my_levels2)
Markers_Telo_Ana <- FindAllMarkers(Telo_Ana, test.use = "wilcox", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
Markers_Telo_Ana %>% group_by(cluster) %>% top_n(5)
write.table(Markers_Telo_Ana, "PATH/FILENAME.txt", sep="\t")

# Heatmap

DoHeatmap(subset(Total_Telo_Ana, downsample = 100), assay = "RNA", features = features, size = 5) + scale_fill_gradientn(colors = c("blue", "white", "red"))

avgexp.Total_Telo_Ana = AverageExpression(Total_Telo_Ana, assays = "RNA", return.seurat = T)

DoHeatmap(subset(avgexp.Total_Telo_Ana), features = features, size = 4, draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.y = element_text(size = 10) + coord_flip())


#stacked violin plot for visualizing single-cell data in Seurat (Check split.by function if you need)

library(Seurat)
library(patchwork)
library(ggplot2)


Idents(Total_Telo_Ana) <- "metadata_column_name"

Idents(Total_Telo_Ana) <- factor(Idents(Total_Telo_Ana))

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, flip = TRUE, 
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, split.by = "stage", pt.size = pt.size, cols = c('blue', "red"))+
    ylab(feature) + 
    theme(legend.position = "none",
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(0.8), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

StackedVlnPlot(obj = Total_Telo_Ana, features = features, assay = "RNA", pt.size = 0)