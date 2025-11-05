# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                          PCB Visium data                           -----
# -----                                                                    -----
# -----                           Chatterjee Lab                           -----
# -----                         University of Iowa                         -----
# -----                                                                    -----
# ------------------------------------------------------------------------------

# Summary: PCB Spatial Transcriptomics data
# Written by: Budhaditya Basu
# Date: 06/10/2024
# Revision: 10/01/2025

library(Seurat)
library(circlize)
library(ggthemes)
library(openxlsx)
library(patchwork)
library(tidyverse)
library(enrichplot)
library(scCustomize)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(clusterProfiler)

################################################################################
##############################PCB 1 ############################################
pcb.1.data.dir <- "/home/bbasu/LSS/lss_schatterj/PCB_visium_v1/PCB_1/outs/"
pcb.1 <- Load10X_Spatial(
  pcb.1.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "PCB_1")
head(pcb.1@meta.data)
pcb.1$orig.ident <- "PCB_1"

plot1 <- SpatialFeaturePlot(pcb.1, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("PCB 1")
VlnPlot(pcb.1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#Subset the biological replicate 1a
pcb.1a.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.1a.coord.csv") %>%
  pull(Barcode)
pcb.1a <- subset(pcb.1, cells = pcb.1a.cells)
pcb.1a$orig.ident <- "PCB_1a"
head(pcb.1a@meta.data)
plot2 <- SpatialFeaturePlot(pcb.1a, features = "nCount_Spatial")+ theme(legend.position = "right")

#Subset the biological replicate 1b
pcb.1b.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.1b.coord.csv") %>%
  pull(Barcode)
pcb.1b <- subset(pcb.1, cells = pcb.1b.cells)
pcb.1b$orig.ident <- "PCB_1b"
head(pcb.1b@meta.data)
plot3 <- SpatialFeaturePlot(pcb.1b, features = "nCount_Spatial")+ theme(legend.position = "right")

pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/1_PCB.1_Spatial_plot.pdf", height = 8, width = 12)
plot3/plot2|plot1
dev.off()
################################################################################
##############################PCB 2 ############################################
pcb.2.data.dir <- "/home/bbasu/LSS/lss_schatterj/PCB_visium_v1/PCB_2/outs/"
pcb.2 <- Load10X_Spatial(
  pcb.2.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "PCB_2")
head(pcb.2@meta.data)
plot4 <- SpatialFeaturePlot(pcb.2, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("PCB 2")
VlnPlot(pcb.2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#Subset the biological replicate 2a
pcb.2a.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.2a.coord.csv") %>%
  pull(Barcode)
pcb.2a <- subset(pcb.2, cells = pcb.2a.cells)
pcb.2a$orig.ident <- "PCB_2a"
head(pcb.2a@meta.data)
plot5 <- SpatialFeaturePlot(pcb.2a, features = "nCount_Spatial")+ theme(legend.position = "right")

#Subset the biological replicate 2b
pcb.2b.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.2b.coord.csv") %>%
  pull(Barcode)
pcb.2b <- subset(pcb.2, cells = pcb.2b.cells)
pcb.2b$orig.ident <- "PCB_2b"
head(pcb.2b@meta.data)
plot6 <- SpatialFeaturePlot(pcb.2b, features = "nCount_Spatial")+ theme(legend.position = "right")

pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/2_PCB.2_Spatial_plot.pdf", height = 8, width = 12)
plot5/plot6 | plot4
dev.off()
################################################################################
##############################Control 1 ########################################

control.1.data.dir <- "/home/bbasu/LSS/lss_schatterj/PCB_visium_v1/Control_1/outs/"
control.1 <- Load10X_Spatial(
  control.1.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Control_1")
head(control.1@meta.data)
plot7 <- SpatialFeaturePlot(control.1, features = "nCount_Spatial") + 
  theme(legend.position = "right")+ ggtitle("Control 1")
VlnPlot(control.1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#Subset the biological replicate 1a
control.1a.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_data/control.1a.coord.csv") %>%
  pull(Barcode)
control.1a <- subset(control.1, cells = control.1a.cells)
control.1a$orig.ident <- "Control_1a"
head(control.1a@meta.data)
plot8 <- SpatialFeaturePlot(control.1a, features = "nCount_Spatial")+ theme(legend.position = "right")

#Subset the biological replicate 1b
control.1b.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_data/control.1b.coord.csv") %>%
  pull(Barcode)
control.1b <- subset(control.1, cells = control.1b.cells)
control.1b$orig.ident <- "Control_1b"
head(control.1b@meta.data)
plot9 <- SpatialFeaturePlot(control.1b, features = "nCount_Spatial")+ theme(legend.position = "right")

pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/3_Control.1_Spatial_plot.pdf", height = 9, width = 12)
plot8/plot9|plot7
dev.off()
################################################################################
##############################Control 2 ########################################

control.2.data.dir <- "/home/bbasu/LSS/lss_schatterj/PCB_visium_v1/Control_2/outs/"
control.2 <- Load10X_Spatial(
  control.2.data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "Control_2")
head(control.2@meta.data)

plot10 <- SpatialFeaturePlot(control.2, features = "nCount_Spatial")  + 
  theme(legend.position = "right")+ ggtitle("Control 2")
VlnPlot(control.2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

#Subset the biological replicate 2a
control.2a.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_data/control.2a.coord.csv") %>%
  pull(Barcode)
control.2a <- subset(control.2, cells = control.2a.cells)
control.2a$orig.ident <- "Control_2a"
head(control.2a@meta.data)
plot11 <- SpatialFeaturePlot(control.2a, features = "nCount_Spatial")+ theme(legend.position = "right")

#Subset the biological replicate 2b
control.2b.cells <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_data/control.2b.coord.csv") %>%
  pull(Barcode)
control.2b <- subset(control.2, cells = control.2b.cells)
control.2b$orig.ident <- "Control_2b"
head(control.2b@meta.data)
plot12 <- SpatialFeaturePlot(control.2b, features = "nCount_Spatial")+ theme(legend.position = "right")

pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/4_Control.2_Spatial_plot.pdf", height = 9, width = 12)
plot11/plot12 | plot10
dev.off()
#===============================================================================
# Save the seurat objects
#===============================================================================
saveRDS(control.1a, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/control.1a.rds")
saveRDS(control.1b, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/control.1b.rds")
saveRDS(control.2a, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/control.2a.rds")
saveRDS(control.2b, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/control.2b.rds")
saveRDS(pcb.1a, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.1a.rds")
saveRDS(pcb.1b, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.1b.rds")
saveRDS(pcb.2a, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.2a.rds")
saveRDS(pcb.2b, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.2b.rds")


#===============================================================================
# Import the biological replicates
control.1a <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/control.1a.rds")
control.1b <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/control.1b.rds") 
control.2a <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/control.2a.rds")
control.2b <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/control.2b.rds")
pcb.1a <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.1a.rds")
pcb.1b <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.1b.rds")
pcb.2a <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.2a.rds")
pcb.2b <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/pcb.2b.rds")
#===============================================================================
# Normalization and Integration 
#===============================================================================

#Assign condition before integration
control.1a$treatment <- "Control"
control.1b$treatment <- "Control"
control.2a$treatment <- "Control"
control.2b$treatment <- "Control"
pcb.1a$treatment <- "PCB"
pcb.1b$treatment <- "PCB"
pcb.2a$treatment <- "PCB"
pcb.2b$treatment <- "PCB"

# Make a list
pcb.list <- list(control.1a,
                 control.1b,
                 control.2a,
                 control.2b,
                 pcb.1a,
                 pcb.1b,
                 pcb.2a,
                 pcb.2b)

# perform SCTransform normalization
for (i in 1:length(pcb.list)) {
  pcb.list[[i]] <- SCTransform(pcb.list[[i]], assay = 'Spatial')
}


# Create RNA assay for all spatial objects
for (i in 1:length(pcb.list)) {
  pcb.list[[i]][["RNA"]] <- pcb.list[[i]][["Spatial"]]
}


anchor_features <- SelectIntegrationFeatures(object.list = pcb.list, nfeatures = 3000)

pcb.list <- PrepSCTIntegration(object.list = pcb.list, anchor.features = anchor_features)

pcb.anchors <- FindIntegrationAnchors(object.list = pcb.list, normalization.method = "SCT",
                                      anchor.features = anchor_features)

merged <- IntegrateData(anchorset = pcb.anchors, normalization.method = "SCT")
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
# merged <- FindClusters(merged, resolution = 0.5)
merged <- FindClusters(merged, resolution = 0.1)

SpatialDimPlot(merged, ncol = 4,
               label = TRUE, label.size = 3)&NoLegend()


#===============================================================================
head(merged@meta.data)
merged$celltype <- as.character(merged$seurat_clusters)
merged$celltype[merged$celltype %in% c('1', '4')] <- 'Neocortex'
merged$celltype[merged$celltype %in% c('3', '7')] <- 'Hippocampal region'
merged$celltype[merged$celltype %in% c('2', '6')] <- 'Thalamus'
merged$celltype[merged$celltype %in% c('0')] <- 'Fiber tracts'
merged$celltype[merged$celltype %in% c('5')] <- 'Caudoputamen'
# merged$celltype[merged$celltype %in% c('6')] <- 'Thalamus'

SpatialDimPlot(merged, ncol = 4,
               label = TRUE, label.size = 3, group.by = "celltype")&NoLegend()

polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
colors <- c(`Neocortex`="#1CBE4FFF", 
            `Hippocampal region`="#782AB6FF",
            `Thalamus`="#90AD1CFF", 
            `Fiber tracts`= "#C4451CFF",
            `Caudoputamen`="#1C7F93FF")


pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/1_spatial_dimplot.pdf",
    height = 10, width = 20)
SpatialDimPlot(merged, ncol = 4,
               label = TRUE, label.size = 3, cols = colors, group.by = "celltype",
               image.alpha = 0.5,
               pt.size.factor = 3.5,
               alpha = 0.7)&NoLegend()
dev.off()

SpatialDimPlot(merged, ncol = 4,
               label = TRUE, label.size = 3, group.by = "seurat_clusters",
               image.alpha = 0.5,
               pt.size.factor = 3.5,
               alpha = 0.7)&NoLegend()
merged <- FindClusters(merged, resolution = 0.5)
#===============================================================================
# Calculate DEGs across all anatomical regions
#===============================================================================

head(merged@meta.data)
Idents(merged) <- merged$celltype
merged$celltype.treatment <- paste(Idents(merged), merged$treatment, sep = "_")

# Prior to performing differential expression, we first run PrepSCTFindMarkers, which ensures that the fixed value is set properly. 
# Then we use FindMarkers(assay="SCT") to find differentially expressed genes. 

merged <- PrepSCTFindMarkers(merged)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Save the integrated object
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveRDS(merged, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/integrated_spatial_seurat.rds")
merged <- readRDS(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/integrated_spatial_seurat.rds")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform DEG across celltype
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEG <- list()
Idents(merged) <- merged$celltype
for(i in levels(merged)){
  cluster = i
  message("Calculating DEG for ",cluster)
  Idents(merged) <- merged$celltype.treatment
  cluster.marker <- FindMarkers(merged, ident.1 = paste0(cluster,"_PCB"), 
                                ident.2 = paste0(cluster, "_Control"), 
                                min.pct=0.2, logfc.threshold=0, 
                                test.use = "wilcox",
                                assay = "SCT") 
  DEG[[cluster]] <- cluster.marker
}

saveRDS(DEG, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/DEG_regions_PCB.rds")
DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/DEG_regions_PCB.rds")

write.xlsx(DEG, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/DEG_complete_PCB_all_regions.xlsx", 
           rowNames = TRUE)

# Filter significant DEGs across all anatomical regions
DEG_sig <- list()
for (df in names(DEG)) {
  tryCatch({
    data = DEG[[df]]
    message("Filtering significant genes for ", df)
    gene.list <- data %>%
      dplyr::filter(abs(avg_log2FC) > 0.20 & p_val_adj < 0.05)%>%
      tibble::rownames_to_column(var = "gene")
    DEG_sig[[df]] <- gene.list
  })
}
write.xlsx(DEG_sig, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/DEG_filtered_PCB_all_regions.xlsx", 
           rowNames = FALSE)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Data Visualization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
merged <- readRDS(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/integrated_spatial_seurat.rds")
head(merged@meta.data)
merged@assays
Idents(merged) <- merged$celltype

DefaultAssay(merged) <- "Spatial"

colors <- c(`Neocortex`="#1CBE4FFF", 
            `Hippocampal region`="#782AB6FF",
            `Thalamus`="#90AD1CFF", 
            `Fiber tracts`= "#C4451CFF",
            `Caudoputamen`="#1C7F93FF")

# UMAP plot showing clusters
#===============================================================================
pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/2_UMAP_plot_clusters.pdf",
    height = 8, width = 10)
DimPlot_scCustom(merged, colors_use = colors, figure_plot = TRUE)
dev.off()

# UMAP plot showing different treatment groups
#===============================================================================
pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/3_UMAP_plot_Cntrl_PCB.pdf",
    height = 8, width = 10)
DimPlot_scCustom(merged, group.by = "treatment", colors_use = c("lightblue", "maroon"),
                 figure_plot = TRUE)
dev.off()

# Violin Plot of selected genes
#===============================================================================
pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/4_Violin_plot_selected_genes.pdf",
    height = 8, width = 8)
Stacked_VlnPlot(merged, features = c("Tcf4", "Mef2a", "Spock1", "Dpysl2", "Morf4l1", "Gstp1"),
                split.by = "treatment", x_lab_rotate = 45,
                colors_use = c("lightblue", "maroon"),
                plot_legend = TRUE)
dev.off()


pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/4a_Violin_plot_selected_genes_datapoints.pdf",
    height = 8, width = 8)
Stacked_VlnPlot(merged, features = c("Tcf4", "Mef2a", "Spock1", "Dpysl2", "Morf4l1", "Gstp1"),
                split.by = "treatment", x_lab_rotate = 45,
                colors_use = c("lightblue", "maroon"), pt.size = 0.01,
                plot_legend = TRUE)
dev.off()

# Box Plot of selected genes
#===============================================================================
pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Box_plot_Gapdh.pdf",
    height = 8, width = 8)
do_BoxPlot(merged, feature = "Gapdh",
           split.by = "treatment",
           slot = "counts",
           xlab = "")+
  scale_fill_manual(values = c("Control"="lightblue",
                               "PCB"="maroon"))
dev.off()

p1 <- do_BoxPlot(merged, feature = "Spock1", 
                 split.by = "treatment", 
                 slot = "counts",
                 xlab = "",
                 legend.position = "none")+
  scale_fill_manual(values = c("Control"="lightblue", 
                               "PCB"="maroon"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p2 <- do_BoxPlot(merged, feature = "Dpysl2", 
                 split.by = "treatment", 
                 slot = "counts",
                 xlab = "",
                 legend.position = "none")+
  scale_fill_manual(values = c("Control"="lightblue", 
                               "PCB"="maroon"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p3 <- do_BoxPlot(merged, feature = "Gstp1", 
                 split.by = "treatment", 
                 slot = "counts",
                 xlab = "",
                 legend.position = "none")+
  scale_fill_manual(values = c("Control"="lightblue", 
                               "PCB"="maroon"))+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p4 <- do_BoxPlot(merged, feature = "Morf4l1", 
                 split.by = "treatment", 
                 slot = "counts",
                 xlab = "",
                 legend.position = "right")+
  scale_fill_manual(values = c("Control"="lightblue", 
                               "PCB"="maroon"))

p1/p2/p3/p4

pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Box_plot_stacked.pdf",
    height = 12, width = 10)
p1/p2/p3/p4
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export QC data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
control1_summary <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_visium_v1/Control_1/outs/metrics_summary.csv")
control2_summary <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_visium_v1/Control_2/outs/metrics_summary.csv")
pcb1_summary <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_visium_v1/PCB_1/outs/metrics_summary.csv")
pcb2_summary <- read.csv("/home/bbasu/LSS/lss_schatterj/PCB_visium_v1/PCB_2/outs/metrics_summary.csv")


combined_qc <- rbind(control1_summary, control2_summary,
                     pcb1_summary, pcb2_summary)
combined_qc$`Age (months)` <- 3
combined_qc$Sex <- "Male"
combined_qc$condition <- c("Vehicle", "Vehicle", "HR-PCB", "HR-PCB")
combined_qc$slide_number <- "V13Y08-074"
combined_qc$area_number <- c("A1", "C1", "B1", "D1")
combined_qc$`#replicate_sections_mounted` <- 2

combined_qc <- as.data.frame(t(combined_qc))


write.xlsx(combined_qc, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/QC_data_combined.xlsx", 
           rowNames = TRUE)

# Marker genes for all clusters
#===============================================================================
Idents(merged) <- merged$celltype
cluster.markers <- FindAllMarkers(merged, assay = "SCT", only.pos = TRUE, logfc.threshold = 0.20,
                                  test.use = "wilcox")

# Retrieve the top 5 marker genes per cluster
# Use whichever genes have the highest values under the AVG_LOG column
top5 <- cluster.markers %>% group_by(cluster) %>%
  dplyr::slice_max(get(grep("^avg_log", colnames(cluster.markers), value = TRUE)),
                   n = 5)
# top10 markers
top10 <- cluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Create the dot plot
pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/dot_plot_top5_marker.pdf", 
    height = 5, width = 10)
Seurat::DotPlot(merged, features = unique(top10$gene)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,
                                                     size = 8, hjust = 1)) +
  Seurat::NoLegend()+
  ggpubr::labs_pubr()
dev.off()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Search for marker gene expression
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Neocortex: Cux1, Cux2
# Fiber Tracts: Sox10, Olig2 
# Thalamus: Lef1, Gbx2
# Hippocampal region: Lct, Crlf1
# Caudoputamen: Adora2a, Gpr6

DefaultAssay(merged) <- "Spatial"
Stacked_VlnPlot(merged, features = c("Cux1", "Cux2", 
                                     "Nme9", "F5",
                                     "Lef1", "Gbx2", 
                                     "Lct", "Crlf1",
                                     "Adora2a", "Gpr6"),
                x_lab_rotate = 45,
                colors_use = colors,
                plot_legend = TRUE)

pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/5_Marker_genes_PCB.pdf", 
    height = 6, width =10)
DotPlot_scCustom(merged, features = c("Cux1", "Cux2", 
                                      "Sox10", "Olig2",
                                      "Lef1", "Gbx2", 
                                      "Lct", "Crlf1",
                                      "Adora2a", "Gpr6"))
dev.off()

# Neocortex marker
p1 <- SpatialFeaturePlot(merged, features = "Myl4", pt.size.factor = 3, images = "Control_1.2")
p2 <- SpatialFeaturePlot(merged, features = "Trbc2", pt.size.factor = 3, images = "Control_1.2")
p3 <- SpatialFeaturePlot(merged, features = "Hs3st2", pt.size.factor = 3, images = "Control_1.2")

# Fiber tracts marker
p4 <- SpatialFeaturePlot(merged, features = "Sox10", pt.size.factor = 3, images = "Control_1.2")
p5 <- SpatialFeaturePlot(merged, features = "Olig2", pt.size.factor = 3, images = "Control_1.2")
p6 <- SpatialFeaturePlot(merged, features = "Opalin", pt.size.factor = 3, images = "Control_1.2")

# Hippocampus marker
p7 <- SpatialFeaturePlot(merged, features = "C1ql2", pt.size.factor = 3, images = "Control_1.2")
p8 <- SpatialFeaturePlot(merged, features = "Lct", pt.size.factor = 3, images = "Control_1.2")
p9 <- SpatialFeaturePlot(merged, features = "Npy2r", pt.size.factor = 3, images = "Control_1.2")

# Thalamus marker
p10 <- SpatialFeaturePlot(merged, features = "Smpx", pt.size.factor = 3, images = "Control_1.2")
p11 <- SpatialFeaturePlot(merged, features = "Atp2a1", pt.size.factor = 3, images = "Control_1.2")
p12 <- SpatialFeaturePlot(merged, features = "Lef1", pt.size.factor = 3, images = "Control_1.2")

# Caudoputamen marker
p13 <- SpatialFeaturePlot(merged, features = "Adora2a", pt.size.factor = 3, images = "Control_1.2")
p14 <- SpatialFeaturePlot(merged, features = "Gpr6", pt.size.factor = 3, images = "Control_1.2")
p15 <- SpatialFeaturePlot(merged, features = "Sh3rf2", pt.size.factor = 3, images = "Control_1.2") 

pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Canonical_Marker_genes_PCB_neocortex_fibretract.pdf", 
    height = 8, width =10)
(p1+p2+p3)/(p4+p5+p6)
dev.off()

pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Canonical_Marker_genes_PCB_hippo_thalamus.pdf", 
    height = 8, width =10)
(p7+p8+p9)/(p10+p11+p12)
dev.off()

pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Canonical_Marker_genes_PCB_caudoputamen_thalamus.pdf", 
    height = 8, width =10)
(p10+p11+p12)/(p13+p14+p0)
dev.off()

pdf(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Canonical_Marker_genes_PCB_all.pdf", 
    height = 12, width =10)
FeaturePlot_scCustom(merged, features = c("Myl4", "Trbc2", "Hs3st2",
                                          "Sox10", "Olig2", "Opalin",
                                          "Lef1", "Atp2a1", "Smpx",
                                          "Lct", "C1ql2", "Npy2r",
                                          "Adora2a", "Gpr6"), num_columns = 3)
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Volcano plot for all clusters
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DefaultAssay(merged) <- "SCT"
DEG <- list()
Idents(merged) <- merged$celltype
for(i in levels(merged)){
  cluster = i
  message("Calculating DEG for ",cluster)
  Idents(merged) <- merged$celltype.treatment
  cluster.marker <- FindMarkers(merged, ident.1 = paste0(cluster,"_PCB"), 
                                ident.2 = paste0(cluster, "_Control"), 
                                min.pct=0.2, logfc.threshold=0, 
                                test.use = "wilcox",
                                assay = "SCT") 
  DEG[[cluster]] <- cluster.marker
}


colors <- c("#1CBE4FFF",
            "#C4451CFF",
            "#90AD1CFF", 
            "#782AB6FF",
            "#1C7F93FF")
plots <- list()
for(i in names(DEG)){
  cluster = i
  message("Plotting for ",cluster)
  cluster_index = which(names(DEG) == cluster)
  #Data wrangling
  cluster.marker <- DEG[[cluster]]
  data <- data.frame(gene = row.names(cluster.marker),
                     pval = -log10(cluster.marker$p_val_adj), 
                     lfc = cluster.marker$avg_log2FC)
  
  data <- mutate(data, color = case_when(data$lfc > 0.2 & data$pval > 1.3 ~ "Increased",
                                         data$lfc < -0.2 & data$pval > 1.3 ~ "Decreased",
                                         data$lfc >= -0.2 & data$lfc <= 0.2 & data$pval > 1.3 ~ "nonsignificant",
                                         data$pval < 1.3 ~ "nonsignificant"))
  
  n_deg_up <- data.frame(table(data$color))%>%
    dplyr::filter(Var1 == "Increased")
  n_deg_down <- data.frame(table(data$color))%>%
    dplyr::filter(Var1 == "Decreased")
  
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(data, aes(x = lfc, y = pval, color = color))+
    geom_point(size = 1, na.rm = T) +
    scale_color_manual(name = "Directionality",
                       values = c(Increased = colors[cluster_index], 
                                  Decreased = colors[cluster_index], 
                                  nonsignificant = "gray90")) +
    theme_base() + # change overall theme
    theme(legend.position = "none") + # change the legend
    xlim(-2,2)+
    scale_y_continuous(breaks = seq(0, 120, by=20), limits=c(0,120))+
    ggtitle(cluster)+
    theme(plot.title = element_text(size = 10))+
    annotate("text", x = -1.5, y = 120, label = capture.output(writeLines(paste0("Up :", n_deg_up$Freq))),
             size = 4,
             fontface =1.2)+
    annotate("text", x = -1.5, y = 110, label = capture.output(writeLines(paste0("Down :", n_deg_down$Freq))),
             size = 4,
             fontface =1.2)
  
  if(cluster %in% c("Neocortex")){
    # Add ggplot2 layers
    p1 <- vol+
      xlab("")+
      ylab(expression(-log[10]~"(Adj P Value)"))
  }
  else if(cluster %in% c("Thalamus", "Caudoputamen")){
    # Add ggplot2 layers
    p1 <- vol +
      xlab(expression(log[2]~"(Fold Change)"))+
      ylab("")
  }
  else if(cluster %in% c("Hippocampal region")){
    # Add ggplot2 layers
    p1 <- vol +
      xlab(expression(log[2]~"(Fold Change)"))+
      ylab(expression(-log[10]~"(Adj P Value)"))
  }
  else{
    # Add ggplot2 layers
    p1 <- vol +
      xlab("")+
      ylab("")
    
  }
  plots[[cluster]] <- p1
}


pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/6_volcano_plot_PCB.pdf", height = 8, width = 10)
plots$Neocortex+plots$`Fiber tracts`+plots$Thalamus+plots$`Hippocampal region`+plots$Caudoputamen
dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Cell type proportion
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Idents(merged) <- merged$celltype
head(merged@meta.data)
##Calculate the proportion of each cell type comparing the conditions
# How many cells are in each cluster
table(Idents(merged))

## How many cells are in each replicate?
table(merged$orig.ident)

## What proportion of cells are in each cluster?
prop.table(table(Idents(merged)))

### How does cluster membership vary by condition?
table(Idents(merged), merged$orig.ident)
table(Idents(merged), merged$treatment)


# Proportion animal level
#===============================================================================
prop_df <- merged@meta.data %>%
  group_by(orig.ident, treatment, celltype) %>%
  summarise(n = n()) %>%
  group_by(orig.ident, treatment) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write.xlsx(prop_df, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/CelltypeProportion_data.xlsx", rowNames = FALSE)

wilcox_results <- prop_df %>%
  group_by(celltype) %>%
  summarise(
    p_value = wilcox.test(prop ~ treatment)$p.value,
    mean_prop_ctrl = mean(prop[treatment == "Control"]),
    mean_prop_trt  = mean(prop[treatment == "PCB"]),
    diff = mean_prop_trt - mean_prop_ctrl
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) # Benjamini & Hochberg /FDR
?p.adjust

write.xlsx(wilcox_results, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Wilcox_result_proportion_analysis.xlsx",
           rowNames = FALSE)

###Keep the proportion in a data frame for plotting
df <- as.data.frame(prop.table(table(Idents(merged), 
                                     merged$orig.ident), margin = 2))
df$Freq <- round(df$Freq, 2)

head(df)

library(ggpubr)
colors <- c("#1CBE4FFF",
            "#C4451CFF",
            "#90AD1CFF", 
            "#782AB6FF",
            "#1C7F93FF")

plot <- ggbarplot(df, x = "Var2", y = "Freq", fill = "Var1", 
                  palette = alpha(colors, 0.9), 
                  xlab = "", 
                  ylab = "Proportion", label = T)+
  labs_pubr()+
  theme(legend.position = "right",
        legend.title = element_blank())+
  rotate_x_text(45)


plot
ggsave(filename = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/CelltypeProportion_animal_level.pdf",
       plot,
       height = 6,
       width = 8,
       units = "in",
       dpi = 600)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)

DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/DEG_regions_PCB.rds")
# Filter significant genes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEG_sig <- list()

for(i in names(DEG)){
  cluster = i
  message("Doing analysis for ", cluster)
  cellType.DEGs = DEG[[cluster]]
  gene.list <- cellType.DEGs %>%
    dplyr::filter(abs(avg_log2FC) > 0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")%>%
    dplyr::mutate(Cluster = cluster)
  gene.list$entrez <- mapIds(x = org.Mm.eg.db,
                             keys = gene.list$gene,
                             column = "ENTREZID",
                             keytype = "SYMBOL",
                             multiVals = "first")
  gene.list$group <- "Upregulated"
  gene.list$group[gene.list$avg_log2FC < 0] <- "Downregulated"
  DEG_sig[[cluster]] <- gene.list
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform KEGG enrichment across all data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mydf <- bind_rows(DEG_sig)

enrich_kegg <- compareCluster(data=mydf, 
                              entrez~group+Cluster,  
                              fun="enrichKEGG",
                              organism = 'mmu')
# Convert EntrezID to Gene Symbol
kegg <- setReadable(enrich_kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
Kegg_data <- kegg@compareClusterResult 
# Store the data
write.xlsx(Kegg_data, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Compare_Cluster_KEGG_data.xlsx", rowNames = FALSE)


pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Compare_Cluster_KEGG.pdf", height = 12, width = 10)
dotplot(enrich_kegg, x="group")+facet_grid(~Cluster.1)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Perform GO enrichment across all data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
enrich_GO <- compareCluster(data=mydf, 
                            entrez~group+Cluster,  
                            fun="enrichGO",
                            OrgDb = org.Mm.eg.db)


enrich_GO_data <- enrich_GO@compareClusterResult %>% distinct(ID, Cluster, .keep_all = TRUE)

enrich_GO@compareClusterResult <- enrich_GO_data

enrichGO_simplify <- simplify(enrich_GO)

GO_data <- enrichGO_simplify@compareClusterResult
# Convert EntrezID to Gene Symbol
GO_df <- setReadable(enrichGO_simplify, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
GO_data <- GO_df@compareClusterResult 
# Store the data
write.xlsx(GO_data, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Compare_Cluster_enrichGO_data.xlsx", rowNames = FALSE)


pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Compare_Cluster_enrichGO.pdf", height = 12, width = 10)
dotplot(enrichGO_simplify, x="group")+facet_grid(~Cluster.1)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# genes detected in <10% spots per replicate
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# proportion of cell types in each biological replicate
head(merged@meta.data)
Assays(merged)
Idents(merged) <- merged$class

prop.table(table(Idents(merged), merged$orig.ident), margin = 2)

# Get replicate info
replicates <- unique(merged$orig.ident)

# Initialize a list to store genes detected in <10% cells per replicate
low_detected_genes <- list()

DefaultAssay(merged) <- "RNA"
# rep = "Control_1a"
# 
for (rep in replicates) {
  # Subset Seurat object by replicate
  rep_obj <- subset(merged, subset = orig.ident == rep)
  
  # expression matrix
  mat <- GetAssayData(rep_obj, layer = "counts")
  
  # fraction of cells where gene is detected (count > 0)
  detection_frac <- Matrix::rowSums(mat > 0) / ncol(mat)
  
  # Genes detected in less than 10% of cells
  low_detected_genes[[rep]] <- names(detection_frac[detection_frac < 0.1])
}

pcb_low_detected_genes <- low_detected_genes[c( "PCB_1a","PCB_1b","PCB_2a","PCB_2b")]
control_low_detected_genes <- low_detected_genes[c("Control_1a","Control_1b","Control_2a","Control_2b")]

# Find genes that are lowly detected in ALL replicates
genes_low_all_reps_pcb <- Reduce(intersect, pcb_low_detected_genes)
genes_low_all_reps_control <- Reduce(intersect, control_low_detected_genes)

# check how many genes
length(genes_low_all_reps_pcb)
length(genes_low_all_reps_control)

# Combine the genes
genes_low_all_reps <- union(genes_low_all_reps_control, genes_low_all_reps_pcb)


# check if these genes were present in the significant DEG list
DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/DEG_regions_PCB.rds")

# Store significant DEGs across regions
DEG_sig <- list()
for(i in names(DEG)){
  cluster = i
  message("Doing analysis for ", cluster)
  cellType.DEGs = DEG[[cluster]]
  gene.list <- cellType.DEGs %>%
    dplyr::filter(abs(avg_log2FC) > 0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")
  DEG_sig[[cluster]] <- gene.list
}


# lowly-expressed DEGs per cluster
low_expr_deg_by_cluster <- list()

# Loop over clusters
for (clust in names(DEG_sig)) {
  
  low_expr_genes <- genes_low_all_reps
  deg_genes <- DEG_sig[[clust]]
  
  # Find overlap
  matched_genes <- intersect(low_expr_genes, deg_genes$gene)
  
  # Store
  low_expr_deg_by_cluster[[clust]] <- matched_genes
  
}

low_expr_deg_by_cluster$Neocortex


# No gene was found in the DEG list of major brain regions which was detected in less than 10% of spots of each replicate !!

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export normalized gene expression data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

merged <- readRDS(file = "/home/bbasu/LSS/lss_schatterj/PCB_data/integrated_spatial_seurat.rds")
Idents(merged) <- merged$celltype

DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/DEG_regions_PCB.rds")

# Filter significant DEGs across all anatomical regions
DEG_sig <- list()
for (df in names(DEG)) {
  tryCatch({
    data = DEG[[df]]
    message("Filtering significant genes for ", df)
    gene.list <- data %>%
      dplyr::filter(abs(avg_log2FC) > 0.20 & p_val_adj < 0.05)%>%
      tibble::rownames_to_column(var = "gene")
    DEG_sig[[df]] <- gene.list
  })
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Extract pseudobulk normalized expression
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_path <- "/home/bbasu/LSS/lss_schatterj/PCB_data/Pseudobulk_data/"

for (df in levels(merged)){
  cluster = df
  message("Pseudobulking for ",cluster)
  subset_data <- subset(merged, idents = cluster)
  Sig_genes <- DEG_sig[[cluster]]$gene
  
  pseudobulk <- AggregateExpression(subset_data, assays = "SCT", return.seurat = T,
                                    group.by = c("orig.ident", "treatment"))
  
  pseudobulk.normalized <- as.data.frame(pseudobulk[["SCT"]]$data)
  
  filename.complete <- paste0(output_path, cluster, "_pseudobulk_completegeneList.xlsx")
  #Save this normalized data for all genes
  write.xlsx(pseudobulk.normalized, file = filename.complete, 
             rowNames = T)
  
  # Filter the normalized expression for DEGs of that cluster
  pseudobulk.normalized.DEGs <- pseudobulk.normalized %>%
    dplyr::filter(rownames(.) %in% Sig_genes)
  
  filename.DEGs <- paste0(output_path, cluster, "_pseudobulk_DEGsOnly.xlsx")
  #Save this normalized data for DEGs only in that cluster
  write.xlsx(pseudobulk.normalized.DEGs, file = filename.DEGs, 
             rowNames = T)
}
