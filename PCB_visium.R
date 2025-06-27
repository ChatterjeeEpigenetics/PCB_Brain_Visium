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


###Proportion of each cluster by condition
prop.table(table(Idents(merged), merged$treatment), margin = 2)

##Fisher's test to evaluate the power of cell composition analysis
library(rstatix)
test <- rstatix::row_wise_fisher_test(as.matrix(table(Idents(merged), 
                                                      merged$treatment)),
                                      p.adjust.method = "BH")

test
write.csv(test, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/Fisher_test_for_cell_composition_analysis.csv",
          row.names = F)
###Keep the proportion in a data frame for plotting
df <- as.data.frame(prop.table(table(Idents(merged), 
                                     merged$treatment), margin = 2))
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
ggsave(filename = "/home/bbasu/LSS/lss_schatterj/PCB_data/plots/7_CelltypeProportion.pdf",
       plot,
       height = 6,
       width = 6,
       units = "in",
       dpi = 600)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Make a GO enrichment heatmap
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEG <- readRDS("/home/bbasu/LSS/lss_schatterj/PCB_data/DEG_regions_PCB.rds")

# Filter significant upregulated genes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DEG_sig_up <- list()
for(i in names(DEG)){
  cluster = i
  message("Doing analysis for ", cluster)
  cellType.DEGs = DEG[[cluster]]
  gene.list <- cellType.DEGs %>%
    dplyr::filter(avg_log2FC > 0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")%>%
    dplyr::mutate(Cluster = cluster)
  DEG_sig_up[[cluster]] <- gene.list
}

# Make a heatmap of the upregulated genes across the clusters
#===============================================================================
# Extract logFC matrix (per cluster)
expr_mat <- DEG_sig_up %>%
  bind_rows()%>%
  dplyr::select(gene, avg_log2FC, Cluster)%>%
  spread(key = "gene", value = "avg_log2FC")%>%
  tibble::column_to_rownames(var = "Cluster")%>%
  t()

# Replace NA value with 0
expr_mat[is.na(expr_mat)] <- 0

# Replace all positive values with 1 to show one color with red
expr_mat[expr_mat > 0] <- 1

expr_mat <- as.data.frame(expr_mat)
mat <- expr_mat
colnames(mat)
#===============================================================================
# Order the column names of the matrix
custom_order <- c("Neocortex", "Fiber tracts", "Thalamus", "Hippocampal region", "Caudoputamen")
mat <- mat[, c("Neocortex", "Fiber tracts", "Thalamus", "Hippocampal region", "Caudoputamen")]

#===============================================================================
# Reorder the rows to first show unique genes, then common genes
common_genes <- row.names(mat[rowSums(mat)> 1, ])
unique_genes <- row.names(mat[rowSums(mat)== 1, ])

# filter the matrix that shows unique genes
unique_gene_mat <- mat[unique_genes, ]

# Now filter unique genes from each cluster
neocortex.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Neocortex"] == 1, ])
fiber.tracts.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Fiber tracts"] == 1, ])
thalamus.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Thalamus"] == 1, ])
hippocampus.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Hippocampal region"] == 1, ])
caudoputamen.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Caudoputamen"] == 1, ])

# Gene order [First the unique genes from each cluster and then the common genes]  
ordered_genes <- c(neocortex.genes, fiber.tracts.genes,
                   thalamus.genes, hippocampus.genes,
                   caudoputamen.genes, common_genes)

#===============================================================================
# Split the heatmap at desired rows (at each cluster's unique genes) so that it would be easy to show GO terms
length(neocortex.genes) #22
length(fiber.tracts.genes) #66
length(thalamus.genes) #49
length(hippocampus.genes) #67
length(caudoputamen.genes) #15
length(common_genes) # 64


row_split <- rep(c("Neocortex", "Fiber tracts", 
                   "Thalamus", "Hippocampal region", 
                   "Caudoputamen", "common"), 
                 times = c(22, 66, 49, 67, 15, 64))

# Turn into a factor with custom order: 
row_split <- factor(row_split, levels = c("Neocortex", "Fiber tracts", "Thalamus", 
                                          "Hippocampal region", "Caudoputamen", "common"))

length(row_split)
nrow(mat)
length(unique_genes)
length(ordered_genes)
length(common_genes)
#===============================================================================
# Annotate the GO terms associated with unique and common genes

# Perform GO enrichment on the unique upregulated genes
gene_list <- list("Neocortex" = neocortex.genes, 
                  "Fiber tracts" = fiber.tracts.genes,
                  "Thalamus" = thalamus.genes, 
                  "Hippocampal region" = hippocampus.genes,
                  "Caudoputamen" = caudoputamen.genes,
                  "common.genes" = common_genes)

# Trycatch was used to ignore errors associated with any cluster
go_results <- lapply(gene_list, function(genes){
  tryCatch(ego <- enrichGO(gene          = genes,
                           OrgDb         = org.Mm.eg.db, 
                           ont           = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           keyType = "SYMBOL",
                           readable      = TRUE), error=function(e) NULL)
  tryCatch(ego@result, error=function(e) NULL)
})

# Check the pathways and the gene names (genes to annotate!)
neocortex_GO <- go_results$Neocortex
fiber_tracts_GO <- go_results$`Fiber tracts`
thalamus_GO <- go_results$Thalamus
hippo_GO <- go_results$`Hippocampal region`
caudoputamen_GO <- go_results$Caudoputamen
common_GO <- go_results$common.genes


# GO_terms = lapply(go_results, function(df) {
#   terms <- df$Description %>% head(3)
# })

# Save the MF enrichment result
write.xlsx(go_results, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/Upregulated_GO_MF_all_regions_unique_common.xlsx", 
           rowNames = F)

#===============================================================================
# colors_celltype <- paletteer::paletteer_d("ggthemes::Tableau_10")

# Annotate celltypes
ha = HeatmapAnnotation(df = data.frame(Celltypes = c("Neocortex", "Fiber tracts", "Thalamus", 
                                                     "Hippocampal region", "Caudoputamen")),
                       col = list(Celltypes = c("Neocortex"="#1CBE4FFF",
                                                "Fiber tracts"= "#C4451CFF",
                                                "Thalamus"="#90AD1CFF",
                                                "Hippocampal region"="#782AB6FF",
                                                "Caudoputamen"="#1C7F93FF")),
                       annotation_legend_param = list(title = "Celltype",
                                                      title_gp = gpar(fontsize = 10,
                                                                      fontface = "bold")),
                       annotation_name_side = "left",
                       simple_anno_size = unit(2, "mm"))

#===============================================================================
# Annotate gene names
mat <- mat[ordered_genes, custom_order]

ann <- mat %>%
  dplyr::filter(rownames(mat) %in% c("Gstm1", "Gstm5",
                                     "Arpp19", "Cnih3", "Fkbp1a", "Fxyd7", "Kcnip3", "Kcnmb4", "Ywhah",
                                     "Fau", "Mrpl30", "Rpl11", "Rpl18", "Rpl18a", "Rpl24", "Rpl27a", "Rps15", "Rps20",
                                     "Rpl10", "Rpl10a", "Rpl12", "Rpl13", "Rpl13a", "Rpl14", "Rpl15", "Rpl17", "Rpl19", "Rps7", "Rps14"))

vrn <- rownames(mat) %in% rownames(ann)

row_anno = rowAnnotation(link = anno_mark(at = which(vrn),
                                          side = "left",
                                          labels = row.names(mat)[vrn],
                                          labels_gp = gpar(fontsize = 12),
                                          padding = unit(1, "mm")))

#===============================================================================
# Right annotation with celltype

right_anno <- rowAnnotation(
  Genes = rep(c("Neocortex", "Fiber tracts", 
                "Thalamus", "Hippocampal region", 
                "Caudoputamen", "common"), 
              times = c(22, 66, 49, 67, 15, 64)),
  col = list(Genes = c("Neocortex"="#1CBE4FFF",
                       "Fiber tracts"= "#C4451CFF",
                       "Thalamus"="#90AD1CFF",
                       "Hippocampal region"="#782AB6FF",
                       "Caudoputamen"="#1C7F93FF",
                       "common"= "#B2A84B"))
)

length(row_split)
nrow(mat)
#===============================================================================
# Make the matrix
mat <- mat[ordered_genes, custom_order]
mat <- as.matrix(mat)
#===============================================================================
# Make the heatmap
ht1 <- Heatmap(mat,cluster_columns = F, cluster_rows = F,
               row_split = row_split,
               width = unit(4, "cm"), 
               height = unit(25, "cm"),
               column_order = custom_order,
               show_column_names = F,
               show_row_names = F,
               col = c("#a9dfde", "red"),
               bottom_annotation = ha,
               left_annotation = row_anno,
               right_annotation = right_anno,
               name = "DEG",
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10,
                                 fontface = "bold")
               ))

ht1

pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Upregulated_genes_GO_enrichment_heatmap.pdf",
    height = 12, width = 8)
ht1
dev.off()
#===============================================================================

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Go enrichment heatmap with significant downregulated genes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Filter significant downregulated genes
DEG_sig_down <- list()
for(i in names(DEG)){
  cluster = i
  message("Doing analysis for ", cluster)
  cellType.DEGs = DEG[[cluster]]
  gene.list <- cellType.DEGs %>%
    dplyr::filter(avg_log2FC < -0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")%>%
    dplyr::mutate(Cluster = cluster)
  DEG_sig_down[[cluster]] <- gene.list
}

# Make a heatmap of the downregulated genes across the clusters
#===============================================================================
# Extract logFC matrix (per cluster)
expr_mat <- DEG_sig_down %>%
  bind_rows()%>%
  dplyr::select(gene, avg_log2FC, Cluster)%>%
  spread(key = "gene", value = "avg_log2FC")%>%
  tibble::column_to_rownames(var = "Cluster")%>%
  t()

# Replace NA value with 0
expr_mat[is.na(expr_mat)] <- 0

# Replace all negative values with 1 to show one color with red
expr_mat[expr_mat < 0] <- 1

expr_mat <- as.data.frame(expr_mat)
mat <- expr_mat
colnames(mat)
#===============================================================================
# Order the column names of the matrix
custom_order <- c("Neocortex", "Fiber tracts", "Thalamus", "Hippocampal region", "Caudoputamen")
mat <- mat[, custom_order]

#===============================================================================
# Reorder the rows to first show unique genes, then common genes
common_genes <- row.names(mat[rowSums(mat)> 1, ])
unique_genes <- row.names(mat[rowSums(mat)== 1, ])

# filter the matrix that shows unique genes
unique_gene_mat <- mat[unique_genes, ]

# Now filter unique genes from each cluster
neocortex.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Neocortex"] == 1, ])
fiber.tracts.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Fiber tracts"] == 1, ])
thalamus.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Thalamus"] == 1, ])
hippocampus.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Hippocampal region"] == 1, ])
caudoputamen.genes <- rownames(unique_gene_mat[unique_gene_mat[, "Caudoputamen"] == 1, ])

# Gene order [First the unique genes from each cluster and then the common genes]  
ordered_genes <- c(neocortex.genes, fiber.tracts.genes,
                   thalamus.genes, hippocampus.genes,
                   caudoputamen.genes, common_genes)

#===============================================================================
# Split the heatmap at desired rows (at each cluster's unique genes) so that it would be easy to show GO terms
length(neocortex.genes) #70
length(fiber.tracts.genes) #46
length(thalamus.genes) #74
length(hippocampus.genes) #25
length(caudoputamen.genes) #18
length(common_genes) #49


row_split <- rep(c("Neocortex", "Fiber tracts", 
                   "Thalamus", "Hippocampal region", 
                   "Caudoputamen", "common"), 
                 times = c(70, 46, 74, 25, 18, 49))

# Turn into a factor with custom order: 
row_split <- factor(row_split, levels = c("Neocortex", "Fiber tracts", "Thalamus", 
                                          "Hippocampal region", "Caudoputamen", "common"))


#===============================================================================
# Annotate the GO terms associated with unique and common genes

# Perform GO enrichment on the unique upregulated genes
gene_list <- list("Neocortex" = neocortex.genes, 
                  "Fiber tracts" = fiber.tracts.genes,
                  "Thalamus" = thalamus.genes, 
                  "Hippocampal region" = hippocampus.genes,
                  "Caudoputamen" = caudoputamen.genes,
                  "common.genes" = common_genes)

# Trycatch was used to ignore errors associated with any cluster
go_results <- lapply(gene_list, function(genes){
  tryCatch(ego <- enrichGO(gene          = genes,
                           OrgDb         = org.Mm.eg.db,
                           ont           = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05,
                           keyType = "SYMBOL",
                           readable      = TRUE), error=function(e) NULL)
  tryCatch(ego@result, error=function(e) NULL)
})

# Check the pathways and the gene names (genes to annotate!)
neocortex_GO <- go_results$Neocortex
fiber_tracts_GO <- go_results$`Fiber tracts`
thalamus_GO <- go_results$Thalamus
hippo_GO <- go_results$`Hippocampal region`
caudoputamen_GO <- go_results$Caudoputamen
common_GO <- go_results$common.genes


# GO_terms = lapply(go_results, function(df) {
#   terms <- df$Description %>% head(3)
# })

# Save the MF enrichment result
write.xlsx(go_results, file = "/home/bbasu/LSS/lss_schatterj/PCB_data/Downregulated_GO_MF_all_regions_unique_common.xlsx", 
           rowNames = F)

#===============================================================================
# colors_celltype <- paletteer::paletteer_d("ggthemes::Tableau_10")

# Annotate celltypes
ha = HeatmapAnnotation(df = data.frame(Celltypes = c("Neocortex", "Fiber tracts", "Thalamus", 
                                                     "Hippocampal region", "Caudoputamen")),
                       col = list(Celltypes = c("Neocortex"="#1CBE4FFF",
                                                "Fiber tracts"= "#C4451CFF",
                                                "Thalamus"="#90AD1CFF",
                                                "Hippocampal region"="#782AB6FF",
                                                "Caudoputamen"="#1C7F93FF")),
                       annotation_legend_param = list(title = "Celltype",
                                                      title_gp = gpar(fontsize = 10,
                                                                      fontface = "bold")),
                       annotation_name_side = "left",
                       simple_anno_size = unit(2, "mm"))

#===============================================================================
# Annotate gene names
mat <- mat[ordered_genes, custom_order]

ann <- mat %>%
  dplyr::filter(rownames(mat) %in% c("Add2","Camk2d","Ewsr1","Map2", "Marcks", "Syt7", "Syt17", "Syt4",
                                     "Chrna4","Ghitm","Hcn2","Kcna2","Kcnab2","Kcnc1","Scn1a","Tmbim6","Ywhag","Ywhaz","Mef2a","Syt1","Tcf4","Ywhah",
                                     "Pde10a","Pde1b",
                                     "Mal","Mbp","Mobp","Plp1", "Ap2a2","App","Cltc","Dnm1","Actb","Ank2","Camk2b","Tuba1a","Tubb2a"))

vrn <- rownames(mat) %in% rownames(ann)

row_anno = rowAnnotation(link = anno_mark(at = which(vrn),
                                          side = "left",
                                          labels = row.names(mat)[vrn],
                                          labels_gp = gpar(fontsize = 12),
                                          padding = unit(1, "mm")))

#===============================================================================
# Right annotation with celltype

right_anno <- rowAnnotation(
  Genes = rep(c("Neocortex", "Fiber tracts", 
                "Thalamus", "Hippocampal region", 
                "Caudoputamen", "common"), 
              times = c(70, 46, 74, 25, 18, 49)),
  col = list(Genes = c("Neocortex"="#1CBE4FFF",
                       "Fiber tracts"= "#C4451CFF",
                       "Thalamus"="#90AD1CFF",
                       "Hippocampal region"="#782AB6FF",
                       "Caudoputamen"="#1C7F93FF",
                       "common"= "#B2A84B"))
)

length(row_split)
nrow(mat)
#===============================================================================
# Make the matrix
mat <- mat[ordered_genes, custom_order]
mat <- as.matrix(mat)
#===============================================================================
# Make the heatmap
ht1 <- Heatmap(mat,cluster_columns = F, cluster_rows = F,
               row_split = row_split,
               width = unit(4, "cm"), 
               height = unit(25, "cm"),
               column_order = custom_order,
               show_column_names = F,
               show_row_names = F,
               col = c("#a9dfde", "#274580"),
               bottom_annotation = ha,
               left_annotation = row_anno,
               right_annotation = right_anno,
               name = "DEG",
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10,
                                 fontface = "bold")
               ))

ht1

pdf("/home/bbasu/LSS/lss_schatterj/PCB_data/plots/Downregulated_genes_GO_enrichment_heatmap.pdf",
    height = 12, width = 8)
ht1
dev.off()
#===============================================================================

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
