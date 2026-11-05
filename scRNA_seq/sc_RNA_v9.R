# load libraries
library(htmltools)
library("sctransform");
library("dplyr");
library("RColorBrewer");
library("ggthemes");
library("cowplot"); 
library("data.table");
library(gridExtra)
library(readxl)
library(celldex)
library(pheatmap)
library(glmGamPoi)
library(gplots)
library('hdf5r')
library('Seurat')
library('Signac')
library(EnsDb.Mmusculus.v79)
library(stringr)
library(patchwork)
library(ggplot2)
library(scales)
library(scCustomize)
library(tidyverse)
library(harmony)
library(future)
library(promises)
library(openxlsx)
library(NMF)
library(CellChat)
library(ggalluvial)
library(ComplexHeatmap)
library(scDblFinder)
library(circlize)
library(biomaRt)

sample_identity <- function(seu_obj,n_color){
  
  seu_obj <- AddMetaData(seu_obj, metadata = seu_obj@active.ident, col.name = 'Assign_name')
  
  cluster_identity = paste(seu_obj$Assign_name,seu_obj@meta.data$orig.ident,sep="_")
  table(cluster_identity)
  
  cluster_identity = cbind(seu_obj$Assign_name,seu_obj@meta.data$orig.ident)
  
  frequency_table = matrix(0,length(unique(seu_obj$orig.ident)),length(levels(seu_obj$Assign_name)))
  rownames(frequency_table) = c("W1","W2","KO1","KO2")
  
  
  for ( i in 1:length(levels(seu_obj$Assign_name))){
    subset = subset(cluster_identity,cluster_identity[,1] ==i)
    frequency = table(subset[,2])
    frequency_table[names(frequency),(i)] = frequency
  }
  
  for ( i in 1:nrow(frequency_table)){
    frequency_table[i,] = frequency_table[i,] / sum(frequency_table[i,])
  }
  
  colnames(frequency_table) = levels(seu_obj$Assign_name)
  WT_m <- apply(frequency_table[1:2,], MARGIN = 2, FUN = mean)
  WT_stdev <- apply(frequency_table[1:2,], MARGIN = 2, FUN = sd)
  
  KO_m <- apply(frequency_table[3:4,], MARGIN = 2, FUN = mean)
  KO_stdev <- apply(frequency_table[3:4,], MARGIN = 2, FUN = sd)
  for (i in 1:dim(frequency_table)[2]){
    t.res = t.test(frequency_table[1:2,i], frequency_table[3:4,i])
    print(colnames(frequency_table)[i])
    print(t.res)
  }
  frequency_table_mean <- rbind(WT_m,KO_m)
  frequency_table_sd <- rbind(WT_stdev,KO_stdev)
  
  frequency_table_mean <- frequency_table_mean %>% as.data.frame %>% mutate(Samples = c("WT","KO")) %>% pivot_longer(cols = -Samples) 
  frequency_table_sd <- frequency_table_sd %>% as.data.frame %>% mutate(Samples = c("WT","KO")) %>% pivot_longer(cols = -Samples) 
  frequency_table <- cbind(frequency_table_mean,frequency_table_sd[,3])
  colnames(frequency_table) = c("Samples","Cell type","Mean","Stdev")
  
  frequency_table$`Cell type` <- factor(frequency_table$`Cell type`,levels = levels(seu_obj$Assign_name))
  
  
  frequency_table %>% mutate(Samples = factor(Samples, levels = c("WT", "KO"))) %>%
    group_by(Samples) %>% mutate(SDPos = 1 - cumsum(Mean) + Mean) %>% 
    ggplot(aes(x=Samples,y = Mean, fill = `Cell type`)) + geom_bar(stat = "identity",color = "black") + theme_classic() +
    labs(title = "Percent(%)", x = "", y = "") +
    geom_errorbar(aes(ymin = SDPos-Stdev, ymax = SDPos+Stdev), width = 0.1, position = "identity") +
    ylim(0, 1.2) + scale_fill_manual(values = n_color)
  
}

mod_SV <- function(SV){
  
  for (ct in seq_along(SV$patches$plots)) {
    max_y <- max(SV[[ct]]$data[1])
    max_y_rounded <- ceiling(max_y * 10) / 10  # Round up to 1 decimal places
    
    SV[[ct]] <- SV[[ct]] + 
      scale_y_continuous(breaks = c(max_y_rounded),
                         labels = scientific_format()) +
      theme(axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_text(),
            axis.title.y = element_text(angle = 0, hjust = 0, vjust = 0.5)
      )
  }
  return(SV)
}

doAddMeta <- function(obj) {
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA)/log10(obj$nCount_RNA)
  obj$mitoRatio <- PercentageFeatureSet(object=obj, pattern = "^mt-")
  obj$mitoRatio <- obj@meta.data$mitoRatio/100
  keepfeat <- grep('^mt-', dimnames(obj)[[1]], value=TRUE, invert=TRUE)
  obj <- subset(obj, features=keepfeat)
  
  return(obj)
}

setwd("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scRNA-seq\\Dec1KO\\RNA-seq_10X\\")

load(".\\Input\\ndW12K12.rda")

ndK1 = UpdateSeuratObject(ndK1)
ndK2 = UpdateSeuratObject(ndK2)
ndW1 = UpdateSeuratObject(ndW1)
ndW2 = UpdateSeuratObject(ndW2)

ndW1$orig.ident = "W1"
ndW2$orig.ident = "W2"
ndK1$orig.ident = "KO1"
ndK2$orig.ident = "KO2"


### Mai QC ###
scrna <- merge(ndK1, y = c(ndK2, ndW1,ndW2))
scrna[["CellIdent"]] = rownames(scrna@meta.data)

## Convert object into singlecellexperiment
scrna.sce <- as.SingleCellExperiment(scrna)
scrna.sce <- scDblFinder(scrna.sce, samples="orig.ident", clusters=FALSE)
sing = rownames(scrna.sce@colData)[which(scrna.sce$scDblFinder.class == "singlet")]
scrna <- subset(scrna,subset = CellIdent %in% sing)
scrna.lst <- SplitObject(scrna, split.by = "orig.ident")
scrna.lst = lapply(scrna.lst, doAddMeta) 

scrna <- merge(scrna.lst[[1]], y = c(scrna.lst[[2]], scrna.lst[[3]],scrna.lst[[4]]))

#save(scrna,file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\RNA_NoDoublets_20240229.rda")

load("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\RNA_NoDoublets_20240229.rda")



data.frame(scrna@meta.data) %>%
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 3000)

data.frame(scrna@meta.data) %>%
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 826)

data.frame(scrna@meta.data) %>%
  ggplot(aes(color=orig.ident, x=log10GenesPerUMI, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.7)

data.frame(scrna@meta.data) %>%
  ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.1)

ml = 826
scrna <- subset(x = scrna, subset = nFeature_RNA > ml)

### 20240305 ###
scrna[["RNA"]] <- as(object = scrna[["RNA"]], Class = "Assay5")
scrna[["RNA"]] <- split(scrna[["RNA"]], f = scrna$orig.ident)
DefaultAssay(scrna) <- 'RNA'
### Log normalization ###

scrna <- NormalizeData(scrna)
scrna <- FindVariableFeatures(scrna)
scrna <- ScaleData(scrna)
scrna <- RunPCA(scrna,seed.use = 465)

### CCA Integration ###
scrna <- IntegrateLayers(object = scrna, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)

scrna[["RNA"]] <- JoinLayers(scrna[["RNA"]])

scrna <- FindNeighbors(scrna, reduction = "integrated.cca", dims = 1:30)
scrna <- FindClusters(scrna, resolution = 0.5, cluster.name = "integrated.cca_clusters",random.seed = 201)
scrna <- RunUMAP(scrna, reduction = "integrated.cca", dims = 1:30,seed.use = 465)

DimPlot(scrna,reduction = "umap",split.by = "orig.ident")
DimPlot(scrna,reduction = "umap")


### 20240418 ###



gene_list_plot <- c(
  "Ins1","Ins2","Iapp","Ucn3", #Beta
  "Pdx1","G6pc2","Nkx6-1",
  "Gcg","Irx2", "Ttr","Arx","Slc7a2",#alpha
  "Sst","Hhex", "Rbp4",#delta
  "Ppy", #PP cell/Gamma 
  "Myl9","Pdgfrb","Rgs5","Col3a1","Acta2", #Stellate
  "Pecam1","Kdr", #Endothelial
  "Krt19", "Pnlip", #Duct
  "Cd86", # Macrophage
  "Cd19", # B cells
  "Cd4","Cd5", # T cells
  "Cd74", 
  "Mki67" # Proliferating cells
) 


Stacked_VlnPlot(seurat_object = scrna, features = gene_list_plot, x_lab_rotate = TRUE)

DimPlot(scrna,label = T,split.by = "orig.ident") + NoLegend()
DimPlot(scrna,label = T) + NoLegend()
DimPlot_scCustom(scrna)

SV_RNA = Stacked_VlnPlot(seurat_object = scrna, features = gene_list_plot, x_lab_rotate = TRUE)
mod_SV(SV_RNA)


### remove doublets LogNormal Seurat Integ###
scrna <- RenameIdents(object = scrna,
                      '2' = "Beta 1",
                      '4' = "Beta 1",
                      '0' = "Beta 2",
                      '1' = "Beta 2",
                      '18' = "Ductual cells",
                      '3' = "Alpha 1",
                      '7' = "Alpha 2",
                      "16" = "Gamma",
                      "5" = "Delta 1",
                      "9" = "Delta 2",
                      "8" = "Immune B cells",
                      "10" = "Immune T cells",
                      "12" = "Macrophages",
                      "14" = "Macrophages",
                      "15" = "Other immune cells",
                      "13" = "Proliferating cells",
                      "6" = "Vascular cells",
                      "11" = "Stellate 1",
                      "17" = "Stellate 2"
)
levels(scrna) = c("Beta 1","Beta 2",
                  "Alpha 1","Alpha 2",
                  "Delta 1","Delta 2",
                  "Gamma","Stellate 1","Stellate 2","Vascular cells","Ductual cells",
                  "Macrophages","Immune B cells","Immune T cells","Other immune cells",
                  "Proliferating cells")
DimPlot_scCustom(scrna)
Stacked_VlnPlot(seurat_object = scrna, features = gene_list_plot, x_lab_rotate = TRUE)
scrna$type = scrna@active.ident


levels(scrna) = rev(levels(scrna))
DotPlot(scrna,features = gene_list_plot) + 
  theme_minimal() +  # Change the overall theme to a minimal style
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels vertically
    axis.title.y = element_blank(),  # Remove y-axis title
    strip.background = element_blank(),  # Remove background of facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  scale_color_gradient(low = "yellow", high = "purple3")

avgexp = AverageExpression(scrna,return.seurat = T)
avgexp = data.frame(GetAssayData(avgexp,slot = "scale.data"))

avgexp = avgexp[rownames(avgexp) %in% gene_list_plot,]
colnames(avgexp) = sub("^RNA\\.", "", colnames(avgexp))

colnames(avgexp) =  gsub("\\.", " ", colnames(avgexp))
avgexp <- rownames_to_column(avgexp, var = "Gene")
# Reshape data into long format
avgexp_long <- gather(avgexp, key = "CellType", value = "AverageExpression",-Gene)


# Define order for x and y axes
x_order <- gene_list_plot
y_order <- levels(scrna)

# Convert variables to factors with the desired order
avgexp_long$CellType <- factor(avgexp_long$CellType, levels = y_order)
avgexp_long$Gene <- factor(avgexp_long$Gene, levels = x_order)


# Create heatmap using ggplot2
ggplot(avgexp_long, aes(x = Gene, y = CellType, fill = AverageExpression)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "purple3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


levels(scrna) = rev(levels(scrna))

save(scrna,file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-RNA\\Outputs\\V7\\Seurat\\Log_integrated_RNA.rda")
load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-RNA\\Outputs\\V7\\Seurat\\Log_integrated_RNA.rda")


#### 0229 ####

merge_cond = 4
DefaultAssay(scrna) <- 'RNA'
cts <- AggregateExpression(scrna,
                           group.by = c("orig.ident","type"),
                           slot = "counts",
                           return.seurat = F)
cts <- cts$RNA
cts.t = t(cts)
cts.t <- as.data.frame(cts.t)
splitRows <- gsub(".*_","",rownames(cts.t))
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))
cts.split.mod <- lapply(cts.split,function(x){
  rownames(x) <- gsub('(.*)_.*','\\1',rownames(x))
  x <- t(x)
})

###############################################
### Pesudo Bulk RNA seq analysis with edgeR ###
###############################################
library(edgeR)
Mth = "TMM"
# remove 0 counts
cts.split.mod <- lapply(cts.split.mod, function(x) {
  x <- x[rowSums(x > 0) == merge_cond, ]
  return(x)
})

#sink("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\tmp.txt")

filtered_files = list()
desq_tot = list()
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

for (i in names(cts.split.mod)){
  if (!is.null(dim(cts.split.mod[[i]])[2])){
    if (dim(cts.split.mod[[i]])[2] == merge_cond){
      counts_type <- cts.split.mod[[i]]
      sampleType<- rep("WT", ncol(counts_type))
      sampleType[grep("K", colnames(counts_type))] <- "KO"
      sampleType = factor(sampleType,levels = c("WT","KO"))


      dgList <- DGEList(counts=counts_type,group = sampleType)
      keep <- filterByExpr(y = dgList)
      dgList <- dgList[keep, , keep.lib.sizes=FALSE]
      dgList <- normLibSizes(dgList,method = Mth)
      design <- model.matrix(~sampleType)
      dgList <- estimateDisp(y = dgList,design)
      ## w/ rep ##
      fit <- glmQLFit(dgList,design)
      qlf <- glmQLFTest(fit)
      tT = as.data.frame(topTags(qlf,n = "Inf"))


      print(i)
      print(tT[which(rownames(tT) == "Bhlhe40"),])
      print(tT[which(rownames(tT) == "Bhlhe41"),])
      print(tT[which(rownames(tT) == "Glul"),])
      print(tT[which(rownames(tT) == "Zfp36l1"),])
      print(tT[which(rownames(tT) == "Stc2"),])
      print(tT[which(rownames(tT) == "Cyp26b1"),])
      print(tT[which(rownames(tT) == "Hip1r"),])
      print(tT[which(rownames(tT) == "Abr"),])
      print(tT[which(rownames(tT) == "Oaz1"),])
      print(tT[which(rownames(tT) == "Mafa"),])
      print(tT[which(rownames(tT) == "Hif1a"),])
      print(tT[which(rownames(tT) == "Gnao1"),])
      print(tT[which(rownames(tT) == "Washc1"),])
      print(tT[which(rownames(tT) == "Scrn1"),])
      print(tT[which(rownames(tT) == "Socs2"),])


      print(dgList$counts[which(rownames(dgList$counts) == "Bhlhe40"),])
      print(dgList$counts[which(rownames(dgList$counts) == "Bhlhe41"),])


      edgeR_result <- (tT)

      #filtered_data <- edgeR_result %>% dplyr::filter(PValue < 0.05)
      filtered_data = edgeR_result
      filtered_files[[i]] <- dim(filtered_data)[1]
      
      

      
      desq_tot[[i]] <- filtered_data
     write.xlsx(filtered_data,paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\",i,"_",Mth,"_nofilter.xlsx"),
              rowNames = T)
    }
  }
}
#sink()

edgeR_count = data.frame(filtered_files)
edgeR_count = t(edgeR_count)
colnames(edgeR_count) = "Number of DEGs"
edgeR_count = data.frame(edgeR_count)
edgeR_count <- rownames_to_column(edgeR_count, var = "Cell.Types")
edgeR_count <- edgeR_count %>% arrange(desc(Number.of.DEGs))

# write.xlsx(edgeR_count,paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\","edgeR_DEG_",Mth,".xlsx"),
#       rowNames = F)

ggplot(edgeR_count, aes(x = reorder(Cell.Types, -Number.of.DEGs), y = Number.of.DEGs)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Number of DEGs per Cell Type", x = "Cell Type", y = "Number of DEGs") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

for (n in names(desq_tot)) {
  
filtered_data = read.xlsx(paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\",n,"_",Mth,"_nofilter.xlsx")
                          ,rowNames = T)

gene_info <- getBM(
  attributes = c("external_gene_name", "ensembl_gene_id"), # Specify the desired attributes
  filters = "external_gene_name",                         # Filter by gene name
  values = rownames(filtered_data),                              # Provide the list of gene names
  mart = ensembl                                           # Use the mouse dataset
)

filtered_data$external_gene_name <- rownames(filtered_data)

# Merge the two datasets on external_gene_name
m_filtered_data <- merge(filtered_data, gene_info, by = "external_gene_name", all.x = TRUE)
m_filtered_data <- m_filtered_data %>% dplyr::select(external_gene_name,PValue,logFC,FDR)

colnames(m_filtered_data) <- c("Gene Symbol Mouse","Expr p-value","Expr Log Ratio","Expr False Discovery Rate")
write.xlsx(m_filtered_data,paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\",n,"_ensembl_",Mth,"_nofilter.xlsx"),
                    rowNames = F)

}



####### DEG comparing ##### 20240418
Idents(scrna) = scrna$type
scrna$condition = scrna$orig.ident
scrna$condition[grep("^K", scrna$condition)] <- "KO"
scrna$condition[grep("^W", scrna$condition)] <- "WT"



dif_gene_betas <- FindMarkers(object = scrna, ident.1 = "KO",ident.2 = "WT",subset.ident=c("Beta 1","Beta 2"),
                              group.by = 'condition',
                              test.use = "wilcox")

dif_gene_beta1 <- FindMarkers(object = scrna, ident.1 = "KO",ident.2 = "WT",subset.ident="Beta 1",
                              group.by = 'condition',
                              test.use = "wilcox")
dif_gene_beta2 <- FindMarkers(object = scrna, ident.1 = "KO",ident.2 = "WT",subset.ident="Beta 2",
                              group.by = 'condition',
                              test.use = "wilcox")
# write.xlsx(dif_gene_beta1,paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\","findmarker_KO_WT_Beta1.xlsx"),
#            rowNames = T)
# write.xlsx(dif_gene_beta2,paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\","findmarker_KO_WT_Beta2.xlsx"),
#            rowNames = T)
dif_gene_beta1 = rownames(dif_gene_beta1 %>% dplyr::filter(p_val < 0.05))
dif_gene_beta2 = rownames(dif_gene_beta2 %>% dplyr::filter(p_val < 0.05))
final_beta_genes = data.frame(union(dif_gene_beta1,dif_gene_beta2))



for (n in levels(scrna$type)[levels(scrna$type) != "Other immune cells"]) {
  dif_gene <- FindMarkers(object = scrna, ident.1 = "KO",ident.2 = "WT",subset.ident=n,
                          group.by = 'condition',
                          test.use = "wilcox")
  
  gene_info <- getBM(
    attributes = c("external_gene_name", "ensembl_gene_id"), # Specify the desired attributes
    filters = "external_gene_name",                         # Filter by gene name
    values = rownames(dif_gene),                              # Provide the list of gene names
    mart = ensembl                                           # Use the mouse dataset
  )
  
  dif_gene$external_gene_name <- rownames(dif_gene)
  
  # Merge the two datasets on external_gene_name
  m_dif_gene <- merge(dif_gene, gene_info, by = "external_gene_name", all.x = TRUE)
  m_dif_gene <- m_dif_gene %>% dplyr::select(external_gene_name,p_val,avg_log2FC,p_val_adj)
  
  colnames(m_filtered_data) <- c("Gene Symbol Mouse","Expr p-value","Expr Log Ratio","Expr False Discovery Rate")
  
  write.xlsx(m_dif_gene,paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\","findmarker_KO_WT_ensembl_",n,".xlsx"),
             rowNames = F)
}

# ### IPA ###
# dif_gene = dif_gene_betas
# gene_info <- getBM(
#   attributes = c("external_gene_name", "ensembl_gene_id"), # Specify the desired attributes
#   filters = "external_gene_name",                         # Filter by gene name
#   values = rownames(dif_gene),                              # Provide the list of gene names
#   mart = ensembl                                           # Use the mouse dataset
# )
# 
# dif_gene$external_gene_name <- rownames(dif_gene)
# 
# # Merge the two datasets on external_gene_name
# m_dif_gene <- merge(dif_gene, gene_info, by = "external_gene_name", all.x = TRUE)
# m_dif_gene <- m_dif_gene %>% dplyr::select(external_gene_name,p_val,avg_log2FC,p_val_adj)
# 
# colnames(m_filtered_data) <- c("Gene Symbol Mouse","Expr p-value","Expr Log Ratio","Expr False Discovery Rate")
# 
# write.xlsx(m_dif_gene,paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\","findmarker_KO_WT_ensembl_Betas",".xlsx"),
#            rowNames = F)
# #write.xlsx(final_beta_genes,
# #           "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\Union_betas.xlsx")
