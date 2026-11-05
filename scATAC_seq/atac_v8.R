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
library(promises)
library(openxlsx)
library(RColorBrewer)
library(ggrepel)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

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

frag = function(dir,combine_peak,md){
  setwd(dir)
  res = CreateFragmentObject(path = "./fragments.tsv.gz",cells = rownames(md))
  
  res.FM = FeatureMatrix(
    fragments = res,
    features = combine_peak,
    cells = rownames(res)
  )
  
  res = CreateChromatinAssay(res.FM, fragments = res)
  
  res = CreateSeuratObject(res, assay = "ATAC", meta.data=md)
  
  return(res)
}

import_atac <- function(dir,count_path, meta_path, fragment_path){
  setwd(dir)
  
  counts <- Read10X_h5(filename = count_path)
  
  meta <- read.csv(
    file = meta_path,
    header = TRUE,
    row.names = 1)
  

  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = fragment_path,
    min.cells = 10,
    min.features = 200
  )
  
  data <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = meta
  )
  data = UpdateSeuratObject(data)
  ## get gene annotation for mm10
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "mm10"
  Annotation(data) <- annotations
  
  ## QCs
  
  # compute nucleosome signal score per cell
  data <- NucleosomeSignal(object = data) #fragment ratio 147-294: <147  ---  mononucleosome:nucleosome-free
  # compute TSS enrichment score per cell
  data <- TSSEnrichment(object = data, fast = FALSE)
  # add blacklist ratio and fraction of reads in peaks
  data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments
  data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100 
  gc()
  
  
  data$high.tss <- ifelse(data$TSS.enrichment > 2, 'High', 'Low')
  TSSPlot(data, group.by = 'high.tss') + NoLegend()
  
  data$nucleosome_group <- ifelse(data$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  FragmentHistogram(object = data, group.by = 'nucleosome_group')
  
  
  return(data)
}

W1 <- import_atac("C:/Users/maili/Box/Data Club/Mai Liu/Computational analysis/scATAC-seq/Inputs/ATAC_seurat/W1",
                  "filtered_peak_bc_matrix.h5",
                     'singlecell.csv',
                     "./fragments.tsv.gz")

W2 <- import_atac("C:/Users/maili/Box/Data Club/Mai Liu/Computational analysis/scATAC-seq/Inputs/ATAC_seurat/W2/",
                  "filtered_peak_bc_matrix.h5",
                  'singlecell.csv',
                  "./fragments.tsv.gz")
KO1 <- import_atac("C:/Users/maili/Box/Data Club/Mai Liu/Computational analysis/scATAC-seq/Inputs/ATAC_seurat/KO1/",
                   "filtered_peak_bc_matrix.h5",
                    'singlecell.csv',
                    "./fragments.tsv.gz")
KO2 <- import_atac("C:/Users/maili/Box/Data Club/Mai Liu/Computational analysis/scATAC-seq/Inputs/ATAC_seurat/KO2/",
                   "filtered_peak_bc_matrix.h5",
                   'singlecell.csv',
                   "./fragments.tsv.gz")

W1$dataset <- "W1"
W2$dataset <- "W2"
KO1$dataset <- "KO1"
KO2$dataset <- "KO2"




# ### sp 1
# save(W1, file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V7\\W1_sp1.rda")
# save(W2, file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V7\\W2_sp1.rda")
# save(KO1, file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V7\\KO1_sp1.rda")
# save(KO2, file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V7\\KO2_sp1.rda")
load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V7\\W1_sp1.rda")
load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V7\\W2_sp1.rda")
load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V7\\KO1_sp1.rda")
load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V7\\KO2_sp1.rda")

# combined <- merge(
#   x = W1,
#   y = list(W2, KO1, KO2),
#   add.cell.ids = c("W1", "W2", "KO1", "KO2")
# )
# 
# data.frame(combined@meta.data) %>%
#   ggplot(aes(color=dataset, x=peak_region_fragments, fill= dataset)) +
#   geom_density(alpha = 0.2) +
#   scale_x_log10() +
#   theme_classic() +
#   geom_vline(xintercept = 50000)
# 
# data.frame(combined@meta.data) %>%
#   ggplot(aes(color=dataset, x=pct_reads_in_peaks, fill= dataset)) +
#   geom_density(alpha = 0.2) +
#   theme_classic() +
#   scale_x_log10() +
#   geom_vline(xintercept = 35)
# 
# data.frame(combined@meta.data) %>%
#   ggplot(aes(color=dataset, x=blacklist_ratio, fill=dataset)) +
#   geom_density(alpha = 0.2) +
#   scale_x_log10() +
#   theme_classic() +
#   geom_vline(xintercept = 0.05)
# 
# data.frame(combined@meta.data) %>%
#   ggplot(aes(color=dataset, x=nucleosome_signal, fill=dataset)) +
#   geom_density(alpha = 0.2) +
#   scale_x_log10() +
#   theme_classic() +
#   geom_vline(xintercept = 2)
# 
# data.frame(combined@meta.data) %>%
#   ggplot(aes(color=dataset, x=TSS.enrichment, fill=dataset)) +
#   geom_density(alpha = 0.2) +
#   scale_x_log10() +
#   theme_classic() +
#   geom_vline(xintercept = 2.5)
# 
# combined <- subset(
#   x = combined,
#   subset = peak_region_fragments > 3000 &
#     peak_region_fragments < 50000 &
#     pct_reads_in_peaks > 35 &
#     blacklist_ratio < 0.05 &
#     nucleosome_signal < 2 &
#     TSS.enrichment > 2.5
# )

W1 = subset(
  x = W1,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 50000 &
    pct_reads_in_peaks > 35 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2.5
)
W2 = subset(
  x = W2,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 50000 &
    pct_reads_in_peaks > 35 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2.5
)
KO1 = subset(
  x = KO1,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 50000 &
    pct_reads_in_peaks > 35 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2.5
)
KO2 = subset(
  x = KO2,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 50000 &
    pct_reads_in_peaks > 35 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2.5
)



#Make combined file, removing overlapping peaks
combined.peaks <- GenomicRanges::reduce(x = c(granges(W1), granges(W2), granges(KO1), granges(KO2)))


#Filter out v small and large peaks
combined.peaks <- combined.peaks[width(combined.peaks) <10000 & width(combined.peaks) >20]




W1 = frag("C:/Users/maili/Box/Data Club/Mai Liu/Computational analysis/scATAC-seq/Inputs/ATAC_seurat/W1",
          combined.peaks,W1@meta.data)
          
W2 = frag("C:/Users/maili/Box/Data Club/Mai Liu/Computational analysis/scATAC-seq/Inputs/ATAC_seurat/W2",
          combined.peaks,W2@meta.data)
KO1 = frag("C:/Users/maili/Box/Data Club/Mai Liu/Computational analysis/scATAC-seq/Inputs/ATAC_seurat/KO1",
           combined.peaks,KO1@meta.data)
KO2 = frag("C:/Users/maili/Box/Data Club/Mai Liu/Computational analysis/scATAC-seq/Inputs/ATAC_seurat/KO2",
           combined.peaks,KO2@meta.data)

### sp 2
# save(W1, file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V5\\Stop_point\\W1_sp2.rda")
# save(W2, file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V5\\Stop_point\\W2_sp2.rda")
# save(KO1, file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V5\\Stop_point\\KO1_sp2.rda")
# save(KO2, file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V5\\Stop_point\\KO2_sp2.rda")

# load("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V5\\Stop_point\\W1_sp2.rda")
# load("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V5\\Stop_point\\W2_sp2.rda")
# load("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V5\\Stop_point\\KO1_sp2.rda")
# load("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\scATAC-seq\\Outputs\\V5\\Stop_point\\KO2_sp2.rda")
 

#########################20240319##############################

### Seurat Inte ###
W1 <- RenameCells(object = W1, add.cell.id = "W1")
W2 <- RenameCells(object = W2, add.cell.id = "W2")
KO1 <- RenameCells(object = KO1, add.cell.id = "KO1")
KO2 <- RenameCells(object = KO2, add.cell.id = "KO2")


W1 <- RunTFIDF(W1)
W1 <- FindTopFeatures(W1, min.cutoff = 20)
W1 <- RunSVD(W1)
W2 <- RunTFIDF(W2)
W2 <- FindTopFeatures(W2, min.cutoff = 20)
W2 <- RunSVD(W2)
KO1 <- RunTFIDF(KO1)
KO1 <- FindTopFeatures(KO1, min.cutoff = 20)
KO1 <- RunSVD(KO1)
KO2 <- RunTFIDF(KO2)
KO2 <- FindTopFeatures(KO2, min.cutoff = 20)
KO2 <- RunSVD(KO2)

unintegrated <- merge(
  x = W1,
  y = list(W2, KO1, KO2)
  
)


DefaultAssay(unintegrated) <- "ATAC"
unintegrated <- RunTFIDF(unintegrated)
unintegrated <- FindTopFeatures(unintegrated, min.cutoff = 50)
unintegrated <- RunSVD(unintegrated)
# unintegrated <- RunUMAP(unintegrated, reduction = 'lsi', dims = 2:30,seed.use = 42)
p1 <- DimPlot(unintegrated, group.by = 'dataset', pt.size = 0.1) + ggplot2::ggtitle("Unintegrated")
p1

integration.anchors <- FindIntegrationAnchors(
  object.list = list(W1,W2,KO1,KO2),
  anchor.features = 20000,
  reduction = "rlsi",
  dims = 2:30
)

integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = unintegrated[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)



# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30,seed.use = 42 )
DimPlot(integrated, group.by = "dataset")



annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) = 'UCSC'
Annotation(integrated) <- annotations
gene.activities = GeneActivity(integrated)
rownames(gene.activities) = str_to_title(tolower(rownames(gene.activities)))
integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(integrated) <- "ACTIVITY"




save(integrated,file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\integration_no_trans_no_embed.rda")
load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\integration_no_trans_no_embed.rda")


### 050624 ###
DefaultAssay(integrated) <- "ACTIVITY"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)
### 050624 ###



transfer.anchors_atac <- FindTransferAnchors(reference = scrna,
                                        query = integrated,
                                        features = (VariableFeatures(object = scrna)),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",normalization.method = "LogNormalize")



predicted.labels_atac <- TransferData(
  anchorset = transfer.anchors_atac,
  refdata = scrna@active.ident,
  weight.reduction = integrated[['integrated_lsi']],
  dims = 2:30
)







integrated <- AddMetaData(integrated, metadata = predicted.labels_atac)


Idents(integrated) = integrated$predicted.id
integrated$orig.ident = integrated$dataset
levels(integrated) = levels(scrna)

DimPlot_scCustom(integrated)

levels(integrated) = rev(levels(integrated))

DefaultAssay(integrated) <- "ACTIVITY"
DotPlot(integrated,features = gene_list_plot) + 
  theme_minimal() +  # Change the overall theme to a minimal style
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels vertically
    axis.title.y = element_blank(),  # Remove y-axis title
    strip.background = element_blank(),  # Remove background of facet labels
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  scale_color_gradient(low = "yellow", high = "purple3")
levels(integrated) = rev(levels(integrated))

#save(integrated,file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\integration_NoFilter.rda")
#load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\integration_NoFilter.rda")

#########################
### filter < 50 cells ###
#########################
filt_thres = table(integrated$predicted.id) > 50
remain_clusters = names(filt_thres[filt_thres])
integrated = subset(x = integrated,subset = (predicted.id %in% remain_clusters))
levels(integrated) = levels(scrna)
save(integrated,file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\integration.rda")
#load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\integration.rda")

DimPlot(integrated,label = T,split.by = "orig.ident") + NoLegend()
DimPlot(integrated,label = T,split.by = "predicted.id") + NoLegend()

DimPlot(integrated,label = T) + NoLegend()
color_select = DiscretePalette_scCustomize(num_colors = 13, palette = "polychrome")

p1 = DimPlot_scCustom(integrated,colors_use = color_select) + NoLegend()
p1 = DimPlot_scCustom(integrated) + NoLegend()

p2 = DimPlot_scCustom(scrna) + NoLegend()
p1 | p2




#########################20240319##############################

integrated$tech <- "ATAC"
scrna$tech <- "RNA"
genes.use <- VariableFeatures(scrna)
refdata <- GetAssayData(scrna,assay = "RNA",slot = "data")[genes.use,]
imputation <- TransferData(anchorset = transfer.anchors_atac,refdata = refdata,
                           weight.reduction = integrated[["integrated_lsi"]],dims = 2:30)
integrated[["RNA"]] <- imputation
coembed <- merge(x = scrna,y = integrated)
coembed <- ScaleData(coembed,features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = F)
coembed <- RunUMAP(coembed,dims = 1:30,seed.use = 42)
coembed$type <- ifelse(!is.na(coembed$type),coembed$type,coembed$predicted.id)

DimPlot_scCustom(coembed)
coembed$condition = coembed$orig.ident
coembed$condition[grep("^K", coembed$condition)] <- "KO"
coembed$condition[grep("^W", coembed$condition)] <- "WT"
coembed[["RNA"]] <- JoinLayers(coembed[["RNA"]])
save(coembed,file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\RNA_ATAC_embed.rda")

load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\RNA_ATAC_embed.rda")


DefaultAssay(coembed) <- 'RNA'

Stacked_VlnPlot(seurat_object = coembed, features = gene_list_plot, x_lab_rotate = TRUE)
r = subset(x = coembed,subset = (tech == "RNA"))
a = subset(x = coembed,subset = (tech == "ATAC"))

DimPlot(coembed,label = T) + NoLegend()
DimPlot(coembed,split.by = "orig.ident",label = T) + NoLegend()
DimPlot(coembed,split.by = "tech",label = T) + NoLegend()
DefaultAssay(coembed) <- 'RNA'
Stacked_VlnPlot(seurat_object = coembed, features = gene_list_plot, x_lab_rotate = TRUE)

n_color = hue_pal()(length(levels(r)))

p1 = sample_identity(r,n_color)

#n_color = n_color[-which(!levels(r) %in% levels(a))]

p2 = sample_identity(a,n_color)

p1 | p2

library(dplyr)
library(tidyr)
library(openxlsx)

# =========================
# Export individual population percentages
# =========================

get_individual_population_percent <- function(seu_obj, tech_name, celltype_col = "type") {
  
  md <- seu_obj@meta.data %>% as.data.frame()
  
  # Use cell type from metadata if available; otherwise use active identity
  if (celltype_col %in% colnames(md)) {
    cell_type_vec <- md[[celltype_col]]
    cell_type_levels <- if (is.factor(cell_type_vec)) {
      levels(cell_type_vec)
    } else {
      unique(as.character(cell_type_vec))
    }
    md$Cell_type <- factor(as.character(cell_type_vec), levels = cell_type_levels)
  } else {
    cell_type_levels <- levels(seu_obj)
    md$Cell_type <- factor(as.character(Idents(seu_obj)), levels = cell_type_levels)
  }
  
  # Keep sample order
  sample_levels <- c("W1", "W2", "KO1", "KO2")
  sample_levels <- c(
    sample_levels[sample_levels %in% unique(md$orig.ident)],
    setdiff(unique(md$orig.ident), sample_levels)
  )
  
  # Count each cell type in each individual sample
  count_long <- md %>%
    mutate(
      Sample = as.character(orig.ident),
      Condition = ifelse(grepl("^K", Sample), "KO", "WT")
    ) %>%
    count(Sample, Condition, Cell_type, name = "Cell_count") %>%
    complete(
      Sample = sample_levels,
      Cell_type = cell_type_levels,
      fill = list(Cell_count = 0)
    ) %>%
    mutate(
      Condition = ifelse(grepl("^K", Sample), "KO", "WT"),
      Sample = factor(Sample, levels = sample_levels),
      Cell_type = factor(Cell_type, levels = cell_type_levels)
    ) %>%
    group_by(Sample) %>%
    mutate(
      Total_cells = sum(Cell_count),
      Fraction = Cell_count / Total_cells,
      Percent = Fraction * 100
    ) %>%
    ungroup() %>%
    mutate(Tech = tech_name, .before = 1)
  
  # Wide percent table: rows = W1/W2/KO1/KO2, columns = cell types
  percent_wide <- count_long %>%
    select(Tech, Sample, Condition, Total_cells, Cell_type, Percent) %>%
    pivot_wider(
      names_from = Cell_type,
      values_from = Percent,
      values_fill = 0
    )
  
  # Wide count table
  count_wide <- count_long %>%
    select(Tech, Sample, Condition, Total_cells, Cell_type, Cell_count) %>%
    pivot_wider(
      names_from = Cell_type,
      values_from = Cell_count,
      values_fill = 0
    )
  
  # Mean and SD used for WT/KO stacked barplot
  mean_sd <- count_long %>%
    group_by(Tech, Condition, Cell_type) %>%
    summarise(
      Mean_percent = mean(Percent, na.rm = TRUE),
      SD_percent = sd(Percent, na.rm = TRUE),
      Mean_fraction = mean(Fraction, na.rm = TRUE),
      SD_fraction = sd(Fraction, na.rm = TRUE),
      n_samples = n_distinct(Sample),
      .groups = "drop"
    )
  
  return(list(
    long = count_long,
    percent_wide = percent_wide,
    count_wide = count_wide,
    mean_sd = mean_sd
  ))
}

# Make sure identities are cell types
Idents(r) <- r$type
Idents(a) <- a$type

# Calculate RNA and ATAC percentages
rna_pop <- get_individual_population_percent(r, tech_name = "RNA")
atac_pop <- get_individual_population_percent(a, tech_name = "ATAC")

# Combine long-format and mean/SD tables
all_long <- bind_rows(rna_pop$long, atac_pop$long)
all_mean_sd <- bind_rows(rna_pop$mean_sd, atac_pop$mean_sd)

# Output folder
out_dir <- "C:/Users/maili/Box/Data Club/Mai Liu/Computational analysis/sc_Mai_outputs/20260423/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "RNA_ATAC_population_percent_by_individual.xlsx")

# Write Excel
wb <- createWorkbook()

addWorksheet(wb, "RNA_percent_wide")
writeData(wb, "RNA_percent_wide", rna_pop$percent_wide)

addWorksheet(wb, "ATAC_percent_wide")
writeData(wb, "ATAC_percent_wide", atac_pop$percent_wide)

addWorksheet(wb, "RNA_counts_wide")
writeData(wb, "RNA_counts_wide", rna_pop$count_wide)

addWorksheet(wb, "ATAC_counts_wide")
writeData(wb, "ATAC_counts_wide", atac_pop$count_wide)

addWorksheet(wb, "All_long")
writeData(wb, "All_long", all_long)

addWorksheet(wb, "Mean_SD_for_barplot")
writeData(wb, "Mean_SD_for_barplot", all_mean_sd)

saveWorkbook(wb, out_file, overwrite = TRUE)

print(paste("Saved to:", out_file))

####################### HERE ########################### 20240418
DefaultAssay(coembed) <- 'ATAC'



# da_peaks <- FindMarkers(object = coembed, ident.1 = "WT",ident.2 = "KO",subset.ident="Beta 1",
#                         group.by = 'condition',
#                         test.use = 'LR', latent.vars = 'peak_region_fragments')
# 
# head(da_peaks)
# save(da_peaks,file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-ATAC\\Outputs\\V7\\DiffPeaks.rda")
# #load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-ATAC\\Outputs\\V7\\DiffPeaks.rda")
# 
# open_WT <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
# open_KO <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])
# 
# closest_genes_open_WT <- ClosestFeature(coembed, regions = open_WT)
# closest_genes_open_KO <- ClosestFeature(coembed, regions = open_KO)
# 
# CoveragePlot(
#   object = coembed,
#   region = rownames(da_peaks)[1],
#   extend.upstream = 40000,
#   extend.downstream = 20000
# )




# # Get a list of motif position frequency matrices from the JASPAR database
# pfm <- getMatrixSet(
#   x = JASPAR2020,
#   opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
# )
# 
# # add motif information
# coembed <- AddMotifs(
#   object = coembed,
#   genome = BSgenome.Mmusculus.UCSC.mm10,
#   pfm = pfm
# )
# 
da_peaks_motif <- FindMarkers(object = coembed, ident.1 = "KO",ident.2 = "WT",subset.ident="Beta 1",
                        group.by = 'condition',
                        test.use = 'LR', min.pct = 0.05,latent.vars = 'nCount_ATAC')
# save(da_peaks_motif,file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-ATAC\\Outputs\\V7\\DiffPeaks_Motif.rda")
# #load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-ATAC\\Outputs\\V7\\DiffPeaks_Motif.rda")
# 
# 
# 
top.da.peak <- rownames(da_peaks_motif[da_peaks_motif$p_val < 0.005 & da_peaks_motif$pct.1 > 0.2, ])
enriched.motifs <- FindMotifs(
  object = coembed,
  features = top.da.peak
)
# 
# 
# MotifPlot(
#   object = coembed,
#   motifs = head(rownames(enriched.motifs))
# )

### ChromVAR ###
library(BiocParallel)
library(chromVAR)
register(SerialParam())

coembed <- RunChromVAR(
  object = coembed,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(coembed) <- 'chromvar'
save(coembed,file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-ATAC\\Outputs\\V7\\chromVAR.rda")
#load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-ATAC\\Outputs\\V7\\chromVAR.rda")
motif.names = GetMotifData(object = coembed, assay = "ATAC", slot = "motif.names")



# differential.activity <- FindMarkers(object = coembed, ident.1 = "KO",ident.2 = "WT",subset.ident=c("Beta 1"),
#                                      group.by = 'condition',
#                                      test.use = 'LR',
#                                      latent.vars = 'nCount_ATAC')
# 
# differential.activity <- FindMarkers(object = coembed, ident.1 = "KO",ident.2 = "WT",subset.ident=c("Beta 1"),
#                                      group.by = 'condition',
#                                      test.use = 'LR', min.pct = 0.05)
## signac ##
differential.activity <- FindMarkers(object = coembed, ident.1 = "KO",ident.2 = "WT",subset.ident=c("Beta 1"),
                                     group.by = 'condition',
                                     mean.fxn = rowMeans,
                                     fc.name = "avg_diff")

motif_lookup <- rownames(differential.activity)
motif_names <- motif.names[motif_lookup]
differential.activity = cbind(differential.activity, gene = unlist(motif_names))
view(differential.activity)



#write.xlsx(differential.activity,paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\Motif_KOvsWT_analysis_Beta1.xlsx"),
#           rowNames = T)

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
differential.activity$delabel <- NA
differential.activity$diffexpressed <- "Not significant"
differential.activity$diffexpressed[differential.activity$avg_diff   > 0.25 & differential.activity$p_val < 0.05] <- "Upregulated"
differential.activity$diffexpressed[differential.activity$avg_diff   < -0.25 & differential.activity$p_val < 0.05] <- "Downregulated"
differential.activity$delabel[which(differential.activity$gene == "BHLHE40")] = "BHLHE40"
differential.activity$delabel[which(differential.activity$gene == "MAFA")] = "MAFA"
#differential.activity$chrom = differential.activity$gene

ggplot(differential.activity, aes(x = avg_diff  , y = -log10(p_val), col = diffexpressed, label = delabel)) +
  geom_point(size = 0.1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "black") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", col = "black") +
  labs(title = "", x = "avg_diff ", y = "-log10(p-value)") +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  theme_classic()+ # Plot title   
  geom_text_repel(max.overlaps = Inf,show.legend =F) # To show all labels 

# write.xlsx(differential.activity %>% dplyr::filter(p_val<= 0.05 & avg_log2FC < 0.25) %>% dplyr::select(gene),paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\DR_KO_analysis_Beta1.xlsx"),
#            rowNames = T)






### Chi-square significance of cell identity
cell_info <- coembed@meta.data %>% dplyr::select(c(tech,type,condition,orig.ident))
for (j in unique(cell_info$tech)) {
  print(j)
  total_cell = colSums(as.matrix(table((cell_info %>% dplyr::filter(tech == j) %>% dplyr::select(c(type,condition))))))
  
  for (i in levels(coembed)) {
    cell_info2 <- cell_info %>% dplyr::filter(tech == j) %>% dplyr::select(c(type,condition)) %>% dplyr::filter(type == i)
    cell_info2 <- table(cell_info2)
    cell_info2 = rbind(cell_info2,total_cell - cell_info2)
    rownames(cell_info2) = c(i,paste("No",i))
    print(i)
    print(chisq.test(cell_info2))  # Prints test summary
  }
}

p_value_RNA = c(2.2e-16,0.001716,0.0008559,0.08321,0.01507,0.9264,0.01936,0.5184,1.363e-05,0.5797,0.2729,2.2e-16,2.2e-16,2.2e-16,2.2e-16,0.1415)
p_value_RNA = scale(-log10(p_value_RNA))
p_value_RNA = data.frame(cbind(p_value_RNA,p_value_RNA))
rownames(p_value_RNA) = levels(coembed)

p_value_ATAC = c(2.2e-16,0.001963,1.549e-06,0.103,0.002846,0.9547,0.00233,0.1845,0.002652,0.0925,0.08376,2.2e-16,2.2e-16,2.2e-16,2.2e-16,0.88)
p_value_ATAC = scale(-log10(p_value_ATAC))
p_value_ATAC = data.frame(cbind(p_value_ATAC,p_value_ATAC))
rownames(p_value_ATAC) = levels(coembed)
# Create a matrix with two rows (one for RNA) and the number of columns equal to the length of p_value_RNA

p1 = pheatmap(as.matrix(p_value_RNA), 
         scale = "none",  # No scaling
         cluster_rows = FALSE,  # Do not cluster rows
         cluster_cols = FALSE,  # Do not cluster columns
         #annotation_col = NULL,  # Remove column annotation
         #annotation_row = NULL,  # Remove row annotation
         show_colnames = F,
         border_color = NA  # Remove border
)
p2 = pheatmap(as.matrix(p_value_ATAC), 
         scale = "none",  # No scaling
         cluster_rows = FALSE,  # Do not cluster rows
         cluster_cols = FALSE,  # Do not cluster columns
         #annotation_col = NULL,  # Remove column annotation
         #annotation_row = NULL,  # Remove row annotation
         show_colnames = F,
         border_color = NA  # Remove border
)



n_color = hue_pal()(length(levels(scrna)))

sample_identity(scrna,n_color)

n_color = n_color[-which(!levels(scrna) %in% levels(combined))]

sample_identity(combined,n_color)

### RNA/ATAC integration ###

