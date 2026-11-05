################
### CellChat ###
################
load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-RNA\\Outputs\\V7\\Seurat\\Log_integrated_RNA.rda")
scrna[["RNA3"]] <- as(object = scrna[["RNA"]], Class = "Assay")
DefaultAssay(scrna) <- 'RNA3'

# Merging #
# Idents(scrna)[which(Idents(scrna) == "Beta 2")] <- "Beta 1"
# scrna$type[which(scrna$type == "Beta 2")] = "Beta 1"

#select_type = c("Beta 1","Beta 2","Alpha 1","Alpha 2","Delta 1", "Delta 2", "Gamma")
select_type = c("Beta 1","Beta 2","Alpha 1","Alpha 2","Delta 1", "Delta 2", "Gamma","Stellate 1","Stellate 2","Vascular cells","Proliferating cells")
#select_type = c("Beta 1","Alpha 1","Alpha 2","Delta 1", "Delta 2", "Gamma","Stellate 1","Stellate 2","Vascular cells","Proliferating cells")

## WT ##
scrna_WT = subset(scrna,(orig.ident == "W1" | orig.ident == "W2") & type %in% select_type)
data.input <- GetAssayData(scrna_WT,assay = "RNA3", slot = "data")
labels <- Idents(scrna_WT)
meta <- data.frame(group = labels,row.names = names(labels))
cellchat_WT = createCellChat(object = data.input, meta = meta, group.by = c("group"))



CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)
### Specific pathway ###
# CellChatDB <- subsetDB(CellChatDB.mouse,search = "Secreted Signaling")
# CellChatDB <- subsetDB(CellChatDB.mouse,search = "ECM-Receptor")
# CellChatDB <- subsetDB(CellChatDB.mouse,search = "Cell-Cell Contact")
CellChatDB.use <- CellChatDB
cellchat_WT@DB <- CellChatDB.use
cellchat_WT <- subsetData(cellchat_WT)
cellchat_WT <- identifyOverExpressedGenes(cellchat_WT)
cellchat_WT <- identifyOverExpressedInteractions(cellchat_WT)
cellchat_WT <- smoothData(cellchat_WT, adj = PPI.mouse)
cellchat_WT <- computeCommunProb(cellchat_WT, raw.use = T)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_WT <- filterCommunication(cellchat_WT, min.cells = 10)
cellchat_WT <- computeCommunProbPathway(cellchat_WT)
cellchat_WT <- aggregateNet(cellchat_WT)

# # 1. Compute the network centrality scores
cellchat_WT <- netAnalysis_computeCentrality(cellchat_WT, slot.name = "netP")


# 1.
## KO ##
scrna_KO = subset(scrna,(orig.ident == "KO1" | orig.ident == "KO2") & type %in% select_type)
data.input <- GetAssayData(scrna_KO,assay = "RNA3", slot = "data")
labels <- Idents(scrna_KO)
meta <- data.frame(group = labels,row.names = names(labels))
cellchat_KO = createCellChat(object = data.input, meta = meta, group.by = c("group"))

# 2. Add the Secreted Signaling database in the CellChat object
cellchat_KO@DB <- CellChatDB.use

# 3.Subset and pre-processing the expression data
cellchat_KO <- subsetData(cellchat_KO)
cellchat_KO <- identifyOverExpressedGenes(cellchat_KO)
cellchat_KO <- identifyOverExpressedInteractions(cellchat_KO)

# 4. project gene expression data onto protein-protein interaction (PPI)
cellchat_KO <- smoothData(cellchat_KO, adj = PPI.mouse) # PPI.mouse for mouse samples

# 5. Compute the communication probability and infer cellular communication network
cellchat_KO <- computeCommunProb(cellchat_KO, raw.use = FALSE)

# 6. Filter out the cell-cell communication min.cells = 10
cellchat_KO <- filterCommunication(cellchat_KO, min.cells = 10)

# 7. Infer the cell-cell communication at a signaling pathway level
cellchat_KO <- computeCommunProbPathway(cellchat_KO)

# 8. Calculate the aggregated cell-cell communication network
cellchat_KO <- aggregateNet(cellchat_KO)

# # 9. Compute the network centrality scores
cellchat_KO <- netAnalysis_computeCentrality(cellchat_KO, slot.name = "netP")

# # 14. Save  cellchat object
saveRDS(cellchat_WT,
        file="C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\Cellchat\\cellchat_non_immune_WT.rds")
saveRDS(cellchat_KO,
        file="C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\Cellchat\\cellchat_non_immune_KO.rds")


















#15. load both NL and KO objects
cellchat_WT <- readRDS(
  "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\Cellchat\\cellchat_non_immune_WT.rds")

cellchat_KO <- readRDS(
  "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-RNA\\Cellchat\\cellchat_non_immune_KO.rds")

object.list <- list(WT = cellchat_WT, KO = cellchat_KO)

names(object.list)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# pdf(paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-RNA\\Cellchat\\WT_KO_NumWeig.pdf"),width = 20,
#     height = 15)



gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
# dev.off()

# pdf(paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-RNA\\Cellchat\\WT_KO_NumWeig_cir.pdf"),width = 20,
#     height = 10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
# dev.off()

# pdf(paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-RNA\\Cellchat\\WT_KO_NumWeig_heatmap.pdf"),width = 20,
#     height = 10)
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
# dev.off()


# 1. Compare the overall information flow of each signaling pathway
# pdf(paste0("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240220\\sc-RNA\\Cellchat\\WT_KO_Pathways.pdf"),width = 10,
#     height = 15)
p1 = rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) # 42 & 47
# dev.off()
p2 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p1 + p2


netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count",label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "weight",label.edge = T)



pathways.show.all <- unique((p1$data %>% dplyr::filter(pvalues < 0.01) %>% arrange(desc(contribution)))$name)[1:10]


# i = 1
# pathways.show.all <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

# Pvalue
cellTypeSelect = 1 # Beta 1 = 1, Beta 2 = 2

# Initialize an empty list to store the computed p-value ratios
dif_pval = list()
ct = 1  # Counter to keep track of list indices
# Loop through top pathways pair in the WT dataset
for (i in cellchat_WT@LR$LRsig$interaction_name[which(cellchat_WT@LR$LRsig$pathway_name %in% pathways.show.all)]) {
  if ((i %in% cellchat_WT@net$LRs & !(i %in% cellchat_KO@net$LRs))) {
    if (any(cellchat_WT@net$pval[,,i] < 0.05)) {
      a <- cellchat_WT@net$pval[,,i] < 0.05
      ratio <- ifelse(a, -Inf, NaN)
      dif_pval[[i]] = ratio
    }
  }
  # Check if the current ligand-receptor pair is also in the both dataset
  else if (i %in% cellchat_KO@net$LRs & i %in% cellchat_WT@net$LRs) {
    a_both_sig <- (cellchat_WT@net$pval[,,i] < 0.05 & cellchat_KO@net$pval[,,i] < 0.05)
    a_wt_sig_ko_nonsig <- (cellchat_WT@net$pval[,,i] < 0.05 & cellchat_KO@net$pval[,,i] >= 0.05)
    a_wt_nonsig_ko_sig <- (cellchat_WT@net$pval[,,i] >= 0.05 & cellchat_KO@net$pval[,,i] < 0.05)
    a_both_nonsig <- (cellchat_WT@net$pval[,,i] >= 0.05 & cellchat_KO@net$pval[,,i] >= 0.05)
    ratio <- matrix(NaN, nrow = dim(cellchat_WT@net$pval)[1], ncol = dim(cellchat_WT@net$pval)[2])
    # Convert the matrix to a dataframe
    
    colnames(ratio) <- colnames(cellchat_WT@net$pval)
    rownames(ratio) <- rownames(cellchat_WT@net$pval)
    
    ratio[a_both_sig] <- -log2(cellchat_KO@net$pval[,,i] / cellchat_WT@net$pval[,,i])
    ratio[a_wt_sig_ko_nonsig] <- -Inf
    ratio[a_wt_nonsig_ko_sig] <- Inf
    ratio[a_both_nonsig] <- NaN
    
    dif_pval[[i]] = ratio
  }
  # Increment the counter
  ct = ct + 1
}


# Prepare data for the heatmap
# Assuming each list element is a matrix of the same size
# Convert the list of matrices to a single data frame
df_list <- lapply(dif_pval, function(x) as.data.frame(x))

# Extract the first row of each data frame in the list and combine them into a new data frame
first_rows <- do.call(rbind, lapply(names(df_list), function(name) {
  df <- df_list[[name]]
  first_row <- df[cellTypeSelect, , drop = FALSE]
  rownames(first_row) <- name
  return(first_row)
}))
# Convert to a data frame
df_combined_first_rows <- as.data.frame(first_rows)
# Find the maximum and minimum finite values
max_value <- max(df_combined_first_rows[is.finite(as.matrix(df_combined_first_rows))], na.rm = TRUE)
min_value <- min(df_combined_first_rows[is.finite(as.matrix(df_combined_first_rows))], na.rm = TRUE)
if (any(df_combined_first_rows == Inf) & max_value == 0 & min_value != 0) {
  max_value = -min_value
} else if (any(df_combined_first_rows == -Inf) & min_value == 0 & max_value != 0) {
  min_value = -max_value
  
} else if (any(df_combined_first_rows == -Inf) & any(df_combined_first_rows == Inf) & min_value == 0 & max_value == 0) {
  min_value = -1
  max_value = 1
}

# Replace Inf with min_value, -Inf with max_value, and NaN with 0
df_combined_first_rows[is.infinite(as.matrix(df_combined_first_rows)) & as.matrix(df_combined_first_rows) < 0] <- 1.2*min_value
df_combined_first_rows[is.infinite(as.matrix(df_combined_first_rows)) & as.matrix(df_combined_first_rows) > 0] <- 1.2*max_value
df_combined_first_rows[is.nan(as.matrix(df_combined_first_rows))] <- 0

## more negative -> smaller p in WT
data_range <- range((df_combined_first_rows), na.rm = TRUE)
range_length <- diff(data_range)
num_colors <- ceiling(range_length / 0.1)

# Melt the data frame to long format
df_combined_first_rows = cbind(rownames(df_combined_first_rows),df_combined_first_rows)
names(df_combined_first_rows)[1] = "Pairs"
df_long <- melt(df_combined_first_rows)

pa = lapply(df_long$Pairs, function(x) {
  
  nu = which(cellchat_WT@LR$LRsig$interaction_name == x)
  return(cellchat_WT@LR$LRsig$pathway_name[nu])
  
})
df_long <- cbind(data.frame(Pathway = unlist(pa)), data.frame(Pairs = df_long$Pairs),data.frame(variable = df_long$variable, value = df_long$value))


# Ensure Pathway and Pairs are factors
df_long$Pathway <- factor(df_long$Pathway, levels = unique(df_long$Pathway))
df_long$Pairs <- factor(df_long$Pairs, levels = unique(df_long$Pairs))

df_long <- df_long %>%
  group_by(Pairs) %>%
  dplyr::filter(any(value != 0)) %>%
  ungroup()


# Create the heatmap with ggplot2
ggplot(df_long, aes(x = variable, y = Pairs, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(name = "-log2(KO/WT)",low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "P value",
       x = "Samples",
       y = "Ligand-Receptor Pairs",
       fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(Pathway~., scales = "free", space = "free_y", switch = "y")



# Strength (Pvalue)
# Initialize an empty list to store the results
dif_pval = list()

# Loop through top pathways pair in the WT dataset
for (i in cellchat_WT@LR$LRsig$interaction_name[which(cellchat_WT@LR$LRsig$pathway_name %in% pathways.show.all)]) {
  if (i %in% cellchat_WT@net$LRs & !(i %in% cellchat_KO@net$LRs)) {
    if (any(cellchat_WT@net$pval[,,i] < 0.05)) {
      a <- (cellchat_WT@net$pval[,,i] < 0.05 & cellchat_WT@net$prob[,,i] > 0)
      ratio <- ifelse(a, -Inf, NaN)
      dif_pval[[i]] = ratio
    }
  } else if (i %in% cellchat_KO@net$LRs & i %in% cellchat_WT@net$LRs) {
    a_both_sig <- (cellchat_WT@net$pval[,,i] < 0.05 & cellchat_KO@net$pval[,,i] < 0.05)
    a_wt_sig_ko_nonsig <- (cellchat_WT@net$pval[,,i] < 0.05 & cellchat_KO@net$pval[,,i] >= 0.05)
    a_wt_nonsig_ko_sig <- (cellchat_WT@net$pval[,,i] >= 0.05 & cellchat_KO@net$pval[,,i] < 0.05)
    a_both_nonsig <- (cellchat_WT@net$pval[,,i] >= 0.05 & cellchat_KO@net$pval[,,i] >= 0.05)
    
    ratio <- matrix(NaN, nrow = dim(cellchat_WT@net$pval)[1], ncol = dim(cellchat_WT@net$pval)[2])
    # Convert the matrix to a dataframe
    colnames(ratio) <- colnames(cellchat_WT@net$pval)
    rownames(ratio) <- rownames(cellchat_WT@net$pval)
    
    ratio[a_both_sig] <- log2(cellchat_KO@net$prob[,,i] / cellchat_WT@net$prob[,,i])
    ratio[a_wt_sig_ko_nonsig] <- -Inf
    ratio[a_wt_nonsig_ko_sig] <- Inf
    ratio[a_both_nonsig] <- NaN
    
    dif_pval[[i]] = ratio
  }
}
# Prepare data for the heatmap
# Assuming each list element is a matrix of the same size
# Convert the list of matrices to a single data frame
df_list <- lapply(dif_pval, function(x) as.data.frame(x))

# Extract the first row of each data frame in the list and combine them into a new data frame
first_rows <- do.call(rbind, lapply(names(df_list), function(name) {
  df <- df_list[[name]]
  first_row <- df[cellTypeSelect, , drop = FALSE]
  rownames(first_row) <- name
  return(first_row)
}))
# Convert to a data frame
df_combined_first_rows <- as.data.frame(first_rows)
# Find the maximum and minimum finite values
max_value <- max(df_combined_first_rows[is.finite(as.matrix(df_combined_first_rows))], na.rm = TRUE)
min_value <- min(df_combined_first_rows[is.finite(as.matrix(df_combined_first_rows))], na.rm = TRUE)
if (any(df_combined_first_rows == Inf) & max_value == 0 & min_value != 0) {
  max_value = -min_value
} else if (any(df_combined_first_rows == -Inf) & min_value == 0 & max_value != 0) {
  min_value = -max_value
  
} else if (any(df_combined_first_rows == -Inf) & any(df_combined_first_rows == Inf) & min_value == 0 & max_value == 0) {
  min_value = -1
  max_value = 1
}

# Replace Inf with min_value, -Inf with max_value, and NaN with 0
df_combined_first_rows[is.infinite(as.matrix(df_combined_first_rows)) & as.matrix(df_combined_first_rows) < 0] <- 1.2*min_value
df_combined_first_rows[is.infinite(as.matrix(df_combined_first_rows)) & as.matrix(df_combined_first_rows) > 0] <- 1.2*max_value
df_combined_first_rows[is.nan(as.matrix(df_combined_first_rows))] <- 0

## more negative -> smaller p in WT
data_range <- range((df_combined_first_rows), na.rm = TRUE)
range_length <- diff(data_range)
num_colors <- ceiling(range_length / 0.1)

# Melt the data frame to long format
df_combined_first_rows = cbind(rownames(df_combined_first_rows),df_combined_first_rows)
names(df_combined_first_rows)[1] = "Pairs"
df_long <- melt(df_combined_first_rows)

pa = lapply(df_long$Pairs, function(x) {
  
  nu = which(cellchat_WT@LR$LRsig$interaction_name == x)
  return(cellchat_WT@LR$LRsig$pathway_name[nu])
  
})
df_long <- cbind(data.frame(Pathway = unlist(pa)), data.frame(Pairs = df_long$Pairs),data.frame(variable = df_long$variable, value = df_long$value))


# Ensure Pathway and Pairs are factors
df_long$Pathway <- factor(df_long$Pathway, levels = unique(df_long$Pathway))
df_long$Pairs <- factor(df_long$Pairs, levels = unique(df_long$Pairs))


df_long <- df_long %>%
  group_by(Pairs) %>%
  dplyr::filter(any(value != 0)) %>%
  ungroup()
# Create the heatmap with ggplot2
ggplot(df_long, aes(x = variable, y = Pairs, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(name = "log2(KO/WT)",low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Strength Probability",
       x = "Samples",
       y = "Ligand-Receptor Pairs",
       fill = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(Pathway~., scales = "free", space = "free_y", switch = "y")












# Sample data (replace this with your actual data reading method)
data <- df_combined_first_rows

data <- data %>%
  dplyr::filter(substring(Pairs, 1, regexpr("_", Pairs) - 1) != substring(Pairs, regexpr("_", Pairs) + 1))


data <- data %>%
  dplyr::filter(rowSums(dplyr::select(., -Pairs)) != 0)

# Extract ligand and receptors
data_long <- data %>%
  pivot_longer(cols = -Pairs, names_to = "To", values_to = "Value") %>%
  separate(Pairs, into = c("Ligand", "Receptor1", "Receptor2"), sep = "_", extra = "merge", fill = "right") %>%
  unite("Receptors", Receptor1:Receptor2, na.rm = TRUE, remove = FALSE, sep = "_") %>%
  dplyr::select(Ligand, Receptors, To, Value)

# Define colors for the cell types
cell_colors <- hue_pal()(length(levels((scrna$type)))) # Assuming scrna$type defines the cell types
names(cell_colors) <- levels((scrna$type))

data_long <- data_long %>% dplyr::filter(Value != 0)
data_long <- data_long %>% arrange(Ligand, To)

color_gradient <- colorRamp2(c(min(data_long$Value), 0, max(data_long$Value)), c("blue", "white", "red"))

# Add color information to the data
data_long <- data_long %>%
  mutate(LigandColor = cell_colors["Beta 1"],  # All ligands from Beta 1 are red
         ReceptorColor = cell_colors[To],      # Receptors color based on cell type
         ConnectionColor = color_gradient(Value))

data_long <- data_long %>%
  arrange(ReceptorColor) %>%
  unite("To_Receptors", To, Receptors, sep = "_") %>%
  arrange(ReceptorColor)

colors <- c(data_long$LigandColor, data_long$ReceptorColor)
names(colors) <- c(data_long$Ligand, data_long$To_Receptors)

# Create a new column with shortened target labels
data_long$Short_Receptors <- data_long$To_Receptors
data_long$Short_Receptors <- gsub(".*_", "", data_long$Short_Receptors)

# Clear previous circos plot
par(mfrow = c(1,1), xpd = TRUE)
circos.clear()

# Adjust the plotting area size
circos.par(canvas.xlim = c(-1.4, 1.4), canvas.ylim = c(-1.4, 1.4))

# Plot the chord diagram using the original labels for Ligand and shortened labels for To_Receptors
chordDiagram(data_long, grid.col = colors, col = data_long$ConnectionColor, directional = TRUE, 
             annotationTrack = "grid", big.gap = 10, small.gap = 0,
             preAllocateTracks = list(track.height = 0.06), link.target.prop = FALSE, 
             transparency = 0.5)

# Adding sector labels for clarity
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_index <- CELL_META$sector.index
  
  if (sector_index %in% data_long$To_Receptors) {
    # For targets, use shortened names
    label <- data_long$Short_Receptors[which(data_long$To_Receptors == sector_index)][1]
  } else {
    # For sources, use original names
    label <- sector_index
  }
  
  circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), label,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# Map ReceptorColor to Cell_Type
color_to_cell_type <- c(
  "#F8766D" = "Beta 1",
  "#E68613" = "Beta 2",
  "#CD9600" = "Alpha 1",
  "#ABA300" = "Alpha 2",
  "#7CAE00" = "Delta 1",
  "#0BB702" = "Delta 2",
  "#00BE67" = "Gamma",
  "#00C19A" = "Stellate 1",
  "#00BFC4" = "Stellate 2",
  "#00B8E7" = "Vascular cells",
  "#00A9FF" = "Ductual cells",
  "#8494FF" = "Macrophages",
  "#C77CFF" = "Immune B cells",
  "#ED68ED" = "Immune T cells",
  "#FF61CC" = "Other immune cells",
  "#FF68A1" = "Proliferating cells"
)

data_long <- data_long %>%
  mutate(Cell_Type = color_to_cell_type[ReceptorColor])

# Add custom boundaries around unique ReceptorColor groups
unique_colors <- unique(data_long$ReceptorColor)

for (color in unique_colors) {
  sectors <- data_long$To_Receptors[data_long$ReceptorColor == color]
  sectors <- sectors[sectors %in% get.all.sector.index()]
  highlight.sector(sectors, track.index = 2, col = NA, border = "black", lwd = 2)
}

unique_colors <- unique(data_long$Ligand)

for (color in unique_colors) {
  sectors <- data_long$Ligand[data_long$Ligand == color]
  sectors <- sectors[sectors %in% get.all.sector.index()]
  highlight.sector(sectors, track.index = 2, col = NA, border = "black", lwd = 2)
}

# Add legends for links and cells
lgd_links <- Legend(at = c(min(data_long$Value), 0, max(data_long$Value)), 
                    col_fun = color_gradient, title = "log2(KO/WT)")
lgd_cells <- Legend(labels = names(cell_colors), 
                    legend_gp = gpar(fill = cell_colors), title = "Cell Types")

draw(lgd_links, x = unit(8, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
draw(lgd_cells, x = unit(8, "mm"), y = unit(50, "mm"), just = c("left", "bottom"))


# Initialize variables to track labeled cell types
labeled_senders <- FALSE
labeled_receptors <- vector()
# Adding a new track to label the cell types for senders and receptors
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector_index <- CELL_META$sector.index
  
  # Label senders (Ligands) with "Beta 1" only once
  if (sector_index %in% unique(data_long$Ligand)) {
    if (!labeled_senders) {
      circos.text(CELL_META$xcenter, mean(CELL_META$ylim), "Beta 1", 
                  facing = "clockwise", niceFacing = TRUE, adj = c(0.5, 0.5), col = "black")
      labeled_senders <<- TRUE  # Use global assignment to update flag
    }
  }
  
  # Label receptors with their Cell_Type only once per unique type
  if (sector_index %in% data_long$To_Receptors) {
    cell_type <- unique(data_long$Cell_Type[data_long$To_Receptors == sector_index])
    if (!(cell_type %in% labeled_receptors)) {
      circos.text(CELL_META$xcenter, mean(CELL_META$ylim), cell_type, 
                  facing = "clockwise", niceFacing = TRUE, adj = c(0.5, 0.5), col = "black")
      labeled_receptors <<- c(labeled_receptors, cell_type)  # Update the list of labeled cell types
    }
  }
}, bg.border = NA)


### IPA ###
# 
# IPA = read.xlsx("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\IPA\\Z-Scores 14 cell types.xlsx",sheet = "p05 pathway zscores",
#                 rowNames = F,colNames = T)
# colnames(IPA) <- IPA[1, ]  
# IPA <- IPA[-1, ]
# colnames(IPA) <- sub("_.*", "", colnames(IPA))
# for (i in 2:ncol(IPA)) {
#   IPA[[i]] <- as.numeric(ifelse(IPA[[i]] == "N/A", NA, IPA[[i]]))
# }
# str(IPA)
# IPA = reshape2::melt(IPA, id.vars = "Canonical Pathways")
# 
# IPA <- rename(IPA, variable = 'Cell_Type' , value =  'Z_Score')
# 
# matches <- IPA$`Canonical Pathways`[grep("Insulin Secretion Signaling Pathway|Epithelial Adherens Junction Signaling|Neuregulin Signaling|WNT/SHH Axonal Guidance Signaling Pathway", 
#                 IPA$`Canonical Pathways`)]
# 
# IPA = IPA %>% dplyr::filter(`Canonical Pathways` %in% matches)
# IPA = IPA %>% dplyr::filter( Cell_Type %in% c("Beta 1","Alpha 1","Delta 1", "Delta 2", "Gamma","Stellate 1","Stellate 2","Vascular cells","Proliferating cells"))
# 
# desired_order <- c("Beta 1", "Alpha 1", 
#                    "Delta 1", "Delta 2", "Gamma", 
#                    "Stellate 1", "Stellate 2", 
#                    "Vascular cells", "Proliferating cells")
# 
# # Reorder Cell_Type by setting it as a factor with the desired order
# IPA$Cell_Type <- factor(IPA$Cell_Type, levels = desired_order)
# 
# # Create the heatmap with ggplot2
# ggplot(IPA, aes(x = Cell_Type, y = `Canonical Pathways`, fill = Z_Score)) +
#   geom_tile() +
#   scale_fill_gradient2(name = "Z-Score",low = "blue", mid = "white", high = "red", midpoint = 0) +
#   theme_minimal() +
#   labs(title = "IPA",
#        x = "Cell types",
#        y = "Canonical Pathways") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
