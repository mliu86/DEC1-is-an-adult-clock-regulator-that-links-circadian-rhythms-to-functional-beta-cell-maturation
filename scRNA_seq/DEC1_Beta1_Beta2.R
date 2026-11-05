library(Seurat)
library(ggplot2)
library(gridExtra)
library(patchwork)
# install vioplot if needed
if (!requireNamespace("vioplot", quietly = TRUE)) {
  install.packages("vioplot", repos = "https://cloud.r-project.org", dependencies = TRUE)
}
library(vioplot)
library(Seurat)

load(file = "C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\RNA_ATAC_embed.rda")

# Make sure RNA assay is active
DefaultAssay(coembed) <- "RNA"

# Subset to Beta 1 and Beta 2
coembed_beta <- subset(coembed, subset = type %in% c("Beta 1", "Beta 2"))

# Fetch expression and metadata
df_beta <- FetchData(coembed_beta, vars = c("Bhlhe40", "type"))
df_beta$type <- factor(df_beta$type, levels = c("Beta 1", "Beta 2"))
df_beta$Bhlhe40_log2 <- log2(df_beta$Bhlhe40 + 1)

# Non-zero subset
df_beta_nonzero <- subset(df_beta, Bhlhe40 > 0)
df_beta_nonzero$type <- factor(df_beta_nonzero$type, levels = c("Beta 1", "Beta 2"))

# P-values FIRST so we can annotate once
p_all     <- wilcox.test(Bhlhe40_log2 ~ type, data = df_beta, exact = FALSE)$p.value
p_nonzero <- wilcox.test(Bhlhe40_log2 ~ type, data = df_beta_nonzero, exact = FALSE)$p.value
fmt_p <- function(p) if (is.na(p)) "p = NA" else if (p < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.3g", p)
lab_all     <- paste0("Wilcoxon: ", fmt_p(p_all))
lab_nonzero <- paste0("Wilcoxon (non-zero): ", fmt_p(p_nonzero))

# Colors & y-limit
pal  <- c("Beta 1" = "skyblue", "Beta 2" = "salmon")
y_max <- 0.5

# --- All cells (with optional dots) ---
p1 <- ggplot(df_beta, aes(x = type, y = Bhlhe40_log2, fill = type, color = type)) +
  geom_violin(trim = TRUE, linewidth = 0.4, alpha = 0.35) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "black", linewidth = 0.4, alpha = 0.95) +
  # uncomment next line to show all dots:
  # geom_jitter(size = 0.4, alpha = 0.25, width = 0.15, height = 0) +
  coord_cartesian(ylim = c(0, y_max), clip = "off") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.margin = margin(10, 20, 10, 10)) +
  labs(title = "Bhlhe40 log2(Expression) (all cells, y≤0.5)",
       x = "Cell Type", y = "log2(Expression + 1)") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = c("Beta 1" = "skyblue", "Beta 2" = "salmon")) +
  annotate("text", x = 1.5, y = y_max * 0.98, label = lab_all, size = 5)

# --- Non-zero only (with optional dots) ---
p2 <- ggplot(df_beta_nonzero, aes(x = type, y = Bhlhe40_log2, fill = type, color = type)) +
  geom_violin(trim = TRUE, linewidth = 0.4, alpha = 0.35) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "black", linewidth = 0.4, alpha = 0.95) +
  # uncomment next line to show all dots:
  # geom_jitter(size = 0.4, alpha = 0.25, width = 0.15, height = 0) +
  coord_cartesian(ylim = c(0, y_max), clip = "off") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none", plot.margin = margin(10, 20, 10, 10)) +
  labs(title = "Bhlhe40 log2(Expression) (non-zero only, y≤0.5)",
       x = "Cell Type", y = "log2(Expression + 1)") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = c("Beta 1" = "skyblue", "Beta 2" = "salmon")) +
  annotate("text", x = 1.5, y = y_max * 0.98, label = lab_nonzero, size = 5)

# Side-by-side
p1 | p2

# install.packages("vioplot")  # run once if needed
library(vioplot)

# Prepare vectors
b1_all <- df_beta$Bhlhe40_log2[df_beta$type == "Beta 1"]
b2_all <- df_beta$Bhlhe40_log2[df_beta$type == "Beta 2"]

nz <- df_beta$Bhlhe40 > 0
b1_nz <- df_beta$Bhlhe40_log2[df_beta$type == "Beta 1" & nz]
b2_nz <- df_beta$Bhlhe40_log2[df_beta$type == "Beta 2" & nz]


# P-values
fmt_p <- function(p) if (is.na(p)) "p = NA" else if (p < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.3g", p)
p_all     <- wilcox.test(b1_all, b2_all)$p.value
p_nonzero <- wilcox.test(b1_nz,  b2_nz)$p.value
lab_all     <- paste0("Wilcoxon: ", fmt_p(p_all))
lab_nonzero <- paste0("Wilcoxon (non-zero): ", fmt_p(p_nonzero))

# Y limit
y_max <- 0.5  # <- change to 2.5 etc. if you like

# Layout
op <- par(mfrow = c(1, 2), mar = c(4.1, 4.6, 3.0, 0.6), mgp = c(2.2, 0.7, 0))
y_max <- 0.5

# Panel 1
vioplot::vioplot(b1_all, b2_all, names = c("Beta 1","Beta 2"),
                 col = c("skyblue","salmon"), border = "grey25",
                 drawRect = FALSE, ylim = c(0, y_max), xlab = "Cell Type", ylab = "", main = "Bhlhe40 (all cells)")
mtext("log2(Expression + 1)", side = 2, line = 2.6, cex = 1.1)  # y-title for panel 1
boxplot(list(b1_all, b2_all), at = c(1,2), add = TRUE, outline = FALSE,
        col = NA, border = "black", boxwex = 0.22, axes = FALSE,
        whisklty = 1, whisklwd = 1.1, staplelty = 1, staplelwd = 1.1, medlwd = 2)

# Panel 2
vioplot::vioplot(b1_nz, b2_nz, names = c("Beta 1","Beta 2"),
                 col = c("skyblue","salmon"), border = "grey25",
                 drawRect = FALSE, ylim = c(0, y_max), xlab = "Cell Type", ylab = "", main = "Bhlhe40 (non-zero only)")
mtext("log2(Expression + 1)", side = 2, line = 2.6, cex = 1.1)  # y-title for panel 2
boxplot(list(b1_nz, b2_nz), at = c(1,2), add = TRUE, outline = FALSE,
        col = NA, border = "black", boxwex = 0.22, axes = FALSE,
        whisklty = 1, whisklwd = 1.1, staplelty = 1, staplelwd = 1.1, medlwd = 2)

par(op)




# vectors for vioplot
b1_all <- df_beta$Bhlhe40_log2[df_beta$type == "Beta 1"]
b2_all <- df_beta$Bhlhe40_log2[df_beta$type == "Beta 2"]
nz <- df_beta$Bhlhe40 > 0
b1_nz <- df_beta$Bhlhe40_log2[df_beta$type == "Beta 1" & nz]
b2_nz <- df_beta$Bhlhe40_log2[df_beta$type == "Beta 2" & nz]

## ======= Stats text =======
fmt_p <- function(p) if (is.na(p)) "p = NA" else if (p < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.3g", p)
p_all     <- wilcox.test(b1_all, b2_all)$p.value
p_nonzero <- wilcox.test(b1_nz,  b2_nz)$p.value
lab_all     <- paste0("Wilcoxon: ", fmt_p(p_all))
lab_nonzero <- paste0("Wilcoxon (non-zero): ", fmt_p(p_nonzero))

## ======= Plot settings =======
y_max <- 0.5
cols  <- c("skyblue", "salmon")

# Optional: save to file
# png("Bhlhe40_vioplot_all_vs_nonzero.png", width = 1400, height = 850, res = 140)

op <- par(mfrow=c(1,2), mar=c(4.1,6.0,3.2,0.8), mgp=c(2.3,0.7,0))


## ---- Panel 1: All cells ----
vioplot::vioplot(b1_all,b2_all, names=c("Beta 1","Beta 2"),
                 col=c("skyblue","salmon"), border="grey25",
                 drawRect=FALSE, ylim=c(0,y_max), xlab="Cell Type", ylab="", main="Bhlhe40 (all cells)")
mtext(expression(log[2]~"(Expression + 1)"), side=2, line=3.2, cex=1.0)


# overlay clean boxplots with whiskers (error bars)
boxplot(list(b1_all, b2_all), at = c(1, 2), add = TRUE, outline = FALSE,
        col = NA, border = "black", boxwex = 0.22, axes = FALSE,
        whisklty = 1, whisklwd = 1.1, staplelty = 1, staplelwd = 1.1, medlwd = 2)
# p-value label near top-left
text(x = 0.6, y = y_max * 0.975, labels = lab_all, adj = c(0, 1))

## ---- Panel 2: Non-zero only ----
vioplot::vioplot(b1_nz,b2_nz, names=c("Beta 1","Beta 2"),
                 col=c("skyblue","salmon"), border="grey25",
                 drawRect=FALSE, ylim=c(0,y_max), xlab="Cell Type", ylab="", main="Bhlhe40 (non-zero only)")
mtext(expression(log[2]~"(Expression + 1)"), side=2, line=3.2, cex=1.0)
boxplot(list(b1_nz, b2_nz), at = c(1, 2), add = TRUE, outline = FALSE,
        col = NA, border = "black", boxwex = 0.22, axes = FALSE,
        whisklty = 1, whisklwd = 1.1, staplelty = 1, staplelwd = 1.1, medlwd = 2)
text(x = 0.6, y = y_max * 0.975, labels = lab_nonzero, adj = c(0, 1))

par(op)
# dev.off()  # uncomment if you used png()








#### 10-28-2025
## ============================== ##
##  Libraries
## ============================== ##
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)    # log10_trans, log1p_trans, label_number
})

## ============================== ##
##  Load object and subset
## ============================== ##
load("C:\\Users\\maili\\Box\\Data Club\\Mai Liu\\Computational analysis\\sc_Mai_outputs\\20240418\\sc-ATAC\\RNA_ATAC_embed.rda")
DefaultAssay(coembed) <- "RNA"

# Beta only
coembed_beta <- subset(coembed, subset = type %in% c("Beta 1","Beta 2"))

## ============================== ##
##  Fetch raw expression + meta
## ============================== ##
df_beta <- FetchData(coembed_beta, vars = c("Bhlhe40","type"))
df_beta <- df_beta[!is.na(df_beta$Bhlhe40) & !is.na(df_beta$type), , drop = FALSE]
df_beta$type <- factor(df_beta$type, levels = c("Beta 1","Beta 2"))
df_beta_nonzero <- subset(df_beta, Bhlhe40 > 0)

## ============================== ##
##  Stats on TRUE values
## ============================== ##
fmt_p <- function(p) if (is.na(p)) "p = NA" else if (p < 2.2e-16) "p < 2.2e-16" else sprintf("p = %.3g", p)
lab_all     <- paste0("Wilcoxon: ",        fmt_p(wilcox.test(Bhlhe40 ~ type, data = df_beta, exact = FALSE)$p.value))
lab_nonzero <- paste0("Wilcoxon (non-zero): ", fmt_p(wilcox.test(Bhlhe40 ~ type, data = df_beta_nonzero, exact = FALSE)$p.value))

## ============================== ##
##  Choose SEMILOG mode
##    "log10"  -> true semilog (needs ε for zeros)
##    "log1p"  -> semilog-style with log1p (zeros OK)
## ============================== ##
semilog_mode <- "log10"   # <- change to "log1p" if you prefer

# Aesthetic colors
pal <- c("Beta 1" = "skyblue", "Beta 2" = "salmon")

## ============================== ##
##  Build plotting variables + scales
## ============================== ##
if (semilog_mode == "log10") {
  # ---- True semilog (log10 y). Replace exact zeros by a tiny epsilon for plotting only.
  pos   <- df_beta$Bhlhe40[df_beta$Bhlhe40 > 0]
  eps   <- if (length(pos)) max(1e-8, as.numeric(quantile(pos, 0.01, na.rm = TRUE)) * 0.3) else 1e-6
  df_beta$y_plot        <- ifelse(df_beta$Bhlhe40 <= 0, eps, df_beta$Bhlhe40)
  df_beta_nonzero$y_plot <- df_beta_nonzero$Bhlhe40
  
  y_trans  <- scales::log10_trans()
  # nice log10 breaks within data range
  min_b <- floor(log10(min(df_beta$y_plot[df_beta$y_plot > 0], na.rm = TRUE)))
  max_b <- ceiling(log10(quantile(df_beta$y_plot, 0.995, na.rm = TRUE)))
  y_breaks <- 10^(min_b:max_b)
  y_labels <- scales::label_number()(y_breaks)
  
  y_title <- "Expression (raw; log10 y)"
} else {
  # ---- Semilog via log1p: zeros allowed, no ε needed
  df_beta$y_plot         <- df_beta$Bhlhe40
  df_beta_nonzero$y_plot <- df_beta_nonzero$Bhlhe40
  
  y_trans  <- scales::log1p_trans()
  # breaks in original units (not on log1p scale): pick pretty positives incl. 0
  upper <- as.numeric(quantile(df_beta$y_plot, 0.995, na.rm = TRUE))
  base_vals <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
  y_breaks  <- base_vals[base_vals <= max(1, upper)]
  if (length(y_breaks) < 5) y_breaks <- pretty(c(0, upper), n = 5)
  y_labels  <- y_breaks
  
  y_title <- "Expression (raw; log1p y)"
}

## ============================== ##
##  ggplot violins on semi-log y
## ============================== ##
common_violin <- list(
  geom_violin(trim = TRUE, adjust = 1.4, linewidth = 0.4, alpha = 0.35, scale = "width"),
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "black", linewidth = 0.4),
  scale_fill_manual(values = pal),
  scale_color_manual(values = pal),
  theme_minimal(base_size = 12),
  theme(legend.position = "none")
)

p1 <- ggplot(df_beta, aes(x = type, y = y_plot, fill = type, color = type)) +
  common_violin +
  scale_y_continuous(
    trans = y_trans,
    breaks = y_breaks, labels = y_labels,
    expand = expansion(mult = c(0.06, 0.10))
  ) +
  labs(title = "Bhlhe40 (all cells, semilog y)", x = "Cell Type", y = y_title) +
  annotate("text", x = 1.5, y = max(y_breaks), label = lab_all, vjust = -0.2, size = 4.5)

p2 <- ggplot(df_beta_nonzero, aes(x = type, y = y_plot, fill = type, color = type)) +
  common_violin +
  scale_y_continuous(
    trans = y_trans,
    breaks = y_breaks, labels = y_labels,
    expand = expansion(mult = c(0.06, 0.10))
  ) +
  labs(title = "Bhlhe40 (non-zero, semilog y)", x = "Cell Type", y = y_title) +
  annotate("text", x = 1.5, y = max(y_breaks), label = lab_nonzero, vjust = -0.2, size = 4.5)

p1 | p2

