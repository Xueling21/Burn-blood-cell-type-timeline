###############################################################
#   MFuzz Time-Series Clustering and Visualization Pipeline
#   - Reads multi-sheet DEG results
#   - Extracts logFC matrices
#   - Performs MFuzz soft clustering
#   - Visualizes temporal gene expression trends
###############################################################


library(readxl)
library(dplyr)
library(Biobase)
library(Mfuzz)
library(ggplot2)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(enrichplot)
library(openxlsx)


##### 1. Read all sheets from the input Excel file####
file <- "high_941celltype_filtered_DEG.xlsx"
sheets <- excel_sheets(file)

# Each sheet is assumed to contain a DEG table for one cell type
data_list <- lapply(sheets, function(s){
  read_excel(file, sheet = s) |> as.data.frame()
})
names(data_list) <- sheets


##### 2. Extract logFC columns from each sheet and standardize format####
logfc_list <- lapply(names(data_list), function(s){
  df <- data_list[[s]]
  gene_col <- colnames(df)[1]
  
  # Identify logFC columns by prefix "logFC"
  logfc_cols <- grep("^logFC", colnames(df), value = TRUE)
  if (length(logfc_cols) == 0){
    warning(paste("Sheet", s, "contains no logFC columns."))
    return(NULL)
  }
  
  # Subset gene + logFC columns
  df_out <- df[, c(gene_col, logfc_cols), drop = FALSE]
  colnames(df_out)[1] <- "Gene"
  
  # Remove missing rows
  df_out <- df_out[complete.cases(df_out), ]
  return(df_out)
})

names(logfc_list) <- names(data_list)
logfc_list <- lapply(logfc_list, function(df){
  if (!is.null(df))
    colnames(df) <- gsub("^logFC_", "", colnames(df))
  return(df)
})

logfc_list <- lapply(logfc_list, function(df){
  if (!is.null(df) && ncol(df) > 1){
    rownames(df) <- df$Gene
    df <- df[, -1, drop = FALSE]
  }
  return(df)
})

##### 3. MFuzz Clustering (Example: Monocyte dataset)####
mat <- as.matrix(logfc_list$Monocytes)
mat_z <- t(scale(t(mat)))

# Construct ExpressionSet object required by MFuzz
eset <- ExpressionSet(assayData = mat_z)

# Filter low-quality or low-variance genes
gene.f <- filter.NA(eset, thres = 0.25)
gene.f <- filter.std(gene.f, min.std = 0)

# Standardize expression values for MFuzz analysis
gene.s <- standardise(gene.f)

# Estimate fuzzification parameter m
m <- mestimate(gene.s)

# Perform soft clustering
cl <- mfuzz(gene.s, c = 6, m = m)

# Retrieve cluster "core" genes above membership threshold
mem_thres <- 0.6
core_gene <- acore(gene.s, cl, min.acore = mem_thres)

##### 4. Visualization Settings (Example: Monocytes cell dataset)####
mat_raw <- as.matrix(logfc_list$Monocytes)
time_labels <- colnames(mat_raw)
xpos <- seq_along(time_labels)

# Plot style parameters
label_top <- 5
curve_lwd <- 0.7
mean_lwd  <- 3
bg_col <- rgb(0.2,0.2,0.2,0.25)
hi_cols <- c("#D55E00","#0072B2","#009E73",
             "#CC79A7","#E69F00","#56B4E9")
iqr_fill <- rgb(0.8,0.1,0.1,0.18)
# Space allocated on the right side of each panel for labels
x_pad <- 2.2
##### 5. Helper function: draw high-membership gene names####
draw_labels_in_zone <- function(labels, cols,
                                x0, x1, y0, y1,
                                ncol = 1,
                                line_len_frac = 0.22,
                                cex_text = 0.85,
                                lwd_line = 2,
                                box_fill = NA,
                                box_border = NA,
                                pad_frac = c(0.10, 0.12),
                                row_spacing = 1.10){
  
  if (length(labels) == 0) return(invisible(NULL))
  
  n <- length(labels)
  ncol <- min(ncol, n)
  nrow <- ceiling(n / ncol)
  
  # Optional background rectangle
  rect(x0, y0, x1, y1, col = box_fill, border = box_border)
  
  xr <- x1 - x0
  yr <- y1 - y0
  
  # Padding within the label panel
  xi <- x0 + xr * pad_frac[1]
  xa <- x1 - xr * pad_frac[1]
  y_top <- y1 - yr * pad_frac[2]
  y_bot <- y0 + yr * pad_frac[2]
  
  xw <- (xa - xi) / ncol
  yh <- (y_top - y_bot) / nrow * row_spacing
  line_len <- xr * line_len_frac
  
  # Draw label rows and colored marker lines
  idx <- 1
  for (c in 1:ncol){
    for (r in 1:nrow){
      if (idx > n) break
      
      yy <- y_top - (r - 1) * yh - yh/2
      xx <- xi + (c - 1) * xw
      
      segments(xx, yy, xx + line_len, yy, col = cols[idx], lwd = lwd_line)
      text(xx + line_len + strwidth("  ", cex = cex_text), yy,
           labels[idx], adj = c(0, 0.5),
           cex = cex_text, col = cols[idx])
      
      idx <- idx + 1
    }
  }
}


##### 6. Plot MFuzz clusters#####
pdf("Mon_MFuzz_finalLayout.pdf", width = 15, height = 9)
op <- par(mfrow = c(3, 2), mar = c(5.6, 5, 3.2, 3.6))

for (i in 1:6){
  # Extract core genes for cluster i
  genes_i_core <- core_gene[[i]]$NAME
  
  if (is.null(genes_i_core)){
    plot.new()
    title(main = sprintf("Cluster %d (N = 0)", i))
    next
  }
  
  # Subset expression matrix
  submat <- mat_raw[intersect(genes_i_core, rownames(mat_raw)), , drop = FALSE]
  
  if (nrow(submat) == 0){
    plot.new()
    title(main = sprintf("Cluster %d (empty)", i))
    next
  }
  
  # Select top membership genes for labeling
  mem_vec <- cl$membership[rownames(submat), i]
  ord <- order(mem_vec, decreasing = TRUE)
  lab_idx <- head(ord, label_top)
  lab_genes <- rownames(submat)[lab_idx]
  lab_cols <- hi_cols[seq_along(lab_genes)]
  
  # Initialize plot
  plot(xpos, submat[1, ], type = "n", xaxt = "n",
       xlab = "", ylab = "logFC",
       xlim = c(1, length(xpos) + x_pad),
       ylim = range(submat, na.rm = TRUE),
       main = sprintf("Cluster %d (core n = %d, m = %.2f)",
                      i, nrow(submat), mem_thres))
  
  axis(1, at = xpos, labels = time_labels, las = 2)
  abline(h = 0, lty = 2)
  
  # Interquartile range shading
  q25 <- apply(submat, 2, quantile, 0.25, na.rm = TRUE)
  q75 <- apply(submat, 2, quantile, 0.75, na.rm = TRUE)
  polygon(c(xpos, rev(xpos)), c(q25, rev(q75)),
          col = iqr_fill, border = NA)
  
  # Background gene curves
  others <- setdiff(rownames(submat), lab_genes)
  if (length(others) > 0){
    apply(submat[others, , drop = FALSE], 1,
          function(y) lines(xpos, y, col = bg_col, lwd = curve_lwd))
  }
  
  # Highlight top membership genes
  for (j in seq_along(lab_genes)){
    g <- lab_genes[j]
    lines(xpos, submat[g, ], col = lab_cols[j], lwd = curve_lwd + 0.2)
  }
  
  # Draw average trend line
  lines(xpos, colMeans(submat), col = "firebrick", lwd = mean_lwd)
  

  # Add labeled gene panel on the right
  usr <- par("usr")
  xr <- usr[2] - usr[1]
  yr <- usr[4] - usr[3]
  
  x0 <- length(xpos) + 0.1
  x1 <- length(xpos) + x_pad - 0.1
  y0 <- usr[3] + 0.1 * yr
  y1 <- usr[3] + 0.6 * yr
  
  par(xpd = NA)
  draw_labels_in_zone(lab_genes, lab_cols,
                      x0, x1, y0, y1,
                      ncol = 1,
                      cex_text = 1,
                      lwd_line = 1.5)
  par(xpd = FALSE)
}

par(op)
dev.off()