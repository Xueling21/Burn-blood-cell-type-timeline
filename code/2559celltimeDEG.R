
library(limma)
library(openxlsx)

####----1. Load expression/phenotype data----
load("/home/sda1/burn_JWJ/data/SCENIC/high_dec/high_dec_data.RData")
Phenotype <- read.table(
  "/home/sda1/burn_JWJ/data/SCENIC/GSE37069+GSE19743_phenotype.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)
phenodata767 <- Phenotype[, 2:5]

####----2. Define time groups based on hours.post.injury----
#    Time = 0 for controls
#    Time = 1 for 0-C3 days (<= 72 h)
#    Time = 2 for 4-C10 days (<= 240 h)
#    Time = 3 for 11-C50 days (<= 1200 h)
#    Time = 4 for 50-C180+ days (> 1200 h)
phenodata767$Time <- ifelse(
  phenodata767$type == "control", 0,
  ifelse(
    phenodata767$hours.post.injury <= 72, 1,
    ifelse(
      phenodata767$hours.post.injury <= 240, 2,
      ifelse(
        phenodata767$hours.post.injury <= 1200, 3,
        4
      )
    )
  )
)

table(phenodata767$Time)

####----3. Function: run limma for one cell type across time windows----
run_limma_by_time <- function(expr_matrix,
                              cell_type,
                              pheno,
                              time_col = "Time",
                              type_col = "type",
                              age_col = "age",
                              sex_col = "sex",
                              log_transform = TRUE,
                              out_prefix = getwd()) {
  message("Running limma for cell type: ", cell_type)
  
  # -------- 3.1 Match samples between expression and phenotype --------
  common_samples <- intersect(colnames(expr_matrix), rownames(pheno))
  expr_matrix <- expr_matrix[, common_samples]
  pheno <- pheno[common_samples, ]
  
  if (log_transform) {
    expr_matrix <- log2(expr_matrix + 1)
  }
  
  gene_time <- list()
  
  # -------- 3.2 Loop over time windows: compare Time j vs control (Time 0) --------
  for (j in 1:4) {
    message("  Processing time group: ", j)
    
    # Keep samples in time group j and controls (Time 0)
    sub_pheno <- subset(pheno, pheno[[time_col]] %in% c(0, j))
    sub_expr  <- expr_matrix[, rownames(sub_pheno)]
    
    # Ensure age is numeric
    sub_pheno[[age_col]] <- as.numeric(sub_pheno[[age_col]])
    
    # Design matrix: burn vs control, adjusting for age and sex
    design <- model.matrix(
      as.formula(paste0("~0 + ", type_col, " + ", age_col, " + ", sex_col)),
      data = sub_pheno
    )
    colnames(design)[1:2] <- c("burn", "control")
    
    # Contrast: burn - control
    contrast.matrix <- makeContrasts(burn - control, levels = design)
    
    fit  <- lmFit(sub_expr, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    # All genes with limma statistics
    deg_all <- topTable(fit2, coef = 1, number = Inf, adjust = "BH")
    
    # -------- 3.3 Compute group means (burn/control) --------
    group_list <- sub_pheno[[type_col]]
    mean_df <- aggregate(t(sub_expr), by = list(group_list), FUN = mean)
    mean_df <- t(mean_df)[-1, ] 
    colnames(mean_df) <- paste0(levels(factor(group_list)), "_mean")
    mean_df <- data.frame(ID = rownames(sub_expr), mean_df)
    
    # Merge limma stats with group means
    deg_all <- data.frame(ID = rownames(deg_all), deg_all)
    deg_all <- merge(deg_all, mean_df, by = "ID", all.x = TRUE)
    rownames(deg_all) <- deg_all$ID
    deg_all$ID <- NULL
    
    # Significant DEGs: adj.P.Val < 0.05 & |logFC| > 1
    deg_sig <- subset(deg_all, adj.P.Val < 0.05 & abs(logFC) > 1)
    
    gene_time[[j]] <- list(
      limma    = deg_all,
      limmaDEG = deg_sig
    )
  }
  
  names(gene_time) <- c("0-3d", "4-10d", "11-50d", "50-180+d")
  
  # -------- 3.4 Save all results to an Excel workbook --------
  wb <- createWorkbook()
  
  for (nm in names(gene_time)) {
    # All genes
    limma_df <- gene_time[[nm]]$limma
    if (!is.null(limma_df) && nrow(limma_df) > 0) {
      limma_df <- data.frame(
        Gene = rownames(limma_df),
        limma_df,
        row.names = NULL
      )
      addWorksheet(wb, paste0(nm, "_all"))
      writeData(wb, sheet = paste0(nm, "_all"), limma_df)
    }
    
    # Significant DEGs
    deg_df <- gene_time[[nm]]$limmaDEG
    if (!is.null(deg_df) && nrow(deg_df) > 0) {
      deg_df <- data.frame(
        Gene = rownames(deg_df),
        deg_df,
        row.names = NULL
      )
      addWorksheet(wb, paste0(nm, "_DEG(logFC>1)"))
      writeData(wb, sheet = paste0(nm, "_DEG(logFC>1)"), deg_df)
    }
  }
  
  file_out <- file.path(out_prefix,
                        paste0(cell_type, "_limma_and_DEG_by_Time.xlsx"))
  saveWorkbook(wb, file_out, overwrite = TRUE)
  
  message("Finished limma for ", cell_type,
          ". Results saved to: ", file_out)
  
  return(gene_time)
}


####----4. Run limma for each cell type----
B_result <- run_limma_by_time(
  expr_matrix = Bcells,
  cell_type   = "Bcells",
  pheno       = phenodata767
)

T_result <- run_limma_by_time(
  expr_matrix = Tcells,
  cell_type   = "Tcells",
  pheno       = phenodata767
)

EOS_result <- run_limma_by_time(
  expr_matrix = Eosinophils,
  cell_type   = "Eosinophils",
  pheno       = phenodata767
)

Mon_result <- run_limma_by_time(
  expr_matrix = Monocytes,
  cell_type   = "Monocytes",
  pheno       = phenodata767
)

Neu_result <- run_limma_by_time(
  expr_matrix = Neutrophils,
  cell_type   = "Neutrophils",
  pheno       = phenodata767
)

NK_result <- run_limma_by_time(
  expr_matrix = NKcells,
  cell_type   = "NKcells",
  pheno       = phenodata767
)


####----5. Summarize logFC and p-values across all time windows----
make_time_table <- function(res,
                            cell_type = "Bcells",
                            out_file = paste0(cell_type, "_logFC_P_byTime.xlsx")) {
  time_names <- names(res)   # e.g. c("0-3d","4-10d","11-50d","50-180+d")
  
  # Extract logFC + P.Value + adj.P.Val from each time window
  gene_tables <- lapply(time_names, function(tn) {
    df <- res[[tn]]$limma          # full limma result for that time window
    df_sub <- df[, c("logFC", "P.Value", "adj.P.Val")]
    
    # Rename columns to indicate time window
    colnames(df_sub) <- c(
      paste0("logFC_", tn),
      paste0("P_", tn),
      paste0("FDR_", tn)
    )
    
    data.frame(
      Gene = rownames(df),
      df_sub,
      row.names = NULL
    )
  })
  
  # Merge all time windows by gene
  merged <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE),
                   gene_tables)
  
  # Write to Excel
  write.xlsx(merged, file = out_file, rowNames = FALSE)
  message("Summary table saved to: ", out_file)
  return(merged)
}

B_time_table   <- make_time_table(B_result,   "Bcells")
T_time_table   <- make_time_table(T_result,   "Tcells")
EOS_time_table <- make_time_table(EOS_result, "Eosinophils")
Mon_time_table <- make_time_table(Mon_result, "Monocytes")
Neu_time_table <- make_time_table(Neu_result, "Neutrophils")
NK_time_table  <- make_time_table(NK_result,  "NKcells")
