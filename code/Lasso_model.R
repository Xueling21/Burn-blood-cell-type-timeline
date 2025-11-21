# Load necessary libraries
library(glmnet)
library(foreach)
library(doParallel)
library(writexl)

# Parallel computation setup
n_cores <- parallel::detectCores() - 10
cl <- makeCluster(n_cores)
registerDoParallel(cl)

####### Seq Burn #########
fraction <- read.table("/home/sda1/burn_JWJ/Burnsequencing/fraction_data.txt", header = TRUE, row.names = 1, sep = "\t")
exp <- read.table("/home/sda1/burn_JWJ/Burnsequencing/sequencingdata.txt", header = TRUE, row.names = 1, sep = "\t")

# Filter burn samples and log-transform expression data
burn_samples <- rownames(fraction)[grep("Bn", rownames(fraction))]
burn_fraction <- fraction[burn_samples, ]
exp <- log2(exp + 1)
burn_exp <- exp[, burn_samples]
burn_fraction <- burn_fraction[match(colnames(burn_exp), rownames(burn_fraction)), ]

# Filter genes with non-zero variance
gene_variance <- apply(burn_exp, 1, var)
burn_exp_select <- burn_exp[gene_variance > 0, ]

# Prepare phenotype data
r_1 <- burn_fraction[, 1:14]
r_1 <- r_1[, apply(r_1, 2, function(x) any(x != 0))]
burn_phen <- read.table("/home/sda1/burn_JWJ/Burnsequencing/burn_phe.txt", header = TRUE, row.names = 1, sep = "\t")
burn_phen$id <- paste0("fpkm_", rownames(burn_phen))
r_1$id <- rownames(r_1)
r_1 <- merge(burn_phen, r_1, by = "id")
rownames(r_1) <- r_1[, 1]
r_1 <- r_1[, 3:18]
r_1$gender <- as.factor(r_1$gender)
# Transpose expression data
burn_exp_select_t <- t(burn_exp_select)

# Filter genes with < 20% non-zero values
non_zero_count <- apply(burn_exp_select_t, 2, function(x) sum(x != 0))
non_zero_ratio <- non_zero_count / nrow(burn_exp_select_t)
genes_to_keep <- names(non_zero_ratio[non_zero_ratio > 0.2])
burn_exp_select_t <- burn_exp_select_t[, genes_to_keep]

# Bootstrap setup
n_bootstrap <- 100
X <- as.matrix(r_1)
Y <- burn_exp_select_t
genes <- colnames(Y)
factors <- colnames(X)

# Step 1: Cross-validation for lambda
global_lambda <- sapply(genes, function(g) {
  cv <- cv.glmnet(X, Y[, g], alpha = 1)
  cv$lambda.min
})

# Step 2: Bootstrap and Lasso regression
burn_coefs <- array(NA, dim = c(length(genes), length(factors), n_bootstrap))
dimnames(burn_coefs) <- list(genes, factors, NULL)

for (b in 1:n_bootstrap) {
  set.seed(b)
  idx <- sample(1:nrow(X), replace = TRUE)
  X_boot <- X[idx, ]
  Y_boot <- Y[idx, ]
  
  gene_results <- foreach(g = seq_along(genes), .combine = rbind, .packages = "glmnet") %dopar% {
    gene <- genes[g]
    y <- Y_boot[, gene]
    model <- glmnet(X_boot, y, alpha = 1, lambda = global_lambda[gene])
    coef_b <- coef(model)
    beta <- as.numeric(coef_b[-1])  # Exclude intercept
    names(beta) <- factors
    return(beta)
  }
  
  burn_coefs[, , b] <- gene_results
  cat("Finished bootstrap", b, "\n")
}

# Step 3: Calculate mean and SD of coefficients
final_results <- list()
for (f in factors) {
  coef_mat <- burn_coefs[, f, ]
  coef_mat[coef_mat < 0] <- 0
  coef_mean <- apply(coef_mat, 1, mean)
  coef_sd <- apply(coef_mat, 1, sd)
  
  final_results[[f]] <- data.frame(
    Gene = rownames(coef_mat),
    coef_mat,
    Mean = coef_mean,
    SD = coef_sd
  )
}


####### Seq Control #########
con_samples <- rownames(fraction)[grep("C", rownames(fraction))]
con_fraction <- fraction[con_samples, ]
con_exp <- exp[, con_samples]
con_fraction <- con_fraction[match(colnames(con_exp), rownames(con_fraction)), ]

# Filter genes with non-zero variance
gene_variance <- apply(con_exp, 1, var)
con_exp_select <- con_exp[gene_variance > 0, ]

# Prepare phenotype data for control samples
r_2 <- con_fraction[, 1:14]
r_2 <- r_2[, apply(r_2, 2, function(x) any(x != 0))]
con_phen <- read.table("/home/sda1/burn_JWJ/Burnsequencing/con_phe.txt", header = TRUE, row.names = 1, sep = "\t")
con_phen$id <- paste0("fpkm_", rownames(con_phen))
r_2$id <- rownames(r_2)
r_2 <- merge(con_phen, r_2, by = "id")
rownames(r_2) <- r_2[, 1]
r_2 <- r_2[, 3:18]
r_2$gender <- as.factor(r_2$gender)
# Transpose the expression data for control samples
con_exp_select_t <- t(con_exp_select)
# Filter genes with < 30% non-zero values
non_zero_count <- apply(con_exp_select_t, 2, function(x) sum(x != 0))
non_zero_ratio <- non_zero_count / nrow(con_exp_select_t)
genes_to_keep <- names(non_zero_ratio[non_zero_ratio > 0.3])
con_exp_select_t <- con_exp_select_t[, genes_to_keep]

# Bootstrap setup for control samples
X <- as.matrix(r_2)
Y <- con_exp_select_t
genes <- colnames(Y)
factors <- colnames(X)

# Step 1: Cross-validation for lambda (control)
global_lambda <- sapply(genes, function(g) {
  cv <- cv.glmnet(X, Y[, g], alpha = 1)
  cv$lambda.min
})

# Step 2: Bootstrap and Lasso regression (control)
con_coefs <- array(NA, dim = c(length(genes), length(factors), n_bootstrap))
dimnames(con_coefs) <- list(genes, factors, NULL)

for (b in 1:n_bootstrap) {
  set.seed(b)
  idx <- sample(1:nrow(X), replace = TRUE)
  X_boot <- X[idx, ]
  Y_boot <- Y[idx, ]
  
  gene_results <- foreach(g = seq_along(genes), .combine = rbind, .packages = "glmnet") %dopar% {
    gene <- genes[g]
    y <- Y_boot[, gene]
    model <- glmnet(X_boot, y, alpha = 1, lambda = global_lambda[gene])
    coef_b <- coef(model)
    beta <- as.numeric(coef_b[-1])  # Exclude intercept
    names(beta) <- factors
    return(beta)
  }
  
  con_coefs[, , b] <- gene_results
  cat("Finished bootstrap", b, "\n")
}

# Step 3: Calculate mean and SD of coefficients (control)
con_final_results <- list()
for (f in factors) {
  coef_mat <- con_coefs[, f, ]
  coef_mat[coef_mat < 0] <- 0
  coef_mean <- apply(coef_mat, 1, mean)
  coef_sd <- apply(coef_mat, 1, sd)
  
  con_final_results[[f]] <- data.frame(
    Gene = rownames(coef_mat),
    coef_mat,
    Mean = coef_mean,
    SD = coef_sd
  )
}

# Stop parallel processing for control bootstrap
stopCluster(cl)

####### Differential Expression Analysis #########
seq_diff_results <- list()
for (cell_type in names(con_final_results)) {
  burn_data <- final_results[[cell_type]]
  con_data <- con_final_results[[cell_type]]
  common_genes <- intersect(burn_data$Gene, con_data$Gene)
  burn_data <- burn_data[burn_data$Gene %in% common_genes, ]
  con_data <- con_data[con_data$Gene %in% common_genes, ]
  
  # Compute Z-scores and p-values
  delta_mean <- burn_data$Mean - con_data$Mean
  se_burn <- burn_data$SD / sqrt(n_bootstrap)
  se_con  <- con_data$SD  / sqrt(n_bootstrap)
  se_combined <- sqrt(se_burn^2 + se_con^2)
  z_scores <- delta_mean / se_combined
  p_values <- 2 * (1 - pnorm(abs(z_scores)))
  fdr <- p.adjust(p_values, method = "fdr")
  
  result_df <- data.frame(
    Gene = burn_data$Gene,
    BurnMean = burn_data$Mean,
    ConMean = con_data$Mean,
    Z = z_scores,
    p_value = p_values,
    FDR = fdr
  )
  
  seq_diff_results[[cell_type]] <- result_df
}

# Filter for significant genes (p-value < 0.05)
seq_significant_genes <- list()
for (cell_type in names(seq_diff_results)) {
  result_df <- seq_diff_results[[cell_type]]
  sig_df <- result_df[result_df$p_value < 0.05, ]
  seq_significant_genes[[cell_type]] <- sig_df
  cat(paste("Significant genes for", cell_type, ":", nrow(sig_df), "\n"))
}

####### Hypergeometric Test #########
hypergeometric_results <- list()
for (cell_type in names(seq_significant_genes)) {
  seq_genes <- seq_significant_genes[[cell_type]]$Gene
  gse_genes <- GSE19743_significant_genes[[cell_type]]$Gene
  intersection_genes <- intersect(seq_genes, gse_genes)
  p_value <- phyper(length(intersection_genes) - 1, length(gse_genes), total_genes - length(gse_genes), length(seq_genes), lower.tail = FALSE)
  hypergeometric_results[[cell_type]] <- list(
    intersection_count = length(intersection_genes),
    p_value = p_value
  )
  cat(paste("Hypergeometric test for", cell_type, ":", length(intersection_genes), "genes overlap, p-value =", p_value, "\n"))
}

