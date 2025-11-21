library('limma')       
library('readxl')      
library('openxlsx')    
setwd("/home/sda1/burn_JWJ/文章/fxy/GSE19743GSE3706920211217/bulk_time_DEGs")
# Load phenotype data
GSE_1937<-read.xlsx("GSE_1937.xlsx", rowNames = TRUE, colNames = TRUE)
phenodata767 <- read.xlsx("phenodata767.xlsx", rowNames = TRUE, colNames = TRUE)
# Create grouped time variable based on time after burn
phenodata767$Time <- ifelse(phenodata767$type == "control", 0,
                            ifelse(phenodata767$time < 12, 1,
                                   ifelse(phenodata767$time < 24, 2,
                                          ifelse(phenodata767$time < 72, 3,
                                                 ifelse(phenodata767$time < 168, 4,
                                                        ifelse(phenodata767$time < 240, 5,
                                                               ifelse(phenodata767$time < 720, 6,
                                                                      ifelse(phenodata767$time < 1200, 7,
                                                                             ifelse(phenodata767$time < 4320, 8, 9)))))))))

#################### limma Analysis ####################
gene_time <- list()
time <- list()
coef_value <- list()
for (j in 1:9) {
  phenodata <- subset(phenodata767, phenodata767$Time == j | phenodata767$Time == 0)
  expdat <-GSE_1937 
  expdat <- expdat[,rownames(phenodata) , drop = FALSE]
  boxplot(expdat)
  # Extract phenotype data
  age <- phenodata$age
  group_list <- phenodata[, 1]
  phenodata <- apply(phenodata, 2, as.factor)
  type <- phenodata[, 1]
  gender <- phenodata[, 3]
  # Build design matrix
  design <- model.matrix(~0 + type + gender + age)
  colnames(design)[1:2] <- c("burn", "control")
  # Build contrast: burn vs control
  matrix <- makeContrasts(burn - control, levels = design)
  # Run limma differential expression pipeline
  fit <- lmFit(expdat, design)
  coef_value[[j]] <- coef(fit)
  fit2 <- contrasts.fit(fit, matrix)
  fit3 <- eBayes(fit2)
  output <- topTable(fit3, coef = 1, n = Inf, adjust = "BH")
  DEG <- na.omit(output)
  # Compute mean and standard deviation per group
  meanmean <- aggregate(data.frame(t(expdat)), list(group_list), mean)
  meanmean <- t(meanmean)
  mean_names <- meanmean[1, ]
  meanmean <- data.frame(apply(meanmean, 2, as.numeric))  
  var <- aggregate(data.frame(t(expdat)), list(group_list), sd)
  var <- t(var)
  var <- data.frame(apply(var, 2, as.numeric))  
  me_var <- cbind(meanmean, var)
  mean_names <- c("burn_mean", "control_mean", "burn_var", "control_var")
  colnames(me_var) <- mean_names
  me_var <- me_var[-1, ]
  rownames(me_var) <- rownames(expdat)
  me_var <- data.frame(ID = rownames(me_var), me_var)
  DEG <- data.frame(ID = rownames(DEG), DEG)
  DEG <- merge(DEG, me_var, by = "ID")
  rownames(DEG) <- DEG[, 1]
  DEG <- DEG[, -1]
  # Filter DEGs with adjusted p-value < 0.05 and |logFC| > 1
  M <- subset(DEG, DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > 1)
  # Save results for each time point
  gene_time[[j]] <- list(limma = DEG, limmaDEG = M)
}

# Set readable names for time points
names(gene_time) <- c("12h", "1d", "2_3d", "4_7d", "8_10d", "11_30d", "31_50d", "51_180d", "181d")

# Load external DEG data (optional)
load("/home/sda1/burn_JWJ/文章/fxy/bulk_time_DEGs.RData")
#################### Collect all time point DEGs ####################
timeDEGs <- NULL
for (i in 1:length(gene_time)) {
  limmaDEG <- gene_time[[i]][[1]]
  selected_genes <- subset(limmaDEG, abs(limmaDEG$logFC) > 1.5 & limmaDEG$adj.P.Val < 0.05)
  timeDEGs <- union(timeDEGs, rownames(selected_genes))
}

#################### Build logFC matrix for all DEGs ####################
DEG_logFC <- gene_time[[1]][[1]]$logFC
for (i in 2:9) {
  DEG_logFC <- cbind(DEG_logFC, gene_time[[i]][[1]]$logFC)
}
colnames(DEG_logFC) <- names(gene_time)
rownames(DEG_logFC) <- rownames(gene_time[[1]][[1]])
genes_dSig <- DEG_logFC[timeDEGs, ]
boxplot(genes_dSig[, 1:9])
genes_dSig <- data.frame(t(genes_dSig))  # Transpose for time-series clustering

#################### Mfuzz Clustering ####################
gene_tpm <- data.matrix(genes_dSig)
temp <- data.frame(gene_tpm)
DEGs_exp_averp <- genes_dSig
eset <- new("ExpressionSet", exprs = DEGs_exp_averp)
# Filter and impute missing data
gene.r <- filter.NA(eset, thres = 0.25)
gene.f <- fill.NA(gene.r, mode = "knn")
tmp <- filter.std(gene.f, min.std = 0)
gene.s <- standardise(tmp)

# Load Mfuzz plotting script (custom version)
source('mfuzz.plot2.R')
# Estimate fuzzifier (m) and set number of clusters (c)
c <- 9
m <- mestimate(gene.s)
# Run Mfuzz clustering
cl <- mfuzz(gene.s, c = c, m = m)
# Plot clustered expression trends
mfuzz.plot(gene.s, cl, mfrow = c(3, 3), time.labels = colnames(DEGs_exp_averp))
# Cluster sizes and membership matrix
cl$size
membership <- cl$membership
head(membership)
# Extract highly representative genes in clusters
tmp <- acore(gene.s, cl, min.acore = 0.5)
SIGs_id2 <- rbind(tmp[[2]])
SDGs_id2 <- rbind(tmp[[4]], tmp[[5]], tmp[[6]])
# Visualize again (in new window)
mfuzz.plot(gene.s, cl, mfrow = c(3, 3), new.window = TRUE)
