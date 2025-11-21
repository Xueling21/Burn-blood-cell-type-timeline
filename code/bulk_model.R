library(WGCNA)
library(dplyr)
library(clusterProfiler)  
library(org.Hs.eg.db)     
library(ggplot2)         
library(ggpubr)    
library(readxl)
library(caret)
library(pROC)
library(glmnet)
library(survival)
library(survminer)
library(DMwR)
library(randomForest)
####----1.GSE182616 WGCNA Analysis----
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
# Load the expression and clinical data
load("/home/jwj/skindata/input/GSE182616/GPL17077_GSE182616exp.RData")

# Mortality: Alive = 1, Death = 0
GSE182616_clinical$mortality <- ifelse(GSE182616_clinical$mortality == "Death", 0, 
                                       ifelse(GSE182616_clinical$mortality == "Alive", 1, NA))

# Transpose expression data and select high variance genes
datExpr0 = as.data.frame(t(GSE182616_cy5))
m.vars = apply(datExpr0, 2, var)
expro.upper = datExpr0[, which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
datExpr1 <- data.matrix(expro.upper)

# Perform hierarchical clustering to detect outliers
sampleTree = hclust(dist(datExpr1), method = "average")
abline(h = 150, col = "red")  # Add a cutoff line for outlier detection

# Select good samples
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
keepSamples = (clust == 1)
datExpr = datExpr1[keepSamples, ]

# Match samples to clinical data
femaleSamples = rownames(datExpr)
traitRows = match(femaleSamples, rownames(GSE182616_clinical))
datTraits = GSE182616_clinical[traitRows, ]
datTraits$age <- as.numeric(datTraits$age)
datTraits$gender <- as.numeric(datTraits$gender)
datTraits$mortality <- as.numeric(datTraits$mortality)
datTraits$time <- as.numeric(datTraits$time)

# Plot sample dendrogram with clinical traits
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# Soft threshold selection for network construction
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale-free topology fit index and mean connectivity
pdf("1Threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n")
abline(h = 0.85, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n")
dev.off()

# Construct gene modules using blockwiseModules
net = blockwiseModules(datExpr, power = 7, TOMType = "unsigned", minModuleSize = 30, 
                       reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, 
                       pamRespectsDendro = FALSE, saveTOMs = TRUE, verbose = 3)

# Plot the dendrogram of modules
mergedColors = labels2colors(net$colors)
pdf("2module.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

# Save module eigengenes and other results
MEs = net$MEs
save(MEs, file = "FemaleLiver-02-networkConstruction-auto.RData")

# Correlate module eigengenes with clinical traits
MEsWW = orderMEs(MEs)
modTraitCor = cor(MEsWW, datTraits, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)

# Plot the module-trait relationship heatmap
pdf("3Module-trait.pdf", width = 6, height = 6)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(datTraits), yLabels = names(MEsWW),
               colors = blueWhiteRed(50), textMatrix = paste(signif(modTraitCor, 2), "\n(", 
                                                             signif(modTraitP, 1), ")", sep = ""), setStdMargins = FALSE, cex.text = 0.5, 
               zlim = c(-1, 1), main = "Module-trait relationships")
dev.off()

# Identify hub genes in a specific module (green module)
module = "green"
column = match(module, modNames)
moduleGenes = moduleColors == module
green_hub <- names(net$colors)[which(moduleColors == module)]
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneTraitSignificance[moduleGenes, 1])
write.table(green_hub, file = "green_hubgenes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

####----2.Green Model genes GO enrichment----
ensem <- bitr(green_hub$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
# Perform GO enrichment analysis (ALL, MF, BP, CC)
go_ALL <- enrichGO(gene = ensem$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", 
                   ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

go_mf <- enrichGO(ensem$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "fdr", 
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = "ENTREZID", readable = TRUE)

go_bp <- enrichGO(ensem$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", 
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = "ENTREZID", readable = TRUE)

go_cc <- enrichGO(ensem$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "fdr", 
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = "ENTREZID", readable = TRUE)

# Store enrichment results
bp <- as.data.frame(go_bp)
mf <- go_mf@result
cc <- go_cc@result
write.csv(bp, file = "ego_result_BP.csv", row.names = T)
write.csv(cc, file = "ego_result_CC.csv", row.names = T)
write.csv(mf, file = "ego_result_MF.csv", row.names = T)

# Select top N pathways for each GO category
display_number = c(8, 8, 8)  
ego_result_BP <- as.data.frame(bp)[1:display_number[1], ]
ego_result_CC <- as.data.frame(cc)[1:display_number[2], ]
ego_result_MF <- as.data.frame(mf)[1:display_number[3], ]

# Combine results for plotting
go_enrich_df <- data.frame(
  ID = c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
  Description = c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
  GeneNumber = c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type = factor(c(rep("BP", display_number[1]), 
                  rep("CC", display_number[2]), 
                  rep("MF", display_number[3])), 
                levels = c("BP", "CC", "MF"))
)

# Create ordered factor for type (BP, CC, MF)
go_enrich_df$type_order <- factor(rev(as.integer(rownames(go_enrich_df))), labels = rev(go_enrich_df$Description))

# Set color palette for plot
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")  

# Create and save the GO enrichment plot
pdf("GO_Terms_of_green_hubgenes.pdf", width = 12, height = 8)
ggplot(data = go_enrich_df, aes(x = type_order, y = GeneNumber, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +  # Bar width
  scale_fill_manual(values = COLS) +  # Set colors
  coord_flip() +  # Flip the bars
  xlab("GO term") + 
  ylab("Gene Number") + 
  labs(title = "The Most Enriched GO Terms of green_hubgenes") +
  theme_bw() +  # Clean theme
  theme(
    axis.title = element_text(size = 14, face = "bold"),  
    axis.text = element_text(size = 14),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  
  )
dev.off()

####----3.RF Model ----
library(pROC)
library(glmnet)
library(DMwR)
library(randomForest)
# Prepare GSE182616 phenotype data (filtering by TBSA, age, and time)
GSE182616_phe <- GSE182616_phenotype %>%
  filter(tbsa > 20, age >= 18 & age <= 60, time <= (11 * 24))
# Filter GSE19743 based on criteria
GSE19743_phe <- GSE19743_phenotype %>%
  filter(tbsa > 20, age >= 18 & age <= 60, hours.post.injury <= (11 * 24))
# Subset expression data for matching phenotypes
a <- rownames(GSE19743_phe)
GSE19743 <- GSE19743_exp[a, ]
GSE19743_survival <- GSE19743_phe[, c("survival", "tbsa")]
# Combine expression data with survival information
GSE19743_gene <- colnames(GSE19743)
GSE19743_survival <- cbind(GSE19743_survival, GSE19743)
# Prepare mortality variable from GSE182616
mortality <- as.factor(GSE182616_phe$mortality)  
mm <- model.matrix(~0 + mortality)  
fit <- lmFit(GSE182616, mm)
# Define contrast for mortality comparison
contr <- makeContrasts(mortality1 - mortality0, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
# Filter differentially expressed genes (DEGs)
DEG_burn <- subset(top.table, P.Value < 0.05 & abs(logFC) > 0.5)
DEG_gene <- rownames(DEG_burn)
# Find common genes with the validation dataset
common_genes <- intersect(DEG_gene, GSE19743_gene)
# Prepare expression data for GSE182616 with common genes
GSE182616 <- as.data.frame(t(GSE182616))
GSE182616 <- GSE182616[, common_genes]
GSE182616_survial <- cbind(survival_data, GSE182616)

# Train-test split
set.seed(123)
train_indices <- createDataPartition(GSE182616_survial$mortality, p = 0.8, list = FALSE)
train_data <- GSE182616_survial[train_indices, ]
test_data <- GSE182616_survial[-train_indices, ]
train_data$mortality <- as.factor(train_data$mortality)
test_data$mortality <- as.factor(test_data$mortality)

####----3.1 LASSO Feature Selection ----
x <- as.matrix(train_data[, -1])  # Exclude label column
y <- as.factor(train_data$mortality)
# LASSO model with cross-validation
set.seed(42)
lasso_model <- cv.glmnet(x, y, alpha = 1, family = "binomial")
best_lambda <- lasso_model$lambda.min
lasso_coef <- coef(lasso_model, s = best_lambda)
# Select non-zero LASSO features
selected_features <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
selected_features <- selected_features[selected_features != "(Intercept)"]
# AUC filtering for LASSO-selected features
selected_auc <- c()
for (feature in selected_features) {
  roc_curve <- roc(y, x[, feature])
  auc_value <- auc(roc_curve)
  selected_auc <- c(selected_auc, auc_value)
  print(paste("feature:", feature, "AUC:", round(auc_value, 3)))
}
# Filter features based on AUC > 0.6
final_features <- selected_features[selected_auc > 0.6]

####----3.2 Random Forest with SMOTE ----
# Apply SMOTE to balance the data
train_data_balanced <- SMOTE(mortality ~ ., data = train_data, perc.over = 200, perc.under = 150)
# Train Random Forest model
rf_model <- randomForest(mortality ~ ., data = train_data_balanced, importance = TRUE)
# Plot Random Forest error vs trees
plot(rf_model, main = "Random Forest Error vs Trees", lwd = 2)
# Variable importance plot (Top 30 features)
varImpPlot(rf_model, n.var = min(30, nrow(rf_model$importance)), main = 'Top 30 Variable Importance')
# Cross-validation error curve
train.cv <- replicate(5, rfcv(train_data_balanced[-1], train_data_balanced$mortality, cv.fold = 10, step = 1.5), simplify = FALSE)
train.cv <- data.frame(sapply(train.cv, '[[', 'error.cv'))
train.cv$otus <- rownames(train.cv)
train.cv <- melt(train.cv, id = 'otus')
train.cv$otus <- as.numeric(as.character(train.cv$otus))
# Cross-validation error plot
p <- ggplot(train.cv, aes(otus, value)) +
  geom_smooth(se = FALSE, method = 'glm', formula = y ~ ns(x, 6)) +
  theme_minimal() +
  labs(title = 'Random Forest Cross-validation', x = 'Number of Features', y = 'CV Error') +
  geom_vline(xintercept = 30, linetype = "dashed")
print(p)

####----3.3 Feature Selection ----
# Select Top 30 important features based on Random Forest importance scores
importance_scores <- importance(rf_model)
selected_features <- rownames(importance_scores)[order(importance_scores[, "MeanDecreaseAccuracy"], decreasing = TRUE)][1:30]
# Final intersected feature set (LASSO + RF)
D <- intersect(selected_features, final_features)

####----3.4 Train Model with Selected Features ----
# Train Random Forest with selected features
train_data_balanced <- train_data_balanced[, c("mortality", D)]
model_rf4 <- randomForest(mortality ~ ., data = train_data_balanced, ntree = 50)
plot(model_rf4, main = "Random Forest Accuracy Curve", lwd = 2)

####----3.5 Test Model Performance ----
# Prediction and performance evaluation on the test set
predictions <- predict(model_rf4, test_data)
conf_matrix <- confusionMatrix(predictions, test_data$mortality, positive = "0")
print(conf_matrix)
# ROC and AUC for the test set
pre_prob <- predict(model_rf4, test_data, type = "prob")[, "0"]
pre_roc <- roc(test_data$mortality, pre_prob)
plot(pre_roc, 
     legacy.axes = TRUE,
     print.auc = TRUE, 
     grid = FALSE, 
     auc.polygon = TRUE,
     auc.polygon.col = "skyblue",
     max.auc.polygon = TRUE,
     main = "Test Set ROC Curve")
# Evaluation metrics: Accuracy, Precision, Recall, F1 Score
accuracy <- conf_matrix$overall["Accuracy"]
precision <- conf_matrix$byClass["Pos Pred Value"]
recall <- conf_matrix$byClass["Sensitivity"]
f1_score <- conf_matrix$byClass["F1"]

####----3.6 Validation Set ----
# Prepare validation data
vid <- GSE19743_survival
vid_features <- vid[, !colnames(vid) %in% c("survival")]
vid_features <- vid_features[, D]  # Select only selected features
vid_label <- vid$survival
# Predictions on validation set
vid_predictions <- predict(model_rf4, vid_features)
vid_confusion <- confusionMatrix(vid_predictions, vid_label, positive = "0")
print(vid_confusion)
# ROC on validation set
GSE19743_prob <- predict(model_rf4, vid, type = "prob")[, "0"]
roc_curve <- roc(vid$survival, GSE19743_prob)
plot(roc_curve, 
     legacy.axes = TRUE,
     print.auc = TRUE, 
     grid = FALSE, 
     auc.polygon = TRUE,
     auc.polygon.col = "skyblue",
     max.auc.polygon = TRUE,
     main = "Validation Set ROC Curve")

