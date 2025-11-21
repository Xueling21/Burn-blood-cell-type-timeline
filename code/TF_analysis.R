
# TF network analysis per cell type using limma + GENIE3

## -------------------- Libraries --------------------
suppressPackageStartupMessages({
  library(ggraph)
  library(ggplot2)
  library(doRNG)
  library(GENIE3)
  library(foreach)
  library(dplyr)
  library(readr)
  library(igraph)
  library(writexl)         
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(limma)
  library(grid)             
})

## -------------------- Read & Prepare Expression Matrices --------------------
# Folders containing *_Window24.txt files
folders <- c(
  "/home/sda1/burn_JWJ/data/SCENIC/high_dec/g1",
  "/home/sda1/burn_JWJ/data/SCENIC/high_dec/g2",
  "/home/sda1/burn_JWJ/data/SCENIC/high_dec/g3"
)

cell_types <- c("Bcells","Eosinophils","Monocytes","Neutrophils","NKcells","Tcells")
cell_type_data <- list()
for (cell_type in cell_types) {
  merged_data <- data.frame()
  for (folder in folders) {
    files <- list.files(folder, full.names = TRUE)
    matching_files <- files[grep(paste0(cell_type, "_Window24.txt$"), files)]
    for (file in matching_files) {
      dat <- read.csv(file, row.names = 1, header = TRUE, sep = "\t")
      merged_data <- rbind(merged_data, dat)  # row-bind replicates across folders
    }
  }
  cell_type_data[[cell_type]] <- merged_data
}

Bcells       <- cell_type_data$Bcells
Tcells       <- cell_type_data$Tcells
Eosinophils  <- cell_type_data$Eosinophils
Monocytes    <- cell_type_data$Monocytes
Neutrophils  <- cell_type_data$Neutrophils
NKcells      <- cell_type_data$NKcells

## -------------------- Phenotype & Design Matrix --------------------
Phenotype <- read.table("/home/sda1/burn_JWJ/data/SCENIC/GSE37069+GSE19743_phenotype.txt",
                        header = TRUE, row.names = 1, sep = "\t")
# Cast types; enforce level order to ensure contrast "burn - control" works reliably
Phenotype$age  <- as.numeric(Phenotype$age)
Phenotype$sex  <- factor(Phenotype$sex)
Phenotype$type <- factor(Phenotype$type, levels = c("control","burn"))
Phenotype <- Phenotype[, 2:5] 
age    <- Phenotype$age
type   <- Phenotype$type
gender <- Phenotype$sex
mm <- model.matrix(~0 + type + gender + age)

## -------------------- load TF list --------------------
# TRRUST human TFs (first column); used to restrict regulators when desired
tf_list <- read.table("/home/sda1/burn_JWJ/data/SCENIC/trrust_rawdata.human.tsv",
                      sep = "\t", header = FALSE)$V1
tf_list <- unique(tf_list)


## ============Neutrophils=================================================
Neutrophils <- log2(Neutrophils + 1)

fit  <- lmFit(Neutrophils, mm)
contr <- makeContrasts(typeburn - typecontrol, levels = colnames(coef(fit)))
tmp  <- contrasts.fit(fit, contr)
tmp  <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

DEG_NEU <- subset(top.table, P.Value < 0.05 & abs(logFC) > 0)
DEG_NEU$gene <- rownames(DEG_NEU)
NEU_genes    <- rownames(DEG_NEU)

# Burn IDs
burn_id <- rownames(Phenotype[Phenotype$type == "burn", ])
expr_burn_NEU <- Neutrophils[, burn_id, drop = FALSE]

# Regulators available in this matrix
available_tfs_NEU <- intersect(tf_list, rownames(expr_burn_NEU))
expr_burn_matrix  <- as.matrix(expr_burn_NEU)

set.seed(123)
weightMatrix_burn <- GENIE3(expr_burn_matrix,
                            regulators   = available_tfs_NEU,
                            nCores       = 8,
                            returnMatrix = TRUE)

# Keep top 1% edges
threshold_burn <- quantile(weightMatrix_burn, probs = 0.99)
links_burn     <- getLinkList(weightMatrix_burn, threshold = threshold_burn)

# TF activity summary
tf_activity <- links_burn %>%
  group_by(regulatoryGene) %>%
  summarise(n_targets = n(),
            max_weight = max(weight),
            mean_weight = mean(weight),
            .groups = "drop") %>%
  arrange(desc(n_targets))

# Build graph (all TFs found)
top_tfs <- tf_activity$regulatoryGene
sub_links <- links_burn[links_burn$regulatoryGene %in% top_tfs, ]
g <- graph_from_data_frame(sub_links, directed = TRUE)
V(g)$degree <- degree(g, mode = "all")
V(g)$type   <- ifelse(V(g)$name %in% top_tfs, "Top TFs", "Target genes")

set.seed(123)
ggraph(g, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = type), alpha = 0.7) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
  scale_color_manual(values = c("Top TFs" = "#FFA07A", "Target genes" = "#87CEFA")) +
  scale_size_continuous(range = c(4, 12)) +
  theme_void() +
  labs(title = "Neutrophil TF Regulatory Network in Burns", color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10)) +
  guides(size = "none",
         color = guide_legend(title.theme = element_text(face = "bold", size = 14),
                              label.theme = element_text(face = "bold", size = 12)))

save(DEG_NEU, weightMatrix_burn, file = "NEU_TF.RData")

# Per-TF GO enrichment (targets by TF)
tf_target_genes_list <- sub_links %>%
  group_by(regulatoryGene) %>%
  summarise(targetGenes = list(unique(targetGene)), .groups = "drop")

run_enrichment_for_tf <- function(tf_name, target_symbols) {
  # Map SYMBOL -> ENTREZ
  target_symbols <- as.character(target_symbols)
  valid_symbols  <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  filtered_syms  <- intersect(target_symbols, valid_symbols)
  if (length(filtered_syms) < 3) return(NULL)
  
  symbol2entrez <- tryCatch(
    bitr(filtered_syms, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (is.null(symbol2entrez) || nrow(symbol2entrez) == 0) return(NULL)
  
  entrez_ids <- symbol2entrez$ENTREZID[!is.na(symbol2entrez$ENTREZID)]
  if (length(entrez_ids) < 3) return(NULL)
  
  enr <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(enr) && nrow(enr@result) > 0) {
    df <- enr@result; df$TF <- tf_name; return(df)
  }
  NULL
}

all_enrichment_results <- lapply(seq_len(nrow(tf_target_genes_list)), function(i) {
  tf <- tf_target_genes_list$regulatoryGene[i]
  tg <- tf_target_genes_list$targetGenes[[i]]
  run_enrichment_for_tf(tf, tg)
})
names(all_enrichment_results) <- tf_target_genes_list$regulatoryGene
go_combined <- do.call(rbind, all_enrichment_results)
write.csv(go_combined, "NEU_TF_GO_enrichment_results.csv", row.names = FALSE)

# Highlight edges to up-/down-regulated DE genes
up_genes   <- DEG_NEU$gene[DEG_NEU$logFC > 0]
down_genes <- DEG_NEU$gene[DEG_NEU$logFC < 0]
edges <- links_burn[links_burn$targetGene %in% c(up_genes, down_genes), ]

all_nodes <- unique(c(edges$regulatoryGene, edges$targetGene))
vertices  <- data.frame(name = all_nodes, stringsAsFactors = FALSE)
vertices$node_class <- dplyr::case_when(
  vertices$name %in% edges$regulatoryGene ~ "TF",
  vertices$name %in% up_genes             ~ "Up",
  vertices$name %in% down_genes           ~ "Down",
  TRUE ~ NA_character_
)
vertices$node_class <- factor(vertices$node_class, levels = c("TF","Up","Down"))

g2 <- graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
V(g2)$degree <- degree(g2, mode = "all")

set.seed(123)
ggraph(g2, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = node_class), alpha = 0.85) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.8, fontface = "bold",
                 max.overlaps = 100, point.padding = unit(0.4, "lines")) +
  scale_color_manual(name = "Node Type",
                     values = c("TF"="#FFA07A","Up"="red","Down"="darkblue"),
                     breaks = c("TF","Up","Down"),
                     labels = c("TF","Up","Down")) +
  scale_size_continuous(range = c(3.5, 10)) +
  theme_void() +
  labs(title = "Neutrophil TF Regulatory Network in Burns") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 12, 10, 12)) +
  guides(size = "none")

## =========================T cells==========================================
Tcells <- log2(Tcells + 1)
fit  <- lmFit(Tcells, mm)
contr <- makeContrasts(typeburn - typecontrol, levels = colnames(coef(fit)))
tmp  <- contrasts.fit(fit, contr)
tmp  <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

DEG_T <- subset(top.table, P.Value < 0.05 & abs(logFC) > 0)
DEG_T$gene <- rownames(DEG_T)
T_genes    <- rownames(DEG_T)

burn_id <- rownames(Phenotype[Phenotype$type == "burn", ])
expr_burn_T <- Tcells[, burn_id, drop = FALSE]

available_tfs_T <- intersect(tf_list, T_genes)
expr_burn_matrix <- as.matrix(expr_burn_T)

set.seed(123)
weightMatrix_burn <- GENIE3(expr_burn_matrix,
                            regulators   = available_tfs_T,
                            nCores       = 8,
                            returnMatrix = TRUE)

threshold_burn <- quantile(weightMatrix_burn, probs = 0.95)
links_burn <- getLinkList(weightMatrix_burn, threshold = threshold_burn)

tf_activity <- links_burn %>%
  group_by(regulatoryGene) %>%
  summarise(n_targets = n(),
            max_weight = max(weight),
            mean_weight = mean(weight),
            .groups = "drop") %>%
  arrange(desc(n_targets))

top_tfs <- tf_activity$regulatoryGene
sub_links <- links_burn[links_burn$regulatoryGene %in% top_tfs, ]
g <- graph_from_data_frame(sub_links, directed = TRUE)
V(g)$degree <- degree(g, mode = "all")
V(g)$type   <- ifelse(V(g)$name %in% top_tfs, "Top TFs", "Target genes")

set.seed(123)
ggraph(g, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = type), alpha = 0.7) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
  scale_color_manual(values = c("Top TFs" = "#FFA07A", "Target genes" = "#87CEFA")) +
  scale_size_continuous(range = c(4, 12)) +
  theme_void() +
  labs(title = "T cells TF Regulatory Network in Burns", color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10)) +
  guides(size = "none",
         color = guide_legend(title.theme = element_text(face = "bold", size = 14),
                              label.theme = element_text(face = "bold", size = 12)))

save(DEG_T, weightMatrix_burn, file = "T_TF.RData")

tf_target_genes_list <- sub_links %>%
  group_by(regulatoryGene) %>%
  summarise(targetGenes = list(unique(targetGene)), .groups = "drop")

tf_target_genes_export <- tf_target_genes_list %>%
  mutate(targetGenes = sapply(targetGenes, function(x) paste(x, collapse = ", ")))
write_xlsx(tf_target_genes_export, path = "T_tf_target.xlsx")
write_xlsx(sub_links, path = "T_weight.xlsx")

# Reuse run_enrichment_for_tf() defined above
all_enrichment_results <- lapply(seq_len(nrow(tf_target_genes_list)), function(i) {
  tf <- tf_target_genes_list$regulatoryGene[i]
  tg <- tf_target_genes_list$targetGenes[[i]]
  run_enrichment_for_tf(tf, tg)
})
names(all_enrichment_results) <- tf_target_genes_list$regulatoryGene
go_combined <- do.call(rbind, all_enrichment_results)
write.csv(go_combined, "T_TF_GO_enrichment_results.csv", row.names = FALSE)

# Highlight up/down targets
up_genes   <- DEG_T$gene[DEG_T$logFC > 0]
down_genes <- DEG_T$gene[DEG_T$logFC < 0]
edges <- links_burn[links_burn$targetGene %in% c(up_genes, down_genes), ]

all_nodes <- unique(c(edges$regulatoryGene, edges$targetGene))
vertices  <- data.frame(name = all_nodes, stringsAsFactors = FALSE)
vertices$node_class <- dplyr::case_when(
  vertices$name %in% edges$regulatoryGene ~ "TF",
  vertices$name %in% up_genes             ~ "Up",
  vertices$name %in% down_genes           ~ "Down",
  TRUE ~ NA_character_
)
vertices$node_class <- factor(vertices$node_class, levels = c("TF","Up","Down"))

g2 <- graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
V(g2)$degree <- degree(g2, mode = "all")

set.seed(123)
ggraph(g2, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = node_class), alpha = 0.85) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.8, fontface = "bold",
                 max.overlaps = 100, point.padding = unit(0.4, "lines")) +
  scale_color_manual(name = "Node Type",
                     values = c("TF"="#FFA07A","Up"="red","Down"="darkblue"),
                     breaks = c("TF","Up","Down")) +
  scale_size_continuous(range = c(3.5, 10)) +
  theme_void() +
  labs(title = "T cells TF Regulatory Network in Burns") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 12, 10, 12)) +
  guides(size = "none")


## =======================B cells============================================
Bcells <- log2(Bcells + 1)
fit  <- lmFit(Bcells, mm)
contr <- makeContrasts(typeburn - typecontrol, levels = colnames(coef(fit)))
tmp  <- contrasts.fit(fit, contr)
tmp  <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

DEG_B <- subset(top.table, P.Value < 0.05 & abs(logFC) > 0)
DEG_B$gene <- rownames(DEG_B)
B_genes    <- rownames(DEG_B)

burn_id <- rownames(Phenotype[Phenotype$type == "burn", ])
expr_burn_B <- Bcells[, burn_id, drop = FALSE]

available_tfs_B <- intersect(tf_list, B_genes)
expr_burn_matrix <- as.matrix(expr_burn_B)

set.seed(123)
weightMatrix_burn <- GENIE3(expr_burn_matrix,
                            regulators   = available_tfs_B,
                            nCores       = 8,
                            returnMatrix = TRUE)

threshold_burn <- quantile(weightMatrix_burn, probs = 0.95)
links_burn <- getLinkList(weightMatrix_burn, threshold = threshold_burn)

tf_activity <- links_burn %>%
  group_by(regulatoryGene) %>%
  summarise(n_targets = n(),
            max_weight = max(weight),
            mean_weight = mean(weight),
            .groups = "drop") %>%
  arrange(desc(n_targets))

top_tfs <- tf_activity$regulatoryGene
sub_links <- links_burn[links_burn$regulatoryGene %in% top_tfs, ]
g <- graph_from_data_frame(sub_links, directed = TRUE)
V(g)$degree <- degree(g, mode = "all")
V(g)$type   <- ifelse(V(g)$name %in% top_tfs, "Top TFs", "Target genes")

set.seed(123)
ggraph(g, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = type), alpha = 0.7) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
  scale_color_manual(values = c("Top TFs" = "#FFA07A", "Target genes" = "#87CEFA")) +
  scale_size_continuous(range = c(4, 12)) +
  theme_void() +
  labs(title = "B cells TF Regulatory Network in Burns", color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10)) +
  guides(size = "none",
         color = guide_legend(title.theme = element_text(face = "bold", size = 14),
                              label.theme = element_text(face = "bold", size = 12)))

save(DEG_B, weightMatrix_burn, file = "B_T.RData")

tf_target_genes_list <- sub_links %>%
  group_by(regulatoryGene) %>%
  summarise(targetGenes = list(unique(targetGene)), .groups = "drop")

tf_target_genes_export <- tf_target_genes_list %>%
  mutate(targetGenes = sapply(targetGenes, function(x) paste(x, collapse = ", ")))
write_xlsx(tf_target_genes_export, path = "B_tf_target.xlsx")
write_xlsx(sub_links, path = "B_weight.xlsx")

# Reuse enrichment
all_enrichment_results <- lapply(seq_len(nrow(tf_target_genes_list)), function(i) {
  tf <- tf_target_genes_list$regulatoryGene[i]
  tg <- tf_target_genes_list$targetGenes[[i]]
  run_enrichment_for_tf(tf, tg)
})
names(all_enrichment_results) <- tf_target_genes_list$regulatoryGene
go_combined <- do.call(rbind, all_enrichment_results)
write.csv(go_combined, "B_TF_GO_enrichment_results.csv", row.names = FALSE)

# Highlight up/down targets
up_genes   <- DEG_B$gene[DEG_B$logFC > 0]
down_genes <- DEG_B$gene[DEG_B$logFC < 0]
edges <- links_burn[links_burn$targetGene %in% c(up_genes, down_genes), ]

all_nodes <- unique(c(edges$regulatoryGene, edges$targetGene))
vertices  <- data.frame(name = all_nodes, stringsAsFactors = FALSE)
vertices$node_class <- dplyr::case_when(
  vertices$name %in% edges$regulatoryGene ~ "TF",
  vertices$name %in% up_genes             ~ "Up",
  vertices$name %in% down_genes           ~ "Down",
  TRUE ~ NA_character_
)
vertices$node_class <- factor(vertices$node_class, levels = c("TF","Up","Down"))

g2 <- graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
V(g2)$degree <- degree(g2, mode = "all")

set.seed(123)
ggraph(g2, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = node_class), alpha = 0.85) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.8, fontface = "bold",
                 max.overlaps = 100, point.padding = unit(0.4, "lines")) +
  scale_color_manual(name = "Node Type",
                     values = c("TF"="#FFA07A","Up"="red","Down"="darkblue"),
                     breaks = c("TF","Up","Down")) +
  scale_size_continuous(range = c(3.5, 10)) +
  theme_void() +
  labs(title = "B cells TF Regulatory Network in Burns") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 12, 10, 12)) +
  guides(size = "none")

## =========================NK cells==========================================
NKcells <- log2(NKcells + 1)
fit  <- lmFit(NKcells, mm)
contr <- makeContrasts(typeburn - typecontrol, levels = colnames(coef(fit)))
tmp  <- contrasts.fit(fit, contr)
tmp  <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

DEG_NK <- subset(top.table, P.Value < 0.05 & abs(logFC) > 0)
DEG_NK$gene <- rownames(DEG_NK)
NK_genes    <- rownames(DEG_NK)

burn_id <- rownames(Phenotype[Phenotype$type == "burn", ])
expr_burn_NK <- NKcells[, burn_id, drop = FALSE]

available_tfs_NK <- intersect(tf_list, NK_genes)
expr_burn_matrix <- as.matrix(expr_burn_NK)

set.seed(123)
weightMatrix_burn <- GENIE3(expr_burn_matrix,
                            regulators   = available_tfs_NK,
                            nCores       = 8,
                            returnMatrix = TRUE)

threshold_burn <- quantile(weightMatrix_burn, probs = 0.95)
links_burn <- getLinkList(weightMatrix_burn, threshold = threshold_burn)

tf_activity <- links_burn %>%
  group_by(regulatoryGene) %>%
  summarise(n_targets = n(),
            max_weight = max(weight),
            mean_weight = mean(weight),
            .groups = "drop") %>%
  arrange(desc(n_targets))

top_tfs <- tf_activity$regulatoryGene
sub_links <- links_burn[links_burn$regulatoryGene %in% top_tfs, ]
g <- graph_from_data_frame(sub_links, directed = TRUE)
V(g)$degree <- degree(g, mode = "all")
V(g)$type   <- ifelse(V(g)$name %in% top_tfs, "Top TFs", "Target genes")

set.seed(123)
ggraph(g, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = type), alpha = 0.7) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
  scale_color_manual(values = c("Top TFs" = "#FFA07A", "Target genes" = "#87CEFA")) +
  scale_size_continuous(range = c(4, 12)) +
  theme_void() +
  labs(title = "NK cells TF Regulatory Network in Burns", color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10)) +
  guides(size = "none",
         color = guide_legend(title.theme = element_text(face = "bold", size = 14),
                              label.theme = element_text(face = "bold", size = 12)))

save(DEG_NK, weightMatrix_burn, file = "NK_TF.RData")

tf_target_genes_list <- sub_links %>%
  group_by(regulatoryGene) %>%
  summarise(targetGenes = list(unique(targetGene)), .groups = "drop")
tf_target_genes_export <- tf_target_genes_list %>%
  mutate(targetGenes = sapply(targetGenes, function(x) paste(x, collapse = ", ")))
write_xlsx(tf_target_genes_export, path = "NK_tf_target.xlsx")
write_xlsx(sub_links, path = "NK_weight.xlsx")

all_enrichment_results <- lapply(seq_len(nrow(tf_target_genes_list)), function(i) {
  tf <- tf_target_genes_list$regulatoryGene[i]
  tg <- tf_target_genes_list$targetGenes[[i]]
  run_enrichment_for_tf(tf, tg)
})
names(all_enrichment_results) <- tf_target_genes_list$regulatoryGene
go_combined <- do.call(rbind, all_enrichment_results)
write.csv(go_combined, "NK_TF_GO_enrichment_results.csv", row.names = FALSE)

# Highlight up/down targets
up_genes   <- DEG_NK$gene[DEG_NK$logFC > 0]
down_genes <- DEG_NK$gene[DEG_NK$logFC < 0]
edges <- links_burn[links_burn$targetGene %in% c(up_genes, down_genes), ]

all_nodes <- unique(c(edges$regulatoryGene, edges$targetGene))
vertices  <- data.frame(name = all_nodes, stringsAsFactors = FALSE)
vertices$node_class <- dplyr::case_when(
  vertices$name %in% edges$regulatoryGene ~ "TF",
  vertices$name %in% up_genes             ~ "Up",
  vertices$name %in% down_genes           ~ "Down",
  TRUE ~ NA_character_
)
vertices$node_class <- factor(vertices$node_class, levels = c("TF","Up","Down"))

g2 <- graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
V(g2)$degree <- degree(g2, mode = "all")

set.seed(123)
ggraph(g2, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = node_class), alpha = 0.85) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.8, fontface = "bold",
                 max.overlaps = 100, point.padding = unit(0.4, "lines")) +
  scale_color_manual(name = "Node Type",
                     values = c("TF"="#FFA07A","Up"="red","Down"="darkblue"),
                     breaks = c("TF","Up","Down")) +
  scale_size_continuous(range = c(3.5, 10)) +
  theme_void() +
  labs(title = "NK cells TF Regulatory Network in Burns") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 12, 10, 12)) +
  guides(size = "none")


## ========================Monocytes===========================================
Monocytes <- log2(Monocytes + 1)
fit  <- lmFit(Monocytes, mm)
contr <- makeContrasts(typeburn - typecontrol, levels = colnames(coef(fit)))
tmp  <- contrasts.fit(fit, contr)
tmp  <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

DEG_Mon <- subset(top.table, P.Value < 0.05 & abs(logFC) > 0)
DEG_Mon$gene <- rownames(DEG_Mon)
Mon_genes    <- rownames(DEG_Mon)

burn_id <- rownames(Phenotype[Phenotype$type == "burn", ])
expr_burn_Mon <- Monocytes[, burn_id, drop = FALSE]

available_tfs_Mon <- intersect(tf_list, Mon_genes)
expr_burn_matrix  <- as.matrix(expr_burn_Mon)

set.seed(123)
weightMatrix_burn <- GENIE3(expr_burn_matrix,
                            regulators   = available_tfs_Mon,
                            nCores       = 8,
                            returnMatrix = TRUE)

threshold_burn <- quantile(weightMatrix_burn, probs = 0.95)
links_burn <- getLinkList(weightMatrix_burn, threshold = threshold_burn)

tf_activity <- links_burn %>%
  group_by(regulatoryGene) %>%
  summarise(n_targets = n(),
            max_weight = max(weight),
            mean_weight = mean(weight),
            .groups = "drop") %>%
  arrange(desc(n_targets))

top_tfs <- tf_activity$regulatoryGene
sub_links <- links_burn[links_burn$regulatoryGene %in% top_tfs, ]
g <- graph_from_data_frame(sub_links, directed = TRUE)
V(g)$degree <- degree(g, mode = "all")
V(g)$type   <- ifelse(V(g)$name %in% top_tfs, "Top TFs", "Target genes")

set.seed(123)
ggraph(g, layout = "kk") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = type), alpha = 0.7) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
  scale_color_manual(values = c("Top TFs" = "#FFA07A", "Target genes" = "#87CEFA")) +
  scale_size_continuous(range = c(4, 12)) +
  theme_void() +
  labs(title = "Monocyte cells TF Regulatory Network in Burns", color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10)) +
  guides(size = "none",
         color = guide_legend(title.theme = element_text(face = "bold", size = 14),
                              label.theme = element_text(face = "bold", size = 12)))

save(DEG_Mon, weightMatrix_burn, file = "Mon_TF.RData")

tf_target_genes_list <- sub_links %>%
  group_by(regulatoryGene) %>%
  summarise(targetGenes = list(unique(targetGene)), .groups = "drop")
tf_target_genes_export <- tf_target_genes_list %>%
  mutate(targetGenes = sapply(targetGenes, function(x) paste(x, collapse = ", ")))
write_xlsx(tf_target_genes_export, path = "Mon_tf_target.xlsx")
write_xlsx(sub_links, path = "Mon_weight.xlsx")

all_enrichment_results <- lapply(seq_len(nrow(tf_target_genes_list)), function(i) {
  tf <- tf_target_genes_list$regulatoryGene[i]
  tg <- tf_target_genes_list$targetGenes[[i]]
  run_enrichment_for_tf(tf, tg)
})
names(all_enrichment_results) <- tf_target_genes_list$regulatoryGene
go_combined <- do.call(rbind, all_enrichment_results)
write.csv(go_combined, "Mon_TF_GO_enrichment_results.csv", row.names = FALSE)

# Highlight up/down targets
up_genes   <- DEG_Mon$gene[DEG_Mon$logFC > 0]
down_genes <- DEG_Mon$gene[DEG_Mon$logFC < 0]
edges <- links_burn[links_burn$targetGene %in% c(up_genes, down_genes), ]

all_nodes <- unique(c(edges$regulatoryGene, edges$targetGene))
vertices  <- data.frame(name = all_nodes, stringsAsFactors = FALSE)
vertices$node_class <- dplyr::case_when(
  vertices$name %in% edges$regulatoryGene ~ "TF",
  vertices$name %in% up_genes             ~ "Up",
  vertices$name %in% down_genes           ~ "Down",
  TRUE ~ NA_character_
)
vertices$node_class <- factor(vertices$node_class, levels = c("TF","Up","Down"))

g2 <- graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
V(g2)$degree <- degree(g2, mode = "all")

set.seed(123)
ggraph(g2, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = node_class), alpha = 0.85) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.8, fontface = "bold",
                 max.overlaps = 100, point.padding = unit(0.4, "lines")) +
  scale_color_manual(name = "Node Type",
                     values = c("TF"="#FFA07A","Up"="red","Down"="darkblue"),
                     breaks = c("TF","Up","Down")) +
  scale_size_continuous(range = c(3.5, 10)) +
  theme_void() +
  labs(title = "Monocytes TF Regulatory Network in Burns") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 12, 10, 12)) +
  guides(size = "none")


## ==========================Eos=========================================
Eosinophils <- log2(Eosinophils + 1)
fit  <- lmFit(Eosinophils, mm)
contr <- makeContrasts(typeburn - typecontrol, levels = colnames(coef(fit)))
tmp  <- contrasts.fit(fit, contr)
tmp  <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)

DEG_Eos <- subset(top.table, P.Value < 0.05 & abs(logFC) > 0)
DEG_Eos$gene <- rownames(DEG_Eos)
Eos_genes    <- rownames(DEG_Eos)

Eos_samp <- Eosinophils[Eos_genes, , drop = FALSE]
burn_id  <- rownames(Phenotype[Phenotype$type == "burn", ])
expr_burn_Eos <- Eos_samp[, burn_id, drop = FALSE]

# You originally ran GENIE3 without restricting regulators for Eos (kept as-is)
expr_burn_matrix <- as.matrix(expr_burn_Eos)

set.seed(123)
weightMatrix_burn <- GENIE3(expr_burn_matrix,
                            # regulators = intersect(tf_list, rownames(expr_burn_Eos)), # optional
                            nCores       = 8,
                            returnMatrix = TRUE)

threshold_burn <- quantile(weightMatrix_burn, probs = 0.99)
links_burn <- getLinkList(weightMatrix_burn, threshold = threshold_burn)

tf_activity <- links_burn %>%
  group_by(regulatoryGene) %>%
  summarise(n_targets = n(),
            max_weight = max(weight),
            mean_weight = mean(weight),
            .groups = "drop") %>%
  arrange(desc(n_targets))

top_tfs <- head(tf_activity$regulatoryGene, 6)
sub_links <- links_burn[links_burn$regulatoryGene %in% top_tfs, ]
g <- graph_from_data_frame(sub_links, directed = TRUE)
V(g)$degree <- degree(g, mode = "all")
V(g)$type   <- ifelse(V(g)$name %in% top_tfs, "Top TFs", "Target genes")

set.seed(123)
ggraph(g, layout = "fr") +
  geom_edge_link(color = "#D3D3D3",
                 arrow = arrow(length = unit(2, "mm")),
                 end_cap = circle(3, 'mm'),
                 start_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = degree, color = type), alpha = 0.7) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, fontface = "bold") +
  scale_color_manual(values = c("Top TFs" = "#FFA07A", "Target genes" = "#87CEFA")) +
  scale_size_continuous(range = c(4, 12)) +
  theme_void() +
  labs(title = "Eosinophils TF Regulatory Network in Burns", color = "Node Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10)) +
  guides(size = "none",
         color = guide_legend(title.theme = element_text(face = "bold", size = 14),
                              label.theme = element_text(face = "bold", size = 12)))

save(DEG_Eos, weightMatrix_burn, file = "Eos_TF.RData")



