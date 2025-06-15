#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Faltan argumentos: --features y --output")
}

features_file <- args[which(args == "--features") + 1]
output_dir <- args[which(args == "--output") + 1]

# Comprobar si el archivo existe y no está vacío
if (!file.exists(features_file) || file.info(features_file)$size == 0) {
  cat("El archivo features.txt está vacío. Finalizando sin análisis.\n")
  quit(save = "no", status = 0)
}

library(tidyverse)
library(optparse)
library(factoextra) # para clustering y PCA
library(ggpubr)

option_list <- list(
  make_option(c("-f", "--features"), type="character", help="Archivo de características (txt o csv)"),
  make_option(c("-o", "--output"), type="character", default="results/analysis", help="Directorio de salida")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Leer características
features <- read.delim(opt$features, stringsAsFactors = FALSE)

# Eliminar columna id para análisis
feat_data <- features %>% select(-id)

# Escalar datos
feat_scaled <- scale(feat_data)

# PCA
pca_res <- prcomp(feat_scaled, center = TRUE, scale. = TRUE)

# Graficar PCA
pca_plot <- fviz_pca_ind(pca_res, geom.ind = "point", col.ind = "cos2",
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE) +
  ggtitle("PCA de características proteicas")

ggsave(filename=file.path(opt$output, "pca_features.png"), plot=pca_plot, width=8, height=6)

# Clustering jerárquico
dist_mat <- dist(feat_scaled)
hc <- hclust(dist_mat, method = "ward.D2")

# Dendrograma
dend_plot <- fviz_dend(hc, k = 3, # número de clusters arbitrario, puedes cambiar
                      rect = TRUE, rect_fill = TRUE,
                      main = "Dendrograma Clustering Jerárquico")

ggsave(filename=file.path(opt$output, "dendrogram.png"), plot=dend_plot, width=8, height=6)

# K-means clustering con k=3 (puedes ajustar)
set.seed(123)
km_res <- kmeans(feat_scaled, centers = 3, nstart = 25)

features$cluster <- as.factor(km_res$cluster)

# Gráfico de distribución (violin) de una característica relevante por cluster
violin_plot <- ggplot(features, aes(x=cluster, y=mw, fill=cluster)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  theme_minimal() +
  labs(title="Distribución del peso molecular por cluster", y="Peso molecular")

ggsave(filename=file.path(opt$output, "violin_mw_by_cluster.png"), plot=violin_plot, width=7, height=5)

# Guardar clusters con IDs para referencia
write.csv(features[, c("id", "cluster")], file=file.path(opt$output, "protein_clusters.csv"), row.names=FALSE)

cat("Análisis de características y clustering completado. Resultados en:", opt$output, "\n")
