#!/usr/bin/env Rscript

library(Biostrings)
library(optparse)
library(tidyverse)
library(ggpubr)

# Función para leer un FASTA y devolver número de secuencias o NULL si no existe o está vacío
safe_read_fasta <- function(file) {
  if (!file.exists(file)) {
    message(paste("Warning:", file, "no encontrado, contando como 0 secuencias."))
    return(NULL)
  }
  fasta <- tryCatch({
    readAAStringSet(file)
  }, error = function(e) {
    message(paste("Warning: Error al leer", file, "- contando como 0 secuencias."))
    return(NULL)
  })
  if (is.null(fasta) || length(fasta) == 0) {
    return(NULL)
  }
  return(fasta)
}

# Definir opciones para línea de comandos
option_list <- list(
  make_option(c("-a", "--annotated"), type="character", help="Archivo annotated.fasta"),
  make_option(c("-u", "--unannotated"), type="character", help="Archivo unannotated.fasta"),
  make_option(c("-p", "--annotated_prot"), type="character", help="Archivo annotated_proteins.fasta"),
  make_option(c("-q", "--unannotated_prot"), type="character", help="Archivo unannotated_proteins.fasta"),
  make_option(c("-o", "--output"), type="character", default="results/analysis", help="Directorio de salida")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Crear directorio de salida si no existe
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Leer cada archivo de forma segura
annotated <- safe_read_fasta(opt$annotated)
unannotated <- safe_read_fasta(opt$unannotated)
annotated_prot <- safe_read_fasta(opt$annotated_prot)
unannotated_prot <- safe_read_fasta(opt$unannotated_prot)

# Contar secuencias por grupo
summary_fasta <- function(fasta, name) {
  n <- ifelse(is.null(fasta), 0, length(fasta))
  message(paste("Número de secuencias en", name, ":", n))
  return(n)
}

n_annotated <- summary_fasta(annotated, "annotated.fasta")
n_unannotated <- summary_fasta(unannotated, "unannotated.fasta")
n_annotated_prot <- summary_fasta(annotated_prot, "annotated_proteins.fasta")
n_unannotated_prot <- summary_fasta(unannotated_prot, "unannotated_proteins.fasta")

# Crear dataframe para gráficos de conteo
counts <- tibble(
  Category = factor(c("Annotated (nuc)", "Unannotated (nuc)", "Annotated (prot)", "Unannotated (prot)"),
                    levels = c("Annotated (nuc)", "Unannotated (nuc)", "Annotated (prot)", "Unannotated (prot)")),
  Count = c(n_annotated, n_unannotated, n_annotated_prot, n_unannotated_prot)
)

# Función para tema con fondo blanco
theme_white <- function() {
  theme_minimal(base_family = "") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA)
    )
}

# Gráfico de barras: Conteo de secuencias
p_counts <- ggplot(counts, aes(x=Category, y=Count, fill=Category)) +
  geom_bar(stat="identity", color="black") +
  labs(title="Conteo de secuencias por categoría", y="Número de secuencias", x=NULL) +
  theme_white() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")

# Función para obtener longitud de secuencias y generar dataframe
get_lengths <- function(fasta, group_name) {
  if (is.null(fasta)) return(tibble(Length=integer(), Group=character()))
  lengths <- width(fasta)
  tibble(Length=lengths, Group=group_name)
}

lengths_df <- bind_rows(
  get_lengths(annotated, "Annotated (nuc)"),
  get_lengths(unannotated, "Unannotated (nuc)"),
  get_lengths(annotated_prot, "Annotated (prot)"),
  get_lengths(unannotated_prot, "Unannotated (prot)")
)

# Boxplot de longitud de secuencias
p_lengths <- ggplot(lengths_df, aes(x=Group, y=Length, fill=Group)) +
  geom_boxplot() +
  labs(title="Distribución de longitud de secuencias", y="Longitud (nt o aa)", x=NULL) +
  theme_white() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")

# Histograma de longitudes
p_hist <- ggplot(lengths_df, aes(x=Length, fill=Group)) +
  geom_histogram(bins=30, alpha=0.7, position="identity") +
  labs(title="Histograma de longitudes de secuencias", x="Longitud", y="Frecuencia") +
  theme_white() +
  theme(legend.position = "bottom")

# Densidad de longitudes
p_density <- ggplot(lengths_df, aes(x=Length, color=Group, fill=Group)) +
  geom_density(alpha=0.3) +
  labs(title="Densidad de longitudes de secuencias", x="Longitud", y="Densidad") +
  theme_white() +
  theme(legend.position = "bottom")

# Boxplot logarítmico para mejor visualización (longitudes grandes)
p_lengths_log <- ggplot(lengths_df, aes(x=Group, y=Length, fill=Group)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title="Distribución logarítmica de longitud de secuencias", y="Longitud (log10)", x=NULL) +
  theme_white() +
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")

# Guardar resultados y gráficos
write_csv(counts, file.path(opt$output, "sequence_counts.csv"))
ggsave(filename = file.path(opt$output, "sequence_counts.png"), plot = p_counts, width = 7, height = 5)
ggsave(filename = file.path(opt$output, "sequence_length_distribution.png"), plot = p_lengths, width = 8, height = 6)
ggsave(filename = file.path(opt$output, "sequence_length_histogram.png"), plot = p_hist, width = 8, height = 6)
ggsave(filename = file.path(opt$output, "sequence_length_density.png"), plot = p_density, width = 8, height = 6)
ggsave(filename = file.path(opt$output, "sequence_length_log_boxplot.png"), plot = p_lengths_log, width = 8, height = 6)

message(paste("Gráficos guardados en:", opt$output))