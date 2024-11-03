library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(ggpubr)
library(poolr)
library(reshape2)
library(rstatix)
library(gginnards)

## Create directories
cidr <- getwd()
dir_depmap <- "/datasets/depmap/"
dir.create(file.path(cidr,dir_depmap), recursive = TRUE)

dir_out <- "/src/results/"
dir_resouces <-"/src/resources/"

dir.create(file.path(cidr,dir_out), recursive = TRUE)
dir.create(file.path(cidr,dir_resouces), recursive = TRUE)
depmap <- readRDS("./datasets/depmap/depmap_data.rds")

## Function to subset data 
subset_object <- function(
    data_object, 
    meta_selection, 
    meta_id="DepMap_ID"){
  
  new_data_object <- list()
  new_data_object$genes <- data_object$genes
  new_data_object$meta_data <- data_object$meta_data[data_object$meta_data[[meta_id]] %in% meta_selection, ]
  new_data_object$mrna <- data_object$mrna[, meta_selection]
  new_data_object$scna <- data_object$snca[, meta_selection]
  new_data_object$gene_fitness_crispr_depmap_q3 <- data_object$gene_fitness_crispr_depmap_q3[, meta_selection]
  new_data_object$gene_fitness_rnai_depmap <- data_object$gene_fitness_rnai_depmap[, meta_selection]
  
  return(new_data_object)
  
}

## Genes of Interest

Genes = c("A", "B")

## Process AML data
l_cell_lines <- rownames(depmap$meta_data[grep("Leukemia", depmap$meta_data$cancer),])
l_depmap <- subset_object(depmap, meta_selection = l_cell_lines)
l_crispr_fitness = l_depmap$gene_fitness_crispr_depmap_q3 %>% as.data.frame()


l_crispr_fitness = l_crispr_fitness %>%
  dplyr::mutate(Gene = row.names(l_crispr_fitness))%>%
  dplyr::filter(Gene %in% Genes)

l_crispr_fitness$Gene = NULL

l_reshape = reshape2::melt(t(l_crispr_fitness))
colnames(l_reshape) <- c("cell_line", "gene_name", "gene_fitness_score")
l_reshape$L_type <- rep("Leukemia", length(unique(l_reshape$gene_name)))

## Process Non AML data
non_l_cell_lines <- rownames(depmap$meta_data[-grep("Leukemia", depmap$meta_data$cancer),])
non_l_depmap <- subset_object(depmap, meta_selection = non_l_cell_lines)

non_l_crispr_fitness = non_l_depmap$gene_fitness_crispr_depmap_q3 %>% as.data.frame()
head(non_l_crispr_fitness)

non_l_crispr_fitness = non_l_crispr_fitness %>%
  dplyr::mutate(Gene = row.names(non_l_crispr_fitness))%>%
  dplyr::filter(Gene %in% Genes)

non_l_crispr_fitness$Gene = NULL
non_l_reshape = reshape2::melt(t(non_l_crispr_fitness))
colnames(non_l_reshape) <- c("cell_line", "gene_name", "gene_fitness_score")
non_l_reshape$L_type <- rep("non_Leukemia", length(unique(non_l_reshape$gene_name)))
l_gene_format = rbind(l_reshape, non_l_reshape)

## Generate boxplots
list = lapply(Genes, function(gene){

  l_gene_format_M = l_gene_format %>%
    dplyr::filter(gene_name == gene) 
  print(l_gene_format_M)
  
  my_comparisons <- list(c("Leukemia", "non_Leukemia"))
  
  p <- ggboxplot(l_gene_format_M, x = "L_type", y = "gene_fitness_score",
                 color = "L_type", palette =c("#1034A6", "gray"),
                 add = "jitter")
  p = p + stat_compare_means(label.y = 1.3, size = 5)
  p$layers[[2]]$aes_params$textsize <- 5
  p <- p + ylim(c(-1.5, 1.5))
  g <- p + stat_compare_means(
    comparisons = my_comparisons, 
    method = "wilcox.test",
    paired = FALSE,
    aes(group = AML_type))

  pg <- ggplot_build(g)
  
  data.frame(gene = gene,
             test = pg$data[[3]]$method,
             Pvalue = pg$data[[3]]$p.adj,
             log10_Pvalue = -log10(pg$data[[3]]$p.adj))
  
}
) 
Pvalues_wilcoxon = as.data.frame(do.call(rbind, list))
Pvalues_wilcoxon = Pvalues_wilcoxon %>%
  dplyr::mutate(FDR = p.adjust(Pvalues_wilcoxon$Pvalue, method = "BH")) %>%
  dplyr::mutate(log10_FDR = -log10(FDR))

write.csv(Pvalues_wilcoxon, file.path(".", dir_out, "Pvalues_wilcoxon_A_nonA.csv"))  
