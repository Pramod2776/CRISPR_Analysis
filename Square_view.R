library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(ggpubr)
library(poolr)
library(reshape2)
library(rstatix)
library(gginnards)
library(ggthemes)
library(ggeasy)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(gridExtra)
library(ggplotify)
library(cowplot)

## Data Preprocessing
depmap <- readRDS("./datasets/depmap/depmap_screen_data.rds")
metadata = depmap$meta_data
aml_cell_lines <-
  rownames(depmap$meta_data[grep("Melanoma", depmap$meta_data$subtype), ])
non_aml_cl <-
  rownames(depmap$meta_data[-grep("Melanoma", depmap$meta_data$subtype), ])

Genes = readxl::read_excel("./src/aml-2019/resources/enriched.xlsx", sheet =3) %>%
  as.data.frame()
Genes = Genes$`Gene/Compound
dependency_score = depmap$gene_fitness_crispr_depmap_q3 %>%
  as.data.frame()

ds_df = dependency_score %>%
  dplyr::filter(rownames(depmap$gene_fitness_crispr_depmap_q3) %in% Genes)

ds_df_aml = ds_df %>%
  as.data.frame() %>%
  dplyr::select(all_of(aml_cell_lines))

col_fun = colorRamp2(c(-2, 0, 2), rev(c("red", "white", "#4682B4")))

ds_df_aml[is.na(ds_df_aml)] <- 0   
ds_df_aml = ds_df_aml %>% 
  dplyr::filter(!if_all(everything(), ~ . == 0))

ds_df_aml = ds_df_aml %>%
  select(where(~ any(. != 0)))


## Generate heatmap
ht2 = Heatmap(as.matrix(ds_df_aml),  cluster_rows = TRUE,
              cluster_columns = TRUE,
              row_labels = rownames(as.matrix(ds_df_aml)),
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              column_title = "Melanoma_cell_lines",
              row_names_gp = gpar(fontsize = 6),
              col = col_fun)

ds_df_non_aml = ds_df %>%
  as.data.frame() %>%
  dplyr::select(all_of(non_aml_cl))

ds_df_non_aml[is.na(ds_df_non_aml)] <- 0   

ds_df_non_aml = ds_df_non_aml %>%
  select(where(~ any(. != 0)))

ds_df_non_aml_t =t(apply(ds_df_non_aml,1,sort))

clusters= draw(ht2)
gene_cluster = row_order(clusters)                       
rownames(ds_df_non_aml_t)[gene_cluster[1]]     

clu_df <- lapply(1:length(gene_cluster), function(i){
  print(i)
  out <- data.frame(coordinates = rownames(ds_df_non_aml_t)[gene_cluster[[i]]],
                    Cluster = paste0("cluster", i), stringsAsFactors = FALSE)
  return(out)
}) %>%
  do.call(rbind, .)

ds_df_non_aml_t_df = ds_df_non_aml[match(clu_df$coordinates, rownames(ds_df_non_aml)),]

ds_df_non_aml_t_S =t(apply(ds_df_non_aml_t_df,1,sort))

dat1 = ds_df_non_aml_t_S
dat1 = t(dat1)

dat1m <- melt(cbind(dat1, ind = rownames(dat1)), id.vars = c('ind'))
levels(dat1m$Var2)

dd =dat1m %>% group_by(Var2) %>% 
  arrange(desc(Var1))%>% 
  mutate(pct = value- mean(value)/sd(value), right = cumsum(pct), left=lag(right, default=-1471), max = 10704)
g = ggplot(dd) + 
  geom_rect(aes(xmin=max, xmax=left, ymin=as.numeric(Var2)-.4, ymax=as.numeric(Var2)+.4, fill= ((value)))) + 
  ##scale_y_continuous(labels=levels(dd$Var2), breaks=1:nlevels(dd$Var2))+
  # theme_classic2()+
  scale_fill_distiller(palette = "RdBu")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ##theme_minimal()+
  ##theme_classic()+
  theme_void()

grob = grid.grabExpr(draw(Heatmap(as.matrix(ds_df_aml),  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
                                  show_row_names=F ,
                                  show_column_names = F,
                                  show_column_dend = FALSE,
                                  show_row_dend = FALSE,
                                  row_names_gp = gpar(fontsize = 6),
                                  col = col_fun))) 

p1 = plot_grid(g, grob)
file <- paste0("m_crispr", fileext = ".pdf")
save_plot(file, p1, ncol = 2, base_asp = 0.5)

#### DE Analyis
library(edgeR)
depmap_mrna <- depmap[["mrna"]]

# Annotate which cell lines have expression and crispr data
dp_cls <- colnames(depmap$mrna)
expr_cls <- dp_cls[colSums(depmap$mrna, na.rm=TRUE) != 0]
crispr_cls <- dp_cls[colSums(depmap$gene_fitness_crispr_depmap_21q3, na.rm=TRUE) != 0]

expr_cls_crispr_cls = intersect(expr_cls, crispr_cls)

depmap_mrna = depmap_mrna %>%
  as.data.frame()

depmap_mrna[is.na(depmap_mrna)] <- 0   
dim(depmap_mrna)

depmap_expr_cls_crispr_clsl = depmap_mrna %>%
  as.data.frame() %>%
  dplyr::select(all_of(expr_cls_crispr_cls))

dim(depmap_expr_cls_crispr_clsl)
aml_interesect = intersect(as.character(colnames(depmap_expr_cls_crispr_clsl)), as.character(aml_cell_lines))
print(length(aml_interesect))

depmap_mrna_aml = depmap_expr_cls_crispr_clsl %>%
  as.data.frame() %>%
  dplyr::select(aml_interesect)
    
dim(depmap_mrna_aml)

non_aml_interesect = intersect(as.character(colnames(depmap_expr_cls_crispr_clsl)), as.character(non_aml_cl))

depmap_mrna_non_aml = depmap_expr_cls_crispr_clsl %>%
  as.data.frame() %>%
  dplyr::select(all_of(non_aml_interesect))

dim(depmap_mrna_non_aml)

depmap_mrna_non_aml = depmap_mrna_non_aml %>%
  select(where(~ any(. != 0)))

depmap_mrna_aml_non_aml = cbind(depmap_mrna_aml, depmap_mrna_non_aml)
dim(depmap_mrna_aml_non_aml)

Groups <- rep(c("Melanoma", "non_Melanoma"), c("65", "924"))
#########
d <- DGEList(counts=depmap_mrna_aml_non_aml,group=factor(Groups))
dgeFull = d
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
head(dgeFull$counts)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples
y =dgeFull
design <- model.matrix(~0+Groups)
colnames(design) = c("Melanoma", "non_Melanoma")
BvsA = makeContrasts(Melanoma-non_Melanoma, levels=colnames(design))

y <- estimateDisp(y, design)

fit <- glmFit(y, contrast=BvsA)
lrt <- glmLRT(fit, contrast=BvsA)
TT = topTags(lrt, n=nrow(dgeFull$counts))

table_aml_vs_nonaml= topTags(lrt, n=nrow(dgeFull$counts))
write.csv(table_aml_vs_nonaml, "table_melanoma_vs_nonmelanoma_redo.csv")

##enhanced valcono plot

head(TT)
TT = TT$table
de = TT
dim(de)

# Remove genes with fold-change = NA
de <- de[is.na(de$logFC) == F, ]

## --------------------------
## Fix inputs for volcano plot
## --------------------------
## fold-change
fc <- de$logFC
head(fc)
max(fc)
min(fc)
fc[fc < -5] <- -5
fc[fc > 5] <- 5
table(is.na(fc))

# adjusted p-value
qval <- de$FDR
table(is.na(qval))
qval[is.na(qval)] <- 1
table(is.na(qval))

# gene names
de$gene_name = rownames(de)
genes <- make.unique(de$gene_name)

de$Sigin = NULL

de_results = de %>%
  mutate(Sigin = ifelse((logFC >= 2 & FDR <= 0.05), 1,
                    ifelse((logFC <= -2 & FDR <= 0.05), -1, NA)))
de = de_results               


# significance
sig <- de_results$Sigin
table(is.na(sig) == T)
sig[is.na(sig)] <- 0

labels = Genes[1:20]

# Create a data frame using gene names, fold changes, and q-values
df <- data.frame(gene_name = as.character(genes), lfc = fc, q = qval, stringsAsFactors = F)
rownames(df) <- genes
head(df)   

## Plot volcano
# y-axis upper limit
max_pval <- -log10(min(qval)) + 0.1


# Create a named vector of custom colors
col_scheme <- c("Up"="red2", "Down"="blue", "N.S."="grey")
cols <- rep(col_scheme['N.S.'], length(genes))
cols[sig == 1] <- col_scheme['Up']
cols[sig == -1] <- col_scheme['Down']
names(cols)[cols == col_scheme['Up']] <- "Up"
names(cols)[cols == col_scheme['Down']] <- "Down"
names(cols)[cols == col_scheme['N.S.']] <- "N.S."


df_intersect = df %>%
  dplyr::filter(q <= 0.05)


labels = intersect (df_intersect$gene_name, Genes)

df_Genes = de %>%
  dplyr::filter(gene_name %in%  Genes)

# Set it globally:
options(ggrepel.max.overlaps = Inf)
plt <- EnhancedVolcano(df, x = 'lfc', y = 'q', lab = df$gene_name,
                       pCutoff = 0.05, FCcutoff = 2.0, 
                       gridlines.major = FALSE, gridlines.minor = FALSE, 
                       drawConnectors = T, legendLabSize = 12,
                       cutoffLineCol = "red", colAlpha = 0.75, 
                       cutoffLineType = "dashed", border = "full",
                       colCustom = cols, legendPosition = "right",
                       pointSize = 2, cutoffLineWidth = 0.4,
                       labFace = "plain", subtitle = "",
                       ylim = c(0, max_pval), xlim = c(-5, 5),
                       axisLabSize = 12, captionLabSize = 12, 
                       xlab = "Log2 Fold-change", ylab = "-Log10 adjusted p-value", title = "", 
                       caption = paste0('Total = ', nrow(df), ' genes'),
                       typeConnectors = "closed", legendIconSize = 2,
                       selectLab = labels, borderWidth = 1.5, boxedLabels = TRUE,
                       maxoverlapsConnectors = Inf
)

print(plt)
ggsave(plt, filename = file.path(".","mRNA_Volcano_select.pdf"), height = 12, width = 12)


### CERES scores lmfit

expr_cls_crispr_cls = intersect(expr_cls, crispr_cls)
depmap_crispr = depmap$gene_fitness_crispr_depmap_21q3
crispr_genes = rownames(table_aml_vs_nonaml)
dim(depmap_crispr)

depmap_crispr_extract = depmap_crispr %>%
  as.data.frame() %>%
  dplyr::select(expr_cls_crispr_cls) %>%
  dplyr::filter(rownames(.) %in% crispr_genes)




aml_interesect_crispr = intersect(as.character(colnames(depmap_crispr_extract)), as.character(aml_cell_lines))
depmap_crispr_extract_aml = depmap_crispr_extract %>%
  dplyr::select(all_of(aml_interesect_crispr))

non_aml_interesect_crispr = intersect(as.character(colnames(depmap_crispr_extract)), as.character(non_aml_cl))
depmap_crispr_extract_non_aml = depmap_crispr_extract %>%
  dplyr::select(all_of(non_aml_interesect_crispr))
dim(depmap_crispr_extract_aml)
depmap_crispr_extract_aml_row =rowMeans(depmap_crispr_extract_aml, na.rm = F)
dim(depmap_crispr_extract_non_aml)
depmap_crispr_extract_non_aml_row =rowMeans(depmap_crispr_extract_non_aml, na.rm = F)
crispr_row_means= cbind(depmap_crispr_extract_aml_row, depmap_crispr_extract_non_aml_row)
colnames(crispr_row_means)

crispr_row_means[is.na(crispr_row_means)] = 0

crispr_row_means = crispr_row_means %>%
  as.data.frame() %>%
  dplyr::mutate(average_CERES = depmap_crispr_extract_aml_row - depmap_crispr_extract_non_aml_row)

depmap_crispr_aml_non_aml = cbind(depmap_crispr_extract_aml, depmap_crispr_extract_non_aml)
dim(depmap_crispr_aml_non_aml)
Groups
mod1 <- model.matrix(~ factor(Groups) -1)
colnames(mod1) <- c("Melanoma", "non_Melanoma")
fit2 <- lmFit(depmap_crispr_aml_non_aml, mod1)
contrast.matrix <- makeContrasts("Melanoma-non_Melanoma", levels = mod1)
contrast.matrix
fit2C <- contrasts.fit(fit2, contrast.matrix)
fit2C <- eBayes(fit2C)
topTable(fit2C)

res <- decideTests(fit2C, p.value=0.10)
summary(res)


tt <- topTable(fit2C, n=Inf)


##enhanced valcono plot

TT1 = tt
head(TT1)

de1 = TT1
dim(de1)

# Remove genes with fold-change = NA
de1 <- de1[is.na(de1$logFC) == F, ]

## --------------------------
## Fix inputs for volcano plot
## --------------------------
## fold-change
fc <- de1$logFC
head(fc)
max(fc)
min(fc)
fc[fc < -0.9] <- -0.9
fc[fc > 0.5] <- 0.5
table(is.na(fc))

# adjusted p-value
qval <- de1$adj.P.Val
table(is.na(qval))
qval[is.na(qval)] <- 1
table(is.na(qval))

# gene names
de1$gene_name = rownames(de1)
genes <- make.unique(de1$gene_name)

de1_results = de1 %>%
  mutate(Sigin = ifelse((logFC >= 0.2 & adj.P.Val <= 0.05), 1,
                        ifelse((logFC <= -0.2 & adj.P.Val <= 0.05), -1, NA)))

de1 = de1_results               


# significance
sig <- de1_results$Sigin
table(is.na(sig) == T)
sig[is.na(sig)] <- 0

labels = Genes[1:20]

# Create a data frame using gene names, fold changes, and q-values
df <- data.frame(gene_name = as.character(genes), lfc = fc, q = qval, stringsAsFactors = F)
rownames(df) <- genes
head(df)   

## Plot volcano
# y-axis upper limit
max_pval <- -log10(min(qval)) + 0.1


# Create a named vector of custom colors
col_scheme <- c("Up"="red2", "Down"="blue", "N.S."="grey")
cols <- rep(col_scheme['N.S.'], length(genes))
cols[sig == 1] <- col_scheme['Up']
cols[sig == -1] <- col_scheme['Down']
names(cols)[cols == col_scheme['Up']] <- "Up"
names(cols)[cols == col_scheme['Down']] <- "Down"
names(cols)[cols == col_scheme['N.S.']] <- "N.S."

# Render the volcano plot

##labels = Genes[1:10]

df_intersect = df %>%
  dplyr::filter(q <= 0.05)
head(df_intersect)

dim(de1_results)

de1_results_pos = de1_results %>%
  dplyr::filter(logFC >= 0)

labels_pos = intersect (de1_results_pos$gene_name, Genes)


labels = intersect (df_intersect$gene_name, Genes)

df_Genes = de1 %>%
  dplyr::filter(gene_name %in%  Genes)

write.csv(df_Genes, "df_Genes_vol_crispr.csv")

##labels = c("SOX10","PAX3", "IGSF11", "SOX5", "BCL2A1", "PRSS1")


labels = labels[1:30]


# Set it globally:
options(ggrepel.max.overlaps = Inf)
plt <- EnhancedVolcano(df, x = 'lfc', y = 'q', lab = df$gene_name,
                       pCutoff = 0.05, FCcutoff = 0.2, 
                       gridlines.major = FALSE, gridlines.minor = FALSE, 
                       drawConnectors = T, legendLabSize = 12,
                       cutoffLineCol = "red", colAlpha = 0.75, 
                       cutoffLineType = "dashed", border = "full",
                       colCustom = cols, legendPosition = "right",
                       pointSize = 2, cutoffLineWidth = 0.4,
                       labFace = "plain", subtitle = "",
                       ylim = c(0, max_pval), xlim = c(-1, 1),
                       axisLabSize = 12, captionLabSize = 12, 
                       xlab = "Log2 Fold-change", ylab = "-Log10 adjusted p-value", title = "", 
                       caption = paste0('Total = ', nrow(df), ' genes'),
                       typeConnectors = "closed", legendIconSize = 2,
                       selectLab = labels, borderWidth = 1.5, boxedLabels = TRUE,
                       maxoverlapsConnectors = Inf
)

print(plt)

ggsave(plt, filename = file.path(".","Ceres_volcano_select.pdf"), height = 12, width = 12)





######### Square plot

dim(table_aml_vs_nonaml)
dim(crispr_row_means)

df_lfc = data.frame(Gene = rownames(table_aml_vs_nonaml$table),
           logfc = table_aml_vs_nonaml$table$logFC)

df_avg = data.frame(Gene = rownames(crispr_row_means),
           CERRS = crispr_row_means$average_CERES)

df_squareview = df_lfc %>%
  dplyr::left_join(df_avg, by = "Gene")
library("MAGeCKFlute")
require(ggExtra)
suppressPackageStartupMessages(library("ggfortify"))
suppressPackageStartupMessages(library("ggforce"))

p = SquareView(df_squareview, ctrlname = "logfc", treatname = "CERRS", label = "Gene", genelist = Genes1)

##max(df_squareview$logfc)
p =p + labs(y = "Difference in Average CERES score cell lines", x = "log2FC diff expression in cell lines") 
p3 = ggMarginal(p, type="histogram", binwidth = 0.2, size=10, fill = "#bdbdbd")
ggsave('crispr_square.png', p3, height = 8, width = 8)



