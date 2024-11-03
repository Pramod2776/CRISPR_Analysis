library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(ggpubr)
library(poolr)
library(reshape2)
library(rstatix)
library(gginnards)
library(GSVA)
library(edgeR)
library(limma)
library(Biobase)
library(GSEABase)
library(GSVAdata)
library(EnhancedVolcano)
library(limma)

## Create directories
cidr <- getwd()
dir_depmap <- "/datasets/depmap/"
dir.create(file.path(cidr,dir_depmap), recursive = TRUE)

dir_out <- "/src/SSGSEA/results/"
dir_resouces <-"/src/SSGSEA/resources/"

dir.create(file.path(cidr,dir_out), recursive = TRUE)
dir.create(file.path(cidr,dir_resouces), recursive = TRUE)
depmap <- readRDS("./datasets/depmap/depmap_Q3_screen_data.rds")
depmap_crispr_aml_nonaml = depmap$gene_fitness_crispr_depmap_q3

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
  new_data_object$gene_fitness_crispr_depmap_21q3 <- data_object$gene_fitness_crispr_depmap_q3[, meta_selection]
  new_data_object$gene_fitness_rnai_depmap <- data_object$gene_fitness_rnai_depmap[, meta_selection]
  
  return(new_data_object)
  
}


## Genes to be filtered
common_ess = read.csv(file.path(".",dir_resouces,"CommonEssentials.csv"))

Genes = common_ess$gene

## Corum complexes
corum_complexes = readxl::read_excel(file.path(".",dir_resouces,"corum.xlsx"), sheet =1) %>% as.data.frame()

## Process AML data
aml_cell_lines <- rownames(depmap$meta_data[grep("AML", depmap$meta_data$subtype),])
aml_depmap <- subset_object(depmap, meta_selection = aml_cell_lines)
aml_crispr_fitness = aml_depmap$gene_fitness_crispr_depmap_q3 %>% as.data.frame()
aml_crispr_fitness = aml_crispr_fitness %>%
  dplyr::mutate(Gene = row.names(aml_crispr_fitness))%>%
  dplyr::filter(!(Gene %in% Genes))

aml_crispr_fitness$Gene = NULL

## Process Non AML data
non_aml_cell_lines <- rownames(depmap$meta_data[-grep("AML", depmap$meta_data$subtype),])
non_aml_depmap <- subset_object(depmap, meta_selection = non_aml_cell_lines)
non_aml_crispr_fitness = non_aml_depmap$gene_fitness_crispr_depmap_21q3 %>% as.data.frame()
non_aml_crispr_fitness = non_aml_crispr_fitness %>%
  dplyr::mutate(Gene = row.names(non_aml_crispr_fitness))%>%
  dplyr::filter(!(Gene %in% Genes))

non_aml_crispr_fitness$Gene = NULL

aml_gene_format = cbind(aml_crispr_fitness, non_aml_crispr_fitness)
rownames(aml_gene_format) = rownames(aml_crispr_fitness)

##GSVA create expression set
aml_gene_format_exprs = as.matrix(aml_gene_format)
aml_gene_format_exprs[is.na(aml_gene_format_exprs)] <- 0
pData = data.frame("cancer_type" = c(rep("AML",53), rep("nonAML", 1743)))
rownames(pData) = colnames(aml_gene_format_exprs)

phenoData <- new("AnnotatedDataFrame",
                 data=pData)

aml_noaml_exprs_eset = ExpressionSet(assayData=aml_gene_format_exprs,
               phenoData=phenoData)

##GSVA create genesets

corum_complexes_Subset = corum_complexes %>%
  dplyr::select("ComplexID", "ComplexName", "Gene_name") %>%
  as.data.frame()

list <- strsplit(as.character(corum_complexes_Subset$Gene_name), ";")
names(list) = corum_complexes_Subset$ComplexID

l2 = lapply(names(list), function(name){
  l2 = list[[name]] %>%
    as.vector()
  l2[!(l2 %in% Genes)]
  
})


l1 <- list[sapply(l2, function(x) length(x) >= 5)]


aml_noaml_exprs_es <- gsva(aml_noaml_exprs_eset, l1, method = "ssgsea")
corum_complexes_names_matched = (corum_complexes_Subset[match(rownames(aml_noaml_exprs_es), corum_complexes_Subset$ComplexID),])$ComplexName

mod1 <- model.matrix(~ factor(aml_noaml_exprs_es$cancer_type) -1)
colnames(mod1) <- c("AML", "NonAML")
fit2 <- lmFit(aml_noaml_exprs_es, mod1)
contrast.matrix <- makeContrasts("AML-NonAML", levels = mod1)
contrast.matrix
fit2C <- contrasts.fit(fit2, contrast.matrix)
fit2C <- eBayes(fit2C)
topTable(fit2C)
res <- decideTests(fit2C, p.value=0.10)
summary(res)

tt <- topTable(fit2C, n=Inf)

##enhanced valcono plot

corum_complexes_names_matched = (corum_complexes_Subset[match(rownames(tt), corum_complexes_Subset$ComplexID),])$ComplexName
tt$ComplexName = corum_complexes_names_matched
tt$ComplexID = rownames(tt)

## --------------------------
## Fix inputs for volcano plot
## --------------------------
## fold-change
fc1 = tt$logFC
head(fc1)
max(fc1)
min(fc1)
fc1[fc1 < -0.2] <- -0.2
fc1[fc1 > 0.2] <- 0.2

# adjusted p-value
qval1 = tt$adj.P.Val
max(qval1)
min(qval1)
qval1[qval1 < 0.01] <- 0.01

# gene names
genes1 <- make.unique(rownames(tt))

# Create a data frame using gene names, fold changes, and q-values

df1 <- data.frame(gene_name = as.character(genes1), lfc = fc1, q = qval1, complexName = tt$ComplexName, stringsAsFactors = F)
rownames(df1) <- genes1
head(df1, 20)

## Plot volcano
# y-axis upper limit
max_pval1 <- -log10(min(qval1)) + 0.1

labels1 = (tt$ComplexName)[1:15]
options(ggrepel.max.overlaps = Inf)

plt1 <- EnhancedVolcano(df1, x = 'lfc', y = 'q', lab = df1$complexName,
                        pCutoff = 0.1, FCcutoff = 1.0,
                        gridlines.major = FALSE, gridlines.minor = FALSE, 
                        drawConnectors = TRUE, legendLabSize = 12, colAlpha = 1.0,
                        cutoffLineCol = "black", cutoffLineType = "dashed",
                        border = "full",
                        ##colCustom = cols1,
                        legendPosition = "right",
                        pointSize = 4, cutoffLineWidth = 0.4,
                        labFace = "plain", subtitle = "",
                        ylim = c(0, max_pval1), 
                        xlim = c(-0.1, 0.1),
                        axisLabSize = 12, captionLabSize = 12,
                        xlab = "SSGSEA enrichment score difference", ylab = "-Log10 adjusted p-value", title = "",
                        caption = paste0('Total = ', nrow(df1), ' genes'),typeConnectors = "closed", legendIconSize = 2,
                        selectLab = labels1,
                        borderWidth = 1.5, boxedLabels = FALSE,
                        widthConnectors = 0.75)
print(plt1)
ggsave(plt1, filename = file.path(".",dir_out,"SSGSEA_aml_nonmal_Volcano_redo.pdf"), height = 9, width = 16)

### enhanced heatmap colored up and down regulated genes
DEAdf = df1
DEAdf$q[360:369] = tt$adj.P.Val[360:369]
keyvals1 <- ifelse(DEAdf$lfc < 0.0 & DEAdf$q <= 0.1, 'red2', ifelse(DEAdf$lfc > 0.0 & DEAdf$q <= 0.1, 'royalblue', 'grey50'))

table(is.na(keyvals1)) 
levels(as.factor(keyvals1))
names(keyvals1)[keyvals1 == 'red2'] <- 'Down-Regulated'
names(keyvals1)[keyvals1 == 'grey50'] <- 'NS'
names(keyvals1)[keyvals1 == 'royalblue'] <- 'Up-Regulated'


plt2 <- EnhancedVolcano(DEAdf,
                x = 'lfc', y = 'q',
                lab = DEAdf$complexName,
                pCutoff = 0.1, FCcutoff = 1.0,
                #xlim = c(min(DEAdf$lfc, na.rm=TRUE), max(DEAdf$lfc, na.rm=TRUE)),
                ylim=c(0, max(-log10(DEAdf$q), na.rm=TRUE) + 0.1),
                pointSize = 5.0,
                labSize = 6.0,
                colCustom = keyvals1,
                selectLab = DEAdf$complexName[which(names(keyvals1) %in% c('Down-Regulated', 'Up-Regulated'))],
                drawConnectors = TRUE, 
                legendLabSize = 16, 
                colAlpha = 1.0,
                gridlines.major = FALSE, gridlines.minor = FALSE,
                legendPosition = "right",
                labFace = "plain", subtitle = "",
                captionLabSize = 12,
                xlim = c(-0.1, 0.1),
                xlab = "SSGSEA enrichment score difference", ylab = "-Log10 adjusted p-value", title = "",
                cutoffLineWidth = 0.8,
                cutoffLineCol = "red", cutoffLineType = "dashed"
                
                )
                

ggsave(plt2, filename = file.path(".",dir_out,"SSGSEA_Volcano_mod_connect.pdf"), height = 16, width = 24)

plt2 <- EnhancedVolcano(DEAdf,
                        x = 'lfc', y = 'q',
                        lab = DEAdf$complexName,
                        pCutoff = 0.1, FCcutoff = 1.0,
                        #xlim = c(min(DEAdf$lfc, na.rm=TRUE), max(DEAdf$lfc, na.rm=TRUE)),
                        ylim=c(0, max(-log10(DEAdf$q), na.rm=TRUE) + 0.1),
                        pointSize = 5.0,
                        labSize = 6.0,
                        colCustom = keyvals1,
                        selectLab = labels1,
                        drawConnectors = TRUE, 
                        legendLabSize = 16, 
                        colAlpha = 1.0,
                        gridlines.major = FALSE, gridlines.minor = FALSE,
                        legendPosition = "right",
                        labFace = "plain", subtitle = "",
                        captionLabSize = 12,
                        xlim = c(-0.1, 0.1),
                        xlab = "SSGSEA enrichment score difference", ylab = "-Log10 adjusted p-value", title = "",
                        cutoffLineWidth = 0.8,
                        cutoffLineCol = "red", cutoffLineType = "dashed",
                       
                        
)


ggsave(plt2, filename = file.path(".",dir_out,"SSGSEA_Volcano_mod_connect_labels.pdf"), height = 16, width = 24)

#####################

## select names

DEAdf_1.5 = DEAdf
# adjusted p-value
qval1.5 = DEAdf_1.5$q
max(qval1.5)
min(qval1.5)
qval1.5[qval1.5 < 0.03] <- 0.03

DEAdf_1.5$q = qval1.5

labels2 = c(labels1[9], labels1[10], labels1[14])


plt3 <- EnhancedVolcano(DEAdf_1.5,
                        x = 'lfc', y = 'q',
                        lab = DEAdf$complexName,
                        pCutoff = 0.1, FCcutoff = 1.0,
                        #xlim = c(min(DEAdf$lfc, na.rm=TRUE), max(DEAdf$lfc, na.rm=TRUE)),
                        ylim=c(0, max(-log10(DEAdf_1.5$q), na.rm=TRUE) + 0.1),
                        pointSize = 6.0,
                        labSize = 6.0,
                        colCustom = keyvals1,
                        selectLab = labels2,
                        drawConnectors = TRUE, 
                        legendLabSize = 16, 
                        colAlpha = 1.0,
                        gridlines.major = FALSE, gridlines.minor = FALSE,
                        legendPosition = "right",
                        labFace = "plain", subtitle = "",
                        captionLabSize = 12,
                        xlim = c(-0.1, 0.1),
                        xlab = "SSGSEA enrichment score difference", ylab = "-Log10 adjusted p-value", title = "",
                        cutoffLineWidth = 0.8,
                        cutoffLineCol = "red", cutoffLineType = "dashed",
                        
                        
)

print(plt3)

ggsave(plt3, filename = file.path(".",dir_out,"SSGSEA_Volcano_mod_connect_labels_1.5.pdf"), height = 8, width = 12)


