setwd("~/Documents/Projects/PamelaLab/DepMap_Data")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(poolr)
library(reshape2)
library(rstatix)
library(gginnards)
library(ggthemes)
library(ggeasy)
library(readxl)
library(dplyr)


## read Depmap Metadata
sample_info = read.csv("./Model.csv")
dim(sample_info)
head(sample_info)

## read CRISPR data from DepMap 24Q2
data_crispr = read.csv("./CRISPR_(DepMap_Public_24Q2+Score,_Chronos).csv")
dim(data_crispr)
data_crispr[1:6, 1:6]
head(colnames(data_crispr))
colnames(data_crispr)[1] = "ModelID"

# colnames(data_crispr) = stringr::str_extract(colnames(data_crispr), "[^..]+")
# head(colnames(data_crispr))

colnames_df = data.frame(names = colnames(data_crispr) )
##writexl::write_xlsx(colnames_df, path = "colnames_df.xlsx")

## read RNAi data from DepMap 24Q2
data_RNAi = read.csv("./RNAi_(Achilles+DRIVE+Marcotte,_DEMETER2).csv")
dim(data_RNAi)
data_RNAi[1:6, 1:6]
head(colnames(data_RNAi))
colnames(data_RNAi)[1] = "ModelID"

##match colnames of RNAi with colnames CRISPR
data_RNAi_crispr = data_RNAi %>%
  dplyr::select(ModelID, colnames(data_crispr)[colnames(data_crispr) %in% colnames(data_RNAi)])

## read mRNA data from DepMap 24Q2
data_mRNA = read.csv("./Batch_corrected_Expression_Public_24Q2.csv")
dim(data_mRNA)
head(colnames(data_mRNA))
colnames(data_mRNA)[1] = "ModelID"

##match colnames of RNAi with colnames CRISPR
data_mRNA_crispr = data_mRNA %>%
  dplyr::select(ModelID, colnames(data_crispr)[colnames(data_crispr) %in% colnames(data_mRNA)])


## sample info for PAAD
PAAD_sample_info = sample_info %>%
dplyr::filter(OncotreeCode %in% "PAAD")

writexl::write_xlsx(PAAD_sample_info, "PAAD_sample_info.xlsx")

PAAD_Primary_sample_info = PAAD_sample_info %>%
  dplyr::filter(PrimaryOrMetastasis %in% "Primary")

PAAD_Metastatic_sample_info = PAAD_sample_info %>%
  dplyr::filter(PrimaryOrMetastasis %in% "Metastatic")

PAAD_sample_info$ModelID

##Extract PAAD Crispr data 
class(data_crispr)
data_crispr_PAAD = data_crispr %>%
  dplyr::filter(ModelID %in% PAAD_sample_info$ModelID)

dim(data_crispr_PAAD)
write.csv(data_crispr_PAAD, "data_crispr_PAAD.csv")
df1 = data_crispr_PAAD[, 2:18436]
median(df1[, 60], na.rm = TRUE)

data_crispr_PAAD_median = rbind(df1, Medians = apply(df1, 2, median, na.rm = TRUE))
as.numeric(data_crispr_PAAD_median[46,])
data_crispr_PAAD_median_extract = data.frame("Gene" = colnames(data_crispr_PAAD_median),
           "Median" = as.numeric(data_crispr_PAAD_median[46,]))
data_crispr_PAAD_median_extract = data_crispr_PAAD_median_extract %>%
  arrange(Median)
data_crispr_PAAD_median_extract$Rank = seq_along(1:dim(data_crispr_PAAD_median_extract)[1])

write.csv(data_crispr_PAAD_median_extract, "data_crispr_PAAD_median_extract.csv") 



## read example data with beta score difference and rank
example_data = data_crispr_PAAD_median_extract
## Genes tp label 
Top5_neg = data_crispr_PAAD_median_extract$Gene[1:5]
Top5_pos = data_crispr_PAAD_median_extract$Gene[18431:18435]
Top5_neg = c( Top5_neg, "ID1")

example_data$grouping <- "Other"
example_data[example_data$Gene %in% Top5_neg, ]$grouping <- "TopNegative"
example_data[example_data$Gene %in% Top5_pos, ]$grouping <- "TopPositive"
class(example_data)
## plot
options(ggrepel.max.overlaps = Inf)
g <-
  ggplot(example_data, aes(x = Median, y = Rank, fill = grouping, label = Gene))+
  geom_point(colour="white", shape=21, size = 4, 
             aes(fill = factor(grouping))) + 
  scale_fill_manual(values=c( "gray", "yellow", "orange"))+
  theme_classic()+
  xlim(-5, 2.5)

g <- g + geom_point(
  pch = 21,
  colour = "black",
  fill = "grey",
  size = 4,
  alpha = 0.1
)

g <- g + geom_point(
  data = example_data[example_data$grouping == "TopNegative",],
  aes(x = Median, y = Rank),
  pch = 21,
  colour = "black",
  fill = "yellow",
  size = 4,
  alpha = 1.0,
  stroke = 1.0
)

g <- g + geom_point(
  data = example_data[example_data$grouping == "TopPositive",],
  aes(x = Median, y = Rank),
  pch = 21,
  colour = "black",
  fill = "orange",
  size = 4,
  alpha = 1.0,
  stroke = 1.0
)

g <- g +
  geom_text_repel(
    data = example_data[example_data$grouping == "TopNegative",],
    # force_pull   = 0, # do not pull toward data points
    # nudge_y      = 0.1,
    # direction    = "y",
    # ##angle        = 90,
    # hjust        = 0,
    # segment.size = 0.5,
    ##max.iter = 1e4, max.time = 1,
    ##min.segment.length = 1,
    ##size = 3.5,
    box.padding = 2.0,
    # point.padding = 0.5
    
  )

g <- g +
  geom_text_repel(
    data = example_data[example_data$grouping == "TopPositive",],
    # force_pull   = 1, # do not pull toward data points
    # nudge_y      = 0.05,
    # direction    = "x",
    # ##angle        = 90,
    # hjust        = 0,
    # segment.size = 0.5,
    # max.iter = 1e4, max.time = 1,
    # min.segment.length = 1,
    box.padding = 1
  )

g <- g + xlab("PAAD Median CRISPR Score")
g <- g + ylab("Rank")
g = g+theme(axis.text.x = element_text(face="bold", color="black", 
                                       size=14, angle=0),
            axis.text.y = element_text(face="bold", color="black", 
                                       size=14, angle=0))+
  theme(text = element_text(size = 20)) +
  theme(axis.ticks.length=unit(.30, "cm"))+
  theme(legend.position = c(0.3, 0.7),
        legend.direction = "vertical")
g = g+
  ggtitle("")+
  ggeasy::easy_center_title()


ggsave(plot=g, filename="./PAAD Median CRISPR Score_gene_ranking.pdf", units = "in", width=8, height= 8, dpi=600)

## Correlation plot mRNA expression vs CRISPR score for ID1 gene

## read mRNA data from DepMap 24Q2
data_mRNA = read.csv("./Batch_corrected_Expression_Public_24Q2.csv")
dim(data_mRNA)
head(colnames(data_mRNA))
colnames(data_mRNA)[1] = "ModelID"

##Extract PAAD mRNA data 
class(data_mRNA)
data_mRNA_PAAD = data_mRNA %>%
  dplyr::filter(ModelID %in% PAAD_sample_info$ModelID)

dim(data_crispr_PAAD)
dim(data_mRNA_PAAD)
write.csv(data_mRNA_PAAD, "data_mRNA_PAAD.csv")

## match data_mRNA_PAAD with data_crispr_PAAD
data_mRNA_PAAD_overlap_crispr = data_mRNA_PAAD %>%
dplyr::select(ModelID, colnames(data_crispr_PAAD)[colnames(data_crispr_PAAD) %in% colnames(data_mRNA_PAAD)])


data_mRNA_PAAD_overlap_crispr = data_mRNA_PAAD_overlap_crispr %>%
  dplyr::filter(ModelID %in% data_crispr_PAAD$ModelID)
dim(data_mRNA_PAAD_overlap_crispr)
dim(data_crispr_PAAD)

data_crispr_PAAD_overlap_mRNA = data_crispr_PAAD %>%
  dplyr::filter(ModelID %in% data_mRNA_PAAD_overlap_crispr$ModelID)
dim(data_crispr_PAAD_overlap_mRNA)

data_mRNA_PAAD_overlap_crispr = data_mRNA_PAAD_overlap_crispr %>%
  arrange(ModelID)

data_mRNA_PAAD_overlap_crispr_ID1 = data_mRNA_PAAD_overlap_crispr %>%
  dplyr::select(ModelID, ID1)

colnames(data_mRNA_PAAD_overlap_crispr_ID1)[2] = "ID1_mRNA"

data_crispr_PAAD_overlap_mRNA_ID1 = data_crispr_PAAD_overlap_mRNA %>%
  dplyr::select(ModelID, ID1)
   
colnames(data_crispr_PAAD_overlap_mRNA_ID1)[2] = "ID1_crispr"

corr_ID1_PAAD_mRNA_crispr = data_crispr_PAAD_overlap_mRNA_ID1 %>%
  left_join(data_mRNA_PAAD_overlap_crispr_ID1, by = "ModelID")

cbind(data_crispr_PAAD_overlap_mRNA$ModelID, data_mRNA_PAAD_overlap_crispr$ModelID)

writexl::write_xlsx(as.data.frame(corr_ID1_PAAD_mRNA_crispr), "corr_ID1_PAAD_mRNA_crispr.xlsx")

## correlation plot

# ## Read input data with LFC or Gene score or Z score for two comparisons or replicates
# subset_edl = readxl::read_excel("./data/subset_edl.xlsx", sheet = 1) %>%
#   as.data.frame()
# rownames(subset_edl) = subset_edl$Gene

## subset data of interest for labeling and highlighting genes
ID1_crispr_high = corr_ID1_PAAD_mRNA_crispr %>%
  arrange(desc(ID1_crispr)) %>%
  slice(1:5) 

ID1_crispr_low = corr_ID1_PAAD_mRNA_crispr %>%
  arrange(desc(ID1_crispr)) %>%
  slice(39:43)

ID1_crispr_map = rbind(ID1_crispr_high, ID1_crispr_low)

ID1_mRNA_high = corr_ID1_PAAD_mRNA_crispr %>%
  arrange(desc(ID1_mRNA)) %>%
  slice(1:5) 

ID1_mRNA_low = corr_ID1_PAAD_mRNA_crispr %>%
  arrange(desc(ID1_mRNA)) %>%
  slice(39:43)

ID1_mRNA_map = rbind(ID1_mRNA_high, ID1_mRNA_low)



## code to generate the plots in slide 2 and 3
options(ggrepel.max.overlaps = Inf)
p = ggplot(corr_ID1_PAAD_mRNA_crispr, aes(ID1_crispr, ID1_mRNA, label = ModelID)) +
  geom_text_repel(
    data          = ID1_crispr_map,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    ##direction    = "y",
    ##angle        = 90,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    min.segment.length = 0,
    box.padding = 0.5
  ) +
  geom_point(size = ifelse(corr_ID1_PAAD_mRNA_crispr$ModelID %in% ID1_crispr_map$ModelID, 2, 0.5),
             color = ifelse(corr_ID1_PAAD_mRNA_crispr$ModelID %in% ID1_crispr_map$ModelID, "red", "black"),
             alpha = ifelse(corr_ID1_PAAD_mRNA_crispr$ModelID %in% ID1_crispr_map$ModelID, 1, 0.1), stroke = 1)+
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 3)+
  geom_hline(yintercept = 2.5, linetype = 3)+
  geom_text_repel(
    data          = ID1_mRNA_map,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    ##direction    = "y",
    ##angle        = 90,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    min.segment.length = 0,
    box.padding = 1.0
  ) +
  geom_point(size = ifelse(corr_ID1_PAAD_mRNA_crispr$ModelID %in% ID1_mRNA_map$ModelID, 2, 0.5),
             color = ifelse(corr_ID1_PAAD_mRNA_crispr$ModelID %in% ID1_mRNA_map$ModelID, "#007575", "black"),
             alpha = ifelse(corr_ID1_PAAD_mRNA_crispr$ModelID %in% ID1_mRNA_map$ModelID, 1, 0.1), stroke = 1)+
  theme_classic()+
  xlab("PAAD_ID1_CRISPR") +
  ylab("PAAD_ID1_mRNA") +
  ggtitle("Correlation of PAAD_ID1_CRISPR and PAAD_ID1_mRNA")
 

ggsave(plot=p, filename="PAAD_ID1_CRISPR_mRNA_Correlation.pdf", units = "in", width=12, height=8)  

##Ranked Status of ID1 dependency in cancer 

dim(sample_info)
cancer_types = unique(sample_info$OncotreePrimaryDisease)
length(cancer_types)
cancer_types[53] = "Undifferentiated Pleomorphic Sarcoma"

ID1_2024_dep = data_crispr %>%
  dplyr::select(ModelID, ID1)
dim(ID1_2024_dep)
head(ID1_2024_dep)

list = lapply(1:89, function(l){
  ##l =1
  print(l)
  print(cancer_types[l])
  
  
  sample_info_cell_lines = sample_info %>%
    dplyr::filter(OncotreePrimaryDisease == cancer_types[l]) %>%
    pull(ModelID)
  print(sample_info_cell_lines)
  
  ID1_2024_dep_stats_cell_lines = ID1_2024_dep %>%
    dplyr::filter(ModelID %in% sample_info_cell_lines) %>%
    pull(ID1)
  print(ID1_2024_dep_stats_cell_lines)
  
  if (length(ID1_2024_dep_stats_cell_lines) == 0) {
    score = NA
  } else {
    score = ID1_2024_dep_stats_cell_lines
  }
  print(score)
  
  data.frame(score = score,
             cancer = print(cancer_types[l]))
  
}) %>%
  set_names(cancer_types) %>%
  bind_rows()
list

df = list
dim(df)
dim(data_crispr)

df = df %>% drop_na()
dim(df)
head(df)

new_order <- with(df, reorder(cancer , score, median , na.rm=T))

head(new_order)

pdf(file = "ID1_complete_box_Plot_primary_cancer.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 10)
par( mar=c(25,10, 2.3, 2))
# boxplot(score ~ cancer, data = df, col = "white",
#         main = "SAGAComplex score CRISPR Dataset",
#         ##xlab = "Cancer type Clinical STATUS",
#         ylab = "Avana gene CERES 24Q2 (SAGAComplexScore)",
#         yaxt = "n", las = 2)
# axis(2, at= seq(-1.25, 1, 0.25))
# mtext(side=1, text="Cancer type Clinical STATUS", line=9)

boxplot(score ~ new_order, data = df, col = "white",
        main = "ID1 DepMap CRISPR score Cancer Status",
        ##xlab = "Cancer type Clinical STATUS",
        ylab = "ID1 DepMap CRISPR score",
        yaxt = "n", las = 2,
        ylim = c(-1.25, 1))
axis(2, at= seq(-1.25, 1, 0.25))
mtext(side=1, text="Cancer type Clinical STATUS", line=9)




manualcolors<-c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")

length(manualcolors)
stripchart(score ~ new_order,
           data = df,
           method = "jitter",
           pch = 19,
           col = manualcolors,
           vertical = TRUE,
           add = TRUE)


dev.off()


##Ranked Status of ID3 dependency in cancer 

dim(sample_info)
cancer_types = unique(sample_info$OncotreePrimaryDisease)
length(cancer_types)
cancer_types[53] = "Undifferentiated Pleomorphic Sarcoma"

ID3_2024_dep = data_crispr %>%
  dplyr::select(ModelID, ID3)
dim(ID3_2024_dep)
head(ID3_2024_dep)

list = lapply(1:89, function(l){
  ##l =1
  print(l)
  print(cancer_types[l])
  
  
  sample_info_cell_lines = sample_info %>%
    dplyr::filter(OncotreePrimaryDisease == cancer_types[l]) %>%
    pull(ModelID)
  print(sample_info_cell_lines)
  
  ID3_2024_dep_stats_cell_lines = ID3_2024_dep %>%
    dplyr::filter(ModelID %in% sample_info_cell_lines) %>%
    pull(ID3)
  print(ID3_2024_dep_stats_cell_lines)
  
  if (length(ID3_2024_dep_stats_cell_lines) == 0) {
    score = NA
  } else {
    score = ID3_2024_dep_stats_cell_lines
  }
  print(score)
  
  data.frame(score = score,
             cancer = print(cancer_types[l]))
  
}) %>%
  set_names(cancer_types) %>%
  bind_rows()
list

df = list
dim(df)
dim(data_crispr)

df = df %>% drop_na()
dim(df)
head(df)

new_order <- with(df, reorder(cancer , score, median , na.rm=T))

head(new_order)

pdf(file = "ID3_complete_box_Plot_primary_cancer.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 10)
par( mar=c(25,10, 2.3, 2))
# boxplot(score ~ cancer, data = df, col = "white",
#         main = "SAGAComplex score CRISPR Dataset",
#         ##xlab = "Cancer type Clinical STATUS",
#         ylab = "Avana gene CERES 24Q2 (SAGAComplexScore)",
#         yaxt = "n", las = 2)
# axis(2, at= seq(-1.25, 1, 0.25))
# mtext(side=1, text="Cancer type Clinical STATUS", line=9)

boxplot(score ~ new_order, data = df, col = "white",
        main = "ID3 DepMap CRISPR score Cancer Status",
        ##xlab = "Cancer type Clinical STATUS",
        ylab = "ID3 DepMap CRISPR score",
        yaxt = "n", las = 2,
        ylim = c(-1.25, 1))
axis(2, at= seq(-1.25, 1, 0.25))
mtext(side=1, text="Cancer type Clinical STATUS", line=9)




manualcolors<-c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                'magenta', 'darkolivegreen4', 'indianred1', 'tan4', 'darkblue', 
                'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 
                'seagreen1', 'moccasin', 'mediumvioletred', 'seagreen','cadetblue1',
                "darkolivegreen1" ,"tan2" ,   "tomato3" , "#7CE3D8","gainsboro")

length(manualcolors)
stripchart(score ~ new_order,
           data = df,
           method = "jitter",
           pch = 19,
           col = manualcolors,
           vertical = TRUE,
           add = TRUE)


dev.off()

## correlation plot ID3 mRNA and CRISPR

data_mRNA_PAAD_overlap_crispr_ID3 = data_mRNA_PAAD_overlap_crispr %>%
  dplyr::select(ModelID, ID3)

colnames(data_mRNA_PAAD_overlap_crispr_ID3)[2] = "ID3_mRNA"

data_crispr_PAAD_overlap_mRNA_ID3 = data_crispr_PAAD_overlap_mRNA %>%
  dplyr::select(ModelID, ID3)

colnames(data_crispr_PAAD_overlap_mRNA_ID3)[2] = "ID3_crispr"

corr_ID3_PAAD_mRNA_crispr = data_crispr_PAAD_overlap_mRNA_ID3 %>%
  left_join(data_mRNA_PAAD_overlap_crispr_ID3, by = "ModelID")



writexl::write_xlsx(as.data.frame(corr_ID3_PAAD_mRNA_crispr), "corr_ID3_PAAD_mRNA_crispr.xlsx")

## correlation plot

# ## Read input data with LFC or Gene score or Z score for two comparisons or replicates
# subset_edl = readxl::read_excel("./data/subset_edl.xlsx", sheet = 1) %>%
#   as.data.frame()
# rownames(subset_edl) = subset_edl$Gene

## subset data of interest for labeling and highlighting genes
ID3_crispr_high = corr_ID3_PAAD_mRNA_crispr %>%
  arrange(desc(ID3_crispr)) %>%
  slice(1:5) 

ID3_crispr_low = corr_ID3_PAAD_mRNA_crispr %>%
  arrange(desc(ID3_crispr)) %>%
  slice(39:43)

ID3_crispr_map = rbind(ID3_crispr_high, ID3_crispr_low)

ID3_mRNA_high = corr_ID3_PAAD_mRNA_crispr %>%
  arrange(desc(ID3_mRNA)) %>%
  slice(1:5) 

ID3_mRNA_low = corr_ID3_PAAD_mRNA_crispr %>%
  arrange(desc(ID3_mRNA)) %>%
  slice(39:43)

ID3_mRNA_map = rbind(ID3_mRNA_high, ID3_mRNA_low)



## code to generate the plots in slide 2 and 3
options(ggrepel.max.overlaps = Inf)
p = ggplot(corr_ID3_PAAD_mRNA_crispr, aes(ID3_crispr, ID3_mRNA, label = ModelID)) +
  geom_text_repel(
    data          = ID3_crispr_map,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    ##direction    = "y",
    ##angle        = 90,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    min.segment.length = 0,
    box.padding = 0.5
  ) +
  geom_point(size = ifelse(corr_ID3_PAAD_mRNA_crispr$ModelID %in% ID3_crispr_map$ModelID, 2, 0.5),
             color = ifelse(corr_ID3_PAAD_mRNA_crispr$ModelID %in% ID3_crispr_map$ModelID, "red", "black"),
             alpha = ifelse(corr_ID3_PAAD_mRNA_crispr$ModelID %in% ID3_crispr_map$ModelID, 1, 0.1), stroke = 1)+
  theme_classic()+
  geom_vline(xintercept = 0, linetype = 3)+
  geom_hline(yintercept = 2.5, linetype = 3)+
  geom_text_repel(
    data          = ID3_mRNA_map,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    ##direction    = "y",
    ##angle        = 90,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    min.segment.length = 0,
    box.padding = 1.0
  ) +
  geom_point(size = ifelse(corr_ID3_PAAD_mRNA_crispr$ModelID %in% ID3_mRNA_map$ModelID, 2, 0.5),
             color = ifelse(corr_ID3_PAAD_mRNA_crispr$ModelID %in% ID3_mRNA_map$ModelID, "#007575", "black"),
             alpha = ifelse(corr_ID3_PAAD_mRNA_crispr$ModelID %in% ID3_mRNA_map$ModelID, 1, 0.1), stroke = 1)+
  theme_classic()+
  xlab("PAAD_ID3_CRISPR") +
  ylab("PAAD_ID3_mRNA") +
  ggtitle("Correlation of PAAD_ID3_CRISPR and PAAD_ID3_mRNA")


ggsave(plot=p, filename="PAAD_ID3_CRISPR_mRNA_Correlation.pdf", units = "in", width=12, height=8)  
