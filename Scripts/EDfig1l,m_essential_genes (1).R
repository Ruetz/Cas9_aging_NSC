library(VennDiagram)
library(readxl)
library(tidyr)
library(ggplot2)

#generate a list of D14 depleted genes from screens that have an FDR value <0.1 in 2/3 screens. Generate for both young and old, and then combine into one list. 
library(readxl)
library(xlsx)
library(tidyverse)

# Import data from excel sheet
d14_scores_from_supp1 <- read_excel('/Users/angiepogson/Dropbox/Screen\ paper\ Tyson/Nature_revisions/Figures/Older\ figures/Angie-figures/EDfig1l,m/Supplementary_Table1_Gene_Scores_Invitro_Screens.xlsx',
                         sheet = "Gene_Scores_d14")

#generate dataframe of list of genes for young1 negative Young1_casTLE_Effect_d14
d14_scores_y1_depleted <- d14_scores_from_supp1[d14_scores_from_supp1$`Young1 casTLE Effect` < 0,]
d14_scores_y1_depleted <- d14_scores_y1_depleted[d14_scores_y1_depleted$`Young1 FDR Q-value` < 0.1,][2]

#generate dataframe of list of genes for young2 negative Young2_casTLE_Effect_d14
d14_scores_y2_depleted <- d14_scores_from_supp1[d14_scores_from_supp1$`Young2 casTLE Effect` < 0,]
d14_scores_y2_depleted <- d14_scores_y2_depleted[d14_scores_y2_depleted$`Young2 FDR Q-value` < 0.1,][2]

#generate dataframe of list of genes for young3 negative Young3_casTLE_Effect_d14
d14_scores_y3_depleted <- d14_scores_from_supp1[d14_scores_from_supp1$`Young3 casTLE Effect` < 0,]
d14_scores_y3_depleted <- d14_scores_y3_depleted[d14_scores_y3_depleted$`Young3 FDR  Q-value` < 0.1,][2]

#keep genes that appear in both d14_scores_y1_depleted and d14_scores_y2_depleted
comb_d14_y1_y2_depeleted_list <- merge(d14_scores_y1_depleted, d14_scores_y2_depleted, by = 'Symbol')

#keep genes that appear in both d14_scores_y2_depleted and d14_scores_y3_depleted
comb_d14_y2_y3_depeleted_list <- merge(d14_scores_y2_depleted, d14_scores_y3_depleted, by = 'Symbol')

#keep genes that appear in both d14_scores_y1_depleted and d14_scores_y3_depleted
comb_d14_y1_y3_depeleted_list <- merge(d14_scores_y1_depleted, d14_scores_y3_depleted, by = 'Symbol')

#combine all 3 sets of dataframes. I want to append all of them to each other and then remove duplicates
depleted_y_d14 = rbind(comb_d14_y1_y2_depeleted_list,comb_d14_y1_y3_depeleted_list)
depleted_y_d14 = rbind(depleted_y_d14,comb_d14_y2_y3_depeleted_list)
depleted_y_d14 = unique(depleted_y_d14)

#Now, I want to do the same thing for the old screens. 
#generate dataframe of list of genes for old1 negative Old1_casTLE_Effect_d14
d14_scores_o1_depleted <- d14_scores_from_supp1[d14_scores_from_supp1$`Old1 casTLE Effect` < 0,]
d14_scores_o1_depleted <- d14_scores_o1_depleted[d14_scores_o1_depleted$`Old1 FDR  Q-value` < 0.1,][2]

#generate dataframe of list of genes for young2 negative Young2_casTLE_Effect_d14
d14_scores_o2_depleted <- d14_scores_from_supp1[d14_scores_from_supp1$`Old2 casTLE Effect` < 0,]
d14_scores_o2_depleted <- d14_scores_o2_depleted[d14_scores_o2_depleted$`Old2 FDR Q-value` < 0.1,][2]

#generate dataframe of list of genes for young3 negative Young3_casTLE_Effect_d14
d14_scores_o3_depleted <- d14_scores_from_supp1[d14_scores_from_supp1$`Old3 casTLE Effect` < 0,]
d14_scores_o3_depleted <- d14_scores_o3_depleted[d14_scores_o3_depleted$`Old3 FDR  Q-value` < 0.1,][2]

#keep genes that appear in both d14_scores_y1_depleted and d14_scores_y2_depleted
comb_d14_o1_o2_depeleted_list <- merge(d14_scores_o1_depleted, d14_scores_o2_depleted, by = 'Symbol')

#keep genes that appear in both d14_scores_y2_depleted and d14_scores_y3_depleted
comb_d14_o2_o3_depeleted_list <- merge(d14_scores_o2_depleted, d14_scores_o3_depleted, by = 'Symbol')

#keep genes that appear in both d14_scores_y1_depleted and d14_scores_y3_depleted
comb_d14_o1_o3_depeleted_list <- merge(d14_scores_o1_depleted, d14_scores_o3_depleted, by = 'Symbol')

#Concatenate all 3 sets of dataframes for the old d14 depleted genes. Remove duplicates
depleted_o_d14 = rbind(comb_d14_o1_o2_depeleted_list,comb_d14_o1_o3_depeleted_list)
depleted_o_d14 = rbind(depleted_o_d14,comb_d14_o2_o3_depeleted_list)
depleted_o_d14 = unique(depleted_o_d14)

#combine young and old essential gene lists
depleted_d14_list = rbind(depleted_o_d14, depleted_y_d14)
depleted_d14_list = unique(depleted_d14_list)

write.xlsx(depleted_d14_list, file = "/Users/angiepogson/Downloads/depeleted_d14.xlsx")

#List of core essential genes 2 from Hart, et al 2017 (Supplementary Table 2)
CEG2 <- read_excel('/Users/brunetlab/Dropbox/Screen\ paper\ Tyson/Nature_revisions/Figures/Older\ figures/Angie-figures/EDfig1l,m/CEG2.xlsx')

#List of essential genes downloaded from Online GEne Essentiality Database from https://v3.ogee.info/#/home (Gurumayum, et al 2021)
OGEE_CSEG_CEG <- read_excel('/Users/brunetlab/Dropbox/Screen\ paper\ Tyson/Nature_revisions/Figures/Older\ figures/Angie-figures/EDfig1l,m/OGEE.xlsx', sheet = "CEG and CSEG")

#Total depleted, overlap with CEG2
venn.diagram(
  x = list(unique(depleted_d14_list$Symbol), unique(CEG2$Symbol)),
  category.names = c("Depeleted Day 14 Genes","CEG2, Hart 2017"),
  filename = 'total_depleted_CEG2_venn-102823.jpg',
  output=TRUE,
  
  lwd = 2,
  fill = c(alpha("seagreen", 0.9),alpha("darkorange1", 0.8)),
  
  cat.pos = c(-30, 30),
  cat.dist = c(0.055,0.085),
  
  cex = 1.6,
  fontface = "bold",
  fontfamily = "sans"
)

#Total depleted, overlap with OGEE - from the CEGs and CSEGs database
venn.diagram(
  x = list(unique(depleted_d14_list$Symbol), unique(OGEE_CSEG_CEG$gene)),
  category.names = c("Total depleted genes","OGEE, Gurumayam 2020"),
  filename = 'total_depleted_OGEE_venn.jpg',
  output=TRUE,
  
  lwd = 2,
  fill = c(alpha("seagreen", 0.6),alpha("gray32", 0.6)),
  
  cat.pos = c(-30, 30),
  cat.dist = c(0.055,0.085),
  
  cex = 1.6,
  fontface = "bold",
  fontfamily = "sans"
)

##########################################################################################
#Is the overlap between these two datasets significant?
##########################################################################################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GeneOverlap")

library(GeneOverlap)

#also need to find the total number of genes included in the dataset. In supp table 1, this is 22906 (the total number of genes listed in the d14 scores tab)
#For OGEE
data(GeneOverlap)
obj.OGEE <- newGeneOverlap(depleted_d14_list$Symbol,
                         OGEE_CSEG_CEG$gene,
                         22905)
obj.OGEE <-testGeneOverlap(obj.OGEE)
print(obj.OGEE)
getContbl(obj.OGEE)

#Same thing for CEG2
obj.CEG2 <- newGeneOverlap(depleted_d14_list$Symbol,
                           CEG2$Symbol,
                           22905)
obj.CEG2 <-testGeneOverlap(obj.CEG2)
print(obj.CEG2)
getContbl(obj.CEG2)

