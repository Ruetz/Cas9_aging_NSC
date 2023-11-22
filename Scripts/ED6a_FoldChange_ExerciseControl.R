library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)

# Load Data
svz <- readRDS("data/Exercise_seurat.SVZ.annotated.2020-04-27.rds")
svz <- subset(svz, AgeCond == "O_Control"  | AgeCond == "Y_Control")
meta <- svz[[]]
d <- t(as.matrix(svz[['RNA']]@counts))

rm(svz)
gc()

# Normalize counts to expression values
d <- sweep(d, MARGIN = 1, FUN = "/", STATS = rowSums(d))
d <- log1p(d * 10000)
d <- as.data.frame(d)

# Subset data to signature genes
colnames(d) <- tolower(colnames(d))

# Reorder factors
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Neuronal",
           "OPC", "Oligodendro", "Endothelial", "Ependymal", "Mural",
           "Microglia", "Macrophage", "Doublet", "T_cell", "Vascular_Leptomeningeal", "Epithelial")
d$celltype <- factor(meta$Celltype.LowRes,  levels=CELLS, ordered=T)
d$age <- factor(meta$AgeCond, levels=c("Y_Control", "O_Control"), ordered=T)

# ADD ZEROS FOR SLC2A2 (since not found)
#d$slc2a2 <- 0


# subset genes and age
genes <- c("slc2a1", "slc2a3", "slc2a4", "slc2a6", "slc2a8", "slc2a10", "slc2a12", "slc2a13") # "slc2a2" not found

d_old <- d %>% dplyr::filter(celltype=="Astrocyte_qNSC") %>% dplyr::filter(age=="O_Control")
d_young <- d %>% dplyr::filter(celltype=="Astrocyte_qNSC") %>% dplyr::filter(age=="Y_Control")

# compute mean expression
old_exp <- colMeans(d_old %>% dplyr::select(genes))
young_exp <- colMeans(d_young %>% dplyr::select(genes))

# compute log2 FC
log2fcs <- log2(old_exp/young_exp)
df <- data.frame(gene=genes,
                 log2fc=log2fcs)
df$gene <- factor(df$gene, levels=genes, ordered=T)


# plot
p <- ggplot(df, aes(x=gene, y=log2fc))+
  geom_bar(stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p
ggsave("plots/Ast_qNSC_log2FC_slc2ax_Exercise.pdf",p, height=6, width=6)
