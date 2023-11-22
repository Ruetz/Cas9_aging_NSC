library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)


filename <- "data/fao_signature_aiderus2018_19genes"
gene_list_GO <- readLines(paste0(filename,".csv"))
genelist <- tolower(c(gene_list_GO))




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
adhesion_data <- d[, colnames(d) %in% genelist]
meta <- meta[, colnames(meta) %in% c("S.Score", "G2M.Score", "Phase","AgeCond", "Replicate", "Celltype.LowRes")]
adhesion_data$adhesion_response <- rowSums(adhesion_data, na.rm=TRUE)
adhesion_data <- cbind(meta, adhesion_data)

# Reorder factors
print(unique(adhesion_data$Celltype.LowRes))
print(unique(adhesion_data$AgeCond))
CELLS <- c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast",
           "OPC", "Oligodendro", "Endothelial")
adhesion_data <- adhesion_data %>% filter(Celltype.LowRes %in% CELLS)
adhesion_data$celltype <- factor(adhesion_data$Celltype.LowRes,  levels=CELLS, ordered=T)
adhesion_data$age <- factor(adhesion_data$AgeCond, levels=c("Y_Control", "O_Control"), ordered=T)



s=0.4
a1=0.4

# Plot parameters
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                    symbols = c("****", "***", "**", "*", "ns"))

ageColors <- c("lightskyblue", "dodgerblue")

p <- ggplot(data=dplyr::filter(adhesion_data, celltype == "Astrocyte_qNSC"), aes(x=age, y=adhesion_response)) +
  geom_violin(aes(fill=age), trim=T,scale="width", draw_quantiles = c(.5)) +
  geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.35, h = 0)) +
  scale_fill_manual(values=ageColors) +
  facet_wrap(~celltype, scales = "free", nrow =2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(legend.position="none") +
  theme_bw() + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(method="wilcox.test")
p

ggsave("plots/Ast_qNSC_violinplot_FAO_Wilcoxon_Exercise.pdf",p, height=6, width=6)



p <- ggplot(data=dplyr::filter(adhesion_data, celltype == "Astrocyte_qNSC"), aes(x=age, y=adhesion_response)) +
  geom_violin(aes(fill=age), trim=T,scale="width", draw_quantiles = c(.5)) +
  geom_point(size=s, alpha=a1, color="black", position = position_jitter(w = 0.35, h = 0)) +
  scale_fill_manual(values=ageColors) +
  facet_wrap(~celltype, scales = "free", nrow =2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(legend.position="none") +
  theme_bw() + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(method="t.test")
p

ggsave("plots/Ast_qNSC_violinplot_FAO_Welch_t_test_Exercise.pdf",p, height=6, width=6)

