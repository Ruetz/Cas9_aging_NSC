library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)

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



### Combined Plot

# Map to "Other cells"

d$celltype_other <- plyr::mapvalues(d$celltype, c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Neuronal",
                                                          "OPC", "Oligodendro", "Endothelial", "Ependymal", "Mural",
                                                          "Microglia", "Macrophage", "Doublet", "T_cell", "Vascular_Leptomeningeal", "Epithelial"),
                                       c("Astrocyte_qNSC", "aNSC_NPC", "Neuroblast", "Other",
                                         "Other", "Other", "Other", "Other", "Other",
                                         "Other", "Other", "Other", "Other", "Other", "Other")) 
d <- d[!is.na(d$celltype_other),]


############################################################
## 1. Create your workhorse function for draw_group
## basically a copy paste of .subset2(GeomViolin, "draw_group") 
## where all stuff related to the violin is removed
draw_group_violin_quantiles <- function(self, data, quantiles, ..., linesize = NULL, 
                                        flipped_aes = FALSE) {
  data <- flip_data(data, flipped_aes)
  data <- transform(data, xminv = x - violinwidth * (x - xmin), 
                    xmaxv = x + violinwidth * (xmax - x))
  if (!(all(quantiles >= 0) && all(quantiles <= 1))) {
    cli::cli_abort("{.arg quantiles} must be between 0 and 1")
  }
  quantiles <- ggplot2:::create_quantile_segment_frame(data, quantiles)
  aesthetics <- data[rep(1, nrow(quantiles)), 
                     setdiff(names(data), 
                             c("x", "y", "group")), drop = FALSE]
  aesthetics$alpha <- rep(1, nrow(quantiles))
  if (!is.null(linesize)) {
    aesthetics$linewidth <- rep(linesize, each = 2)
  }
  both <- vctrs::vec_cbind(quantiles, aesthetics)
  both <- both[!is.na(both$group), , drop = FALSE]
  both <- flip_data(both, flipped_aes)
  quantile_grob <- GeomPath$draw_panel(both, ...)
  ggplot2:::ggname("geom_violin_quantiles", quantile_grob)
}

## 2. Create your Geom proto object
## here we just inherit everything from `GeomViolin` and just replace the draw_group
GeomViolinQuantiles <- ggproto(
  "GeomViolinQuantiles",
  GeomViolin,
  draw_group = draw_group_violin_quantiles
)

## 3. Create the user facting geom_violin_qunatiles function
geom_violin_quantiles <- function(mapping = NULL,
                                  data = NULL,
                                  stat = "ydensity",
                                  position = "dodge",
                                  quantiles = NULL,
                                  ...,
                                  linesize = NULL,
                                  trim = TRUE, scale = "area", 
                                  na.rm = FALSE, orientation = NA, 
                                  show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomViolinQuantiles, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, quantiles = quantiles, 
                      na.rm = na.rm, orientation = orientation, linesize = linesize, ...))
}

############################################################


s=0.4
a1=0.4

# Plot parameters
ageColors <- c("deepskyblue2", "red2")
#ageColors <- c("lightblue", "lightcoral")
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                    symbols = c("****", "***", "**", "*", "ns"))

p <- ggplot(data=d, aes(x=celltype_other, y=slc2a4, fill=age)) +
  geom_point(size=s, alpha=a1, color="grey50", position = position_jitterdodge(dodge=0.9)) +
  geom_violin(aes(fill=age), trim=T,scale="width", draw_quantiles = c(.5), alpha=0.5, linewidth=0.3) +
  geom_violin_quantiles(quantiles = 0.5, trim=T,scale="width", linesize = 1) +
  scale_fill_manual(values=ageColors) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(legend.position="none") +
  theme_bw() + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(method="wilcox.test", label="p.format")
p

ggsave("plots/Violinplot_Slc2a4_Wilcoxon_Exercise.pdf",p, height=4, width=8)




p <- ggplot(data=d, aes(x=celltype_other, y=slc2a4, fill=age)) +
  geom_point(size=s, alpha=a1, color="grey50", position = position_jitterdodge(dodge=0.9)) +
  geom_violin(aes(fill=age), trim=T,scale="width", draw_quantiles = c(.5), alpha=0.5, linewidth=0.3) +
  geom_violin_quantiles(quantiles = 0.5, trim=T,scale="width", linesize = 1) +
  scale_fill_manual(values=ageColors) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(legend.position="none") +
  theme_bw() + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(method="t.test", label="p.format")
p

ggsave("plots/Violinplot_Slc2a4_Welch_t_test_Exercise.pdf",p, height=4, width=8)