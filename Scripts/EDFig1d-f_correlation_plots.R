# Load necessary packages
library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(gghighlight)
library(ggpubr)

#Purpose is to generate 3 scatterplots that compare the directionality and size of the casTLE effects per each screen replicate in the Young Day 14 screen.

# Import day 14 CasTLE scores from Supplementary Table 1
d14_scores_from_supp1 <- read_excel('/Users/brunetlab/Dropbox/Screen\ paper\ Tyson/Nature_revisions/Figures/Older\ figures/Angie-figures/EDfig1d-f/Supplementary_Table1_Gene_Scores_Invitro_Screens.xlsx',
                         sheet = "Gene_Scores_d14")

#CasTLE score is the magnitude of the gene knockout, whereas casTLE Effect is the directionality of the effect (positive value is enriching in the screen, negative is depleting). 
#Calculate the casTLE score multiplied by the direction for each gene knockout for Young Day 14 replicate 1
d14_scores_from_supp1$young1_direction_score_d14 <- d14_scores_from_supp1$`Young1 casTLE Score` * sign(d14_scores_from_supp1$`Young1 casTLE Effect`)
#Calculate the casTLE score multiplied by the direction for each gene knockout for Young Day 14 replicate 2
d14_scores_from_supp1$young2_direction_score_d14 <- d14_scores_from_supp1$`Young2 casTLE Score` * sign(d14_scores_from_supp1$`Young2 casTLE Effect`)
#Calculate the casTLE score multiplied by the direction for each gene knockout for Young Day 14 replicate 3
d14_scores_from_supp1$young3_direction_score_d14 <- d14_scores_from_supp1$`Young3 casTLE Score` * sign(d14_scores_from_supp1$`Young3 casTLE Effect`)

#Generate a plot that compares the direction and magnitude of the score of the gene knockout for pairs of each of the screens. 
#Additionally, calculate a Pearson's correlation coefficient (r) and add a regression line
d14_y1_y2_score <- ggplot(d14_scores_from_supp1, aes(x = young1_direction_score_d14, y = young2_direction_score_d14)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, level = 0.99) +
  stat_cor(aes(label=..rr.label..), label.x=70, label.y=78) +
  xlab("Young_rep1 casTLE Score") +
  ylab("Young_rep2 casTLE Score") +
  ggtitle("Day 14: Young_rep1 v Young_rep2 Scores")+
  scale_x_continuous(limits = c(-60, 80)) +
  scale_y_continuous(limits = c(-60, 80)) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) 
d14_y1_y2_score


d14_y2_y3_score <- ggplot(d14_scores, aes(x = young2_direction_score_d14, y = young3_direction_score_d14)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, level = 0.99) +
  stat_cor(aes(label=..rr.label..), label.x=70, label.y=78) +
  xlab("Young_rep2 casTLE Score") +
  ylab("Young_rep3 casTLE Score") +
  ggtitle("Day 14: Young_rep2 v Young_rep3 Scores")  +
  scale_x_continuous(limits = c(-60, 80)) +
  scale_y_continuous(limits = c(-60, 80)) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) 
d14_y2_y3_score

d14_y1_y3_score <- ggplot(d14_scores, aes(x = young1_direction_score_d14, y = young3_direction_score_d14)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", se = TRUE, level = 0.99) +
  stat_cor(aes(label=..rr.label..), label.x=70, label.y=78, cor.coef.name = c("rho")) +
  xlab("Young_rep1 casTLE Score") +
  ylab("Young_rep3 casTLE Score") +
  ggtitle("Day 14: Young_rep1 v Young_rep3 Scores") +
  scale_x_continuous(limits = c(-60, 80)) +
  scale_y_continuous(limits = c(-60, 80)) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) 
d14_y1_y3_score
