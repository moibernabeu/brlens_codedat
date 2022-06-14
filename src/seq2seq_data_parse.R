# Organising output data
# Moisès Bernabeu
# Barcelona, April 2022

library(tidyr)
library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

# Definitions ----
get_other <- function(x, ref) {
  if (x['from_sp'] != ref) {
    y <- x['from_sp']
  } else {
    y <- x['to_sp']
  }
  return(y)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Import and save image of data ----

# Read raw data
ryedat <- read.csv('../../02_get_distances/outputs/0005_dist.csv')
rhudat <- read.csv('../../02_get_distances/outputs/0076_dist.csv')

yespt <- read.csv('../../02_get_distances/outputs/0005_sptree_dist.csv')
huspt <- read.csv('../../02_get_distances/outputs/0076_sptree_dist.csv')

# Get species to
yedat <- ryedat
yedat$sp_to <- apply(ryedat, 1, get_other, ref = 'YEAST')
hudat <- rhudat
hudat$sp_to <- apply(rhudat, 1, get_other, ref = 'HUMAN')

yespt$sp_to <- apply(yespt, 1, get_other, ref = 'YEAST')
huspt$sp_to <- apply(huspt, 1, get_other, ref = 'HUMAN')

# Numbers description
count_df <- data.frame(Yeast = c(Total = length(which(yedat$sp_to != 'YEAST' &
                                            (yedat$to_sp == 'YEAST' |
                                               yedat$from_sp == 'YEAST') &
                                            yedat$dist != 0))),
           Human = c(total = length(which(hudat$sp_to != 'HUMAN' &
                                            (hudat$to_sp == 'HUMAN' |
                                               hudat$from_sp == 'HUMAN') &
                                            hudat$dist != 0))))

# Filter data
yedat <- yedat[which(yedat$mrca_type == 'S' & yedat$sp_to != 'YEAST' &
                       (yedat$to_sp == 'YEAST' | yedat$from_sp == 'YEAST') &
                       yedat$dist != 0), ]
hudat <- hudat[which(hudat$mrca_type == 'S' & hudat$sp_to != 'HUMAN' &
                       (hudat$to_sp == 'HUMAN' | hudat$from_sp == 'HUMAN') &
                       hudat$dist != 0), ]

count_df <- rbind(count_df, Speciation = c(dim(yedat)[1], dim(hudat)[1]))

yespt <- yespt[which(yespt$mrca_type == 'S' & yespt$sp_to != 'YEAST' &
                       (yespt$sp_to == 'YEAST' | yespt$from_sp == 'YEAST')), ]
huspt <- huspt[which(huspt$mrca_type == 'S' & huspt$sp_to != 'HUMAN' &
                       (huspt$sp_to == 'HUMAN' | huspt$from_sp == 'HUMAN')), ]

# Filter bad trees
yewdthratio <- yedat$whole_width / yedat$norm_width
yeq099 <- quantile(yewdthratio, 0.9)

huwdthratio <- hudat$whole_width / hudat$norm_width
huq099 <- quantile(huwdthratio, 0.9)

par(mfrow = c(1, 2))
hist(yewdthratio)
abline(v = yeq099, col = 'red', lty = 4)
hist(huwdthratio)
abline(v = huq099, col = 'red', lty = 4)
par(mfrow = c(1, 1))

length(which(yedat$whole_width / yedat$norm_width >= yeq099))
yedat <- yedat[-which(yedat$whole_width / yedat$norm_width >= yeq099), ]

length(which(hudat$whole_width / hudat$norm_width >= huq099))
hudat <- hudat[-which(hudat$whole_width / hudat$norm_width >= huq099), ]

count_df <- rbind(count_df, Width_ratio = c(dim(yedat)[1], dim(hudat)[1]))

# save(hudat, yedat, yespt, huspt, file = '../data/seed2sp_dist.Rdata')

# Descriptive plots
count_df$data <- row.names(count_df)
count_df <- gather(count_df, value = 'val', key = 'process', -data)

yedat_nms_sort <- names(sort(by(yedat$ndist_A, yedat$sp_to, median)))
yesortdf <- data.frame(table(yedat$sp_to)[yedat_nms_sort])

hudat_nms_sort <- names(sort(by(hudat$ndist_A, hudat$sp_to, median)))
husortdf <- data.frame(table(hudat$sp_to)[hudat_nms_sort])

filterpl <- ggplot(count_df, aes(x = reorder(data, -val), y = val, fill = process)) +
  geom_col(position = 'stack') +
  labs(fill = 'Phylome') +
  ylab('Number of trees') +
  xlab('') +
  geom_text(aes(label = val), position = position_stack(.5)) +
  theme(legend.position = 'left')

yespto_tr <- ggplot(yesortdf, aes(x = reorder(Var1, -Var1), y = Freq)) +
  geom_point(colour = gg_color_hue(2)[2]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  xlab('Species to') +
  ylab('Number of trees')

huspto_tr <- ggplot(husortdf, aes(x = reorder(Var1, -Var1), y = Freq)) +
  geom_point(colour = gg_color_hue(1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  xlab('Species to') +
  ylab('Number of trees')

pdf('../msct_plots/sumstats.pdf', width = 14, height = 6)
ggarrange(filterpl,
          ggarrange(yespto_tr, huspto_tr, align = 'hv', ncol = 1, labels = c('b', 'c')),
          widths = c(2, 1), labels = c('a', ''))
dev.off()


summary(yesortdf$Freq)
summary(husortdf$Freq)
