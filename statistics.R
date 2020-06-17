#Code for statistics in Hammerschmith, Woodruff, and Phillips 2020 "Opposing directions of stage-specific body length change in a close relative of C. elegans"


library(reshape2)

#renamed "all_species_dauer_size.csv" to "all_species_dauer_size.tsv" for data sharing
dat <- read.table("all_species_dauer_size.tsv", sep="\t", header=T)



#convert to microns
dat$length <- dat$length*1000
dat$width <- dat$width*1000

#get l:w ratio

dat$l_over_w <- dat$length/dat$width

#melt data
dat_melt <- melt(dat,measure.vars=2:4)

#reorder and rename species levels

#re-order and re-name species levels
dat_melt$species <- factor(dat_melt$species, levels = c("briggsae","tropicalis","elegans","inopinata"))

levels(dat_melt$species)[levels(dat_melt$species)=="briggsae"] <- "C. briggsae"
levels(dat_melt$species)[levels(dat_melt$species)=="elegans"] <- "C. elegans"
levels(dat_melt$species)[levels(dat_melt$species)=="inopinata"] <- "C. inopinata"
levels(dat_melt$species)[levels(dat_melt$species)=="tropicalis"] <- "C. tropicalis"


dat$species <- factor(dat$species, levels = c("briggsae","tropicalis","elegans","inopinata"))

levels(dat$species)[levels(dat$species)=="briggsae"] <- "C. briggsae"
levels(dat$species)[levels(dat$species)=="elegans"] <- "C. elegans"
levels(dat$species)[levels(dat$species)=="inopinata"] <- "C. inopinata"
levels(dat$species)[levels(dat$species)=="tropicalis"] <- "C. tropicalis"

#stats

ino <- dat[dat$species == "C. inopinata",]
ele <- dat[dat$species == "C. elegans",]

summary(ino$length)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  326.0   399.0   432.0   429.5   462.0   524.0 
#

sd(ino$length)
#[1] 44.79935


summary(ele$length)

#> summary(ele$length)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  333.0   399.0   444.0   455.8   496.0   662.0 

sd(ele$length)

#[1] 65.21081

wilcox.test(ino$length,ele$length)

#> wilcox.test(ino$length,ele$length)
#
#	Wilcoxon rank sum test with continuity correction
#
#data:  ino$length and ele$length
#W = 4337, p-value = 0.009222
#alternative hypothesis: true location shift is not equal to 0





summary(ino$width)

#> summary(ino$width)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  19.00   24.00   26.00   26.98   29.00   44.00 
#

sd(ino$width)

#[1] 4.198163

summary(ele$width)

#> summary(ele$width)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   9.00   16.00   17.00   17.88   21.00   28.00 


sd(ele$width)

#[1] 3.534908


wilcox.test(ino$width,ele$width)

#> wilcox.test(ino$width,ele$width)
#
#	Wilcoxon rank sum test with continuity correction
#
#data:  ino$width and ele$width
#W = 10530, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0
#

#Tukey tests

library(multcompView)

#remove inopinata

dat_no_ino <- dat[dat$species != "C. inopinata",]

dat_no_ino$species <- droplevels(dat_no_ino$species)

library(nparcomp)

data.mctp <- mctp(length ~ species, data = dat_no_ino, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)

summary(data.mctp)

# #----------------Nonparametric Multiple Comparisons for relative effects---------------# 
 #
# - Alternative Hypothesis:  True differences of relative effects are less or equal than 0 
# - Estimation Method: Global Pseudo ranks 
# - Type of Contrast : Tukey 
# - Confidence Level: 95 % 
# - Method = Fisher with 182 DF 
 #
# #--------------------------------------------------------------------------------------# 
 #
# #----Data Info-------------------------------------------------------------------------# 
#                     Sample Size    Effect     Lower     Upper
#C. briggsae     C. briggsae   85 0.4097920 0.3760789 0.4443747
#C. tropicalis C. tropicalis   93 0.6638383 0.6309568 0.6952053
#C. elegans       C. elegans  113 0.4263697 0.3887656 0.4648450
#
# #----Contrast--------------------------------------------------------------------------# 
#                            C. briggsae C. tropicalis C. elegans
#C. tropicalis - C. briggsae          -1             1          0
#C. elegans - C. briggsae             -1             0          1
#C. elegans - C. tropicalis            0            -1          1
#
# #----Analysis--------------------------------------------------------------------------# 
#                            Estimator  Lower  Upper Statistic      p.Value
#C. tropicalis - C. briggsae     0.254  0.175  0.330     7.351 1.399281e-11
#C. elegans - C. briggsae        0.017 -0.076  0.109     0.420 9.069100e-01
#C. elegans - C. tropicalis     -0.237 -0.324 -0.147    -6.092 2.478822e-08
#
# #----Overall---------------------------------------------------------------------------# 
#  Quantile      p.Value
#1 2.360246 1.399281e-11
#
# #--------------------------------------------------------------------------------------# 



data.mctp <- mctp(width ~ species, data = dat_no_ino, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)

summary(data.mctp)

# #----------------Nonparametric Multiple Comparisons for relative effects---------------# 
 #
# - Alternative Hypothesis:  True differences of relative effects are less or equal than 0 
# - Estimation Method: Global Pseudo ranks 
# - Type of Contrast : Tukey 
# - Confidence Level: 95 % 
# - Method = Fisher with 174 DF 
 #
# #--------------------------------------------------------------------------------------# 
 #
# #----Data Info-------------------------------------------------------------------------# 
#                     Sample Size    Effect     Lower     Upper
#C. briggsae     C. briggsae   85 0.4493459 0.4099474 0.4893904
#C. tropicalis C. tropicalis   93 0.5621659 0.5235360 0.6000566
#C. elegans       C. elegans  113 0.4884881 0.4501493 0.5269628
#
# #----Contrast--------------------------------------------------------------------------# 
#                            C. briggsae C. tropicalis C. elegans
#C. tropicalis - C. briggsae          -1             1          0
#C. elegans - C. briggsae             -1             0          1
#C. elegans - C. tropicalis            0            -1          1
#
# #----Analysis--------------------------------------------------------------------------# 
#                            Estimator  Lower Upper Statistic   p.Value
#C. tropicalis - C. briggsae     0.113  0.014 0.209     2.705 0.0203603
#C. elegans - C. briggsae        0.039 -0.059 0.137     0.942 0.6143795
#C. elegans - C. tropicalis     -0.074 -0.167 0.021    -1.838 0.1603457
#
# #----Overall---------------------------------------------------------------------------# 
#  Quantile   p.Value
#1 2.363675 0.0203603


data.mctp <- mctp(l_over_w ~ species, data = dat_no_ino, type = "Tukey", conf.level = 0.95, asy.method = "fisher", info = FALSE)

summary(data.mctp)

# #----------------Nonparametric Multiple Comparisons for relative effects---------------# 
 #
# - Alternative Hypothesis:  True differences of relative effects are less or equal than 0 
# - Estimation Method: Global Pseudo ranks 
# - Type of Contrast : Tukey 
# - Confidence Level: 95 % 
# - Method = Fisher with 163 DF 
 #
# #--------------------------------------------------------------------------------------# 
 #
# #----Data Info-------------------------------------------------------------------------# 
#                     Sample Size    Effect     Lower     Upper
#C. briggsae     C. briggsae   85 0.4908359 0.4496513 0.5321452
#C. tropicalis C. tropicalis   93 0.5387729 0.5005745 0.5765214
#C. elegans       C. elegans  113 0.4703912 0.4313360 0.5098123
#
# #----Contrast--------------------------------------------------------------------------# 
#                            C. briggsae C. tropicalis C. elegans
#C. tropicalis - C. briggsae          -1             1          0
#C. elegans - C. briggsae             -1             0          1
#C. elegans - C. tropicalis            0            -1          1
#
# #----Analysis--------------------------------------------------------------------------# 
#                            Estimator  Lower Upper Statistic   p.Value
#C. tropicalis - C. briggsae     0.048 -0.051 0.146     1.140 0.4904356
#C. elegans - C. briggsae       -0.020 -0.122 0.082    -0.473 0.8842002
#C. elegans - C. tropicalis     -0.068 -0.162 0.026    -1.714 0.2026504
#
# #----Overall---------------------------------------------------------------------------# 
#  Quantile   p.Value
#1 2.364683 0.2026504
#
# #--------------------------------------------------------------------------------------# 





####Variation in SDS resistance among C. inopinata lines


#formerly "inopinata_sds.csv" , now "inopinata_sds.tsv"
dat <- read.table("inopinata_sds.tsv", sep="\t", header=T)

dat$total_worms <- dat$alive + dat$dead

dat$fraction_alive <- dat$alive/dat$total_worms

dat$perc_alive <- dat$fraction_alive*100

summary(dat$perc_alive)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.2074  0.1832  2.0000 



sd(dat$perc_alive)

#[1] 0.4327949

> tapply(dat$perc_alive,dat$strain,mean)
#      NKZ2      NKZ22      NKZ27      NKZ43      NKZ44      NKZ45      NKZ46 
#1.04159577 0.00000000 0.51117442 0.39682540 0.00000000 0.00000000 0.08082194 
#     NKZ47      NKZ49      NKZ50      NKZ51      NKZ52      NKZ54      NKZ55 
#0.06451613 0.00000000 0.00000000 0.30674227 0.34965035 0.01594896 0.29411765 
#     NKZ56      NKZ57      NKZ59      NKZ60      NKZ63      NKZ64      NKZ66 
#0.35114046 0.00000000 0.00000000 0.21588915 0.03906250 0.09847395 0.00000000 
#     NKZ67      NKZ68      NKZ69      NKZ70      NKZ72      NKZ73      NKZ75 
#0.06647231 0.00000000 0.06110402 0.00000000 0.00000000 0.35134663 0.23615185 
#     NKZ88 
#1.13211067 

tapply(dat$perc_alive,dat$strain,length)

#> tapply(dat$perc_alive,dat$strain,length)
# NKZ2 NKZ22 NKZ27 NKZ43 NKZ44 NKZ45 NKZ46 NKZ47 NKZ49 NKZ50 NKZ51 NKZ52 NKZ54
#    8     5     2     4     5    10    10     5     5     4     5     2    10
#NKZ55 NKZ56 NKZ57 NKZ59 NKZ60 NKZ63 NKZ64 NKZ66 NKZ67 NKZ68 NKZ69 NKZ70 NKZ72
#    5     5     5     5     5    10    10     5    10     5     5     5     4
#NKZ73 NKZ75 NKZ88
#   10     5    10

summary(dat$total_worms)

#> summary(dat$total_worms)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   17.0   233.0   397.0   449.8   624.5  1557.0

irio <- dat[dat$island == "Iriomote ",]

ishi <- dat[dat$island == "Ishigaki",]

wilcox.test(irio$perc_alive,ishi$perc_alive)

#> wilcox.test(irio$perc_alive,ishi$perc_alive)
#
#	Wilcoxon rank sum test with continuity correction
#
#data:  irio$perc_alive and ishi$perc_alive
#W = 4067.5, p-value = 0.4656
#alternative hypothesis: true location shift is not equal to 0



