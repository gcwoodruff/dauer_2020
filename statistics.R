#this is the revised code for analyzing data associated with manuscript, 
#"Opposing directions of stage-specific body length change in a close relative of C. elegans"
#by Hammerschmith, Woodruff, Johnson and Phillips
#2021

#email Gavin at gcwoodruff@ou.edu if you have any questions

#load libraries

library(dplyr)
library(rstatix)
library(lsmeans)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(cowplot)
library(lemon)
library(ggforce)
library(reshape2)
library(caret)
library(MASS)

#get data in there

dat_o <- read.table("original_length_width_data.tsv", sep="\t", header=T)

#get factor levels right
dat_o$species <- factor(dat_o$species, levels = c("C. elegans","C. inopinata", "C. briggsae", "C. tropicalis"))
dat_o$stage <- factor(dat_o$stage, levels = c("Dauer","L1", "L2", "L3", "L4", "Young adult", "Gravid adult"))
dat_o$species.stage <- factor(dat_o$species.stage, levels = c("C. inopinata Dauer", "C. elegans Dauer", "C. briggsae Dauer", "C. tropicalis Dauer", "C. elegans L1", "C. elegans L2", "C. elegans L3", "C. elegans L4", "C. elegans Young adult", "C. elegans Gravid adult", "C. inopinata L1", "C. inopinata L2", "C. inopinata L3", "C. inopinata L4", "C. inopinata Young adult", "C. inopinata Gravid adult"))


#sample sizes for all species.stage groups

N_df <- count(dat_o,species.stage)

#summary statistics per group: length

sum_df <- aggregate(dat_o$length,list(dat_o$species.stage), summary)
#include sd
sd_df<- aggregate(dat_o$length,list(dat_o$species.stage), sd)
#make new df
stat_df<- data.frame(species.stage=sum_df$Group.1,N=N_df$n,round(sum_df[,-1]),sd=round(sd_df[,-1]))
#this is Supplemental Table Sheet 1
write.table(stat_df,"sheet_1_length_summary_statistics.tsv",sep="\t", quote=FALSE, row.names=FALSE)

stat_df

               #species.stage   N Min. X1st.Qu. Median Mean X3rd.Qu. Max.  sd
#1         C. inopinata Dauer  97  326      399    432  429      462  524  45
#2           C. elegans Dauer 113  333      399    444  456      496  662  65
#3          C. briggsae Dauer  85  319      423    452  451      480  574  45
#4        C. tropicalis Dauer  93  340      464    496  499      534  604  48
#5              C. elegans L1  40  201      229    248  248      264  314  25
#6              C. elegans L2  40  322      358    377  384      401  484  39
#7              C. elegans L3  36  403      502    520  518      538  642  47
#8              C. elegans L4  41  632      714    749  748      774  870  47
#9     C. elegans Young adult  30  686      812    826  832      874 1022  69
#10   C. elegans Gravid adult  40  842      921    962  982     1052 1143  78
#11           C. inopinata L1  20  266      311    323  320      336  345  19
#12           C. inopinata L2  22  328      378    446  438      479  552  67
#13           C. inopinata L3  40  512      582    710  689      782  860 104
#14           C. inopinata L4  33  759      890    960  947      994 1064  72
#15  C. inopinata Young adult  16 1127     1276   1372 1363     1456 1582 127
#16 C. inopinata Gravid adult  32 1319     1510   1596 1608     1702 1821 120

#summary statistics per group: width

sum_df <- aggregate(dat_o$width,list(dat_o$species.stage), summary)
#include sd
sd_df<- aggregate(dat_o$width,list(dat_o$species.stage), sd)
#make new df
stat_df<- data.frame(species.stage=sum_df$Group.1,N=N_df$n,round(sum_df[,-1]),sd=round(sd_df[,-1]))
#this is Supplemental Table Sheet 2
write.table(stat_df,"sheet_2_width_summary_statistics.tsv",sep="\t", quote=FALSE, row.names=FALSE)
stat_df

#               species.stage   N Min. X1st.Qu. Median Mean X3rd.Qu. Max. sd
#1         C. inopinata Dauer  97   19       24     26   27       29   44  4
#2           C. elegans Dauer 113    9       16     17   18       21   28  4
#3          C. briggsae Dauer  85   12       15     17   17       19   26  3
#4        C. tropicalis Dauer  93   12       17     18   19       21   31  3
#5              C. elegans L1  40   11       14     16   16       16   21  2
#6              C. elegans L2  40   18       22     24   23       24   28  2
#7              C. elegans L3  36   28       29     31   31       32   39  3
#8              C. elegans L4  41   37       41     44   44       48   56  5
#9     C. elegans Young adult  30   37       46     50   51       54   65  6
#10   C. elegans Gravid adult  40   48       59     66   65       71   82  8
#11           C. inopinata L1  20   14       16     17   17       18   20  2
#12           C. inopinata L2  22   17       19     22   23       25   35  5
#13           C. inopinata L3  40   22       31     32   34       38   44  5
#14           C. inopinata L4  33   33       42     46   45       48   52  5
#15  C. inopinata Young adult  16   39       50     54   55       60   67  8
#16 C. inopinata Gravid adult  32   59       65     70   69       73   88  6



#percent differences in mean length among groups

#get mean lengths per group
lmean_df <- aggregate(dat_o$length,list(dat_o$species.stage), mean)
#get pairwise percent differences among those means
lmean_perc_diffm <- ((outer(lmean_df$x, lmean_df$x, `-`))/lmean_df$x)*-100
#round
lmean_perc_diffm <- round(lmean_perc_diffm)
#get column and row ids right
colnames(lmean_perc_diffm) <- lmean_df$Group.1
rownames(lmean_perc_diffm) <- lmean_df$Group.1

#this is Supplemental Table Sheet 3
write.table(lmean_perc_diffm,"sheet_3_length_pairwise_percent_differences.tsv",sep="\t", quote=FALSE,row.names=TRUE)



#percent differences in mean width among groups

#get mean widths per group
wmean_df <- aggregate(dat_o$width,list(dat_o$species.stage), mean)
#get pairwise percent differences among those means
wmean_perc_diffm <- ((outer(wmean_df$x, wmean_df$x, `-`))/wmean_df$x)*-100
#round
wmean_perc_diffm <- round(wmean_perc_diffm)
#get column and row ids right
colnames(wmean_perc_diffm) <- wmean_df$Group.1
rownames(wmean_perc_diffm) <- wmean_df$Group.1

#this is Supplemental Table Sheet 4
write.table(wmean_perc_diffm,"sheet_4_width_pairwise_percent_differences.tsv",sep="\t", quote=FALSE,row.names=TRUE)



#Cohen's d effect sizes for length among all group pairwise comparisons
#this is Supplemental Table Sheet 5

write.table(as.data.frame(dat_o %>% cohens_d(length ~ species.stage,ci = TRUE)), "sheet_5_length_pairwise_effect_sizes.tsv",sep="\t", quote=FALSE,row.names=FALSE)


#Cohen's d effect sizes for width among all group pairwise comparisons
#this is Supplemental Table Sheet 6

write.table(as.data.frame(dat_o %>% cohens_d(width ~ species.stage,ci = TRUE)), "sheet_6_width_pairwise_effect_sizes.tsv",sep="\t", quote=FALSE,row.names=FALSE)




#global Kruskal-Wallis rank sum test test for length and width differences and pairwise comparisons among groups


#length
kwt_dato <- kruskal.test(length ~ species.stage, data = dat_o)
kwt_dato

#	Kruskal-Wallis rank sum test
#
#data:  length by species.stage
#Kruskal-Wallis chi-squared = 643.57, df = 15, p-value < 2.2e-16

kwt_dato$p.value
#[1] 1.933859e-127

#width
kwt_dato <- kruskal.test(width ~ species.stage, data = dat_o)
kwt_dato

#	Kruskal-Wallis rank sum test
#
#data:  width by species.stage
#Kruskal-Wallis chi-squared = 666.64, df = 15, p-value < 2.2e-16

kwt_dato$p.value
#[1] 2.380006e-132



#Wilcoxon tests for length among all group pairwise comparisons
#this is Supplemental Table Sheet 7

write.table(as.data.frame(dat_o %>% wilcox_test(length ~ species.stage,p.adjust.method="BH")), "sheet_7_length_all_pairwise_wilcox_p_bh_adjusted.tsv",sep="\t", quote=FALSE,row.names=FALSE)

pair_wilc_len <- dat_o %>% wilcox_test(length ~ species.stage,p.adjust.method="BH")


#how does the C. inopinata dauer length compare to all other species stage groups?
#get C. inopinata dauer comparison index
pair_wilc_len.ci_d_index <- grep("C. inopinata Dauer", pair_wilc_len$group1)
#grep("C. inopinata Dauer", pair_wilc_len$group2)
#integer(0)
#all pairwise comparisons with C. inopinata Dauer have that group in column "group1"

#get those comparisons
pair_wilc_len.ci_d <- pair_wilc_len[pair_wilc_len.ci_d_index,]
#which are not significant
pair_wilc_len.ci_d_ns <- pair_wilc_len.ci_d[which(pair_wilc_len.ci_d$p.adj.signif == "ns"),]

#count fraction not significant
num.sig <- length(pair_wilc_len.ci_d_ns$p.adj.signif)
tot.comp <- length(pair_wilc_len.ci_d$p.adj.signif)
fra.sig <- num.sig/tot.comp

paste(num.sig,tot.comp,fra.sig)
#[1] "1 15 0.0666666666666667"

#which comparison is not significant?
paste(pair_wilc_len.ci_d_ns$group1, pair_wilc_len.ci_d_ns$group2)
#[1] "C. inopinata Dauer C. inopinata L2"

#only one comparison is ns, inopinata dauer and inopinata L2

pair_wilc_wid <- dat_o %>% wilcox_test(width ~ species.stage,p.adjust.method="BH")

#how does the C. inopinata dauer width compare to all other species stage groups?
#get C. inopinata dauer comparison index
pair_wilc_wid.ci_d_index <- grep("C. inopinata Dauer", pair_wilc_wid$group1)
#grep("C. inopinata Dauer", pair_wilc_wid$group2)
#integer(0)
#all pairwise comparisons with C. inopinata Dauer have that group in column "group1"

#get those comparisons
pair_wilc_wid.ci_d <- pair_wilc_wid[pair_wilc_wid.ci_d_index,]
#which are not significant
pair_wilc_wid.ci_d_ns <- pair_wilc_wid.ci_d[which(pair_wilc_wid.ci_d$p.adj.signif == "ns"),]

#count fraction not significant
num.sig <- length(pair_wilc_wid.ci_d_ns$p.adj.signif)
tot.comp <- length(pair_wilc_wid.ci_d$p.adj.signif)
fra.sig <- num.sig/tot.comp

paste(num.sig,tot.comp,fra.sig)
#[1] "0 15 0"

#C. inopinata dauer has significant difference in width among all pairwise comparisons with wilcoxon test


#Wilcoxon tests for width among all group pairwise comparisons
#this is Supplemental Table Sheet 8

write.table(as.data.frame(dat_o %>% wilcox_test(width ~ species.stage,p.adjust.method="BH")), "sheet_8_width_all_pairwise_wilcox_p_bh_adjusted.tsv",sep="\t", quote=FALSE,row.names=FALSE)






##new data shows reproducibility and that the tail only modestly contributes to length difference
#sep elegans and inopinata

#load additional inopinata and elegans dauer data generated for revisions
dat_n <- read.table("new_data_for_revisions_dauer_elegans_inopinata_tails_sheath.tsv", sep="\t", header=T)

#sample sizes

count(dat_n,species)

#       species  n
#1   C. elegans 51
#2 C. inopinata 34


#factor levels with inopinata first for effect size
dat_n$species <- factor(dat_n$species, levels = c("C. inopinata", "C. elegans"))
#set aside species
ci_n<- dat_n[dat_n$species == "C. inopinata",]
ce_n <- dat_n[dat_n$species == "C. elegans",]


# a function to get the percent difference in something
perc_diff = function(w,z){
	((w-z)/z)*100
}

#percent difference in mean dauer length between elegans and inopinata with new data
perc_diff(mean(ci_n$length_include_tail_micron), mean(ce_n$length_include_tail_micron))
#[1] -15.81731

#cohen's d effect size

dat_n %>% cohens_d(length_include_tail_micron ~ species,ci = TRUE)

## A tibble: 1 x 9
#  .y.          group1   group2  effsize    n1    n2 conf.low conf.high magnitude
#* <chr>        <chr>    <chr>     <dbl> <int> <int>    <dbl>     <dbl> <ord>
#1 length_incl… C. inop… C. ele…   -2.14    34    51    -2.74     -1.79 large

#wilcoxon rank sum test

dat_n %>% wilcox_test(length_include_tail_micron ~ species)
#
## A tibble: 1 x 7
#  .y.                      group1       group2       n1    n2 statistic        p
#* <chr>                    <chr>        <chr>     <int> <int>     <dbl>    <dbl>
#1 length_include_tail_mic… C. inopinata C. elega…    34    51        55 3.35e-13

#tail length

dat_n$tail_length_micron <- round(dat_n$tail_length_micron)

ci_n<- dat_n[dat_n$species == "C. inopinata",]
ce_n <- dat_n[dat_n$species == "C. elegans",]


perc_diff(mean(na.omit(ci_n$tail_length_micron)), mean(na.omit(ce_n$tail_length_micron)))
#[1] -46.66335

#cohen's d effect size

dat_n %>% cohens_d(tail_length_micron ~ species,ci = TRUE)

#  .y.         group1    group2  effsize    n1    n2 conf.low conf.high magnitude
#* <chr>       <chr>     <chr>     <dbl> <int> <int>    <dbl>     <dbl> <ord>
#1 tail_lengt… C. inopi… C. ele…   -3.10    34    51    -3.81     -2.62 large



dat_n %>% cohens_d(tail_length_micron ~ species,ci = TRUE)

#wilcoxon rank sum test


dat_n %>% wilcox_test(tail_length_micron ~ species)

## A tibble: 1 x 7
#  .y.                group1       group2        n1    n2 statistic        p
#* <chr>              <chr>        <chr>      <int> <int>     <dbl>    <dbl>
#1 tail_length_micron C. inopinata C. elegans    34    51        14 3.42e-14

#length without tail

dat_n$length_no_tail <- round(dat_n$length_include_tail_micron - dat_n$tail_length_micron)

ci_n<- dat_n[dat_n$species == "C. inopinata",]
ce_n <- dat_n[dat_n$species == "C. elegans",]


#percent difference among mean length without tail
perc_diff(mean(na.omit(ci_n$length_no_tail)), mean(na.omit(ce_n$length_no_tail)))
#[1] -13.84804
#cohen's d effect size


dat_n %>% cohens_d(length_no_tail ~ species,ci = TRUE)

#  .y.        group1    group2   effsize    n1    n2 conf.low conf.high magnitude
#* <chr>      <chr>     <chr>      <dbl> <int> <int>    <dbl>     <dbl> <ord>
#1 length_no… C. inopi… C. eleg…   -1.78    34    51    -2.33      -1.4 large


#wilcoxon rank sum test


dat_n %>% wilcox_test(length_no_tail ~ species)

## A tibble: 1 x 7
#  .y.            group1       group2        n1    n2 statistic        p
#* <chr>          <chr>        <chr>      <int> <int>     <dbl>    <dbl>
#1 length_no_tail C. inopinata C. elegans    34    51      132. 8.09e-11



#length-width relationship across non-dauer stages


#elegans non dauer stages
ce_repr <- dat_o[dat_o$species == "C. elegans" & dat_o$stage != "Dauer",]
#elegans linear fit
ce_repr_fit <- lm(ce_repr$length~ce_repr$width)
#elegans dauer stage
ce_dauer <- dat_o[dat_o$species == "C. elegans" & dat_o$stage == "Dauer",]
#do dauers fall above or below fit?
ce_dauer$is_above_regression <- ifelse(ce_dauer$length > (ce_repr_fit$coefficients[2]*ce_dauer$width) + ce_repr_fit$coefficients[1],'Dauer above','Dauer below')
#summary
summary(ce_repr_fit)

#Call:
#lm(formula = ce_repr$length ~ ce_repr$width)
#
#Residuals:
#     Min       1Q   Median       3Q      Max
#-283.572  -40.897   -0.135   36.963  194.144
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)    62.2214    10.6027   5.868 1.56e-08 ***
#ce_repr$width  14.4799     0.2532  57.194  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 67.65 on 225 degrees of freedom
#Multiple R-squared:  0.9356,	Adjusted R-squared:  0.9354
#F-statistic:  3271 on 1 and 225 DF,  p-value: < 2.2e-16
summary(ce_repr_fit)$coefficients
#              Estimate Std. Error   t value      Pr(>|t|)
#(Intercept)   62.22140 10.6027330  5.868431  1.562474e-08
#ce_repr$width 14.47988  0.2531729 57.193649 5.091459e-136


#inopinata non dauer stages
ci_repr <- dat_o[dat_o$species == "C. inopinata" & dat_o$stage != "Dauer",]
#inopinata linear fit
ci_repr_fit <- lm(ci_repr$length~ci_repr$width)
#inopinata dauer stage
ci_dauer <- dat_o[dat_o$species == "C. inopinata" & dat_o$stage == "Dauer",]
#do dauers fall above or below fit?
ci_dauer$is_above_regression <- ifelse(ci_dauer$length > (ci_repr_fit$coefficients[2]*ci_dauer$width) + ci_repr_fit$coefficients[1],'Dauer above','Dauer below')
#summary
summary(ci_repr_fit)
#Call:
#lm(formula = ci_repr$length ~ ci_repr$width)
#
#Residuals:
#    Min      1Q  Median      3Q     Max
#-577.42  -65.44   -8.30   44.63  433.87
#
#Coefficients:
#              Estimate Std. Error t value Pr(>|t|)
#(Intercept)   -80.9462    26.1831  -3.092  0.00235 **
#ci_repr$width  23.8224     0.5764  41.327  < 2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 135.1 on 161 degrees of freedom
#Multiple R-squared:  0.9139,	Adjusted R-squared:  0.9133
#F-statistic:  1708 on 1 and 161 DF,  p-value: < 2.2e-16
#coefficients
summary(ci_repr_fit)$coefficients
#               Estimate Std. Error  t value     Pr(>|t|)
#(Intercept)   -80.94624 26.1831486 -3.09154 2.347217e-03
#ci_repr$width  23.82239  0.5764307 41.32742 1.267441e-87



#chi-square classification of dauers above/below fit

dauers <- (rbind(ci_dauer,ce_dauer))
dauers$species<-droplevels(dauers$species)
dFrequencies <- xtabs( ~ is_above_regression + species, data = dauers)
xtabs( formula = ~ is_above_regression + species, data = dauers)

#is_above_regression C. elegans C. inopinata
#        Dauer above        111            7
#        Dauer below          2           90

chisq.test(dFrequencies)

#	Pearson's Chi-squared test with Yates' continuity correction
#
#data:  dFrequencies
#X-squared = 171.96, df = 1, p-value < 2.2e-16



#is there an interaction between species and width on length?

repr_dat <- rbind(ce_repr,ci_repr)

m.interaction <- lm(length~width*species,data=repr_dat)

anova(m.interaction)
#Analysis of Variance Table
#
#Response: length
#               Df   Sum Sq  Mean Sq F value    Pr(>F)
#width           1 46847194 46847194 4555.08 < 2.2e-16 ***
#species         1  4984021  4984021  484.61 < 2.2e-16 ***
#width:species   1  2710655  2710655  263.56 < 2.2e-16 ***
#Residuals     386  3969860    10285
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#coefficients
summary(m.interaction)$coefficients
#                            Estimate Std. Error   t value      Pr(>|t|)
#(Intercept)                 62.22140 15.8936529  3.914858  1.069269e-04
#width                       14.47988  0.3795099 38.154161 4.802043e-133
#speciesC. inopinata       -143.16764 25.2727388 -5.664904  2.879424e-08
#width:speciesC. inopinata    9.34251  0.5754669 16.234660  1.506259e-45

#get coffecients for each species
m.lst <- lstrends(m.interaction, "species", var="width")
m.lst
# species      width.trend    SE  df lower.CL upper.CL
# C. elegans          14.5 0.380 386     13.7     15.2
# C. inopinata        23.8 0.433 386     23.0     24.7
#Confidence level used: 0.95

# Compare slopes
pairs(m.lst)
# contrast                  estimate    SE  df t.ratio p.value
# C. elegans - C. inopinata    -9.34 0.575 386 -16.235 <.0001




#permanova, size among species.stage groups

#make size matrix
size_matrix <- dat_o[, c("length", "width")]

## Perform the permanova:
pmnva <- adonis(size_matrix ~ dat_o$species + dat_o$stage + dat_o$species:dat_o$stage, method = "euclidean")
pmnva

#Call:
#adonis(formula = size_matrix ~ dat_o$species + dat_o$stage +      dat_o$species:dat_o$stage, method = "euclidean")
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#                           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
#dat_o$species               3   7819986 2606662  647.67 0.10430  0.001 ***
#dat_o$stage                 6  56672987 9445498 2346.91 0.75589  0.001 ***
#dat_o$species:dat_o$stage   6   7415890 1235982  307.10 0.09891  0.001 ***
#Residuals                 762   3066788    4025         0.04090
#Total                     777  74975651                 1.00000
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#do post-hoc pairwise tests to compare all groups
post_hoc_pairwise_permanova <- pairwise.adonis(size_matrix,dat_o$species.stage,p.adjust.m="BH")

#this is Supplemental Table Sheet 9 
write.table(post_hoc_pairwise_permanova,file="sheet_9_post-hoc_pairwise_permanova_original_data.tsv",quote=FALSE,row.names = FALSE,sep="\t")




#how does the C. inopinata dauer length compare to all other species stage groups?

#get C. inopinata dauer index
post_hoc_pairwise_permanova.ci_d_index <- grep("C. inopinata Dauer", post_hoc_pairwise_permanova$pairs)

#get the C. inopinata dauer comparisons
post_hoc_pairwise_permanova.ci_d <- post_hoc_pairwise_permanova[post_hoc_pairwise_permanova.ci_d_index,]

#get those that are not significant
ci.dau.pairs.pmva.ns <- post_hoc_pairwise_permanova.ci_d[which(post_hoc_pairwise_permanova.ci_d$sig == ""),]

ns.ci.dau.pmnva.comp <- as.character(droplevels(ci.dau.pairs.pmva.ns$pairs))

#get the fraction not significant
num.sig <- length(ci.dau.pairs.pmva.ns$sig)
tot.comp <- length(post_hoc_pairwise_permanova[post_hoc_pairwise_permanova.ci_d_index,]$sig)
fra.sig <- num.sig/tot.comp

paste(num.sig,tot.comp,fra.sig)
#[1] "1 15 0.0666666666666667"
ns.ci.dau.pmnva.comp 
#[1] "C. inopinata Dauer vs C. inopinata L2"

#in this framework,  C. inopinata dauer is significantly different from all comparisons aside from C. inopinata L2

#k means clustering (for subsampling see subsample.R)


dat_lw <- data.frame(length=dat_o$length,width=dat_o$width)

dat_lws <- scale(dat_lw)



#trying out various values of k

k1 <- kmeans(dat_lws, centers = 1, iter.max = 20, nstart = 25)
k2 <- kmeans(dat_lws, centers = 2,iter.max = 20, nstart = 25)
k3 <- kmeans(dat_lws, centers = 3,iter.max = 20, nstart = 25)
k4 <- kmeans(dat_lws, centers = 4,iter.max = 20, nstart = 25)
k5 <- kmeans(dat_lws, centers = 5,iter.max = 20, nstart = 25)
k6 <- kmeans(dat_lws, centers = 6,iter.max = 20, nstart = 25)
k7 <- kmeans(dat_lws, centers = 7,iter.max = 20, nstart = 25)
k8 <- kmeans(dat_lws, centers = 8,iter.max = 20, nstart = 25)
k9 <- kmeans(dat_lws, centers = 9,iter.max = 20, nstart = 25)
k10 <- kmeans(dat_lws, centers = 10,iter.max = 20, nstart = 25)
k11 <- kmeans(dat_lws, centers = 11,iter.max = 20, nstart = 25)
k12 <- kmeans(dat_lws, centers = 12,iter.max = 20, nstart = 25)
k13 <- kmeans(dat_lws, centers = 13,iter.max = 20, nstart = 25)
k14 <- kmeans(dat_lws, centers = 14,iter.max = 20, nstart = 25)
k15 <- kmeans(dat_lws, centers = 15,iter.max = 20, nstart = 25)
k16 <- kmeans(dat_lws, centers = 16,iter.max = 20, nstart = 25)
k17 <- kmeans(dat_lws, centers = 17,iter.max = 20, nstart = 25)
k18 <- kmeans(dat_lws, centers = 18,iter.max = 20, nstart = 25)
k19 <- kmeans(dat_lws, centers = 19,iter.max = 20, nstart = 25)
k20 <- kmeans(dat_lws, centers = 20,iter.max = 20, nstart = 25)

#add classifications to data
dat_o <- cbind(dat_o, k1 = k1$cluster,k2 = k2$cluster,k3 = k3$cluster,k4 = k4$cluster,k5 = k5$cluster,k6 = k6$cluster,k7 = k7$cluster,k8 = k8$cluster,k9 = k9$cluster,k10 = k10$cluster,k11 = k11$cluster,k12 = k12$cluster,k13 = k13$cluster,k14 = k14$cluster,k15 = k15$cluster,k16 = k16$cluster,k17 = k17$cluster,k18 = k18$cluster,k19 = k19$cluster,k20 = k20$cluster)

#what is the best value of k?
#try elbow plot approach

#get within-cluster sum of squares for all values of k 
clust_var <-data.frame(k=c(1:20), wss= c(k1$tot.withinss,k2$tot.withinss,k3$tot.withinss,k4$tot.withinss,k5$tot.withinss,k6$tot.withinss,k7$tot.withinss,k8$tot.withinss,k9$tot.withinss,k10$tot.withinss,k11$tot.withinss,k12$tot.withinss,k13$tot.withinss,k14$tot.withinss,k15$tot.withinss,k16$tot.withinss,k17$tot.withinss,k18$tot.withinss,k19$tot.withinss,k20$tot.withinss))

ggplot(clust_var,aes(x=k,y=wss)) + geom_line() + theme_cowplot() + scale_x_continuous(breaks=c(1:20)) + ylab("Total within-cluster sum of squares")


#best k is 3-5 by elbow plot, what does BIC say?
# http://sherrytowers.com/2013/10/24/k-means-clustering/ , https://stackoverflow.com/a/25557162

#define function to get AIC and BIC
kmeansAIC = function(fit){

m = ncol(fit$centers)
n = length(fit$cluster)
k = nrow(fit$centers)
D = fit$tot.withinss
return(data.frame(AIC = D + 2*m*k, BIC = D + log(n)*m*k))
}

#add BIC to data
clust_var <-cbind(clust_var, BIC = c(kmeansAIC(k1)$BIC,kmeansAIC(k2)$BIC,kmeansAIC(k3)$BIC,kmeansAIC(k4)$BIC,kmeansAIC(k5)$BIC,kmeansAIC(k6)$BIC,kmeansAIC(k7)$BIC,kmeansAIC(k8)$BIC,kmeansAIC(k9)$BIC,kmeansAIC(k10)$BIC,kmeansAIC(k11)$BIC,kmeansAIC(k12)$BIC,kmeansAIC(k13)$BIC,kmeansAIC(k14)$BIC,kmeansAIC(k15)$BIC,kmeansAIC(k16)$BIC,kmeansAIC(k17)$BIC,kmeansAIC(k18)$BIC,kmeansAIC(k19)$BIC,kmeansAIC(k20)$BIC), AIC = c(kmeansAIC(k1)$AIC,kmeansAIC(k2)$AIC,kmeansAIC(k3)$AIC,kmeansAIC(k4)$AIC,kmeansAIC(k5)$AIC,kmeansAIC(k6)$AIC,kmeansAIC(k7)$AIC,kmeansAIC(k8)$AIC,kmeansAIC(k9)$AIC,kmeansAIC(k10)$AIC,kmeansAIC(k11)$AIC,kmeansAIC(k12)$AIC,kmeansAIC(k13)$AIC,kmeansAIC(k14)$AIC,kmeansAIC(k15)$AIC,kmeansAIC(k16)$AIC,kmeansAIC(k17)$AIC,kmeansAIC(k18)$AIC,kmeansAIC(k19)$AIC,kmeansAIC(k20)$AIC))


ggplot(clust_var,aes(x=k,y=BIC)) + geom_point() + theme_cowplot() + scale_x_continuous(breaks=c(1:20)) 

#get best k that minimizes BIC
clust_var[which.min(clust_var$BIC),]$k
#[1] 7
min(clust_var$BIC)
#[1] 179.832

#min BIC is that of k = 7

#so elbow test best k = 3-5, BIC = 7, and there are seven stages, four species, and 16 species.stage groups. Let's look at k=7.


#set aside species
ci <- dat_o[dat_o$species=="C. inopinata",]
ce <- dat_o[dat_o$species=="C. elegans",]
cb <- dat_o[dat_o$species=="C. briggsae",]
ct <- dat_o[dat_o$species=="C. tropicalis",]


#visualize k = 7

dat_o$k7 <- as.factor(dat_o$k7)
dat_o$Stage <- dat_o$stage

ggplot(dat_o, aes(x=width, y=length, colour=k7)) + geom_point(aes(shape=Stage),alpha=0.9) + facet_rep_wrap(~species,nrow=2) + scale_shape_manual(values=1:7) + theme_cowplot() + scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","black")) + theme(strip.text.x = element_text(face = "italic"), strip.background = element_blank()) + xlab("Width (microns)") + ylab("Length (microns)") + xlim(0,100) + scale_y_continuous(limits=c(0,2000),breaks=c(0,250,500,750,1000,1250,1500,1750,2000)) + labs(colour="Cluster")



#this is supplemental figure 6
#supplemental figure 7 below, other supplemental figures in figures.R

#get chi square stats for all pairwise comparisons


#get all species.stage pairs
species.stage.pairs <- t(combn(unique(dat_o$species.stage),2))

species.stage.pairs.df <- data.frame(grp1=species.stage.pairs[,1], grp2=species.stage.pairs[,2])

species.stage.pairs.df$grp1 <- as.character(species.stage.pairs.df$grp1)

species.stage.pairs.df$grp2 <- as.character(species.stage.pairs.df$grp2)

#fix data
dat_o_km <- dat_o

dat_o_km$k7 <- as.factor(dat_o_km$k7)
dat_o_km$species.stage <- as.factor(dat_o_km$species.stage)

p.df <- NULL


#get pairwise comparisons, chisquare and fishers exact test for classification with k means clustering

CS_func_ii = function(x){

	grp1 <- x[1]

	grp2 <- x[2]

	k7_df <- dat_o_km[dat_o_km$species.stage == grp1 | dat_o_km$species.stage == grp2,]
	
	k7_df$species.stage <- droplevels(k7_df$species.stage)
	k7_df$k7 <- droplevels(k7_df$k7)
	
	k7Frequencies <- xtabs( ~ k7 + species.stage, data = k7_df)

	if(nrow(k7Frequencies)<=1){

	p.df <- rbind(p.df,data.frame(grp1 = grp1, grp2 = grp2, chisq_stat=NA,chisq_p=NA,fisher_p=NA))

	print(p.df)

	} else { 
	
	p.df <- rbind(p.df,data.frame(grp1 = grp1, grp2 = grp2, chisq_stat=as.numeric(chisq.test(k7Frequencies)$statistic),chisq_p=as.numeric(chisq.test(k7Frequencies)$p.value),fisher_p=fisher.test(k7Frequencies)$p.value))

	print(p.df)}
}


CS_out <- apply(species.stage.pairs.df, 1, CS_func_ii)


CS_out_df <- bind_rows(CS_out)

#adjust p for multiple comparisons
CS_out_df$chisq_p_adjust <- p.adjust(CS_out_df$chisq_p,method='BH')
CS_out_df$fisher_p_adjust <- p.adjust(CS_out_df$fisher_p,method='BH')

#how many pairwise comparisons are significant among C. inopinata dauers??
pair_chsq_fet_km.ci_d_index <- grep("C. inopinata Dauer", CS_out_df$grp1)

#grep("C. inopinata Dauer", CS_out_df$grp2)
#integer(0)
#C. inopinata Dauer is always in first group

pair_chsq_fet_km.ci_d <- CS_out_df[pair_chsq_fet_km.ci_d_index,]

pair_chsq_fet_km.ns <- pair_chsq_fet_km.ci_d[which(pair_chsq_fet_km.ci_d$chisq_p_adjust >= 0.05),]

num.sig <- length(pair_chsq_fet_km.ns$chisq_p_adjust)
tot.comp <- length(pair_chsq_fet_km.ci_d$chisq_p_adjust)
fra.sig <- num.sig/tot.comp

paste(num.sig,tot.comp,fra.sig)
#[1] "0 15 0"

#in this framework, C. inopinata dauers are uniquely classified compared to all other species-stage groups

#this is Supplemental Table Sheet 10
write.table(CS_out_df,file="sheet_10_k_means_chisq_fet_pairwise_comparisons.tsv",quote=FALSE,row.names = FALSE,sep="\t")

#make supplemental figure 7
#set aside species
ci <- dat_o[dat_o$species=="C. inopinata",]
ce <- dat_o[dat_o$species=="C. elegans",]
cb <- dat_o[dat_o$species=="C. briggsae",]
ct <- dat_o[dat_o$species=="C. tropicalis",]



#get counts of all classifications for all species and stages
ci_stage_count <- count(ci,stage)
ci_L1_k_count <- count(ci[ci$stage== "L1",],k7)
ci_L2_k_count <- count(ci[ci$stage== "L2",],k7)
ci_L3_k_count <- count(ci[ci$stage== "L3",],k7)
ci_L4_k_count <- count(ci[ci$stage== "L4",],k7)
ci_young_adult_k_count <- count(ci[ci$stage== "Young adult",],k7)
ci_gravid_adult_k_count <- count(ci[ci$stage== "Gravid adult",],k7)
ci_dauer_k_count <- count(ci[ci$stage== "Dauer",],k7)

ci_L1_k_count$fra_k <- ci_L1_k_count$n/ci_stage_count[ci_stage_count$stage=="L1",]$n
ci_L2_k_count$fra_k <- ci_L2_k_count$n/ci_stage_count[ci_stage_count$stage=="L2",]$n
ci_L3_k_count$fra_k <- ci_L3_k_count$n/ci_stage_count[ci_stage_count$stage=="L3",]$n
ci_L4_k_count$fra_k <- ci_L4_k_count$n/ci_stage_count[ci_stage_count$stage=="L4",]$n
ci_young_adult_k_count$fra_k <- ci_young_adult_k_count$n/ci_stage_count[ci_stage_count$stage=="Young adult",]$n
ci_gravid_adult_k_count$fra_k <- ci_gravid_adult_k_count$n/ci_stage_count[ci_stage_count$stage=="Gravid adult",]$n
ci_dauer_k_count$fra_k <- ci_dauer_k_count$n/ci_stage_count[ci_stage_count$stage=="Dauer",]$n

ce_stage_count <- count(ce,stage)
ce_L1_k_count <- count(ce[ce$stage== "L1",],k7)
ce_L2_k_count <- count(ce[ce$stage== "L2",],k7)
ce_L3_k_count <- count(ce[ce$stage== "L3",],k7)
ce_L4_k_count <- count(ce[ce$stage== "L4",],k7)
ce_young_adult_k_count <- count(ce[ce$stage== "Young adult",],k7)
ce_gravid_adult_k_count <- count(ce[ce$stage== "Gravid adult",],k7)
ce_dauer_k_count <- count(ce[ce$stage== "Dauer",],k7)

ce_L1_k_count$fra_k <- ce_L1_k_count$n/ce_stage_count[ce_stage_count$stage=="L1",]$n
ce_L2_k_count$fra_k <- ce_L2_k_count$n/ce_stage_count[ce_stage_count$stage=="L2",]$n
ce_L3_k_count$fra_k <- ce_L3_k_count$n/ce_stage_count[ce_stage_count$stage=="L3",]$n
ce_L4_k_count$fra_k <- ce_L4_k_count$n/ce_stage_count[ce_stage_count$stage=="L4",]$n
ce_young_adult_k_count$fra_k <- ce_young_adult_k_count$n/ce_stage_count[ce_stage_count$stage=="Young adult",]$n
ce_gravid_adult_k_count$fra_k <- ce_gravid_adult_k_count$n/ce_stage_count[ce_stage_count$stage=="Gravid adult",]$n
ce_dauer_k_count$fra_k <- ce_dauer_k_count$n/ce_stage_count[ce_stage_count$stage=="Dauer",]$n

ci_L1_k_count$stage <- "L1"
ci_L2_k_count$stage <- "L2"
ci_L3_k_count$stage <- "L3"
ci_L4_k_count$stage <- "L4"
ci_young_adult_k_count$stage <- "Young adult"
ci_gravid_adult_k_count$stage <- "Gravid adult"
ci_dauer_k_count$stage <- "Dauer"

ce_L1_k_count$stage <- "L1"
ce_L2_k_count$stage <- "L2"
ce_L3_k_count$stage <- "L3"
ce_L4_k_count$stage <- "L4"
ce_young_adult_k_count$stage <- "Young adult"
ce_gravid_adult_k_count$stage <- "Gravid adult"
ce_dauer_k_count$stage <- "Dauer"

ci2<- rbind(ci_L1_k_count,ci_L2_k_count,ci_L3_k_count,ci_L4_k_count,ci_young_adult_k_count,ci_gravid_adult_k_count,ci_dauer_k_count)


ce2<- rbind(ce_L1_k_count,ce_L2_k_count,ce_L3_k_count,ce_L4_k_count,ce_young_adult_k_count,ce_gravid_adult_k_count,ce_dauer_k_count)

ci2$species <- "C. inopinata"
ce2$species <- "C. elegans"


cb_stage_count <- count(cb,stage)
cb_dauer_k_count <- count(cb[cb$stage== "Dauer",],k7)
cb_dauer_k_count$fra_k <- cb_dauer_k_count$n/cb_stage_count[cb_stage_count$stage=="Dauer",]$n
cb_dauer_k_count$stage <- "Dauer"
cb_dauer_k_count$species <- "C. briggsae"

ct_stage_count <- count(ct,stage)
ct_dauer_k_count <- count(ct[ct$stage== "Dauer",],k7)
ct_dauer_k_count$fra_k <- ct_dauer_k_count$n/ct_stage_count[ct_stage_count$stage=="Dauer",]$n
ct_dauer_k_count$stage <- "Dauer"
ct_dauer_k_count$species <- "C. tropicalis"

dat_k7 <- rbind(ci2,ce2,cb_dauer_k_count,ct_dauer_k_count)

dat_k7$k7 <- as.factor(dat_k7$k7)


dat_k7$stage <- factor(dat_k7$stage, levels = c("Dauer","L1", "L2", "L3", "L4", "Young adult", "Gravid adult"))


dat_k7$species <- factor(dat_k7$species, levels = c("C. elegans","C. inopinata", "C. briggsae", "C. tropicalis"))



ggplot(dat_k7, aes(x = stage, y = n, fill = k7)) + geom_col() + facet_rep_wrap(~species,nrow=2) + scale_fill_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","black")) + theme_cowplot() + theme(axis.text.x = element_text(angle = 45,hjust=1), strip.text.x = element_text(face = "italic"),, strip.background = element_blank()) + ylim(0,170) + xlab("Stage") + ylab("Number of animals") + labs(fill="Cluster")
	#this is supplemental figure 7



#for subsampling see subsample.R





#linear discriminant analysis

	#http://www.sthda.com/english/articles/36-classification-methods-essentials/146-discriminant-analysis-essentials-in-r/


library(tidyverse)
library(caret)
library(MASS)

#reload data to get rid of k means crap
dat_o <- read.table("original_length_width_data.tsv", sep="\t", header=T)
dat_o$species <- factor(dat_o$species, levels = c("C. elegans","C. inopinata", "C. briggsae", "C. tropicalis"))
dat_o$stage <- factor(dat_o$stage, levels = c("Dauer","L1", "L2", "L3", "L4", "Young adult", "Gravid adult"))


o_d <- dat_o[dat_o$stage == "Dauer",]

#training and test data for predictions
set.seed(123)
training.samples <- o_d$species %>% createDataPartition(p = 0.8, list = FALSE)
train.data <- o_d[training.samples, ]
test.data <- o_d[-training.samples, ]

# Estimate preprocessing parameters
preproc.param <- train.data %>% preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)

# Fit the model
model <- lda(species~length+width, data = train.transformed)

# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy (for all species)
mean(predictions$class==test.transformed$species)
#[1] 0.5657895

#is the prediction correct?
predictions$correct.prediction <- predictions$class==test.transformed$species
pred_df <- data.frame(class=predictions$class,correct.prediction=predictions$correct.prediction)


#model accuracy-- what fraction of points predicted as C. inopinata are C. inopinata??
count(pred_df[pred_df$class=="C. inopinata",],correct.prediction)

#  correct.prediction  n
#1              FALSE  1
#2               TRUE 18

#0.9473684 accuracy for C. inopinata

#now, do all of the data

# Estimate preprocessing parameters
preproc.param <- o_d %>% preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
o_d.transformed <- preproc.param %>% predict(o_d)

# Fit the model
model <- lda(species~length+width, data = o_d.transformed)
model

#Call:
#lda(species ~ length + width, data = o_d.transformed)
#
#Prior probabilities of groups:
#   C. elegans  C. inopinata   C. briggsae C. tropicalis
#    0.2912371     0.2500000     0.2190722     0.2396907
#
#Group means:
#                   length      width
#C. elegans    -0.04455042 -0.4451998
#C. inopinata  -0.50102791  1.2648863
#C. briggsae   -0.13645842 -0.5408137
#C. tropicalis  0.70142871 -0.2840562
#
#Coefficients of linear discriminants:
#              LD1        LD2
#length  0.7118029 -0.9520852
#width  -1.5317487 -0.4107314
#
#Proportion of trace:
#   LD1    LD2
#0.9387 0.0613

#visualize
lda.data <- cbind(o_d.transformed, predict(model)$x)

#this is supplemental figure 9
ggplot(lda.data, aes(LD1, LD2)) + geom_point(aes(color = species)) + scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442")) + theme_cowplot() + theme(legend.text = element_text(face = "italic"))




#dauer formation frequency

sds_dat <- read.table("dauer_formation_frequency.tsv", sep="\t", header=T)



#get species factor levels right
levels(sds_dat$Species)[levels(sds_dat$Species)=="Elegans"] <- "C. elegans"
levels(sds_dat$Species)[levels(sds_dat$Species)=="sp_34"] <- "C. inopinata"
levels(sds_dat$Species)[levels(sds_dat$Species)=="Briggsae"] <- "C. briggsae"
levels(sds_dat$Species)[levels(sds_dat$Species)=="Tropicalis "] <- "C. tropicalis"

sds_dat$Species <- factor(sds_dat$Species, levels = c("C. elegans","C. inopinata", "C. briggsae", "C. tropicalis"))


#get estimate of dauer survival

#total worms

sds_dat$Control.Total <- sds_dat$Control.Alive + sds_dat$Control.Dead

sds_dat$SDS.Total <- sds_dat$SDS.Alive + sds_dat$SDS.Dead

#fraction alive
sds_dat$Control.fraction_alive <- sds_dat$Control.Alive/sds_dat$Control.Total

sds_dat$SDS.fraction_alive <- sds_dat$SDS.Alive/sds_dat$SDS.Total

#fraction dauer


sds_dat$fraction_dauer <- sds_dat$Control.fraction_alive * sds_dat$SDS.fraction_alive 

#inopinata variation in dauer formation frequency


dat <- sds_dat[sds_dat$Species == "C. inopinata",]

dat$Species <- droplevels(dat$Species)
dat$Strain <- droplevels(dat$Strain)

summary(dat$fraction_dauer)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000000 0.000000 0.000000 0.001776 0.001673 0.017215

sd(dat$fraction_dauer)

#[1] 0.003691608

> tapply(dat$fraction_dauer,dat$Strain,mean)
#        NKZ2        NKZ22        NKZ27        NKZ43        NKZ44        NKZ45
#0.0083323946 0.0000000000 0.0045784186 0.0038285267 0.0000000000 0.0000000000
#       NKZ46        NKZ47        NKZ49        NKZ50        NKZ51        NKZ52
#0.0007565893 0.0005972533 0.0000000000 0.0000000000 0.0030524620 0.0031347962
#       NKZ54        NKZ55        NKZ56        NKZ57        NKZ59        NKZ60
#0.0001415653 0.0029107505 0.0030186725 0.0000000000 0.0000000000 0.0018537911
#       NKZ63        NKZ64        NKZ66        NKZ67        NKZ68        NKZ69
#0.0003557922 0.0009128790 0.0000000000 0.0006093522 0.0000000000 0.0005918483
#       NKZ70        NKZ72        NKZ73        NKZ75        NKZ88
#0.0000000000 0.0000000000 0.0029297565 0.0017353076 0.0094646542

#number of observations per strain
tapply(dat$fraction_dauer,dat$Strain,length)

#> tapply(dat$perc_alive,dat$strain,length)
# NKZ2 NKZ22 NKZ27 NKZ43 NKZ44 NKZ45 NKZ46 NKZ47 NKZ49 NKZ50 NKZ51 NKZ52 NKZ54
#    8     5     2     4     5    10    10     5     5     4     5     2    10
#NKZ55 NKZ56 NKZ57 NKZ59 NKZ60 NKZ63 NKZ64 NKZ66 NKZ67 NKZ68 NKZ69 NKZ70 NKZ72
#    5     5     5     5     5    10    10     5    10     5     5     5     4
#NKZ73 NKZ75 NKZ88
#   10     5    10

summary(dat$SDS.Total)

#> summary(dat$total_worms)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   17.0   233.0   397.0   449.8   624.5  1557.0

#set aside islands
irio <- dat[dat$Island == "Iriomote",]

ishi <- dat[dat$Island == "Ishigaki",]

wilcox.test(irio$fraction_dauer,ishi$fraction_dauer)

#	Wilcoxon rank sum test with continuity correction
#
#data:  irio$fraction_dauer and ishi$fraction_dauer
#W = 4062.5, p-value = 0.4769
#alternative hypothesis: true location shift is not equal to 0







#new data, dauer formation frequency in C. inopinata at 30 deg.


sds_new_dat <- read.table("dauer_formation_frequency_new_data.tsv", sep="\t", header=T)

#get estimate of dauer survival

#total worms

sds_new_dat$Control.Total <- sds_new_dat$Control.Alive + sds_new_dat$Control.Dead

sds_new_dat$SDS.Total <- sds_new_dat$SDS.Alive + sds_new_dat$SDS.Dead

#fraction alive
sds_new_dat$Control.fraction_alive <- sds_new_dat$Control.Alive/sds_new_dat$Control.Total

sds_new_dat$SDS.fraction_alive <- sds_new_dat$SDS.Alive/sds_new_dat$SDS.Total

#fraction dauer

sds_new_dat$fraction_dauer <- sds_new_dat$Control.fraction_alive * sds_new_dat$SDS.fraction_alive 

#groups

sds_new_dat$temp.strain <- paste(sds_new_dat$temperature,sds_new_dat$strain)


kruskal.test(fraction_dauer ~ temp.strain, data = sds_new_dat)

#are there significant differences among groups?

#	Kruskal-Wallis rank sum test
#
#data:  fraction_dauer by temp.strain
#Kruskal-Wallis chi-squared = 7.6435, df = 3, p-value = 0.05398

#pairwise wilcox

sds_new_dat %>% wilcox_test(fraction_dauer ~ temp.strain,p.adjust.method="BH")

## A tibble: 6 x 9
#  .y.            group1   group2     n1    n2 statistic     p p.adj p.adj.signif
#* <chr>          <chr>    <chr>   <int> <int>     <dbl> <dbl> <dbl> <chr>
#1 fraction_dauer 25 NKZ2  25 PX7…     4     3         5 0.854 0.854 ns
#2 fraction_dauer 25 NKZ2  30 NKZ2     4     4        12 0.186 0.275 ns
#3 fraction_dauer 25 NKZ2  30 PX7…     4     4         3 0.191 0.275 ns
#4 fraction_dauer 25 PX723 30 NKZ2     3     4        10 0.123 0.275 ns
#5 fraction_dauer 25 PX723 30 PX7…     3     4         2 0.229 0.275 ns
#6 fraction_dauer 30 NKZ2  30 PX7…     4     4         0 0.021 0.127 ns

#and, does temperature matter?

kruskal.test(fraction_dauer ~ temperature, data = sds_new_dat)

#	Kruskal-Wallis rank sum test
#
#data:  fraction_dauer by temperature
#Kruskal-Wallis chi-squared = 0.014881, df = 1, p-value = 0.9029

#frustratingly, no it does not
