#this is the revised code for making figures associated with manuscript, 
#"Opposing directions of stage-specific body length change in a close relative of C. elegans"
#by Hammerschmith, Woodruff, Johnson and Phillips
#2021

#email Gavin at gcwoodruff@ou.edu if you have any questions

#load libraries

library(tidyverse)
library(ggplot2)
library(cowplot)
library(lemon)
library(ggforce)
library(patchwork)
library(reshape2)
library(ggtree)
library(ape)
library(ggmap)

#get data in there

dat_o <- read.table("original_length_width_data.tsv", sep="\t", header=T)

#get factor levels right
dat_o$species <- factor(dat_o$species, levels = c("C. elegans","C. inopinata", "C. briggsae", "C. tropicalis"))
dat_o$stage <- factor(dat_o$stage, levels = c("Dauer","L1", "L2", "L3", "L4", "Young adult", "Gravid adult"))
dat_o$species.stage <- factor(dat_o$species.stage, levels = c("C. inopinata Dauer", "C. elegans Dauer", "C. briggsae Dauer", "C. tropicalis Dauer", "C. elegans L1", "C. elegans L2", "C. elegans L3", "C. elegans L4", "C. elegans Young adult", "C. elegans Gravid adult", "C. inopinata L1", "C. inopinata L2", "C. inopinata L3", "C. inopinata L4", "C. inopinata Young adult", "C. inopinata Gravid adult"))


#figure 2 -- phylogeny


	#this is the bayesian tree from Stevens et al. 2019 ; retrieved from (https://zenodo.org/record/1402254)
tree <- read.tree("PhyloBayes_species_tree.nwk")

tree_1 <- keep.tip(tree,c("CKAMA","CELEG","CSP34","CWALL","CTROP","CDOUG","CBREN","CREMA","CLATE","CSINI","CSP40","CSP26","CBRIG","CNIGO"))

#root the tree

tree_2 <- root(tree_1, outgroup= "CKAMA")

#new tip labels

tree_2$tip.label <- c("C. kamaaina","C. inopinata","C. elegans","C. remanei","C. latens", "C. tribulationis", "C. zanzibari", "C. sinica","C. nigoni","C. briggsae","C. wallacei","C. tropicalis","C. doughertyi","C. brenneri")


ggtree(tree_2,branch.length=0.75) + geom_tiplab() + geom_treescale() + xlim(0,1)



#figure 2 -- scatterplot

ce_ci <- dat_o[dat_o$species == "C. elegans" | dat_o$species == "C. inopinata",]

ce_ci$species <- droplevels(ce_ci$species)

ggplot(ce_ci, aes(x=width, y=length, colour=stage)) + geom_point(size=0.75,alpha=0.9) + facet_rep_wrap(~species,nrow=1) + scale_colour_manual(name="Developmental\nstage", values=c("red3", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5","#084594")) + theme_cowplot() + theme(strip.text.x = element_text(size=12, face = "italic"),strip.background = element_blank()) + xlab("Width (microns)") + ylab("Length (microns)") + xlim(0,100) + scale_y_continuous(limits=c(0,2000),breaks=c(0,250,500,750,1000,1250,1500,1750,2000)) + guides(colour = guide_legend(override.aes = list(size=1,alpha=1))) 


#figure 2 -- sina

dauer_dat <- dat_o[dat_o$stage=="Dauer",]
dauer_dat$species <- factor(dauer_dat$species, levels = c("C. briggsae","C. tropicalis","C. elegans","C. inopinata"))


	#length sina
ggplot(dauer_dat, aes(x=species, y=length)) + geom_sina(size=0.5,alpha=0.5) + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme_cowplot() + theme(axis.text.x=element_text(size=12, colour="black",face="italic", angle = 45, hjust = 1)) + xlab("Species") + ylab("Length (microns)") + scale_y_continuous(limits=c(0,700),breaks=c(0,100,200,300,400,500,600,700))

	#width sina
ggplot(dauer_dat, aes(x=species, y=width)) + geom_sina(size=0.5,alpha=0.5) + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme_cowplot() + theme(axis.text.x=element_text(size=12, colour="black",face="italic", angle = 45, hjust = 1)) + xlab("Species") + ylab("Width (microns)") + scale_y_continuous(limits=c(0,50),breaks=c(0,10,20,30,40,50))



#figure 3 -- inopinata variation in dauer formation frequency

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

#inopinata only
ino_sds_dat <- sds_dat[sds_dat$Species == "C. inopinata",]

ino_sds_dat$Species <- droplevels(ino_sds_dat$Species)
ino_sds_dat$Strain <- droplevels(ino_sds_dat$Strain)

ino_sds_dat$perc_alive <- ino_sds_dat$fraction_dauer*100

ggplot(ino_sds_dat, aes(x= reorder(Strain, -fraction_dauer, FUN=mean),y=fraction_dauer)) + geom_dotplot(binaxis="y", binwidth=.00030, stackdir="center", aes(fill=Island, colour=Island),alpha=1) + scale_fill_brewer(name="Island", palette="Set1") + scale_colour_brewer(name="Island", palette="Set1") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5,size=0.25) + theme_cowplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Strain") + ylab("Fraction dauer") + guides(colour = guide_legend(override.aes = list(size=1,alpha=1))) + scale_y_continuous(limits=c(0,0.02),breaks=c(0,0.005,0.01,0.015,0.02))




#supplemental figure 1
#length sina
ggplot(dat_o, aes(x = stage, y = length)) + geom_sina(aes(colour=species),size=0.75, alpha=0.75,scale="width") + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442")) + theme_cowplot() + xlab("Stage") + ylab("Length (microns)") + labs(colour="Species") + theme(legend.text = element_text(face = "italic")) + scale_y_continuous(limits=c(0,2000),breaks=c(0,250,500,750,1000,1250,1500,1750,2000)) + scale_x_discrete(labels=c(c("Dauer","L1", "L2", "L3", "L4", "Young\nadult", "Gravid\nadult"))) + guides(colour = guide_legend(override.aes = list(size=2,alpha=1))) 



#supplemental figure 2 effect sizes (see statistics.R for how to generate these df's)

length_effsize_df <- read.table("sheet_5_length_pairwise_effect_sizes.tsv", sep="\t", header=T)
width_effsize_df <- read.table("sheet_6_width_pairwise_effect_sizes.tsv", sep="\t", header=T)

#organize df to plot
length_effsize_df[length_effsize_df$group1 == "C. inopinata L2" & length_effsize_df$group2 == "C. elegans L2",]

dauer_length_effsize <- data.frame(stage= "Dauer",effsize= length_effsize_df[length_effsize_df$group1 == "C. inopinata Dauer" & length_effsize_df$group2 == "C. elegans Dauer",]$effsize, upper=length_effsize_df[length_effsize_df$group1 == "C. inopinata Dauer" & length_effsize_df$group2 == "C. elegans Dauer",]$conf.high,lower=length_effsize_df[length_effsize_df$group1 == "C. inopinata Dauer" & length_effsize_df$group2 == "C. elegans Dauer",]$conf.low,dimension="Length")
	#because direction in df is different in pairwise comparisons we need to multiply by negative 1 to make things comparable among groups!!! Also flip upper and lower CI!!
L1_length_effsize <- data.frame(stage= "L1",effsize= length_effsize_df[length_effsize_df$group1 == "C. elegans L1" & length_effsize_df$group2 == "C. inopinata L1",]$effsize*-1, lower=length_effsize_df[length_effsize_df$group1 == "C. elegans L1" & length_effsize_df$group2 == "C. inopinata L1",]$conf.high*-1,upper=length_effsize_df[length_effsize_df$group1 == "C. elegans L1" & length_effsize_df$group2 == "C. inopinata L1",]$conf.low*-1,dimension="Length")
L2_length_effsize <- data.frame(stage= "L2",effsize= length_effsize_df[length_effsize_df$group1 == "C. elegans L2" & length_effsize_df$group2 == "C. inopinata L2",]$effsize*-1, lower=length_effsize_df[length_effsize_df$group1 == "C. elegans L2" & length_effsize_df$group2 == "C. inopinata L2",]$conf.high*-1,upper=length_effsize_df[length_effsize_df$group1 == "C. elegans L2" & length_effsize_df$group2 == "C. inopinata L2",]$conf.low*-1,dimension="Length")
L3_length_effsize <- data.frame(stage= "L3",effsize= length_effsize_df[length_effsize_df$group1 == "C. elegans L3" & length_effsize_df$group2 == "C. inopinata L3",]$effsize*-1, lower=length_effsize_df[length_effsize_df$group1 == "C. elegans L3" & length_effsize_df$group2 == "C. inopinata L3",]$conf.high*-1,upper=length_effsize_df[length_effsize_df$group1 == "C. elegans L3" & length_effsize_df$group2 == "C. inopinata L3",]$conf.low*-1,dimension="Length")
L4_length_effsize <- data.frame(stage= "L4",effsize= length_effsize_df[length_effsize_df$group1 == "C. elegans L4" & length_effsize_df$group2 == "C. inopinata L4",]$effsize*-1, lower=length_effsize_df[length_effsize_df$group1 == "C. elegans L4" & length_effsize_df$group2 == "C. inopinata L4",]$conf.high*-1,upper=length_effsize_df[length_effsize_df$group1 == "C. elegans L4" & length_effsize_df$group2 == "C. inopinata L4",]$conf.low*-1,dimension="Length")
young_adult_length_effsize <- data.frame(stage= "Young adult",effsize= length_effsize_df[length_effsize_df$group1 == "C. elegans Young adult" & length_effsize_df$group2 == "C. inopinata Young adult",]$effsize*-1, lower=length_effsize_df[length_effsize_df$group1 == "C. elegans Young adult" & length_effsize_df$group2 == "C. inopinata Young adult",]$conf.high*-1,upper=length_effsize_df[length_effsize_df$group1 == "C. elegans Young adult" & length_effsize_df$group2 == "C. inopinata Young adult",]$conf.low*-1,dimension="Length")
gravid_adult_length_effsize <- data.frame(stage= "Gravid adult",effsize= length_effsize_df[length_effsize_df$group1 == "C. elegans Gravid adult" & length_effsize_df$group2 == "C. inopinata Gravid adult",]$effsize*-1, lower=length_effsize_df[length_effsize_df$group1 == "C. elegans Gravid adult" & length_effsize_df$group2 == "C. inopinata Gravid adult",]$conf.high*-1,upper=length_effsize_df[length_effsize_df$group1 == "C. elegans Gravid adult" & length_effsize_df$group2 == "C. inopinata Gravid adult",]$conf.low*-1,dimension="Length")


len_effsize_plot_df <- rbind.data.frame(dauer_length_effsize,L1_length_effsize,L2_length_effsize,L3_length_effsize,L4_length_effsize,young_adult_length_effsize,gravid_adult_length_effsize)

#width
dauer_width_effsize <- data.frame(stage= "Dauer",effsize= width_effsize_df[width_effsize_df$group1 == "C. inopinata Dauer" & width_effsize_df$group2 == "C. elegans Dauer",]$effsize, upper=width_effsize_df[width_effsize_df$group1 == "C. inopinata Dauer" & width_effsize_df$group2 == "C. elegans Dauer",]$conf.high,lower=width_effsize_df[width_effsize_df$group1 == "C. inopinata Dauer" & width_effsize_df$group2 == "C. elegans Dauer",]$conf.low,dimension="Width")
	#because direction in df is different in pairwise comparisons we need to multiply by negative 1 to make things comparable among groups!!! Also flip upper and lower!!
L1_width_effsize <- data.frame(stage= "L1",effsize= width_effsize_df[width_effsize_df$group1 == "C. elegans L1" & width_effsize_df$group2 == "C. inopinata L1",]$effsize*-1, lower=width_effsize_df[width_effsize_df$group1 == "C. elegans L1" & width_effsize_df$group2 == "C. inopinata L1",]$conf.high*-1,upper=width_effsize_df[width_effsize_df$group1 == "C. elegans L1" & width_effsize_df$group2 == "C. inopinata L1",]$conf.low*-1,dimension="Width")
L2_width_effsize <- data.frame(stage= "L2",effsize= width_effsize_df[width_effsize_df$group1 == "C. elegans L2" & width_effsize_df$group2 == "C. inopinata L2",]$effsize*-1, lower=width_effsize_df[width_effsize_df$group1 == "C. elegans L2" & width_effsize_df$group2 == "C. inopinata L2",]$conf.high*-1,upper=width_effsize_df[width_effsize_df$group1 == "C. elegans L2" & width_effsize_df$group2 == "C. inopinata L2",]$conf.low*-1,dimension="Width")
L3_width_effsize <- data.frame(stage= "L3",effsize= width_effsize_df[width_effsize_df$group1 == "C. elegans L3" & width_effsize_df$group2 == "C. inopinata L3",]$effsize*-1, lower=width_effsize_df[width_effsize_df$group1 == "C. elegans L3" & width_effsize_df$group2 == "C. inopinata L3",]$conf.high*-1,upper=width_effsize_df[width_effsize_df$group1 == "C. elegans L3" & width_effsize_df$group2 == "C. inopinata L3",]$conf.low*-1,dimension="Width")
L4_width_effsize <- data.frame(stage= "L4",effsize= width_effsize_df[width_effsize_df$group1 == "C. elegans L4" & width_effsize_df$group2 == "C. inopinata L4",]$effsize*-1, lower=width_effsize_df[width_effsize_df$group1 == "C. elegans L4" & width_effsize_df$group2 == "C. inopinata L4",]$conf.high*-1,upper=width_effsize_df[width_effsize_df$group1 == "C. elegans L4" & width_effsize_df$group2 == "C. inopinata L4",]$conf.low*-1,dimension="Width")
young_adult_width_effsize <- data.frame(stage= "Young adult",effsize= width_effsize_df[width_effsize_df$group1 == "C. elegans Young adult" & width_effsize_df$group2 == "C. inopinata Young adult",]$effsize*-1, lower=width_effsize_df[width_effsize_df$group1 == "C. elegans Young adult" & width_effsize_df$group2 == "C. inopinata Young adult",]$conf.high*-1,upper=width_effsize_df[width_effsize_df$group1 == "C. elegans Young adult" & width_effsize_df$group2 == "C. inopinata Young adult",]$conf.low*-1,dimension="Width")
gravid_adult_width_effsize <- data.frame(stage= "Gravid adult",effsize= width_effsize_df[width_effsize_df$group1 == "C. elegans Gravid adult" & width_effsize_df$group2 == "C. inopinata Gravid adult",]$effsize*-1, lower=width_effsize_df[width_effsize_df$group1 == "C. elegans Gravid adult" & width_effsize_df$group2 == "C. inopinata Gravid adult",]$conf.high*-1,upper=width_effsize_df[width_effsize_df$group1 == "C. elegans Gravid adult" & width_effsize_df$group2 == "C. inopinata Gravid adult",]$conf.low*-1,dimension="Width")


wid_effsize_plot_df <- rbind.data.frame(dauer_width_effsize,L1_width_effsize,L2_width_effsize,L3_width_effsize,L4_width_effsize,young_adult_width_effsize,gravid_adult_width_effsize)


effsize_plot_df <- rbind.data.frame(len_effsize_plot_df,wid_effsize_plot_df)

#make a y axis title
y_title <- expression(paste("Cohen's ", italic("d"), " effect size"))

#plot!
ggplot(effsize_plot_df, aes(x = stage, y = effsize)) + geom_col(fill = "#56B4E9") + geom_errorbar(aes(ymin = lower, ymax = upper), width = .1) + facet_rep_wrap(~dimension,ncol=1) + theme_cowplot() + xlab("Stage")+ theme(strip.background=element_blank()) + geom_hline(yintercept=0,linetype="dashed") + scale_x_discrete(labels=c(c("Dauer","L1", "L2", "L3", "L4", "Young\nadult", "Gravid\nadult"))) + labs(y=y_title) + scale_y_continuous(breaks=c(-1:8),limits=c(-1,8))




##supplemental figure 3, new data, tails

#get data in there
dat_n <- read.table("new_data_for_revisions_dauer_elegans_inopinata_tails_sheath.tsv", sep="\t", header=T)

#round and remove tail
dat_n$length_include_tail_micron <- round(dat_n$length_include_tail_micron)
dat_n$tail_length_micron <- round(dat_n$tail_length_micron)
dat_n$length_no_tail <- dat_n$length_include_tail_micron - dat_n$tail_length_micron


a <- ggplot(dat_n, aes(x = species, y = length_include_tail_micron)) + geom_sina(colour="#56B4E9",size=1, alpha=1,scale="width") + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", colour="black") + theme_cowplot() + xlab("Species") + ylab("Total dauer length (microns)") + xlab("Species") + theme(axis.text.x = element_text(face = "italic"),axis.title.x=element_blank()) + scale_y_continuous(limits=c(0,700),breaks=c(0,100,200,300,400,500,600,700)) 

b <- ggplot(dat_n, aes(x = species, y = length_no_tail)) + geom_sina(colour="#56B4E9",size=1, alpha=1,scale="width") + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", colour="black") + theme_cowplot() + xlab("Species") + ylab("Dauer length\nwith tail excluded (microns)") + xlab("Species") + theme(axis.text.x = element_text(face = "italic")) + scale_y_continuous(limits=c(0,700),breaks=c(0,100,200,300,400,500,600,700)) 

c <- ggplot(dat_n, aes(x = species, y = tail_length_micron)) + geom_sina(colour="#56B4E9",size=1, alpha=1,scale="width") + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", colour="black") + theme_cowplot() + xlab("Species") + ylab("Dauer tail length (microns)") + xlab("Species") + theme(axis.text.x = element_text(face = "italic"),axis.title.x=element_blank()) + scale_y_continuous(limits=c(0,60),breaks=c(0,10,20,30,40,50,60)) 

d <- ggplot(dat_n, aes(x = species, y = width_exclude_sheath_micron)) + geom_sina(colour="#56B4E9",size=1, alpha=1,scale="width") + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", colour="black") + theme_cowplot() + xlab("Species") + ylab("Width (microns)") + xlab("Species") + theme(axis.text.x = element_text(face = "italic"),axis.title.x=element_blank()) + scale_y_continuous(limits=c(0,30),breaks=c(0,5,10,15,20,25,30)) 


(a+b)/(c+d)




#supplemental figure 4
#width sina
ggplot(dat_o, aes(x = stage, y = width)) + geom_sina(aes(colour=species),size=0.75, alpha=0.75,scale="width") + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442")) + theme_cowplot() + xlab("Stage") + ylab("Width (microns)") + labs(colour="Species") + theme(legend.text = element_text(face = "italic")) + scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100)) + scale_x_discrete(labels=c(c("Dauer","L1", "L2", "L3", "L4", "Young\nadult", "Gravid\nadult"))) + guides(colour = guide_legend(override.aes = list(size=2,alpha=1))) 




#supplemental figure 5, relationship of dauers to non-dauer stages in length-width space


#get elegans reproductive stages
ce_repr <- dat_o[dat_o$species == "C. elegans" & dat_o$stage != "Dauer",]
#get length width fit
ce_repr_fit <- lm(ce_repr$length~ce_repr$width)
#get elegans dauers
ce_dauer <- dat_o[dat_o$species == "C. elegans" & dat_o$stage == "Dauer",]
#are elegans dauers above fit?
ce_dauer$is_above_regression <- ifelse(ce_dauer$length > (ce_repr_fit$coefficients[2]*ce_dauer$width) + ce_repr_fit$coefficients[1],'Dauer above','Dauer below')


#get inopinata reproductive stages
ci_repr <- dat_o[dat_o$species == "C. inopinata" & dat_o$stage != "Dauer",]
#get length width fit
ci_repr_fit <- lm(ci_repr$length~ci_repr$width)
#get inopinata dauers
ci_dauer <- dat_o[dat_o$species == "C. inopinata" & dat_o$stage == "Dauer",]
#are inopinata dauers above fit?
ci_dauer$is_above_regression <- ifelse(ci_dauer$length > (ci_repr_fit$coefficients[2]*ci_dauer$width) + ci_repr_fit$coefficients[1],'Dauer above','Dauer below')
#get labels right
ce_repr$is_above_regression <- 'Reproductive cycle stages'
ci_repr$is_above_regression <- 'Reproductive cycle stages'

#cat data for plots
fdat <- rbind(ce_repr,ce_dauer,ci_repr,ci_dauer)
#scatterplot colored by stage and relationship to fit
a <- ggplot(fdat, aes(x=width, y=length, colour=is_above_regression)) + geom_point(alpha=1,size=0.75) + geom_smooth(data = subset(dat_o, stage != "Dauer"),se=FALSE,method=lm,linetype="dotted",colour="black") + facet_rep_wrap(~species,nrow=1)  + theme_cowplot() + scale_colour_manual(values=c("red1","rosybrown","steelblue3")) + theme(strip.text.x = element_text(face = "italic",size=14),axis.title=element_text(size=16),axis.text=element_text(size=14),legend.title=element_text(size=16), legend.text=element_text(size=14), strip.background = element_blank()) + xlab("Width (microns)") + ylab("Length (microns)") + xlim(0,100) + scale_y_continuous(limits=c(0,2000),breaks=c(0,250,500,750,1000,1250,1500,1750,2000)) + labs(colour="Dauer above fit?") + guides(colour = guide_legend(override.aes = list(size=2)))

#get bar plot data for above/below fit classifcation of dauers
ci_is_above_count <- count(ci_dauer,is_above_regression)
ci_is_above_count$species <- "C. inopinata"

ce_is_above_count <- count(ce_dauer,is_above_regression)
ce_is_above_count$species <- "C. elegans"

dat_is_above <- rbind(ce_is_above_count,ci_is_above_count)

#bar plot
b <- ggplot(dat_is_above, aes(x = species, y = n, fill = is_above_regression)) + geom_col() + scale_fill_manual(values=c("red1","rosybrown")) + theme_cowplot() + theme(strip.text.x = element_text(face = "italic",size=14),axis.title=element_text(size=16),axis.text=element_text(size=14),legend.title=element_text(size=16), legend.text=element_text(size=14), strip.background = element_blank()) + ylim(0,120) + xlab("Species") + ylab("Number of animals") + labs(fill="Dauer above fit?")
#combine plots
(a/b) +  plot_layout(heights = c(3, 1))
#stats and fit equation added in adobe illustrator


#supplemental figure 6 and supplemental figure 7 code is in statistics.R

#supplemental figure 8 code is in kmeans_subsample.R

#supplemental figure 9 is in statistics.R


#supplemental figure 10 


#data
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


ggplot(sds_new_dat, aes(x=temperature, y=fraction_dauer)) + geom_sina(size=0.5,alpha=0.5) + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme_cowplot() + xlab("Species") + ylab("Fraction dauer") 

sds_new_dat$temperature <- as.factor(sds_new_dat$temperature)

geom_point(aes(colour=strain))

ggplot(sds_new_dat, aes(x = temperature, y = fraction_dauer)) + geom_jitter(aes(colour=strain,width = 0.9),size=1, alpha=1,scale="width") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("skyblue","tan1")) + theme_cowplot() + xlab("Temperature (°C)") + ylab("Fraction dauer") + labs(colour="Strain") + scale_y_continuous(limits=c(-0.00002,0.007)) + guides(colour = guide_legend(override.aes = list(size=2,alpha=1))) 



#Supplemental figure 11 was made with an academic subscription to google maps in 2019.


# register_google(key = "")


loc_dat <- read.table("strains_dauer_paper_gps.tsv", header=TRUE, sep="\t")


sbbox <- make_bbox(lon = loc_dat$Longitude, lat = loc_dat$Latitude, f = 1)




sq_map <- get_map(location = sbbox, maptype = "satellite", source = "google")

ggmap(sq_map) + geom_point(data = loc_dat, mapping = aes(x = Longitude, y = Latitude), size=1, colour="white") + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())




sbbox <- make_bbox(lon = loc_dat$Longitude, lat = loc_dat$Latitude, f = 20)


sq_map <- get_map(location = sbbox, maptype = "satellite", source = "google")

ggmap(sq_map) + geom_point(data = loc_dat, mapping = aes(x = Longitude, y = Latitude), size=1, colour="white") + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())



#supplemental figure 12


ggplot(sds_dat, aes(x=reorder(Species, -fraction_dauer, FUN=mean), y=fraction_dauer)) + geom_sina(size=0.5,alpha=0.75) + stat_summary(aes(group=Species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25,size=0.25, colour="red",position = position_dodge(width = 0.9)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.y=element_text(size=12, colour="black"), axis.text.x=element_text(size=12, colour="black",face="italic", angle = 45, hjust = 1),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=14),legend.title = element_text(colour="black", size=13), legend.text = element_text(colour="black", size = 12),legend.key=element_blank(),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank()) + xlab("Species") + ylab("Fraction dauer") 



#supplemental figure 13, length-width ratio

dat_o$length.width.ratio <- dat_o$length/dat_o$width


ggplot(dat_o, aes(x = stage, y = length.width.ratio)) + geom_sina(aes(colour=species),size=0.75, alpha=0.75,scale="width") + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", colour="black",position = position_dodge(width = 0.9)) + scale_colour_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442")) + theme_cowplot() + xlab("Stage") + ylab("Length:width ratio") + labs(colour="Species") + theme(legend.text = element_text(face = "italic")) + scale_y_continuous(limits=c(0,50),breaks=c(0,10,20,30,40,50)) + scale_x_discrete(labels=c(c("Dauer","L1", "L2", "L3", "L4", "Young\nadult", "Gravid\nadult"))) + guides(colour = guide_legend(override.aes = list(size=1.5,alpha=1)))

