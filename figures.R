
#Code for figures in Hammerschmith, Woodruff, and Phillips 2020 "Opposing directions of stage-specific body length change in a close relative of C. elegans"

#Figure 2A

library(ggplot2)
library(lemon)

dat <- read.table("body_size_developmental_stage.tsv", sep="\t", header=T)

levels(dat$species)[levels(dat$species)=="elegans"] <- "C. elegans"
levels(dat$species)[levels(dat$species)=="inopinata"] <- "C. inopinata"

dat$stage <- factor(dat$stage, levels = c("dauer","L1", "L2", "L3", "L4", "young_adult", "gravid_adult"))

ggplot(dat, aes(x=width, y=length, colour=stage)) + geom_point(size=0.75,alpha=0.9) + facet_rep_wrap(~species,nrow=1) + scale_colour_manual(name="Developmental\nstage", breaks=c("dauer","L1", "L2", "L3", "L4", "young_adult", "gravid_adult"),labels=c("Dauer","L1", "L2", "L3", "L4", "Young adult", "Gravid adult"),values=c("red3", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5","#084594")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12, colour="black"),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=14),legend.title = element_text(colour="black", size=13), legend.text = element_text(colour="black", size = 12),legend.key=element_blank(),strip.text.x = element_text(size=12, face = "italic"),strip.text.y = element_text(size=12), strip.background = element_blank()) + xlab("Width (microns)") + ylab("Length (microns)") + xlim(0,100) + scale_y_continuous(limits=c(0,2000),breaks=c(0,250,500,750,1000,1250,1500,1750,2000)) + guides(colour = guide_legend(override.aes = list(size=1,alpha=1)))




#Figure 2B

library(ggplot2)
library(ggforce)

#renamed "all_species_dauer_size.csv" to "all_species_dauer_size.tsv" for data sharing
dat <- read.table("all_species_dauer_size.tsv", sep="\t", header=T)



#convert to microns
dat$length <- dat$length*1000
dat$width <- dat$width*1000

#get l:w ratio

dat$l_over_w <- dat$length/dat$width

#length
ggplot(dat, aes(x=species, y=length)) + geom_sina(size=0.5,alpha=0.5) + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.y=element_text(size=12, colour="black"), axis.text.x=element_text(size=12, colour="black",face="italic", angle = 45, hjust = 1),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=14),legend.title = element_text(colour="black", size=13), legend.text = element_text(colour="black", size = 12),legend.key=element_blank(),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank()) + xlab("Species") + ylab("Length (microns)") + scale_y_continuous(limits=c(0,700),breaks=c(0,100,200,300,400,500,600,700))


#width
ggplot(dat, aes(x=species, y=width)) + geom_sina(size=0.5,alpha=0.5) + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.y=element_text(size=12, colour="black"), axis.text.x=element_text(size=12, colour="black",face="italic", angle = 45, hjust = 1),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=14),legend.title = element_text(colour="black", size=13), legend.text = element_text(colour="black", size = 12),legend.key=element_blank(),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank()) + xlab("Species") + ylab("Width (microns)") + scale_y_continuous(limits=c(0,50),breaks=c(0,10,20,30,40,50))


#length to width ratio
ggplot(dat, aes(x=species, y=l_over_w)) + geom_sina(size=0.5,alpha=0.5) + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, colour="red",position = position_dodge(width = 0.9)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.y=element_text(size=12, colour="black"), axis.text.x=element_text(size=12, colour="black",face="italic", angle = 45, hjust = 1),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=14),legend.title = element_text(colour="black", size=13), legend.text = element_text(colour="black", size = 12),legend.key=element_blank(),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank()) + xlab("Species") + ylab("Length:width ratio") + scale_y_continuous(limits=c(0,50),breaks=c(0,10,20,30,40,50))



#Figure 3

library(ggplot2)

#formerly "inopinata_sds.csv" , now "inopinata_sds.tsv" for data sharing
dat <- read.table("inopinata_sds.tsv", sep="\t", header=T)

dat$total_worms <- dat$alive + dat$dead

dat$fraction_alive <- dat$alive/dat$total_worms

dat$perc_alive <- dat$fraction_alive*100


ggplot(dat, aes(x= reorder(strain, -fraction_alive, FUN=mean),y=fraction_alive)) + geom_dotplot(binaxis="y", binwidth=.00030, stackdir="center", aes(fill=island, colour=island),alpha=1) + scale_fill_brewer(name="Island", palette="Set1") + scale_colour_brewer(name="Island", palette="Set1") + stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.5,size=0.25) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=11, colour="black"),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=13),legend.title = element_text(colour="black", size=13), legend.text = element_text(colour="black", size = 12),legend.key=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Strain") + ylab("Fraction alive") + guides(colour = guide_legend(override.aes = list(size=1,alpha=1)))




#Supplemental figure 1, sampling map

library(ggplot2)
library(ggmap)

#used google maps will need api key, uncomment line below with a working key
#register_google(key = "")


loc_dat <- read.table("strains_dauer_paper_gps.tsv", header=TRUE, sep="\t")

#Supplemental figure 1B
sbbox <- make_bbox(lon = loc_dat$Longitude, lat = loc_dat$Latitude, f = 1)

sq_map <- get_map(location = sbbox, maptype = "satellite", source = "google")

ggmap(sq_map) + geom_point(data = loc_dat, mapping = aes(x = Longitude, y = Latitude), size=1, colour="white") + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())


#Supplemental figure 1A
sbbox <- make_bbox(lon = loc_dat$Longitude, lat = loc_dat$Latitude, f = 20)


sq_map <- get_map(location = sbbox, maptype = "satellite", source = "google")

ggmap(sq_map) + geom_point(data = loc_dat, mapping = aes(x = Longitude, y = Latitude), size=1, colour="white") + theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())




#Supplemental figure 2


library(ggplot2)
library(ggforce)
library(lemon)
library(reshape2)

#renamed "sds_fraction_alive.csv" to "sds_all_species.tsv"
dat <- read.table("sds_fraction_alive.csv", sep="\t", header=T)

dat$species <- factor(dat$species, levels = c("briggsae","tropicalis","elegans","inopinata"))


levels(dat$species)[levels(dat$species)=="briggsae"] <- "C. briggsae"
levels(dat$species)[levels(dat$species)=="elegans"] <- "C. elegans"
levels(dat$species)[levels(dat$species)=="inopinata"] <- "C. inopinata"
levels(dat$species)[levels(dat$species)=="tropicalis"] <- "C. tropicalis"

ggplot(dat, aes(x=reorder(species, -fraction_alive, FUN=mean), y=fraction_alive)) + geom_sina(size=0.5,alpha=0.75) + stat_summary(aes(group=species),fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25,size=0.25, colour="red",position = position_dodge(width = 0.9)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.y=element_text(size=12, colour="black"), axis.text.x=element_text(size=12, colour="black",face="italic", angle = 45, hjust = 1),axis.ticks = element_line(colour = "black"), axis.title=element_text(size=14),legend.title = element_text(colour="black", size=13), legend.text = element_text(colour="black", size = 12),legend.key=element_blank(),strip.text.x = element_text(size=12),strip.text.y = element_text(size=12), strip.background = element_blank()) + xlab("Species") + ylab("Fraction alive") 
