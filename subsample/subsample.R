#subsampling data to look at the effects of unequal sample sizes. 

library(dplyr)
library(ggplot2)
library(lemon)
library(cowplot)
library(effsize)
library(ggforce)
library(reshape2)

#get data in there
dat_o <- read.table("original_length_width_data.tsv", sep="\t", header=T)

#get factor levels right
dat_o$species <- factor(dat_o$species, levels = c("C. elegans","C. inopinata", "C. briggsae", "C. tropicalis"))
dat_o$stage <- factor(dat_o$stage, levels = c("Dauer","L1", "L2", "L3", "L4", "Young adult", "Gravid adult"))
dat_o$species.stage <- factor(dat_o$species.stage, levels = c("C. inopinata Dauer", "C. elegans Dauer", "C. briggsae Dauer", "C. tropicalis Dauer", "C. elegans L1", "C. elegans L2", "C. elegans L3", "C. elegans L4", "C. elegans Young adult", "C. elegans Gravid adult", "C. inopinata L1", "C. inopinata L2", "C. inopinata L3", "C. inopinata L4", "C. inopinata Young adult", "C. inopinata Gravid adult"))


#combine adults to get N >= 20 for each developmental stage among elegans and inopinata
levels(dat_o$stage)[levels(dat_o$stage)=="Young adult"] <- "Adult"
levels(dat_o$stage)[levels(dat_o$stage)=="Gravid adult"] <- "Adult"
dat_o$stage <- droplevels(dat_o$stage)
dat_o$species.stage <- paste(dat_o$species,dat_o$stage)
dat_o$species.stage <- as.factor(dat_o$species.stage)
dat_o$stage <- factor(dat_o$stage, levels = c("Dauer","L1", "L2", "L3", "L4", "Adult"))

#reproducible random samples

set.seed(54321)
#get sample sizes for all groups

count(dat_o,species.stage)

#         species.stage   n
#1    C. briggsae Dauer  85
#2     C. elegans Adult  70
#3     C. elegans Dauer 113
#4        C. elegans L1  40
#5        C. elegans L2  40
#6        C. elegans L3  36
#7        C. elegans L4  41
#8   C. inopinata Adult  48
#9   C. inopinata Dauer  97
#10     C. inopinata L1  20
#11     C. inopinata L2  22
#12     C. inopinata L3  40
#13     C. inopinata L4  33
#14 C. tropicalis Dauer  93

#randomly subsample 20 data points per group 100 times
sampled_dat_list <- lapply(1:100, function(i) dat_o %>% group_by(species.stage) %>% slice_sample(n=20))

#before, we did best k is 7.

#define functions that get the cluster classification for all species-stage groups

#extract classifications for each stage and put in a df
KCOUNT = function(species1,data){
	ci <- data[data$species == species1,]
	#get classification of all stages
	ci_stage_count <- count(ci,stage)
	ci_L1_k_count <- count(ci[ci$stage== "L1",],k7)
	ci_L2_k_count <- count(ci[ci$stage== "L2",],k7)
	ci_L3_k_count <- count(ci[ci$stage== "L3",],k7)
	ci_L4_k_count <- count(ci[ci$stage== "L4",],k7)
	ci_adult_k_count <- count(ci[ci$stage== "Adult",],k7)
	ci_dauer_k_count <- count(ci[ci$stage== "Dauer",],k7)
	#name stages in df
	ci_L1_k_count$stage <- "L1"
	ci_L2_k_count$stage <- "L2"
	ci_L3_k_count$stage <- "L3"
	ci_L4_k_count$stage <- "L4"
	ci_adult_k_count$stage <- "Adult"
	ci_dauer_k_count$stage <- "Dauer"
	#cat df's
	ci2<- rbind(ci_L1_k_count,ci_L2_k_count,ci_L3_k_count,ci_L4_k_count,ci_adult_k_count,ci_dauer_k_count)
	#add species id
	ci2$species <- species1

	return(ci2)
}

#do pair-wise chi-square and fisher exact tests of k means classifications among all species.stage groups, get stats and put in a df
CS_func_ii = function(x,y){

	grp1 <- x[1]

	grp2 <- x[2]

	k7_df <- y[y$species.stage == grp1 | y$species.stage == grp2,]
	
	k7_df$species.stage <- droplevels(k7_df$species.stage)
	k7_df$k7 <- droplevels(k7_df$k7)
	
	k7Frequencies <- xtabs(n ~ k7 + species.stage, data = k7_df)

	print(k7Frequencies)

	if(nrow(k7Frequencies)==1){

	p.df <- rbind(p.df,data.frame(grp1 = grp1, grp2 = grp2, chisq_stat=NA,chisq_p=NA,fisher_p=NA))

	print(p.df)

	} else { 
	
	p.df <- rbind(p.df,data.frame(grp1 = grp1, grp2 = grp2, chisq_stat=as.numeric(chisq.test(k7Frequencies)$statistic),chisq_p=as.numeric(chisq.test(k7Frequencies)$p.value),fisher_p=fisher.test(k7Frequencies)$p.value))

	print(p.df)}
}

#get a df with the species.stage group pairs
species.stage.pairs <- t(combn(unique(dat_o$species.stage),2))

species.stage.pairs.df <- data.frame(grp1=species.stage.pairs[,1], grp2=species.stage.pairs[,2])

species.stage.pairs.df$grp1 <- as.character(species.stage.pairs.df$grp1)

species.stage.pairs.df$grp2 <- as.character(species.stage.pairs.df$grp2)


#the massive function that does permanova, wilcoxon rank sum tests, k means clustering, chisquare/fishers exact tests for all subsamples
get_clusters = function(df){

	#permanova
	size_matrix <- df[, c("length", "width")]
	pmnva <- adonis(size_matrix ~ df$species + df$stage + df$species:df$stage, method = "euclidean")
	pmnva.df <- tidy(pmnva$aov.tab)
		#above is part one of list
	post_hoc_pairwise_permanova <- pairwise.adonis(size_matrix,df$species.stage,p.adjust.m="BH")
		#above is part two of list
	post_hoc_pairwise_permanova.ci_d_index <- grep("C. inopinata Dauer", post_hoc_pairwise_permanova$pairs)
	post_hoc_pairwise_permanova.ci_d <- post_hoc_pairwise_permanova[post_hoc_pairwise_permanova.ci_d_index,]
	ci.dau.pairs.pmva.ns <- post_hoc_pairwise_permanova.ci_d[which(post_hoc_pairwise_permanova.ci_d$sig == ""),]
	ns.ci.dau.pmnva.comp <- as.character(droplevels(ci.dau.pairs.pmva.ns$pairs))
		#above is part three of list
	num.nsig <- sum(post_hoc_pairwise_permanova[post_hoc_pairwise_permanova.ci_d_index,]$sig == "", na.rm = TRUE)
	tot.comp <- length(post_hoc_pairwise_permanova[post_hoc_pairwise_permanova.ci_d_index,]$sig)
	fra.sig <- 1- (num.nsig/tot.comp)
	ci.dau.pmva.comp.df <- data.frame(num.sig = num.sig, tot.comp=tot.comp, fra.sig=fra.sig)
		#above is part four of list

	#pairwise wilcoxon

	l.wil.df <-  as.data.frame(as.data.frame(df) %>% wilcox_test(length ~ species.stage,p.adjust.method="BH"))
		#above is part five  of list

	w.wil.df <- as.data.frame(as.data.frame(df) %>% wilcox_test(width ~ species.stage,p.adjust.method="BH"))
		#above is part six  of list



	#get just length and width
	dat_lw <- data.frame(length=df$length,width=df$width)
	#normalize
	dat_lws <- scale(dat_lw)
	#do k means clustering
	k7 <- kmeans(dat_lws, centers = 7,iter.max = 20, nstart = 25)
	#call KCOUNT function to get classifications
	kdat <- cbind(data.frame(df, k7 = k7$cluster))
	ci.k.count <- KCOUNT("C. inopinata", kdat)
	ce.k.count <- KCOUNT("C. elegans", kdat)
	#do for dauer for tropicalis and briggsae
	cb <- kdat[kdat$species=="C. briggsae",]
	ct <- kdat[kdat$species=="C. tropicalis",]
	
	cb_stage_count <- count(cb,stage)
	cb_dauer_k_count <- count(cb[cb$stage== "Dauer",],k7)
	cb_dauer_k_count$stage <- "Dauer"
	cb_dauer_k_count$species <- "C. briggsae"
	
	ct_stage_count <- count(ct,stage)
	ct_dauer_k_count <- count(ct[ct$stage== "Dauer",],k7)
	ct_dauer_k_count$stage <- "Dauer"
	ct_dauer_k_count$species <- "C. tropicalis"
	
	#combine df's
	dat_k7 <- rbind(ci.k.count,ce.k.count,cb_dauer_k_count,ct_dauer_k_count)
		#this is part seven of list
	#get data right, factor levels, data types
	dat_k7$cluster <- as.factor(dat_k7$k7)
	dat_k7$k7 <- as.factor(dat_k7$k7)

	
	dat_k7$stage <- factor(dat_k7$stage, levels = c("Dauer","L1", "L2", "L3", "L4", "Adult"))
	
	
	dat_k7$species <- factor(dat_k7$species, levels = c("C. elegans","C. inopinata", "C. briggsae", "C. tropicalis"))

	dat_k7$species.stage <- paste(dat_k7$species,dat_k7$stage)
	dat_k7$species.stage <- factor(dat_k7$species.stage, levels = c("C. inopinata Dauer", "C. elegans Dauer", "C. briggsae Dauer", "C. tropicalis Dauer", "C. elegans L1", "C. elegans L2", "C. elegans L3", "C. elegans L4", "C. elegans Adult", "C. inopinata L1", "C. inopinata L2", "C. inopinata L3", "C. inopinata L4", "C. inopinata Adult"))
	
	#pairwise chisquare, fishers exact tests for classifications

	CS_out <- apply(species.stage.pairs.df, 1, CS_func_ii,y=dat_k7)


	CS_out_df <- bind_rows(CS_out)
		#this is part eight
	a_list <- list(pmnva.df,post_hoc_pairwise_permanova,ns.ci.dau.pmnva.comp,ci.dau.pmva.comp.df,l.wil.df,w.wil.df,dat_k7,CS_out_df)
	
	subs_cluster_list <- append(subs_cluster_list, a_list)

}


#make an object to put the classifications of subsamples in
subs_cluster_list <- NULL
#call the functions for all subsamples!!
subs_cluster_list <- lapply(sampled_dat_list,get_clusters)

#extract global permanova stats
glob_pmva_list <- lapply(subs_cluster_list, `[`, 1)
glob_pmva_df <- bind_rows(glob_pmva_list, .id = "column_label")

#extract posthoc pairwise permanova statistics

pair_pmva_list <- lapply(subs_cluster_list, `[`, 2)

pair_pmva_df <- bind_rows(pair_pmva_list, .id = "column_label")

#extract non-significant C. inopinata Dauer permanova pairwise comparisons 

Ci_dau_pair_ns_pmva_list <- lapply(subs_cluster_list, `[`, 3)
names(Ci_dau_pair_ns_pmva_list) <- 1:100
Ci_dau_pair_ns_pmva_df <- do.call(rbind.data.frame, Ci_dau_pair_ns_pmva_list)
colnames(Ci_dau_pair_ns_pmva_df)[1] <- "comparison"
Ci_dau_pair_ns_pmva_df$subsample.index <- rownames(Ci_dau_pair_ns_pmva_df)



#fraction significant permanova 

fra_sig_pmva_list <- lapply(subs_cluster_list, `[`, 4)

fra_sig_pmva_df <- bind_rows(fra_sig_pmva_list, .id = "column_label")



#wilcoxon length

wilcox_length_list <- lapply(subs_cluster_list, `[`, 5)

wilcox_length_list_df <- bind_rows(wilcox_length_list, .id = "column_label")


#wilcoxon width

wilcox_width_list <- lapply(subs_cluster_list, `[`, 6)

wilcox_width_list_df <- bind_rows(wilcox_width_list, .id = "column_label")


#k means classifications

k_means_classification_list <- lapply(subs_cluster_list, `[`, 7)

k_means_classification_df <- bind_rows(k_means_classification_list, .id = "column_label")


#k means chi square


k_means_chisq_list <- lapply(subs_cluster_list, `[`, 8)

k_means_chisq_df <- bind_rows(k_means_chisq_list, .id = "column_label")


#adjust p values (others adjusted)
k_means_chisq_df$chisq_p_adjust <- p.adjust(k_means_chisq_df$chisq_p,method='BH')
k_means_chisq_df$fisher_p_adjust <- p.adjust(k_means_chisq_df$fisher_p,method='BH')
glob_pmva_df$p_adjust <- p.adjust(glob_pmva_df$p.value,method='BH')



sampled_dat_df <- bind_rows(sampled_dat_list, .id = "column_label")

#prep subsampled data and stats for output


write.table(sampled_dat_df,"subsample_20_input_data.tsv",sep='\t',row.names = FALSE, quote=FALSE)
write.table(glob_pmva_df,"subsample_20_permanova.tsv",sep='\t',row.names = FALSE, quote=FALSE)
write.table(pair_pmva_df,"subsample_20_pairwise_posthoc_permanova.tsv",sep='\t',row.names = FALSE, quote=FALSE)
write.table(Ci_dau_pair_ns_pmva_df,"subsample_20_Ci_dauer_non-significant_pairwise_posthoc_permanova.tsv",sep='\t',row.names = FALSE, quote=FALSE)
write.table(fra_sig_pmva_df,"subsample_20_fraction_Ci_dauer_significant_pairwise_posthoc_permanova.tsv",sep='\t',row.names = FALSE, quote=FALSE)
write.table(wilcox_length_list_df,"subsample_20_pairwise_wilcoxon_length.tsv",sep='\t',row.names = FALSE, quote=FALSE)
write.table(wilcox_width_list_df,"subsample_20_pairwise_wilcoxon_width.tsv",sep='\t',row.names = FALSE, quote=FALSE)
write.table(k_means_classification_df,"subsample_20_kmeans_classifications.tsv",sep='\t',row.names = FALSE, quote=FALSE)
write.table(k_means_chisq_df,"subsample_20_kmeans_pairwise_chisq_fet.tsv",sep='\t',row.names = FALSE, quote=FALSE)



#how often do we see a significant difference in k means classification among inopinata dauers and all other species.stage groups?

#get the non-significant comparisons with C. inopinata dauers, the chisquare/k means framework
get_sig_chisq = function(dat){

	pair_chsq_fet_km.ci_d_index <- grep("C. inopinata Dauer", dat$grp1)
	pair_chsq_fet_km.ci_d <- dat[pair_chsq_fet_km.ci_d_index,]
	pair_chsq_fet_km.ns <- pair_chsq_fet_km.ci_d[which(pair_chsq_fet_km.ci_d$chisq_p_adjust >= 0.05),]
	num.sig <- length(pair_chsq_fet_km.ns$chisq_p_adjust)
	tot.comp <- length(pair_chsq_fet_km.ci_d$chisq_p_adjust)
	fra.sig <- num.sig/tot.comp
	return(fra.sig)

}

fra_ns_chisq_sub <- by(k_means_chisq_df, k_means_chisq_df$column_label, get_sig_chisq)

table(fra_ns_chisq_sub)


#                 0 0.0769230769230769  0.153846153846154  0.230769230769231
#                19                 64                 13                  4

#for 19% of subsamples, C. inopinata dauer can be distinguished from all others.
# for 64%, there is one species.stage group that cannot be distinguished from C. inopinata dauer
# 13%, two groups
# 4%, three groups


#which groups cannot be distinguished from C. inopinata dauers?


get_sig_chisq = function(dat){

	pair_chsq_fet_km.ci_d_index <- grep("C. inopinata Dauer", dat$grp1)
	pair_chsq_fet_km.ci_d <- dat[pair_chsq_fet_km.ci_d_index,]
	pair_chsq_fet_km.ns <- pair_chsq_fet_km.ci_d[which(pair_chsq_fet_km.ci_d$chisq_p_adjust >= 0.05),]
	num.sig <- length(pair_chsq_fet_km.ns$chisq_p_adjust)
	tot.comp <- length(pair_chsq_fet_km.ci_d$chisq_p_adjust)
	fra.sig <- num.sig/tot.comp
	return(pair_chsq_fet_km.ns)

}


ns_chisq_sub <- by(k_means_chisq_df, k_means_chisq_df$column_label, get_sig_chisq)
names(ns_chisq_sub) <- 1:100
ns_chisq_sub_df <- do.call(rbind.data.frame, ns_chisq_sub)

ns_chisq_sub_df$grp1.grp2 <- paste(ns_chisq_sub_df$grp1, ns_chisq_sub_df$grp2)

table(ns_chisq_sub_df$grp1.grp2)
#   C. inopinata Dauer C. elegans Dauer       C. inopinata Dauer C. elegans L2
#                                     3                                      7
#      C. inopinata Dauer C. elegans L3     C. inopinata Dauer C. inopinata L2
#                                    24                                     61
#    C. inopinata Dauer C. inopinata L3 C. inopinata Dauer C. tropicalis Dauer
#                                     1                                      6

sum(table(ns_chisq_sub_df$grp1.grp2))
#102

#mostly, it is the C. inopinata L2 that cannot be distinguished from the C. inopinata dauer
#next most common, C. elegans L3 and C. inopinata dauer

#inopinata dauer cannot be distinguished from elegans dauer in 3/100 subsamples.
#inopinata dauer cannot be distinguished from tropicalis dauer in 6/100 subsamples.
#inopinata dauer is distinguishable from C. briggsae dauer in all subsamples.



#what about the same as above, but the pairwise permanova framework? how often can C. inopinata dauers not be distinguished from other groups?

get_pair_perm = function(dat){

	post_hoc_pairwise_permanova.ci_d_index <- grep("C. inopinata Dauer", dat$pairs)

	post_hoc_pairwise_permanova.ci_d <- dat[post_hoc_pairwise_permanova.ci_d_index,]

	ci.dau.pairs.pmva.ns <- post_hoc_pairwise_permanova.ci_d[which(post_hoc_pairwise_permanova.ci_d$sig == ""),]

	ns.ci.dau.pmnva.comp <- as.character(droplevels(ci.dau.pairs.pmva.ns$pairs))

	num.nsig <- length(ci.dau.pairs.pmva.ns$sig)
	tot.comp <- length(post_hoc_pairwise_permanova[post_hoc_pairwise_permanova.ci_d_index,]$sig)
	fra.nsig <- num.nsig/tot.comp

	return(fra.nsig)

}

fra_ns_pmva_sub <- by(pair_pmva_df, pair_pmva_df$column_label, get_pair_perm)

table(fra_ns_pmva_sub)

#fra_ns_pmva_sub
#0.0769230769230769  0.153846153846154  0.230769230769231  0.307692307692308
#                31                 49                 19                  1

#in 49% of cases, 2/13 pairs are not significant; 31% 1/13 pair; 19% 3/13 pairs,  1% 4/13 pairs


#which groups cannot be distinguished from C. inopinata dauers?

get_pair_perm = function(dat){

	post_hoc_pairwise_permanova.ci_d_index <- grep("C. inopinata Dauer", dat$pairs)

	post_hoc_pairwise_permanova.ci_d <- dat[post_hoc_pairwise_permanova.ci_d_index,]

	ci.dau.pairs.pmva.ns <- post_hoc_pairwise_permanova.ci_d[which(post_hoc_pairwise_permanova.ci_d$sig == ""),]

	ns.ci.dau.pmnva.comp <- as.character(droplevels(ci.dau.pairs.pmva.ns$pairs))

	num.nsig <- length(ci.dau.pairs.pmva.ns$sig)
	tot.comp <- length(post_hoc_pairwise_permanova.ci_d$sig)
	fra.nsig <- num.nsig/tot.comp

	return(ci.dau.pairs.pmva.ns)

}

ns_pmva_sub <- by(pair_pmva_df, pair_pmva_df$column_label, get_pair_perm)

names(ns_pmva_sub) <- 1:100
ns_pmva_sub_df <- do.call(rbind.data.frame, ns_pmva_sub)

ns_pmva_sub_df$pairs <- droplevels(ns_pmva_sub_df$pairs)

table(ns_pmva_sub_df$pairs)

#C. briggsae Dauer vs C. inopinata Dauer  C. elegans Dauer vs C. inopinata Dauer
#                                     32                                      55
#    C. elegans L2 vs C. inopinata Dauer   C. inopinata Dauer vs C. inopinata L2
#                                      3                                     100

#inopinata dauer continues to be indistinguishable from inopinata L2 in this framework

#most of time, 55%, indistinguishable from elegans dauer
# a substantial fraction of the time, 32%, indistinguishable from briggsae dauer.

#all of this is expected if effect sizes are small and we are subsampling.



#plot dauer classification between elegans and inopinata in all 100 subsamples
#get just dauer classifications for all subsamples
dauer_clust_stage_df <- k_means_classification_df[k_means_classification_df$stage=="Dauer",]
#just get elegans and inopinata dauers
ele_ino_dauer_clust_stage_df <- dauer_clust_stage_df[dauer_clust_stage_df$species == "C. elegans" | dauer_clust_stage_df$species == "C. inopinata",]
#get shorter species label for plotting
ele_ino_dauer_clust_stage_df$species.short <- ele_ino_dauer_clust_stage_df$species

levels(ele_ino_dauer_clust_stage_df$species.short)[levels(ele_ino_dauer_clust_stage_df$species.short)=="C. elegans"] <- "Ce"
levels(ele_ino_dauer_clust_stage_df$species.short)[levels(ele_ino_dauer_clust_stage_df$species.short)=="C. inopinata"] <- "Ci"
#plot it!
ggplot(ele_ino_dauer_clust_stage_df, aes(x = species.short, y = n, fill = cluster)) + geom_col(position="fill",nrow=5) + facet_rep_wrap(~column_label) + scale_fill_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")) + theme_cowplot() + theme(strip.background = element_blank(),strip.text.x = element_blank(),axis.text.x = element_text(face = "italic"))  + xlab("Species") + ylab("Fraction of animals") + scale_y_continuous(breaks=c(0,0.5,1))
#this is supplemental figure 8
