######## AA Figure plots


library(dplyr)
library(ggplot2)
library(stringr)
options(scipen=999)
library(data.table)
library(RcppRoll)
library(TTR)
library(scales)

setwd("~/ERB/PDD-2023-02-07_EXP-046/")

## Import subset fits
sig_fits <- read.csv("./final_fits/SUBSET_NLSfits_SIRWT.csv", row.names = 1, stringsAsFactors = F)


pdf("./plots/restime_density.pdf")
ggplot(sig_fits, aes(x=avg_restime, col=""))+
  geom_density()+
  theme_classic()+
  ylab("Peak Density")+
  scale_color_manual(values=c("#006666"))+
  scale_x_continuous(trans = log_trans(), breaks = c(5, 10, 20, 30, 45, 60), name="Residence Time (min)")
dev.off()


pdf("./plots/restime_histogram.pdf")
ggplot(sig_fits, aes(x=avg_restime, fill=""))+
  geom_histogram(bins=100)+
  theme_classic()+
  ylab("Number of peaks")+
  scale_fill_manual(values=c("#006666", "#99CCCC","#666699"))+
  scale_x_continuous(name="Residence Time (min)", limits=c(0,34))
dev.off()


pdf("./plots/restime_histogram_categories.pdf")
ggplot(sig_fits, aes(x=avg_restime, fill=cat))+
  geom_histogram(bins=100)+
  theme_classic()+
  ylab("Number of peaks")+
  scale_fill_manual(values=c("#006666", "#99CCCC","#666699"))+
  scale_x_continuous(name="Residence Time (min)", limits=c(0,34))
dev.off()


########################################################################
#check replicate fits for correlation


pdf("./plots/koff_replicates.pdf")
ggplot(sig_fits, aes(x=log_koff.237, y=log_koff.238))+
  geom_point(alpha=0.2, col="#006666")+
  geom_abline(slope=1)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_x_continuous("koff (log) Replicate 1", limits=c(-4,-1.5))+
  scale_y_continuous("koff (log) Replicate 2", limits=c(-4,-1.5))+
  annotate("text", x=-3.75, y=-1.75, label="r = 0.87")

dev.off()


cor.test(sig_fits$log_koff.237, sig_fits$log_koff.238)



#### Quantiles
sig_fits <- sig_fits %>%
  mutate(quantile = ntile(avg_restime, 4))
sig_fits$group[sig_fits$quantile==1] <- "Shortest"
sig_fits$group[sig_fits$quantile==2] <- "Short"
sig_fits$group[sig_fits$quantile==3] <- "Long"
sig_fits$group[sig_fits$quantile==4] <- "Longest"

sig_fits$group <- factor(sig_fits$group, levels = c("Shortest", "Short", "Long", "Longest"))

write.csv(sig_fits, "./subset_peaks/PDD-2023_02_10-SigFits.csv")


### Counts
counts <- read.table("./normalized_counts/SUBSET_normalized.counts.txt", stringsAsFactors = F, check.names = F)
counts$peaknames <- row.names(counts)

occupancy <- merge(sig_fits, counts, by="peaknames") %>% dplyr::select(c("peaknames", "avg_logkoff", "avg_restime", "group", "237_00_IP", "238_00_IP"))

occupancy$avg_occ <- (occupancy$`237_00_IP` + occupancy$`238_00_IP`)/2


## Plot

pdf("./plots/ResTime_vs_Occupancy.pdf")
ggplot(occupancy, aes(x=avg_restime, y=avg_occ, col=group))+
  geom_point()+
  scale_y_log10(name="Rap1 Occupancy (IP time=0)")+
  xlab("Residence Time (min)")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  annotate("text", x=25, y=50, label="r = 0.45")
dev.off()  


cor.test(x=occupancy$avg_restime, y=occupancy$avg_occ)



#######################################################################

#Ribosomal protein genes

rp <- sig_fits
rp$rp <- "Other"
rp$rp[str_detect(rp$watson_genename, "RPS") | str_detect(rp$watson_genename, "RPL")] <- "RP"
rp$rp[str_detect(rp$crick_genename, "RPS") | str_detect(rp$crick_genename, "RPL")] <- "RP"
rp$rp <- factor(rp$rp, levels = c("RP", "Other"))

### plot

pdf("./plots/RP_boxplot.pdf")
ggplot(rp, aes(x=rp, y=avg_restime, fill=rp))+
  geom_boxplot(col="black")+
  theme_classic()+
  scale_fill_manual(values=c( "#006666","#99CCCC"))+
  xlab("Gene Category")+
  ylab("Residence Time (min)")+
  annotate("text", x="Other", y=30, label="p < 2.2e-16 ")+
  annotate("text", x="RP", y=19.5, label="n = 104")+
  annotate("text", x="Other", y=10.5, label="n = 273")
dev.off()


sum(rp$rp!="RP")
wilcox.test(avg_restime ~ rp, data=rp) 

pdf("./plots/RP_boxplot_OCC.pdf")
ggplot(rp, aes(x=rp, y=X237_00_IP, fill=rp))+
  geom_boxplot(col="black")+
  theme_classic()+
  scale_fill_manual(values=c( "#006666","#99CCCC"))+
  xlab("Gene Category")+
  scale_y_log10(name="Rap1 Occupancy (AU)")+
  annotate("text", x="RP", y=10000, label="p < 2.2e-16")
dev.off()
wilcox.test(X237_00_IP ~ rp, data=rp) 





#######################################################################

#Sub-telometric genes

telo <- sig_fits
telo$subtelo <- factor(telo$subtelo, levels=c(TRUE, FALSE))

### plot

pdf("./plots/Telomere_boxplot.pdf")
ggplot(telo, aes(x=subtelo, y=avg_restime, fill=subtelo))+
  geom_boxplot(col="black")+
  theme_classic()+
  scale_fill_manual(values=c( "#006666","#99CCCC"))+
  xlab("Gene Category")+
  ylab("Residence Time (min)")+
  annotate("text", x=TRUE, y=30, label="p = 0.02 ")+
  annotate("text", x=TRUE, y=19.5, label="n = 33")+
  annotate("text", x=FALSE, y=12, label="n = 344")
dev.off()


sum(telo$subtelo==TRUE)
wilcox.test(avg_restime ~ subtelo, data=telo) 


pdf("./plots/Telomere_boxplot_OCC.pdf")
ggplot(telo, aes(x=subtelo, y=X237_00_IP, fill=subtelo))+
  geom_boxplot(col="black")+
  theme_classic()+
  scale_fill_manual(values=c( "#006666","#99CCCC"))+
  xlab("Gene Category")+
  scale_y_log10(name="Rap1 Occupancy (AU)")+
  annotate("text", x=TRUE, y=10000, label="p = 0.04 ")+
  annotate("text", x=TRUE, y=12, label="n = 33")+
  annotate("text", x=FALSE, y=12, label="n = 344")
dev.off()

wilcox.test(X237_00_IP ~ subtelo, data=telo) 


telo_only <- telo[telo$subtelo == TRUE,]
telo_only$distance[telo_only$summit < 15000] <- telo_only$summit[telo_only$summit < 15000]
telo_only$distance[telo_only$summit > 15000] <- telo_only$V2[telo_only$summit > 15000] - telo_only$summit[telo_only$summit > 15000]

pdf("./plots/Distance_Telomere.pdf")
ggplot(telo_only, aes(x=distance, y=avg_restime, col=""))+
  geom_point()+
  theme_classic()+
  xlab("Distance to telomere (bp)")+
  ylab("Average Res. time (min)")+
  annotate("text", x=5000, y=20, label="p = 0.91")+
  annotate("text", x=5000, y=18, label="r = 0.02")+
  scale_color_manual(values=c("#006666"))+
  theme(legend.position = "none")
dev.off()
cor.test(telo_only$distance, telo_only$avg_restime, method="spearman")


pdf("./plots/Distance_Telomere_OCC.pdf")
ggplot(telo_only, aes(x=distance, y=X237_00_IP, col=""))+
  geom_point()+
  theme_classic()+
  xlab("Distance to telomere (bp)")+
  scale_y_log10(name="Rap1 Occupancy (AU)")+
  annotate("text", x=5000, y=30, label="p = 0.06")+
  annotate("text", x=5000, y=20, label="r = 0.33")+
  scale_color_manual(values=c("#006666"))+
  theme(legend.position = "none")
dev.off()
cor.test(telo_only$distance, telo_only$X237_00_IP)

##############################################################

#### H3 coverage


#### For this I have a bedgraph of H3 occupance
h3 <- read.table("/mnt/ingolialab/paige_diamond/ERB/EXP_046/other_data/H3_occupancy.bedgraph", stringsAsFactors = F)


##### I am going to make a bed file of the positions I want to get coverage of (1000 bp up and downstream of summit) 
width=1000
h3_coverage <- data.frame(matrix(NA, nrow = nrow(sig_fits)*width, ncol = 8))
colnames(h3_coverage)<- c("peakname", "chr", "start", "end", "residence", "group", "relative_pos", "h3_cov")
h3_coverage$peakname <- rep(sig_fits$peaknames,1000)
h3_coverage$chr <- rep(sig_fits$chr,1000)
h3_coverage$start <- rep(sig_fits$start,1000)
h3_coverage$end  <- rep(sig_fits$end,1000)
h3_coverage$residence  <- rep(sig_fits$avg_restime,1000)
h3_coverage$group  <- rep(sig_fits$group,1000)
h3_coverage$relative_pos <- rep(-500:499, each=nrow(sig_fits))
h3_coverage$h3_pos <- ifelse(h3_coverage$end -150 + h3_coverage$relative_pos > 0, h3_coverage$end -150 + h3_coverage$relative_pos, NA)

h3_coverage <- h3_coverage[!is.na(h3_coverage$h3_pos),]

################################
##### In the terminal, I am going to use bedtools genomecov to get coverage at each of those positions from the H3 occupancy bedgraph
h3_bed <- data.frame(chr=paste("chr", gsub("_cer", "", h3_coverage$chr), sep=""), start=h3_coverage$h3_pos, end=h3_coverage$h3_pos+1)


write.table(h3_bed, "./other_data/H3_coverage.bed", quote=F, row.names=F, col.names = F, sep="\t")


## bedtools sort -i H3_coverage.bed > H3_coverage_sorted.bed
## bedtools map -a H3_coverage_sorted.bed -b H3_occupancy.bedgraph -c 4 -o mean > H3_coverage_new.bed 

h3_values <- read.table("./other_data/H3_coverage_new.bed")
colnames(h3_values) <- c("chr", "start", "end", "h3_cov")
h3_values$chr <- gsub("chr", "", h3_values$chr)


dt_coverage = data.table(h3_coverage, key=c("chr", "h3_pos"))
dt_h3values = data.table(h3_values, key=c("chr", "start"))

merge_h3 <- merge(dt_coverage, dt_h3values, by.x=c("chr", "h3_pos"), by.y=c("chr", "start"))
merge_h3 <- data.frame(merge_h3)
merge_h3 <- merge_h3[!duplicated(merge_h3[c("chr", "h3_pos")]),]


h3_collapse <- merge_h3 %>% group_by(relative_pos, group) %>% summarise(mean_coverage = mean(h3_cov.y), sd=sd(h3_cov.y), n=n()) %>%  mutate(se = sd / sqrt(n),
                                                                                                                 lower.ci = mean_coverage - qt(1 - (0.05 / 2), n - 1) * se,
                                                                                                                 upper.ci = mean_coverage + qt(1 - (0.05 / 2), n - 1) * se)
pdf("./plots/H3_metagene.pdf")
ggplot(h3_collapse, aes(x=relative_pos, y=mean_coverage))+
  geom_line(aes(col=group))+
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill=group), alpha = 0.1)+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_fill_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_y_continuous("H3 signal (AU)")+
  scale_x_continuous("Position Relative to Peak")
dev.off()


#####################################################

##### Remake H3 plot based on occupancy 

occ_group <- sig_fits

#### Quantiles
occ_group <- occ_group %>%
  mutate(quantile = ntile(X237_00_IP, 4))
occ_group$group <- NA
occ_group$group[occ_group$quantile==1] <- "Least"
occ_group$group[occ_group$quantile==2] <- "Less"
occ_group$group[occ_group$quantile==3] <- "More"
occ_group$group[occ_group$quantile==4] <- "Most"

occ_group$group <- factor(occ_group$group, levels = c("Least", "Less", "More", "Most"))


width=1000
h3_occ_cov <- data.frame(matrix(NA, nrow = nrow(occ_group)*width, ncol = 8))
colnames(h3_occ_cov)<- c("peakname", "chr", "start", "end", "residence", "group", "relative_pos", "h3_cov")
h3_occ_cov$peakname <- rep(occ_group$peaknames,1000)
h3_occ_cov$chr <- rep(occ_group$chr,1000)
h3_occ_cov$start <- rep(occ_group$start,1000)
h3_occ_cov$end  <- rep(occ_group$end,1000)
h3_occ_cov$residence  <- rep(occ_group$avg_restime,1000)
h3_occ_cov$group  <- rep(occ_group$group,1000)
h3_occ_cov$relative_pos <- rep(-500:499, each=nrow(occ_group))
h3_occ_cov$h3_pos <- ifelse(h3_occ_cov$end -150 + h3_occ_cov$relative_pos > 0, h3_occ_cov$end -150 + h3_occ_cov$relative_pos, NA)


dt_coverage = data.table(h3_occ_cov, key=c("chr", "h3_pos"))
dt_h3values = data.table(h3_values, key=c("chr", "start"))

merge_h3 <- merge(dt_coverage, dt_h3values, by.x=c("chr", "h3_pos"), by.y=c("chr", "start"))
merge_h3 <- data.frame(merge_h3)
merge_h3 <- merge_h3[!duplicated(merge_h3[c("chr", "h3_pos")]),]


h3_collapse <- merge_h3 %>% group_by(relative_pos, group) %>% summarise(mean_coverage = mean(h3_cov.y), sd=sd(h3_cov.y), n=n()) %>%  mutate(se = sd / sqrt(n),
                                                                                                                                            lower.ci = mean_coverage - qt(1 - (0.05 / 2), n - 1) * se,
                                                                                                                                            upper.ci = mean_coverage + qt(1 - (0.05 / 2), n - 1) * se)
pdf("./plots/H3_metagene_OCC.pdf")
ggplot(h3_collapse, aes(x=relative_pos, y=mean_coverage))+
  geom_line(aes(col=group))+
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill=group), alpha = 0.1)+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_fill_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_y_continuous("H3 signal (AU)")+
  scale_x_continuous("Position Relative to Peak")
dev.off()



######################################################
##### Distance to ORF

### For divergent promoters, each ORF is considered separately 

corf_repeat <- sig_fits[!is.na(sig_fits$crick),]
corf_repeat$ORF <- corf_repeat$crick

worf_repeat <- sig_fits[!is.na(sig_fits$watson),]
worf_repeat$ORF <- worf_repeat$watson

orf_repeat <- rbind(corf_repeat, worf_repeat)

gtf <- read.table("~/Genomes/Saccharomyces_cerevisiae.R64-1-1.85.gtf", skip=5, sep = '\t')
gtf <- gtf[,c(1,3,4,5,7,9)]
colnames(gtf) <- c("chr", "type", "start", "end", "strand", "name")
gtf <- gtf[gtf$type=="gene",]
gtf$orf <- gsub("\\;.*", "", gtf$name)
gtf$orf <- gsub("gene_id ", "", gtf$orf)

orf_repeat <- orf_repeat %>% select(c("peaknames", "avg_restime", "ORF", "X237_00_IP", "summit", "group")) %>% merge(gtf, by.x = "ORF", by.y="orf")
orf_repeat$orf_distance[orf_repeat$strand=="+"] <- abs(orf_repeat$start[orf_repeat$strand=="+"] - orf_repeat$summit[orf_repeat$strand=="+"])
orf_repeat$orf_distance[orf_repeat$strand=="-"] <- abs(orf_repeat$end[orf_repeat$strand=="-"] - orf_repeat$summit[orf_repeat$strand=="-"])

pdf("ORF_startdistance.pdf")
ggplot(orf_repeat, aes(x=orf_distance, y=avg_restime, col=group))+
  geom_point()+
  ylab("Residence Time (min)")+
  xlab("Distance to ORF start (bp)")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  annotate("text", x=100, y=35, label="r = 0.09")+
  annotate("text", x=100, y=33, label="p = 0.04")
dev.off()

cor.test(orf_repeat$orf_distance, orf_repeat$avg_restime)



##########################################################

### TSS distance 
## Import major transcript isoform coordinates from Pelechano et al. 2013
## This file has TSS and polyA sites for 4928 ORFs 
tss <- read.csv("./other_data/Pelechano_mTIF_anno_file.csv", row.names = 1)
orf_repeat <- rbind(corf_repeat, worf_repeat)


tss_repeat <- orf_repeat %>% select(c("peaknames", "avg_restime", "ORF", "X237_00_IP", "summit", "group")) %>% merge(tss, by.x = "ORF", by.y="orf")
tss_repeat$tss_distance <- abs(tss_repeat$summit - tss_repeat$tss)

pdf("TSS_startdistance.pdf")
ggplot(tss_repeat, aes(x=tss_distance, y=avg_restime, col=group))+
  geom_point()+
  ylab("Residence Time (min)")+
  scale_x_log10(name="Distance to TSS start (bp)")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  annotate("text", x=50, y=35, label="rho = 0.19")+
  annotate("text", x=50, y=33, label="p = 3.5e-4")
dev.off()

cor.test(tss_repeat$tss_distance, tss_repeat$avg_restime, method="spearman")



####################################################

### Net-seq data
#pos_data <- read.table("./other_data/netseq_pos.bedgraph")
#pos_data$cov <- pos_data$V4 + pos_data$V5 + pos_data$V6 + pos_data$V7
#pos_data <- pos_data %>% select(c("V1", "V2", "V3", "cov"))
#write.table(pos_data, "./other_data/netseq_pos.bedgraph", row.names = F, col.names = F, quote=F, sep = '\t')

#neg_data <- read.table("./other_data/netseq_neg.bedgraph")
#neg_data$cov <- neg_data$V4 + neg_data$V5 + neg_data$V6 + neg_data$V7
#neg_data <- neg_data %>% select(c("V1", "V2", "V3", "cov"))
#write.table(neg_data, "./other_data/netseq_neg.bedgraph", row.names = F, col.names = F, quote=F, sep = '\t')


### TSS data from above
wtss_repeat <- tss_repeat[tss_repeat$strand == '+',]
wtss_repeat$chr <- as.character(as.roman(wtss_repeat$chr))
ctss_repeat <- tss_repeat[tss_repeat$strand == '-',]
ctss_repeat$chr <- as.character(as.roman(ctss_repeat$chr))

###### Bed files to count over

##### 500 bp downstream of TSS
width=500

bed <- data.frame(peaknames=NA, orf=NA, avg_restime=NA, group=NA, tss=NA, polya=NA, relative_pos=NA, netseq_pos=NA, chr=NA, ip=NA)

for(peaks in c(1:nrow(wtss_repeat))){
  if(wtss_repeat$polya[peaks] - wtss_repeat$tss[peaks] > width){
    tmp <- data.frame(peaknames=rep(wtss_repeat$peaknames[peaks],width),
                      orf=rep(wtss_repeat$ORF[peaks],width), 
                      avg_restime=rep(wtss_repeat$avg_restime[peaks],width), 
                      group=rep(wtss_repeat$group[peaks],width), 
                      tss=rep(wtss_repeat$tss[peaks],width), 
                      polya=rep(wtss_repeat$polya[peaks],width),
                      relative_pos=c(0:(width-1)),
                      netseq_pos=wtss_repeat$tss[peaks]+c(0:(width-1)),
                      chr=rep(wtss_repeat$chr[peaks],width), 
                      ip=rep(wtss_repeat$X237_00_IP[peaks],width))
  }
  else {
    gene_width <- wtss_repeat$polya[peaks] - wtss_repeat$tss[peaks] 
    tmp <- data.frame(peaknames=rep(wtss_repeat$peaknames[peaks],gene_width),
                      orf=rep(wtss_repeat$ORF[peaks],gene_width), 
                      avg_restime=rep(wtss_repeat$avg_restime[peaks],gene_width), 
                      group=rep(wtss_repeat$group[peaks],gene_width), 
                      tss=rep(wtss_repeat$tss[peaks],gene_width), 
                      polya=rep(wtss_repeat$polya[peaks],gene_width),
                      relative_pos=c(0:(gene_width-1)),
                      netseq_pos=wtss_repeat$tss[peaks]+c(0:(gene_width-1)),
                      chr=rep(wtss_repeat$chr[peaks],gene_width),
                      ip=rep(wtss_repeat$X237_00_IP[peaks],gene_width))
  }
  bed <- rbind(bed,tmp)
}

negative_bed <- data.frame(peaknames=NA, orf=NA, avg_restime=NA, group=NA, tss=NA, polya=NA, relative_pos=NA, netseq_pos=NA, chr=NA, ip=NA)

for(peaks in c(1:nrow(ctss_repeat))){
  if(ctss_repeat$tss[peaks] - ctss_repeat$polya[peaks] > width){
    tmp <- data.frame(peaknames=rep(ctss_repeat$peaknames[peaks],width),
                      orf=rep(ctss_repeat$ORF[peaks],width), 
                      avg_restime=rep(ctss_repeat$avg_restime[peaks],width), 
                      group=rep(ctss_repeat$group[peaks],width), 
                      tss=rep(ctss_repeat$tss[peaks],width), 
                      polya=rep(ctss_repeat$polya[peaks],width),
                      relative_pos=c(0:(width-1)),
                      netseq_pos=ctss_repeat$tss[peaks]-c(0:(width-1)),
                      chr=rep(ctss_repeat$chr[peaks],width), 
                      ip=rep(ctss_repeat$X237_00_IP[peaks],width))
  }
  else {
    gene_width <- ctss_repeat$tss[peaks] - ctss_repeat$polya[peaks] 
    tmp <- data.frame(peaknames=rep(ctss_repeat$peaknames[peaks],gene_width),
                      orf=rep(ctss_repeat$ORF[peaks],gene_width), 
                      avg_restime=rep(ctss_repeat$avg_restime[peaks],gene_width),
                      group=rep(ctss_repeat$group[peaks],gene_width), 
                      tss=rep(ctss_repeat$tss[peaks],gene_width), 
                      polya=rep(ctss_repeat$polya[peaks],gene_width),
                      relative_pos=c(0:(gene_width-1)),
                      netseq_pos=ctss_repeat$tss[peaks]-c(0:(gene_width-1)),
                      chr=rep(ctss_repeat$chr[peaks],gene_width),
                      ip=rep(ctss_repeat$X237_00_IP[peaks],gene_width))
  }
  negative_bed <- rbind(negative_bed,tmp)
}


################################
##### In the terminal, I am going to use bedtools genomecov to get coverage at each of those positions from the H3 occupancy bedgraph
bed <- bed[!is.na(bed$orf),]
negative_bed <- negative_bed[!is.na(negative_bed$orf),]
netseq_bed <- data.frame(chr=paste("chr", gsub("_cer", "", bed$chr), sep=""), start=bed$netseq_pos, end=bed$netseq_pos+1)
negnet_bed <- data.frame(chr=paste("chr", gsub("_cer", "", negative_bed$chr), sep=""), start=negative_bed$netseq_pos, end=negative_bed$netseq_pos+1)


write.table(netseq_bed, "./other_data/netseq_positive_positions.bed", quote=F, row.names=F, col.names = F, sep="\t")
write.table(negnet_bed, "./other_data/netseq_negative_positions.bed", quote=F, row.names=F, col.names = F, sep="\t")

## bedtools sort -i netseq_positive_positions.bed > netseq_pos_positions_sorted.bed
## bedtools map -a netseq_pos_positions_sorted.bed -b netseq_pos.bedgraph -c 4 -o mean > netseq_pos_coverage.bed

## bedtools sort -i netseq_negative_positions.bed > netseq_neg_positions_sorted.bed
## bedtools map -a netseq_neg_positions_sorted.bed -b netseq_neg.bedgraph -c 4 -o mean > netseq_neg_coverage.bed





#####################################################


pos_values <- read.table("./other_data/netseq_pos_coverage.bed", stringsAsFactors = F)
neg_values <- read.table("./other_data/netseq_neg_coverage.bed", stringsAsFactors = F)

netseq_values <- rbind(pos_values, neg_values)
colnames(netseq_values) <- c("chr", "start", "end", "cov")
netseq_values$chr <- gsub("chr", "", netseq_values$chr)

bed <- rbind(bed, negative_bed)
dt_coverage = data.table(bed, key=c("chr", "netseq_pos"))
dt_netseq_values = data.table(netseq_values, key=c("chr", "start"))

merge_netseq <- merge(dt_coverage, dt_netseq_values, by.x=c("chr", "netseq_pos"), by.y=c("chr", "start"))
merge_netseq <- data.frame(merge_netseq)
merge_netseq <- merge_netseq[!duplicated(merge_netseq[c("chr", "netseq_pos")]),]
merge_netseq$cov[merge_netseq$cov == "."] <- 0
merge_netseq$cov <- as.numeric(merge_netseq$cov)
merge_netseq$group <- factor(merge_netseq$group, levels = c("Shortest", "Short", "Long", "Longest"))


#### Quantiles
merge_netseq <- merge_netseq %>%
  mutate(quantile = ntile(ip, 4))
merge_netseq$ip_group[merge_netseq$quantile==1] <- "Least"
merge_netseq$ip_group[merge_netseq$quantile==2] <- "Less"
merge_netseq$ip_group[merge_netseq$quantile==3] <- "More"
merge_netseq$ip_group[merge_netseq$quantile==4] <- "Most"

merge_netseq$ip_group <- factor(merge_netseq$ip_group, levels = c("Least", "Less", "More", "Most"))


netseq_collapse <- merge_netseq %>% group_by(relative_pos, group) %>% summarise(mean_coverage = mean(cov), sd=sd(cov), n=n()) %>%  mutate(se = sd / sqrt(n),
                                                                                                                                            lower.ci = mean_coverage - qt(1 - (0.05 / 2), n - 1) * se,
                                                                                                                                            upper.ci = mean_coverage + qt(1 - (0.05 / 2), n - 1) * se)
pdf("./plots/netseq_metagene_500bp.pdf")
ggplot(netseq_collapse, aes(x=relative_pos, y=mean_coverage))+
  geom_line(aes(col=group))+
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill=group), alpha = 0.1)+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_fill_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_y_continuous("Net-seq signal (AU)", limits = c(0,100))+
  scale_x_continuous("Position Relative to TSS")
dev.off()


rolling_netseq <- netseq_collapse %>% group_by(group) %>% mutate(rolling_mean = runMean(mean_coverage, 20))

pdf("./plots/netseq_metagene_500bp_rolling.pdf")
ggplot(rolling_netseq, aes(x=relative_pos, y=rolling_mean, col=group))+
  geom_line()+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_y_continuous("Net-seq signal (AU)")+
  scale_x_continuous("Position relative to TSS")
dev.off()


################################################################
netseq_collapse_occupancy <- merge_netseq %>% group_by(relative_pos, ip_group) %>% summarise(mean_coverage = mean(cov), sd=sd(cov), n=n()) %>%  mutate(se = sd / sqrt(n),
                                                                                                                                          lower.ci = mean_coverage - qt(1 - (0.05 / 2), n - 1) * se,
                                                                                                                                          upper.ci = mean_coverage + qt(1 - (0.05 / 2), n - 1) * se)
pdf("./plots/netseq_metagene_500bp_OCC.pdf")
ggplot(netseq_collapse_occupancy, aes(x=relative_pos, y=mean_coverage))+
  geom_line(aes(col=ip_group))+
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill=ip_group), alpha = 0.1)+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_fill_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_y_continuous("Net-seq signal (AU)", limits = c(0,100))+
  scale_x_continuous("Position Relative to TSS")
dev.off()


rolling_netseq_occ <- netseq_collapse_occupancy %>% group_by(ip_group) %>% mutate(rolling_mean = runMean(mean_coverage, 20))

pdf("./plots/netseq_metagene_500bp_rolling_OCC.pdf")
ggplot(rolling_netseq_occ, aes(x=relative_pos, y=rolling_mean, col=ip_group))+
  geom_line()+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_y_continuous("Net-seq signal (AU)")+
  scale_x_continuous("Position relative to TSS")
dev.off()


#################################################################

###### Bed files to count over

##### Length normalized ORF

bed <- data.frame(peaknames=NA, orf=NA, avg_restime=NA, group=NA, tss=NA, polya=NA, relative_pos=NA, netseq_pos=NA, chr=NA, ip=NA)

for(peaks in c(1:nrow(wtss_repeat))){
    gene_width <- wtss_repeat$polya[peaks] - wtss_repeat$tss[peaks] 
    tmp <- data.frame(peaknames=rep(wtss_repeat$peaknames[peaks],gene_width),
                      orf=rep(wtss_repeat$ORF[peaks],gene_width), 
                      avg_restime=rep(wtss_repeat$avg_restime[peaks],gene_width), 
                      group=rep(wtss_repeat$group[peaks],gene_width), 
                      tss=rep(wtss_repeat$tss[peaks],gene_width), 
                      polya=rep(wtss_repeat$polya[peaks],gene_width),
                      relative_pos=c(0:(gene_width-1))/gene_width,
                      netseq_pos=wtss_repeat$tss[peaks]+c(0:(gene_width-1)),
                      chr=rep(wtss_repeat$chr[peaks],gene_width),
                      ip=rep(wtss_repeat$X237_00_IP[peaks],gene_width))
  bed <- rbind(bed,tmp)
}

negative_bed <- data.frame(peaknames=NA, orf=NA, avg_restime=NA, group=NA, tss=NA, polya=NA, relative_pos=NA, netseq_pos=NA, chr=NA, ip=NA)

for(peaks in c(1:nrow(ctss_repeat))){
    gene_width <- ctss_repeat$tss[peaks] - ctss_repeat$polya[peaks] 
    tmp <- data.frame(peaknames=rep(ctss_repeat$peaknames[peaks],gene_width),
                      orf=rep(ctss_repeat$ORF[peaks],gene_width), 
                      avg_restime=rep(ctss_repeat$avg_restime[peaks],gene_width),
                      group=rep(ctss_repeat$group[peaks],gene_width), 
                      tss=rep(ctss_repeat$tss[peaks],gene_width), 
                      polya=rep(ctss_repeat$polya[peaks],gene_width),
                      relative_pos=c(0:(gene_width-1))/gene_width,
                      netseq_pos=ctss_repeat$tss[peaks]-c(0:(gene_width-1)),
                      chr=rep(ctss_repeat$chr[peaks],gene_width),
                      ip=rep(ctss_repeat$X237_00_IP[peaks],gene_width))
  negative_bed <- rbind(negative_bed,tmp)
}


################################
##### In the terminal, I am going to use bedtools genomecov to get coverage at each of those positions from the H3 occupancy bedgraph
bed <- bed[!is.na(bed$orf),]
negative_bed <- negative_bed[!is.na(negative_bed$orf),]
netseq_bed <- data.frame(chr=paste("chr", gsub("_cer", "", bed$chr), sep=""), start=bed$netseq_pos, end=bed$netseq_pos+1)
negnet_bed <- data.frame(chr=paste("chr", gsub("_cer", "", negative_bed$chr), sep=""), start=negative_bed$netseq_pos, end=negative_bed$netseq_pos+1)


write.table(netseq_bed, "./other_data/netseq_positive_positions.bed", quote=F, row.names=F, col.names = F, sep="\t")
write.table(negnet_bed, "./other_data/netseq_negative_positions.bed", quote=F, row.names=F, col.names = F, sep="\t")

## bedtools sort -i netseq_positive_positions.bed > netseq_pos_positions_sorted.bed
## bedtools map -a netseq_pos_positions_sorted.bed -b netseq_pos.bedgraph -c 4 -o mean > netseq_pos_coverage.bed

## bedtools sort -i netseq_negative_positions.bed > netseq_neg_positions_sorted.bed
## bedtools map -a netseq_neg_positions_sorted.bed -b netseq_neg.bedgraph -c 4 -o mean > netseq_neg_coverage.bed





#####################################################


pos_values <- read.table("./other_data/netseq_pos_coverage.bed", stringsAsFactors = F)
neg_values <- read.table("./other_data/netseq_neg_coverage.bed", stringsAsFactors = F)

netseq_values <- rbind(pos_values, neg_values)
colnames(netseq_values) <- c("chr", "start", "end", "cov")
netseq_values$chr <- gsub("chr", "", netseq_values$chr)

bed <- rbind(bed, negative_bed)
dt_coverage = data.table(bed, key=c("chr", "netseq_pos"))
dt_netseq_values = data.table(netseq_values, key=c("chr", "start"))

merge_netseq <- merge(dt_coverage, dt_netseq_values, by.x=c("chr", "netseq_pos"), by.y=c("chr", "start"))
merge_netseq <- data.frame(merge_netseq)
merge_netseq <- merge_netseq[!duplicated(merge_netseq[c("chr", "netseq_pos")]),]
merge_netseq$cov[merge_netseq$cov == "."] <- 0
merge_netseq$cov <- as.numeric(merge_netseq$cov)
merge_netseq$group <- factor(merge_netseq$group, levels = c("Shortest", "Short", "Long", "Longest"))
merge_netseq$relative_pos <- round(merge_netseq$relative_pos,2)



#### Quantiles
merge_netseq <- merge_netseq %>%
  mutate(quantile = ntile(ip, 4))
merge_netseq$ip_group[merge_netseq$quantile==1] <- "Least"
merge_netseq$ip_group[merge_netseq$quantile==2] <- "Less"
merge_netseq$ip_group[merge_netseq$quantile==3] <- "More"
merge_netseq$ip_group[merge_netseq$quantile==4] <- "Most"

merge_netseq$ip_group <- factor(merge_netseq$ip_group, levels = c("Least", "Less", "More", "Most"))
                                                                                                                                     
netseq_collapse <- merge_netseq %>% group_by(relative_pos, group) %>% summarise(mean_coverage = mean(cov), sd=sd(cov), n=n()) %>%  mutate(se = sd / sqrt(n),
                                                                                                                                          lower.ci = mean_coverage - qt(1 - (0.05 / 2), n - 1) * se,
                                                                                                                                          upper.ci = mean_coverage + qt(1 - (0.05 / 2), n - 1) * se)
pdf("./plots/netseq_metagene_entire-transcript.pdf")
ggplot(netseq_collapse, aes(x=relative_pos, y=mean_coverage))+
  geom_line(aes(col=group))+
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill=group), alpha = 0.1)+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_fill_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_y_continuous("Net-seq signal (AU)", limits = c(0,50))+
  scale_x_continuous("Position Relative within transcript (TSS to pA)")
dev.off()

netseq_collapse_occupancy <- merge_netseq %>% group_by(relative_pos, ip_group) %>% summarise(mean_coverage = mean(cov), sd=sd(cov), n=n()) %>%  mutate(se = sd / sqrt(n),
                                                                                                                                          lower.ci = mean_coverage - qt(1 - (0.05 / 2), n - 1) * se,
                                                                                                                                          upper.ci = mean_coverage + qt(1 - (0.05 / 2), n - 1) * se)
pdf("./plots/netseq_metagene_entire-transcript_OCC.pdf")
ggplot(netseq_collapse_occupancy, aes(x=relative_pos, y=mean_coverage))+
  geom_line(aes(col=ip_group))+
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill=ip_group), alpha = 0.1)+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_fill_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_y_continuous("Net-seq signal (AU)", limits = c(0,50))+
  scale_x_continuous("Position Relative within transcript (TSS to pA)")
dev.off()

#################################################################

####################################################



###### Bed files to count over
watson_gtf <- read.csv("watson_gtf.csv")
worf_repeat <- merge(worf_repeat, watson_gtf, by.x="watson", by.y="orf")

crick_gtf <- read.csv("crick_gtf.csv")
corf_repeat <- merge(corf_repeat, crick_gtf, by.x="crick", by.y="orf")

##### 500 bp downstream of ORF
width=500

bed <- data.frame(peaknames=NA, orf=NA, avg_restime=NA, group=NA, start=NA, end=NA, relative_pos=NA, netseq_pos=NA, chr=NA)

for(peaks in c(1:nrow(worf_repeat))){
  if(worf_repeat$end.y[peaks] - worf_repeat$start.y[peaks] > width){
    tmp <- data.frame(peaknames=rep(worf_repeat$peaknames[peaks],width),
                      orf=rep(worf_repeat$watson[peaks],width), 
                      avg_restime=rep(worf_repeat$avg_restime[peaks],width), 
                      group=rep(worf_repeat$group[peaks],width), 
                      start=rep(worf_repeat$start.y[peaks],width), 
                      end=rep(worf_repeat$end.y[peaks],width),
                      relative_pos=c(0:(width-1)),
                      netseq_pos=worf_repeat$start.y[peaks]+c(0:(width-1)),
                      chr=rep(worf_repeat$chr.x[peaks],width))
  }
  else {
    gene_width <- worf_repeat$end.y[peaks] - worf_repeat$start.y[peaks] 
    tmp <- data.frame(peaknames=rep(worf_repeat$peaknames[peaks],gene_width),
                      orf=rep(worf_repeat$watson[peaks],gene_width), 
                      avg_restime=rep(worf_repeat$avg_restime[peaks],gene_width), 
                      group=rep(worf_repeat$group[peaks],gene_width), 
                      start=rep(worf_repeat$start.y[peaks],gene_width), 
                      end=rep(worf_repeat$end.y[peaks],gene_width),
                      relative_pos=c(0:(gene_width-1)),
                      netseq_pos=worf_repeat$start.y[peaks]+c(0:(gene_width-1)),
                      chr=rep(worf_repeat$chr.x[peaks],gene_width))
  }
  bed <- rbind(bed,tmp)
}

negative_bed <- data.frame(peaknames=NA, orf=NA, avg_restime=NA, group=NA, start=NA, end=NA, relative_pos=NA, netseq_pos=NA, chr=NA)

for(peaks in c(1:nrow(corf_repeat))){
  if(corf_repeat$end.y[peaks] - corf_repeat$start.y[peaks] > width){
    tmp <- data.frame(peaknames=rep(corf_repeat$peaknames[peaks],width),
                      orf=rep(corf_repeat$crick[peaks],width), 
                      avg_restime=rep(corf_repeat$avg_restime[peaks],width), 
                      group=rep(corf_repeat$group[peaks],width), 
                      start=rep(corf_repeat$start.y[peaks],width), 
                      end=rep(corf_repeat$end.y[peaks],width),
                      relative_pos=c(0:(width-1)),
                      netseq_pos=corf_repeat$end.y[peaks]-c(0:(width-1)),
                      chr=rep(corf_repeat$chr.x[peaks],width))
  }
  else {
    gene_width <- corf_repeat$end.y[peaks] - corf_repeat$start.y[peaks] 
    tmp <- data.frame(peaknames=rep(corf_repeat$peaknames[peaks],gene_width),
                      orf=rep(corf_repeat$crick[peaks],gene_width), 
                      avg_restime=rep(corf_repeat$avg_restime[peaks],gene_width),
                      group=rep(corf_repeat$group[peaks],gene_width), 
                      start=rep(corf_repeat$start.y[peaks],gene_width), 
                      end=rep(corf_repeat$end.y[peaks],gene_width),
                      relative_pos=c(0:(gene_width-1)),
                      netseq_pos=corf_repeat$end.y[peaks]-c(0:(gene_width-1)),
                      chr=rep(corf_repeat$chr.x[peaks],gene_width))
  }
  negative_bed <- rbind(negative_bed,tmp)
}


################################
##### In the terminal, I am going to use bedtools genomecov to get coverage at each of those positions from the H3 occupancy bedgraph
bed <- bed[!is.na(bed$orf),]
negative_bed <- negative_bed[!is.na(negative_bed$orf),]
netseq_bed <- data.frame(chr=paste("chr", gsub("_cer", "", bed$chr), sep=""), start=bed$netseq_pos, end=bed$netseq_pos+1)
negnet_bed <- data.frame(chr=paste("chr", gsub("_cer", "", negative_bed$chr), sep=""), start=negative_bed$netseq_pos, end=negative_bed$netseq_pos+1)


write.table(netseq_bed, "./other_data/netseq_positive_positions.bed", quote=F, row.names=F, col.names = F, sep="\t")
write.table(negnet_bed, "./other_data/netseq_negative_positions.bed", quote=F, row.names=F, col.names = F, sep="\t")

## bedtools sort -i netseq_positive_positions.bed > netseq_pos_positions_sorted.bed
## bedtools map -a netseq_pos_positions_sorted.bed -b netseq_pos.bedgraph -c 4 -o mean > netseq_pos_coverage.bed

## bedtools sort -i netseq_negative_positions.bed > netseq_neg_positions_sorted.bed
## bedtools map -a netseq_neg_positions_sorted.bed -b netseq_neg.bedgraph -c 4 -o mean > netseq_neg_coverage.bed





#####################################################


pos_values <- read.table("./other_data/netseq_pos_coverage.bed", stringsAsFactors = F)
neg_values <- read.table("./other_data/netseq_neg_coverage.bed", stringsAsFactors = F)

netseq_values <- rbind(pos_values, neg_values)
colnames(netseq_values) <- c("chr", "start", "end", "cov")
netseq_values$chr <- gsub("chr", "", netseq_values$chr)

bed <- rbind(bed, negative_bed)
dt_coverage = data.table(bed, key=c("chr", "netseq_pos"))
dt_netseq_values = data.table(netseq_values, key=c("chr", "start"))

merge_netseq <- merge(dt_coverage, dt_netseq_values, by.x=c("chr", "netseq_pos"), by.y=c("chr", "start"))
merge_netseq <- data.frame(merge_netseq)
merge_netseq <- merge_netseq[!duplicated(merge_netseq[c("chr", "netseq_pos")]),]
merge_netseq$cov[merge_netseq$cov == "."] <- 0
merge_netseq$cov <- as.numeric(merge_netseq$cov)
merge_netseq$group <- factor(merge_netseq$group, levels = c("Shortest", "Short", "Long", "Longest"))


netseq_collapse <- merge_netseq %>% group_by(relative_pos, group) %>% summarise(mean_coverage = mean(cov), sd=sd(cov), n=n()) %>%  mutate(se = sd / sqrt(n),
                                                                                                                                          lower.ci = mean_coverage - qt(1 - (0.05 / 2), n - 1) * se,
                                                                                                                                          upper.ci = mean_coverage + qt(1 - (0.05 / 2), n - 1) * se)
pdf("./plots/netseq_metagene_500bpORF.pdf")
ggplot(netseq_collapse, aes(x=relative_pos, y=mean_coverage))+
  geom_line(aes(col=group))+
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill=group), alpha = 0.1)+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_fill_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  scale_y_continuous("Net-seq signal (AU)", limits = c(0,100))+
  scale_x_continuous("Position Relative to ORF (bp)")
dev.off()




#################################################################

## TAF1 occupancy

taf1_bed <- sig_fits[,c("peaknames", "chr", "start", "end", "X237_00_IP", "avg_restime", "group", "watson", "crick")]


## I am going to manually fix chrIII, which has been adjusted in the genome the data was aligned to
taf1_bed$new_start <- taf1_bed$start
taf1_bed$new_end <- taf1_bed$end

taf1_bed$new_start[taf1_bed$chr=="III" & taf1_bed$start > 201212] <- taf1_bed$start[taf1_bed$chr=="III" & taf1_bed$start > 201212] + 69
taf1_bed$new_end[taf1_bed$chr=="III" & taf1_bed$start > 201212] <- taf1_bed$end[taf1_bed$chr=="III" & taf1_bed$start > 201212] + 69

taf1_bed$new_start[taf1_bed$chr=="III" & taf1_bed$start < 201212 & taf1_bed$start > 13020] <- taf1_bed$start[taf1_bed$chr=="III" & taf1_bed$start < 201212 & taf1_bed$start > 13020] + 397
taf1_bed$new_end[taf1_bed$chr=="III" & taf1_bed$start < 201212 & taf1_bed$start > 13020] <- taf1_bed$end[taf1_bed$chr=="III" & taf1_bed$start < 201212 & taf1_bed$start > 13020] + 397


taf1 <- data.frame(chr=taf1_bed$chr, start=taf1_bed$new_start, end=taf1_bed$new_end)

write.table(taf1, "./other_data/taf1_peakpositions.bed", quote=F, row.names=F, col.names = F, sep="\t")

### ON BOWIE IN ~/EXP_046/rmdups/taf1_data/bamfiles/
#### bedtools sort -i taf1_peakpositions.bed > taf1_peakpositions_sorted.bed
#### bedtools multicov -bed taf1_peakpositions_sorted.bed -bams 438_IP-sorted.bam 439_IP-sorted.bam  > taf1_ip_coverage.bed


####### Import taf1 count data 
taf1_cov <- read.table("./other_data/taf1_ip_coverage.bed", stringsAsFactors = F)
taf1_input <- read.table("./other_data/taf1_input_coverage.bed", stringsAsFactors = F)

colnames(taf1_cov) <- c("chr", "new_start", "new_end", "rep1", "rep2")
taf1_cov$cov <- (taf1_cov$rep1 + taf1_cov$rep2)/2
colnames(taf1_input) <- c("chr", "new_start", "new_end", "rep1", "rep2")
taf1_input$cov <- (taf1_input$rep1 + taf1_input$rep2)/2

taf1_merge <- merge(taf1_cov, taf1_bed, by=c("chr", "new_start"))

## Plot TAF1 counts over Rap1 peaks versus res_time
pdf("./plots/TAF1_Restime_scatter.pdf")
ggplot(taf1_merge, aes(x=avg_restime, y=cov, col=group))+
  geom_point()+
  scale_y_log10(name="Taf1 Occupancy (reads)")+
  xlab("Residence Time (min)")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  annotate("text", x=30, y=10, label="rho = 0.43")+
  annotate("text", x=30, y=30, label="p = 2.2e-16")
dev.off()

cor.test(x=taf1_merge$avg_restime, y=taf1_merge$cov, method="spearman")

pdf("./plots/TAF1_Restime_boxplot.pdf")
ggplot(taf1_merge, aes(x=group, y=cov, col=group))+
  geom_boxplot()+
  scale_y_log10(name="Taf1 Occupancy")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))
dev.off()

one.way <- aov(cov ~ group, data = taf1_merge)
TukeyHSD(one.way)
summary(one.way)


################ 
### plot Taf1 data versus Occupancy

pdf("./plots/TAF1_Rap1Occ_scatter.pdf")
ggplot(taf1_merge, aes(x=X237_00_IP, y=cov, col=""))+
  geom_point()+
  scale_y_log10(name="Taf1 Occupancy")+
  scale_x_log10("Rap1 Occupancy (time=0)")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  annotate("text", x=1000, y=10, label="rho = 0.38")+
  annotate("text", x=1000, y=30, label="p < 1.5e-15")
dev.off()
cor.test(x=taf1_merge$X237_00_IP, y=taf1_merge$cov, method="spearman")


taf1_input <- merge(taf1_input, taf1_bed, by=c("chr", "new_start"))

## Plot input TAF1 data versus restime 

pdf("./plots/TAF1_input_Restime_scatter.pdf")
ggplot(taf1_input, aes(x=avg_restime, y=cov, col=group))+
  geom_point()+
  scale_y_log10(name="Taf1 Occupancy")+
  xlab("Residence Time (min)")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))
dev.off()

cor.test(x=taf1_input$avg_restime, y=taf1_input$cov, method="spearman")




#################################################

#### Sumo data

#### ON BOWIE IN ~/EXP_061/rmdup_bam/
#### bedtools multicov -bed rap1_peaks_sorted.bed -bams 258-IP-sorted.bam 259-IP-sorted.bam  > sumo_ip_rap1peak_coverage.bed

#### Sumo counts in Rap1 peaks 

####### Import sumo count data 
sumo_cov <- read.table("./other_data/sumo_ip_rap1peak_coverage.bed", stringsAsFactors = F)

colnames(sumo_cov) <- c("chr", "new_start", "new_end", "rep1", "rep2")
sumo_cov$cov <- (sumo_cov$rep1 + sumo_cov$rep2)/2


sumo_merge <- merge(sumo_cov, taf1_bed, by=c("chr", "new_start"))

## Plot TAF1 counts over Rap1 peaks versus res_time
pdf("./plots/SUMO_Restime_scatter.pdf")
ggplot(sumo_merge, aes(x=avg_restime, y=cov, col=group))+
  geom_point()+
  scale_y_log10(name="Sumo Occupancy (reads)")+
  xlab("Residence Time (min)")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  annotate("text", x=30, y=20, label="r = 0.11")+
  annotate("text", x=30, y=30, label="p = 0.021")
dev.off()

cor.test(x=sumo_merge$avg_restime, y=sumo_merge$cov)

pdf("./plots/SUMO_Restime_boxplot.pdf")
ggplot(sumo_merge, aes(x=group, y=cov, col=group))+
  geom_boxplot()+
  scale_y_log10(name="SUMO Occupancy")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))
dev.off()



#### Sumo peak intersect with Rap1 peak


#### ON BOWIE IN ~/EXP_061/rmdup_bam/
#### bedtools intersect -a rap1_peaks_sorted.bed -b sumo_peak_regions.bed -wa > intersect_rap1_sumo.bed
sumo_peak <- read.table("./other_data/intersect_rap1_sumo.bed", stringsAsFactors = F)

sumo_peak$peak <- TRUE

colnames(sumo_peak) <- c("chr", "new_start", "new_end", "peak")

sumo_merge <- merge(sumo_peak, sumo_merge, by=c("chr", "new_start"), all.y=TRUE)
sumo_merge$peak[is.na(sumo_merge$peak)] <- FALSE

sumo_only <- sumo_merge[sumo_merge$peak == TRUE,]

pdf("./plots/SUMO_peaks_Restime_scatter.pdf")
ggplot(sumo_only, aes(x=avg_restime, y=cov, col=group))+
  geom_point()+
  scale_y_log10(name="Sumo Occupancy (reads)")+
  xlab("Residence Time (min)")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  annotate("text", x=30, y=20, label="rho = 0.10")+
  annotate("text", x=30, y=30, label="p = 0.19")
dev.off()


cor.test(x=sumo_only$avg_restime, y=sumo_only$cov, method="spearman")


pdf("./plots/SUMO_peak_Restime_boxplot.pdf")
ggplot(sumo_only, aes(x=group, y=cov, col=group))+
  geom_boxplot()+
  scale_y_log10(name="SUMO Occupancy")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))
dev.off()



pdf("./plots/SUMO_peaks_Rap1OCC_scatter.pdf")
ggplot(sumo_only, aes(x=X237_00_IP, y=cov, col=""))+
  geom_point()+
  scale_y_log10(name="Sumo Occupancy (reads)")+
  scale_x_log10(name="Rap1 Occupancy (min)")+
  theme_classic()+
  scale_color_manual(values=c("#666699", "#ADADF2", "#99CCCC", "#006666"))+
  annotate("text", x=1000, y=100, label="rho = 0.53")+
  annotate("text", x=1000, y=150, label="p = 3.2e-14")
dev.off()


cor.test(x=sumo_only$X237_00_IP, y=sumo_only$cov, method="spearman")


