####################################################
library(dplyr)
library(reshape2)
library(ggplot2)
library(scales)
library(RColorBrewer)



setwd("~/ERB/PDD-2023-02-07_EXP-046/")

#### I am going to subset data here to include only peaks nearby genes and sub-telometric
#### First count data correlations 

norm_counts <- read.table("./normalized_counts/Normalized_counts.txt", check.names = F)
write.csv(row.names(norm_counts), "./subset_peaks/peaknames.csv", row.names=F)
peak_assignment <- read.csv("./subset_peaks/peak_assignment.csv")

keeps <- peak_assignment$peakname[(!is.na(peak_assignment$watson) | !is.na(peak_assignment$crick) | peak_assignment$subtelo) & !peak_assignment$rDNA & is.na(peak_assignment$trna)]


ip_norm <- norm_counts[row.names(norm_counts) %in% keeps,str_detect(colnames(norm_counts), "IP")]
ip_237 <- ip_norm[,str_detect(colnames(ip_norm), "237")]



###### Heat maps over the time course
matrix <- as.matrix(ip_237)

pdf("./plots/Heatmap_IP_237_subset-counts.pdf")
heatmap(matrix, scale="row", Colv=NA, labRow = FALSE, col= colorRampPalette(brewer.pal(8, "Reds"))(25))
dev.off()

ip_238 <- ip_norm[,str_detect(colnames(ip_norm), "238")]
matrix <- as.matrix(ip_238)

pdf("./plots/Heatmap_IP_238_subset-counts.pdf")
heatmap(matrix, scale="row", Colv=NA, labRow = FALSE, col= colorRampPalette(brewer.pal(8, "Reds"))(25))
dev.off()


###################################################################################

### Correlation heat map based on subset count data
ip_norm_wt <- ip_norm[,(str_detect(colnames(ip_norm), "237")|str_detect(colnames(ip_norm), "238"))]
cor_matrix <- round(cor(ip_norm_wt),2)

get_upper_tri<-function(cor_matrix){
  cor_matrix[lower.tri(cor_matrix)] <- NA
  return(cor_matrix)
}
cor_matrix <- get_upper_tri(cor_matrix)
cor_matrix <- melt(cor_matrix, na.rm = TRUE)

pdf("./plots/Correlation_heatmap.pdf")
ggheatmap <- ggplot(cor_matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "#F9EBEA", high = "#922B21", mid = "#CD6155", 
                       midpoint = 0.9, limit = c(0.8,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
ggheatmap
dev.off()


################################################################
### A few individual replicate plots

pdf("./plots/Time0_replicates.pdf")
ggplot(ip_norm_wt, aes(x=`237_00_IP`, y=`238_00_IP`))+
  geom_point()+
  scale_x_log10(name="IP peak counts t=0 rep_1")+
  scale_y_log10(name="IP peak counts t=0 rep_2")+
  theme_classic()
dev.off()

pdf("./plots/Time30_replicates.pdf")
ggplot(ip_norm_wt, aes(x=`237_30_IP`, y=`238_30_IP`))+
  geom_point()+
  scale_x_log10(name="IP peak counts t=30 rep_1")+
  scale_y_log10(name="IP peak counts t=30 rep_2")+
  theme_classic()
dev.off()



###########################################################################

### Justification for choosing to proceed with the peaks upstream of ORFS and telomere genes
sir_wt_allpeak <- read.csv('./final_fits/NLS_fit_WT.csv', stringsAsFactors = F)

peak_info <- read.csv("./subset_peaks/peak_assignment.csv", stringsAsFactors = F)
anno <- peak_info  %>% merge(sir_wt_allpeak, by.x="peakname", by.y="peaknames")
anno$gene_by <- FALSE
anno$gene_by[!is.na(anno$watson)|!is.na(anno$crick)] <- TRUE
anno <- norm_counts %>% select(c(4)) %>% merge(anno, by.x=0, by.y="peakname")

anno$cat <- "Neither"
anno$cat[anno$gene_by == TRUE] <- "Gene Proximal"
anno$cat[anno$subtelo == TRUE] <- "Sub-telometric"
anno$cat[anno$gene_by == TRUE & anno$subtelo == TRUE] <- "Both"
anno$cat[!is.na(anno$trna)] <- "tRNA"
anno$cat[anno$summit < 500 | anno$V2-anno$summit < 500] <- "Telometric"


#### plot IP t=0 coverage for different peak classifications 

pdf("./plots/Peak_category_coveragehist.pdf")
ggplot(anno, aes(x=`237_00_IP`, fill=cat))+
  geom_histogram(bins=150)+
  scale_x_log10(name="Peak IP signal (time=0) AU")+
  ylab("Number of peaks")+
  theme_classic()+
  scale_fill_manual(values=c("#006666", "#99CCCC", "grey", "#666699", "#ADADF2", "red"))
dev.off()


#### plot "error" in fit for different peak classifications  
pdf("./plots/Peak_category_sigmahist.pdf")
ggplot(anno, aes(x=sigma.237, fill=cat))+
  geom_histogram(bins=150)+
  scale_x_log10(name="Sigma (residual st. dev.)")+
  ylab("Number of peaks")+
  theme_classic()+
  scale_fill_manual(values=c("#006666", "#99CCCC", "grey", "#666699", "#ADADF2", "red"))
dev.off()



###### Pie charts for peak classifications
peak_count <- anno %>% count(cat)

pdf("./plots/Peak_category_pie.pdf")
ggplot(peak_count, aes(x = "", y = n, fill = cat)) +
  geom_col(color = "black") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  scale_fill_manual(values=c("#006666", "#99CCCC", "grey", "#666699", "#ADADF2"))+
  theme_void()
dev.off()



##### Pie chart for peaks that don't fit 
no_fit_peak_count <- anno %>% subset(is.na(avg_restime)) %>% count(cat)

pdf("./plots/Peak_category_pie_nofit.pdf")
ggplot(no_fit_peak_count, aes(x = "", y = n, fill = cat)) +
  geom_col(color = "black") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y")+
  scale_fill_manual(values=c("#006666", "#99CCCC", "grey", "#666699", "#ADADF2"))+
  theme_void()
dev.off()


##########################################################
#Subset Peaks

subset_counts <- norm_counts[row.names(norm_counts) %in% keeps,]
write.table(subset_counts, "./normalized_counts/SUBSET_normalized.counts.txt")
colnames(anno)[1] <- "peaknames"
subset_fits <- anno[anno$peaknames %in% keeps,]

subset_fits <- subset_fits[!is.na(subset_fits$avg_logkoff),]
subset_fits <- subset_fits[subset_fits$pval_koff.237 < 0.05 & subset_fits$pval_koff.238 < 0.05,]

if (!file.exists("SGD_features.tab")) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile="SGD_features.tab")
}
sgd <- read.delim("SGD_features.tab", header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))

sgd <- sgd %>% select(c("name", "gene"))

subset_fits <- subset_fits %>% merge(sgd, by.x="watson", by.y="name", all.x=TRUE) 
colnames(subset_fits)[colnames(subset_fits)=="gene"] <- "watson_genename"

subset_fits <- subset_fits %>% merge(sgd, by.x="crick", by.y="name", all.x=TRUE) 
colnames(subset_fits)[colnames(subset_fits)=="gene"] <- "crick_genename"

write.csv(subset_fits, "./final_fits/SUBSET_NLSfits_SIRWT.csv")

