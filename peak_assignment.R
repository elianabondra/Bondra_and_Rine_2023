library(stringr)
library(data.table)
library(dplyr)

#### Assignment of peaks to ORFs 

all_peaks <- read.csv("./subset_peaks/peaknames.csv", stringsAsFactors = F)
colnames(all_peaks) <- c("peakname")
all_peaks$species[str_detect(all_peaks$peakname, "cer")] <- "cer"
all_peaks$species[str_detect(all_peaks$peakname, "para")] <- "para"


## get peak summits and chr from peaknames
all_peaks$start <- gsub("(?:[^.]+\\.){1}([^.]+).*", "\\1", all_peaks$peakname)
all_peaks$end <- gsub("(?:[^.]+\\.){2}([^.]+).*", "\\1", all_peaks$peakname)
all_peaks$summit <- (as.numeric(all_peaks$start) + as.numeric(all_peaks$end)-1)/2
all_peaks$chr <- gsub("\\_cer.*", "", all_peaks$peakname)


## I am going to manually fix chrIII, which has been adjusted in the genome the data was aligned to
all_peaks$summit[all_peaks$chr=="III" & all_peaks$summit > 296146] <- all_peaks$summit[all_peaks$chr=="III" & all_peaks$summit > 296146] - 109


## We define peaks proximal to genes using 300 bp upstream of ORF start
gtf <- read.table("~/Genomes/Saccharomyces_cerevisiae.R64-1-1.85.gtf", skip=5, sep = '\t')
gtf <- gtf[,c(1,3,4,5,7,9)]
colnames(gtf) <- c("chr", "type", "start", "end", "strand", "name")
gtf <- gtf[gtf$type=="gene",]
gtf$orf <- gsub("\\;.*", "", gtf$name)
gtf$orf <- gsub("gene_id ", "", gtf$orf)
trna <- gtf[str_detect(gtf$orf, "t"),]
gtf <- gtf[str_detect(gtf$orf, "Y"),]
gtf$min_pos[gtf$strand=='+'] <- gtf$start[gtf$strand=='+'] - 500
gtf$min_pos[gtf$strand=='-'] <- gtf$end[gtf$strand=='-']

gtf$max_pos[gtf$strand=='+'] <- gtf$start[gtf$strand=='+']
gtf$max_pos[gtf$strand=='-'] <- gtf$end[gtf$strand=='-'] + 500

watson_gtf <- gtf[gtf$strand == "+",]
write.csv(watson_gtf, "watson_gtf.csv")
crick_gtf <- gtf[gtf$strand == "-",]
write.csv(crick_gtf, "crick_gtf.csv")

merge_orfs <- setDT(all_peaks)[setDT(watson_gtf), on = c("chr", "summit>=min_pos", "summit<=max_pos"), watson:= orf]
merge_orfs <- setDT(merge_orfs)[setDT(crick_gtf), on = c("chr", "summit>=min_pos", "summit<=max_pos"), crick:= orf]
merge_orfs <- setDT(merge_orfs)[setDT(trna), on = c("chr", "summit>=start", "summit<=end"), trna:= orf]

sum(!is.na(merge_orfs$watson) | !is.na(merge_orfs$crick))

merge_orfs$rDNA <- ifelse(merge_orfs$chr == 'XII'& merge_orfs$summit < 491000 & merge_orfs$summit > 450000, TRUE, FALSE)

## From Challal 2018. molec cel
rap1_targets2 <- read.csv("./other_data/rap1_targets_Challal.csv", col.names = 0)

sum(merge_orfs$watson %in% rap1_targets2$X0 | merge_orfs$crick %in% rap1_targets2$X0)


## These are the exact same matches

yeast_sizes <- read.table("./yeast-sizes.txt")
yeast_sizes$species[str_detect(yeast_sizes$V1, "cer")] <- "cer"
yeast_sizes <- yeast_sizes[yeast_sizes$species == "cer",]
yeast_sizes$chr <- gsub("\\_cer.*", "", yeast_sizes$V1)


merge_orfs <- yeast_sizes %>% select(c("V2", "chr")) %>% merge(merge_orfs, by="chr")

merge_orfs$subtelo <- FALSE
merge_orfs$subtelo[merge_orfs$summit < 15000 | (merge_orfs$summit > (merge_orfs$V2-15000))] <- TRUE

## I am going to move forward for now with these peaks. I am pretty confident that these are better assignments than we had before. 
write.csv(merge_orfs, "./subset_peaks/peak_assignment.csv")





####################################################
#
## Import major transcript isoform coordinates from Pelechano et al. 2013
## This file has TSS and polyA sites for 4928 ORFs 
anno_file <- read.csv("../other_data/Pelechano_mTIF_anno_file.csv")
anno_file$chr <- as.roman(anno_file$chr)
anno_file$chr <- as.character(anno_file$chr)
anno_file <- na.omit(anno_file)

## Set boundaries for Rap1 summit distance from TSS, i.e. if Rap1 summit is in the X bp upstream of TSS those should be matched
## Lets try 500
anno_file$min_tss[anno_file$strand=='+'] <- anno_file$tss[anno_file$strand=='+'] - 500
anno_file$min_tss[anno_file$strand=='-'] <- anno_file$tss[anno_file$strand=='-']

anno_file$max_tss[anno_file$strand=='+'] <- anno_file$tss[anno_file$strand=='+']
anno_file$max_tss[anno_file$strand=='-'] <- anno_file$tss[anno_file$strand=='-'] + 500

watson_anno_file <- anno_file[anno_file$strand=='+',]
crick_anno_file <- anno_file[anno_file$strand=='-',]


merge_peaks <- setDT(all_peaks)[setDT(watson_anno_file), on = c("chr", "summit>=min_tss", "summit<=max_tss"), watson := orf, mult="all"]
merge_peaks <- setDT(merge_peaks)[setDT(crick_anno_file), on = c("chr", "summit>=min_tss", "summit<=max_tss"), crick := orf, mult="all"]

