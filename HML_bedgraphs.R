#in this case I am importing just HML from bedgraphs that are not normalized, and will use para_counts as a normalization factor

library(dplyr)
library(data.table)
library(ggplot2)


## Import bedgraphs - ONLY CHROMOSOME III-10000:16000
setwd("~/ERB/EXP_046/HML_bedgraphs")
options(stringsAsFactors=FALSE)
bedgraphDir <- "~/ERB/EXP_046/HML_bedgraphs/" #Tells where to find the bedgraph files. 
samples <- sub("_HML.bedgraph", "", list.files(bedgraphDir, pattern="*_HML.bedgraph$"))


#Double check that your .bedgraph files were read correctly. 
samples 


#Create a dataframe from the bedgraph files 
## Be prepared. This is a doooozy 
#Skips first row - make sure your files look like this
bedgraphList <- lapply(samples, function(s) {
  fread(paste(bedgraphDir, s, "_HML.bedgraph" , sep=""), header=FALSE, skip=1, col.names = c("chr", "start", "end", "cov"))
})
names(bedgraphList) <- samples

big_data <- bind_rows(bedgraphList, .id='id')

#sample norm factor <- input real csv here #### IN THIS CASE USING PARA READCOUNTS
norm <- read.csv("~/ERB/EXP_046/para_sums.csv", header = FALSE)
colnames(norm) <- c("sample","coverage")
norm <- norm %>% select(sample, coverage)
norm2 <- data.frame(t(norm[-1]))
colnames(norm2) <- norm[, 1]
norm2 <- mutate_all(norm2, function(x) as.numeric(as.character(x)))
rownames(norm2) <- c("norm_value")
df <- c(norm$sample)
norm2 <- rbind(norm2,df)
row.names(norm2) <- c("norm_value","sample")
norm_factor <- data.frame(id=norm$sample, normalize = norm$coverage)

dt_bigdata = data.table(big_data, key="id")
dt_normfactor = data.table(norm_factor, key="id")

#this is also a big call 
normalized_bd <- merge(dt_bigdata, dt_normfactor, all.x=TRUE)

normalized_bd$adj_cov <- normalized_bd$cov / normalized_bd$normalize

setwd("~/ERB/EXP_046/final_fits/")
fwrite(normalized_bd, "merged_normalized_HMLbedgraph.txt")


normalized_bd <- fread("merged_normalized_HMLbedgraph.txt")

############################# START HERE

## build list of which samples you would like to plot
input237 <- c("237_00_IN", "237_05_IN", "237_10_IN", "237_15_IN", "237_20_IN", "237_30_IN", "237_45_IN","237_60_IN")
ip237 <- c("237_00_IP", "237_05_IP", "237_10_IP", "237_15_IP", "237_20_IP", "237_30_IP", "237_45_IP","237_60_IP")


input238 <- c("238_00_IN", "238_05_IN", "238_10_IN", "238_15_IN", "238_20_IN", "238_30_IN", "238_45_IN","238_60_IN")
ip238 <- c("238_00_IP", "238_05_IP", "238_10_IP", "238_15_IP", "238_20_IP", "238_30_IP", "238_45_IP","238_60_IP")

input361 <- c("361_00_IN", "361_05_IN", "361_10_IN", "361_15_IN", "361_20_IN", "361_30_IN", "361_45_IN","361_60_IN")
ip361 <- c("361_00_IP", "361_05_IP", "361_10_IP", "361_15_IP", "361_20_IP", "361_30_IP", "361_45_IP","361_60_IP")

input261 <- c("261_00_IN", "261_05_IN", "261_10_IN", "261_15_IN", "261_20_IN", "261_30_IN", "261_45_IN","261_60_IN")
ip261 <- c("261_00_IP", "261_05_IP", "261_10_IP", "261_15_IP", "261_20_IP", "261_30_IP", "261_45_IP","261_60_IP")

input262 <- c("262_00_IN", "262_05_IN", "262_10_IN", "262_15_IN", "262_20_IN", "262_30_IN", "262_45_IN","262_60_IN")
ip262 <- c("262_00_IP", "262_05_IP", "262_10_IP", "262_15_IP", "262_20_IP", "262_30_IP", "262_45_IP","262_60_IP")

#subset data frames to include only those samples
df_input237 <- normalized_bd[(normalized_bd$id %in% input237),]
df_ip237 <- normalized_bd[(normalized_bd$id %in% ip237),]

df_input238 <- normalized_bd[(normalized_bd$id %in% input238),]
df_ip238 <- normalized_bd[(normalized_bd$id %in% ip238),]

df_input361 <- normalized_bd[(normalized_bd$id %in% input361),]
df_ip361 <- normalized_bd[(normalized_bd$id %in% ip361),]

df_input261 <- normalized_bd[(normalized_bd$id %in% input261),]
df_ip261 <- normalized_bd[(normalized_bd$id %in% ip261),]

df_input262 <- normalized_bd[(normalized_bd$id %in% input262),]
df_ip262 <- normalized_bd[(normalized_bd$id %in% ip262),]


df_ip237$time <- c(rep("00",6000),rep("05",6000),rep("10",6000),rep("15",6000),rep("20",6000),rep("30",6000),rep("45",6000),rep("90",6000))
df_ip238$time <- c(rep("00",6000),rep("05",6000),rep("10",6000),rep("15",6000),rep("20",6000),rep("30",6000),rep("45",6000),rep("90",6000))
df_sir_wt_IP <- merge(df_ip237,df_ip238, by.x = c('time','start'), by.y=c('time','start'), all.x = TRUE)
df_sir_wt_IP$avg_cov <- (df_sir_wt_IP$adj_cov.x + df_sir_wt_IP$adj_cov.y)/2

df_ip261$time <- c(rep("00",6000),rep("05",6000),rep("10",6000),rep("15",6000),rep("20",6000),rep("30",6000),rep("45",6000),rep("90",6000))
df_ip262$time <- c(rep("00",6000),rep("05",6000),rep("10",6000),rep("15",6000),rep("20",6000),rep("30",6000),rep("45",6000),rep("90",6000))
df_sir4del_IP <- merge(df_ip261, df_ip262, by.x = c('time','start'), by.y=c('time','start'), all.x = TRUE)
df_sir4del_IP$avg_cov <- (df_sir4del_IP$adj_cov.x + df_sir4del_IP$adj_cov.y) / 2

#Edit roi function to adjust what the plot looks like
roi <- function(chromosome, region, dat){
  new_dat <- dat[dat$chr == chromosome & between(dat$start, region[1], region[2]),]
  ggplot(new_dat, aes(x= start, y=cov, col=id))+
    geom_line()+
    theme_classic()
}



avg_roi <- function(chromosome, region, dat,dat2){
  new_dat <- dat[dat$chr.x == chromosome & between(dat$start, region[1], region[2]),]
  new_dat2 <- dat2[dat2$chr.x == chromosome & between(dat2$start, region[1], region[2]),]
  ggplot(new_dat, aes(x= start, y=avg_cov, col=time))+
    geom_line()+
    theme_classic()
}


#Edit roi function to adjust what the plot looks like
wrap_avg_roi <- function(chromosome, region, dat){
  new_dat <- dat[dat$chr.x == chromosome & between(dat$start, region[1], region[2]),]
  ggplot(new_dat, aes(x= start, y= avg_cov))+
    geom_line(size = 1.5, color="#006666")+
    theme_classic()+
    #geom_segment(aes(x=11147,xend=11348,y=-2,yend=-2), size=1.5, col="black")+
    #geom_segment(aes(x=14950,xend=15151,y=-2,yend=-2), size=1.5, col="black")+
    #geom_segment(aes(x=13440,xend=13641,y=-2,yend=-2), size=1.5, col="black")+
    #geom_segment(aes(x=12378,xend=13413,y=-4, yend=-4), size=2, col = "#2E3192")+
    #geom_segment(aes(x=13670,xend=14517,y=-4, yend=-4), size=2, col = "#2E3192")+
    facet_wrap(~time, ncol=1)+
    theme(strip.background = element_blank(),strip.text.x = element_blank())+
    theme(axis.text.y = element_blank())
}



#Run the function of your region of interest with the subsetted df 
roi(chromosome="III_cer", region=c(10000,16000), dat=df_input237)
roi(chromosome="III_cer", region=c(10000,16000), dat=df_ip237)
roi(chromosome="III_cer", region=c(10000,16000), dat=df_ip238)
roi(chromosome="III_cer", region=c(10000,16000), dat=df_ip361)

roi(chromosome="III_cer", region=c(10000,16000), dat=df_ip261)
roi(chromosome="III_cer", region=c(10000,16000), dat=df_ip262)

avg_roi(chromosome="III_cer", region=c(10000,16000), dat=df_sir_wt_IP)



roi(chromosome="III_cer", region=c(199000,201000), dat=df_ip237)

wrap_avg_roi(chromosome="III_cer", region=c(10000,16000), dat=df_sir_wt_IP)
wrap_avg_roi(chromosome="III_cer", region= c(10000,16000), dat=df_sir4del_IP)

