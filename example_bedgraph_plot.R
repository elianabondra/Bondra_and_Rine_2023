library(ggplot2)
library(dplyr)

setwd("/home/ebondra/EXP_041/part2/chromIII/")
getwd()
#read in tables as dataframes, skip first 3 lines of metadata 

in246 <- read.table("246_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"), stringsAsFactors = FALSE)
in246$sample <- "in246"
ip246 <- read.table("246_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip246$sample <- "ip246"
in247 <- read.table("247_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip247 <- read.table("247_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))

in334 <- read.table("334_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip334 <- read.table("334_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
in335 <- read.table("335_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip335 <- read.table("335_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))

in337 <- read.table("337_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip337 <- read.table("337_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
in338 <- read.table("338_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip338 <- read.table("338_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))

in340 <- read.table("340_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip340 <- read.table("340_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
in341 <- read.table("341_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip341 <- read.table("341_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))

in343 <- read.table("343_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip343 <- read.table("343_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
in344 <- read.table("344_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip344 <- read.table("344_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))

in346 <- read.table("346_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip346 <- read.table("346_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
in347 <- read.table("347_IN_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))
ip347 <- read.table("347_IP_III.bedgraph", skip=3, col.names = c("chr", "start", "end", "cov"))

ipV5notag <- read.table("/home/ebondra/manuscript/supp_notag_V5/bedgraph/10790_IP_med.bedgraph", col.names = c("chr", "start", "end", "cov"), stringsAsFactors = FALSE)
ipV5aa <- read.table("/home/ebondra/manuscript/supp_notag_aa/12923_med.bedgraph", col.names = c("chr", "start", "end", "cov"), stringsAsFactors = FALSE)


in246$sample <- "in246"
in247$sample <- "in247"
ip246$sample <- "ip246"
ip247$sample <- "ip247"
in334$sample <- "in334"
in335$sample <- "in335"
ip334$sample <- "ip334"
ip335$sample <- "ip335"

in337$sample <- "in337"
in338$sample <- "in338"
ip337$sample <- "ip337"
ip338$sample <- "ip338"

in340$sample <- "in340"
ip340$sample <- "ip340"
in341$sample <- "in341"
ip341$sample <- "ip341"

in343$sample <- "in343"
ip343$sample <- "ip343"
in344$sample <- "in344"
ip344$sample <- "ip344"

in346$sample <- "in346"
ip346$sample <- "ip346"
in347$sample <- "in347"
ip347$sample <- "ip347"

ipV5notag$sample <- "notag"
ipV5aa$sample <- "notagaa"

notagV5 <- averagecov(ipV5aa, ipV5notag, chromosome="III", region=c(10000,16000), "notag")


# averaging duplicates and making them into a new dataframe. The idea is to create a new df using the columns
#from the samples above. But then making a new "coverage" column that is (a+b)/2
#merge by chromosome and start position instead because they are different lengths

averagecov <- function(sample.a, sample.b, chromosome, region, samplename){
  newdat.a <- sample.a[sample.a$chr == chromosome & between(sample.a$start, region[1], region[2]),]
  newdat.b <- sample.b[sample.b$chr == chromosome & between(sample.b$start, region[1], region[2]),]
  newdat <- merge(newdat.a, newdat.b, by =c("chr", "start"), sort = FALSE)
  newdat$cov <- ((newdat$cov.x+newdat$cov.y)/2)
  newdat$sample <- samplename 
  keeps <- c("chr", "start", "cov.x", "cov.y", "cov", "sample")
  newdat <- newdat[,keeps]
  return(newdat)
}

SIRwtIP <- averagecov(ip246, ip247, chromosome="III", region=c(0,220000), "SIRwtIP")
SIRwtIN <- averagecov(in246, in247, chromosome="III", region=c(0,220000), "SIRwtIN")

SIRwt <- rbind(SIRwtIP, SIRwtIN, notagV5)
#SIRwt_ind <-rbind(ip246,ip247,in247,in246)

SIRmutIP <- averagecov(ip334, ip335, chromosome="III", region=c(10000,220000), "SIRmutIP")
SIRmutIN <- averagecov(in334, in335, chromosome="III", region=c(10000,220000), "SIRmutIN")

SIRmut <- rbind(SIRmutIP, SIRmutIN)
#SIRmut_ind <-rbind(ip334,ip335,in334,in335)

SIRtotal <- rbind(SIRwt, SIRmut, notagV5)
SIRtotal_ind <-rbind(SIRwt_ind,SIRmut_ind)

sir4delwtIP <- averagecov(ip340, ip341, chromosome="III", region=c(10000,220000), "sir4delwtIP")
sir4delwtIN <- averagecov(in340, in341, chromosome="III", region=c(10000,220000), "sir4delwtIN")
sir4delwt <- rbind(sir4delwtIP, sir4delwtIN)
#sir4delwt_ind <-rbind(ip340,ip341,in340,in341)

sir4delmutIP <- averagecov(ip343, ip344, chromosome="III", region=c(10000,220000), "sir4delmutIP")
sir4delmutIN <- averagecov(in343, in344, chromosome="III", region=c(10000,220000), "sir4delmutIN")
sir4delmut <- rbind(sir4delmutIP, sir4delmutIN)
#sir4delmut_ind <- rbind(ip343,ip344,in343,in344)

sir4deltotal <- rbind(sir4delmut,sir4delwt)


#### FINAL FIG 2A
roiA <- function(chromosome, region, dat){
  new_dat <- dat[dat$chr == chromosome & between(dat$start, region[1], region[2]),]
  ggplot(new_dat, aes(x= start, y=cov, col=sample))+
    geom_line(size=1.5)+
    scale_color_manual(values=c("#58595B","#CCCCCC", "#99CCCC","#CCCCCC", "#006666"))+
    theme_classic()+
    scale_y_continuous(limits= c(-10,75), name="Normalized Coverage")+
    xlab("Position")+
    geom_segment(aes(x=11147,xend=11348,y=-3,yend=-3), size=1.5, col="black")+
    geom_segment(aes(x=14950,xend=15151,y=-3,yend=-3), size=1.5, col="black")+
    geom_segment(aes(x=13440,xend=13641,y=-3,yend=-3), size=1.5, col="black")+
    geom_segment(aes(x=12378,xend=13405,y=-5, yend=-5), size=2, col = "#2E3192")+
    geom_segment(aes(x=13675,xend=14517,y=-5, yend=-5), size=2, col = "#2E3192")+
    theme(legend.position = "none")
  # add things here to make your ggplot prettier
}

roiA(chromosome="III", region=c(10000,16000), dat=SIRtotal)