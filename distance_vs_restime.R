peaks <- read.csv("/mnt/ingolialab/paige_diamond/ERB/EXP_046/final_targets/sig_fits.csv")
library(dplyr)
library(stringr)
library(tidyr)

peaks <- select(sig_fits, c("peaknames", "avg_restime", "chr", "start","end"))
peaks <- peaks[!duplicated(peaks$peaknames),]
#my_data[!duplicated(my_data$Sepal.Width), ]
shortest <- peaks %>% slice(1:112)
longest <- peaks %>% slice(336:448)

shortest <- select(shortest, c("chr","start","end"))
shortest <- shortest %>% distinct()
longest <- select(longest, c("chr","start","end"))
longest <- longest %>% distinct()
write.table(shortest, "/mnt/ingolialab/paige_diamond/ERB/EXP_046/final_targets/shortest.bed", row.names = FALSE, col.names = FALSE)
write.table(longest, "/mnt/ingolialab/paige_diamond/ERB/EXP_046/final_targets/longest.bed", row.names = FALSE, col.names = FALSE)


lengths <- read.table("/mnt/ingolialab/paige_diamond/ERB/EXP_046/final_targets/yeast-sizes.txt", stringsAsFactors = FALSE)
colnames(lengths) <- c("chr", "length")
peaks2 <- merge(peaks,lengths, by = "chr")
peaks2$summit <- (peaks2$start + peaks2$end +1)/2

dat <- c(peaks2$peaknames)
dat <- c(peaks2$summit)
ifelse(peaks2$summit < (peaks2$length/2)) {
  print((peaks2$summit))
} ifelse {
  print("NA")
}
f <- function(t) {
  ifelse(t < 0, 0, (2*t)/((1+t^2)^2)
}
This should work if you run f(t) where t is a numeric vector. The syntax is:
  
  ifelse(condition, do_if_true, do_if_false)
f<- function(peaks2){
  ifelse(peaks2$summit < (peaks2$length/2), print(peaks2$summit), print(peaks2$length-peaks2$summit))
}
peaks2$distance <- f(peaks2)

ggplot(peaks2, aes(x = avg_restime, y = distance))+
  geom_point()
cor.test(x=peaks2$avg_restime, y=peaks2$distance)
