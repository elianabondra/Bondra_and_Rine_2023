#Here I am going to try to normalize with just paradoxus read counts, because basd on "normalized_11042022", normalizing is causing 
#a greater difference in correlation across samples, so I'm going to see if I can better maintain the correlation over time 
#by using a different normalization factor. I will see whether the correlation of countTable 00 v 00, etc. looks more or less correlated
#than the norm.counttable
#usingALL featurecounts now, with rmdups

library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

## Import count data from featureCounts

options(stringsAsFactors=FALSE)
countDir <- "~/ERB/EXP_046/2023-02-16_counts/" #Tells where to find the count files. 
samples <- sub("-counts.txt$", "", list.files(countDir, pattern="-counts.txt$"))

#Double check that your _count.txt files were read correctly. 
samples 

#Create a dataframe from the text files from fpcount program. 
countList <- lapply(samples, function(s) {
  read.delim(sprintf("%s/%s-counts.txt", countDir, s), header=FALSE, skip =2, row.names=1, na.strings=c("N/A"))
})
names(countList) <- samples

countTable <- do.call(cbind.data.frame, lapply(countList, function(qe) { qe[[6]] }))
colnames(countTable) <- samples 

setwd("~/ERB/PDD-2023-02-07_EXP-046")

#Make sure they were all read in as numbers
countTable <- mutate_all(countTable, function(x) as.numeric(as.character(x)))
row.names(countTable) <- row.names(countList[[1]])
ggplot(countTable, aes (x = `238_30_IP`, y = `237_30_IP`))+
  geom_point()+
  geom_abline(slope=1)

## read in normalization table and mutate it to be able to combine with count table 
norm <- read.csv("~/ERB/EXP_046/para_sums.csv", header = FALSE)
colnames(norm) <- c("sample","coverage")
norm <- norm %>% select(sample, coverage)
norm$coverage <- norm$coverage/1000000

#Transpose and make first row column names
norm2 <- data.frame(t(norm[-1]))
colnames(norm2) <- norm[, 1]
norm2 <- mutate_all(norm2, function(x) as.numeric(as.character(x)))
rownames(norm2) <- c("norm_value")

norm.countTable <- countTable
#norm.countTable <- subset (norm.countTable, select = -c(1,2)) #remove 12923 samples because they are missing from the norm2 table and then I can't merge


norm.countTable <- bind_rows(norm.countTable, norm2)

norm2.countTable <- data.frame(lapply(norm.countTable, function(X) X/X[nrow(norm.countTable)]))
rownames(norm2.countTable) <- rownames(norm.countTable) #renaming rows to have the full informational name with position etc
colnames(norm2.countTable) <- colnames(norm.countTable) #renaming columns to not have an X in front of them


write.table(norm2.countTable, sep='\t', "./normalized_counts/Normalized_counts.txt")
subset_norm2 <- read.table("./normalized_counts/Normalized_counts.txt", sep='\t', check.names = F)



######################################################################

input_norm <- subset_norm2[,str_detect(colnames(subset_norm2), "IN")]
ip_norm <- subset_norm2[,str_detect(colnames(subset_norm2), "IP")]

ip_237 <- ip_norm[,str_detect(colnames(ip_norm), "237")]
write.csv(row.names(ip_237), "./subset_peaks/peaknames.csv")

write.csv(ip_237,"./normalized_counts/ip.norm_237.csv")
ip_238 <- ip_norm[,str_detect(colnames(ip_norm), "238")]
write.csv(ip_238,"./normalized_counts/ip.norm_238.csv")
ip_261 <- ip_norm[,str_detect(colnames(ip_norm), "261")]
write.csv(ip_261,"./normalized_counts/ip.norm_261.csv")
ip_262 <- ip_norm[,str_detect(colnames(ip_norm), "262")]
write.csv(ip_262,"./normalized_counts/ip.norm_262.csv")
ip_361 <- ip_norm[,str_detect(colnames(ip_norm), "361")]


##### makes a function that takes in a dataframe and time variable 

time <- c(0,5,10,15,20,30,45,60)


nls_fit <- function(dat, time){
  df_nls <- data.frame(matrix(NA, nrow = nrow(dat), ncol = 13))
  colnames(df_nls)<- c("coef_yf", "coef_y0", "log_koff", "ste_yf", "ste_y0", "ste_log_koff", "tval_yf", "tval_y0", "tval_log_koff", "pval_yf", "pval_y0", "pval_koff", "sigma")
  row.names(df_nls)<- row.names(dat)
  
  
  for (row in 1:nrow(dat)) {
    ip_0 <- as.numeric(dat[row,1])
    ip <- as.numeric(dat[row,])/ip_0
    fit <- tryCatch(nls(ip ~ SSasymp(time, yf, y0, log_koff)), error=function(e) NA)
    if (!is.na(fit)) {
      df_nls[row,] <- c(coef(summary(fit), complete=FALSE), sigma(fit))
    }
    else {
      df_nls[row,] <- NA
    }
  }
  df_nls$peaknames <- row.names(df_nls)
  return(df_nls)
}


nlsfit_237 <- nls_fit(ip_237, time)
write.csv(nlsfit_237, "./final_fits/NLS_Fit_237_final.csv")
nlsfit_238 <- nls_fit(ip_238, time)
write.csv(nlsfit_238, "./final_fits/NLS_Fit_238_final.csv")
nlsfit_261 <- nls_fit(ip_261, time)
write.csv(nlsfit_261, "./final_fits/NLS_Fit_261_final.csv")
nlsfit_262 <- nls_fit(ip_262, time)
write.csv(nlsfit_262, "./final_fits/NLS_Fit_262_final.csv")
nlsfit_361 <- nls_fit(ip_361, time)
write.csv(nlsfit_361, "./final_fits/NLS_Fit_361_final.csv")


## Some statistcs on how well this worked
1-mean(is.na(nlsfit_237))
#0.961 (all peaks)




1-mean(is.na(nlsfit_238))
#0.966
sum(na.omit(nlsfit_238$pval_koff) < 0.05)
#1504



1-mean(is.na(nlsfit_361))
#0.951
sum(na.omit(nlsfit_361$pval_koff) < 0.05)
#1386



###################################################
#now I am going to merge wildtype and sir samples into the same dataframe 


sir_wt <- merge(nlsfit_237, nlsfit_238, by="peaknames") %>% select("peaknames", "log_koff.x", "ste_log_koff.x", "pval_koff.x", "sigma.x", "log_koff.y", "ste_log_koff.y", "pval_koff.y", "sigma.y")
colnames(sir_wt) <- c("peaknames", "log_koff.237", "ste_log_koff.237", "pval_koff.237", "sigma.237", "log_koff.238", "ste_log_koff.238", "pval_koff.238", "sigma.238")

sir_wt <- merge(sir_wt, nlsfit_361, by="peaknames") %>% select("peaknames", "log_koff.237", "ste_log_koff.237", "pval_koff.237", "sigma.237", "log_koff.238", "ste_log_koff.238", "pval_koff.238", "sigma.238", "log_koff", "ste_log_koff", "pval_koff", "sigma")
colnames(sir_wt)[c(10,11,12, 13)] <- c("log_koff.361", "ste_log_koff.361", "pval_koff.361", "sigma.361")

sir_wt$avg_logkoff <- (sir_wt$log_koff.237 + sir_wt$log_koff.238)/2
sir_wt$avg_restime <- 1/exp(sir_wt$avg_logkoff)


write.csv(sir_wt, "./final_fits/NLS_fit_WT.csv")

sir4_del <- merge(nlsfit_261, nlsfit_262, by="peaknames") %>% select("peaknames", "log_koff.x", "ste_log_koff.x", "pval_koff.x", "sigma.x", "log_koff.y", "ste_log_koff.y", "pval_koff.y", "sigma.y")
colnames(sir4_del) <- c("peaknames", "log_koff.261", "ste_log_koff.261", "pval_koff.261", "sigma.261", "log_koff.262", "ste_log_koff.262", "pval_koff.262", "sigma.262")

sir4_del$avg_logkoff <- (sir4_del$log_koff.261 + sir4_del$log_koff.262)/2
sir4_del$avg_restime <- 1/exp(sir4_del$avg_logkoff)

write.csv(sir4_del, "./final_fits/NLS_fit_sir4.csv")



###############################################################################
