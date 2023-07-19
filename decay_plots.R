results237 <- read.csv("~/ERB/EXP_046/final_fits/ip.norm_237.csv")
results238 <- read.csv("~/ERB/EXP_046/final_fits/ip.norm_238.csv")
results261 <- read.csv("~/ERB/EXP_046/final_fits/ip.norm_261.csv")
results262 <- read.csv("~/ERB/EXP_046/final_fits/ip.norm_262.csv")

colnames(results237) <- c("peaknames", "time00", "time05","time10","time15","time20","time30","time45","time60")
colnames(results238) <- c("peaknames", "time00", "time05","time10","time15","time20","time30","time45","time60")
colnames(results261) <- c("peaknames", "time00", "time05","time10","time15","time20","time30","time45","time60")
colnames(results262) <- c("peaknames", "time00", "time05","time10","time15","time20","time30","time45","time60")
results237$sample <- "IP237"
results238$sample <- "IP238"
results261$sample <- "IP261"
results262$sample <- "IP262"
SIR_results <- rbind(results237,results238)
sir4del_results <- rbind(results261, results262)

full_results <- rbind(sir4del_results,SIR_results)

hml_e="III_cer.11091.11392"
hml_p="III_cer.12992.13293"
hml_i="III_cer.14503.14804"
mat_p="III_cer.200150.200451"



row=hml_p
poi <- full_results[full_results$peaknames==row,]
poi <- melt(poi)
ip237_norm <- poi$value[poi$sample == "IP237" & poi$variable == "time00"] /100
ip238_norm <- poi$value[poi$sample == "IP238" & poi$variable == "time00"] /100
#ip361_norm <- poi$value[poi$sample == "IP361" & poi$variable == "time00"] / 100
poi$norm_value[poi$sample == "IP237"] <- poi$value[poi$sample == "IP237"] / ip237_norm
poi$norm_value[poi$sample == "IP238"] <- poi$value[poi$sample == "IP238"] / ip238_norm

row=mat_p
poi2 <- full_results[full_results$peaknames==row,]
poi2 <- melt(poi2)
ip237_norm2 <- poi2$value[poi2$sample == "IP237" & poi2$variable == "time00"] /100
ip238_norm2 <- poi2$value[poi2$sample == "IP238" & poi2$variable == "time00"] /100
poi2$norm_value[poi2$sample == "IP237"] <- poi2$value[poi2$sample == "IP237"] / ip237_norm2
poi2$norm_value[poi2$sample == "IP238"] <- poi2$value[poi2$sample == "IP238"] / ip238_norm2

poi <- rbind(poi,poi2)


poi$time[poi$variable == "time00"] <- 0
poi$time[poi$variable == "time05"] <- 5
poi$time[poi$variable == "time10"] <- 10
poi$time[poi$variable == "time15"] <- 15
poi$time[poi$variable == "time20"] <- 20
poi$time[poi$variable == "time30"] <- 30
poi$time[poi$variable == "time45"] <- 45
poi$time[poi$variable == "time60"] <- 60

poi$color[poi$peaknames == "III_cer.12992.13293" & poi$sample=="IP237"] <- "#006666"
poi$color[poi$peaknames == "III_cer.12992.13293" & poi$sample=="IP238"] <- "#006666"

poi$color[poi$peaknames == "III_cer.200150.200451" & poi$sample=="IP237"] <- "#99CCCC"
poi$color[poi$peaknames == "III_cer.200150.200451" & poi$sample=="IP238"] <- "#99CCCC"

#row <- hml_i
#poi <- data.frame(time = c(0,5,10, 15, 20, 30, 45, 60), ip= as.numeric(SIR_results[SIR_results$peaknames == row,ip.data])*scalefactor, sample="SIR_results")

ggplot(poi, aes(x = time, y = norm_value, col = color))+
  geom_point(size=2.5)+
  geom_smooth(data=poi[poi$sample=="IP237" & poi$peaknames==hml_p,], aes(), color="#006666", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  geom_smooth(data=poi[poi$sample=="IP238" & poi$peaknames==hml_p,], aes(), color="#006666", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  geom_smooth(data=poi[poi$sample=="IP237" & poi$peaknames=="III_cer.200150.200451",], aes(), color="#99CCCC", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  geom_smooth(data=poi[poi$sample=="IP238"& poi$peaknames=="III_cer.200150.200451",], aes(), color="#99CCCC", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #geom_smooth(data=poi[poi$sample=="IP361",], aes(), color="green", linetype=2, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #geom_smooth(data=poi[poi$sample=="IP261",], aes(), color="#99CCCC", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #geom_smooth(data=poi[poi$sample=="IP262",], aes(), color="#99CCCC", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #stat_smooth(method = "nls", formula=y~SSasymp(x, yf, y0, log_koff),  method.args = list(koff=coef(summary(fit))[3]), se = FALSE)+
  theme_classic()+
  scale_color_manual(values=c("#006666","#99CCCC","#006666","#99CCCC"))+
  #scale_color_manual(values=c("#339999","#339999"))+
  #scale_color_manual(values=c("#99CCCC","#99CCCC"))+
  ylab("Normalized IP counts (AU)")+
  xlab("Time (min)")+
  theme(legend.position = "none")

#hml_p only plot
ggplot(poi, aes(x = time, y = norm_value, col = sample))+
  geom_point(size=2.5)+
  geom_smooth(data=poi[poi$sample=="IP237" & poi$peaknames=="III_cer.12992.13293",], aes(), color="#006666", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  geom_smooth(data=poi[poi$sample=="IP238" & poi$peaknames=="III_cer.12992.13293",], aes(), color="#006666", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #geom_smooth(data=poi[poi$sample=="IP237",], aes(), color="#99CCCC", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #geom_smooth(data=poi[poi$sample=="IP238",], aes(), color="#99CCCC", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #geom_smooth(data=poi[poi$sample=="IP361",], aes(), color="green", linetype=2, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #geom_smooth(data=poi[poi$sample=="IP261",], aes(), color="#99CCCC", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #geom_smooth(data=poi[poi$sample=="IP262",], aes(), color="#99CCCC", linetype=5, size = 1, method='nls', formula=y~SSasymp(x, yf, y0, log_koff), se=F) +
  #stat_smooth(method = "nls", formula=y~SSasymp(x, yf, y0, log_koff),  method.args = list(koff=coef(summary(fit))[3]), se = FALSE)+
  theme_classic()+
  scale_color_manual(values=c("#006666","#006666","#99CCCC","#99CCCC"))+
  #scale_color_manual(values=c("#339999","#339999"))+
  #scale_color_manual(values=c("#666699","#666699"))+
  ylab("Normalized IP counts (AU)")+
  xlab("Time (min)")+
  theme(legend.position = "none")

nls(formula=norm_value~SSasymp(time, yf, y0, log_koff), data=poi[poi$sample=="IP237" & poi$peaknames==hml_p,])
