#rm(list=ls())

#install.packages("remotes")
#remotes::install_github("ha0ye/rEDM", dependencies=TRUE)
#remotes::install_github("canankarakoc/EDMhelper", dependencies=TRUE)

#require(rEDMhelper)
#some will be loaded with the package 
require(rEDM)
library(plyr)
library(zoo)
library(viridis)
library(gridExtra)
library(scales)
library(ggplot2)
library(vegan)
library(reshape)
library(ggpubr)


#ggplot theme
mytheme<-theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        legend.text = element_text(size=12),
        legend.background = element_blank(),
        legend.title= element_text(size=12,face="bold"),
        plot.title = element_text(size=14, face="bold", hjust = 0.5),
        strip.text =element_text(size=14, face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#WD
#setwd("~/Documents/Everyfthing/r_package_EDMhelper/TEST")
setwd("C:/Users/karakoc/Documents/Everyfthing/r_package_EDMhelper/TEST_21.10")

#easy summary statistics
source("Summary_statistics_functions.R")

#load data
p_data <- read.table("p_combined.csv", header=T, sep=";",dec=",", stringsAsFactors=F)

#fill missing values with linear intrapolation in total ~3% of the all data used for the analysis
p_data <- as.data.frame(cbind(day=p_data$day,
                              ddply(p_data[,2:7],
                                    .(replicate=as.factor(p_data$replicate)),
                                    colwise(na.approx))))

#get segments (all observations are in one single segment)
#here replicates are the segments
sgm <- segment_data(p_data[,2])

#select columns for the tests
dnames <- colnames(p_data[,3:8])

#make standardized data
#here I use truezero=T, so that the real 0 measurements don't change 
p_data_norm <- normalize_data(p_data[,dnames], truezero = T)$xnorm

#run simplex on all columns
maxE_test <- round(sqrt(min(sgm[,2]-sgm[,1]))) #maximum E to test is sqrt(time series length)

simplex_out <- lapply(dnames, function(dnames){
  simplex(p_data_norm[,c(dnames)], E=2:maxE_test, lib=sgm, pred=sgm, tau=1, silent = TRUE)
})
names(simplex_out)<-dnames

#plot simplex results
simplex_p <- plot_simplex(simplex_out, predtype = "mae")+
  ylab("Predictability, MAE")+
  mytheme

#find best E for all time series
best_E <- find_max_E(simplex_out, predtype = "mae") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#run CCM for communities with predator
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(2210)

#*******************
#uncomment to run!!!
#*******************

#ccm_out <- ccm_easy(df=p_data_norm, best_E=best_E, lib_segments=sgm, pred_segments=sgm, num_samples=10000)
#saveRDS(ccm_out, file = "ccm_out.rds")

#read saved analysis
ccm_out <- readRDS("ccm_out.rds")

#summarize results
summaryout <- ccm_summary(ccm_output = ccm_out, predtype = "mae")

#run significance tests
sigout <- significance_test(ccm_output = ccm_out, predtype = "mae")

#this is just for a nice looking plot for the paper
summaryout$plot <-mapvalues(summaryout$direction, from = c("A causes Am", "Am causes A",
                                                           "A causes K", "K causes A",
                                                           "A causes Km", "Km causes A",
                                                           "A causes S", "S causes A",
                                                           "A causes T", "T causes A",
                                                           "Am causes K", "K causes Am",
                                                           "Am causes Km", "Km causes Am",
                                                           "Am causes S", "S causes Am",
                                                           "Am causes T", "T causes Am",
                                                           "K causes Km", "Km causes K",
                                                           "K causes S", "S causes K",
                                                           "K causes T", "T causes K",
                                                           "Km causes S", "S causes Km",
                                                           "Km causes T", "T causes Km",
                                                           "S causes T", "T causes S"),
                            to=c("A causes Am\nAm causes A", "A causes Am\nAm causes A",
                                 "A causes K\nK causes A", "A causes K\nK causes A",
                                 "A causes Km\nKm causes A", "A causes Km\nKm causes A",
                                 "A causes S\nS causes A", "A causes S\nS causes A",
                                 "A causes T\nT causes A", "A causes T\nT causes A",
                                 "Am causes K\nK causes Am", "Am causes K\nK causes Am",
                                 "Am causes Km\nKm causes Am", "Am causes Km\nKm causes Am",
                                 "Am causes S\nS causes Am", "Am causes S\nS causes Am",
                                 "Am causes T\nT causes Am", "Am causes T\nT causes Am",
                                 "K causes Km\nKm causes K", "K causes Km\nKm causes K",
                                 "K causes S\nS causes K", "K causes S\nS causes K",
                                 "K causes T\nT causes K", "K causes T\nT causes K",
                                 "Km causes S\nS causes Km", "Km causes S\nS causes Km",
                                 "Km causes T\nT causes Km", "Km causes T\nT causes Km",
                                 "S causes T\nT causes S", "S causes T\nT causes S"))


summaryout$dummy <-mapvalues(summaryout$direction, from = c("A causes Am", "Am causes A",
                                                            "A causes K", "K causes A",
                                                            "A causes Km", "Km causes A",
                                                            "A causes S", "S causes A",
                                                            "A causes T", "T causes A",
                                                            "Am causes K", "K causes Am",
                                                            "Am causes Km", "Km causes Am",
                                                            "Am causes S", "S causes Am",
                                                            "Am causes T", "T causes Am",
                                                            "K causes Km", "Km causes K",
                                                            "K causes S", "S causes K",
                                                            "K causes T", "T causes K",
                                                            "Km causes S", "S causes Km",
                                                            "Km causes T", "T causes Km",
                                                            "S causes T", "T causes S"),
                                                    to=rep(c(1,2), times=15))



ccm_fig <- ggplot(summaryout, aes(x=lib_size, y=Q50., group=direction))+
  geom_ribbon(aes(ymin=Q15.86553., ymax=Q84.13447., x=lib_size),alpha=0.1)+
  geom_line(aes(color=dummy))+
  facet_wrap(~plot)+
  scale_color_manual(values = c("darkred","midnightblue"))+
  scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1))+
  labs(y="Predictability, MAE", x="Library Size")+
  mytheme
ccm_fig

#~~~~~~~~~~~~~~~~~~
#run time-delay CCM
#~~~~~~~~~~~~~~~~~~

#uncomment to run
#ccm_del_out <- ccm_delay_easy(df=p_data_norm, best_E = best_E, lib_segments=sgm,
#                            pred_segments=sgm, tp = -5:0, num_samples=10000)
#saveRDS(ccm_del_out, file = "ccm_del_out.rds")

#read the saved analysis
ccm_del_out<-readRDS("ccm_del_out.rds")

#summarize results
summaryout_del <- ccm_summary(ccm_output = ccm_del_out, predtype = "mae")

#run significance tests
sigout_del <- significance_test(ccm_output = ccm_del_out, predtype = "mae")

#supplementary table
tableS <- subset(sigout_del, sigout_del$tp<=0 & sigout_del$pval<=0.05)
capture.output(tableS, file="tableS.txt")

#find best prediction lag
sigout_tp <- find_tp(sigout_del, cutoff = 0.05)

#plot result
#ccm_easy_plot(summaryout_del)

#plot prediction skill at the maximum library size as a function of tp
#ccm_delay_plot(summaryout_del)

#this is just for a fancy plot for the paper
summaryout_del$plot <- mapvalues(summaryout_del$direction, from = c("A causes Am", "Am causes A",
                                                                   "A causes K", "K causes A",
                                                                   "A causes Km", "Km causes A",
                                                                   "A causes S", "S causes A",
                                                                   "A causes T", "T causes A",
                                                                   "Am causes K", "K causes Am",
                                                                   "Am causes Km", "Km causes Am",
                                                                   "Am causes S", "S causes Am",
                                                                   "Am causes T", "T causes Am",
                                                                   "K causes Km", "Km causes K",
                                                                   "K causes S", "S causes K",
                                                                   "K causes T", "T causes K",
                                                                   "Km causes S", "S causes Km",
                                                                   "Km causes T", "T causes Km",
                                                                   "S causes T", "T causes S"),
                                to=c("A causes Am\nAm causes A", "A causes Am\nAm causes A",
                                     "A causes K\nK causes A", "A causes K\nK causes A",
                                     "A causes Km\nKm causes A", "A causes Km\nKm causes A",
                                     "A causes S\nS causes A", "A causes S\nS causes A",
                                     "A causes T\nT causes A", "A causes T\nT causes A",
                                     "Am causes K\nK causes Am", "Am causes K\nK causes Am",
                                     "Am causes Km\nKm causes Am", "Am causes Km\nKm causes Am",
                                     "Am causes S\nS causes Am", "Am causes S\nS causes Am",
                                     "Am causes T\nT causes Am", "Am causes T\nT causes Am",
                                     "K causes Km\nKm causes K", "K causes Km\nKm causes K",
                                     "K causes S\nS causes K", "K causes S\nS causes K",
                                     "K causes T\nT causes K", "K causes T\nT causes K",
                                     "Km causes S\nS causes Km", "Km causes S\nS causes Km",
                                     "Km causes T\nT causes Km", "Km causes T\nT causes Km",
                                     "S causes T\nT causes S", "S causes T\nT causes S"))


summaryout_del$dummy <-mapvalues(summaryout_del$direction, from = c("A causes Am", "Am causes A",
                                                                    "A causes K", "K causes A",
                                                                    "A causes Km", "Km causes A",
                                                                    "A causes S", "S causes A",
                                                                    "A causes T", "T causes A",
                                                                    "Am causes K", "K causes Am",
                                                                    "Am causes Km", "Km causes Am",
                                                                    "Am causes S", "S causes Am",
                                                                    "Am causes T", "T causes Am",
                                                                    "K causes Km", "Km causes K",
                                                                    "K causes S", "S causes K",
                                                                    "K causes T", "T causes K",
                                                                    "Km causes S", "S causes Km",
                                                                    "Km causes T", "T causes Km",
                                                                    "S causes T", "T causes S"),
                                                                to=rep(c(1,2), times=15))
#plot
#here I just plot the max. library size 
#summaryout_del_sub <- subset(summaryout_del, summaryout_del$lib_size==10 | summaryout_del$lib_size==100)
summaryout_del_sub <- subset(summaryout_del, summaryout_del$lib_size==100)
as.numeric(summaryout_del_sub$tp)

time_delay_ccm_fig <- ggplot(summaryout_del_sub, aes(x=tp, y=Q50.,
                                                     color=dummy))+
  geom_vline(xintercept=0, linetype="dashed", color="grey")+
  geom_ribbon(aes(ymin=Q15.86553., ymax=Q84.13447., x=tp),alpha=0.1, linetype=0)+
  geom_line()+
  facet_wrap(~plot)+
  scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1))+
  scale_color_manual(values = c("darkred", "midnightblue"))+
  labs(x="Prediction lag", y="Predictability, MAE")+
  mytheme
time_delay_ccm_fig

#~~~~~~~~~~~~~~~~~~~
#s-maps coefficients
#~~~~~~~~~~~~~~~~~~~

smap_out <- lapply(dnames, function(dnames){
  s_map(p_data_norm[,c(dnames)], E=best_E[dnames], lib=sgm, pred=sgm, tau=1, silent = TRUE)
})
names(smap_out) <- dnames

#find best theta for all time series
best_theta <- find_max_theta(smap_out, predtype = "mae")
names(best_theta) <- dnames

#plot predictability as a function of theta
theta_p <- plot_smap(smap_out, predtype = "mae")+
  ylab("Predictability, MAE")+
  mytheme

#calculate s-map coefficients
smap_coef_out_FULL_tp <- get_smap_coef(df=p_data_norm, lib_segments = sgm,
                                      sigout = sigout_tp, best_E = best_E, best_theta = best_theta, selfref = TRUE)

#plot coefficients
#here I select the significant interactions manually using sigout_tp and smap results,
#and and create a data set with all the coefficients

#K on A
a <- mean(smap_coef_out_FULL_tp$smap_out[[11]]$smap_coefficients[[1]]$K_t2, na.rm=T)
#K on Km
b <- mean(smap_coef_out_FULL_tp$smap_out[[13]]$smap_coefficients[[1]]$K_t1, na.rm=T)
#Km on A
c <- mean(smap_coef_out_FULL_tp$smap_out[[16]]$smap_coefficients[[1]]$Km_t0, na.rm=T)
#Km on T
d <- mean(smap_coef_out_FULL_tp$smap_out[[20]]$smap_coefficients[[1]]$Km_t0, na.rm=T)
#T on A
e <- mean(smap_coef_out_FULL_tp$smap_out[[26]]$smap_coefficients[[1]]$T_t1, na.rm=T)
#T on Km
f <- mean(smap_coef_out_FULL_tp$smap_out[[29]]$smap_coefficients[[1]]$T_t1, na.rm=T)
#T on S
g <- mean(smap_coef_out_FULL_tp$smap_out[[30]]$smap_coefficients[[1]]$T_t0, na.rm=T)

#function plots prediction segments side by side
#i just get the data from the plot with a small ggplot hack and save for later custom plotting 
#K on A
p1 <- smap_coef_plot(index_sp=11, smap_coef_out=smap_coef_out_FULL_tp,
                     covar="time", ycol="K_t2", xlab="time")
KA <- as.data.frame(cbind((ggplot_build(p1)$data[[2]])[,2:3], y=(ggplot_build(p1)$data[[3]])[,2]))

#K on Km
p2 <- smap_coef_plot(index_sp=13, smap_coef_out=smap_coef_out_FULL_tp,
                     covar="time", ycol="K_t1", xlab="time")
KKm <- as.data.frame(cbind((ggplot_build(p2)$data[[2]])[,2:3], y=(ggplot_build(p2)$data[[3]])[,2]))

#Km on A
p3 <- smap_coef_plot(index_sp=16, smap_coef_out=smap_coef_out_FULL_tp,
                     covar="time", ycol="Km_t0", xlab="time")
KmA <- as.data.frame(cbind((ggplot_build(p3)$data[[2]])[,2:3], y=(ggplot_build(p3)$data[[3]])[,2]))

#Km on T
p4 <- smap_coef_plot(index_sp=20, smap_coef_out=smap_coef_out_FULL_tp,
                     covar="time", ycol="Km_t0", xlab="time")
KmT <- as.data.frame(cbind((ggplot_build(p4)$data[[2]])[,2:3], y=(ggplot_build(p4)$data[[3]])[,2]))

#T on A
p5 <- smap_coef_plot(index_sp=26, smap_coef_out=smap_coef_out_FULL_tp,
                     covar="time", ycol="T_t1", xlab="time")
TA <- as.data.frame(cbind((ggplot_build(p5)$data[[2]])[,2:3], y=(ggplot_build(p5)$data[[3]])[,2]))

#T on Km
p6 <- smap_coef_plot(index_sp=29, smap_coef_out=smap_coef_out_FULL_tp,
                     covar="time", ycol="T_t1", xlab="time")
TKm <- as.data.frame(cbind((ggplot_build(p6)$data[[2]])[,2:3], y=(ggplot_build(p6)$data[[3]])[,2]))

#T on S
p7 <- smap_coef_plot(index_sp=30, smap_coef_out=smap_coef_out_FULL_tp,
                     covar="time", ycol="T_t0", xlab="time")
TS <- as.data.frame(cbind((ggplot_build(p7)$data[[2]])[,2:3], y=(ggplot_build(p7)$data[[3]])[,2]))

allcoef <- rbind(KA, KKm, KmA, KmT, TA, TKm, TS)
allcoef$day <- rep(c(3:34), times=21)
allcoef$replicates <- as.factor(rep(c(1,2,3), each=32, times=7))
allcoef$var <- rep(c("A/K", "Km/K", "A/Km", "T/Km", "A/T", "Km/T", "S/T"), each=96)

values_deriv <- c(expression(~~partialdiff~A/~partialdiff~K), 
                  expression(~~partialdiff~Km/~partialdiff~K),
                  expression(~~partialdiff~A/~partialdiff~Km),
                  expression(~~partialdiff~T/~partialdiff~Km),
                  expression(~~partialdiff~A/~partialdiff~T),
                  expression(~~partialdiff~Km/~partialdiff~T),
                  expression(~~partialdiff~S/~partialdiff~T))

#Plot
aggdata <- aggregate(allcoef, by=list(allcoef$var,allcoef$day), FUN=mean, na.rm=TRUE)

aggcoefplot <- ggplot(aggdata, aes(x=day, y=y, color=Group.1, linetype=Group.1, group=Group.1))+
  geom_hline(yintercept=0, color="grey")+
  geom_ribbon(aes(ymin=ymin, ymax=ymax, x=day),alpha=0.1, linetype=0)+
  geom_line()+
  scale_color_viridis(discrete=TRUE, labels=values_deriv)+
  scale_linetype_discrete(labels=values_deriv)+
  scale_y_continuous(limits=c(-2.15,0.75))+
  scale_x_continuous(limits=c(3,34))+
  labs(x="Time (days)", y="S-map coefficients")+
  mytheme+
  theme(legend.position = c(.8,0.25), legend.title = element_blank())+
  ggtitle("Predator present")
aggcoefplot

#~~~~~~~~~~~~~~~~~~~~~~~~
#local stability analysis
#~~~~~~~~~~~~~~~~~~~~~~~~

stability_out <- stability_analysis(smap_coef_out_FULL_tp)

p.stab.data <- stability_out$all_indices
p.stab.data$replicates <- as.factor(rep(c(1,2,3), each=32))
p.stab.data$day <- rep(3:34, times=3)

#averaged over replicates
aggdata.stab <-aggregate(p.stab.data, by=list(p.stab.data$day), FUN=mean, na.rm=TRUE)

#for the plot with error bars
p.dyn.stab <- summarySE(p.stab.data, measurevar="local_stability", groupvars="day", na.rm=T)
p.mean.int <- summarySE(p.stab.data, measurevar="mean_interactions", groupvars="day", na.rm=T)
p.weak.int <- summarySE(p.stab.data, measurevar="weak_interactions", groupvars="day", na.rm=T)

#*****************************************************************
#*******same analysis with the communities without predator*******
#*****************************************************************

#load data
wp_data <- read.table("wp_combined.csv", header=T,sep=";",dec=",", stringsAsFactors=F)

#Fill missing values with linear intrapolation 
wp_data <- as.data.frame(cbind(day=wp_data$day,
                               ddply(wp_data[,c(2,3,5)],
                                     .(replicate=as.factor(wp_data$replicate)),
                                     colwise(na.approx))))

#get segments (all observations are in one single segment)
sgm_wp <- segment_data(wp_data[,2])

#select the columns for the analysis
dnames_wp <- colnames(wp_data[,3:5])

#make standardized data
wp_data_norm <- normalize_data(wp_data[,dnames_wp], truezero = T)$xnorm

#run simplex on all columns
maxE_test_wp <- round(sqrt(min(sgm_wp[,2]-sgm_wp[,1]))) #maximum E to test is sqrt(time series length)

simplex_out_wp <- lapply(dnames_wp, function(dnames_wp){
  simplex(wp_data_norm[,c(dnames_wp)], E=2:maxE_test_wp, lib=sgm_wp, pred=sgm_wp, tau=1, silent = TRUE)
})
names(simplex_out_wp) <- dnames_wp

simplex_wp <- plot_simplex(simplex_out_wp, predtype = "mae")+
  ylab("Predictability, MAE")+
  mytheme

#find best E for all time series
best_E_wp <- find_max_E(simplex_out_wp, predtype = "mae")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#supplementary figure
#simplex results combined
ggarrange(simplex_wp, simplex_p, ncol = 1, nrow = 2, heights = c(0.55,1), labels = c("a", "b"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~
#Run CCM
#~~~~~~~

#*******************
#uncomment to run!!!
#*******************

#ccm_out_wp <- ccm_easy(df=wp_data_norm, best_E=best_E_wp, lib_segments=sgm_wp, pred_segments=sgm_wp, num_samples=10000)
#saveRDS(ccm_out_wp, file = "ccm_out_wp.rds")

#read the saved analysis
ccm_out_wp <- readRDS("ccm_out_wp.rds")

#summarize results
summaryout_wp <- ccm_summary(ccm_output = ccm_out_wp, predtype = "mae")

#run significance tests
sigout_wp <- significance_test(ccm_output = ccm_out_wp, predtype = "mae")

#plot result
#ccm_easy_plot(summaryout_wp)

#this is just for fancy plot for the paper 
summaryout_wp$plot <- mapvalues(summaryout_wp$direction, from = c("A causes Am", "Am causes A",
                                                           "A causes S", "S causes A",
                                                           "Am causes S", "S causes Am"),
                            to=c("A causes Am\nAm causes A", "A causes Am\nAm causes A",
                                 "A causes S\nS causes A", "A causes S\nS causes A",
                                 "Am causes S\nS causes Am", "Am causes S\nS causes Am"))

summaryout_wp$dummy <- mapvalues(summaryout_wp$direction, from = c("A causes Am", "Am causes A",
                                                            "A causes S", "S causes A",
                                                            "Am causes S", "S causes Am"),
                                                          to=rep(c(1,2), times=3))

ccm_fig_wp <- ggplot(summaryout_wp, aes(x=lib_size, y=Q50., group=direction))+
  geom_ribbon(aes(ymin=Q15.86553., ymax=Q84.13447., x=lib_size),alpha=0.1)+
  geom_line(aes(color=dummy))+
  facet_wrap(~plot)+
  scale_color_manual(values = c("darkred","midnightblue"))+
  labs(y="Predictability, MAE", x="Library Size")+
  mytheme

#~~~~~~~~~~~~~~~~~~
#Run time-delay CCM
#~~~~~~~~~~~~~~~~~~

#*******************
#uncomment to run!!!
#*******************

#ccm_del_out_wp <- ccm_delay_easy(df=wp_data_norm, best_E = best_E_wp, lib_segments=sgm_wp,
#                            pred_segments=sgm_wp, tp = -5:0, num_samples=10000)
#saveRDS(ccm_del_out_wp, file = "ccm_del_out_wp.rds")

#read the saved analysis
ccm_del_out_wp <- readRDS("ccm_del_out_wp.rds")

#summarize results
summaryout_del_wp <- ccm_summary(ccm_output = ccm_del_out_wp, predtype = "mae")

#run significance tests
sigout_del_wp <- significance_test(ccm_output = ccm_del_out_wp, predtype = "mae")

#find best prediction lag
sigout_tp_wp <- find_tp(sigout_del_wp, cutoff = 0.05)

#supplementary table
tableS_wp <- subset(sigout_del_wp, sigout_del_wp$tp<=0 & sigout_del_wp$pval<=0.05)
capture.output(tableS_wp, file="tableS_wp.txt")

#plot result
#ccm_easy_plot(summaryout_del_wp)

#plot prediction skill at the maximum library size as a function of tp
#ccm_delay_plot(summaryout_del_wp)

#this is just for fancy plot for the paper 
summaryout_del_wp$plot <- mapvalues(summaryout_del_wp$direction, from = c("A causes Am", "Am causes A",
                                                                   "A causes S", "S causes A",
                                                                   "Am causes S", "S causes Am"),
                                to=c("A causes Am\nAm causes A", "A causes Am\nAm causes A",
                                     "A causes S\nS causes A", "A causes S\nS causes A",
                                     "Am causes S\nS causes Am", "Am causes S\nS causes Am"))


summaryout_del_wp$dummy <- mapvalues(summaryout_del_wp$direction, from = c("A causes Am", "Am causes A",
                                                                    "A causes S", "S causes A",
                                                                    "Am causes S", "S causes Am"),
                                                                   to=rep(c(1,2), times=3))

summaryout_del_wp_sub <- subset(summaryout_del_wp, summaryout_del_wp$lib_size==100)
as.numeric(summaryout_del_wp_sub$tp)

ccm_fig_del_wp <- ggplot(summaryout_del_wp_sub, aes(x=tp, y=Q50., group=direction))+
  geom_ribbon(aes(ymin=Q15.86553., ymax=Q84.13447., x=tp),alpha=0.1)+
  geom_line(aes(color=dummy))+
  facet_wrap(~plot)+
  scale_color_manual(values = c("darkred","midnightblue"))+
  labs(y="Predictability, MAE", x="Prediction lag")+
  mytheme
ccm_fig_del_wp

#~~~~~~~~~~~~~~~~~~~
#s-maps coefficients 
#~~~~~~~~~~~~~~~~~~~

smap_out_wp <- lapply(dnames_wp, function(dnames_wp){
  s_map(wp_data_norm[,c(dnames_wp)], E=best_E_wp[dnames_wp], lib=sgm_wp, pred=sgm_wp, tau=1, silent = TRUE)
})
names(smap_out_wp) <- dnames_wp

#find best theta for all time series
best_theta_wp <- find_max_theta(smap_out_wp, predtype = "mae")
names(best_theta_wp) <- dnames_wp

theta_wp <- plot_smap(smap_out_wp, predtype = "mae")+
  ylab("Predictability, MAE")+
  mytheme

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#supplementary figure
#theta combined
ggarrange(theta_wp, theta_p, ncol = 1, nrow = 2, heights = c(0.55,1), labels = c("a", "b"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#calculate s-map coefficients
smap_coef_out_FULL_tp_wp <- get_smap_coef(df=wp_data_norm, lib_segments = sgm_wp,
                                     sigout = sigout_tp_wp, best_E = best_E_wp, best_theta = best_theta_wp, selfref = TRUE)

#there is only one interaction
#A on Am
h <- mean(smap_coef_out_FULL_tp_wp$smap_out[[1]]$smap_coefficients[[1]]$A_t0, na.rm = T)

#A on Am
p8 <- smap_coef_plot(index_sp=1, smap_coef_out=smap_coef_out_FULL_tp_wp,
                     covar="time", ycol="A_t0", xlab="time")
AAm <- as.data.frame(cbind((ggplot_build(p8)$data[[2]])[,2:3], y=(ggplot_build(p8)$data[[3]])[,2]))

allcoef_wp <- AAm
allcoef_wp$day <- rep(c(3:34), times=3)
allcoef_wp$replicates <- as.factor(rep(c(1,2,3), each=32, times=1))
allcoef_wp$var <- rep("Am/A", times=96)

#Plot
aggdata_wp <- aggregate(allcoef_wp, by=list(allcoef_wp$var,allcoef_wp$day), FUN=mean, na.rm=TRUE)

aggcoefplot_wp <- ggplot(aggdata_wp, aes(x=day, y=y, color=Group.1, group=Group.1))+
  geom_ribbon(aes(ymin=ymin, ymax=ymax, x=day),alpha=0.1, linetype=0)+
  geom_hline(yintercept=0, color="grey")+
  geom_line()+
  scale_y_continuous(limits=c(-2.15,0.75))+
  scale_x_continuous(limits=c(3,34))+
  scale_color_manual(values="black", labels = expression(~partialdiff~Am/~partialdiff~A))+
  labs(x="Time (days)", y="S-map coefficients")+
  mytheme+
  theme(legend.position = c(.8,0.3), legend.title=element_blank())+
  ggtitle("Predator absent")
aggcoefplot_wp

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 3
#combined plot for the paper 
ggarrange(aggcoefplot_wp, aggcoefplot, ncol = 2, nrow = 1, labels = c("a", "b"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 2
#interaction strengths for the network plot
#I drew in Inkscape
round(rescale(abs(c(a, b, c, d, e, f, g, h)), c(0.4,2)), 2)
#0.55, 0.43, 1.48, 2.00, 0.75, 0.45, 0.40, 0.53
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~
#local stability analysis
#~~~~~~~~~~~~~~~~~~~~~~~~

stability_out_wp <- stability_analysis(smap_coef_out_FULL_tp_wp)

wp.stab.data <- stability_out_wp$all_indices
wp.stab.data$replicates <- as.factor(rep(c(1,2,3), each=32))
wp.stab.data$day <- rep(3:34, times=3)

# for the plot with error bars
wp.dyn.stab <- summarySE(wp.stab.data, measurevar="local_stability", groupvars="day", na.rm=T)
wp.mean.int <- summarySE(wp.stab.data, measurevar="mean_interactions", groupvars="day", na.rm=T)
wp.weak.int <- summarySE(wp.stab.data, measurevar="weak_interactions", groupvars="day", na.rm=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#add diversity and error bars
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p_data$div  <- diversity(p_data[,3:8], index="invsimpson")
wp_data$div <- diversity(wp_data[,3:5], index="invsimpson")

#summary diversity
p_div_sum  <- summarySE(p_data, measurevar="div", groupvars="day", na.rm=T)
wp_div_sum <- summarySE(wp_data, measurevar="div", groupvars="day", na.rm=T)
p_div_sum$var <- rep("present", nrow(p_div_sum))
wp_div_sum$var <- rep("absent", nrow(wp_div_sum))

#combined data
div_comb <- as.data.frame(rbind(p_div_sum, wp_div_sum ))
comb.div <- ggplot(div_comb, aes(x=day,y=div, group=var, color=var))+
  geom_errorbar(aes(ymin=div-sd, ymax=div+sd),
                color="grey60", width=0.2, size=0.5, alpha=.4)+
  geom_line(size=1)+
  scale_color_manual(values=c("grey25", "firebrick"))+
  #facet_grid(~var)+
  scale_x_continuous(limits=c(3,35))+
  geom_point(size=2)+
  ylab("Diversity (inv. Simpson)")+
  xlab("Time (days)")+
  mytheme+
  theme(legend.position = "none")
comb.div

#~~~~~~~~~~~~~~~~~~~~~~~
#Interactions error bars
#~~~~~~~~~~~~~~~~~~~~~~~

#stability

p.dyn.stab$var  <- rep("present", nrow(p.dyn.stab))
wp.dyn.stab$var <- rep("absent", nrow(wp.dyn.stab))
stab_combined   <- as.data.frame(rbind(p.dyn.stab, wp.dyn.stab))

comb.stab <- ggplot(stab_combined, aes(x=day,y=local_stability, group=var, color=var))+
  geom_hline(yintercept = 1, linetype="dashed", color="grey")+
  geom_errorbar(aes(ymin=local_stability-sd, ymax=local_stability+sd),
                color="grey60", width=0.2, size=0.5, alpha=.4)+
  geom_line(size=1)+
  geom_point(size=2)+
  annotate("text", size=4, angle = 90, hjust = 0, x = 3, y = 0.65,  label = "stable")+
  annotate("text", size=4, angle = 90, hjust = 0, x = 3, y = 1.1, label = "unstable")+
  scale_color_manual(name= "Predator", values=c("grey25", "firebrick"), labels=c("Absent", "Present"))+
  scale_x_continuous(limits=c(3,35))+
  ylab("Dynamic stability")+
  xlab("Time (days)")+
  mytheme+
  theme(legend.position = c(0.25,0.45))
comb.stab

#mean interactions

p.mean.int$var  <- rep("present", nrow(p.mean.int))
wp.mean.int$var <- rep("absent", nrow(wp.mean.int))
mean_combined   <- as.data.frame(rbind(p.mean.int, wp.mean.int))

comb.mean<- ggplot(mean_combined, aes(x=day,y=mean_interactions, group=var, color=var))+
  geom_errorbar(aes(ymin=mean_interactions-sd, ymax=mean_interactions+sd),
                color="grey60", width=0.2, size=0.5, alpha=.4)+
  geom_line(size=1)+
  geom_point(size=2)+
  scale_color_manual(values=c("grey25", "firebrick"))+
  scale_x_continuous(limits=c(3,35))+
  ylab("Mean interaction strength")+
  xlab("Time (days)")+
  mytheme+
  theme(legend.position = "none")
comb.mean

#weak interactions

p.weak.int$var  <- rep("present", nrow(p.weak.int))
wp.weak.int$var <- rep("absent", nrow(wp.weak.int))
weak_combined   <- as.data.frame(rbind(p.weak.int, wp.weak.int))

comb.weak<- ggplot(weak_combined, aes(x=day,y=weak_interactions, group=var, color=var))+
  geom_errorbar(aes(ymin=weak_interactions-sd, ymax=weak_interactions+sd),
                color="grey60", width=0.2, size=0.5, alpha=.4)+
  geom_line(size=1)+
  geom_point(size=2)+
  annotate("text", size=4, angle = 90, hjust = 0, x = 4, y = 0.03, label = "dominance of weak interactions")+
  geom_segment(aes(x = 3, y = 0.34, xend = 3, yend = 0.009), arrow = arrow(length = unit(0.2, "inches")))+
  scale_color_manual(values=c("grey25", "firebrick"))+
  scale_x_continuous(limits=c(3,35))+
  ylab("Mean:max. int. strength")+
  xlab("Time (days)")+
  mytheme+
  theme(legend.position = "none")
comb.weak

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 4
#all indices
ggarrange(comb.stab, comb.div, comb.mean, comb.weak, ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#test the effects on stability
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

p.test.dat  <- cbind(p.stab.data, diversity=p_data$div)[,c("day", "replicates", "local_stability",
                                                           "mean_interactions", "weak_interactions", "diversity")]
wp.test.dat <- cbind(wp.stab.data, diversity=wp_data$div)[,c("day", "replicates", "local_stability",
                                                             "mean_interactions", "weak_interactions", "diversity")]
p.dnames  <- colnames(p.test.dat[,3:6])
wp.dnames <- colnames(wp.test.dat[,3:6])

#make standardized data
#p_test_data_norm  <- normalize_data(p.test.dat[,p.dnames], truezero = F)$xnorm
#wp_test_data_norm <- normalize_data(wp.test.dat[,wp.dnames], truezero = T)$xnorm

#get segments (all observations are in one single segment)
p.sgm <- segment_data(p.test.dat[,2])
wp.sgm <- segment_data(wp.test.dat[,2])

#run simplex on all columns
maxE_test <-round(sqrt(min(p.sgm[,2]-p.sgm[,1]))) #maximum E to test is sqrt(time series length)

#run simplex for communities with predator
p.simplex_out <- lapply(p.dnames, function(p.dnames){
  simplex(p.test.dat[,c(p.dnames)], E=2:maxE_test, lib=p.sgm, pred=p.sgm, tau=1, silent = TRUE)
})
names(p.simplex_out) <- p.dnames

#plot predictability as a function of E
simplex_ptest <- plot_simplex(p.simplex_out, predtype = "mae")+
  ylab("Predictability, MAE")+
  mytheme

#run simplex for communities without predator
wp.simplex_out <- lapply(wp.dnames, function(wp.dnames){
  simplex(wp.test.dat[,c(wp.dnames)], E=2:maxE_test, lib=wp.sgm, pred=wp.sgm, tau=1, silent = TRUE)
})
names(wp.simplex_out) <-wp.dnames

#plot predictability as a function of E
simplex_wptest <- plot_simplex(wp.simplex_out, predtype = "mae")+
  ylab("Predictability, MAE")+
  mytheme

#find best E for all time series
p.best_E <- find_max_E(p.simplex_out, predtype = "mae") 
wp.best_E <- find_max_E(wp.simplex_out, predtype = "mae")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#supplementary figure
#simplex results combined
ggarrange(simplex_wptest, simplex_ptest, ncol = 1, nrow = 2, heights = c(0.5,0.5), labels = c("a", "b"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#*******
#run CCM
#*******

#*******************
#uncomment to run!!!
#*******************

#p.ccm_out <- ccm_easy(df=p.test.dat[,c(p.dnames)], best_E=p.best_E, lib_segments=p.sgm, pred_segments=p.sgm, num_samples=10000)
#saveRDS(p.ccm_out, file = "test_interaction_ccm_out_p.rds")

#read saved analysis
p.ccm_out <- readRDS("test_interaction_ccm_out_p.rds")

#wp.ccm_out <- ccm_easy(df=wp.test.dat[,c(wp.dnames)], best_E=wp.best_E, lib_segments=wp.sgm, pred_segments=wp.sgm, num_samples=10000)
#saveRDS(wp.ccm_out, file = "test_interaction_ccm_out_wp.rds")

#read saved analysis
wp.ccm_out <- readRDS("test_interaction_ccm_out_wp.rds")

#summarize results
p.summaryout <- ccm_summary(ccm_output = p.ccm_out, predtype = "mae")
wp.summaryout <- ccm_summary(ccm_output = wp.ccm_out, predtype = "mae")

#run significance tests
p.sigout  <- significance_test(ccm_output = p.ccm_out, predtype = "mae")
wp.sigout <- significance_test(ccm_output = wp.ccm_out, predtype = "mae")

#******************
#run time-delay CCM
#******************

#*******************
#uncomment to run!!!
#*******************

#p.ccm_del_out <- ccm_delay_easy(df = p.test.dat[,c(p.dnames)], best_E = p.best_E, lib_segments=p.sgm, 
#                                pred_segments=p.sgm, tp = -5:0, num_samples=10000)
#saveRDS(p.ccm_del_out, file = "ccm_del_out_test_p.rds")

#read saved analysis
p.ccm_del_out <- readRDS("ccm_del_out_test_p.rds")

#wp.ccm_del_out <-ccm_delay_easy(df = wp.test.dat[,c(wp.dnames)], best_E = wp.best_E, lib_segments=wp.sgm,
#                           pred_segments=wp.sgm, tp = -5:0, num_samples=10000)
#saveRDS(wp.ccm_del_out, file = "ccm_del_out_test_wp.rds")

#read saved analysis
wp.ccm_del_out <- readRDS("ccm_del_out_test_wp.rds")

#summarize results
p.summaryout_del <- ccm_summary(ccm_output = p.ccm_del_out, predtype = "mae")
wp.summaryout_del <- ccm_summary(ccm_output = wp.ccm_del_out, predtype = "mae")

#run significance tests
p.sigout_del <- significance_test(ccm_output = p.ccm_del_out, predtype = "mae")
p.sigout_tp <- find_tp(p.sigout_del, cutoff = 0.05)

wp.sigout_del <- significance_test(ccm_output = wp.ccm_del_out, predtype = "mae")
wp.sigout_tp <- find_tp(wp.sigout_del, cutoff = 0.05)


#subset the summaries with the significant causations and prediction lags 
p.summaryout_del_sub_divloc <- subset(p.summaryout_del, p.summaryout_del$direction=="diversity causes local_stability" &
                                 p.summaryout_del$tp==-5)

p.summaryout_del_sub_locdiv <- subset(p.summaryout_del, p.summaryout_del$direction=="local_stability causes diversity" &
                                     p.summaryout_del$tp==-5)

p.summaryout_del_sub_locweak <- subset(p.summaryout_del, p.summaryout_del$direction=="local_stability causes weak_interactions" &
                                     p.summaryout_del$tp==0)

p.summaryout_del_sub_meandiv <- subset(p.summaryout_del, p.summaryout_del$direction=="mean_interactions causes diversity" &
                                      p.summaryout_del$tp==-2)

p.summaryout_del_sub_meanweak <- subset(p.summaryout_del, p.summaryout_del$direction=="mean_interactions causes weak_interactions" &
                                      p.summaryout_del$tp==-5)

comb.p_test <- as.data.frame (rbind(p.summaryout_del_sub_divloc, p.summaryout_del_sub_locdiv, 
                                    p.summaryout_del_sub_locweak, p.summaryout_del_sub_meandiv,
                                    p.summaryout_del_sub_meanweak))
comb.p_test$var <- rep("present", nrow(comb.p_test))

comb.wp_test <-subset(wp.summaryout_del, wp.summaryout_del$direction=="diversity causes mean_interactions" &
                                    wp.summaryout_del$tp==-4)
comb.wp_test$var <- rep("absent", nrow(comb.wp_test))

#combined figure for the paper with the delayed responses

comb.all_plot <- as.data.frame(rbind(comb.p_test, comb.wp_test))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Figure 5
comb.all_plot$var <- factor(comb.all_plot$var, labels = c("Predator absent", "Predator present"))
comb.all_plot$direction <- factor(comb.all_plot$direction, labels = c("Diversity causes local stability", 
                                  "Diversity causes mean int. strength", "Local stability causes diversity", 
                                  "Local stability causes weak int. strength", "Mean int. strength causes diversity", 
                                  "Mean int. strength causes weak int. strength"))

comb.ccm_fig <- ggplot(comb.all_plot, aes(x=lib_size, y=Q50., group=direction))+
  geom_ribbon(aes(ymin=Q15.86553., ymax=Q84.13447., x=lib_size),alpha=0.1)+
  geom_line(aes(color=direction, linetype=direction))+
  facet_wrap(~var)+
  labs(y="Predictability, MAE", x="Library Size")+
  scale_color_manual(values=c("darkred", "black", "darkred", "midnightblue", "grey50", "darkgreen"))+
  scale_linetype_manual(values=c("solid", "solid", "dashed", "solid", "solid", "solid"))+
  mytheme+
  theme(legend.position = c(.25,.3), legend.direction="vertical", legend.title=element_blank()) 
  
comb.ccm_fig
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Time series plot for the paper
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Read the data above again!
p_data
wp_data

p_ts <- melt(p_data, id=c("day","replicate"))
np_ts <- melt(wp_data, id=c("day","replicate"))

dfc_p <- summarySE(p_ts, measurevar="value", groupvars=c("day","variable"),na.rm=T)
dfc_np <- summarySE(np_ts, measurevar="value", groupvars=c("day","variable"),na.rm=T)


#With predator
pred.mean <- ggplot(dfc_p, aes(x=day,y=log10(value)))+
  geom_errorbar(aes(ymin=log10((value-sd)),ymax=log10((value+sd))),
                color="grey50", width=0.2, size=0.5, alpha=0.4)+
  geom_line(size=1, aes(color=variable), data = dfc_p[!is.na(dfc_p$value),])+
  geom_point(size=3, aes(fill=variable, shape=variable))+
  ylab("Abundance")+
  xlab("Time (days)")+
  scale_shape_manual(values=c(21,21,21,21,21,24), name = "Species")+
  scale_fill_manual(values=c("darkred", "lightpink3","midnightblue", "skyblue3",
                             "darkgreen", "grey25"), name = "Species")+
  scale_color_manual(values=c("darkred", "lightpink3","midnightblue", "skyblue3",
                              "darkgreen", "grey25"), name = "Species")+
  scale_y_continuous(limits=c(4, 9))+
  mytheme

pred.mean

#Without predator
comp.mean <- ggplot(dfc_np, aes(x=day,y=log10(value)))+
  geom_errorbar(aes(ymin=log10((value-sd)),ymax=log10((value+sd))),
                color="grey50", width=0.2, size=0.5, alpha=0.4)+
  geom_line(size=1, aes(color=variable), data = dfc_np[!is.na(dfc_np$value),])+
  geom_point(size=3, aes(fill=variable, shape=variable))+
  ylab("Abundance")+
  xlab("Time (days)")+
  scale_shape_manual(values=c(21,21,21,21,21,24), name = "Species")+
  scale_fill_manual(values=c("darkred", "lightpink3",  "midnightblue",
                             "darkgreen"), name = "Species")+
  scale_color_manual(values=c("darkred", "lightpink3", "midnightblue", "darkgreen"), name = "Species")+
  scale_y_continuous(limits=c(4, 9))+
  mytheme

comp.mean

grid.arrange(comp.mean,pred.mean, ncol=2)
