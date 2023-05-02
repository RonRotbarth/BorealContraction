#############################Introduction######################################
# R script containing data preparation, exploration, visualisation and analysis for the paper titled 'North American boreal forests: Northern expansion is not compensating for southern declines' published in Nature Communications 2023
# 
# Authors: Ronny Rotbarth*, Egbert H. van Nes, Marten Scheffer, Jane Uhd Jepsen, 
# Ole Petter Laksforsmo Vindstad, Chi Xu, Milena Holmgren
# 
# *corresponding author
# 
# Publication year: 2023
# 
# R script version dated 1 May 2023
# 
# All analyses were performed in RStudio Version 1.1.463
# 
# R version: 4.2.0
# 
# 
# Paper abstract: Climate change is expected to shift the boreal biome northward through expansion at the northern and contraction at the southern boundary respectively. However, biome-scale evidence of such a shift is rare. Here, we used remotely-sensed tree cover data to quantify temporal changes across the North American boreal biome from 2000 to 2019. We reveal a strong north-south asymmetry in tree cover change, coupled with a range shrinkage of tree cover distributions. We found no evidence for tree cover expansion in the northern biome, while tree cover increased markedly in the core of the biome range. By contrast, tree cover declined along the southern biome boundary, where losses were related largely to wildfires and timber logging.  We show that these contrasting trends are structural indicators for a possible onset of a biome contraction which may lead to long-term carbon declines. 
# 
# 



#### Read libraries and working directory####
library(tidyverse)
library(plyr)
library(corrplot)
library(gam)
library(mgcv)
library(splitstackshape)
library(mgcViz)
library(nlme)
library(PerformanceAnalytics)
library(scales)
library(reshape2)
library(ggpubr)
library(metR)
library(gridExtra)
library(png)
library(cowplot)
library(ggplotify)
library(ggfortify)
library(ggimage)


setwd()


#### Set up functions to be used####
#function to calcualte moving averages within a time series
ma = function(x, n = 5){
  stats::filter(x, rep(1 / n, n), method = "convolution", sides = 2)
}




#### Read in dataset ####
all.boreal = readxl::read_xlsx("SI_data1.xlsx")

str(all.boreal)

all.boreal$GLC = as.factor(all.boreal$GLC) #land cover type as factor

all.boreal$GLC = factor(all.boreal$GLC, levels = c("Non-woody", "Shrubs", "Needleleaf forest", "Mixed forest", "Broadleaf forest", "Unknown forest"))

all.boreal$ecozone = as.factor(all.boreal$ecozone) 
all.boreal$dist_type3 = as.factor(all.boreal$dist_type3)
all.boreal$dist_type3 = factor(all.boreal$dist_type3, levels = c("Disturbed <1985/Undisturbed", "Fire 1985-1999", "Harvest 1985-1999", "Fire 2000-2019", "Harvest 2000-2019"))


#round distance to 0.025
all.boreal$rounddistance = round(all.boreal$distance/0.025)*0.025


#Boreal distribution boundaries as factors
all.boreal$flon = as.factor(all.boreal$lon) #transect longitudes as fact


#Include broader classification of disturbance without disturbance timing
all.boreal$dist_type = "Undisturbed"
all.boreal$dist_type[grepl("Fire", all.boreal$dist_type3, fixed = T)] = "Fire"
all.boreal$dist_type[grepl("Harvest", all.boreal$dist_type3, fixed = T)] = "Harvest"

all.boreal$dist_type = as.factor(all.boreal$dist_type)
all.boreal$dist_type = factor(all.boreal$dist_type, levels = c("undisturbed", "Fire", "Harvest"))


#### Data exploration: Correlation between predictors ####
variables = dplyr::select(all.boreal, c(distance,elev,meantc,temp,meanpre,tchange8019,pchange8019,tchange0019,pchange0019))

var.matrix = cor(variables, use = "pairwise.complete.obs")

corrplot.all = corrplot(var.matrix,
                        type = "lower",
                        method = "number",
                        tl.pos = "lt",
                        tl.cex = 1.5,
                        tl.col = "black",
                        tl.offset = 0.5,
                        tl.srt = 45,
                        cl.pos = F,
                        diag = F)




env_names = c("SBD", "Elevation", "Mean TC","MAT", "MAP", "trend MAT80-19", "trend MAP80-19", "trend MAT00-19",  "trend MAP00-19")

variables$trend.class = "1"
variables$trend.class[all.boreal$trend.yp < 0 & all.boreal$sig.yp <= 0.05] = "3"
variables$trend.class[all.boreal$trend.yp > 0 & all.boreal$sig.yp <= 0.05] = "2"

#variables$trend.yp = variables$trend.yp*100/boreal$meantc

pca = prcomp(dplyr::select(variables, -trend.class), scale. = T)


pcaplot = autoplot(pca,
                   data = variables,
                   loadings = T,
                   loadings.label = T,
                   loadings.label.size = 7,
                   loadings.label.repel = T,
                   loadings.label.colour = "black",
                   loadings.colour = "black",
                   loadings.label.label = env_names,
                   colour = "trend.class",
                   size = 3,
                   alpha = .2)+
  stat_ellipse(aes(col = trend.class), lwd = 1, type = "norm")+
  geom_hline(yintercept = 0, lwd = .4, lty = "dashed")+
  geom_vline(xintercept = 0, lwd = .4, lty = "dashed")+
  scale_color_manual(name = "",
                     values = c("darkgrey", "deepskyblue4", "tomato4"),
                     label = c("No Trend", "Gain", "Loss"))+
  theme_classic()+
  theme(axis.line = element_line(colour = "black"),
        text = element_text(size = 25),
        legend.text = element_text(size = 25),)

pcaplot 




#### 4.1 Generalised additive mixed-effects models####
#Tree cover change and boundary distance by disturbance: As tree cover changes (trend.yp) follow a clear non-linear pattern along the south-north transects, we perform statistical analyses using a generalised additive model. We only include the standardised distance from the boreal boundary as covariate in the model and let it depend on disturbance type. 

#We first check for a potential effect of transects (expressed as unique longitudes) on tree cover changes. .
ggplot(data = all.boreal, aes(x = flon,
                          y = trend.yp)) +
  geom_boxplot(fill = "grey70")+ 
  geom_hline(yintercept = 0, 
             lty = "dashed",
             lwd = 1)+
  facet_wrap(.~dist_type3)+
  xlab("Transect (expressed as longitudes)")+ 
  ylab("Tree cover change")+
  theme_classic()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(size = 10, 
                                   angle = 45,
                                   vjust = .5))

#Some differences in tree cover change between transects exist. We therefore include transect as random effect in a generalised mixed-effects model. 

#GAMM model fit

gamm.distance = gam(data = all.boreal, 
                    trend.yp ~ 1 + s(distance, by = dist_type3) +
                      s(flon, bs = "re"),
                    method = "REML")

summary(gamm.distance)

#draw distance effect in the model
draw(gamm.distance, select = 1)

#extract standard deviations
gam.vcomp(gamm.distance) 

#calculate intra-class correlation ICC using mixed effects variance and overall variance
ICC <- 0.34676822  ^2 / (0.34676822  ^2 + 0.08852528 ^2)
ICC

# THE ICC is almost 1 which indicates low mixed effects per transect. However, we decided to include the mixed effect in the model fit nonetheless. 



#### 4.1.1 Model validation####

#' Get residuals and fitted values:
all.boreal$E1  <- resid(gamm.distance)
all.boreal $mu1 <- fitted(gamm.distance)



#' Pearson residuals versus fitted values.
p1 <- ggplot(all.boreal, aes(x = mu1, y = E1))+
  geom_point(size = 1)+
  geom_hline(yintercept = 0, 
             linetype="dashed",
             color = "red")+
  labs(x = "Fitted values", y = "Residuals")+
  theme(text = element_text(size = 10))

#' Pearson residuals versus disturbance type
p2 <- ggplot(all.boreal, aes(x = dist_type3, y = E1))+
  geom_boxplot()+
  geom_hline(yintercept = 0, color = "red")+
  labs(x = "Disturbance", y = "Residuals")+
  theme(text = element_text(size=15),
        axis.text.x = element_text(color = "#993333", 
                                   size = 6, 
                                   angle = 45))


#' Fitted values versus observed data.
p3 <- ggplot(all.boreal, aes(x = mu1, y = trend.yp))+
  geom_point(size = 0.5, col= grey(0.5))+
  geom_abline(intercept = 0, slope = 1, color = "red")+
  labs(x = "Fitted values", y = "Observed Tree cover change")+
  theme(text = element_text(size = 10))

#' Residuals versus distance
p4 <- ggplot(all.boreal, aes(x = distance, y = E1))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept = 0, color = "red")+
  labs(x = "Distance", y = "Residuals")+
  theme(text = element_text(size = 10))

#' Residuals versus longitude
p5 <- ggplot(all.boreal, aes(x = flon, y = E1))+
  geom_boxplot()+
  geom_hline(yintercept = 0, color = "red")+
  labs(x = "Longitude", y = "Residuals")+
  theme(text = element_text(size=15),
        axis.text.x = element_text(color = "#993333", 
                                   size = 6, 
                                   angle = 45))

plot_grid(p1, p4, 
          p5, p2, p3,  ncol = 2, 
          labels = c('A', 'B','C', 'D'), 
          rel_widths = c(1, 1))

#Model fit seems satisfying. No patterns in residuals. 


#### 4.1.2 Check for spatial correlation of model residuals####
#plot variogram and fitted function
#perform variogram using coordinates for spatial distance calculation
v1 = variogram(object =  E1 ~ x+y,
               locations = ~x+y,
               data = all.boreal,
               cutoff = 250,
               width = 250/100)
#cutoff is the maximum distance to be considered, width the number of points
v1

v1.fit = fit.variogram(v1,vgm(psill = 0.1, model = "Exp", range = 250, nugget = 1))  #fit a function to the variogram. The model choice should be informed by visual fit. The available models can be shown by command "vgm()". range and nugget values are start parameter values and should not be too far from the realistic values. Otherwise the model will not converge.

v1.line = variogramLine(v1.fit, maxdist = 250, n = 100) #fit a line to the model which can then be used in ggplot.


ggplot(v1, aes(dist, gamma))+
  geom_line(data = v1.line, size = 2)+
  geom_point(size = 3, col = "deepskyblue4")+
  geom_vline(xintercept = v1.fit$range[2]*3, lty = "dashed", size = 1)+
  geom_hline(yintercept = c(v1.fit$psill[1], v1.fit$psill[1] + v1.fit$psill[2]), lty = "dashed", size = 1)+
  labs(x = "Lag distance (km)", y = "Variance")+
  theme_classic()+
  theme(text = element_text(size = 30),
        plot.margin = unit(c(0,3,0,0), "line"),
        legend.position = "right")


ggsave(filename = "variogram_gam_tc_distance.png", plot = last_plot(), path = "C:/Users/rotba001/OneDrive - WageningenUR/PhD_project/Scientific_data/Figures/Tree_cover_final/whole_transects", width = 14, height = 10)

#Remaining spatial correlation between sample plots along transects with a lag distance of around 100km. We will use an additive mixed-effects model with a spatial correlation structure to account for this remaining spatial correlation. We will follow this model approach for each predictor separately below. 



#### 4.2 GAMM with spatial correlation structure: model fitting of relative and absolute tree cover trends and environmental variables #### 
#### 4.2.1 Model RELATIVE TREE COVER and.... ####
#...distance by disturbance and vegetation type
all.boreal$dist_GLC = interaction(all.boreal$dist_type3, all.boreal$GLC)



# Model expressions for...

#... vegetation type and boundary distance
M01 = gamm(reltrend ~ 
             s(distance, by = dist_GLC),
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)

sumM01 = summary(M01$gam)
sumM01



#... standard boundary distance
M1 = gamm(reltrend ~ 
            s(distance, by = dist_type3), # smooth term for predictor
          
          random = list(lon =~1), # random effect of transect
          correlation = corExp(1, form = ~ x + y),  #spatial correlation structure
          data = all.boreal)



sumM1 = summary(M1$gam)
sumM1

#...mean temperature 1980-2019
M2 = gamm(reltrend ~ 
            s(temp, by = dist_type3),
          
          random = list(lon =~1),
          correlation = corExp(1, form = ~ x + y), 
          data = all.boreal)


sumM2 = summary(M2$gam)
sumM2

#...mean annual precipitation 1980-2019
M3 = gamm(reltrend ~ 
            s(meanpre, by = dist_type3),
          
          random = list(lon =~1),
          correlation = corExp(1, form = ~ x + y), 
          data = all.boreal)

sumM3 = summary(M3$gam)
sumM3

#...mean tree cover 2000-2019
M4 = gamm(reltrend ~ 
            s(meantc, by = dist_type3),
          
          random = list(lon =~1),
          correlation = corExp(1, form = ~ x + y), 
          data = all.boreal)

sumM4 = summary(M4$gam)
sumM4

#...elevation
M5 = gamm(reltrend ~ 
            s(elev, by = dist_type3),
          
          random = list(lon =~1),
          correlation = corExp(1, form = ~ x + y), 
          data = all.boreal)

sumM5 = summary(M5$gam)
sumM5

#...trend in annual temperature 1980-2019
M6 = gamm(reltrend ~ 
            s(tchange8019, by = dist_type3),
          
          random = list(lon =~1),
          correlation = corExp(1, form = ~ x + y), 
          data = all.boreal)

sumM6 = summary(M6$gam)
sumM6

#...trend in annual precipation 1980-2019
M7 = gamm(reltrend ~ 
            s(pchange8019, by = dist_type3),
          
          random = list(lon =~1),
          correlation = corExp(1, form = ~ x + y), 
          data = all.boreal)

sumM7 = summary(M7$gam)
sumM7

#...trend in annual precipiation 2000-2019
M8 = gamm(reltrend ~ 
            s(pchange0019, by = dist_type3),
          
          random = list(lon =~1),
          correlation = corExp(1, form = ~ x + y), 
          data = all.boreal)

sumM8 = summary(M8$gam)
sumM8

#...trend in annual temperature 2000-2019
M9 = gamm(reltrend ~ 
            s(tchange0019, by = dist_type3),
          
          random = list(lon =~1),
          correlation = corExp(1, form = ~ x + y), 
          data = all.boreal)

sumM9 = summary(M9$gam)
sumM9


#### 4.2.2 Model validation####
#' Get residuals and fitted values:
all.boreal$E01  <- resid(M01$gam)
all.boreal$mu01 <- fitted(M01$gam)
all.boreal$E1  <- resid(M1$gam)
all.boreal$mu1 <- fitted(M1$gam)
all.boreal$E2  <- resid(M2$gam)
all.boreal$mu2 <- fitted(M2$gam)
all.boreal$E3  <- resid(M3$gam)
all.boreal$mu3 <- fitted(M3$gam)
all.boreal$E4  <- resid(M4$gam)
all.boreal$mu4 <- fitted(M4$gam)
all.boreal$E5  <- resid(M5$gam)
all.boreal$mu5 <- fitted(M5$gam)
all.boreal$E6  <- resid(M6$gam)
all.boreal$mu6 <- fitted(M6$gam)
all.boreal$E7  <- resid(M7$gam)
all.boreal$mu7 <- fitted(M7$gam)
all.boreal$E8  <- resid(M8$gam)
all.boreal$mu8 <- fitted(M8$gam)
all.boreal$E9  <- resid(M9$gam)
all.boreal$mu9 <- fitted(M9$gam)


#The code below can be used for validation of each model by changing the fitted and residual values, mux and Ex, as well as the predictor used in the model. 
#' Pearson residuals versus fitted values.

p1 <- ggplot(all.boreal, aes(x = mu1, y = E1))+
  geom_point(size = 1)+
  geom_hline(yintercept = 0, 
             linetype="dashed",
             color = "red")+
  labs(x = "Fitted values", y = "Residuals")+
  theme(text = element_text(size = 10))

#' Pearson residuals versus disturbance type
p2 <- ggplot(all.boreal, aes(x = dist_type3, y = E1))+
  geom_boxplot()+
  geom_hline(yintercept = 0, color = "red")+
  labs(x = "Disturbance", y = "Residuals")+
  theme(text = element_text(size=15),
        axis.text.x = element_text(color = "#993333", 
                                   size = 6, 
                                   angle = 45))


#' Fitted values versus observed data.
p3 <- ggplot(all.boreal, aes(x = mu1, y = reltrend))+
  geom_point(size = 0.5, col= grey(0.5))+
  geom_abline(intercept = 0, slope = 1, color = "red")+
  labs(x = "Fitted values", y = "Observed Tree cover change")+
  theme(text = element_text(size = 10))

#' Residuals versus predictor in the model
p4 <- ggplot(all.boreal, aes(x = distance, y = E1))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept = 0, color = "red")+
  labs(y = "Residuals")+
  theme(text = element_text(size = 10))

#' Residuals versus longitude
p5 <- ggplot(all.boreal, aes(x = flon, y = E1))+
  geom_boxplot()+
  geom_hline(yintercept = 0, color = "red")+
  labs(x = "Longitude", y = "Residuals")+
  theme(text = element_text(size=15),
        axis.text.x = element_text(color = "#993333", 
                                   size = 6, 
                                   angle = 45))

plot_grid(p1, p4, 
          p5, p2, p3,  ncol = 2, 
          labels = c('A', 'B','C', 'D', 'E'), 
          rel_widths = c(1, 1))

#All model fits seem satisfying. No concerning patterns in residuals. 



#### 4.2.3 Model prediction and visualisation####
#Create table with min and max values for each predictor variables in the models and for each disturbance type
glcminmax = ddply(all.boreal, c("GLC", "dist_type3"), summarise,
                  mindist = min(distance),
                  maxdist = max(distance))


climminmax = ddply(all.boreal, "dist_type3", summarise,
                   mint = min(temp),
                   maxt = max(temp),
                   minp = min(meanpre),
                   maxp = max(meanpre),
                   mintc = min(meantc),
                   maxtc = max(meantc),
                   mindist = min(distance),
                   maxdist = max(distance),
                   minelev = min(elev),
                   maxelev = max(elev),
                   mintt = min(tchange8019),
                   maxtt = max(tchange8019),
                   mintp = min(pchange8019),
                   maxtp = max(pchange8019),
                   mintt0019 = min(tchange0019),
                   maxtt0019 = max(tchange0019),
                   mintp0019 = min(pchange0019),
                   maxtp0019 = max(pchange0019))



tempdist = list()
predist = list()
tcdist = list()
ddist = list()
elevdist = list()
ttdist = list()
tpdist = list()
tt0019dist = list()
tp0019dist = list()

distglc = list()


#for each disturbance type create new data as a sequence between min and max predictor values
for (i in 1:length(levels(all.boreal$dist_type3))) {
  print(paste("Loop", i, "/", length(levels(all.boreal$dist_type3))))
  t = seq(climminmax$mint[i], climminmax$maxt[i], length.out = 100)
  tempdist[[i]] = t
  p = seq(climminmax$minp[i], climminmax$maxp[i], length.out = 100)
  predist[[i]] = p
  tc = seq(climminmax$mintc[i], climminmax$maxtc[i], length.out = 100)
  tcdist[[i]] = tc
  dist = seq(climminmax$mindist[i], climminmax$maxdist[i], length.out = 100)
  ddist[[i]] = dist
  elev = seq(climminmax$minelev[i], climminmax$maxelev[i], length.out = 100)
  elevdist[[i]] = elev
  ttemp.yp = seq(climminmax$mintt[i], climminmax$maxtt[i], length.out = 100)
  ttdist[[i]] = ttemp.yp
  tpre.yp = seq(climminmax$mintp[i], climminmax$maxtp[i], length.out = 100)
  tpdist[[i]] = tpre.yp
  ttemp.yp0019 = seq(climminmax$mintt0019[i], climminmax$maxtt0019[i], length.out = 100)
  tt0019dist[[i]] = ttemp.yp0019
  tpre.yp0019 = seq(climminmax$mintp0019[i], climminmax$maxtp0019[i], length.out = 100)
  tp0019dist[[i]] = tpre.yp0019
}


for (i in 1:length(levels(all.boreal$dist_GLC))) {
  distglc[[i]] = seq(glcminmax$mindist[i], glcminmax$maxdist[i], length.out = 100)
}



# combine all values into one dummy table
dummy = data.frame(temp = unlist(tempdist),
                   meanpre = unlist(predist),
                   meantc = unlist(tcdist),
                   distance = unlist(ddist),
                   elev = unlist(elevdist),
                   tchange8019 = unlist(ttdist),
                   pchange8019 = unlist(tpdist),
                   tchange0019 = unlist(tt0019dist),
                   pchange0019 = unlist(tp0019dist),
                   dist_type3 = as.factor(rep(c(levels(all.boreal$dist_type3)), each = 100)))

dummy$dist_type3 = factor(dummy$dist_type3, levels = levels(all.boreal$dist_type3))

dummy$dist_label = dummy$dist_type3
levels(dummy$dist_label) = c("Disturbed <1985/Undisturbed", "Fire 1985-1999", "Harvest 1985-1999", "Fire 2000-2019", "Harvest 2000-2019")


dummyglc = data.frame(distance = unlist(distglc),
                      GLC = as.factor(rep(levels(all.boreal$GLC), each = 500)),
                      dist_type3 = as.factor(rep(c(levels(all.boreal$dist_type3)), each = 100)))

dummyglc$GLC = factor(dummyglc$GLC, levels = levels(all.boreal$GLC))
dummyglc$dist_type3 = factor(dummyglc$dist_type3, levels = levels(all.boreal$dist_type3))
dummyglc$dist_GLC = interaction(dummyglc$dist_type3, dummyglc$GLC)
dummyglc$dist_GLC = factor(dummyglc$dist_GLC, levels = levels(all.boreal$dist_GLC))




#From dummy data predict and visualise...
#distance by vegetation type and disturbance
dummyfit = predict.gam(object = M01$gam, newdata = dummyglc, se.fit = T)

dummyglc$fdistrel = dummyfit$fit
dummyglc$sdistrel = dummyfit$se.fit

dummyglc$dist_type = "Undisturbed"
dummyglc$dist_type[grepl("Fire", dummyglc$dist_type3, fixed = T)] = "Fire"
dummyglc$dist_type[grepl("Harvest", dummyglc$dist_type3, fixed = T)] = "Harvest"
dummyglc$dist_type = as.factor(dummyglc$dist_type)
dummyglc$dist_type = factor(dummyglc$dist_type, levels = levels(all.boreal$dist_type))

distcol = c("#018571", "#FDBE85", "#FD8D3C",  "#DFC27D", "#A6611A")

textdf = data.frame(x = c(1.1,0.1,-0.9),
                    y = -5.5, 
                    text = c("North", "Interior", "South"))


glcrelplot = 
  ggplot(dummyglc, aes(distance, fdistrel))+
  geom_point(data = all.boreal, aes(distance, reltrend, col = dist_type3), pch = 1, alpha = .4)+
  geom_line(aes(col = dist_type3), lwd = 1.3)+
  # geom_line(aes(distance, fdist + sdist), lwd = .5)+
  # geom_line(aes(distance, fdist - sdist), lwd = .5)+
  # 
  scale_color_manual(name = "", 
                     values = distcol,
                     labels = levels(dummyglc$dist_type3))+
  
  scale_x_continuous(limits = c(-1.6, 1.6), breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  
  guides(col = guide_legend(nrow = 2))+
  
  geom_text(data = textdf, aes(x, y, label = text),
            hjust = 0,
            vjust = 0,
            size = 9)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = c(-1, 0, 1), lty = "dashed", lwd = 1)+
  
  
  coord_flip()+
  
  facet_grid(GLC~dist_type,
             scale = "fixed")+
  
  labs(x = "Standardised boundary distance",
       y = expression(paste('Relative tree cover change (% ' ~'%'^-1 ~year^-1, ')')),
       tag = "B")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 25),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "")
glcrelplot





#Distance 
dummyfit = predict.gam(object = M1$gam, newdata = dummy, se.fit = T)

dummy$fdist = dummyfit$fit
dummy$sdist = dummyfit$se.fit

textdf = data.frame(x = c(1.1,0.1,-0.9),
                    y = -5, 
                    text = c("North", "Interior", "South"))



distplot = 
  ggplot(dummy, aes(distance, fdist))+
  geom_point(data = all.boreal, aes(distance, reltrend), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(distance, fdist + sdist), lwd = .5)+
  geom_line(aes(distance, fdist - sdist), lwd = .5)+
  
  scale_x_continuous(limits = c(-1.6, 1.6), breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  
  
  
  geom_text(data = textdf, aes(x, y, label = text),
            hjust = 0,
            size = 9)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = c(-1, 0, 1), lty = "dashed", lwd = 1)+
  
  
  coord_flip()+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Stand. boundary distance",
       y = expression(paste('TC change (% ' ~"%"^-1 ~year^-1, ')')),
       tag = "A")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
distplot





#Temperature 
dummyfit = predict.gam(object = M2$gam, newdata = dummy, se.fit = T)

dummy$ftemp = dummyfit$fit

dummy$stemp = dummyfit$se.fit



matplot = 
  ggplot(dummy, aes(temp, ftemp))+
  geom_point(data = all.boreal, aes(temp, reltrend), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(temp, ftemp + stemp), lwd = .5)+
  geom_line(aes(temp, ftemp - stemp), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Mean annual temperature (°C)",
       y = expression(paste('TC change (% ' ~"%"^-1 ~year^-1, ')')),
       tag = "D")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
matplot



#precipitation
dummyfit = predict.gam(object = M3$gam, newdata = dummy, se.fit = T)

dummy$fpre = dummyfit$fit
dummy$spre = dummyfit$se.fit



mapplot = 
  ggplot(dummy, aes(meanpre, fpre))+
  geom_point(data = all.boreal, aes(meanpre, reltrend), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(meanpre, fpre + spre), lwd = .5)+
  geom_line(aes(meanpre, fpre - spre), lwd = .5)+
  
  scale_x_continuous(limits = c(0,2000))+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Mean annual precipitation (mm)",
       y = expression(paste('TC change (% ' ~"%"^-1 ~year^-1, ')')),
       tag = "E")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        axis.text.x = element_text(size = 10),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
mapplot






#Mean tree cover
dummyfit = predict.gam(object = M4$gam, newdata = dummy, se.fit = T)

dummy$ftc = dummyfit$fit
dummy$stc = dummyfit$se.fit



tcplot = 
  ggplot(dummy, aes(meantc, ftc))+
  geom_point(data = all.boreal, aes(meantc, reltrend), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(meantc, ftc + stc), lwd = .5)+
  geom_line(aes(meantc, ftc - stc), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Mean tree cover (%)",
       y = expression(paste('TC change (% ' ~"%"^-1 ~year^-1, ')')),
       tag = "B")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
tcplot




#Elevation
dummyfit = predict.gam(object = M5$gam, newdata = dummy, se.fit = T)

dummy$felev = dummyfit$fit
dummy$selev = dummyfit$se.fit



elevplot.rel = 
  ggplot(dummy, aes(elev, felev))+
  geom_point(data = all.boreal, aes(elev, reltrend), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(elev, felev + selev), lwd = .5)+
  geom_line(aes(elev, felev - selev), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Elevation (m)",
       y = expression(paste('TC change (% ' ~"%"^-1 ~year^-1, ')')),
       tag = "C")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        axis.text.x = element_text(size = 10),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
elevplot.rel

#Temperature change 1980-2019
dummyfit = predict.gam(object = M6$gam, newdata = dummy, se.fit = T)

dummy$ftt = dummyfit$fit
dummy$stt = dummyfit$se.fit



ttplot = 
  ggplot(dummy, aes(tchange8019, ftt))+
  geom_point(data = all.boreal, aes(tchange8019, reltrend), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(tchange8019, ftt + stt), lwd = .5)+
  geom_line(aes(tchange8019, ftt - stt), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = expression(paste('Mean annual temperature change 1980-2019 (°C ' ~year^-1, ')')),
       y = expression(paste('TC change (% ' ~"%"^-1 ~year^-1, ')')),
       tag = "F")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
ttplot


#Precipitation change 1980-2019
dummyfit = predict.gam(object = M7$gam, newdata = dummy, se.fit = T)

dummy$ftp = dummyfit$fit
dummy$stp = dummyfit$se.fit



tpplot = 
  ggplot(dummy, aes(pchange8019, ftp))+
  geom_point(data = all.boreal, aes(pchange8019, reltrend), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(pchange8019, ftp + stp), lwd = .5)+
  geom_line(aes(pchange8019, ftp - stp), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = expression(paste('Mean annual precipitation change 1980-2019 (mm ' ~year^-1, ')')),
       y = expression(paste('TC change (% ' ~"%"^-1 ~year^-1, ')')),
       tag = "G")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
tpplot


#Temperature change 2000-2019
dummyfit = predict.gam(object = M9$gam, newdata = dummy, se.fit = T)

dummy$ftt0019 = dummyfit$fit
dummy$stt0019 = dummyfit$se.fit



tt0019plot = 
  ggplot(dummy, aes(tchange0019, ftt0019))+
  geom_point(data = all.boreal, aes(tchange0019, reltrend), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(tchange0019, ftt0019 + stt0019), lwd = .5)+
  geom_line(aes(tchange0019, ftt0019 - stt0019), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = expression(paste('Mean annual temperature change 2000-2019 (°C ' ~year^-1, ')')),
       y = expression(paste('TC change (% ' ~"%"^-1 ~year^-1, ')')),
       tag = "H")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
tt0019plot


#Precipitation change 2000-2019
dummyfit = predict.gam(object = M8$gam, newdata = dummy, se.fit = T)

dummy$ftp0019 = dummyfit$fit
dummy$stp0019 = dummyfit$se.fit



tp0019plot = 
  ggplot(dummy, aes(pchange0019, ftp0019))+
  geom_point(data = all.boreal, aes(pchange0019, reltrend), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(pchange0019, ftp0019 + stp0019), lwd = .5)+
  geom_line(aes(pchange0019, ftp0019 - stp0019), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = expression(paste('Mean annual precipitation change 2000-2019 (mm ' ~year^-1, ')')),
       y = expression(paste('TC change (% ' ~"%"^-1 ~year^-1, ')')),
       tag = "I")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
tp0019plot



  #### 4.2.4 Model ABSOLUTE TREE COVER and.... ####
#...distance by disturbance and vegetation type
M02 = gamm(trend.yp ~ 
             s(distance, by = dist_GLC),
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)

sumM02 = summary(M02$gam)
sumM02

saveRDS(M02, "model02.rds")


#...distance 
M11 = gamm(trend.yp ~ 
             s(distance, by = dist_type3),
           
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)



sumM11 = summary(M11$gam)
sumM11

#...Mean annual temperature 1980-2019
M12 = gamm(trend.yp ~ 
             s(temp, by = dist_type3),
           
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)


sumM12 = summary(M12$gam)
sumM12

#...Mean annual precipitation 1980-2019
M13 = gamm(trend.yp ~ 
             s(meanpre, by = dist_type3),
           
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)

sumM13 = summary(M13$gam)
sumM13

#...Meant tree cover 2000-2019
M14 = gamm(trend.yp ~ 
             s(meantc, by = dist_type3),
           
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)

sumM14 = summary(M4$gam)
sumM14

#...elevation
M15 = gamm(trend.yp ~ 
             s(elev, by = dist_type3),
           
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)

sumM15 = summary(M15$gam)
sumM15

#...Trend in annual temperature 1980-2019
M16 = gamm(trend.yp ~ 
             s(tchange8019, by = dist_type3),
           
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)

sumM16 = summary(M16$gam)
sumM16

#...Trend in annual precipiation 1980-2019
M17 = gamm(trend.yp ~ 
             s(pchange8019, by = dist_type3),
           
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)

sumM17 = summary(M17$gam)
sumM17

#...Trend in annual precipiation 2000-2019
M18 = gamm(trend.yp ~ 
             s(pchange0019, by = dist_type3),
           
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)

sumM18 = summary(M18$gam)
sumM18

#...Trend in annual temperature 2000-2019
M19 = gamm(trend.yp ~ 
             s(tchange0019, by = dist_type3),
           
           random = list(lon =~1),
           correlation = corExp(1, form = ~ x + y), 
           data = all.boreal)

sumM19 = summary(M19$gam)
sumM19



#### 4.2.5 Model validation####
#' Get residuals and fitted values:
all.boreal$E02  <- resid(M02$gam)
all.boreal$mu02 <- fitted(M02$gam)
all.boreal$E11  <- resid(M11$gam)
all.boreal$mu11 <- fitted(M11$gam)
all.boreal$E12  <- resid(M12$gam)
all.boreal$mu12 <- fitted(M12$gam)
all.boreal$E13  <- resid(M13$gam)
all.boreal$mu13 <- fitted(M13$gam)
all.boreal$E14  <- resid(M14$gam)
all.boreal$mu14 <- fitted(M14$gam)
all.boreal$E15  <- resid(M15$gam)
all.boreal$mu15 <- fitted(M15$gam)
all.boreal$E16  <- resid(M16$gam)
all.boreal$mu16 <- fitted(M16$gam)
all.boreal$E17  <- resid(M17$gam)
all.boreal$mu17 <- fitted(M17$gam)
all.boreal$E18  <- resid(M18$gam)
all.boreal$mu18 <- fitted(M18$gam)
all.boreal$E19  <- resid(M19$gam)
all.boreal$mu19 <- fitted(M19$gam)


#The code below can be used for validation of each model by changing the fitted and residual values, mux and Ex, as well as the predictor used in the model. 
#' Pearson residuals versus fitted values.

p1 <- ggplot(all.boreal, aes(x = mu11, y = E11))+
  geom_point(size = 1)+
  geom_hline(yintercept = 0, 
             linetype="dashed",
             color = "red")+
  labs(x = "Fitted values", y = "Residuals")+
  theme(text = element_text(size = 10))

#' Pearson residuals versus disturbance type
p2 <- ggplot(all.boreal, aes(x = dist_type3, y = E11))+
  geom_boxplot()+
  geom_hline(yintercept = 0, color = "red")+
  labs(x = "Disturbance", y = "Residuals")+
  theme(text = element_text(size=15),
        axis.text.x = element_text(color = "#993333", 
                                   size = 6, 
                                   angle = 45))


#' Fitted values versus observed data.
p3 <- ggplot(all.boreal, aes(x = mu11, y = trend.yp))+
  geom_point(size = 0.5, col= grey(0.5))+
  geom_abline(intercept = 0, slope = 1, color = "red")+
  labs(x = "Fitted values", y = "Observed Tree cover change")+
  theme(text = element_text(size = 10))

#' Residuals versus predictor in the model
p4 <- ggplot(all.boreal, aes(x = distance, y = E11))+
  geom_point()+
  geom_smooth()+
  geom_hline(yintercept = 0, color = "red")+
  labs(y = "Residuals")+
  theme(text = element_text(size = 10))

#' Residuals versus longitude
p5 <- ggplot(all.boreal, aes(x = flon, y = E11))+
  geom_boxplot()+
  geom_hline(yintercept = 0, color = "red")+
  labs(x = "Longitude", y = "Residuals")+
  theme(text = element_text(size=15),
        axis.text.x = element_text(color = "#993333", 
                                   size = 6, 
                                   angle = 45))

plot_grid(p1, p4, 
          p5, p2, p3,  ncol = 2, 
          labels = c('A', 'B','C', 'D', 'E'), 
          rel_widths = c(1, 1))

#All model fits seem satisfying. No concerning patterns in residuals. 



#### 4.2.6 Model prediction and visualisation####
#From dummy data predict and visualise...
#distance by vegetation type and disturbance
dummyglcfit = predict.gam(object = M02$gam, newdata = dummyglc, se.fit = T)

dummyglc$fdist = dummyglcfit$fit
dummyglc$sdist = dummyglcfit$se.fit


textdf = data.frame(x = c(1.1,0.1,-0.9),
                    y = -2, 
                    text = c("North", "Interior", "South"))

distcol = c("#018571", "#FDBE85", "#FD8D3C",  "#DFC27D", "#A6611A")

glcabsplot = 
  ggplot(dummyglc, aes(distance, fdist))+
  geom_point(data = all.boreal, aes(distance, trend.yp, col = dist_type3), pch = 1, alpha = .4)+
  geom_line(aes(col = dist_type3), lwd = 1.3)+
  # geom_line(aes(distance, fdist + sdist), lwd = .5)+
  # geom_line(aes(distance, fdist - sdist), lwd = .5)+
  # 
  scale_color_manual(name = "",
                     values = distcol, 
                     labels = levels(dummyglc$dist_type3))+
  
  scale_x_continuous(limits = c(-1.6, 1.6), breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  
  geom_text(data = textdf, aes(x, y, label = text),
            hjust = 0,
            size = 9)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = c(-1, 0, 1), lty = "dashed", lwd = 1)+
  
  
  coord_flip()+
  
  facet_grid(GLC~dist_type, 
             scale = "fixed")+
  
  labs(x = "Standardised boundary distance",
       y = expression(paste('Tree cover change (% ' ~year^-1, ')')),
       tag = "A")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 25),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "")
glcabsplot




#Distance 
dummyfit = predict.gam(object = M11$gam, newdata = dummy, se.fit = T)

dummy$fdist = dummyfit$fit
dummy$sdist = dummyfit$se.fit

textdf = data.frame(x = c(1.1,0.1,-0.9),
                    y = -1.5, 
                    text = c("North", "Interior", "South"))



distplot = 
  ggplot(dummy, aes(distance, fdist))+
  geom_point(data = all.boreal, aes(distance, trend.yp), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(distance, fdist + sdist), lwd = .5)+
  geom_line(aes(distance, fdist - sdist), lwd = .5)+
  
  scale_x_continuous(limits = c(-1.6, 1.6), breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))+
  
  
  
  geom_text(data = textdf, aes(x, y, label = text),
            hjust = 0,
            size = 9)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = c(-1, 0, 1), lty = "dashed", lwd = 1)+
  
  
  coord_flip()+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Stand. boundary distance",
       y = expression(paste('Tree cover change (% ' ~year^-1, ')')),
       tag = "A")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
distplot





#Mean annual temperatures 1980-2019 
dummyfit = predict.gam(object = M12$gam, newdata = dummy, se.fit = T)

dummy$ftemp = dummyfit$fit

dummy$stemp = dummyfit$se.fit



matplot = 
  ggplot(dummy, aes(temp, ftemp))+
  geom_point(data = all.boreal, aes(temp, trend.yp), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(temp, ftemp + stemp), lwd = .5)+
  geom_line(aes(temp, ftemp - stemp), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Mean annual temperature (°C)",
       y = expression(paste('TC change (% ' ~year^-1, ')')),
       tag = "D")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
matplot



#Mean annual precipitation 1980-2019
dummyfit = predict.gam(object = M13$gam, newdata = dummy, se.fit = T)

dummy$fpre = dummyfit$fit
dummy$spre = dummyfit$se.fit



mapplot = 
  ggplot(dummy, aes(meanpre, fpre))+
  geom_point(data = all.boreal, aes(meanpre, trend.yp), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(meanpre, fpre + spre), lwd = .5)+
  geom_line(aes(meanpre, fpre - spre), lwd = .5)+
  
  scale_x_continuous(limits = c(0,1700))+
  scale_y_continuous(limits = c(-2.5,2.5))+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Mean annual precipitation (mm)",
       y = expression(paste('TC change (% ' ~year^-1, ')')),
       tag = "E")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
mapplot



#Mean tree cover 2000-2019
dummyfit = predict.gam(object = M14$gam, newdata = dummy, se.fit = T)

dummy$ftc = dummyfit$fit
dummy$stc = dummyfit$se.fit



tcplot = 
  ggplot(dummy, aes(meantc, ftc))+
  geom_point(data = all.boreal, aes(meantc, trend.yp), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(meantc, ftc + stc), lwd = .5)+
  geom_line(aes(meantc, ftc - stc), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Mean tree cover (%)",
       y = expression(paste('TC change (% ' ~year^-1, ')')),
       tag = "B")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
tcplot




#Elevation
dummyfit = predict.gam(object = M15$gam, newdata = dummy, se.fit = T)

dummy$felev = dummyfit$fit
dummy$selev = dummyfit$se.fit



elevplot = 
  ggplot(dummy, aes(elev, felev))+
  geom_point(data = all.boreal, aes(elev, trend.yp), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(elev, felev + selev), lwd = .5)+
  geom_line(aes(elev, felev - selev), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = "Elevation (m)",
       y = expression(paste('TC change (% ' ~year^-1, ')')),
       tag = "C")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
elevplot





#Temperature change 1980-2019
dummyfit = predict.gam(object = M16$gam, newdata = dummy, se.fit = T)

dummy$ftt = dummyfit$fit
dummy$stt = dummyfit$se.fit



ttplot = 
  ggplot(dummy, aes(tchange8019, ftt))+
  geom_point(data = all.boreal, aes(tchange8019, trend.yp), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(tchange8019, ftt + stt), lwd = .5)+
  geom_line(aes(tchange8019, ftt - stt), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = expression(paste('Mean annual temperature change 1980-2019 (°C ' ~year^-1, ')')),
       y = expression(paste('TC change (% ' ~year^-1, ')')),
       tag = "F")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
ttplot





#Precipitation change 1980-2019
dummyfit = predict.gam(object = M17$gam, newdata = dummy, se.fit = T)

dummy$ftp = dummyfit$fit
dummy$stp = dummyfit$se.fit



tpplot = 
  ggplot(dummy, aes(pchange8019, ftp))+
  geom_point(data = all.boreal, aes(pchange8019, trend.yp), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(pchange8019, ftp + stp), lwd = .5)+
  geom_line(aes(pchange8019, ftp - stp), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = expression(paste('Mean annual precipitation change 1980-2019 (mm ' ~year^-1, ')')),
       y = expression(paste('TC change (% ' ~year^-1, ')')),
       tag = "G")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
tpplot








#Temperature change 2000-2019
dummyfit = predict.gam(object = M19$gam, newdata = dummy, se.fit = T)

dummy$ftt0019 = dummyfit$fit
dummy$stt0019 = dummyfit$se.fit



tt0019plot = 
  ggplot(dummy, aes(tchange0019, ftt0019))+
  geom_point(data = all.boreal, aes(tchange0019, trend.yp), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(tchange0019, ftt0019 + stt0019), lwd = .5)+
  geom_line(aes(tchange0019, ftt0019 - stt0019), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = expression(paste('Mean annual temperature change 2000-2019 (°C ' ~year^-1, ')')),
       y = expression(paste('TC change (% ' ~year^-1, ')')),
       tag = "H")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
tt0019plot





#Precipitation change 2000-2019
dummyfit = predict.gam(object = M18$gam, newdata = dummy, se.fit = T)

dummy$ftp0019 = dummyfit$fit
dummy$stp0019 = dummyfit$se.fit



tp0019plot = 
  ggplot(dummy, aes(pchange0019, ftp0019))+
  geom_point(data = all.boreal, aes(pchange0019, trend.yp), pch = 1, col = "grey", alpha = .6)+
  geom_line(lwd = 1.3)+
  geom_line(aes(pchange0019, ftp0019 + stp0019), lwd = .5)+
  geom_line(aes(pchange0019, ftp0019 - stp0019), lwd = .5)+
  
  geom_hline(yintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  geom_vline(xintercept = 0,
             lty = "dashed",
             lwd= 1)+
  
  facet_grid(.~dist_label, scale = "fixed", labeller = label_wrap_gen(width=20))+
  
  labs(x = expression(paste('Mean annual precipitation change 2000-2019 (mm ' ~year^-1, ')')),
       y = expression(paste('TC change (% ' ~year^-1, ')')),
       tag = "I")+
  theme_classic()+
  theme(legend.key.width = unit(2, units = "cm"),
        text = element_text(size = 23),
        panel.spacing = unit(2, "lines"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", face = "bold"),
        legend.position = "bottom")
tp0019plot


