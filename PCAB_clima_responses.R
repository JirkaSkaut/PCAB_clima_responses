setwd("C:/Users/jirka/Desktop/CR_PCAB/")

library(dplR) # tree ring analysis
library(treeclim) # climatic signal
library(ncdf4) # nc loading for climadata
library(raster)
library(rgdal)
library(SPEI) # SPEI calculation
library(ggplot2) # figures drawings
library(cowplot)
library(ggh4x)
library(tidytext)
library(ggtext)
library(scales)
library(Hmisc) # correlation matrix
library(tidyr) # data reorganization
library(reshape2)
library(relaimpo) # relative importance
library(car) # Variance Inflation Factor
library(multcompView)

## List of sites loading
PCAB_loc<- read.table("Data/Loc.txt", header=TRUE, sep="\t", na.strings="NA", dec=",", strip.white=TRUE, encoding = "UTF-8")

##############################################################
#         1. Calculation of TRI chronologies
##############################################################

Chron<- data.frame(matrix(ncol = nrow(PCAB_loc), nrow = 423)); colnames(Chron)<- PCAB_loc$CODE; rownames(Chron)<- c(1600:2022)

Table_S1<- data.frame(matrix(nrow = nrow(PCAB_loc), ncol = 5)); rownames(Table_S1)<- PCAB_loc$CODE; colnames(Table_S1)<- c("AGE", "AGE_std", "TRW_mean", "N_core","Rbar")

for (i in c(1:nrow(PCAB_loc))){
  
  serie <- read.rwl(paste("C:/Users/jirka/Desktop/CR_PCAB/Data/PCAB_RWL/", PCAB_loc[i, "CODE"], ".rwl", sep = ""))
  
  detrend_serie <- detrend(serie, method = "Spline", nyrs = 30)
  chronology <- chron(detrend_serie, biweight = T)
  
  stat<- rwl.stats(serie)
  info<- rwi.stats(detrend_serie)
  
  Table_S1[i, "AGE"]<- mean(stat$year)
  Table_S1[i, "AGE_std"]<- sd(stat$year)
  Table_S1[i, "TRW_mean"]<- mean(stat$mean)
  Table_S1[i, "N_core"]<- ncol(serie)
  Table_S1[i, "Rbar"]<- info$rbar.tot
  
  for (b in rownames(chronology)) {
    
    Chron[rownames(Chron[b,]==rownames(chronology[b,])), colnames(Chron) == PCAB_loc[i, "CODE"]] <- chronology[b,1]
    
  }
}

TRI<- Chron[rownames(Chron)<2011 & rownames(Chron)>1984,]

write.table(TRI, "Data/TRI.txt", row.names = T, col.names = T, sep = "\t")

##############################################################
#         2. Detrending of NDVI data
##############################################################

NDVI_orig <- read.table("Data/NDVI_orig.txt", check.names=FALSE, row.names = 1, header=TRUE, sep="\t", na.strings="NA", dec=",")

NDVI_det<- data.frame(matrix(ncol = ncol(NDVI_orig), nrow = nrow(NDVI_orig))); colnames(NDVI_det)<- colnames(NDVI_orig); rownames(NDVI_det)<- rownames(NDVI_orig)

for (a in c(1:ncol(NDVI_orig))) {
  
  model<- lm(NDVI_orig[,a] ~as.numeric(rownames(NDVI_orig)))
  
  NDVI_det[,a] <- model$residuals+mean(NDVI_orig[,a])
}

NDVI<- NDVI_det[rownames(NDVI_det)<2011& rownames(NDVI_det)>1984,]

write.table(NDVI, "Data/NDVI.txt", row.names = T, col.names = T, sep = "\t")

##############################################################
#         4. Climatic signal of TRI and NDVI
##############################################################

## upload of clima grids
temp <- nc_open("Data/Clima_grids/t_wgs_month.nc")
prec <- nc_open("Data/Clima_grids/sra_wgs_month.nc")

## tables for correlations
TRI_T_cor<- as.data.frame(matrix(ncol = 18, nrow = nrow(PCAB_loc))); rownames(TRI_T_cor)<- PCAB_loc$CODE
TRI_S_cor<- as.data.frame(matrix(ncol = 18, nrow = nrow(PCAB_loc))); rownames(TRI_S_cor)<- PCAB_loc$CODE

NDVI_T_cor<- as.data.frame(matrix(ncol = 18, nrow = nrow(PCAB_loc))); rownames(NDVI_T_cor)<- PCAB_loc$CODE
NDVI_S_cor<- as.data.frame(matrix(ncol = 18, nrow = nrow(PCAB_loc))); rownames(NDVI_S_cor)<- PCAB_loc$CODE

## tables for correlations significance
TRI_T_sig<- as.data.frame(matrix(ncol = 18, nrow = nrow(PCAB_loc))); rownames(TRI_T_sig)<- PCAB_loc$CODE
TRI_S_sig<- as.data.frame(matrix(ncol = 18, nrow = nrow(PCAB_loc))); rownames(TRI_S_sig)<- PCAB_loc$CODE

NDVI_T_sig<- as.data.frame(matrix(ncol = 18, nrow = nrow(PCAB_loc))); rownames(NDVI_T_sig)<- PCAB_loc$CODE
NDVI_S_sig<- as.data.frame(matrix(ncol = 18, nrow = nrow(PCAB_loc))); rownames(NDVI_S_sig)<- PCAB_loc$CODE

## tables for clima time-series
Temp_data<- as.data.frame(matrix(ncol = ncol(TRI), nrow = nrow(TRI))); colnames(Temp_data)<- colnames(TRI); rownames(Temp_data)<- rownames(TRI)
SPEI_data<- as.data.frame(matrix(ncol = ncol(TRI), nrow = nrow(TRI))); colnames(SPEI_data)<- colnames(TRI); rownames(SPEI_data)<- rownames(TRI)
Prec_data<- as.data.frame(matrix(ncol = ncol(TRI), nrow = nrow(TRI))); colnames(Prec_data)<- colnames(TRI); rownames(Prec_data)<- rownames(TRI)

Temp_data1<- as.data.frame(matrix(ncol = ncol(TRI), nrow = nrow(TRI))); colnames(Temp_data1)<- colnames(TRI); rownames(Temp_data1)<- rownames(TRI)
Prec_data1<- as.data.frame(matrix(ncol = ncol(TRI), nrow = nrow(TRI))); colnames(Prec_data1)<- colnames(TRI); rownames(Prec_data1)<- rownames(TRI)
SPEI_data1<- as.data.frame(matrix(ncol = ncol(TRI), nrow = nrow(TRI))); colnames(SPEI_data1)<- colnames(TRI); rownames(SPEI_data1)<- rownames(TRI)

for (i in c(1:nrow(PCAB_loc))) {
  
  TRI_chrono<- as.data.frame(TRI[,i]); rownames(TRI_chrono)<- rownames(TRI); colnames(TRI_chrono)<- colnames(TRI)[i]
  NDVI_chrono<- as.data.frame(NDVI[,i]); rownames(NDVI_chrono)<- rownames(TRI); colnames(NDVI_chrono)<- colnames(NDVI)[i]
  
  ## Tables for climadata of given site
  STORAGE_T <- data.frame(TIMESTAMP = seq(from = as.Date("19610101", tryFormats = c("%Y%m%d")), to = as.Date("20201201", tryFormats = c("%Y%m%d")), by = "month"))
  STORAGE_P <- data.frame(TIMESTAMP = seq(from = as.Date("19610101", tryFormats = c("%Y%m%d")), to = as.Date("20201201", tryFormats = c("%Y%m%d")), by = "month"))
  
  ## Coordinates of site
  site1x <- PCAB_loc[i, "LON"]; site1y <- PCAB_loc[i, "LAT"]
  
  ### Selecting pixel based on coordinates of site
  Order.X <- data.frame(ORDER = c(1:temp$dim$x$len), GRID = ncvar_get(nc = temp, varid = "x"))
  Order.Y <- data.frame(ORDER = c(1:temp$dim$y$len), GRID = ncvar_get(nc = temp, varid = "y"))
  Order.X$DIFFERENCE <- abs(Order.X$GRID - site1x)
  Order.Y$DIFFERENCE <- abs(Order.Y$GRID - site1y)
  site1x.order <- Order.X[Order.X$DIFFERENCE == min(Order.X$DIFFERENCE),"ORDER"] # Poradi pixelu, pro ktery je rozdil pozice stredu pixelu od zadaneho bodu nejmensi
  site1y.order <- Order.Y[Order.Y$DIFFERENCE == min(Order.Y$DIFFERENCE),"ORDER"]
  
  ### Extraction of climadata from given pixel
  STORAGE_T[,"Temperature"] <- ncvar_get(nc = temp, varid = "T", start = c(site1x.order, site1y.order, 1), count = c(1,1,-1)) 
  STORAGE_P[,"Precipitation"] <- ncvar_get(nc = prec, varid = "SRA", start = c(site1x.order, site1y.order, 1), count = c(1,1,-1)) 
  
  STORAGE_T$YEAR<- lapply(strsplit(as.character(STORAGE_T$TIMESTAMP), "\\-"), "[",1)
  STORAGE_T$MONTH<- lapply(strsplit(as.character(STORAGE_T$TIMESTAMP), "\\-"), "[",2)
  STORAGE_T$TIMESTAMP<- NULL
  Temp<- pivot_wider(STORAGE_T, names_from = MONTH, values_from = Temperature)
  Temp <- as.data.frame(lapply(Temp, unlist))
  Temp$YEAR<- as.numeric(Temp$YEAR)
  Temp_SPEI<- Temp[,2:13]
  rownames(Temp_SPEI)<- Temp$YEAR
  
  STORAGE_P$YEAR<- lapply(strsplit(as.character(STORAGE_P$TIMESTAMP), "\\-"), "[",1)
  STORAGE_P$MONTH<- lapply(strsplit(as.character(STORAGE_P$TIMESTAMP), "\\-"), "[",2)
  STORAGE_P$TIMESTAMP<- NULL
  Prec<- pivot_wider(STORAGE_P, names_from = MONTH, values_from = Precipitation)
  Prec <- as.data.frame(lapply(Prec, unlist))
  Prec$YEAR<- as.numeric(Prec$YEAR)
  Prec_SPEI<- Prec[,2:13]
  rownames(Prec_SPEI)<- Prec$YEAR
  
  ## calculation of SPEI
  PET<-thornthwaite(as.numeric(t(Temp_SPEI)),lat=50)
  
  SPEI_calc<-spei(as.numeric(t(Prec_SPEI))-PET,scale=2)
  
  SPEI<-as.data.frame(matrix(ncol=12,nrow=nrow(Temp_SPEI)))
  
  for(a in c(1:12))  {
    
    SPEI[,a]<-as.numeric(SPEI_calc$fitted)[seq(a,length(as.numeric(SPEI_calc$fitted)),12)]
    
  }
  
  colnames(SPEI)<-paste("S",1:12,sep="")
  rownames(SPEI)<-rownames(Temp_SPEI)
  
  SPEI$YEAR <- as.numeric(rownames(SPEI))
  SPEI<- SPEI[, c(13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]
  rownames(SPEI)<- NULL
  
  ## Climadata saving
  Temp_data[,i]<- rowMeans(Temp[c(25:50),7:9])
  SPEI_data[,i]<- rowMeans(SPEI[c(25:50),7:9])
  Prec_data[,i]<- rowMeans(Prec[c(25:50),7:9])
  
  Temp_data1[,i]<- rowMeans(Temp[c(25:50),2:13])
  SPEI_data1[,i]<- rowMeans(SPEI[c(25:50),2:13])
  Prec_data1[,i]<- rowSums(Prec[c(25:50),2:13])
  
  ## climasignal TRI
  TRI_COR_T <- dcc(TRI_chrono, Temp, selection= .range(-6:9) + .mean(6:8)+ .mean(-6:-8), method="correlation")
  TRI_COR_S <- dcc(TRI_chrono, SPEI, selection= .range(-6:9) + .mean(6:8)+ .mean(-6:-8), method="correlation")
  
  TRI_COR_T_coef<- TRI_COR_T$coef
  TRI_COR_S_coef<- TRI_COR_S$coef
  
  TRI_T_cor[i,]<- TRI_COR_T_coef$coef
  TRI_S_cor[i,]<- TRI_COR_S_coef$coef
  
  TRI_T_sig[i,]<- TRI_COR_T_coef$significant
  TRI_S_sig[i,]<- TRI_COR_S_coef$significant
  
  ## climasignal NDVI
  NDVI_COR_T <- dcc(NDVI_chrono, Temp, selection= .range(-6:9) + .mean(6:8)+ .mean(-6:-8), method="correlation")
  NDVI_COR_S <- dcc(NDVI_chrono, SPEI, selection= .range(-6:9) + .mean(6:8)+ .mean(-6:-8), method="correlation")
  
  NDVI_COR_T_coef<- NDVI_COR_T$coef
  NDVI_COR_S_coef<- NDVI_COR_S$coef
  
  NDVI_T_cor[i,]<- NDVI_COR_T_coef$coef
  NDVI_S_cor[i,]<- NDVI_COR_S_coef$coef
  
  NDVI_T_sig[i,]<- NDVI_COR_T_coef$significant
  NDVI_S_sig[i,]<- NDVI_COR_S_coef$significant
  
}

#### TRI
colnames(TRI_T_cor)<- TRI_COR_T_coef$month
colnames(TRI_S_cor)<- TRI_COR_S_coef$month

colnames(TRI_T_sig)<- TRI_COR_T_coef$month
colnames(TRI_S_sig)<- TRI_COR_T_coef$month

#### NDVI
colnames(NDVI_T_cor)<- NDVI_COR_T_coef$month
colnames(NDVI_S_cor)<- NDVI_COR_S_coef$month

colnames(NDVI_T_sig)<- NDVI_COR_T_coef$month
colnames(NDVI_S_sig)<- NDVI_COR_T_coef$month


TRI_T_cor$VAR<- "TRI"
TRI_S_cor$VAR<- "TRI"

TRI_T_sig$VAR<- "TRI"
TRI_S_sig$VAR<- "TRI"

NDVI_T_cor$VAR<- "NDVI"
NDVI_S_cor$VAR<- "NDVI"

NDVI_T_sig$VAR<- "NDVI"
NDVI_S_sig$VAR<- "NDVI"


#### TRI
TRI_T_cor$SITE<- rownames(TRI_T_cor)
TRI_T_cor<- gather(TRI_T_cor, "MONTH", "COR", 1:18)
TRI_T_cor$CLIMA<- "Temperature"

TRI_S_cor$SITE<- rownames(TRI_S_cor)
TRI_S_cor<- gather(TRI_S_cor, "MONTH", "COR", 1:18)
TRI_S_cor$CLIMA<- "SPEI"

TRI_T_sig$SITE<- rownames(TRI_T_sig)
TRI_T_sig<- gather(TRI_T_sig, "MONTH", "sig", 1:18)
TRI_T_sig$CLIMA<- "Temperature"

TRI_S_sig$SITE<- rownames(TRI_S_sig)
TRI_S_sig<- gather(TRI_S_sig, "MONTH", "sig", 1:18)
TRI_S_sig$CLIMA<- "SPEI"

#### NDVI
NDVI_T_cor$SITE<- rownames(NDVI_T_cor)
NDVI_T_cor<- gather(NDVI_T_cor, "MONTH", "COR", 1:18)
NDVI_T_cor$CLIMA<- "Temperature"

NDVI_S_cor$SITE<- rownames(NDVI_S_cor)
NDVI_S_cor<- gather(NDVI_S_cor, "MONTH", "COR", 1:18)
NDVI_S_cor$CLIMA<- "SPEI"

NDVI_T_sig$SITE<- rownames(NDVI_T_sig)
NDVI_T_sig<- gather(NDVI_T_sig, "MONTH", "sig", 1:18)
NDVI_T_sig$CLIMA<- "Temperature"

NDVI_S_sig$SITE<- rownames(NDVI_S_sig)
NDVI_S_sig<- gather(NDVI_S_sig, "MONTH", "sig", 1:18)
NDVI_S_sig$CLIMA<- "SPEI"

COR<- rbind(TRI_T_cor, TRI_S_cor, NDVI_T_cor, NDVI_S_cor)
SIG<- rbind(TRI_T_sig, TRI_S_sig, NDVI_T_sig, NDVI_S_sig)

COR$SIG<- SIG$sig

### JUN...AUG and Jun–Aug rename
for (i in c(1:nrow(COR))) {
  
  if(COR[i, "MONTH"]=="JUN...AUG"){
    
    COR[i, "MONTH"]<- "JUN–AUG"
    
  }
  
  if(COR[i, "MONTH"]=="jun...aug"){
    
    COR[i, "MONTH"]<- "Jun–Aug"
    
  }
  
  else(
    
    COR[i, "MONTH"]<- COR[i, "MONTH"]
    
  )
  
}

## adding explaining factors
Prec_data1<- as.data.frame(t(Prec_data1))
Temp_data1<- as.data.frame(t(Temp_data1))
SPEI_data1<- as.data.frame(t(SPEI_data1))

PCAB_loc$PREC<- rowMeans(Prec_data1)
PCAB_loc$MAT<- rowMeans(Temp_data1)
PCAB_loc$SPEI<- rowMeans(SPEI_data1)

COR<- merge(COR, PCAB_loc[,c(2, 8, 9, 10, 14, 15, 16, 19)], by.x="SITE", by.y="CODE", all.x=T)

## selecting of JJA correlations
COR_gg<- COR[COR$MONTH=="JUN–AUG",]

##############################################################
#            4. Model of TRI and NDVI climasignal
##############################################################

# intercorrelation of predictors
Intercor <- rcorr(as.matrix(COR_gg[c(7, 8, 9, 10, 13)]))

Cor_mat<- as.data.frame(Intercor$r)
P_mat<- as.data.frame(Intercor$P)

Cor_mat$VARs<- rownames(Cor_mat)
P_mat$VARs<- rownames(P_mat)

Cor_mat <- melt(Cor_mat, na.rm = T)
P_mat <- melt(P_mat, na.rm = T)

colnames(Cor_mat)<- c("Var1", "Var2", "COR")
colnames(P_mat)<- c("Var1", "Var2", "P_val")

Cor_mat$Code<- paste(Cor_mat$Var1, Cor_mat$Var2, sep = "_")
P_mat$Code<- paste(P_mat$Var1, P_mat$Var2, sep = "_")

Cor_mat<- merge(Cor_mat, P_mat[,c(4,3)], by.x="Code", by.y="Code", all.x=T)

order<- c("MAT", "TWI", "SLP", "HLI", "AGE")

Cor_mat<- Cor_mat[-c(2, 3, 4, 5, 8, 9, 10, 18, 20, 23),]

### Figure S2A (correlation matrix)
S3A<- ggplot(Cor_mat, aes(factor(Var2, level=order), factor(Var1, level=order), fill = COR, alpha= P_val<0.05))+ geom_tile()+
  scale_fill_gradient2(low = "#d7191c", mid= "#ffffbf", high = "#2c7bb6",  midpoint = 0, limit = c(-1,1))+
  geom_text(aes(Var2, Var1, label = round(COR, digits = 2)), color = "black", size = 4) +
  scale_alpha_discrete(range = c(0.4, 1))+
  guides(alpha = "none")+
  
  labs(x= "", y= "")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12, color="black"))+
  theme(axis.text.y = element_text(size = 12, color="black"))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 12))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 12))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "none", legend.margin=margin(t=-25))+ ## umístění legendy dolu
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white")) ## odstranění šedého podbaevení legend


### Models
COR_gg$CODE<- paste(COR_gg$VAR, COR_gg$CLIMA, sep = "–")

Levels<- as.data.frame(unique(COR_gg$CODE)); colnames(Levels)<- "CODE"
MAT_models<- as.data.frame(unique(COR_gg$CODE)); colnames(MAT_models)<- "CODE"

for (d in c(1:4)) {
  
  Set<- subset(COR_gg, subset=COR_gg$CODE==Levels[d, "CODE"])
  
  Model_clima<- lm(COR~AGE+MAT+TWI+SLP+HLI+as.factor(SOIL), data = Set)
  mod_summary<-summary(Model_clima)
  Rel_Imp<-calc.relimp(Model_clima)
  Rel_Imp<- as.data.frame(Rel_Imp$lmg)
  
  Levels[d, "R2"]<- summary(Model_clima)$r.squared
  Levels[d, "adj R2"]<- summary(Model_clima)$adj.r.squared
  Levels[d, "p_val"]<- pf(mod_summary$fstatistic[1],mod_summary$fstatistic[2],mod_summary$fstatistic[3],lower.tail = FALSE)
  
  Levels[d, "AGE_p"]<- mod_summary$coefficients["AGE",4]
  Levels[d, "MAT_p"]<- mod_summary$coefficients["MAT",4]
  Levels[d, "TWI_p"]<- mod_summary$coefficients["TWI",4]
  Levels[d, "SLP_p"]<- mod_summary$coefficients["SLP",4]
  Levels[d, "HLI_p"]<- mod_summary$coefficients["HLI",4]
  Levels[d, "SOIL_p"]<- NA
  
  Levels[d, "AGE_s"]<- cor(Set$COR, Set$AGE)
  Levels[d, "MAT_s"]<- cor(Set$COR, Set$MAT)
  Levels[d, "TWI_s"]<- cor(Set$COR, Set$TWI)
  Levels[d, "SLP_s"]<- cor(Set$COR, Set$SLP)
  Levels[d, "HLI_s"]<- cor(Set$COR, Set$HLI)
  Levels[d, "SOIL_s"]<- NA
  
  Levels[d, "AGE_r"]<- Rel_Imp["AGE",1]
  Levels[d, "MAT_r"]<- Rel_Imp["MAT",1]
  Levels[d, "TWI_r"]<- Rel_Imp["TWI",1]
  Levels[d, "SLP_r"]<- Rel_Imp["SLP",1]
  Levels[d, "HLI_r"]<- Rel_Imp["HLI",1]
  Levels[d, "SOIL_r"]<- Rel_Imp["as.factor(SOIL)",1]
  
  MAT_model<- lm(COR~MAT, data=Set)
  MAT_summary<-summary(MAT_model)
  
  MAT_models[d, "R2"]<- MAT_summary$r.squared
  MAT_models[d, "p_val"]<- pf(MAT_summary$fstatistic[1],MAT_summary$fstatistic[2],MAT_summary$fstatistic[3],lower.tail = FALSE)
  
}

### Figure S2B (variance inflation factor)
VIF<- as.data.frame(vif(Model_clima)); colnames(VIF)<- "VIF"
VIF$Variables<- rownames(VIF)
VIF[6, "Variables"]<- "SOIL"
VIF<- VIF[-6,]

S3B<- ggplot(VIF, aes(x=factor(Variables,level=order), y=VIF))+
  geom_bar(stat = "identity", fill="#2b8cbe")+
  
  labs(x= "", y= "Variance inflation factor")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(axis.text.x = element_text(size = 12, color="black"))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.y = element_text(size = 12, color="black"))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 12))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 12))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "right", legend.margin=margin(t=-25))+ ## umístění legendy dolu
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white")) ## odstranění šedého podbaevení legend

plot_grid(S3A, S3B, nrow=2, ncol=1, labels = c("A", "B"), label_size = 16)

ggsave("Graphs/Figure S3 (Factor relations).tiff", height = 180, width = 110, units = "mm", dpi = 300)

### Figure 2 (clima signal)

order1<- c("TRI",	"NDVI")

COR_gg$VAR_order <- factor(COR_gg$VAR, levels=c("TRI",	"NDVI"))

ggplot(COR_gg, aes(x=MAT, y=COR, color=CLIMA, alpha=SIG=="TRUE"))+geom_point()+
  geom_smooth(method=lm, aes(color=CLIMA, fill=CLIMA, alpha=SIG), alpha=0.5, level=0.95)+
  facet_wrap2(~VAR_order, axes = "all", ncol = 2)+

  labs(x= "Mean annual temperature (°C)", y= "Correlation")+
  scale_alpha_discrete(range = c(0.2, 1))+
  scale_x_continuous(breaks = seq(1, 9, 2))+
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3))+
  scale_color_manual(values = c("#f0d524", "#f54052"))+
  scale_fill_manual(values = c("#f0d524", "#f54052"))+
  geom_hline(yintercept=0, color = "black",linetype="dashed")+
  
  theme(strip.text  = element_text(size = 14))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 14, color="black"))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 14, color="black"))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 14, margin = margin(t = 5, r = 20, b = 3, l = 0)))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 14))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom", legend.margin=margin(t=-8))+ ## umístění legendy dolu
  guides(alpha = "none")+
  theme(axis.ticks = element_line(color = "black"))+
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white")) ## odstranění šedého podbaevení legend

ggsave("Graphs/Figure 2 (Clima signal).tiff", height = 100, width = 180, units = "mm", dpi = 300)

##############################################################
#         5. Correlation of TRI and NDVI
##############################################################

COR_TRI_NDVI<- as.data.frame(matrix(ncol = 3, nrow = nrow(PCAB_loc))); colnames(COR_TRI_NDVI)<- c("SITE","COR", "SIG")

for (c in c(1:nrow(PCAB_loc))) {
  
  correlation<- cor.test(TRI[,c], NDVI[,c])
  
  COR_TRI_NDVI[c, "SITE"]<- PCAB_loc[c,"CODE"]
  COR_TRI_NDVI[c,"COR"]<- correlation$estimate
  COR_TRI_NDVI[c,"SIG"]<- correlation$p.value
  
}

COR_TRI_NDVI<- merge(COR_TRI_NDVI, PCAB_loc[,c(2, 8, 9, 10, 14, 15, 16, 19)], by.x="SITE", by.y="CODE", all.x=T)

##############################################################
#            6. Model of TRI and NDVI reliationship
##############################################################

Model_TRI_NDVI<-lm(COR~AGE+MAT+TWI+SLP+HLI+as.factor(SOIL), data = COR_TRI_NDVI)

mod_summary3<-summary(Model_TRI_NDVI)
TRI_NDVI_ri<-calc.relimp(Model_TRI_NDVI)
TRI_NDVI_ri<- as.data.frame(TRI_NDVI_ri$lmg)

Levels[5, 1]<- "COR_TRI_NDVI"
Levels[5, "R2"]<- summary(Model_TRI_NDVI)$r.squared
Levels[5, "adj R2"]<- summary(Model_TRI_NDVI)$adj.r.squared
Levels[5, "p_val"]<- pf(mod_summary3$fstatistic[1],mod_summary3$fstatistic[2],mod_summary3$fstatistic[3],lower.tail = FALSE)

Levels[5, "AGE_p"]<- mod_summary3$coefficients["AGE",4]
Levels[5, "MAT_p"]<- mod_summary3$coefficients["MAT",4]
Levels[5, "TWI_p"]<- mod_summary3$coefficients["TWI",4]
Levels[5, "SLP_p"]<- mod_summary3$coefficients["SLP",4]
Levels[5, "HLI_p"]<- mod_summary3$coefficients["HLI",4]
Levels[5, "SOIL_p"]<- NA

Levels[5, "AGE_s"]<- cor(COR_TRI_NDVI$COR, COR_TRI_NDVI$AGE)
Levels[5, "MAT_s"]<- cor(COR_TRI_NDVI$COR, COR_TRI_NDVI$MAT)
Levels[5, "TWI_s"]<- cor(COR_TRI_NDVI$COR, COR_TRI_NDVI$TWI)
Levels[5, "SLP_s"]<- cor(COR_TRI_NDVI$COR, COR_TRI_NDVI$SLP)
Levels[5, "HLI_s"]<- cor(COR_TRI_NDVI$COR, COR_TRI_NDVI$HLI)
Levels[5, "SOIL_s"]<- NA

Levels[5, "AGE_r"]<- TRI_NDVI_ri["AGE",1]
Levels[5, "MAT_r"]<- TRI_NDVI_ri["MAT",1]
Levels[5, "TWI_r"]<- TRI_NDVI_ri["TWI",1]
Levels[5, "SLP_r"]<- TRI_NDVI_ri["SLP",1]
Levels[5, "HLI_r"]<- TRI_NDVI_ri["HLI",1]
Levels[5, "SOIL_r"]<- TRI_NDVI_ri["as.factor(SOIL)",1]

##############################

GG_models<- Levels[-5,-c(2,3,4)]

P_vals<- gather(GG_models[,c(1:7)], "FACTOR", "P_val", 2:7)
Efect<- gather(GG_models[,c(1,8:13)], "FACTOR", "EFECT", 2:7)
Relaimpo<- gather(GG_models[,c(1,14:19)], "FACTOR", "Relaimpo", 2:7)

GG_models<- cbind(P_vals, Efect$EFECT, Relaimpo$Relaimpo)
GG_models$FACTOR<- lapply(strsplit(as.character(GG_models$FACTOR), "\\_"), "[",1)
GG_models <- as.data.frame(lapply(GG_models, unlist))
colnames(GG_models)<- c("CODE", "FACTOR", "P_val", "EFECT", "Relaimpo")
GG_models$Relaimpo<- GG_models$Relaimpo*100

GG_models$CODE_f <- factor(GG_models$CODE, levels=c("TRI–SPEI", "NDVI–SPEI", "TRI–Temperature", "NDVI–Temperature"))

Model_info<- Levels[1:4,c(1,2,4)]
Model_info$R2<- sprintf("%0.2f", round(Model_info$R2, digits = 2))
Model_info$p_val<- signif(Model_info$p_val, digits = 2)

gg_text1<- data.frame(label = Model_info$R2, CODE = Model_info$CODE, x = c(5, 5, 5, 5),y = c(25, 25, 25, 25))
gg_text1$CODE_f<- factor(gg_text1$CODE, levels=c("TRI–SPEI", "NDVI–SPEI", "TRI–Temperature", "NDVI–Temperature"))

gg_text2<- data.frame(label = 0.05, CODE = Model_info$CODE, x = c(5, 5, 5, 5),y = c(20, 20, 20, 20))
gg_text2$CODE_f<- factor(gg_text2$CODE, levels=c("TRI–SPEI", "NDVI–SPEI", "TRI–Temperature", "NDVI–Temperature"))

## Figure 3 (Models)

ggplot(GG_models, aes(x=reorder_within(FACTOR, -Relaimpo, CODE), y=Relaimpo)) + 
  geom_bar(stat = "identity", aes(fill= EFECT>0, alpha = P_val<0.05, color = P_val<0.05))+ 
  labs(x ="", y= "Explained variability (%)")+
  facet_wrap2(~CODE_f, axes = "all", scales = "free")+
  coord_cartesian(ylim= c(0, 31))+
  scale_fill_manual(values = c("#d7191c", "#2c7bb6"))+
  scale_color_manual(values = c("#ffffff", "#000000"), guide=F)+
  
  theme(legend.position = "none")+
  scale_x_reordered()+
  scale_alpha_discrete(range = c(0.2, 1))+
  geom_richtext(data = gg_text1, mapping = aes(x = x, y = y, label = paste0("R<sup>2</sup> = ", label)), size = 5, fill = NA, label.color = NA)+
  geom_richtext(data = gg_text2, mapping = aes(x = x, y = y, label = paste0("p < ", label)), size = 5, fill = NA, label.color = NA)+
  
  theme(axis.text.x = element_text(size = 12, color="black"))+
  theme(axis.text.y = element_text(size = 12, color="black"))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.y = element_text(size = 15,vjust = +3))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(axis.ticks = element_line(color = "black"))+
  theme(panel.spacing.y = unit(1, "lines"))

ggsave("Graphs/Figure 3 (Models).tiff", height = 160, width = 260, units = "mm", dpi = 300)

### Fig. S4 Soil COR

ggplot(COR_gg, aes(x=MAT, y=COR))+geom_point(aes(color=COR, alpha=SIG=="TRUE"))+
  scale_color_gradient2(low = "#d7191c", mid= "#ffffbf", high = "#2c7bb6",  midpoint = 0, limit = c(-1,1), name="Correlation")+
  
  facet_grid2(SOIL~CODE, axes = "all", remove_labels = "all")+
  labs(x ="Mean annual temperature (°C)", y= "Correlation")+
  scale_alpha_discrete(range = c(0.4, 1))+
  guides(alpha = "none")+
  geom_smooth(method = "lm", color ="black")+
  scale_x_continuous(breaks = seq(1, 9, 2))+
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3))+
  geom_hline(yintercept=0, color = "black",linetype="dashed")+
  
  theme(axis.text.x = element_text(size = 12, color="black"))+
  theme(axis.text.y = element_text(size = 12, color="black"))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.text.x = element_text(size = 12))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.y = element_text(size = 12,vjust = +3))+ ## velikost popisků osy y
  theme(legend.position = "bottom", legend.margin=margin(t=-8))+ ## umístění legendy dolu
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(axis.ticks = element_line(color = "black"))+
  theme(panel.spacing.y = unit(1, "lines"))

ggsave("Graphs/Figure S4 (SOIL_COR).tiff", height = 200, width = 230, units = "mm", dpi = 300)

## Figure 4 (Cor TRI-NDVI vs MAT)

ggplot(COR_TRI_NDVI, aes(x=MAT, y=COR))+geom_point(aes(color = SIG<0.05))+
  geom_smooth(method = "lm", color="black")+
  labs(x ="Mean annual temperature (°C)", y= "Correlation")+
  scale_x_continuous(breaks = seq(1, 9, 2))+
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3))+
  scale_color_manual(values = c("#b0b4b8", "#2b8cbe"))+
  geom_hline(yintercept=0, color = "black",linetype="dashed")+
  
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size = 12, color="black"))+
  theme(axis.text.y = element_text(size = 12, color="black"))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.text.x = element_text(size = 12))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.y = element_text(size = 12, vjust = +3))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(axis.ticks = element_line(color = "black"))+
  theme(panel.spacing.y = unit(1, "lines"))

ggsave("Graphs/Figure 4 (Correlations TRI–NDVI).tiff", height = 90, width = 140, units = "mm", dpi = 300)
ggsave("Graphs/Graf_COR.png", height = 90, width = 120, units = "mm", dpi = 300)

## Figure S5A (Model)

#####
GG_model<- Levels[-c(1:4),-c(2,3,4)]

P_vals<- as.data.frame(t(GG_model[,c(2:7)]))
Efect<- as.data.frame(t(GG_model[,c(8:13)]))
Relaimpo<- as.data.frame(t(GG_model[,c(14:19)]))

GG_model<- cbind(P_vals, Efect[,1], Relaimpo[,1])
GG_model$FACTOR<- lapply(strsplit(as.character(rownames(GG_model)), "\\_"), "[",1)
GG_model <- as.data.frame(lapply(GG_model, unlist))
colnames(GG_model)<- c("P_val", "EFECT", "Relaimpo", "FACTOR")
GG_model$Relaimpo<- GG_model$Relaimpo*100

FS5A<- ggplot(GG_model, aes(x=reorder(FACTOR, -Relaimpo), y=Relaimpo)) + 
  geom_bar(stat = "identity", aes(fill= EFECT>0, alpha = P_val<0.05, color = P_val<0.05))+ 
  labs(x ="", y= "Explained variability (%)")+
  scale_y_continuous(breaks = seq(0, 13, 5))+
  
  scale_fill_manual(values = c("#d7191c", "#2c7bb6"))+
  scale_color_manual(values = c("#ffffff", "#000000"), guide=F)+
  
  theme(legend.position = "none")+
  scale_x_reordered()+
  scale_alpha_discrete(range = c(0.4, 1))+
  annotate("text", x = 5, y = 10, label = paste0("R^2 == ", round(Levels[5,2], digits = 2)), parse = TRUE, color = "black", size=5)+
  annotate("text", x = 5, y = 8, label = "p < 0.05", parse = TRUE, color = "black", size=5)+

  theme(axis.text.x = element_text(size = 12, color="black"))+
  theme(axis.text.y = element_text(size = 12, color="black"))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.text.x = element_text(size = 12))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.y = element_text(size = 12,vjust = +3))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(axis.ticks = element_line(color = "black"))+
  theme(panel.spacing.y = unit(1, "lines"))

FS5B<- ggplot(COR_TRI_NDVI, aes(x=MAT, y=COR, color=SOIL))+geom_point()+
  geom_smooth(method = "lm", aes(color=SOIL, fill=SOIL))+
  scale_color_manual(values = c("#eb5959", "#f7df43", "#31cf2b", "#3f47e0", "#ca7af5"), name = "")+
  scale_fill_manual(values = c("#eb5959", "#f7df43", "#31cf2b", "#3f47e0", "#ca7af5"), name = "")+
  labs(x ="Mean annual temperature (°C)", y= "Correlation")+
  scale_x_continuous(breaks = seq(1, 9, 2))+
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3))+
  geom_hline(yintercept=0, color = "black",linetype="dashed")+
  
  theme(legend.position = "bottom", legend.margin=margin(t=-8, l=-20))+
  theme(axis.text.x = element_text(size = 12, color="black"))+
  theme(axis.text.y = element_text(size = 12, color="black"))+
  theme(axis.title = element_text(size = 12))+
  theme(strip.text.x = element_text(size = 12))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.y = element_text(size = 12,vjust = +3))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(axis.ticks = element_line(color = "black"))+
  theme(panel.spacing.y = unit(1, "lines"))

plot_grid(FS5A, FS5B, nrow=2, ncol=1, labels = c("A", "B"), label_size = 16, rel_heights = c(1, 1, 1.1))

ggsave("Graphs/Figure S5 (Model).tiff", height = 150, width = 130, units = "mm", dpi = 300)

## Figure S6 (Clima correlations pJJA)

COR_pJJA<- COR[COR$MONTH=="Jun–Aug",]

COR_pJJA$VAR_order <- factor(COR_pJJA$VAR, levels=c("TRI",	"NDVI"))

ggplot(COR_pJJA, aes(x=MAT, y=COR, color=CLIMA, alpha=SIG=="TRUE"))+geom_point()+
  geom_smooth(method=lm, aes(color=CLIMA, fill=CLIMA, alpha=SIG), alpha=0.5, level=0.95)+
  facet_wrap2(~VAR_order, axes = "all", ncol = 2)+
  labs(x= "Mean annual temperature (°C)", y= "Correlation")+
  scale_alpha_discrete(range = c(0.2, 1))+
  scale_x_continuous(breaks = seq(1, 9, 2))+
  scale_y_continuous(breaks = seq(-0.6, 0.6, 0.3))+
  scale_color_manual(values = c("#f0d524", "#f54052"))+
  scale_fill_manual(values = c("#f0d524", "#f54052"))+
  geom_hline(yintercept=0, color = "black",linetype="dashed")+
  
  theme(strip.text  = element_text(size = 14))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 14, color="black"))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 14, color="black"))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 14, margin = margin(t = 5, r = 20, b = 3, l = 0)))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 14))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom", legend.margin=margin(t=-8))+ ## umístění legendy dolu
  guides(alpha = "none")+
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(axis.ticks = element_line(color = "black"))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white")) ## odstranění šedého podbaevení legend

ggsave("Graphs/Figure S6 (Clima signal pJJA).tiff", height = 100, width = 180, units = "mm", dpi = 300)

## Figure S2 (Clima correlations)

MAT<- as.data.frame(sort(unique(COR$MAT))); colnames(MAT)<-"MAT"
MAT$NUM<- as.numeric(rownames(MAT))

COR<- merge(COR, MAT, by.x="MAT", by.y="MAT", all.x=T)

order1<- c("Jun",	"Jul",	"Aug",	"Sep",	"Oct",	"Nov",	"Dec",	"JAN",	"FEB",	"MAR",	"APR",	"MAY",	"JUN",	"JUL",	"AUG",	"SEP", "Jun–Aug","JUN–AUG")

COR$VAR_CLIM <- paste(COR$VAR, COR$CLIMA, sep = "–")

COR$VAR_CLIM_order <- factor(COR$VAR_CLIM, levels=c("TRI–SPEI", "NDVI–SPEI", "TRI–Temperature", "NDVI–Temperature"))

ggplot(COR, aes(x=factor(MONTH,level=order1), y=NUM, fill=COR, alpha=SIG=="TRUE"))+geom_tile()+
  facet_wrap2(~VAR_CLIM_order, axes = "all")+
  scale_fill_gradient2(low = "#d7191c", mid= "#ffffbf", high = "#2c7bb6",  midpoint = 0, limit = c(-1,1), name="Correlation")+
  
  scale_y_continuous(breaks = c(1, 17, 46, 88, 121), labels = c(2, 4, 5, 7, 9))+
  
  labs(x= "", y= "Mean annual temperature (°C)")+
  guides(alpha = "none")+
  theme(strip.text  = element_text(size = 15))+ ## velikost nadpisů
  theme(axis.text.y = element_text(size = 12, color="black"))+ ## velikost názvu osy y
  theme(axis.text.x = element_text(size = 12, angle=45, hjust = 1, color="black"))+ ## velikost n?zvu osy x
  theme(axis.title.x = element_text(size = 12))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 15))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom", legend.margin=margin(t=-25))+ ## umístění legendy dolu
  theme(axis.ticks = element_line(color = "black"))+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

ggsave("Graphs/Figure S2 (Clima correlations).tiff", height = 200, width = 300, units = "mm", dpi = 300)

## Figure S1 (Variables)

TRI_model<-TRI
TRI_model$YEAR<- rownames(TRI_model)
TRI_model<- gather(TRI_model, "LOC", "TRI", 1:138)

NDVI_model<- NDVI
NDVI_model$YEAR<- rownames(NDVI_model)
NDVI_model<- gather(NDVI_model, "LOC", "NDVI", 1:138)

Temp_data$YEAR<- rownames(Temp_data)
Temp_data_model<- gather(Temp_data, "LOC", "Temp", 1:138)

Prec_data$YEAR<- rownames(Prec_data)
Prec_data_model<- gather(Prec_data, "LOC", "Prec", 1:138)

SPEI_data$YEAR<- rownames(SPEI_data)
SPEI_data_model<- gather(SPEI_data, "LOC", "SPEI", 1:138)

Model_data<- cbind(TRI_model, NDVI_model$NDVI, Temp_data_model$Temp, SPEI_data_model$SPEI)
colnames(Model_data)[4:6]<- c("NDVI", "June–August temperature", "June–August SPEI")

Model_data_gg<- gather(Model_data, "VAR", "VAL",3:6)

Model_data_gg<- merge(Model_data_gg, PCAB_loc[,c(2,19)], by.x="LOC", by.y="CODE", all.x=T)

Model_data_gg$VAR1 = factor(Model_data_gg$VAR, levels=c("TRI","NDVI", "June–August temperature", "June–August SPEI"))

ggplot(Model_data_gg, aes(x=as.numeric(YEAR), y=VAL))+geom_line(aes(group=LOC, color=MAT), linewidth=0.3)+
  facet_wrap2(~VAR1, scales = "free", axes = "all", nrow = 2)+
  scale_color_gradient2(low = "#0571b0", mid= "#f7f7f7", high = "#ca0020",  midpoint = mean(PCAB_loc$MAT), limit = c(min(PCAB_loc$MAT),max(PCAB_loc$MAT))) +
  labs(x= "", y= "")+
  
  theme(panel.spacing.y = unit(0.5, "lines"))+
  theme(panel.spacing.x = unit(0.5, "lines"))+
  theme(strip.text  = element_text(size = 14))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12, color="black"))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12, color="black"))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 14))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 14))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "none", legend.margin=margin(t= -25))+ ## umístění legendy dolu
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))+ ## odstranění šedého podbaevení legendy
  guides(col = guide_legend(nrow = 1))+ ## počet řádků legendy
  theme(axis.ticks = element_line(color = "black"))+
  theme(panel.spacing = unit(2, "lines"))

ggsave("Graphs/Figure S1 (Variables).tiff", height = 130, width = 200, units = "mm", dpi = 300)

#####################################
#     Temp vs PREC for LOC
#####################################

ggplot(PCAB_loc, aes(x=MAT, y=PREC, color=ELE))+geom_point(size=2)+
  labs(x= "Mean annual temperature (°C)", y= "Mean annual precipitation (mm)")+
  scale_color_gradient2(low = "#729a5a", mid= "#f0cb87", high = "#c28c7c",  midpoint = mean(PCAB_loc$ELE), limit = c(min(PCAB_loc$ELE),max(PCAB_loc$ELE))) +
  scale_x_continuous(breaks = seq(1, 9, 2))+
  
  theme(strip.text  = element_text(size = 14))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12, color="black"))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12, color="black"))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 14))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 14))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "none")+ ## umístění legendy dolu
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))+ ## odstranění šedého podbaevení legendy
  theme(axis.ticks = element_line(color = "black"))+
  theme(legend.margin = margin(t = -8, r = 0, b = 0, l = 0, unit = "pt"))

ggsave("Graphs/Graf_clima.png", height = 100, width = 120, units = "mm", dpi = 300)

### Range and percentage of significant correlations

max(COR_gg[COR_gg$CODE=="TRI–SPEI","COR"])
min(COR_gg[COR_gg$CODE=="TRI–SPEI","COR"])
nrow(COR_gg[COR_gg$CODE=="TRI–SPEI" & COR_gg$SIG=="TRUE",])/nrow(PCAB_loc)*100

max(COR_gg[COR_gg$CODE=="TRI–Temperature","COR"])
min(COR_gg[COR_gg$CODE=="TRI–Temperature","COR"])
nrow(COR_gg[COR_gg$CODE=="TRI–Temperature" & COR_gg$SIG=="TRUE",])/nrow(PCAB_loc)*100

mean(COR_gg[COR_gg$CODE=="NDVI–SPEI","COR"])
nrow(COR_gg[COR_gg$CODE=="NDVI–SPEI" & COR_gg$SIG=="TRUE",])/nrow(PCAB_loc)*100

mean(COR_gg[COR_gg$CODE=="NDVI–Temperature","COR"])
nrow(COR_gg[COR_gg$CODE=="NDVI–Temperature" & COR_gg$SIG=="TRUE",])/nrow(PCAB_loc)*100


nrow(COR_TRI_NDVI[COR_TRI_NDVI$SIG<0.05 & COR_TRI_NDVI$COR<0,])
nrow(COR_TRI_NDVI[COR_TRI_NDVI$SIG<0.05 & COR_TRI_NDVI$COR>0,])
nrow(COR_TRI_NDVI[COR_TRI_NDVI$SIG<0.05,])/nrow(COR_TRI_NDVI)*100


### Signiificance of COR~TEMP
summary(lm(COR_gg[COR_gg$CODE=="TRI–SPEI", "COR"]~ COR_gg[COR_gg$CODE=="TRI–SPEI", "MAT"]))
summary(lm(COR_gg[COR_gg$CODE=="TRI–Temperature", "COR"]~ COR_gg[COR_gg$CODE=="TRI–SPEI", "MAT"]))

summary(lm(COR_gg[COR_gg$CODE=="NDVI–SPEI", "COR"]~ COR_gg[COR_gg$CODE=="TRI–SPEI", "MAT"]))
summary(lm(COR_gg[COR_gg$CODE=="NDVI–Temperature", "COR"]~ COR_gg[COR_gg$CODE=="TRI–SPEI", "MAT"]))


### Autocorrelations

for (i in c(1:ncol(TRI))) {
  
  ACF_TRI<-acf(TRI[,i], plot = FALSE)
  ACF_TRI<-as.data.frame(ACF_TRI$acf); colnames(ACF_TRI)<- "ACF_TRI"
  
  Table_S1[i, "ACF_TRI"]<- ACF_TRI[2, "ACF_TRI"]
  
  ACF_NDVI<-acf(NDVI[,i], plot = FALSE)
  ACF_NDVI<-as.data.frame(ACF_NDVI$acf); colnames(ACF_NDVI)<- "ACF_NDVI"
  
  Table_S1[i, "ACF_NDVI"]<- ACF_NDVI[2, "ACF_NDVI"]
  
  ACF_Temp<-acf(Temp_data[,i], plot = FALSE)
  ACF_Temp<-as.data.frame(ACF_Temp$acf); colnames(ACF_Temp)<- "ACF_Temp"
  
  Table_S1[i, "ACF_Temp"]<- ACF_Temp[2, "ACF_Temp"]
  
  ACF_SPEI<-acf(SPEI_data[,i], plot = FALSE)
  ACF_SPEI<-as.data.frame(ACF_SPEI$acf); colnames(ACF_SPEI)<- "ACF_SPEI"
  
  Table_S1[i, "ACF_SPEI"]<- ACF_SPEI[2, "ACF_SPEI"]
  
}

Table_S1<- merge(Table_S1, PCAB_loc[,c(2,4,5,7)], by.x=0, by.y="CODE", all.x=T)

write.table(Table_S1, "Graphs/Table S1.txt", row.names = T, col.names = T, sep = "\t")

#####################################
#   TRI a NDVI time-serie example
#####################################

Abstract_data<- cbind(TRI_model[TRI_model$LOC=="C004009PCAB",], NDVI_model[NDVI_model$LOC=="C004009PCAB", 3])
colnames(Abstract_data)[4]<- "NDVI"

ylim.prim<- c(0, 1.3)
ylim.sec<- c(0.5, 0.85)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

ggplot(Abstract_data)+
  geom_line(aes(x=as.numeric(YEAR), y=TRI), linewidth=1.5, color="#965e03")+
  geom_line(aes(x=as.numeric(YEAR), y=a+NDVI*b), linewidth=1.5, color="#0cb00c")+
  scale_y_continuous("Tree-ring index", sec.axis = sec_axis(~ (. - a)/b, name = "NDVI"))+
  labs(x= "")+
  
  theme(strip.text  = element_text(size = 14))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12, color="black"))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12, color="black"))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 14))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 14))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 12))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom")+ ## umístění legendy dolu
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))+ ## odstranění šedého podbaevení legendy
  theme(axis.ticks = element_line(color = "black"))

ggsave("Graphs/Graf_TRI_NDVI.png", height = 100, width = 180, units = "mm", dpi = 300)

### Citations

citation("dplR")
citation("treeclim")
citation("ncdf4")
citation("raster")
citation("rgdal")
citation("SPEI")
citation("gstat")
citation("mgcv")
citation()