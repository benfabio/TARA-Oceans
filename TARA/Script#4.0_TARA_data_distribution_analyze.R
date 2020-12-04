
### 16/06/2020: Script to analyze the updated zooplankton imaging data from TARA Oceans (WP2, Bongo & Regent).
### Aims to: 
### - For each net, examine the normality of ESD and abunda data, apply transformations (log, log10, log1p, sqrt) and choose final
###   response vars as a fuxntion of N organisms per stations
### - Also, examine the normality of the hydro data, apply transformations, and choose final explanatory vars, examine their covariations 
### - Perform PCA
### - Examine covariance between ESD a,d abund per group (heatmap with correlation coeff)

### Last update: 15/07/2020

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("vegan")
library("FactoMineR")
library("matrixStats")
library("lubridate")
library("viridis")
library("maps")
library("ggrepel")

world <- map_data("world") # coastline for maps
WD <- getwd()

### Draw the correlation coeff heatmaps for each net
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
# Utiliser la corrélation entre les variables
  # comme mésure de distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}
augment.PCA <- function(x, dims = c(1:4), which="col") {
  .get <- function(x, element, dims) {
    y <- as.data.frame(x[[element]]$coord[,dims])
    if (nrow(y) == 0) {
      y <- NULL
    } else {
      y$type <- element
    }
    return(y)
  }
  if (which == "col") {
    y <- rbind(.get(x, "var", dims), .get(x, "quanti.sup", dims))
  } else {
    y <- rbind(.get(x, "ind", dims), .get(x, "quali.sup", dims))
  }
  y$var <- row.names(y)
  row.names(y) <- NULL
  return(y)
}

# --------------------------------------------------------------------------------------------------------------------------------


### 1°) Examine data distribution before and after various transformations  ------------------------------------------------------
abund <- read.table("table_abund+hydro_allnets_16_06_20.txt", h = T, sep = "\t")
dim(abund)
colnames(abund)
summary(abund$Copepoda)# ok, no NA
unique(abund$net)

esd <- read.table("table_ESD+hydro_allnets_16_06_20.txt", h = T, sep = "\t")
dim(esd)
colnames(esd)
summary(esd$Copepoda)# ok, no NA
unique(esd$net)


### A) Examine nb of NA per hydro data to discard those thta have most NAs (information loss too big)
for(net in unique(abund$net)) {
    
    for(c in colnames(abund)[c(38:90)] ) {
            # c <- "T_sur"
            n <- sum(length(which(is.na(abund[abund$net == net,c]))))
            nrows <- nrow(abund[abund$net == net,])
            message(paste("For net == ",net,"   ||  ",c, " == ",n,"/",nrows, " (",round(n/nrows,3)*100,"%)", sep = ""))
    } # eo 2nd for loop
    
    message(paste(" ", sep = ""))
    message(paste(" ", sep = ""))
    
} # eo 1st for loop - n in nets
### OK 

### Display only those with round(n/nrows,3)*100 < 25%
for(net in unique(abund$net)) {
    
    for(c in colnames(abund)[c(38:90)] ) {
            # c <- "T_sur"
            n <- sum(length(which(is.na(abund[abund$net == net,c]))))
            nrows <- nrow(abund[abund$net == net,])
            perc <- round(n/nrows,3)*100
            if(perc <= 30) {
                message(paste("For net == ",net,"   ||  ",c, " == ",n,"/",nrows, " (",perc,"%)", sep = ""))
            } else {
             message(paste(".", sep = ""))   
            }
    } # eo 2nd for loop
    
    message(paste(" ", sep = ""))
    message(paste(" ", sep = ""))
    
} # eo 1st for loop - n in nets

### Choose vars to keep per type (Temp, Sal, O2, nuts...)
## WP2:
# - Temp: Temp_10m vs Temp_DCM 
# - Sal: Sal_10m vs Sal_DCM 
# - O2: O2_10m vs O2_DCM ? Depth_max_O2/ Depth_min_O2
# - MLD
# - Zeu
# - Biomass: Chl_a_10m vs Chl_a_DCM ? Depth_chl_max
# - Nutrients: NO3_10m vs NO3_DCM
# - Backscattering: Kd_PAR_ vs beta470_sur vs bb470_sur vs fCDOM_sur vs bac660_sur

## Bongo: 
# - Temp: Temp_10m vs Temp_DCM 
# - Sal: Sal_10m vs Sal_DCM 
# - O2: O2_10m vs O2_DCM ? Depth_max_O2/ Depth_min_O2
# - MLD
# - Zeu
# - Biomass: Chl_a_10m vs Chl_a_DCM ? Depth_chl_max
# - Nutrients: NO3_10m vs NO3_DCM vs NO2_sur vs Phosphate_sur vs NO3_NO2_sur NO3_sur 
# - Backscattering: Kd_PAR_ vs Kd490 vs bbp470 vs acCDOM vs beta470_sur vs bb470_sur vs fCDOM_sur vs bac660_sur
# - PAR

## Regent: same as Bongo

### Draw the correlation coeff heatmaps for each net

names <- colnames(abund)[c(41:47,49:58,76,78,86:89)] ; names
mydata <- abund[abund$net == "wp2",names]
head(mydata)
cormat <- round(cor(na.omit(mydata), method = "spearman"),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

#quartz()
heatmap.wp2 <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
            limit = c(-1,1), space = "Lab", name = "Spearman correlation\ncoefficient (WP2)") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
    coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.ticks = element_blank(),
          legend.justification = c(1, 0), legend.position = c(0.6, 0.7), legend.direction = "horizontal") +
   guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))


### Now for Bongo + regent
names <- colnames(abund)[c(38:47,49:58,76,78,86:90)] ; names
mydata <- abund[abund$net == "bongo",names]
cormat <- round(cor(na.omit(mydata), method = "spearman"),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

#quartz()
heatmap.bongo <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
            limit = c(-1,1), space = "Lab", name = "Spearman correlation\ncoefficient (Bongo)") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
    coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.ticks = element_blank(),
          legend.justification = c(1, 0), legend.position = c(0.6, 0.7), legend.direction = "horizontal") +
   guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))


# Régent
names <- colnames(abund)[c(38:47,49:58,76,78,86:90)] ; names
mydata <- abund[abund$net == "regent",names]
cormat <- round(cor(na.omit(mydata), method = "spearman"),2)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

#quartz()
heatmap.regent <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
            limit = c(-1,1), space = "Lab", name = "Spearman correlation\ncoefficient (Regent)") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
    coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.ticks = element_blank(),
          legend.justification = c(1, 0), legend.position = c(0.6, 0.7), legend.direction = "horizontal") +
   guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))

### Save plots
ggsave(plot = heatmap.wp2, filename = "heatmap_spearman_hydro_WP2.jpg", dpi = 300, height = 15, width = 15)
ggsave(plot = heatmap.bongo, filename = "heatmap_spearman_hydro_bongo.jpg", dpi = 300, height = 15, width = 15)
ggsave(plot = heatmap.regent, filename = "heatmap_spearman_hydro_regent.jpg", dpi = 300, height = 15, width = 15)

### Ok so, choice of hydro vars 
# - WP2: Temp_10m, Sal_10m, Chla_10m, 02_10m, NO3_10m, MLD, bbp470 OR beta470 (choose based on distrb)
# - Bongo: Temp_10m, Sal_10m, Chla_10m, 02_10m, NO3_10m, MLD, PAR, bbp470 OR beta470 (choose based on distrb)
# - Regent: Temp_10m, Sal_10m, Chla_10m, 02_10m, NO3_10m, MLD, PAR, bbp470 OR beta470 (choose based on distrb)


### B) Examine distribution of bio and hydro
distrib.analyzer <- function(var, data) {
    
                    # Separate per dataset/net
                    for(net in unique(data$net)) {
                        
                        # Subset data
                        message(paste("Examining ",var," for ",net, sep = ""))
                        subset <- data[which(data$net == net),]
        
                        # Sometimes, the net dataset has no value (or not enough) for the variable of interest -> length(!is.na(subset$var)) == 0    
                        if( sum(!is.na(subset[,var])) < 5 ) {
                        
                            message(paste("\nNo values for ",var," for ",net," || moving to next variable", sep = ""))
                            
                        } else {
                                
                            # perform Shapiro-Wilk normality test
                            norm.test <- shapiro.test(subset[,var]) # if p-val < 0.05, var is not normally distrbuted
                            pval <- norm.test$p.value
                            if( pval < 0.05 ) { norm <- "Not normal" } else { norm <- "Normal" }
                            
                            # Plot original distribution, and then try some transformations: sqrt(), log10(), log(), cube root (^(1/3))
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Distribution_plots/imagery/abund/",net, sep = ""))
                            
                            plot <- ggplot(subset, aes(x = get(var))) + geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) +
                                        geom_density(alpha=.2, fill="#dd3497") + xlab(var) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                        theme_classic()
                                        
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_raw_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                                
                            # sqrt transformation 
                            norm.test <- shapiro.test(sqrt(subset[,var]))
                            pval <- norm.test$p.value
                            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                            plot <- ggplot(subset, aes(x = sqrt(get(var)))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (sqrt)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_sqrt_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                            
                            # log transformation 
                            norm.test <- shapiro.test(log(subset[,var]+0.0000001))
                            pval <- norm.test$p.value
                            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                            plot <- ggplot(subset, aes(x = log(get(var)))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_log_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                           
                            # log10 transformation 
                            norm.test <- shapiro.test(log10(subset[,var]+0.0000001))
                            pval <- norm.test$p.value
                            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                            plot <- ggplot(subset, aes(x = log10(get(var)))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_log10_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                                
                            # cubic transformation 
                            norm.test <- shapiro.test((subset[,var])^(1/3))
                            pval <- norm.test$p.value
                            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                            plot <- ggplot(subset, aes(x = (get(var))^(1/3))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (cubic)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_cubic_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                        
                            
                        } # eo if else loop - if no values 
                        
                        
            } # eo for loop - net types
        
} # eo fun

# Apply fun above in for loop for each bio and hydro var
# Change colnames #76 to PO4
# colnames(esd)[76] <- "PO4_sur"
# colnames(abund)[76] <- "PO4_sur"
# colnames(esd)[78] <- "SiO2_sur"
# colnames(abund)[78] <- "SiO2_sur"
variables <- c(colnames(abund)[c(1:33,91:95)],"Temp_10m","Sal_10m","Chl_a_10m","O2_10m","NO3_10m",
            "PO4_sur","SiO2_sur","MLD","PAR","bb470_sur","beta470_sur")
for(var in variables) {
    distrib.analyzer(var = var, data = abund)   
}

# ?log1p?
# shapiro.test(log1p(esd$Zooplankton))$p.value
# shapiro.test(log1p(esd$Calanoida))$p.value
# shapiro.test(log10(esd$Calanoida))$p.value
# shapiro.test(log(esd$Calanoida))$p.value
#
# shapiro.test(log1p(abund$Zooplankton))$p.value
# shapiro.test(log1p(abund$Calanoida))$p.value
# shapiro.test(log10(abund$Calanoida))$p.value
# shapiro.test(log(abund$Calanoida))$p.value

### Choice of transformations: 
# ESD: log transform
# Abundances: cubic transformation
# Temp, O2, Sal, PAR, MLD = keep raw 
# NO3, PO4, SiO2: log transform
# Chl_a: cubic transformation
# bbp470 (backscattering): log transform

### Add %micro/%nano/%pico from Hydro TARA Final 08_03_20 (and PAR2?)

### Now that you've mdae your choices, filter unwanted columns and apply transformations, and change colnames
colnames(esd)
sizes2 <- esd[,c(34:37,1:33,91:95,49,55,51,90,43,53,87)]
head(sizes2)
### Apply transformations
# Transform ESD data
colnames(sizes2)
sizes2transf <- colnames(sizes2)[c(5:42)] ; sizes2transf
sizes2[sizes2transf] <- lapply(sizes2[sizes2transf], log)
# log transfrom macro-nutrients concentrations + bbp470
sizes2transf <- colnames(sizes2)[c(49)] ; sizes2transf
sizes2[sizes2transf] <- lapply(sizes2[sizes2transf], log)
# cubic transform Chl-a
sizes2$Chl_a_10m <- (sizes2$Chl_a_10m)^(1/3)
# Change colnames
colnames(sizes2)[c(5:42)] <- paste("ESD_",colnames(sizes2)[c(5:42)], sep = "")
colnames(sizes2)[c(43:45)] <- c("Temperature","O2","Salinity")
colnames(sizes2)[c(48:49)] <- c("Chla","bbp470")

### Same for abund but cubic transform 
abund2 <- abund[,c(34:37,33,1:32,91:95,49,55,51,90,43,53,87)]
head(abund2) ; colnames(abund2)
### Apply transformations
# Transform abund data with cubic transf
abund2[,c(5:42)] <- (abund2[,c(5:42)])^(1/3)
# log transfrom macro-nutrients concentrations + bbp470
abund2transf <- colnames(abund2)[c(49)] ; abund2transf
abund2[abund2transf] <- lapply(abund2[abund2transf], log)
# cubic transform Chl-a
abund2$Chl_a_10m <- (abund2$Chl_a_10m)^(1/3)
# Change colnames
colnames(abund2)[c(5:42)] <- paste("Ab_",colnames(abund2)[c(5:42)], sep = "")
colnames(abund2)[c(43:45)] <- c("Temperature","O2","Salinity")
colnames(sizes2)[c(48:49)] <- c("Chla","bbp470")

### Go to main dir and complete abund2+sizes2 data with missing coords (-__-) + % pigments + PAR2 + flux200m
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/")
hydro.final <- read.csv("Hydro_final_TARA_08_03_20.csv", h = T, sep = ",") ; dim(hydro.final)
head(hydro.final)
### For sizes2
commons1 <- intersect(unique(sizes2$Station), unique(hydro.final$Station))
length(commons1)
commons2 <- intersect(unique(abund2$Station), unique(hydro.final$Station))
length(commons2)
# Find those that are lacking coords in sizes2
sizes2[which(sizes2$Station %in% commons1 & is.na(sizes2$Longitude)),c("Latitude","Longitude")] <- hydro.final[hydro.final$Station %in% c(1,2,15),c("Lat","Lon")]
abund2[which(abund2$Stations %in% commons2 & is.na(abund2$Longitude)),c("Latitude","Longitude")] <- hydro.final[hydro.final$Station %in% c(1,2,15),c("Lat","Lon")]

    
### Add pigments, PAR2 and flux200m
sizes2$Micro <- NA
sizes2$Nano <- NA
sizes2$Pico <- NA
sizes2$PAR2 <- NA
sizes2$Flux200 <- NA
# fill with for loop with common stations
for(c in commons1) {
        message(paste(c, sep = ""))
        sizes2[sizes2$Station == c,"Micro"] <- hydro.final[hydro.final$Station == c,"micro"]
        sizes2[sizes2$Station == c,"Nano"] <- hydro.final[hydro.final$Station == c,"nano"]
        sizes2[sizes2$Station == c,"Pico"] <- hydro.final[hydro.final$Station == c,"pico"]
        sizes2[sizes2$Station == c,"PAR2"] <- hydro.final[hydro.final$Station == c,"PAR"]
        sizes2[sizes2$Station == c,"Flux200"] <- hydro.final[hydro.final$Station == c,"flux200m"]
} # eo for loop 

# Same for abundances data
abund2$Micro <- NA
abund2$Nano <- NA
abund2$Pico <- NA
abund2$PAR2 <- NA
abund2$Flux200 <- NA
for(c in commons2) {
        message(paste(c, sep = ""))
        abund2[abund2$Stations == c,"Micro"] <- hydro.final[hydro.final$Station == c,"micro"]
        abund2[abund2$Stations == c,"Nano"] <- hydro.final[hydro.final$Station == c,"nano"]
        abund2[abund2$Stations == c,"Pico"] <- hydro.final[hydro.final$Station == c,"pico"]
        abund2[abund2$Stations == c,"PAR2"] <- hydro.final[hydro.final$Station == c,"PAR"]
        abund2[abund2$Stations == c,"Flux200"] <- hydro.final[hydro.final$Station == c,"flux200m"]
} # eo for loop 

### ANd nutrients
sizes2$NO2NO3 <- NA
sizes2$SiO2 <- NA
sizes2$PO4 <- NA
# fill with for loop with common stations
for(c in commons1) {
        message(paste(c, sep = ""))
        sizes2[sizes2$Station == c,"NO2NO3"] <- hydro.final[hydro.final$Station == c,"NO2NO3"]
        sizes2[sizes2$Station == c,"SiO2"] <- hydro.final[hydro.final$Station == c,"Si"]
        sizes2[sizes2$Station == c,"PO4"] <- hydro.final[hydro.final$Station == c,"PO4"]
} # eo for loop 

# Same for abundances data
abund2$NO2NO3 <- NA
abund2$SiO2 <- NA
abund2$PO4 <- NA
for(c in commons2) {
        message(paste(c, sep = ""))
        abund2[abund2$Station == c,"NO2NO3"] <- hydro.final[hydro.final$Station == c,"NO2NO3"]
        abund2[abund2$Station == c,"SiO2"] <- hydro.final[hydro.final$Station == c,"Si"]
        abund2[abund2$Station == c,"PO4"] <- hydro.final[hydro.final$Station == c,"PO4"]
} # eo for loop 
# Check nNA
summary(sizes2)
summary(abund2)

distrib.analyzer <- function(var, data) {
    
                    # Separate per dataset/net
                    for(net in unique(data$net)) {
                        
                        # Subset data
                        message(paste("Examining ",var," for ",net, sep = ""))
                        subset <- data[which(data$net == net),]
        
                        # Sometimes, the net dataset has no value (or not enough) for the variable of interest -> length(!is.na(subset$var)) == 0    
                        if( sum(!is.na(subset[,var])) < 5 ) {
                        
                            message(paste("\nNo values for ",var," for ",net," || moving to next variable", sep = ""))
                            
                        } else {
                                
                            # perform Shapiro-Wilk normality test
                            norm.test <- shapiro.test(subset[,var]) # if p-val < 0.05, var is not normally distrbuted
                            pval <- norm.test$p.value
                            if( pval < 0.05 ) { norm <- "Not normal" } else { norm <- "Normal" }
                            
                            # Plot original distribution, and then try some transformations: sqrt(), log10(), log(), cube root (^(1/3))
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Distribution_plots/imagery/abund/",net, sep = ""))
                            
                            plot <- ggplot(subset, aes(x = get(var))) + geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) +
                                        geom_density(alpha=.2, fill="#dd3497") + xlab(var) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                        theme_classic()
                                        
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_raw_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                                
                            # sqrt transformation 
                            norm.test <- shapiro.test(sqrt(subset[,var]))
                            pval <- norm.test$p.value
                            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                            plot <- ggplot(subset, aes(x = sqrt(get(var)))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (sqrt)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_sqrt_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                            
                            # log transformation 
                            norm.test <- shapiro.test(log(subset[,var]+0.0000001))
                            pval <- norm.test$p.value
                            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                            plot <- ggplot(subset, aes(x = log(get(var)))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_log_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                           
                            # log10 transformation 
                            norm.test <- shapiro.test(log10(subset[,var]+0.0000001))
                            pval <- norm.test$p.value
                            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                            plot <- ggplot(subset, aes(x = log10(get(var)))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_log10_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                                
                            # cubic transformation 
                            norm.test <- shapiro.test((subset[,var])^(1/3))
                            pval <- norm.test$p.value
                            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                            plot <- ggplot(subset, aes(x = (get(var))^(1/3))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (cubic)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                            ggsave(plot = plot, filename = paste("plot_distrib_",var,"_cubic_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                        
                            
                        } # eo if else loop - if no values 
                        
                        
            } # eo for loop - net types
        
} # eo fun


variables <- c("NO2NO3","SiO2","PO4")
for(var in variables) {
    distrib.analyzer(var = var, data = abund2)   
}

### Nutrients need to be cubic transformed in sizes2 and abund2
sizes2$NO2NO3 <- (sizes2$NO2NO3)^(1/3)
sizes2$SiO2 <- (sizes2$SiO2)^(1/3)
sizes2$PO4 <- (sizes2$PO4)^(1/3)
abund2$NO2NO3 <- (abund2$NO2NO3)^(1/3)
abund2$SiO2 <- (abund2$SiO2)^(1/3)
abund2$PO4 <- (abund2$PO4)^(1/3)

### Save both tables, next step will be plotting spatial patterns :-) 
setwd(WD)
write.table(abund2, "table_abund+hydro_allnets_17_06_20.txt", sep = "\t")
write.table(sizes2, "table_ESD+hydro_allnets_17_06_20.txt", sep = "\t")


### 2°) Re-examine collinearity  --------------------------------------------------------------------------

### Load the transformed data from above
abund <- read.table("table_abund+hydro_allnets_17_06_20.txt", sep = "\t", h = T)
colnames(abund)[2] <- "Station"
esd <- read.table("table_ESD+hydro_allnets_17_06_20.txt", sep = "\t", h = T)
dim(abund) ; dim(esd)
head(abund) ; head(esd)
# net <- "wp2"
for(net in unique(abund$net)) {
    
    names <- colnames(abund)[c(43:57)] #; names
    mydata <- abund[abund$net == net,names]
    cormat <- round(cor(na.omit(mydata), method = "spearman"),2)
    upper_tri <- get_upper_tri(cormat)
    melted_cormat <- melt(upper_tri, na.rm = TRUE)

    #quartz()
    htmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                limit = c(-1,1), space = "Lab", name = "Spearman correlation\ncoefficient (Regent)") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
        coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              panel.grid.major = element_blank(), panel.border = element_blank(),
              panel.background = element_blank(), axis.ticks = element_blank(),
              legend.justification = c(1, 0), legend.position = c(0.6, 0.7), legend.direction = "horizontal") +
       guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
    
    ggsave(plot = htmap, filename = paste("heatmap_spearman_hydro_transf_",net,".jpg", sep = ""), dpi = 300, height = 7.5, width = 7.5)
    
} # eo for loop 

### Conclusion: 
# - choose between PAR and PAR2
quartz()
ggplot(data = abund[abund$net == "wp2",]) + geom_point(aes(x = PAR, y = PAR2)) + theme_classic()
# Relationshop's the same for all 3 nets
# Check distrib of PAR vs PAR2
norm.test <- shapiro.test(abund[abund$net == "wp2","PAR"])
pval <- norm.test$p.value
if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
quartz()
ggplot(abund[abund$net == "wp2",], aes(x = PAR)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
    xlab(paste("PAR",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
    theme_classic()
#
norm.test <- shapiro.test(abund[abund$net == "wp2","PAR2"])
pval <- norm.test$p.value
if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
quartz()
ggplot(abund[abund$net == "wp2",], aes(x = PAR2)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
    xlab(paste("PAR v2",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
    theme_classic()
    
summary(abund[abund$net == "bongo","PAR"]) 
summary(abund[abund$net == "bongo","PAR2"]) 
# PAR2 slightly less normally distrbuted but has way less misisng values
# Map both
ggplot() + geom_point(aes(x = Longitude, y = Latitude, fill = PAR), data = abund[abund$net == "wp2",], pch = 21, colour = "black", size = 2.5) +
	scale_fill_viridis(name = "PAR v1") + 
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "#BAB0AC", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,190),
               labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right")
#
ggplot() + geom_point(aes(x = Longitude, y = Latitude, fill = PAR2), data = abund[abund$net == "wp2",], pch = 21, colour = "black", size = 2.5) +
	scale_fill_viridis(name = "PAR v2") + 
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "#BAB0AC", colour = "black", size = 0.3) +
 	coord_quickmap() + scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,190),
               labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
 	scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90),
 		      labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
   	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "right")    

### --> Choose PARv2

# - choose NO2NO3 over PO4 and SiO2
# - choose %micro over %pico


### 3°) Perform PCA on all predictors (but PARv1)  --------------------------------------------------------------------------
require("FactoMineR")
require("ggrepel")
# net <- "wp2"
esd$PC1 <- NA
esd$PC2 <- NA
for(net in unique(esd$net)) {
    
    message(paste("Performing PCA for ESD data based on ",net, sep = ""))
    names <- colnames(esd)[c(1:4,43:45,47:53,55:57)] # names
    data.pca <- na.omit(esd[esd$net == net,names])
    pca <- PCA(data.pca[,c(5:length(data.pca))], scale.unit = T, ncp = 4, graph = F)
    # summary(pca)
    pc1 <- round(pca$eig[,"percentage of variance"][1], digits = 2)
    pc2 <- round(pca$eig[,"percentage of variance"][2], digits = 2)
    pc3 <- round(pca$eig[,"percentage of variance"][3], digits = 2)
    pc4 <- round(pca$eig[,"percentage of variance"][4], digits = 2)
    message(paste(pc1+pc2+pc2+pc4,"% of variance explained", sep = ""))
    pcad <- augment.PCA(pca)
    plot.pca1 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
        annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
        geom_segment(aes(x=0, xend = Dim.1, y=0, yend = Dim.2, colour = factor(type)), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
        scale_colour_manual(name = "", values = c("#c51b7d"), guide = F) + 
        geom_text_repel(aes(x=Dim.1, y=Dim.2, colour=type, label=var), 
                data=filter(pcad, (Dim.1^2+Dim.2^2) > 0.2^2), segment.alpha=0.5) +
        xlab(paste("PC1 ","(",pc1,"%)", sep = "")) + ylab(paste("PC2 ","(",pc2,"%)", sep = "")) + theme_bw()
    #
    plot.pca2 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
        annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
        geom_segment(aes(x=0, xend = Dim.3, y=0, yend = Dim.4, colour = factor(type)), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
        scale_colour_manual(name = "", values = c("#c51b7d"), guide = F) + 
        geom_text_repel(aes(x=Dim.3, y=Dim.4, colour=type, label=var), 
                data=filter(pcad, (Dim.3^2+Dim.4^2) > 0.2^2), segment.alpha=0.5) +
        xlab(paste("PC3 ","(",pc3,"%)", sep = "")) + ylab(paste("PC4 ","(",pc4,"%)", sep = "")) + theme_bw()
        
    # Save plots
    ggsave(plot = plot.pca1, filename = paste("plot_pca_ESD_",net,"_pc1+pc2.jpg", sep = ""), dpi = 300, width = 5, height = 5)
    ggsave(plot = plot.pca2, filename = paste("plot_pca_ESD_",net,"_pc3+pc4.jpg", sep = ""), dpi = 300, width = 5, height = 5)
    
    # provide PC scores (PC1 and PC2) to data.pca and then 'esd'
    message(paste("Providing PC scores for ESD data based on ",net, sep = ""))
    data.pca$PC1 <- pca$ind$coord[,1]
    data.pca$PC2 <- pca$ind$coord[,2]
    esd[esd$net == net & esd$Station %in% data.pca$Station,"PC1"] <- data.pca$PC1
    esd[esd$net == net & esd$Station %in% data.pca$Station,"PC2"] <- data.pca$PC2

} # eo for loop 
### Check
summary(esd)

### Same for abund
abund$PC1 <- NA
abund$PC2 <- NA
for(net in unique(abund$net)) {
    
    message(paste("Performing PCA for ABUND data based on ",net, sep = ""))
    names <- colnames(abund)[c(1:4,43:45,47:53,55:57)] # names
    data.pca <- na.omit(abund[abund$net == net,names])
    pca <- PCA(data.pca[,c(5:length(data.pca))], scale.unit = T, ncp = 4, graph = F)
    # summary(pca)
    pc1 <- round(pca$eig[,"percentage of variance"][1], digits = 2)
    pc2 <- round(pca$eig[,"percentage of variance"][2], digits = 2)
    pc3 <- round(pca$eig[,"percentage of variance"][3], digits = 2)
    pc4 <- round(pca$eig[,"percentage of variance"][4], digits = 2)
    
    message(paste(pc1+pc2+pc2+pc4,"% of variance explained", sep = ""))
    
    pcad <- augment.PCA(pca)
    plot.pca1 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
        annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
        geom_segment(aes(x=0, xend = Dim.1, y=0, yend = Dim.2, colour = factor(type)), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
        scale_colour_manual(name = "", values = c("#c51b7d"), guide = F) + 
        geom_text_repel(aes(x=Dim.1, y=Dim.2, colour=type, label=var), 
                data=filter(pcad, (Dim.1^2+Dim.2^2) > 0.2^2), segment.alpha=0.5) +
        xlab(paste("PC1 ","(",pc1,"%)", sep = "")) + ylab(paste("PC2 ","(",pc2,"%)", sep = "")) + theme_bw()
    #
    plot.pca2 <- ggplot(pcad) + coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
        annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
        geom_segment(aes(x=0, xend = Dim.3, y=0, yend = Dim.4, colour = factor(type)), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
        scale_colour_manual(name = "", values = c("#c51b7d"), guide = F) + 
        geom_text_repel(aes(x=Dim.3, y=Dim.4, colour=type, label=var), 
                data=filter(pcad, (Dim.3^2+Dim.4^2) > 0.2^2), segment.alpha=0.5) +
        xlab(paste("PC3 ","(",pc3,"%)", sep = "")) + ylab(paste("PC4 ","(",pc4,"%)", sep = "")) + theme_bw()
        
    # Save plots
    ggsave(plot = plot.pca1, filename = paste("plot_pca_abund_",net,"_pc1+pc2.jpg", sep = ""), dpi = 300, width = 5, height = 5)
    ggsave(plot = plot.pca2, filename = paste("plot_pca_abund_",net,"_pc3+pc4.jpg", sep = ""), dpi = 300, width = 5, height = 5)
    
    # provide PC scores (PC1 and PC2) to data.pca and then 'esd'
    message(paste("Providing PC scores for abund data based on ",net, sep = ""))
    data.pca$PC1 <- pca$ind$coord[,1]
    data.pca$PC2 <- pca$ind$coord[,2]
    abund[abund$net == net & abund$Station %in% data.pca$Station,"PC1"] <- data.pca$PC1
    abund[abund$net == net & abund$Station %in% data.pca$Station,"PC2"] <- data.pca$PC2
    
} # eo for loop 
### Gut ! 


### 4°) Examine covariance between the groups' abundance vs. ESD  --------------------------------------------------------------------------

### First, need to get the counts of individuals per group and stations to account for those stations that have sufficent N for decent median ESD
counts.bongo <- read.table("table_counts_groups_bongo.txt", h = T, sep = "\t")
counts.wp2 <- read.table("table_counts_groups_wp2.txt", h = T, sep = "\t")
counts.regent <- read.table("table_counts_groups_regent.txt", h = T, sep = "\t")
#dim(counts.bongo) ; dim(counts.wp2) ; dim(counts.regent)


### In a lapply, per group and net, merge corresponding ESD and abund data. Restrict to stations that have n > 25 for ESD and abund
nets <- c("bongo","wp2","regent")
net <- "wp2"
res <- lapply(nets, function(net) {
    
        message(paste("", sep = ""))
        message(paste("Examining covariance between abundance and size for ",net, sep = ""))
        message(paste("", sep = ""))
        sizes <- esd[esd$net == net,]
        abs <- abund[abund$net == net,]
        
        # 2nd lapply based on the groups 
        names <- colnames(sizes)[c(5:42)]
        groups <- str_replace_all(names,"ESD_","")
        colnames(sizes)[c(5:42)] <- groups
        colnames(abs)[c(5:42)] <- groups
 
        # unique(counts$group)
        if(net == "bongo") {
            counts <- counts.bongo
        } else if(net == "wp2") {
            counts <- counts.wp2
        } else {
            counts <- counts.regent
        } # eo if else loop
        # unique(counts$group)
        ### Change some levels in counts so they match the colnames of abund/esd
        levels(counts$group)[levels(counts$group) == "Calanoida_unid"] <- "Calanoida_unidentified"
        levels(counts$group)[levels(counts$group) == "Copepoda_unid"] <- "Copepoda_unidentified"
        levels(counts$group)[levels(counts$group) == "Crustaceans"] <- "Crust"
        levels(counts$group)[levels(counts$group) == "Grazers"] <- "Small_grazers"
       
        common.groups <- intersect(groups, unique(counts$group))
        
        # gr <- "Copepoda"
        grps <- lapply(common.groups, function(gr) {
            
                    message(paste(gr, sep = ""))            
                    stations2keep <- counts[counts$group == gr & counts$n >= 25 ,"station"]
                    gr.size <- sizes[sizes$Station %in% stations2keep, c("Station",gr)]
                    gr.abund <- abs[abs$Station %in% stations2keep, c("Station",gr)]
                    # na.omit
                    gr.size <- na.omit(gr.size)
                    gr.abund <- na.omit(gr.abund)
                    gr.size <- gr.size[order(gr.size$Station),]
                    gr.abund <- gr.abund[order(gr.abund$Station),]
                  
                    # Restrict co common stations to compute correlation
                    common_stations <- intersect(unique(gr.size$Station), unique(gr.abund$Station))
                    gr.size2 <- gr.size[which(gr.size$Station %in% common_stations),]
                    gr.abund2 <- gr.abund[which(gr.abund$Station %in% common_stations),]
                    # dim(gr.size2) ; dim(gr.abund2)
                    
                    # There can be several abund measurements for abund
                    gr.abund3 <- data.frame(gr.abund2 %>% group_by(Station) %>% summarize(abund = mean(!! sym(gr)) ) ) # noice
                    
                    # If both esd and abund gave same dimensios and there are enough stations (n = 10)
                    # na.omit
                    if(nrow(gr.size2) == nrow(gr.abund3) & nrow(gr.size2) > 10) {
                        
                        rho <- cor(gr.size2[,2], gr.abund3[,2], method = "spearman")
                        pval <- cor.test(x = gr.size2[,2], y = gr.abund3[,2], method = "spearman")$p.value
                        # Put in a ddf
                        res.cor <- data.frame(group = gr, net = net, n = nrow(gr.size2), rho = rho, pval = pval)
                        return(res.cor)
                        
                    } else {
                        
                         message(paste("Not enough stations for ", gr, sep = ""))     
                        
                    }
                    
            } # eo 2nd FUN
            
        ) # eo 2nd LAPPLY
        ### Rbind results
        ddf <- bind_rows(grps)
        rm(grps)
        # ddf
        return(ddf)
        
    } # eo 1st FUN
    
) # eo 1st LAPPLY
ddf <- bind_rows(res)
dim(ddf)
rm(res)
#                     group    net   n          rho         pval
# 1                Copepoda  bongo 117  0.258467321 4.897312e-03
# 2            Chaetognatha  bongo  95  0.238698223 1.983131e-02
# 3                Gel_carn  bongo  49  0.030529117 8.350451e-01
# 4                  Gel_FF  bongo  86  0.033588228 7.588334e-01
# 5                   Crust  bongo  72  0.032904424 7.837841e-01
# 6               Pteropoda  bongo  36  0.380485089 2.206756e-02
# 7                Rhizaria  bongo  74  0.136086444 2.476331e-01
# 8           Small_grazers  bongo 100  0.025961104 7.976495e-01
# 9               Calanoida  bongo 117  0.191051777 3.907461e-02
# 10          Harpacticoida  bongo  33  0.351936911 4.458545e-02
# 11           Augaptilidae  bongo  62 -0.130001876 3.138937e-01
# 12          Centropagidae  bongo  12 -0.074609889 8.177457e-01
# 13            Eucalanidae  bongo  45  0.098446254 5.199778e-01
# 14        Heterorhabdidae  bongo  86 -0.012664613 9.078643e-01
# 15  Copepoda_unidentified  bongo 116  0.264208082 4.158247e-03
# 16               Copepoda    wp2 152  0.461612057 2.155102e-09
# 17           Chaetognatha    wp2  96  0.080835363 4.336689e-01
# 18               Gel_carn    wp2  37  0.073250854 6.665744e-01
# 19                 Gel_FF    wp2  98  0.194668703 5.475509e-02
# 20                  Crust    wp2  31 -0.220806715 2.326027e-01
# 21              Pteropoda    wp2  45  0.288158350 5.491277e-02
# 22               Rhizaria    wp2  37 -0.220662801 1.893796e-01
# 23          Small_grazers    wp2  81  0.079609628 4.798941e-01
# 24              Calanoida    wp2 149  0.364605402 4.843295e-06
# 25          Harpacticoida    wp2  46  0.289669581 5.086030e-02
# 26      Poecilostomatoida    wp2 125  0.345500510 7.938817e-05
# 27 Calanoida_unidentified    wp2  34  0.532645024 1.183073e-03
# 28            Corycaeidae    wp2  59  0.124029959 3.493098e-01
# 29            Eucalanidae    wp2  18 -0.258842611 2.996669e-01
# 30           Metridinidae    wp2  13  0.236458025 4.366955e-01
# 31              Oncaeidae    wp2  17  0.073764754 7.784369e-01
# 32          Sapphirinidae    wp2  27 -0.139878985 4.865085e-01
# 33              Temoridae    wp2  43 -0.070901102 6.514139e-01
# 34               Copepoda regent 117  0.427499863 1.529892e-06
# 35           Chaetognatha regent  64 -0.023439281 8.541361e-01
# 36               Gel_carn regent  61  0.002451182 9.850419e-01
# 37                 Gel_FF regent  14  0.607744082 2.113796e-02
# 38                  Crust regent  34  0.160199892 3.654428e-01
# 39               Rhizaria regent  15  0.103845984 7.126496e-01
# 40          Small_grazers regent  27 -0.208168730 2.974289e-01
# 41              Calanoida regent 114  0.430932654 1.699927e-06
# 42            Candaciidae regent  19 -0.346961859 1.455691e-01
# 43            Euchaetidae regent  47  0.212920130 1.507360e-01
# 44  Copepoda_unidentified regent  45  0.170545466 2.626755e-01


### Check those groups that display signif. positive relationships between abundance and median size
ddf[ddf$pval < 0.05,]
#                      group    net   n       rho         pval
# 1                Copepoda  bongo 117 0.2584673 4.897312e-03
# 2            Chaetognatha  bongo  95 0.2386982 1.983131e-02
# 6               Pteropoda  bongo  36 0.3804851 2.206756e-02
# 9               Calanoida  bongo 117 0.1910518 3.907461e-02
# 10          Harpacticoida  bongo  33 0.3519369 4.458545e-02
# 15  Copepoda_unidentified  bongo 116 0.2642081 4.158247e-03

# 16               Copepoda    wp2 152 0.4616121 2.155102e-09
# 24              Calanoida    wp2 149 0.3646054 4.843295e-06
# 26      Poecilostomatoida    wp2 125 0.3455005 7.938817e-05
# 27 Calanoida_unidentified    wp2  34 0.5326450 1.183073e-03

# 34               Copepoda regent 117 0.4274999 1.529892e-06
# 37                 Gel_FF regent  14 0.6077441 2.113796e-02
# 41              Calanoida regent 114 0.4309327 1.699927e-06



### Save both tables, next step will be plotting spatial patterns :-) 

### 19/06/2020: Modify some colnames for the hydro data
colnames(abund)[c(48,49)] <- c("Chla","bbp470")

setwd(WD)
write.table(abund, "table_abund+hydro_allnets_18_06_20.txt", sep = "\t")
write.table(esd, "table_ESD+hydro_allnets_18_06_20.txt", sep = "\t")


### 5°) Plotting time: map and plot zonal patterns  --------------------------------------------------------------------------
library("ggpubr")

# For tests
net <- "bongo"
data <- esd
var <- "ESD_Poecilostomatoida"

distrib.analyzer <- function(var, data) {
    
                    # Separate per dataset/net
                    for(net in unique(data$net)) {
                        
                        # Subset data
                        message(paste("Plotting zonal gradients of ",var," for ",net, sep = ""))
                        subset <- data[which(data$net == net),]
                        
                        # Drop the NAS if there are any
                        if( sum(is.na(subset[,var])) > 0 ) {
                                subset <- subset %>% drop_na(var)
                        } # eo if loop
                        
                        # unique(counts$group)
                        if(net == "bongo") {
                            counts <- counts.bongo
                        } else if(net == "wp2") {
                            counts <- counts.wp2
                        } else {
                            counts <- counts.regent
                        } # eo if else loop
                        ### Change some levels in counts so they match the colnames of abund/esd
                        levels(counts$group)[levels(counts$group) == "Calanoida_unid"] <- "Calanoida_unidentified"
                        levels(counts$group)[levels(counts$group) == "Copepoda_unid"] <- "Copepoda_unidentified"
                        levels(counts$group)[levels(counts$group) == "Crustaceans"] <- "Crust"
                        levels(counts$group)[levels(counts$group) == "Grazers"] <- "Small_grazers"
   
                        # Change varname for title/ axes caption
                        title <- str_replace(var,"Ab_","")
                        title1 <- paste(title,"\nmedian abundance", sep = "")
                        title2 <- paste("Fitted ",title," median abundance", sep = "")
                   
                        if(var == "Ab_Zooplankton") {
                            
                             stations2keep <- unique(subset$Station)
                             nstations <- length(stations2keep)
                             subset2 <- subset
                             
                        } else {
                            
                             stations2keep <- counts[counts$group == title & counts$n >= 20 ,"station"]
                             nstations <- length(stations2keep)
                             subset2 <- subset[which(subset$Station %in% stations2keep),]
                             
                        }                    
                     
                        ### Do we even have more than 30 stations?
                        if(nstations >= 30) {
                            
                            # Map 'em
                            map1 <- ggplot() + geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey85", colour = "grey50", size = 0.3) +
                                geom_point(aes(x = Longitude, y = Latitude, fill = get(var)), data = subset2, colour = "black", pch = 21, size = 2) + 
                                scale_fill_distiller(name = title1, palette = "Spectral") + 
                            	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(-180,-120,-60,0,60,120,180),
                                       	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                            	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                            	      	labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                            			panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") 
                        
                            # Fit gams if enough stations without NAs
                            if(nrow(subset2) >= 30) {
                                
                                require("mgcv")
                                gam1 <- mgcv::gam(data = subset2, get(var) ~ s(Latitude, bs="tp"), method = "REML")
                                # Extract deviance explained : str(summary(gam1))
                                r2.1 <- round(summary(gam1)$r.sq,3)
                                smooth.pval <- summary(gam1)$s.pv
                                if(smooth.pval < 0.001) {
                                    sign <- "***"
                                } else if(smooth.pval < 0.01) {
                                    sign <- "**"
                                } else if(smooth.pval < 0.05) {
                                    sign <- "*"
                                } else {
                                    sign <- "ns"
                                }
                                # Gather obs, y, fit and conf interval in a ddf
                                pred1 <- data.frame(y = subset2$Latitude, obs = subset2[,var], fit = predict(gam1,se.fit=T)$fit, se = predict(gam1,se.fit=T)$se.fit)
                                pred1 <- pred1[order(pred1$y),]
                                # Make zonal plots
                                plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = fit - se, xmax = fit + se), fill = "black", data = pred1, alpha = 0.25) +
                                            geom_path(aes(y = y, x = fit), data = pred1, colour = "black") +
                                            scale_y_continuous(position = "right", name = "Latitude (°)", limits = c(-65,80)) + 
                                            scale_x_continuous(name = paste(title2,"\n(r2 = ",r2.1,"; n = ",nstations,"; ",sign,")", sep = "")) + 
                                            theme_classic()   
                                
                                # Arrangethem in a panel
                                panel <- ggarrange(map1, plot1, ncol = 2, nrow = 1, widths = c(3,1))
                                setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Latitudinal_gradients/Imagery/",net, sep = ""))
                                ggsave(plot = panel, filename = paste("map_zonal_",var,"_",net,".jpg", sep = ""), dpi = 300, width = 9, height = 3.15)       
                                
                            } else {
                                
                                message(paste("NOT ENOUGH STATIONS", sep = ""))
                                
                            }
         
                            
                        } else {
                            
                            message(paste("NOT ENOUGH STATIONS", sep = ""))
                            
                        }
                                               
            } # eo for loop - net types
        
} # eo fun
# Apply fun above in for loop
variables <- colnames(abund)[c(5:42)]
for(var in variables) {
    distrib.analyzer(var = var, data = esd)   
}


### Modify the script a bit for abundances 
# For tests
net <- "bongo"
data <- abund
var <- "Ab_Zooplankton"
distrib.analyzer <- function(var, data) {
    
                    # Separate per dataset/net
                    for(net in unique(data$net)) {
                        
                        # Subset data
                        message(paste("Plotting zonal gradients of ",var," for ",net, sep = ""))
                        subset <- data[which(data$net == net),]
                        
                        # Drop the NAS if there are any
                        if( sum(is.na(subset[,var])) > 0 ) {
                                subset <- subset %>% drop_na(var)
                        } # eo if loop
                        
                        # Change varname for title/ axes caption
                        title <- str_replace(var,"Ab_","")
                        title1 <- paste(title,"\nabundance", sep = "")
                        title2 <- paste("Fitted ",title," abundance", sep = "")
                   
                        nstations <- nrow(subset)           
                     
                        ### Do we even have more than 30 stations?
                        if(nstations >= 30) {
                            
                            # Map 'em
                            map1 <- ggplot() + geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey85", colour = "grey50", size = 0.3) +
                                geom_point(aes(x = Longitude, y = Latitude, fill = get(var)), data = subset, colour = "black", pch = 21, size = 2) + 
                                scale_fill_distiller(name = title1, palette = "Spectral") + 
                            	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(-180,-120,-60,0,60,120,180),
                                       	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                            	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                            	      	labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                            			panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") 
                        
                            # Fit gams if enough stations without NAs
                            require("mgcv")
                            gam1 <- mgcv::gam(data = subset, get(var) ~ s(Latitude, bs="tp"), method = "REML")
                            # Extract deviance explained : str(summary(gam1))
                            r2.1 <- round(summary(gam1)$r.sq,3)
                            smooth.pval <- summary(gam1)$s.pv
                            if(smooth.pval < 0.001) {
                                    sign <- "***"
                            } else if(smooth.pval < 0.01) {
                                    sign <- "**"
                            } else if(smooth.pval < 0.05) {
                                    sign <- "*"
                            } else {
                                    sign <- "ns"
                            }
                            # Gather obs, y, fit and conf interval in a ddf
                            pred1 <- data.frame(y = subset$Latitude, obs = subset[,var], fit = predict(gam1,se.fit=T)$fit, se = predict(gam1,se.fit=T)$se.fit)
                            pred1 <- pred1[order(pred1$y),]
                            
                            # Make zonal plots
                            plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = fit - se, xmax = fit + se), fill = "black", data = pred1, alpha = 0.25) +
                                        geom_path(aes(y = y, x = fit), data = pred1, colour = "black") +
                                        scale_y_continuous(position = "right", name = "Latitude (°)", limits = c(-65,80)) + 
                                        scale_x_continuous(name = paste(title2,"\n(r2 = ",r2.1,"; n = ",nstations,"; ",sign,")", sep = "")) + 
                                        theme_classic()   
                                
                            # Arrangethem in a panel
                            panel <- ggarrange(map1, plot1, ncol = 2, nrow = 1, widths = c(3,1))
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Latitudinal_gradients/Imagery/",net, sep = ""))
                            ggsave(plot = panel, filename = paste("map_zonal_",var,"_",net,".jpg", sep = ""), dpi = 300, width = 9, height = 3.15)       
                                
                            } else {
                                
                                message(paste("NOT ENOUGH STATIONS", sep = ""))
                                
                            }
                                               
            } # eo for loop - net types
        
} # eo fun
# Apply fun above in for loop
variables <- colnames(abund)[c(5:42)]
for(var in variables) {
    distrib.analyzer(var = var, data = abund)   
}



### 5°) Correlograms: abund and ESD of chosen groups versus hydro data  --------------------------------------------------------------------------

abund <- read.table("table_abund+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
esd <- read.table("table_ESD+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
dim(abund) ; dim(esd)
# colnames(abund) --> 5:45 + 47,48 + 50:54 + 56:59
# colnames(esd) --> 5:45 + 47,48 + 50:54 + 56:59
net <- "wp2"

for(net in unique(abund$net)) {
    
        message(paste("", sep = ""))
        message(paste("Examining covariance between abundance and size for ",net, sep = ""))
        message(paste("", sep = ""))
        names <- colnames(abund)[c(5:21,23:27,29:37,43:45,47:59)] #; names
        abs <- abund[abund$net == net,names]
        
        # 2nd lapply based on the groups 
        groups <- colnames(abs)[c(1:31)]
        groups <- str_replace_all(groups,"Ab_","")
        colnames(abs)[c(1:31)] <- groups
        
        # Remove columns that ONLY have NAs (zoo groups)
        abs <- abs %>% select_if(~sum(!is.na(.)) > 0)
        mydata <- na.omit(abs)

        cormat <- round(cor(mydata, method = "spearman"),2)
        p.mat <- cor_pmat(mydata, method = "spearman", conf.level = 0.95)
        # dim(corr) ; dim(p.mat)
        # Re-order?
        cor_tri <- get_upper_tri(cormat)
        cor_tri <- melt(cor_tri, na.rm = T)
        pval_tri <- get_upper_tri(p.mat)
        pval_tri <- melt(pval_tri, na.rm = T)
        colnames(cor_tri)[3] <- "rho"
        colnames(pval_tri)[3] <- "pval"
        cor_tri$pval <- pval_tri$pval
        rm(pval_tri,p.mat,cormat) ; gc()
        #unique(cor_tri$Var1)
        #unique(cor_tri$Var2)
        
        ### First heatmap: main zoo groups versus hydro
        zoo_groups <- c("Zooplankton","Copepoda","Chaetognatha","Gel_carn","Gel_FF","Crust","Pteropoda","Rhizaria","Small_grazers")
        hydro_vars <- c("Temperature","O2","Salinity","MLD","Chla","bbp470","Micro","Nano","Pico","PAR2","NO2NO3","PO4","SiO2","PC1")
        cormat <- cor_tri[which(cor_tri$Var1 %in% zoo_groups & cor_tri$Var2 %in% hydro_vars),]
        cormat <- cormat[!(cormat$Var1 == cormat$Var2),]
        
        # Change factor levels for beter names
        levels(cormat$Var1)[levels(cormat$Var1) == "Gel_carn"] <- "Cnidaria"
        levels(cormat$Var1)[levels(cormat$Var1) == "Gel_FF"] <- "Tunicata"
        levels(cormat$Var1)[levels(cormat$Var1) == "Crust"] <- "Eumalacostraca"
        levels(cormat$Var1)[levels(cormat$Var1) == "Small_grazers"] <- "Ostracoda+Cladocera"
        # Same for hydrobio data
        levels(cormat$Var2)[levels(cormat$Var2) == "O2"] <- "Oxygen"
        levels(cormat$Var2)[levels(cormat$Var2) == "Chla"] <- "Chlorophylla"
        levels(cormat$Var2)[levels(cormat$Var2) == "Micro"] <- "%Micro"
        levels(cormat$Var2)[levels(cormat$Var2) == "Nano"] <- "%Nano"
        levels(cormat$Var2)[levels(cormat$Var2) == "Pico"] <- "%Pico"
        levels(cormat$Var2)[levels(cormat$Var2) == "PAR2"] <- "PAR"
        levels(cormat$Var2)[levels(cormat$Var2) == "NO2NO3"] <- "Nitrate"
        levels(cormat$Var2)[levels(cormat$Var2) == "PO4"] <- "Phosphate"
        levels(cormat$Var2)[levels(cormat$Var2) == "SiO2"] <- "Silicate"
         
        # dim(cormat)
        ### Add signif label
        cormat$signif <- NA
        cormat[cormat$pval <= 0.001,"signif"] <- "***"
        cormat[cormat$pval > 0.001 & cormat$pval <= 0.01,"signif"] <- "**"
        cormat[cormat$pval > 0.01 & cormat$pval < 0.05,"signif"] <- "*"
        cormat[cormat$pval >= 0.05,"signif"] <- "-"
       
        plot1 <- ggplot(cormat, aes(factor(Var2), factor(Var1), fill = rho)) + geom_tile(color = "white") +
            scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-0.7,0.7), 
                    name = paste("Spearman correlation\ncoefficient (",net,")", sep = "") ) + xlab("") + ylab("") + 
            theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
            axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))+
            coord_fixed() + geom_text(aes(factor(Var2), factor(Var1), label = signif), color = "black", size = 4)
        
        
        ### Second heatmaps: zoo groups versus zoo groups
        zoo_groups2 <- c("Zooplankton","Copepoda","Chaetognatha","Gel_carn","Gel_FF","Crust","Pteropoda","Rhizaria","Small_grazers",
                        "Copepoda_unid","Harpacticoida","Calanidae","Oncaeidae","Corycaeidae","Oithonidae",
                        "Candaciidae","Centropagidae","Eucalanidae","Euchaetidae","Heterorhabdidae","Metridinidae",
                        "Paracalanidae","Rhincalanidae","Sapphirinidae","Temoridae")
                        
        cormat <- cor_tri[which(cor_tri$Var1 %in% zoo_groups2 & cor_tri$Var2 %in% zoo_groups2),]
        cormat <- cormat[!(cormat$Var1 == cormat$Var2),]
        # Change factor levels for beter names
        levels(cormat$Var1)[levels(cormat$Var1) == "Gel_carn"] <- "Cnidaria"
        levels(cormat$Var1)[levels(cormat$Var1) == "Gel_FF"] <- "Tunicata"
        levels(cormat$Var1)[levels(cormat$Var1) == "Crust"] <- "Eumalacostraca"
        levels(cormat$Var2)[levels(cormat$Var2) == "Gel_carn"] <- "Cnidaria"
        levels(cormat$Var2)[levels(cormat$Var2) == "Gel_FF"] <- "Tunicata"
        levels(cormat$Var2)[levels(cormat$Var2) == "Crust"] <- "Eumalacostraca"
        levels(cormat$Var1)[levels(cormat$Var1) == "Small_grazers"] <- "Ostracoda+Cladocera"
        levels(cormat$Var2)[levels(cormat$Var2) == "Small_grazers"] <- "Ostracoda+Cladocera"
        
        cormat$signif <- NA
        cormat[cormat$pval <= 0.001,"signif"] <- "***"
        cormat[cormat$pval > 0.001 & cormat$pval <= 0.01,"signif"] <- "**"
        cormat[cormat$pval > 0.01 & cormat$pval < 0.05,"signif"] <- "*"
        cormat[cormat$pval >= 0.05,"signif"] <- "-"
        
        plot2 <- ggplot(cormat, aes(factor(Var1), factor(Var2), fill = rho)) + geom_tile(color = "white") +
            scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(max(cormat$rho)*-1,max(cormat$rho)), 
                    name = paste("Spearman correlation\ncoefficient (",net,")", sep = "") ) + xlab("") + ylab("") + 
            theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
            axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))+
            coord_fixed() + geom_text(aes(factor(Var1), factor(Var2), label = signif), color = "black", size = 4)
        
        # Save plots
        ggsave(plot = plot1, filename = paste("heatmap_corr_abund_hydro_",net,".jpg", sep = ""), dpi = 300, width = 9, height = 6)
        ggsave(plot = plot2, filename = paste("heatmap_corr_abund_only_",net,".jpg", sep = ""), dpi = 300, width = 11, height = 11)
        
} # eo for loop
    


