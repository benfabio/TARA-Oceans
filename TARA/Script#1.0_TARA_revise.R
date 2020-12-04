
##### 29/11/19: R script to perform the revisions needed for the TARA study on zooplankton size structure/abundance/diversity
### To do list:
#   - Examine the data distribution for the different net types and depth levels --> choose data transformation
#   - Examine spatial distribution of ESD/Abund/Div for the different net types and depth levels
#   - Perform spatial GAMs to model zonal patterns 
#   - Perform multivariate GAMs to model zooplankton and group patterns of ESD

### last update: 10/03/2020

### ---------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("vegan")
library("FactoMineR")
library("matrixStats")
library("lubridate")
library("viridis")
library("maps")
library("marmap")
library("mgcv")
library("nlme")
library("MuMIn")
library("ggrepel")

world <- map_data("world") # coastline for maps

WD <- getwd()

### ---------------------------------------------------------------------------------------------------------------------------

### 1°) Load the data, examine distribution, spatial patterns etc.
imagery <- read.csv("TARA_OCEANS_table_imagery_dataV2.txt", h = T, sep = ";")
metaB <- read.csv("TARA_OCEANS_table_rel_abund.csv", h = T, sep = ";")
shannon <- read.csv("TARA_OCEANS_table_shannon.csv", h = T, sep = ";")
dim(imagery); dim(metaB); dim(shannon)
colnames(imagery)[3] <- "net"
colnames(imagery)[42] <- "ESD_Sapphirinidae"

# str(imagery); str(metaB); str(shannon)
unique(imagery$station) # 166
unique(metaB$station) # 77
unique(shannon$station) # 92
# Check colnames
colnames(imagery)
# ESD: 9-16 and then 25 to 45
# Abund: 17-24 and then 46-66
colnames(metaB) # abund = 6-13
colnames(shannon) # shannon from Ibarlbaz --> 8-17

### 19/12/19: Map the spatial distribution of the stations providing the imagery data for each net 
bathy <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = 90, resolution = 15)
bathy <- fortify(bathy) # convert to ddf
# Plot in a for loop
# n <- "wp2"
for(n in unique(imagery$net)) {
    
        sub <- imagery[imagery$net == n,]
        
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = z), data = bathy[bathy$z <= 0,]) + 
            scale_fill_distiller(name = "Depth (m)", palette = "Greys") + 
            geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey90", colour = "grey75", size = 0.3) +
        	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(-180,-120,-60,0,60,120,180),
                   	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
        	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
        	      	labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
        	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        			panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) + 
            geom_point(aes(x = x, y = y), data = sub, colour = "black", pch = 21, fill = "#fd8d3c")
            
        # Save
        ggsave(plot = map, filename = paste("map_stations_",n,".pdf", sep = ""), dpi = 300, height = 6, width = 8) 
                    
} # eo for loop

### A) Examine the distribution of the variables of interest. Maybe transform.
### A.1) Imagery data (watchout: need to seperate per ner type)
# imagery$net # ok
variables <- colnames(imagery)[9:94]
# For testing fun
# data <- imagery
# var <- "Abund_Protista"
# net <- "wp2"
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
                        
                            # In multinet samples, some station are repeated
                            if( net %in% c("multinet_surf","multinet_meso") ) {
                                
                                subset2 <- data.frame(subset %>% group_by(station) %>% summarize(median = median(get(var),na.rm=T)) )
                               
                                # perform Shapiro-Wilk normality test
                                norm.test <- shapiro.test(subset2$median) # if p-val < 0.05, var is not normally distrbuted
                                pval <- norm.test$p.value
                                if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                                # Plot original distribution, and then try some transformations: sqrt(), log10(), log(), cube root (^(1/3))
                                setwd(paste(WD,"/","Distribution_plots/imagery/",net,sep=""))
                                
                                plot <- ggplot(subset2, aes(x = median)) + geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) +
                                        geom_density(alpha=.2, fill="#dd3497") + xlab(var) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                        theme_classic()
                                        
                                ggsave(plot = plot, filename = paste("plot_distrib_raw_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                                
                                # sqrt transformation 
                                norm.test <- shapiro.test(sqrt(subset2$median))
                                pval <- norm.test$p.value
                                if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                                plot <- ggplot(subset2, aes(x = sqrt(median))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (sqrt)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                                ggsave(plot = plot, filename = paste("plot_distrib_sqrt_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                            
                                # log transformation 
                                norm.test <- shapiro.test(log(subset2$median+0.0000001))
                                pval <- norm.test$p.value
                                if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                                plot <- ggplot(subset2, aes(x = log(median))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                                ggsave(plot = plot, filename = paste("plot_distrib_log_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                           
                                # log10 transformation 
                                norm.test <- shapiro.test(log10(subset2$median+0.0000001))
                                pval <- norm.test$p.value
                                if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                                plot <- ggplot(subset2, aes(x = log10(median))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                                ggsave(plot = plot, filename = paste("plot_distrib_log10_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                                
                                # cubic transformation 
                                norm.test <- shapiro.test((subset2$median)^(1/3))
                                pval <- norm.test$p.value
                                if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                                plot <- ggplot(subset2, aes(x = (median)^(1/3))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (cubic)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                                ggsave(plot = plot, filename = paste("plot_distrib_cubic_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                           
                            
                            } else {
                                
                                # perform Shapiro-Wilk normality test
                                norm.test <- shapiro.test(subset[,var]) # if p-val < 0.05, var is not normally distrbuted
                                pval <- norm.test$p.value
                                if( pval < 0.05 ) { norm <- "Not normal" } else { norm <- "Normal" }
                                # Plot original distribution, and then try some transformations: sqrt(), log10(), log(), cube root (^(1/3))
                                setwd(paste(WD,"/","Distribution_plots/imagery/",net,sep=""))
                                plot <- ggplot(subset, aes(x = get(var))) + geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) +
                                        geom_density(alpha=.2, fill="#dd3497") + xlab(var) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                        theme_classic()
                                        
                                ggsave(plot = plot, filename = paste("plot_distrib_raw_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                                
                                # sqrt transformation 
                                norm.test <- shapiro.test(sqrt(subset[,var]))
                                pval <- norm.test$p.value
                                if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                                plot <- ggplot(subset, aes(x = sqrt(get(var)))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (sqrt)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                                ggsave(plot = plot, filename = paste("plot_distrib_sqrt_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                            
                                # log transformation 
                                norm.test <- shapiro.test(log(subset[,var]+0.0000001))
                                pval <- norm.test$p.value
                                if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                                plot <- ggplot(subset, aes(x = log(get(var)))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                                ggsave(plot = plot, filename = paste("plot_distrib_log_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                           
                                # log10 transformation 
                                norm.test <- shapiro.test(log10(subset[,var]+0.0000001))
                                pval <- norm.test$p.value
                                if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                                plot <- ggplot(subset, aes(x = log10(get(var)))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                                ggsave(plot = plot, filename = paste("plot_distrib_log10_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                                
                                # cubic transformation 
                                norm.test <- shapiro.test((subset[,var])^(1/3))
                                pval <- norm.test$p.value
                                if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
                                plot <- ggplot(subset, aes(x = (get(var))^(1/3))) + 
                                            geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                                            xlab(paste(var, " (cubic)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                                            theme_classic()
                                            
                                ggsave(plot = plot, filename = paste("plot_distrib_cubic_",var,"_",net,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
                                
                            } # eo if else loop - for summarizing values across depths ranges in multinet samples
                            
                            
                        } # eo if else loop - if no values 
                        
                        
            } # eo for loop - net types
        
} # eo fun
# Apply fun above in for loop
for(var in variables) {
    distrib.analyzer(var = var, data = imagery)   
}

distrib.analyzer(var = "PAR2", data = imagery)   

### Conclusion 
### For WP2: 
# - log transform ESD
# - log OR cubic transform abundances (try both)
# - keep RAW % of phytoplankton types
# - min O2 depth -> sqrt or cubic
# - log transform DCM
# - keep RAW z_eu
# - keep raw Temp at 10m, same for T_zeu and T_mld
# - same for Salinity (S), don't transform
# - ignore Density because same as Salinity
# - keep raw O2 for all versions
# - cubic transform Chla at 10m/ Zeu/ MLD 
# - log transform or cubic transform integrated Chla
# - cubic or log transform macronutrients (NO3/PO4/SiO2) and Fluorescence
# - log transform backscattering
# - cubic transform PAR
# - RAW PAR2
imagery[imagery$net == "wp2","PAR"]
imagery[imagery$net == "wp2","PAR2"]

### For Bongo: Apply the same transformations
### !!! issue wit PAR values --> NOT the same scale as in WP2 --> units??
imagery[imagery$net == "bongo","PAR"]
imagery[imagery$net == "bongo","PAR2"]
# Not enough data --> DO NOT CONSIDER PAR, use the values from WP2
# USE PAR2

### For Regent: Apply the same transformations
### !!! issue wit PAR values --> NOT the same scale as in WP2 --> units??
imagery[imagery$net == "regent","PAR"]
imagery[imagery$net == "regent","PAR2"]
# Not enough data --> DO NOT CONSIDER PAR, use the values from WP2
# USE PAR2

### For multinet_surf and multinet_meso --> apply the same transformations

### A.2) metaB data: examine only for the 180-2000 size fraction (no need to use a FUN)  ------------------------------------------
variables <- colnames(metaB)[c(6:13)]
# var <- "abund_Copepoda"
for(var in variables) {
        
        # Useless message
        message(paste("Examining the distribution of ", var, sep = ""))
        subset <- metaB[which(metaB$fraction_size == "180-2000"),]
        
        # Sometimes, the net dataset has no value (or not enough) for the variable of interest -> length(!is.na(subset$var)) == 0    
        if( sum(!is.na(subset[,var])) < 5 ) {
        
            message(paste("\nNo values for ",var," || moving to next variable", sep = ""))
            
        } else {
         
            # perform Shapiro-Wilk normality test
            norm.test <- shapiro.test(subset[,var]) # if p-val < 0.05, var is not normally distrbuted
            pval <- norm.test$p.value
            if( pval < 0.05 ) { norm <- "Not normal" } else { norm <- "Normal" }
            # Plot original distribution, and then try some transformations: sqrt(), log10(), log(), cube root (^(1/3))
            setwd(paste(WD,"/","Distribution_plots/metaB", sep = ""))
            plot <- ggplot(subset, aes(x = get(var))) + geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) +
                    geom_density(alpha=.2, fill="#dd3497") + xlab(var) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                    theme_classic()
                    
            ggsave(plot = plot, filename = paste("plot_distrib_raw_",var,"_","180-2000","_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
            
            # sqrt transformation 
            norm.test <- shapiro.test(sqrt(subset[,var]))
            pval <- norm.test$p.value
            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
            plot <- ggplot(subset, aes(x = sqrt(get(var)))) + 
                        geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                        xlab(paste(var, " (sqrt)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                        theme_classic()
                        
            ggsave(plot = plot, filename = paste("plot_distrib_sqrt_",var,"_","180-2000","_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
        
            # log transformation 
            norm.test <- shapiro.test(log(subset[,var]+0.0000001))
            pval <- norm.test$p.value
            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
            plot <- ggplot(subset, aes(x = log(get(var)))) + 
                        geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                        xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                        theme_classic()
                        
            ggsave(plot = plot, filename = paste("plot_distrib_log_",var,"_","180-2000","_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
       
            # log10 transformation 
            norm.test <- shapiro.test(log10(subset[,var]+0.0000001))
            pval <- norm.test$p.value
            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
            plot <- ggplot(subset, aes(x = log10(get(var)))) + 
                        geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                        xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                        theme_classic()
                        
            ggsave(plot = plot, filename = paste("plot_distrib_log10_",var,"_","180-2000","_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
            
            # cubic transformation 
            norm.test <- shapiro.test((subset[,var])^(1/3))
            pval <- norm.test$p.value
            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
            plot <- ggplot(subset, aes(x = (get(var))^(1/3))) + 
                        geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                        xlab(paste(var, " (cubic)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                        theme_classic()
                        
            ggsave(plot = plot, filename = paste("plot_distrib_cubic_",var,"_","180-2000","_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
            
            
        } # eo if else loop - if no values
    
} # eo for loop - for each variables
### for the relative abund based on metaB, cubic transform seems the best, but remember these are RELATIVE abundances already
### Pteropoda, Chaetognatha, and Decapoda display VERY LOW contribution


### A.3) Div indices from metaB (shannon): need to separate per size fraction and layers  ------------------------------------------
variables <- colnames(shannon)[c(8:17)]
# var <- "shannon_copepoda"
for(var in variables) {
        
        # Useless message
        message(paste("Examining the distribution of ", var, sep = ""))
        
        # Sometimes, the net dataset has no value (or not enough) for the variable of interest -> length(!is.na(subset$var)) == 0    
        if( sum(!is.na(shannon[,var])) < 5 ) {
        
            message(paste("\nNo values for ",var," || moving to next variable", sep = ""))
            
        } else {
         
            # perform Shapiro-Wilk normality test
            norm.test <- shapiro.test(shannon[,var]) # if p-val < 0.05, var is not normally distrbuted
            pval <- norm.test$p.value
            if( pval < 0.05 ) { norm <- "Not normal" } else { norm <- "Normal" }
            # Plot original distribution, and then try some transformations: sqrt(), log10(), log(), cube root (^(1/3))
            setwd(paste(WD,"/","Distribution_plots/Shannons", sep = ""))
            plot <- ggplot(shannon, aes(x = get(var))) + geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) +
                    geom_density(alpha=.2, fill="#dd3497") + xlab(var) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                    theme_classic()
                    
            ggsave(plot = plot, filename = paste("plot_distrib_raw_",var,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
            
            # sqrt transformation 
            norm.test <- shapiro.test(sqrt(shannon[,var]))
            pval <- norm.test$p.value
            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
            plot <- ggplot(shannon, aes(x = sqrt(get(var)))) + 
                        geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                        xlab(paste(var, " (sqrt)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                        theme_classic()
                        
            ggsave(plot = plot, filename = paste("plot_distrib_sqrt_",var,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
        
            # log transformation 
            norm.test <- shapiro.test(log(shannon[,var]+0.0000001))
            pval <- norm.test$p.value
            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
            plot <- ggplot(shannon, aes(x = log(get(var)))) + 
                        geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                        xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                        theme_classic()
                        
            ggsave(plot = plot, filename = paste("plot_distrib_log_",var,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
       
            # log10 transformation 
            norm.test <- shapiro.test(log10(shannon[,var]+0.0000001))
            pval <- norm.test$p.value
            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
            plot <- ggplot(shannon, aes(x = log10(get(var)))) + 
                        geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                        xlab(paste(var, " (log)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                        theme_classic()
                        
            ggsave(plot = plot, filename = paste("plot_distrib_log10_",var,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
            
            # cubic transformation 
            norm.test <- shapiro.test((shannon[,var])^(1/3))
            pval <- norm.test$p.value
            if(pval < 0.05) { norm <- "Not normal" } else { norm <- "Normal" }
            plot <- ggplot(shannon, aes(x = (get(var))^(1/3))) + 
                        geom_histogram(aes(y=..density..), colour="black", fill="grey55", bins = 25) + geom_density(alpha=.2, fill="#dd3497") + 
                        xlab(paste(var, " (cubic)",sep="")) + ylab(paste("Density (",norm,")\n",pval,sep="")) + 
                        theme_classic()
                        
            ggsave(plot = plot, filename = paste("plot_distrib_cubic_",var,"_",norm,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
            
            
        } # eo if else loop - if no values
    
}
### For Shannon indices --> keep the raw values


### B) 19/12/19: Examine the zonal latitudinal gradients of the variables of interest (ESD and Abund) by fitting a GAM
### (like Ibarbalz et al.) rather than the loess fit
variables <- colnames(imagery)[9:66] ; variables
var <- "ESD_Centropagidae"
net <- "bongo"
data <- imagery
### Make function to map the transfromed data and to fit basic GAMs and plot the fit vs. the observations 
distrib.analyzer <- function(var, data) {
    
                    # Separate per dataset/net
                    for(net in unique(data$net)) {
                        
                        # Subset data
                        message(paste("Plotting zonal gradients of ",var," for ",net, sep = ""))
                        subset <- data[which(data$net == net),]
        
                        # Sometimes, the net dataset has no value (or not enough) for the variable of interest -> length(!is.na(subset$var)) == 0    
                        if( sum(!is.na(subset[,var])) < 5 ) {
                        
                            message(paste("\nNo values for ",var," for ",net," || moving to next variable", sep = ""))
                            
                        } else {
                        
                            # In multinet samples, some station are repeated
                            if( net %in% c("multinet_surf","multinet_meso") ) {
                                
                                subset2 <- data.frame(subset %>% group_by(station) %>% summarize(median = median(get(var),na.rm=T), y = mean(y), x = mean(x)) )
                                subset2$logtrans <- log1p(subset2$median)
                                subset2$cubic <- (subset2$median)^(1/3)
                                # Drop the NAS if there are any
                                if( sum(is.na(subset2[,"logtrans"])) > 0 ) {
                                        subset2 <- subset2 %>% drop_na(logtrans,cubic)
                                } # eo if loop
                                
                                # Map 'em
                                map1 <- ggplot() +
                                    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey90", colour = "grey75", size = 0.3) +
                                    geom_point(aes(x = x, y = y, fill = logtrans), data = subset2, colour = "black", pch = 21) + 
                                    scale_fill_distiller(name = var, palette = "RdYlBu") + 
                                	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(-180,-120,-60,0,60,120,180),
                                           	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                                	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                                	      	labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                                	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                                			panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 
                                    
                                #
                                map2 <- ggplot() + 
                                    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey90", colour = "grey75", size = 0.3) +
                                    geom_point(aes(x = x, y = y, fill = cubic), data = subset2, colour = "black", pch = 21) +
                                    scale_fill_distiller(name = var, palette = "RdYlBu") + 
                                	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(-180,-120,-60,0,60,120,180),
                                           	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                                	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                                	      	labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                                	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                                			panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                                
                                # Save 'em
                                setwd(paste(WD,"/","Latitudinal_gradients/imagery/maps/", net, sep = ""))
                                ggsave(plot = map1, filename = paste("map_",var,"_","log","_",net,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
                                ggsave(plot = map2, filename = paste("map_",var,"_","cubic","_",net,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
                                
                                # Fit gams
                                require("mgcv")
                                gam1 <- mgcv::gam(data = subset2, logtrans ~ s(y,bs="tp"), method = "REML")
                                gam2 <- mgcv::gam(data = subset2, cubic ~ s(y,bs="tp"), method = "REML")
                                # Extract deviance explained : str(summary(gam1))
                                r2.1 <- round(summary(gam1)$r.sq,3)
                                r2.2 <- round(summary(gam2)$r.sq,3)
                                # Gather ovs, y, fit and conf interval in a ddf
                                pred1 <- data.frame(y = subset2$y, obs = subset2$logtrans, fit = predict(gam1,se.fit=T)$fit, se = predict(gam1,se.fit=T)$se.fit)
                                pred2 <- data.frame(y = subset2$y, obs = subset2$cubic, fit = predict(gam2,se.fit=T)$fit, se = predict(gam2,se.fit=T)$se.fit)
                                # Make zonal plots
                                plot1 <- ggplot(data = pred1) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
                                        geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") + 
                                        geom_line(aes(x = y, y = fit), colour = "#0868ac") + 
                                        xlab("Latitude") + ylab(paste(var," (log)",sep="")) + theme_classic() +
                                        ggtitle(paste("GAM fit (R2 = ",r2.1,")", sep=""))

                                plot2 <- ggplot(data = pred2) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
                                        geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") + 
                                        geom_line(aes(x = y, y = fit), colour = "#0868ac") + 
                                        xlab("Latitude") + ylab(paste(var," (cubic)",sep="")) + theme_classic() +
                                        ggtitle(paste("GAM fit (R2 = ",r2.2,")", sep=""))
                                        
                                # Save 'em
                                setwd(paste(WD,"/","Latitudinal_gradients/imagery/",net,sep=""))
                                ggsave(plot = plot1, filename = paste("plot_zonal_",var,"_","log","_",net,".pdf", sep = ""), dpi = 300, width = 7, height = 4)
                                ggsave(plot = plot2, filename = paste("plot_zonal_",var,"_","cubic","_",net,".pdf", sep = ""), dpi = 300, width = 7, height = 4)
                            
                            } else {
                                
                                # Log and cubic transform variable
                                subset$logtrans <- log1p(subset[,var])
                                subset$cubic <- (subset[,var])^(1/3)
                                
                                # Map 'em first
                                map1 <- ggplot() + 
                                    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey90", colour = "grey75", size = 0.3) +
                                    geom_point(aes(x = x, y = y, fill = logtrans), data = subset, colour = "black", pch = 21) + 
                                    scale_fill_distiller(name = var, palette = "RdYlBu") + 
                                	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(-180,-120,-60,0,60,120,180),
                                           	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                                	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                                	      	labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                                	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                                			panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
                                    
                                #
                                map2 <- ggplot() + 
                                    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey90", colour = "grey75", size = 0.3) +
                                    geom_point(aes(x = x, y = y, fill = cubic), data = subset, colour = "black", pch = 21) + 
                                    scale_fill_distiller(name = var, palette = "RdYlBu") + 
                                	coord_quickmap() + scale_x_continuous(name = "Longitude", breaks = c(-180,-120,-60,0,60,120,180),
                                           	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                                	scale_y_continuous(name = "Latitude", breaks = c(-90,-60,-30,0,30,60,90),
                                	      	labels = c("-90°N","-60°N","-30°N","0°N","30°N","60°N","90°N"), expand = c(0,0)) +
                                	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                                			panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) 
                                
                               # Save 'em
                               setwd(paste(WD,"/","Latitudinal_gradients/imagery/maps/", net, sep = ""))
                               ggsave(plot = map1, filename = paste("map_",var,"_","log","_",net,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
                               ggsave(plot = map2, filename = paste("map_",var,"_","cubic","_",net,".pdf", sep = ""), dpi = 300, width = 8, height = 5)
                           
                               ### Sometimes there are NAs in the ESD data (no organisms so no size measurements so no data, unlike abundances)
                                if( sum(is.na(subset[,"logtrans"])) > 0 ) {
                                    
                                    # Drop the NAS
                                    subset2 <- subset %>% drop_na(logtrans,cubic)
                                    
                                    # Fit gams with y
                                    require("mgcv")
                                    gam1 <- mgcv::gam(data = subset2, logtrans ~ s(y,bs="tp"), method = "REML")
                                    gam2 <- mgcv::gam(data = subset2, cubic ~ s(y,bs="tp"), method = "REML")
                                    # Extract deviance explained : str(summary(gam1))
                                    r2.1 <- round(summary(gam1)$r.sq,3)
                                    r2.2 <- round(summary(gam2)$r.sq,3)
                                    # Gather ovs, y, fit and conf interval in a ddf
                                    pred1 <- data.frame(y = subset2$y, obs = subset2$logtrans, fit = predict(gam1,se.fit=T)$fit, se = predict(gam1,se.fit=T)$se.fit)
                                    pred2 <- data.frame(y = subset2$y, obs = subset2$cubic, fit = predict(gam2,se.fit=T)$fit, se = predict(gam2,se.fit=T)$se.fit)
                                    # Make zonal plots
                                    plot1 <- ggplot(data = pred1) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
                                            geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") + 
                                            geom_line(aes(x = y, y = fit), colour = "#0868ac") + 
                                            xlab("Latitude") + ylab(paste(var," (log)",sep = "")) + theme_classic() +
                                            ggtitle(paste("GAM fit (R2 = ",r2.1,")", sep = ""))

                                    plot2 <- ggplot(data = pred2) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
                                            geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") + 
                                            geom_line(aes(x = y, y = fit), colour = "#0868ac") + 
                                            xlab("Latitude") + ylab(paste(var," (cubic)", sep = "")) + theme_classic() +
                                            ggtitle(paste("GAM fit (R2 = ",r2.2,")", sep = ""))
                                        
                                    # Save 'em
                                    setwd(paste(WD,"/","Latitudinal_gradients/imagery/",net,sep=""))
                                    ggsave(plot = plot1, filename = paste("plot_zonal_",var,"_","log","_",net,".pdf", sep = ""), dpi = 300, width = 4.5, height = 3)
                                    ggsave(plot = plot2, filename = paste("plot_zonal_",var,"_","cubic","_",net,".pdf", sep = ""), dpi = 300, width = 4.5, height = 3)
                                    
                                } else {
                                    
                                    # Fit gams with y
                                    require("mgcv")
                                    gam1 <- mgcv::gam(data = subset, logtrans ~ s(y,bs="tp"), method = "REML")
                                    gam2 <- mgcv::gam(data = subset, cubic ~ s(y,bs="tp"), method = "REML")
                                    # Extract deviance explained : str(summary(gam1))
                                    r2.1 <- round(summary(gam1)$r.sq,3)
                                    r2.2 <- round(summary(gam2)$r.sq,3)
                                    # Gather ovs, y, fit and conf interval in a ddf
                                    pred1 <- data.frame(y = subset$y, obs = subset$logtrans, fit = predict(gam1,se.fit=T)$fit, se = predict(gam1,se.fit=T)$se.fit)
                                    pred2 <- data.frame(y = subset$y, obs = subset$cubic, fit = predict(gam2,se.fit=T)$fit, se = predict(gam2,se.fit=T)$se.fit)
                                    # Make zonal plots
                                    plot1 <- ggplot(data = pred1) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
                                            geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") + 
                                            geom_line(aes(x = y, y = fit), colour = "#0868ac") + 
                                            xlab("Latitude") + ylab(paste(var," (log)",sep="")) + theme_classic() +
                                            ggtitle(paste("GAM fit (R2 = ",r2.1,")", sep=""))

                                    plot2 <- ggplot(data = pred2) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
                                            geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") + 
                                            geom_line(aes(x = y, y = fit), colour = "#0868ac") + 
                                            xlab("Latitude") + ylab(paste(var," (cubic)",sep="")) + theme_classic() +
                                            ggtitle(paste("GAM fit (R2 = ",r2.2,")", sep=""))
                                        
                                    # Save 'em
                                    setwd(paste(WD,"/","Latitudinal_gradients/imagery/",net,sep=""))
                                    ggsave(plot = plot1, filename = paste("plot_zonal_",var,"_","log","_",net,".pdf", sep = ""), dpi = 300, width = 4.5, height = 3)
                                    ggsave(plot = plot2, filename = paste("plot_zonal_",var,"_","cubic","_",net,".pdf", sep = ""), dpi = 300, width = 4.5, height = 3)
                     
                                } # eo if else loop - if there are NAs in ESD data

                            
                            } # eo if else loop - for summarizing values across depths ranges in multinet samples
                            
                            
                        } # eo if else loop - if no values 
                        
                        
            } # eo for loop - net types
        
} # eo fun
# Apply fun above in for loop
for(var in variables) {
    distrib.analyzer(var = var, data = imagery)   
}

### Conclusions: 
# For Régent net: for abundances and ESD -> log transformation (try cubic for abund later)
# For WP2 net: for abundances and ESD -> log transformation (try cubic for abund later)
# For Bongo net: for abundances and ESD -> log transformation (try cubic for abund later)
### Same for the multinet samples

# ### 19/12/19: Try some random forests out instead of GAMs
# library("ranger")
# # ?ranger
#
# ### First, examine evolutionof r2 with ntrees (for each ntree value, run 100 rf)
# trees <- seq(from = 10, to = 1025, by = 25)
# var <- "ESD_Zooplankton"
# net <- "bongo"
# # nt <- 133
# res <- lapply(trees, function(nt) {
#
#                 # Subset data
#                 message(paste("Training 100 RF with ",nt," trees for ",var," sampled with ",net, sep = ""))
#                 subset <- data[which(data$net == net),]
#
#                 # Log and cubic transform variable
#                 subset$logtrans <- log1p(subset[,var])
#                 subset$cubic <- (subset[,var])^(1/3)
#
#                 # Prepare empty vectors
#                 R2_log <- NA
#                 R2_cubic <- NA
#                 MSE_log <- NA
#                 MSE_cubic <- NA
#
#                 # Fill the empty vectors with a for loop
#                 for(i in c(1:100)) {
#
#                     # Fit random forests
#                     require("ranger")
#                     rf1 <- ranger(data = subset, logtrans ~ y, importance = "none",
#                                 write.forest = F, replace = T, splitrule = "variance",
#                                 oob.error = T, num.trees = nt)
#                     rf2 <- ranger(data = subset, cubic ~ y, importance = "none",
#                                 write.forest = F, replace = T, splitrule = "variance",
#                                 oob.error = T, num.trees = nt)
#
#                     # Extract deviance explained : rf1$r.squared
#                     r2.1 <- round(rf1$r.squared,5)
#                     r2.2 <- round(rf2$r.squared,5)
#                     mse.1 <- round(rf1$prediction.error,5)
#                     mse.2 <- round(rf2$prediction.error,5)
#
#                     # Supply the R2 values to their matching vector
#                     R2_log[i] <- r2.1
#                     R2_cubic[i] <- r2.2
#                     MSE_log[i] <- mse.1
#                     MSE_cubic[i] <- mse.2
#
#                 } # eo for loop - i in 1:100
#                 R2_log[R2_log < 0] <- 0.0000001
#                 R2_cubic[R2_cubic < 0] <- 0.0000001
#                 # Join in a single table that you will return
#                 table <- data.frame(ntree = nt, R2_log = R2_log, R2_cubic = R2_cubic, MSE_log = MSE_log, MSE_cubic = MSE_cubic)
#                 return(table)
#
#         } # eo FUN
#
# ) # eo lapply
#
# ### Bind and examine relationships bewteen thr RF models' R2 and their number of trees
# ddf <- bind_rows(res)
# dim(ddf) ; str(ddf)
# rm(res); gc()
# # Compute mean and sdev R2
# ddf.2 <- data.frame(ddf %>%
#             group_by(ntree) %>%
#             summarize(R2.1 = mean(R2_log), R2.2 = mean(R2_cubic),
#             sd.1 = sd(R2_log), sd.2 = sd(R2_cubic),
#             MSE.1 = mean(MSE_log), MSE.2 = mean(MSE_cubic),
#             sd.3 = sd(MSE_log), sd.4 = sd(MSE_cubic)
#         )
# ) # eo ddf
# ddf.2
#
# # Plot profiles of R2 and MSE
# quartz()
# ggplot(data = ddf.2) + geom_point(aes(x = ntree, y = R2.2), colour = "#4eb3d3") +
#         geom_ribbon(aes(x = ntree, ymin = R2.2-sd.2, ymax = R2.2+sd.2), alpha = 0.2, fill = "black") +
#         geom_line(aes(x = ntree, y = R2.2), colour = "#0868ac") +
#         xlab("Number of trees") + ylab(paste("R2"," (cubic)",sep="")) + theme_classic()
# #
# quartz()
# ggplot(data = ddf.2) + geom_point(aes(x = ntree, y = R2.1), colour = "#4eb3d3") +
#         geom_ribbon(aes(x = ntree, ymin = R2.1-sd.1, ymax = R2.1+sd.1), alpha = 0.2, fill = "black") +
#         geom_line(aes(x = ntree, y = R2.1), colour = "#0868ac") +
#         xlab("Number of trees") + ylab(paste("R2"," (log)",sep="")) + theme_classic()
# #
# quartz()
# ggplot(data = ddf.2) + geom_point(aes(x = ntree, y = MSE.1), colour = "#4eb3d3") +
#         geom_ribbon(aes(x = ntree, ymin = MSE.1-sd.3, ymax = MSE.1+sd.3), alpha = 0.2, fill = "black") +
#         geom_line(aes(x = ntree, y = MSE.1), colour = "#0868ac") +
#         xlab("Number of trees") + ylab(paste("MSE"," (log)",sep="")) + theme_classic()
# #
# quartz()
# ggplot(data = ddf.2) + geom_point(aes(x = ntree, y = MSE.2), colour = "#4eb3d3") +
#         geom_ribbon(aes(x = ntree, ymin = MSE.2-sd.4, ymax = MSE.2+sd.4), alpha = 0.2, fill = "black") +
#         geom_line(aes(x = ntree, y = MSE.2), colour = "#0868ac") +
#         xlab("Number of trees") + ylab(paste("MSE"," (cubic)",sep="")) + theme_classic()
#
# # Choose ntree == 750
#
# ### Write a FUN to fit latitudinal RF models based on the re abobe
# var <- "ESD_Zooplankton"
# net <- "bongo"
# data <- imagery
#
# rf.fitter <- function(var, data) {
#
#                     # Separate per dataset/net
#                     for(net in unique(data$net)) {
#
#                         # Subset data
#                         message(paste("Plotting zonal gradients of ",var," for ",net, sep = ""))
#                         subset <- data[which(data$net == net),]
#
#                         # Sometimes, the net dataset has no value (or not enough) for the variable of interest -> length(!is.na(subset$var)) == 0
#                         if( sum(!is.na(subset[,var])) < 5 ) {
#
#                             message(paste("\nNo values for ",var," for ",net," || moving to next variable", sep = ""))
#
#                         } else {
#
#                             # In multinet samples, some station are repeated
#                             if( net %in% c("multinet_surf","multinet_meso") ) {
#
#                                 subset2 <- data.frame(subset %>% group_by(station) %>% summarize(median = median(get(var),na.rm=T), y = mean(y), x = mean(x)) )
#                                 subset2$logtrans <- log1p(subset2$median)
#                                 subset2$cubic <- (subset2$median)^(1/3)
#                                 # Drop the NAS if there are any
#                                 if( sum(is.na(subset2[,"logtrans"])) > 0 ) {
#                                         subset2 <- subset2 %>% drop_na(logtrans,cubic)
#                                 } # eo if loop
#
#                                 # Fit random forests
#                                 require("ranger")
#                                 rf1 <- ranger(data = subset2, logtrans ~ y, importance = "none", keep.inbag = T,
#                                             write.forest = T, replace = T, splitrule = "variance",
#                                             oob.error = T, num.trees = 750)
#                                 rf2 <- ranger(data = subset2, cubic ~ y, importance = "none", keep.inbag = T,
#                                             write.forest = T, replace = T, splitrule = "variance",
#                                             oob.error = T, num.trees = 750)
#
#                                 # Extract deviance explained : rf1$r.squared
#                                 r2.1 <- round(rf1$r.squared, 3)
#                                 r2.2 <- round(rf2$r.squared, 3)
#
#                                 # Gather obs, y, fit and conf interval in a ddf
#                                 pred1 <- data.frame(y = subset2$y, obs = subset2$logtrans,
#                                             fit = rf1$predictions, se = predict(rf1, data = subset2, type = "se")$se)
#
#                                 pred2 <- data.frame(y = subset2$y, obs = subset2$cubic,
#                                             fit = rf2$predictions, se = predict(rf2, data = subset2, type = "se")$se)
#
#                                 # Make zonal plots
#                                 plot1 <- ggplot(data = pred1) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
#                                         geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") +
#                                         geom_line(aes(x = y, y = fit), colour = "#0868ac") +
#                                         xlab("Latitude") + ylab(paste(var," (log)", sep = "")) + theme_classic() +
#                                         ggtitle(paste("RF fit (R2 = ",r2.1,")", sep = ""))
#
#                                 plot2 <- ggplot(data = pred2) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
#                                         geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") +
#                                         geom_line(aes(x = y, y = fit), colour = "#0868ac") +
#                                         xlab("Latitude") + ylab(paste(var," (cubic)", sep = "")) + theme_classic() +
#                                         ggtitle(paste("RF fit (R2 = ",r2.2,")", sep = ""))
#
#                                 # Save 'em
#                                 setwd(paste(WD,"/","Latitudinal_gradients/Imagery/random.forests/",net,sep=""))
#                                 ggsave(plot = plot1, filename = paste("plot_zonal_RF_",var,"_","log","_",net,".pdf", sep = ""), dpi = 300, width = 7, height = 4)
#                                 ggsave(plot = plot2, filename = paste("plot_zonal_RF_",var,"_","cubic","_",net,".pdf", sep = ""), dpi = 300, width = 7, height = 4)
#
#                             } else {
#
#                                # Log and cubic transform variable
#                                subset$logtrans <- log1p(subset[,var])
#                                subset$cubic <- (subset[,var])^(1/3)
#
#                                ### Sometimes there are NAs in the ESD data (no organisms so no size measurements so no data, unlike abundances)
#                                if( sum(is.na(subset[,"logtrans"])) > 0 ) {
#
#                                     # Drop the NAS
#                                     subset2 <- subset %>% drop_na(logtrans,cubic)
#
#                                     # Fit random forests
#                                     require("ranger")
#                                     rf1 <- ranger(data = subset2, logtrans ~ y, importance = "none", keep.inbag = T,
#                                                 write.forest = T, replace = T, splitrule = "variance",
#                                                 oob.error = T, num.trees = 750)
#                                     rf2 <- ranger(data = subset2, cubic ~ y, importance = "none", keep.inbag = T,
#                                                 write.forest = T, replace = T, splitrule = "variance",
#                                                 oob.error = T, num.trees = 750)
#
#                                     # Extract deviance explained : rf1$r.squared
#                                     r2.1 <- round(rf1$r.squared, 3)
#                                     r2.2 <- round(rf2$r.squared, 3)
#
#                                     # Gather obs, y, fit and conf interval in a ddf
#                                     pred1 <- data.frame(y = subset2$y, obs = subset2$logtrans,
#                                                 fit = rf1$predictions, se = predict(rf1, data = subset2, type = "se")$se)
#
#                                     pred2 <- data.frame(y = subset2$y, obs = subset2$cubic,
#                                                 fit = rf2$predictions, se = predict(rf2, data = subset2, type = "se")$se)
#
#                                     # Make zonal plots
#                                     plot1 <- ggplot(data = pred1) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
#                                             geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") +
#                                             geom_line(aes(x = y, y = fit), colour = "#0868ac") +
#                                             xlab("Latitude") + ylab(paste(var," (log)", sep = "")) + theme_classic() +
#                                             ggtitle(paste("RF fit (R2 = ",r2.1,")", sep = ""))
#
#                                     plot2 <- ggplot(data = pred2) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
#                                             geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") +
#                                             geom_line(aes(x = y, y = fit), colour = "#0868ac") +
#                                             xlab("Latitude") + ylab(paste(var," (cubic)", sep = "")) + theme_classic() +
#                                             ggtitle(paste("RF fit (R2 = ",r2.2,")", sep = ""))
#                                     # Save 'em
#                                     setwd(paste(WD,"/","Latitudinal_gradients/Imagery/random.forests/",net,sep=""))
#                                     ggsave(plot = plot1, filename = paste("plot_zonal_RF_",var,"_","log","_",net,".pdf", sep = ""), dpi = 300, width = 4.5, height = 3)
#                                     ggsave(plot = plot2, filename = paste("plot_zonal_RF_",var,"_","cubic","_",net,".pdf", sep = ""), dpi = 300, width = 4.5, height = 3)
#
#                                 } else {
#
#                                     # Fit random forests
#                                     require("ranger")
#                                     rf1 <- ranger(data = subset, logtrans ~ y, importance = "none", keep.inbag = T,
#                                                 write.forest = T, replace = T, splitrule = "variance",
#                                                 oob.error = T, num.trees = 750)
#                                     rf2 <- ranger(data = subset, cubic ~ y, importance = "none", keep.inbag = T,
#                                                 write.forest = T, replace = T, splitrule = "variance",
#                                                 oob.error = T, num.trees = 750)
#
#                                     # Extract deviance explained : rf1$r.squared
#                                     r2.1 <- round(rf1$r.squared, 3)
#                                     r2.2 <- round(rf2$r.squared, 3)
#
#                                     # Gather obs, y, fit and conf interval in a ddf
#                                     pred1 <- data.frame(y = subset$y, obs = subset$logtrans,
#                                                 fit = rf1$predictions, se = predict(rf1, data = subset, type = "se")$se)
#
#                                     pred2 <- data.frame(y = subset$y, obs = subset$cubic,
#                                                 fit = rf2$predictions, se = predict(rf2, data = subset, type = "se")$se)
#
#                                     # Make zonal plots
#                                     plot1 <- ggplot(data = pred1) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
#                                             geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") +
#                                             geom_line(aes(x = y, y = fit), colour = "#0868ac") +
#                                             xlab("Latitude") + ylab(paste(var," (log)", sep = "")) + theme_classic() +
#                                             ggtitle(paste("RF fit (R2 = ",r2.1,")", sep = ""))
#
#                                     plot2 <- ggplot(data = pred2) + geom_point(aes(x = y, y = obs), colour = "#4eb3d3") +
#                                             geom_ribbon(aes(x = y, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black") +
#                                             geom_line(aes(x = y, y = fit), colour = "#0868ac") +
#                                             xlab("Latitude") + ylab(paste(var," (cubic)", sep = "")) + theme_classic() +
#                                             ggtitle(paste("RF fit (R2 = ",r2.2,")", sep = ""))
#
#                                     # Save 'em
#                                     setwd(paste(WD,"/","Latitudinal_gradients/Imagery/random.forests/",net,sep=""))
#                                     ggsave(plot = plot1, filename = paste("plot_zonal_RF_",var,"_","log","_",net,".pdf", sep = ""), dpi = 300, width = 4.5, height = 3)
#                                     ggsave(plot = plot2, filename = paste("plot_zonal_RF_",var,"_","cubic","_",net,".pdf", sep = ""), dpi = 300, width = 4.5, height = 3)
#
#                                 } # eo if else loop - if there are NAs in ESD data
#
#
#                             } # eo if else loop - for summarizing values across depths ranges in multinet samples
#
#
#                         } # eo if else loop - if no values
#
#
#             } # eo for loop - net types
#
# } # eo fun
# # Apply fun above in for loop
# for(var in variables) {
#     rf.fitter(var = var, data = imagery)
# }

### --> that was a bad idea all in all ^^


### ---------------------------------------------------------------------------------------------------------------------------

### 19/12/19: Prior to the multivariate modelling, examine collinearity/ coeff of correlations between responses
### (logESD and logABund) and predictors (env + shannons)

### First, apply a log1p transformation to all ESD and ABund data
imagery2 <- imagery
variables <- colnames(imagery2)[9:66]
imagery2[variables] <- lapply(imagery2[variables], log1p)

### Second, apply chosen transformations on env predictors (code from 18/12/19)
colnames(imagery2)
# - min O2 depth -> sqrt or cubic
imagery2$min_o2_depth <- sqrt(imagery2$min_o2_depth)
# - log transform DCM
imagery2$DCM <- log(imagery2$DCM)
# - cubic transform Chla at 10m/ Zeu/ MLD 
imagery2$Chla_10m <- (imagery2$Chla_10m)^(1/3)
imagery2$Chla_zeu <- (imagery2$Chla_zeu)^(1/3)
imagery2$Chla_mld <- (imagery2$Chla_mld)^(1/3)
# - log transform or cubic transform integrated Chla
imagery2$Chla <- log(imagery2$Chla)
# - cubic or log transform macronutrients (NO3/PO4/SiO2) and Fluorescence
imagery2$PO4 <- log(imagery2$PO4)
imagery2$NO3 <- log(imagery2$NO3)
imagery2$Si <- log(imagery2$Si)
imagery2$Fluo <- log(imagery2$Fluo)
# - log transform backscattering
imagery2$bac660 <- log1p(imagery2$bac660)
# - cubic transform PAR
imagery2$PAR <- (imagery2$PAR)^(1/3)

### Second, provide the shannon index values for additional predictors
colnames(shannon) # cols 8:17
names <- colnames(shannon)[c(8:17)]
imagery2[names] <- NA
# Fill with for loop
# s <- 7
for(s in unique(shannon$station) ) {
    
    for(v in names) {
        
        message(paste("Supplying ",v," for station ",s, sep = ""))
        value <- shannon[shannon$station == s,v]
        
        if(length(value) > 1) {
            value[!is.na(value)]
        } else {
            imagery2[imagery2$station == s,v] <- value
        } # eo if else loop
        
    } # eo 2nd for loop
    
} # eo 2nd for loop
### Check results
#summary(imagery2[names])
# Cool


### 20/12/19: To examine relationships between responses and predictors, plot biplot with basic gam functions with geom_smooth
# colnames(imagery2)
nets <- unique(imagery2$net) # nets
responses <- colnames(imagery2)[c(9:66)]
preds1 <- colnames(imagery2)[c(67:78,82:93,95)]
preds2 <- colnames(imagery2)[c(96:104)]
# Examine which Temp and Sal, and Chla versions have most NaNs
length(is.na(imagery2$T_10m))
length(is.na(imagery2$T_zeu))
length(is.na(imagery2$T_mld))
# All Temp versions have equal nb of values
length(is.na(imagery2$S_10m))
length(is.na(imagery2$S_zeu))
length(is.na(imagery2$S_mld))
# All Sal versions have equal nb of values
length(is.na(imagery2$Chla_10m))
length(is.na(imagery2$Chla_zeu))
length(is.na(imagery2$Chla_mld))
length(is.na(imagery2$Chla))
# Same 

### For each net, examine resp x preds relationships
# For testing
n <- "wp2"
r <- "ESD_Copepoda"
p <- "shannon_protist_parasit"

for(n in nets) {
    
    subset <- imagery2[imagery2$net == n,]
    #if( n %in% c("bongo","regent","multinet_surf","multinet_meso") ) {
        #preds1 <- preds1[!(preds1 == "PAR")]
        #} # eo if loop - remove PAR if not WP2
    
    for(r in responses) {
        
        #for(p in preds2) {
        
        for(p in "PAR2") {
            
            message(paste("Plotting relationship between ",r, " and ",p," for ",n, sep = ""))
            data.temp <- na.omit(subset[,c(r,p)])
            #if(p %in% c("S_10m","S_zeu","S_mld")) {
                #subset <- subset[subset[,p] > 30,]
                #} # eo if loop - remove obs < 30
            
            #if(p == "bac660") {
                #subset <- subset[subset[,p] < 1,]
                #} # eo if loop - remove bac660 > 1
            
            # Check if enough data 
            if( nrow(data.temp) >= 10 ) {
                
                setwd(paste(WD,"/","Biplots and correlations/",n, sep = ""))
                plot <- ggplot(data = data.temp) + geom_point(aes(x = get(p), y = get(r)), colour = "#4eb3d3") + 
                        geom_smooth(aes(x = get(p), y = get(r)), method = mgcv::gam, colour = "#0868ac") +
                        xlab(p) + ylab(r) + theme_classic()
                # Save 
                ggsave(plot = plot, filename = paste("plot_smooth_gam_",n,"_",r,"_",p,".pdf", sep = ""), dpi = 300, width = 5, height = 4)
                        
            } else {
                
                message(paste("Not enough data for plotting the relationship between ",r," and ",p," || ",nrow(data.temp), sep = ""))
                rm(data.temp); gc()
                
            }
            
        } # eo for loop - pred in preds
        
    } # eo for loop - resp in reponses
    
} # eo for loop - n in nets
# Re-run with preds2 too

### Prior to running the actual GAMs, examine collinearity between env predictors to make a first selection (again, per nets)
# n <- "multinet_meso"
for(n in nets) {
    
        subset <- imagery2[imagery2$net == n,]
        message(paste("Examining predictors collinearity for net = ",n, sep = ""))
    
        #if( n %in% c("bongo","regent","multinet_surf","multinet_meso") ) {
            #preds1 <- preds1[!(preds1 == "PAR")]
            #} # eo if loop - remove PAR if not WP2
    
        data.temp <- na.omit(subset[,preds1])
        cormat <- round(cor(data.temp, method = "spearman"),3)
        melted_cormat <- melt(cormat)
        
        get_lower_tri <- function(cormat) {
          cormat[upper.tri(cormat)] <- NA
          return(cormat)
        }
        get_upper_tri <- function(cormat) {
          cormat[lower.tri(cormat)]<- NA
          return(cormat)
        }
        reorder_cormat <- function(cormat) {
        # Use correlation between variables as distance
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd)
        cormat <-cormat[hc$order, hc$order]
        }
    
        # Create a ggheatmap
        ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
                     scale_fill_gradient2(low = "#4575b4", high = "#d73027", mid = "white", midpoint = 0,
                             limit = c(-1,1), space = "Lab", name = "Spearman\nCorrelation", guide=FALSE) +
                     theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+ 
                     coord_fixed() + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
                     theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
                             panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank(),
                             legend.justification = c(1,0),legend.position = c(0.6,0.7),legend.direction = "horizontal")
        setwd(paste(WD,"/","Biplots and correlations/", sep = ""))
        ggsave(plot = ggheatmap, filename = paste("heatmap_env_preds_",n,"v2.pdf", sep = ""), dpi = 300, width = 15, height = 15)
    
} # eo for loop - n in nets

### Conclusions:
### A) WP2: 
# pick one Temp version out of the 3 versions available (the most normally distributed) -> run alternative models #
# Same with Salinity 
# pick one O2 version out of the 3 (the most normally distributed) -> run alternative models 
# pick one CHLA version out of the 4 (the most normally distributed) -> run alternative models 
# need have to choose O2 and Temp though. Run GAMs with the two, separately
# choose between DCM and Zeu (choose most normally distributed) -> DCM excludes Chla @ 10m...choose Zeu raw. Also, DCM overall more correlated to CHLA than Zeu
# choose between Pico and Micro --> Micro because less correlated to SST and O2
# choose between Fluo and backscaretting --> backscaretting
# choose between the 3 macronutrients (Si, PO4, NO3) -> run alternative models 

### Final choice of predictors: Temperature(x3) or O2(x3) + Sal (x3) + CHLA (x4) + Zeu + Micro + Nano + bac660 + macronutrients (Si or PO4 or NO3) + min_o2_depth 

### B) Bongo: 
# Same as above but do not have to choose between Fluo and bac660
### Final choice of predictors: Temperature(x3) or O2(x3) + Sal(x3) + CHLA (x4) + Zeu + Micro + + Fluo + Nano + bac660 + macronutrients (Si or PO4 or NO3) + min_o2_depth 

### C) Regent: same choices as Bongo !
### Temperature(x3) or O2(x3) + Sal(x3) + CHLA (x4) + Zeu + Micro + + Fluo + Nano + bac660 + macronutrients (Si or PO4 or NO3) + min_o2_depth 

### NOTE: make sure to choose the versions of CHLA/SAL/TEMP accordingly (don"t mix Temp at 10m with CHLA from MLD)


# Write a function that will fit GAMs based on these sets of predictors. Use thin plate regression, default k values,
# REML for smoothing parameter estimation method, select = T to add an extra penalty to each term so that it can be penalized to zero
# preds1

### 10/03/2020: Add PAR2 to all predictors list ! 

### Define lists of predictors set for each of the 3 net types (2 lists because bongo and regent have the same)
list.pred.bongo <- list(
        c("T_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","Si","Fluo","bac660","PAR2"), 
        c("T_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","Fluo","bac660","PAR2"), 
        c("T_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","Fluo","bac660","PAR2"), 
        c("T_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","Si","Fluo","bac660","PAR2"), 
        c("T_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","Fluo","bac660","PAR2"), 
        c("T_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","Fluo","bac660","PAR2"), 
        c("T_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","Si","Fluo","bac660","PAR2"), 
        c("T_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","Fluo","bac660","PAR2"), 
        c("T_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","Fluo","bac660","PAR2"), 
        c("o2_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","Si","Fluo","bac660","PAR2"), 
        c("o2_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","Fluo","bac660","PAR2"), 
        c("o2_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","Fluo","bac660","PAR2"), 
        c("o2_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","Si","Fluo","bac660","PAR2"), 
        c("o2_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","Fluo","bac660","PAR2"), 
        c("o2_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","Fluo","bac660","PAR2"), 
        c("o2_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","Si","Fluo","bac660","PAR2"), 
        c("o2_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","Fluo","bac660","PAR2"), 
        c("o2_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","Fluo","bac660","PAR2")
) # eo list.pred.bongo
list.pred.wp2 <- list(
        c("T_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","Si","bac660","PAR2"), 
        c("T_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","bac660","PAR2"), 
        c("T_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","bac660","PAR2"), 
        c("T_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","Si","bac660","PAR2"), 
        c("T_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","bac660","PAR2"), 
        c("T_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","bac660","PAR2"), 
        c("T_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","Si","bac660","PAR2"), 
        c("T_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","bac660","PAR2"), 
        c("T_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","bac660","PAR2"), 
        c("o2_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","Si","bac660","PAR2"), 
        c("o2_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","bac660","PAR2"), 
        c("o2_10m","S_10m","Chla_10m","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","bac660","PAR2"), 
        c("o2_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","Si","bac660","PAR2"), 
        c("o2_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","bac660","PAR2"), 
        c("o2_zeu","S_zeu","Chla_zeu","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","bac660","PAR2"), 
        c("o2_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","Si","bac660","PAR2"), 
        c("o2_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","PO4","bac660","PAR2"), 
        c("o2_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","bac660","PAR2")
) # eo list.pred.bongo
# v <- "Abund_Corycaeidae"
# net <- "bongo"
# i <- 4
gam.fitter <- function(net, variables) {
    
                setwd(paste(WD,"/","Models/GAMs/",net, sep = ""))
                message(paste("\nTraining GAMs for net = ",net, sep = ""))
                subset <- imagery2[imagery2$net == net,]
                
                # Choose list of predictor sets
                if(net == "wp2") {
                    list.preds <- list.pred.wp2
                    npred <- 10
                } else {
                    list.preds <- list.pred.bongo
                    npred <- 11
                } # eo if else loop
                
                # Run GAMs for each element of this list
                require("parallel")
                # i <- 1
                mclapply(c(1:length(list.preds)), function(i) {
                    
                        preds <- list.preds[[i]]
                        formula <- paste(preds, collapse = "+")
                        
                        # Train GAM model for each variable
                        for(v in variables) {
                            
                            if( npred == 10 ) {
                                
                                subset2 <- na.omit(subset[,c(v,preds)])
                                nobs <- nrow(subset2)
                                
                                # Control for a minimum of 15 obs, or if the resp var is a flat one (abundances always = 0)
                                if( nobs >= 25 & sum(subset2[,v]) > 0 ) {
                                    
                                    # The number of coeff in the model should lower than the nb of obs -> nobs*9
                                    message(paste("Training GAM with list ",i,"  | ",v, " ~ ",formula, sep = ""))
                                    k <- round(nobs/npred, digits = 0)
                                    model <- mgcv::gam(get(v) ~ s(get(preds[1]),bs="tp",k=k)+s(get(preds[2]),bs="tp",k=k)+s(get(preds[3]),bs="tp",k=k)+
                                                s(get(preds[4]),bs="tp",k=k)+s(get(preds[5]),bs="tp",k=k)+s(get(preds[6]),bs="tp",k=k)+
                                                s(get(preds[7]),bs="tp",k=k)+s(get(preds[8]),bs="tp",k=k)+s(get(preds[9]),bs="tp",k=k)+
                                                s(get(preds[10]),bs="tp",k=k), data = subset2, select = T, method = "REML")
                                            
                                    save(model, file = paste("gam_",v,"_",formula,"_",net,".Rdata", sep = "") )      
                                    
                                } else {
                                    
                                    message(paste("Stopped because not enough data || n = ",nobs, sep = ""))
                                    
                                } # eo if else loop - for controlling nobs     
                                
                            } else if( npred == 11 ) {
                                
                                subset2 <- na.omit(subset[,c(v,preds)])
                                nobs <- nrow(subset2)
                                
                                # Control for a minimum of 15 obs
                                if( nobs >= 25 & sum(subset2[,v]) > 0 ) {
                                    
                                    # The number of coeff in the model should lower than the nb of obs -> nobs*10
                                    message(paste("Training GAM with list ",i,"  | ",v, " ~ ",formula, sep = ""))
                                    k <- round(nobs/npred, digits = 0)
                                    model <- mgcv::gam(get(v) ~ s(get(preds[1]),bs="tp",k=k)+s(get(preds[2]),bs="tp",k=k)+s(get(preds[3]),bs="tp",k=k)+
                                                s(get(preds[4]),bs="tp",k=k)+s(get(preds[5]),bs="tp",k=k)+s(get(preds[6]),bs="tp",k=k)+
                                                s(get(preds[7]),bs="tp",k=k)+s(get(preds[8]),bs="tp",k=k)+s(get(preds[9]),bs="tp",k=k)+
                                                s(get(preds[10]),bs="tp",k=k)+s(get(preds[11]),bs="tp",k=k), data = subset2, select = T, method = "REML")
                                                
                                    save(model, file = paste("gam_",v,"_",formula,"_",net,".Rdata", sep = "") )         
                                    
                                } else {
                                    
                                    message(paste("Stopped because not enough data || n = ",nobs, sep = ""))
                                    
                                } # eo if else loop - for controlling nobs
                                        
                            } # eo if else loop - npred
                            
                            
                        } # eo for v in variables
                    
                    }, mc.cores = 2
                
                ) # eo mclapply
    
} # eo gam.fitter FUN
# Apply to each 3 net
gam.fitter(net = "wp2", variables = responses)
gam.fitter(net = "bongo", variables = responses)
gam.fitter(net = "regent", variables = responses)

### Next, examine models' results per group, response var (ESD or abund), examine models' skills etc.

### 21/12/19: Do the same with random forests !
library("ranger")

# First, run some tests based on a couple of variables to choose ntrees
preds <- c("o2_mld","S_mld","Chla_mld","x_Micro","x_Nano","min_o2_depth","z_eu","NO3","bac660")
var <- "Abund_Cnidaria"
net <- "wp2"
ntrees <- seq(from = 10, to = 500, by = 10) ; ntrees
# Train 100 RF models per ntree values
# t <- 50
res <- lapply(ntrees, function(t) {
    
            message(paste("Running 100 RF models with ntrees = ", t, sep = ""))
            subset <- imagery2[imagery2$net == net,]
            subset2 <- na.omit(subset[,c(var,preds)])
            # nobs <- nrow(subset2)
            r2 <- NA
            mse <- NA
            for(i in c(1:100)) {
             
                rf <- ranger(data = subset2, as.formula(paste(var,"~",paste(preds,collapse="+"), sep = "")),
                            importance = "impurity", keep.inbag = F, write.forest = F,
                            replace = T, splitrule = "variance", oob.error = T, num.trees = t)
                # rf$variable.importance
                # Return explained deviance/ variance
                r2[i] <- round(rf$r.squared,3)
                mse[i] <- round(rf$prediction.error,3)
                
            } # eo for loop - i in 1-100
            
            table <- data.frame(ntree = t, R2 = mean(r2), sd1 = sd(r2), MSE = mean(mse), sd2 = sd(mse))
            return(table)
                    
    } # eo lapply
) # eo lapply
# Rbind
ddf <- do.call(rbind, res)
dim(ddf); head(ddf)
rm(res); gc()
quartz()
ggplot(ddf) + geom_ribbon(aes(x = ntree, ymin = R2-sd1, ymax = R2+sd1), alpha = 0.2, fill = "black") +
    geom_line(aes(x = ntree, y = R2), colour = "black") +
    xlab("Number of trees") + ylab("Explained variance (%)") + theme_classic()
#
quartz()
ggplot(ddf) + geom_ribbon(aes(x = ntree, ymin = MSE-sd2, ymax = MSE+sd2), alpha = 0.2, fill = "black") +
    geom_line(aes(x = ntree, y = MSE), colour = "black") +
    xlab("Number of trees") + ylab("MSE") + theme_classic()

### --> choose ntree = 200

v <- "Abund_Cnidaria"
net <- "wp2"

### Write fun to train RF models 
rf.fitter <- function(net, variables) {
    
                setwd(paste(WD,"/","Models/RF/",net, sep = ""))
                message(paste("\nTraining RF models for net = ",net, sep = ""))
                subset <- imagery2[imagery2$net == net,]
                
                # Choose list of predictor sets
                if(net == "wp2") {
                    list.preds <- list.pred.wp2
                    npred <- 10
                } else {
                    list.preds <- list.pred.bongo
                    npred <- 11
                } # eo if else loop
                
                # Run GAMs for each element of this list
                require("parallel")
                # i <- 1
                mclapply(c(1:length(list.preds)), function(i) {
                    
                        preds <- list.preds[[i]]
                        formula <- paste(preds, collapse = "+")
                        
                        # Train GAM model for each variable
                        for(v in variables) {
                            
                            if( npred == 10 ) {
                                
                                subset2 <- na.omit(subset[,c(v,preds)])
                                nobs <- nrow(subset2)
                                
                                # Control for a minimum of 15 obs, or if the resp var is a flat one (abundances always = 0)
                                if( nobs >= 25 & sum(subset2[,v]) > 0 ) {
                                    
                                    # The number of coeff in the model should lower than the nb of obs -> nobs*9
                                    message(paste("Training RF with list ",i,"  | ",v, " ~ ",formula, sep = ""))

                                    model <- ranger(data = subset2, as.formula(paste(v,"~",paste(preds,collapse="+"), sep = "")),
                                            importance = "impurity", keep.inbag = T, write.forest = T,
                                            replace = T, splitrule = "variance", oob.error = T, num.trees = 200)
                                            
                                    save(model, file = paste("rf_",v,"_",formula,"_",net,".Rdata", sep = "") )      
                                    
                                } else {
                                    
                                    message(paste("Stopped because not enough data || n = ",nobs, sep = ""))
                                    
                                } # eo if else loop - for controlling nobs     
                                
                            } else if( npred == 11 ) {
                                
                                subset2 <- na.omit(subset[,c(v,preds)])
                                nobs <- nrow(subset2)
                                
                                # Control for a minimum of 15 obs
                                if( nobs >= 25 & sum(subset2[,v]) > 0 ) {
                                    
                                    # The number of coeff in the model should lower than the nb of obs -> nobs*10
                                    message(paste("Training RF with list ",i,"  | ",v, " ~ ",formula, sep = ""))
                                    
                                    model <- ranger(data = subset2, as.formula(paste(v,"~",paste(preds,collapse="+"), sep = "")),
                                            importance = "impurity", keep.inbag = T, write.forest = T,
                                            replace = T, splitrule = "variance", oob.error = T, num.trees = 200)
                                            
                                    save(model, file = paste("rf_",v,"_",formula,"_",net,".Rdata", sep = "") )         
                                    
                                } else {
                                    
                                    message(paste("Stopped because not enough data || n = ",nobs, sep = ""))
                                    
                                } # eo if else loop - for controlling nobs
                                        
                            } # eo if else loop - npred
                            
                            
                        } # eo for v in variables
                    
                    }, mc.cores = 2
                
                ) # eo mclapply
    
} # eo rf.fitter FUN

# Apply to each 3 net
rf.fitter(net = "wp2", variables = responses)
rf.fitter(net = "bongo", variables = responses)
rf.fitter(net = "regent", variables = responses)


### ---------------------------------------------------------------------------------------------------------------------------

### 06/01/20: Examine outputs of the GAMs and RF models:
# - for which reponse var, which is the best model combination
# - which are the best explaining variables (importance in RF, p-values in GAMs)
# - then, plot response curves

### First, examine the structure of the R objects
# setwd(paste(WD,"/","Models/GAMs", sep = ""))
# net <- "wp2"
# setwd(paste(WD,"/","Models/GAMs/",net, sep = ""))
# model <- get(load( dir()[1] ))
# str(model)
# summary(model)
# str( summary(model) )

### For each net, retrieve the results of the GAMs
nets <- c("wp2","bongo","regent")
n <- "bongo"
require("parallel")
res <- mclapply(nets, function(n) {
    
                setwd(paste(WD,"/","Models/GAMs/",n, sep = ""))
                files <- dir()[grep("gam_",dir())]
                # f <- "gam_ESD_Zooplankton_T_zeu+S_zeu+Chla_zeu+x_Micro+x_Nano+min_o2_depth+z_eu+Si+Fluo+bac660+PAR2_bongo.Rdata"
                gams <- lapply(files, function(f) {
                            
                            message(paste("Retrieving GAM results for ",f, sep = ""))
                            model <- get(load(f))
                            # Extract the terms of the model from the file name
                            terms <- do.call(cbind, strsplit(as.character(f),"\\+"))
                            # 1st element contains response var and first predictor
                            terms[1,] <- str_replace_all(terms[1,], "gam_", "")
                            # unlist(strsplit(terms[1,],"_"))
                            # Watchout, if other_Copepoda is the group, then there's an extra underscore to consider
                            if( grepl("other_Copepoda", terms[1,]) ) {
                                
                                resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2],unlist(strsplit(terms[1,],"_"))[3],sep = "_")
                                pred1 <- paste(unlist(strsplit(terms[1,],"_"))[4], unlist(strsplit(terms[1,],"_"))[5], sep = "_")
                                # And extract the last pred from the last element of terms
                                predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
                                
                            } else {
                                
                                resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2], sep = "_")
                                pred1 <- paste(unlist(strsplit(terms[1,],"_"))[3], unlist(strsplit(terms[1,],"_"))[4], sep = "_")
                                # And extract the last pred from the last element of terms
                                predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
                                
                            } # eo if else loop - grepl("other_Copepoda", terms[1,]) 
                                
                            if(n == "wp2") {
                                preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],terms[9,],predlast)
                            } else {  
                                preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],terms[9,],terms[10,],predlast)    
                            } # eo if else loop
                            
                            # Extract r2 and AIC of the GAM model
                            aic <- round(model$aic, 4)
                            r2 <- round(summary(model)$r.sq, 4)
                            formula <- paste(preds, collapse = "+")
                            pvalues <- summary(model)$s.pv
                            Fstats <- summary(model)$s.table[,3]
                            
                            table <- data.frame(resp = resp, preds = preds, net = n, formula = formula, R2 = r2, AIC = aic, pval = pvalues, F = Fstats)
                            # Rank predictors based on their F scores
                            table$rank <- table[,"F"]/ max(table[,"F"])
                            
                            rownames(table) <- NULL
                            return(table)
                            
                    } # eo FUN
                    
                ) # eo lapply
                
                # Rbind 
                table <- do.call(rbind, gams)
                return(table)
                rm(gams);gc()
    
        }, mc.cores = 3
        
) # eo mclapply
# Rbind 
table <- bind_rows(res)
head(table); dim(table)
summary(table)
rm(res); gc()

### Examine the distribution of R2 scores for each response variable, separate abundances from ESD
sd( table[table$resp == "Abund_Augaptilidae","R2"] )
sd( table[table$resp == "Abund_Augaptilidae","AIC"] )

scores <- data.frame(table %>% group_by(resp,net) %>% summarize(mean.R2 = mean(R2), sd.R2 = sd(R2), mean.AIC = mean(AIC), sd.AIC = sd(AIC)) )
scores$type <- NA
scores[c(1:67),"type"] <- "Abundance"
scores[c(68:117),"type"] <- "ESD"

# Plot distrbution of R2 and AIC by displaying mean R2/AIC and associated sd
quartz()
ggplot(scores[scores$type == "ESD",], aes(x = factor(resp), y = mean.R2)) +
    geom_errorbar(width = 0.1, aes(ymin = mean.R2-sd.R2, ymax = mean.R2+sd.R2) ) +
    geom_point(shape=21, size=3, aes(fill = factor(net))) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("R2") + xlab("") + coord_flip() + 
    facet_wrap(factor(scores[scores$type == "ESD","net"]), nrow = 3)

quartz()
ggplot(scores[scores$type == "Abundance",], aes(x = factor(resp), y = mean.R2)) +
    geom_errorbar(width = 0.1, aes(ymin = mean.R2-sd.R2, ymax = mean.R2+sd.R2) ) +
    geom_point(shape=21, size=3, aes(fill = factor(net))) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("R2") + xlab("") + coord_flip() + 
    facet_wrap(factor(scores[scores$type == "Abundance","net"]), nrow = 3)


# Same with AIC
quartz()
ggplot(scores[scores$type == "ESD",], aes(x = factor(resp), y = mean.AIC)) +
    geom_errorbar(width = 0.1, aes(ymin = mean.AIC-sd.AIC, ymax = mean.AIC+sd.AIC) ) +
    geom_point(shape=21, size=3, aes(fill = factor(net))) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("AIC") + xlab("") + coord_flip() + 
    facet_wrap(factor(scores[scores$type == "ESD","net"]), nrow = 3)

quartz()
ggplot(scores[scores$type == "Abundance",], aes(x = factor(resp), y = mean.AIC)) +
    geom_errorbar(width = 0.1, aes(ymin = mean.AIC-sd.AIC, ymax = mean.AIC+sd.AIC) ) +
    geom_point(shape=21, size=3, aes(fill = factor(net))) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("AIC") + xlab("") + coord_flip() + 
    facet_wrap(factor(scores[scores$type == "Abundance","net"]), nrow = 3)


### Examine, for each response variable, the best models
best.models <- data.frame(table %>% group_by(net,resp,formula) %>% summarize(R2 = unique(R2), AIC = unique(AIC)) )
head(best.models)
quartz()
ggplot(data = best.models[best.models$resp == "ESD_Zooplankton",]) +
    geom_point(aes(x = R2, y = AIC, fill = factor(net)), shape = 21, size = 3) +
    theme_bw() + facet_wrap(factor(best.models[best.models$resp == "ESD_Zooplankton","net"]), nrow = 3)
    
# The relationships between R2 and AIC is unclear: for the Bongo dataset, AIC actually increase with R2. Most of the time, for a similar AIC, you can have models spanning 0.4 units of R2. 
# Is your goal model parsimony or the predictive power of the model? If parsimony, then use AIC, if predictive power then 𝑅2.
# Usually the answer is similar, but if you are comparing models with very similar 𝑅2 or a number of low quality predictors the answers can be different.

### For each response variable, examine which type if variable is associated with the lowest AIC/ highest R2 between: 
# - T_10m vs T_zeu vs T_mld
# - S_10m vs S_zeu vs S_mld
# - o2_10m vs o2_zeu vs o2_mld
# - Chla_10m vs Chla_zeu vs Chla_mld
# - Si vs NO3 vs PO4
# - o2 vs Temp

### A) T_10m vs T_zeu vs T_mld
quartz()
ggplot(aes(x = factor(preds), y = R2), data = table[table$preds %in% c("T_10m","T_zeu","T_mld"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("R2") +
    facet_wrap(factor(table[table$preds %in% c("T_10m","T_zeu","T_mld"),"net"]), nrow = 3, scales = "free")
# Unclear
quartz()
ggplot(aes(x = factor(preds), y = AIC), data = table[table$preds %in% c("T_10m","T_zeu","T_mld"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("AIC") +
    facet_wrap(factor(table[table$preds %in% c("T_10m","T_zeu","T_mld"),"net"]), nrow = 3, scales = "free")
# Still unclear, look at F scores or p-values
quartz()
ggplot(aes(x = factor(preds), y = log(get("pval"))), data = table[table$preds %in% c("T_10m","T_zeu","T_mld"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = F) + theme_classic() + xlab("") + ylab("F scores") +
    facet_wrap(factor(table[table$preds %in% c("T_10m","T_zeu","T_mld"),"net"]), nrow = 3)

### B) S_10m vs S_zeu vs S_mld
quartz()
ggplot(aes(x = factor(preds), y = R2), data = table[table$preds %in% c("S_10m","S_zeu","S_mld"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("R2") +
    facet_wrap(factor(table[table$preds %in% c("S_10m","S_zeu","S_mld"),"net"]), nrow = 3, scales = "free")
# Unclear
quartz()
ggplot(aes(x = factor(preds), y = AIC), data = table[table$preds %in% c("S_10m","S_zeu","S_mld"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("AIC") +
    facet_wrap(factor(table[table$preds %in% c("S_10m","S_zeu","S_mld"),"net"]), nrow = 3, scales = "free")

### C) o2_10m vs o2_zeu vs o2_mld
quartz()
ggplot(aes(x = factor(preds), y = R2), data = table[table$preds %in% c("o2_10m","o2_zeu","o2_mld"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("R2") +
    facet_wrap(factor(table[table$preds %in% c("o2_10m","o2_zeu","o2_mld"),"net"]), nrow = 3, scales = "free")
#
quartz()
ggplot(aes(x = factor(preds), y = AIC), data = table[table$preds %in% c("o2_10m","o2_zeu","o2_mld"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("AIC") +
    facet_wrap(factor(table[table$preds %in% c("o2_10m","o2_zeu","o2_mld"),"net"]), nrow = 3, scales = "free")
    
### D) Chla_10m vs Chla_zeu vs Chla_mld
quartz()
ggplot(aes(x = factor(preds), y = R2), data = table[table$preds %in% c("Chla_10m","Chla_zeu","Chla_mld"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("R2") +
    facet_wrap(factor(table[table$preds %in% c("Chla_10m","Chla_zeu","Chla_mld"),"net"]), nrow = 3, scales = "free")
#
quartz()
ggplot(aes(x = factor(preds), y = AIC), data = table[table$preds %in% c("Chla_10m","Chla_zeu","Chla_mld"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("AIC") +
    facet_wrap(factor(table[table$preds %in% c("Chla_10m","Chla_zeu","Chla_mld"),"net"]), nrow = 3, scales = "free")


### 08/01/2020: also check the range covered by the different versions and the nb of missing values
summary(imagery2[imagery2$net == "regent",c("T_10m","T_zeu","T_mld")])
summary(imagery2[imagery2$net == "regent",c("S_10m","S_zeu","S_mld")])
summary(imagery2[imagery2$net == "regent",c("o2_10m","o2_zeu","o2_mld")])
summary(imagery2[imagery2$net == "regent",c("Chla_10m","Chla_zeu","Chla_mld")])
### --> measurements @ the basis of the euphotic layer present way more missing values. Meanwhile the 2 remaining versions present very similar ranges
### and the models present simialr distrbutions in expalined deviance/ R2 etc. 

### CONCLUSIONS --> use the 10m depth only

### E) Si vs NO3 vs PO4
summary(imagery2[imagery2$net == "wp2",c("Si","NO3","PO4")])
summary(imagery2[imagery2$net == "bongo",c("Si","NO3","PO4")])
summary(imagery2[imagery2$net == "regent",c("Si","NO3","PO4")])

quartz()
ggplot(aes(x = factor(preds), y = R2), data = table[table$preds %in% c("Si","NO3","PO4"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("R2") +
    facet_wrap(factor(table[table$preds %in% c("Si","NO3","PO4"),"net"]), nrow = 3, scales = "free")
# Unclear
quartz()
ggplot(aes(x = factor(preds), y = AIC), data = table[table$preds %in% c("Si","NO3","PO4"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("AIC") +
    facet_wrap(factor(table[table$preds %in% c("Si","NO3","PO4"),"net"]), nrow = 3, scales = "free")
# Still unclear, look at F scores or p-values
quartz()
ggplot(aes(x = factor(preds), y = log(get('F'))), data = table[table$preds %in% c("Si","NO3","PO4"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = F) + theme_classic() + xlab("") + ylab("F scores") +
    facet_wrap(factor(table[table$preds %in% c("Si","NO3","PO4"),"net"]), nrow = 3, scales = "free")
# Still unclear...
quartz()
ggplot(aes(x = factor(preds), y = log(pval)), data = table[table$preds %in% c("Si","NO3","PO4"),]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = F) + theme_classic() + xlab("") + ylab("P") +
    facet_wrap(factor(table[table$preds %in% c("Si","NO3","PO4"),"net"]), nrow = 3, scales = "free")

### --> all 3 macronutrients present similar nb of missing values (~10 to 15 depending on the net)
### All 3 present simialr range of units (~9/10 units)
### Impossible to choose yet


### F) o2 vs Temp
quartz()
ggplot(aes(x = factor(preds), y = R2),
    data = table[table$preds %in% c("T_10m","o2_10m") & table$net == "bongo",]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("R2") +
    facet_wrap(factor(table[table$preds %in% c("T_10m","o2_10m") & table$net == "bongo","net"]), nrow = 3, scales = "free")
#
quartz()
ggplot(aes(x = factor(preds), y = AIC),
    data = table[table$preds %in% c("T_10m","o2_10m") & table$net == "bongo",]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("AIC") +
    facet_wrap(factor(table[table$preds %in% c("T_10m","o2_10m") & table$net == "bongo","net"]), nrow = 3, scales = "free")
#
quartz()
ggplot(aes(x = factor(preds), y = log(get('F')) ),
    data = table[table$preds %in% c("T_10m","o2_10m") & table$net == "bongo",]) + 
    geom_boxplot(aes(fill = factor(preds)), colour = "black", notch = T) + theme_classic() + xlab("") + ylab("F scores (log)") +
    facet_wrap(factor(table[table$preds %in% c("T_10m","o2_10m") & table$net == "bongo","net"]), nrow = 3, scales = "free")

### For WP2, o2 actually seems to lead to slightly better models. For the regent net it's impossible to tell yet. Same for the Bongo. 

### CONCLUSIONS: now, for each net and every response variable, find the best model (highest R2/ lowest AIC) among the ones that contain:
# - 10m depth measurements
# - T_10m   |
# - o2_10m  |
# - all 3 macronutrients to help choose 
require("parallel")
resp <- unique(table$resp)
# r <- "Abund_Calanoida"
# n <- "wp2"
res <- mclapply(resp, function(r) {
    
            # Subset based on r
            subset <- table[table$resp == r,]
            
            res2 <- lapply(nets, function(n) {
                
                        # Useless message
                        message(paste("Finding best GAM model for ",r," sampled with ",n, sep = ""))
                        subset2 <- subset[subset$net == n,]
                        # Pick the best model with T_10m
                        max.r2.1 <- max(subset2[which(subset2$preds == "T_10m"),"R2"])
                        # and the one with o2_10m
                        max.r2.2 <- max(subset2[which(subset2$preds == "o2_10m"),"R2"])
                        # Subset base don that
                        best.temp <- subset2[subset2$R2 == max.r2.1,]
                        best.o2 <- subset2[subset2$R2 == max.r2.2,]
                        # Rbind and return
                        best <- rbind(best.temp, best.o2)
                        return(best)
                        
                    } # eo FUN
                    
            ) # eo lapply - nets 
            # Rbind
            table.best <- do.call(rbind, res2)
            rm(res2); gc()
            return(table.best)
    
    }, mc.cores = 3
    
) # eo mclapply - resp
# Rbind
best.models <- do.call(rbind, res)
head(best.models) ; dim(best.models)

### Examine the best models for each resp*net
best.models[best.models$net == "wp2",]  # dim(best.models[best.models$net == "wp2",])
### Would be better to switch to a wider data.frame format -> dcast
best <- dcast(data = best.models[,c(1:6,8)], formula =  net+resp+formula+R2+AIC ~ preds, fun.aggregate = mean)
head(best) ; dim(best)
best[1:10,]
summary(best)
# Ok, get rid of the ones with _mld and _zeu
models2rm <- c(unique(best$formula)[grepl("_zeu",unique(best$formula))] , unique(best$formula)[grepl("_mld",unique(best$formula))])
models2rm
# Remove them from 'best'
best2 <- best[!(best$formula %in% models2rm),]
dim(best2)
unique(best2$resp) # Ok, no pred got removed
summary(best2)

# Drop the _zeu and _mld columns
drop.cols <- c("T_mld","S_zeu","S_mld","o2_zeu","o2_mld","Chla_zeu","Chla_mld")
best2 <- best2 %>% select (-drop.cols)
summary(best2)

### In a for loop, display the outputs for each resp
# r <- "Abund_Candaciidae"
best2[best2$resp == r,]

for(r in unique(best2$resp)) {
        message(paste("Displaying best models for ",r," for Régent net -------------------------------------", sep = ""))
        sub <- best2[best2$resp == r & best2$net == "regent",]
        message("   ")
        # And then display the ranking of vars n table
        print( sub )
        if(nrow(sub) > 0) {
                for(i in c(1:nrow(sub))) {
                    form <- sub[i,"formula"]
                    message(paste("Displaying predictors ranking for ",r," modelled with ",form, sep = ""))
                    sub.table <- table[which(table$formula == form & table$resp == r & table$net == "regent"),]
                    print( sub.table[order(sub.table$F, decreasing = T),c("preds","pval","F")] )
                } # eo for loop
        } # eo if loop 
        message("   ")
        message("   ")
} # eo for loop - r in resp


### Examine distrbution of ESD and Abund models' R2 per net
resp1 <- unique(best$resp)[grepl("ESD",unique(best$resp))] # ESD response vars
resp2 <- unique(best$resp)[grepl("Abund",unique(best$resp))] # Abund resp vars
# ESD first

quartz()
ggplot(aes(x = factor(net), y = R2), data = best[best$resp %in% resp1,]) + 
    geom_violin(aes(fill = factor(net)), colour = "black",) + 
    geom_boxplot(fill = "white ", colour = "black", width = 0.1) + 
    theme_classic() + xlab("Plankton net type") + ylab("R2 of ESD-based GAMs") + 
    scale_fill_brewer(name = "", palette = "Paired")
#
quartz()
ggplot(aes(x = factor(net), y = R2), data = best[best$resp %in% resp2,]) + 
    geom_violin(aes(fill = factor(net)), colour = "black") + 
    geom_boxplot(fill = "white ", colour = "black", width = 0.1) + 
    theme_classic() + xlab("Plankton net type") + ylab("R2 of Abundance-based GAMs") +
    scale_fill_brewer(name = "", palette = "Paired")
#
quartz()
ggplot(aes(x = factor(net), y = AIC), data = best[best$resp %in% resp1,]) + 
    geom_violin(aes(fill = factor(net)), colour = "black",) + 
    geom_boxplot(fill = "white ", colour = "black", width = 0.1) + 
    theme_classic() + xlab("Plankton net type") + ylab("AIC of ESD-based GAMs") + 
    scale_fill_brewer(name = "", palette = "Paired")
#
quartz()
ggplot(aes(x = factor(net), y = AIC), data = best[best$resp %in% resp2,]) + 
    geom_violin(aes(fill = factor(net)), colour = "black") + 
    geom_boxplot(fill = "white ", colour = "black", width = 0.1) + 
    theme_classic() + xlab("Plankton net type") + ylab("AIC of Abundance-based GAMs") +
    scale_fill_brewer(name = "", palette = "Paired")
    
### --> Bongo's the way to go

### Save datasets so you don't need to re-run the script everytime
setwd("/Users/fabiobenedetti/Desktop/~work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/")
write.table(table, file = "table_GAMs_stats_all.txt", sep = "\t")
write.table(best.models, file = "table_GAMs_stats_best.txt", sep = "\t")
write.table(imagery2, file = "TARA_OCEANS_table_imagery2_transformed.txt", sep = "\t")

### ---------------------------------------------------------------------------------------------------------------------------

### 10/01/2020: Now that you've narrowed down the main net type and those response variables that were adequately model by the GAMs, 
### use the models and the training data to predict the responses and draw the response curves, for the interesting predictors

library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("vegan")
library("FactoMineR")
library("matrixStats")
library("lubridate")
library("viridis")
library("maps")
library("marmap")
library("mgcv")
library("nlme")
library("MuMIn")
library("ggrepel")

world <- map_data("world") # coastline for maps
WD <- getwd()

table <- read.table("table_GAMs_stats_all.txt", sep = "\t", h = T) 
best <- read.table("table_GAMs_stats_best.txt", sep = "\t", h = T) 
# dim(table) ; dim(best) # 23% of all models
imagery2 <- read.table("TARA_OCEANS_table_imagery2_transformed.txt", sep = "\t", h = T) 

### Define the response variables to draw the resp curves for the Bongo dataset
esd.resp <- c("ESD_Augaptilidae","ESD_Calanidae","ESD_Calanoida","ESD_Candaciidae","ESD_Chaetognatha","ESD_Cnidaria",
            "ESD_Copepoda","ESD_Eucalanidae","ESD_Euchaetidae","ESD_Oithonidae","ESD_Oncaeidae","ESD_Protista",
            "ESD_Sapphirinidae","ESD_Sapphirinidae","ESD_Zooplankton")
# same for abundances
abund.resp <- c("Abund_Acartiidae","Abund_Aetideidae","Abund_Augaptilidae","Abund_Calanoida","Abund_Candaciidae",
                "Abund_Centropagidae","Abund_Chaetognatha","Abund_Copepoda","Abund_Cyclopoida","Abund_Decapoda",
                "Abund_Euchaetidae","Abund_Harpacticoida","Abund_Oithonidae","Abund_otherCopepoda","Abund_Protista",
                "Abund_Rhincalanidae","Abund_Sapphirinidae","Abund_Temoridae","Abund_Urochordata","Abund_Zooplankton")

### Develop a code for plotting resp curves of signif terms for one test response var
n <- "bongo"
r <- "ESD_Zooplankton"
best[best$net == n & best$resp == r,]
# Associated unique formulae
forms <- unique(best[best$net == n & best$resp == r,"formula"])
# Get the corresponding model objects
model.names <- paste("gam_",r,"_",forms,"_",n,".Rdata", sep = "")
m <- model.names[3]

for(m in model.names) {
        
        # get model object
        message(paste("Drawing response curves for ", r," based on model = ", m, sep=""))
        setwd(paste(WD,"/","Models/GAMs/",n, sep = ""))
        model <- get(load(m))
        
        # Extract the terms of the model from the filename
        terms <- do.call(cbind, strsplit(as.character(m),"\\+"))
        # 1st element contains response var and first predictor
        terms[1,] <- str_replace_all(terms[1,], "gam_", "")
        
        # Watchout, if other_Copepoda is the group, then there's an extra underscore to consider
        if( grepl("other_Copepoda", terms[1,]) ) {
            resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2],unlist(strsplit(terms[1,],"_"))[3],sep = "_")
            pred1 <- paste(unlist(strsplit(terms[1,],"_"))[4], unlist(strsplit(terms[1,],"_"))[5], sep = "_")
            predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
        } else {
            resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2], sep = "_")
            pred1 <- paste(unlist(strsplit(terms[1,],"_"))[3], unlist(strsplit(terms[1,],"_"))[4], sep = "_")
            predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
        } # eo if else loop - grepl("other_Copepoda", terms[1,]) 
            
        if(n == "wp2") {
            preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],predlast)
        } else {  
            preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],terms[9,],predlast)    
        } # eo if else loop
        
        # Subset 'imagery2' to draw resp curves
        subset <- imagery2[imagery2$net == n,]
        subset <- subset[,c(resp,preds)]
        subset <- na.omit(subset)
        colnames(subset) <- c("y","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10")
        k <- round(nrow(subset)/10, digits = 0)
        
        # Train the modle again but by changing the names of the terms to draw the resp curves
        model.new <- mgcv::gam(formula = y ~ s(x1,bs="tp",k=k) + s(x2,bs="tp",k=k) + s(x3,bs="tp",k=k) + s(x4,bs="tp",k=k) +
                        s(x5,bs="tp",k=k) + s(x6,bs="tp",k=k) + s(x7,bs="tp",k=k) + s(x8,bs="tp",k=k) + s(x9,bs="tp",k=k) +
                        s(x10,bs="tp",k=k), data = subset, select = TRUE, method = "REML")
        
        ## Predict r from model
        # mgcv::predict.gam(object = model.new, newdata = subset[,c(2:11)], type = "response", se.fit = T)
        # Now, for each of the 10 predictors, generate a new data.frame allowing to examine the resp curve of all the signif. terms
        P <- data.frame(term = c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"), pred = preds, pval = summary(model.new)$s.pv)
        
        # Terms to plot 
        signif <- P[P$pval < 0.05,c("pred","term")]
        
        # For each of the 10 predictors, create a new data frame meant to help you draw response curves
        new.data.1 <- data.frame(x1 = seq(from=min(subset$x1),to=max(subset$x1),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
            x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
            x9 = mean(subset$x9), x10 = mean(subset$x10), pred2draw = P[1,"pred"]) # eo ddf
        # Predict
        new.data.1$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.1, type = "response", se.fit = T)$fit
        new.data.1$se <- mgcv::predict.gam(object = model.new, newdata = new.data.1, type = "response", se.fit = T)$se
        
        new.data.2 <- data.frame(x2 = seq(from=min(subset$x2),to=max(subset$x2),length=100), x1 = mean(subset$x1), x3 = mean(subset$x3),
            x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
            x9 = mean(subset$x9), x10 = mean(subset$x10), pred2draw = P[2,"pred"]) # eo ddf
        # Predict
        new.data.2$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.2, type = "response", se.fit = T)$fit
        new.data.2$se <- mgcv::predict.gam(object = model.new, newdata = new.data.2, type = "response", se.fit = T)$se
      
        
        new.data.3 <- data.frame(x3 = seq(from=min(subset$x3),to=max(subset$x3),length=100), x2 = mean(subset$x2), x1 = mean(subset$x1),
            x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
            x9 = mean(subset$x9), x10 = mean(subset$x10), pred2draw = P[3,"pred"]) # eo ddf
        # Predict
        new.data.3$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.3, type = "response", se.fit = T)$fit
        new.data.3$se <- mgcv::predict.gam(object = model.new, newdata = new.data.3, type = "response", se.fit = T)$se
      
        
        new.data.4 <- data.frame(x4 = seq(from=min(subset$x4),to=max(subset$x4),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
            x1 = mean(subset$x1), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
            x9 = mean(subset$x9), x10 = mean(subset$x10), pred2draw = P[4,"pred"]) # eo ddf
        # Predict
        new.data.4$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.4, type = "response", se.fit = T)$fit
        new.data.4$se <- mgcv::predict.gam(object = model.new, newdata = new.data.4, type = "response", se.fit = T)$se
      
        
        new.data.5 <- data.frame(x5 = seq(from=min(subset$x5),to=max(subset$x5),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
            x4 = mean(subset$x4), x1 = mean(subset$x1), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
            x9 = mean(subset$x9), x10 = mean(subset$x10), pred2draw = P[5,"pred"]) # eo ddf
        # Predict
        new.data.5$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.5, type = "response", se.fit = T)$fit
        new.data.5$se <- mgcv::predict.gam(object = model.new, newdata = new.data.5, type = "response", se.fit = T)$se
      
        
        new.data.6 <- data.frame(x6 = seq(from=min(subset$x6),to=max(subset$x6),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
            x4 = mean(subset$x4), x5 = mean(subset$x5), x1 = mean(subset$x1), x7 = mean(subset$x7), x8 = mean(subset$x8),
            x9 = mean(subset$x9), x10 = mean(subset$x10), pred2draw = P[6,"pred"]) # eo ddf
        # Predict
        new.data.6$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.6, type = "response", se.fit = T)$fit
        new.data.6$se <- mgcv::predict.gam(object = model.new, newdata = new.data.6, type = "response", se.fit = T)$se
      
        
        new.data.7 <- data.frame(x7 = seq(from=min(subset$x7),to=max(subset$x7),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
            x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x1 = mean(subset$x1), x8 = mean(subset$x8),
            x9 = mean(subset$x9), x10 = mean(subset$x10), pred2draw = P[7,"pred"]) # eo ddf
        # Predict
        new.data.7$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.7, type = "response", se.fit = T)$fit
        new.data.7$se <- mgcv::predict.gam(object = model.new, newdata = new.data.7, type = "response", se.fit = T)$se
      
        
        new.data.8 <- data.frame(x8 = seq(from=min(subset$x8),to=max(subset$x8),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
            x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x1 = mean(subset$x1),
            x9 = mean(subset$x9), x10 = mean(subset$x10), pred2draw = P[8,"pred"]) # eo ddf
        # Predict
        new.data.8$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.8, type = "response", se.fit = T)$fit
        new.data.8$se <- mgcv::predict.gam(object = model.new, newdata = new.data.8, type = "response", se.fit = T)$se
      
        
        new.data.9 <- data.frame(x9 = seq(from=min(subset$x9),to=max(subset$x9),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
            x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
            x1 = mean(subset$x1), x10 = mean(subset$x10), pred2draw = P[9,"pred"]) # eo ddf
        # Predict
        new.data.9$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.9, type = "response", se.fit = T)$fit
        new.data.9$se <- mgcv::predict.gam(object = model.new, newdata = new.data.9, type = "response", se.fit = T)$se
      
        new.data.10 <- data.frame(x10 = seq(from=min(subset$x10),to=max(subset$x10),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
            x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
            x9 = mean(subset$x9), x1 = mean(subset$x1), pred2draw = P[10,"pred"]) # eo ddf
        # Predict
        new.data.10$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.10, type = "response", se.fit = T)$fit
        new.data.10$se <- mgcv::predict.gam(object = model.new, newdata = new.data.10, type = "response", se.fit = T)$se
      
        ### Make the 10 plots
        setwd(paste(WD,"/","Models/GAMs/",n,"/rcurves2/", sep = ""))
        dir.create(m)
        setwd(paste(WD,"/","Models/GAMs/",n,"/rcurves2/",m, sep = ""))
        
        p1 <- ggplot() + geom_point(aes(x = x1, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x1, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.1) +
                geom_line(aes(x = x1, y = fit), colour = "black", data = new.data.1) +
                xlab(paste(P[P$term=="x1","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p1, filename = paste("plot_rcurve_",r,"_",P[P$term=="x1","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        p2 <- ggplot() + geom_point(aes(x = x2, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x2, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.2) +
                geom_line(aes(x = x2, y = fit), colour = "black", data = new.data.2) +
                xlab(paste(P[P$term=="x2","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p2, filename = paste("plot_rcurve_",r,"_",P[P$term=="x2","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        p3 <- ggplot() + geom_point(aes(x = x3, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x3, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.3) +
                geom_line(aes(x = x3, y = fit), colour = "black", data = new.data.3) +
                xlab(paste(P[P$term=="x3","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p3, filename = paste("plot_rcurve_",r,"_",P[P$term=="x3","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        p4 <- ggplot() + geom_point(aes(x = x4, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x4, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.4) +
                geom_line(aes(x = x4, y = fit), colour = "black", data = new.data.4) +
                xlab(paste(P[P$term=="x4","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p4, filename = paste("plot_rcurve_",r,"_",P[P$term=="x4","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        p5 <- ggplot() + geom_point(aes(x = x5, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x5, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.5) +
                geom_line(aes(x = x5, y = fit), colour = "black", data = new.data.5) +
                xlab(paste(P[P$term=="x5","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p5, filename = paste("plot_rcurve_",r,"_",P[P$term=="x5","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        p6 <- ggplot() + geom_point(aes(x = x6, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x6, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.6) +
                geom_line(aes(x = x6, y = fit), colour = "black", data = new.data.6) +
                xlab(paste(P[P$term=="x6","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p6, filename = paste("plot_rcurve_",r,"_",P[P$term=="x6","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        p7 <- ggplot() + geom_point(aes(x = x7, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x7, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.7) +
                geom_line(aes(x = x7, y = fit), colour = "black", data = new.data.7) +
                xlab(paste(P[P$term=="x7","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p7, filename = paste("plot_rcurve_",r,"_",P[P$term=="x7","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        p8 <- ggplot() + geom_point(aes(x = x8, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x8, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.8) +
                geom_line(aes(x = x8, y = fit), colour = "black", data = new.data.8) +
                xlab(paste(P[P$term=="x8","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p8, filename = paste("plot_rcurve_",r,"_",P[P$term=="x8","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        p9 <- ggplot() + geom_point(aes(x = x9, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x9, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.9) +
                geom_line(aes(x = x9, y = fit), colour = "black", data = new.data.9) +
                xlab(paste(P[P$term=="x9","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p9, filename = paste("plot_rcurve_",r,"_",P[P$term=="x9","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        p10 <- ggplot() + geom_point(aes(x = x10, y = y), data = subset, colour = "black", alpha = 0.2) + 
                geom_ribbon(aes(x = x10, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.10) +
                geom_line(aes(x = x10, y = fit), colour = "black", data = new.data.10) +
                xlab(paste(P[P$term=="x10","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
        ggsave(plot = p10, filename = paste("plot_rcurve_",r,"_",P[P$term=="x10","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
        ### Clean stuff etc.
        rm(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,new.data.1,new.data.2,new.data.3,new.data.4,new.data.5,new.data.6,new.data.7,new.data.8,
            new.data.9,new.data.10,P)
        setwd(paste(WD,"/","Models/GAMs/",n, sep = ""))
        
}

### Now use this for loop for every model
n <- "regent"
# r <- "ESD_Protista" # for testing
for(r in abund.resp) {
    
    message(paste("", sep = ""))
    message(paste("Drawing response curves for ",r, sep = ""))
    message(paste("", sep = ""))
    message(paste("", sep = ""))
    
    # best[best$net == n & best$resp == r,]
    # Associated unique formulae
    forms <- unique(best[best$net == n & best$resp == r,"formula"])
    
    if( length(forms) > 0 ) {
    
        # Get the corresponding model objects
        model.names <- paste("gam_",r,"_",forms,"_",n,".Rdata", sep = "")
        # m <- model.names[1]
        for(m in model.names) {
        
            # get model object
            message(paste("Drawing response curves for ", r," based on model = ", m, sep=""))
            setwd(paste(WD,"/","Models/GAMs/",n, sep = ""))
            model <- get(load(m))   # summary(model)
            # quartz() ; plot(model)
            # Extract the terms of the model from the filename
            terms <- do.call(cbind, strsplit(as.character(m),"\\+"))
            # 1st element contains response var and first predictor
            terms[1,] <- str_replace_all(terms[1,], "gam_", "")
        
            # Watchout, if other_Copepoda is the group, then there's an extra underscore to consider
            if( grepl("other_Copepoda", terms[1,]) ) {
                resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2],unlist(strsplit(terms[1,],"_"))[3],sep = "_")
                pred1 <- paste(unlist(strsplit(terms[1,],"_"))[4], unlist(strsplit(terms[1,],"_"))[5], sep = "_")
                predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
            } else {
                resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2], sep = "_")
                pred1 <- paste(unlist(strsplit(terms[1,],"_"))[3], unlist(strsplit(terms[1,],"_"))[4], sep = "_")
                predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
            } # eo if else loop - grepl("other_Copepoda", terms[1,]) 
            
            if(n == "wp2") {
                preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],terms[9,],predlast)
            } else {  
                preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],terms[9,],terms[10,],predlast)    
            } # eo if else loop
        
            # Subset 'imagery2' to draw resp curves
            subset <- imagery2[imagery2$net == n,]
            subset <- subset[,c(resp,preds)]
            subset <- na.omit(subset)
            colnames(subset) <- c("y","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11")
            k <- round(nrow(subset)/11, digits = 0)
        
            # Train the model again but by changing the names of the terms to draw the resp curves
            model.new <- mgcv::gam(formula = y ~ s(x1,bs="tp",k=k) + s(x2,bs="tp",k=k) + s(x3,bs="tp",k=k) + s(x4,bs="tp",k=k) +
                            s(x5,bs="tp",k=k) + s(x6,bs="tp",k=k) + s(x7,bs="tp",k=k) + s(x8,bs="tp",k=k) + s(x9,bs="tp",k=k) +
                            s(x10,bs="tp",k=k) + s(x11,bs="tp",k=k) 
                            , data = subset, select = T, method = "REML")
        
            ## Predict r from model
            # mgcv::predict.gam(object = model.new, newdata = subset[,c(2:11)], type = "response", se.fit = T)
            # Now, for each of the 10 predictors, generate a new data.frame allowing to examine the resp curve of all the signif. terms
            P <- data.frame(term = c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11")
                , pred = preds, pval = summary(model.new)$s.pv)
        
            # Terms to plot 
            signif <- P[P$pval < 0.05,c("pred","term")]
        
            # For each of the 10 predictors, create a new data frame meant to help you draw response curves
            new.data.1 <- data.frame(x1 = seq(from=min(subset$x1),to=max(subset$x1),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
                x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
                x9 = mean(subset$x9), x10 = mean(subset$x10), x11 = mean(subset$x11),
                pred2draw = P[1,"pred"]) # eo ddf
            # Predict
            new.data.1$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.1, type = "response", se.fit = T)$fit
            new.data.1$se <- mgcv::predict.gam(object = model.new, newdata = new.data.1, type = "response", se.fit = T)$se
        
            new.data.2 <- data.frame(x2 = seq(from=min(subset$x2),to=max(subset$x2),length=100), x1 = mean(subset$x1), x3 = mean(subset$x3),
                x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
                x9 = mean(subset$x9), x10 = mean(subset$x10), x11 = mean(subset$x11),
                pred2draw = P[2,"pred"]) # eo ddf
            # Predict
            new.data.2$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.2, type = "response", se.fit = T)$fit
            new.data.2$se <- mgcv::predict.gam(object = model.new, newdata = new.data.2, type = "response", se.fit = T)$se
      
        
            new.data.3 <- data.frame(x3 = seq(from=min(subset$x3),to=max(subset$x3),length=100), x2 = mean(subset$x2), x1 = mean(subset$x1),
                x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
                x9 = mean(subset$x9), x10 = mean(subset$x10), x11 = mean(subset$x11),
                pred2draw = P[3,"pred"]) # eo ddf
            # Predict
            new.data.3$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.3, type = "response", se.fit = T)$fit
            new.data.3$se <- mgcv::predict.gam(object = model.new, newdata = new.data.3, type = "response", se.fit = T)$se
      
        
            new.data.4 <- data.frame(x4 = seq(from=min(subset$x4),to=max(subset$x4),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
                x1 = mean(subset$x1), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
                x9 = mean(subset$x9), x10 = mean(subset$x10), x11 = mean(subset$x11),
                pred2draw = P[4,"pred"]) # eo ddf
            # Predict
            new.data.4$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.4, type = "response", se.fit = T)$fit
            new.data.4$se <- mgcv::predict.gam(object = model.new, newdata = new.data.4, type = "response", se.fit = T)$se
      
        
            new.data.5 <- data.frame(x5 = seq(from=min(subset$x5),to=max(subset$x5),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
                x4 = mean(subset$x4), x1 = mean(subset$x1), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
                x9 = mean(subset$x9), x10 = mean(subset$x10), x11 = mean(subset$x11),
                pred2draw = P[5,"pred"]) # eo ddf
            # Predict
            new.data.5$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.5, type = "response", se.fit = T)$fit
            new.data.5$se <- mgcv::predict.gam(object = model.new, newdata = new.data.5, type = "response", se.fit = T)$se
      
        
            new.data.6 <- data.frame(x6 = seq(from=min(subset$x6),to=max(subset$x6),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
                x4 = mean(subset$x4), x5 = mean(subset$x5), x1 = mean(subset$x1), x7 = mean(subset$x7), x8 = mean(subset$x8),
                x9 = mean(subset$x9), x10 = mean(subset$x10), x11 = mean(subset$x11),
                pred2draw = P[6,"pred"]) # eo ddf
            # Predict
            new.data.6$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.6, type = "response", se.fit = T)$fit
            new.data.6$se <- mgcv::predict.gam(object = model.new, newdata = new.data.6, type = "response", se.fit = T)$se
      
        
            new.data.7 <- data.frame(x7 = seq(from=min(subset$x7),to=max(subset$x7),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
                x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x1 = mean(subset$x1), x8 = mean(subset$x8),
                x9 = mean(subset$x9), x10 = mean(subset$x10), x11 = mean(subset$x11),
                pred2draw = P[7,"pred"]) # eo ddf
            # Predict
            new.data.7$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.7, type = "response", se.fit = T)$fit
            new.data.7$se <- mgcv::predict.gam(object = model.new, newdata = new.data.7, type = "response", se.fit = T)$se
      
        
            new.data.8 <- data.frame(x8 = seq(from=min(subset$x8),to=max(subset$x8),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
                x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x1 = mean(subset$x1),
                x9 = mean(subset$x9), x10 = mean(subset$x10), x11 = mean(subset$x11),
                pred2draw = P[8,"pred"]) # eo ddf
            # Predict
            new.data.8$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.8, type = "response", se.fit = T)$fit
            new.data.8$se <- mgcv::predict.gam(object = model.new, newdata = new.data.8, type = "response", se.fit = T)$se
      
        
            new.data.9 <- data.frame(x9 = seq(from=min(subset$x9),to=max(subset$x9),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
                x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
                x1 = mean(subset$x1), x10 = mean(subset$x10), x11 = mean(subset$x11),
                pred2draw = P[9,"pred"]) # eo ddf
            # Predict
            new.data.9$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.9, type = "response", se.fit = T)$fit
            new.data.9$se <- mgcv::predict.gam(object = model.new, newdata = new.data.9, type = "response", se.fit = T)$se
      
            new.data.10 <- data.frame(x10 = seq(from=min(subset$x10),to=max(subset$x10),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
                x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
                x9 = mean(subset$x9), x1 = mean(subset$x1), x11 = mean(subset$x11),
                pred2draw = P[10,"pred"]) # eo ddf
            # Predict
            new.data.10$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.10, type = "response", se.fit = T)$fit
            new.data.10$se <- mgcv::predict.gam(object = model.new, newdata = new.data.10, type = "response", se.fit = T)$se
      
             new.data.11 <- data.frame(x11 = seq(from=min(subset$x11),to=max(subset$x11),length=100), x2 = mean(subset$x2), x3 = mean(subset$x3),
                 x4 = mean(subset$x4), x5 = mean(subset$x5), x6 = mean(subset$x6), x7 = mean(subset$x7), x8 = mean(subset$x8),
                 x9 = mean(subset$x9), x1 = mean(subset$x1), x10 = mean(subset$x10),
                 pred2draw = P[11,"pred"]) # eo ddf
             # Predict
             new.data.11$fit <- mgcv::predict.gam(object = model.new, newdata = new.data.11, type = "response", se.fit = T)$fit
             new.data.11$se <- mgcv::predict.gam(object = model.new, newdata = new.data.11, type = "response", se.fit = T)$se

            ### Make the 10 plots
            setwd(paste(WD,"/","Models/GAMs/",n,"/rcurves2/", sep = ""))
            dir.create(m)
            setwd(paste(WD,"/","Models/GAMs/",n,"/rcurves2/",m, sep = ""))
        
            p1 <- ggplot() + #geom_point(aes(x = x1, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x1, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.1) +
                    geom_line(aes(x = x1, y = fit), colour = "black", data = new.data.1) +
                    xlab(paste(P[P$term=="x1","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p1, filename = paste("plot_rcurve_",r,"_",P[P$term=="x1","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p2 <- ggplot() + #geom_point(aes(x = x2, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x2, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.2) +
                    geom_line(aes(x = x2, y = fit), colour = "black", data = new.data.2) +
                    xlab(paste(P[P$term=="x2","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p2, filename = paste("plot_rcurve_",r,"_",P[P$term=="x2","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p3 <- ggplot() + #geom_point(aes(x = x3, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x3, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.3) +
                    geom_line(aes(x = x3, y = fit), colour = "black", data = new.data.3) +
                    xlab(paste(P[P$term=="x3","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p3, filename = paste("plot_rcurve_",r,"_",P[P$term=="x3","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p4 <- ggplot() + #geom_point(aes(x = x4, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x4, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.4) +
                    geom_line(aes(x = x4, y = fit), colour = "black", data = new.data.4) +
                    xlab(paste(P[P$term=="x4","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p4, filename = paste("plot_rcurve_",r,"_",P[P$term=="x4","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p5 <- ggplot() + #geom_point(aes(x = x5, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x5, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.5) +
                    geom_line(aes(x = x5, y = fit), colour = "black", data = new.data.5) +
                    xlab(paste(P[P$term=="x5","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p5, filename = paste("plot_rcurve_",r,"_",P[P$term=="x5","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p6 <- ggplot() + #geom_point(aes(x = x6, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x6, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.6) +
                    geom_line(aes(x = x6, y = fit), colour = "black", data = new.data.6) +
                    xlab(paste(P[P$term=="x6","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p6, filename = paste("plot_rcurve_",r,"_",P[P$term=="x6","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p7 <- ggplot() + #geom_point(aes(x = x7, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x7, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.7) +
                    geom_line(aes(x = x7, y = fit), colour = "black", data = new.data.7) +
                    xlab(paste(P[P$term=="x7","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p7, filename = paste("plot_rcurve_",r,"_",P[P$term=="x7","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p8 <- ggplot() + #geom_point(aes(x = x8, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x8, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.8) +
                    geom_line(aes(x = x8, y = fit), colour = "black", data = new.data.8) +
                    xlab(paste(P[P$term=="x8","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p8, filename = paste("plot_rcurve_",r,"_",P[P$term=="x8","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p9 <- ggplot() + #geom_point(aes(x = x9, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x9, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.9) +
                    geom_line(aes(x = x9, y = fit), colour = "black", data = new.data.9) +
                    xlab(paste(P[P$term=="x9","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p9, filename = paste("plot_rcurve_",r,"_",P[P$term=="x9","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p10 <- ggplot() + #geom_point(aes(x = x10, y = y), data = subset, colour = "black", alpha = 0.2) + 
                    geom_ribbon(aes(x = x10, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.10) +
                    geom_line(aes(x = x10, y = fit), colour = "black", data = new.data.10) +
                    xlab(paste(P[P$term=="x10","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
                
            ggsave(plot = p10, filename = paste("plot_rcurve_",r,"_",P[P$term=="x10","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
        
            p11 <- ggplot() + #geom_point(aes(x = x10, y = y), data = subset, colour = "black", alpha = 0.2) +
                     geom_ribbon(aes(x = x11, ymin = fit-se, ymax = fit+se), alpha = 0.2, fill = "black", data = new.data.11) +
                     geom_line(aes(x = x11, y = fit), colour = "black", data = new.data.11) +
                     xlab(paste(P[P$term=="x11","pred"],sep="")) + ylab(paste(r,sep="")) + theme_classic()
 
            ggsave(plot = p11, filename = paste("plot_rcurve_",r,"_",P[P$term=="x11","pred"],".pdf", sep = ""), dpi = 300, width = 6, height = 5)
    
            ### Clean stuff etc.
            rm(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
                new.data.1,new.data.2,new.data.3,
                new.data.4,new.data.5,new.data.6,
                new.data.7,new.data.8,new.data.9,
                new.data.10,P)
            setwd(paste(WD,"/","Models/GAMs/",n, sep = ""))       
            
            } # eo for loop
        
         } else {
             
             message(paste("Skipping ",r, sep = ""))
        
    }
    
} # eo for loop - r in abund.resp


### ---------------------------------------------------------------------------------------------------------------------------

### 19/02/2020: After talking with Fabien Lombard @ OSM 2020, decided on some stuff:
# - stick mainly to ESD (new stuff)
# - avoid focusing on copepod families or subgroups finer than classes because they have much lower number of objects (unrealiable due to unsufficient data density)
# - use overall ranks of predictors within the GAMs to identify the top 6/7/8 variables 
# - use those in a signle model and extract the r curves from this model to show on a panel of figure (Groups x predictors)

library("tidyverse")
library("reshape2")
library("RColorBrewer")
library("vegan")
library("FactoMineR")
library("matrixStats")
library("lubridate")
library("viridis")
library("maps")
library("marmap")
library("mgcv")
library("nlme")
library("MuMIn")
library("ggrepel")

world <- map_data("world") # coastline for maps
WD <- getwd()

table <- read.table("table_GAMs_stats_all.txt", sep = "\t", h = T) 
best <- read.table("table_GAMs_stats_best.txt", sep = "\t", h = T) 
# dim(table) ; dim(best) # 23% of all models
imagery2 <- read.table("TARA_OCEANS_table_imagery2_transformed.txt", sep = "\t", h = T) 

### 11/03/2020: For these responses variables, check if you have any good models and rank the variables so you can identify the reponse curves to put in new panel figure 
# [16] ESD_Calanoida       ESD_Chaetognatha    ESD_Cnidaria       
# [19] ESD_Copepoda        ESD_Cyclopoida     
# [22] ESD_other_Copepoda 
# [25] ESD_Protista        ESD_Pteropoda      
# [28] ESD_Urochordata     ESD_Zooplankton
unique(best$resp)
esd.resp <- c("ESD_Zooplankton","ESD_Copepoda","ESD_Protista","ESD_Cnidaria","ESD_Chaetognatha","ESD_Pteropoda","ESD_Urochordata","ESD_Calanoida","ESD_Cyclopoida","ESD_other_Copepoda")

table[which(table$resp == "ESD_Cyclopoida" & table$net == "bongo"),] 

### A°) For Bongo
bestbest <- best[which(best$resp %in% esd.resp & best$net == "regent"),] 
dim(bestbest)
# quartz()
# ggplot(data = bestbest, aes(x = factor(preds), y = rank, fill = factor(preds))) +
#     geom_boxplot(colour = "black") + xlab("Predictor") + ylab("Rank in GAM models") +
#     scale_fill_brewer(name = "", palette = "Paired") +
#     theme_classic()

unique(bestbest$preds)

bestbest[which(bestbest$resp == "ESD_Pteropoda"),]
mean(bestbest[which(bestbest$resp == "ESD_Pteropoda"),"R2"]) ; sd(bestbest[which(bestbest$resp == "ESD_Pteropoda"),"R2"])
bestbest[which(bestbest$resp == "ESD_Pteropoda"),c("preds","pval")]
ranks <- data.frame(bestbest[which(bestbest$resp == "ESD_Pteropoda"),] %>% group_by(preds) %>% summarize(pval = mean(pval)) ) # eo ddf
ranks[order(ranks$pval, decreasing = F),]


### Define the response variables to draw the resp curves for the Bongo dataset
esd.resp <- unique(best$resp)[c(16:29)]
# Filter unwanted ones
bestbest <- best[which(best$resp %in% esd.resp),] 
unique(bestbest$preds) # select values at 10m only (1:13)
bestbest <- bestbest[which(bestbest$preds %in% unique(bestbest$preds)[c(1:13)]),]
# And from that plot distribution (violins/ ggridges/ boxplots) of the predictors' ranks/Fvalue
quartz()
ggplot(data = bestbest[bestbest$net == "bongo",], aes(x = factor(preds), y = rank, fill = factor(preds))) + 
    geom_boxplot(colour = "black") + xlab("Predictor") + ylab("Rank in GAM models") + 
    scale_fill_brewer(name = "", palette = "Paired") +
    theme_classic()
#
quartz()
ggplot(data = bestbest[bestbest$net == "regent",], aes(x = factor(preds), y = rank, fill = factor(preds))) + 
    geom_boxplot(colour = "black") + xlab("Predictor") + ylab("Rank in GAM models") + 
    scale_fill_brewer(name = "", palette = "Paired") +
    theme_classic()
#
quartz()
ggplot(data = bestbest[bestbest$net == "wp2",], aes(x = factor(preds), y = rank, fill = factor(preds))) + 
    geom_boxplot(colour = "black") + xlab("Predictor") + ylab("Rank in GAM models") + 
    scale_fill_brewer(name = "", palette = "Paired") +
    theme_classic()
    

