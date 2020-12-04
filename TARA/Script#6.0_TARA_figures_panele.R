
### 03/07/2020: Script to plot the final figures for the revised manuscript
### Aims to plot the following: 
### - Figure 1: panel of spatial patterns in groups' abundance (WP2)
### - Figure 2: --> Manoela's donut plots
### - Figure 3: panel of spatial patterns in groups' ESD (WP2)
### - Figure 4: correlation heatmaps (ESD/ abund?)
### - Figures 5 and 6: Panel of GAMs response curve --> Script # 5.0

### Last update: 04/08/20

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

# --------------------------------------------------------------------------------------------------------------------------------

### Figure 1: Panel of chosen groups abundance patterns (adjust scales, labels etc.)
setwd(WD)
abund <- read.table("table_abund+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
esd <- read.table("table_ESD+hydro_allnets_18_06_20.txt", sep = "\t", h = T)

# List below the groups displaying significant latitudinal terms in the spatial GAMs
# Compare N stations per net
length(unique(abund[abund$net == "wp2","Station"])) # 152
length(unique(abund[abund$net == "bongo","Station"])) # 117
length(unique(abund[abund$net == "regent","Station"])) # 120


### WP2 -----------------------------------------------------------------------------
## Abundances: Acartiidae, Augaptilidae, Calanidae, Calanoida, Candaciidae, Centro, Chaetognatha
##              Copepoda (unid), Copepoda, Corycaeidae, Cyclopoida/Oitho, Eucalanidae, Euchaetidae, 
##              Tunicata, Oncaeidae, Paracalanidae, Poecilostomatoida, Rhizaria, Sapphirinidae, Temo
##              Small grazers and Zooplankton

### Bongo ---------------------------------------------------------------------------
## Abundances: Augaptilidae, Calanidae, Calanoida, Calanoida unid, Copepoda, Corycaeidae, Eucalanidae, 
##              Euchaetidae, Tunicata, Oncaeidae, Paracala, Rhizaria, Zooplankton
## RE-COMPUTE Poecilo (Oncaeidae + Cory)

### Régent --------------------------------------------------------------------------
## Abundances: Augaptilidae, Calanidae, Calanoida, Calanoida (unid), Candaciidae, Chaeto, Copepoda,
##              Eumalaca, Eucalanidae, Euchaetidae, Cnidaria, Metri, Rhinca, Rhizaria, Sapph, Zooplankton

### --> show WP2 for main Figure (because more signif patterns for broader groups and because it has 35 more stations) 
var <- "Ab_Copepoda"
data <- abund

distrib.analyzer <- function(var, data) {
    
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
                        title2 <- paste(title," abundance", sep = "")
                   
                        nstations <- nrow(subset)           
                     
                        ### Do we even have more than 30 stations?
                        if(nstations >= 30) {
                            
                            #lab <- paste0("Abundance<br>(ind.m<sup>2</sup>)")                   
                            # Map 'em
                            map1 <- ggplot() + geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey85", colour = "grey50", size = 0.3) +
                                geom_point(aes(x = Longitude, y = Latitude, fill = get(var)), data = subset, colour = "black", pch = 21, size = 2) + 
                                scale_fill_distiller(name = paste0("Abundance"), palette = "Spectral") + 
                            	coord_quickmap() + scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,180),
                                       	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                            	scale_y_continuous(name = "", labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N"), limits = c(-90,90), 
                            	      	breaks = c(-90,-60,-30,0,30,60,90), expand = c(0,0)) +
                            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                            			panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
                                        ggtitle(title)
                        
                            # Fit gams if enough stations without NAs
                            require("mgcv")
                            gam1 <- mgcv::gam(data = subset, get(var) ~ s(Latitude, bs="tp"), method = "REML")
                            # Extract deviance explained : str(summary(gam1))
                            r2.1 <- round(summary(gam1)$r.sq,2)
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
                            library("ggtext")
                            library("ggpubr")
                            plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = fit - se, xmax = fit + se), fill = "black", data = pred1, alpha = 0.25) +
                                        geom_path(aes(y = y, x = fit), data = pred1, colour = "black") +
                                        scale_y_continuous(position = "right", limits = c(-90,90), name = "Latitude", expand = c(0,0),
                                            breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N")) + 
                                        xlab( expression("Fitted abundance" ~ (ind.m^{2})) ) +
                                        annotate("text", size = 4.5, x = max(pred1$fit)-0.2, y = -65, label = as.expression(bquote(R^2~"="~.(r2.1))) ) + 
                                        annotate("text", size = 4.5, x = max(pred1$fit)-0.2, y = -75, label = paste("n = ",nstations,sep="") ) + 
                                        annotate("text", size = 4.5, x = max(pred1$fit)-0.2, y = -85, label = paste(sign,sep="")) + 
                                        theme_classic()   
                                
                            # Save plots
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Latitudinal_gradients/", sep = ""))
                            #panel <- ggarrange(map1, plot1, ncol = 2, nrow = 1, widths = c(3,1))
                            ggsave(plot = map1, filename = paste("map_",var,"_",net,".pdf", sep = ""), dpi = 300, width = 7, height = 5)       
                            ggsave(plot = plot1, filename = paste("plot_zonal_",var,"_",net,".pdf", sep = ""), dpi = 300, width = 3, height = 4)       

                                
                            } else {
                                
                                message(paste("NOT ENOUGH STATIONS", sep = ""))
                                
                            }
                                                       
} # eo fun

# Apply fun above in for loop
colnames(abund)[c(5:9,12)] <- c("Ab_Zooplankton","Ab_Copepoda","Ab_Chaetognatha","Ab_Cnidaria","Ab_Tunicata","Ab_Rhizaria")
variables <- colnames(abund)[c(5:9,12)] ; variables
net <- "wp2"
for(var in variables) {
    distrib.analyzer(var = var, data = abund)   
}

### --> Adapt same code for the other groups and nets --> SI



# -----------------------------------------------------------------------------------

### Figure 3: Panel of chosen groups median ESD patterns (adjust scales, labels etc.)

### WP2: Copepoda, Calanoida, Poecilo, Cnidaria, Rhizaria, total zoo weakly signf.
### Bongo: Copepoda, Calanoida, Heterorhabdidae, small grazers, Zooplankon
### Régent: Copepoda, Calanoida, Copepoda (unid), Cnidaria, Euchaetidae+Heterorhabdidae, small grazers and Zooplankton

#Need the counts for this one
counts.bongo <- read.table("table_counts_groups_bongo.txt", h = T, sep = "\t")
counts.wp2 <- read.table("table_counts_groups_wp2.txt", h = T, sep = "\t")
counts.regent <- read.table("table_counts_groups_regent.txt", h = T, sep = "\t")

### --> show a panel of median Copepoda ESD only for all 3 nets 
data <- esd
var <- "ESD_Copepoda"
net <- "bongo"

### You want to use the same palette scale this time: check ranges of median ESD across nets
#summary(data[data$net == "wp2","ESD_Copepoda"])
#summary(data[data$net == "bongo","ESD_Copepoda"])
#summary(data[data$net == "regent","ESD_Copepoda"])

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
                        title <- str_replace(var,"ESD_","")
                        #title1 <- paste(title,"\nmedian abundance", sep = "")
                        #title2 <- paste("Fitted ",title," median abundance", sep = "")
                        if(net == "wp2") {
                            titlemap <- "WP2"
                        } else if(net == "bongo") {
                            titlemap <- "Bongo"
                        } else {
                            titlemap <- "Régent"
                        }
                   
                        if(var == "ESD_Zooplankton") {
                            
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
                                scale_fill_distiller(name = "Median ESD\nlog(µm)", palette = "Spectral") + 
                            	coord_quickmap() + scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,180),
                                       	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                            	scale_y_continuous(name = "", labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N"), limits = c(-90,90), 
                            	      	breaks = c(-90,-60,-30,0,30,60,90), expand = c(0,0)) +
                            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                            			panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
                                        ggtitle(titlemap)
                        
                            # Fit gams if enough stations without NAs
                            if( nrow(subset2) >= 30 ) {
                                
                                require("mgcv")
                                gam1 <- mgcv::gam(data = subset2, get(var) ~ s(Latitude, bs="tp"), method = "REML")
                                # Extract deviance explained : str(summary(gam1))
                                r2.1 <- round(summary(gam1)$r.sq,2)
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
                                library("ggtext")
                                library("ggpubr")
                                if(net != "regent") {
                                    
                                    plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = fit - se, xmax = fit + se), fill = "black", data = pred1, alpha = 0.25) +
                                                geom_path(aes(y = y, x = fit), data = pred1, colour = "black") +
                                                scale_y_continuous(position = "right", limits = c(-90,90), name = "Latitude", expand = c(0,0),
                                                    breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N")) + 
                                                xlab( expression("Fitted median ESD" ~ log(µm)) ) +
                                                annotate("text", size = 4.5, x = max(pred1$fit)-0.03, y = -65, label = as.expression(bquote(R^2~"="~.(r2.1))) ) + 
                                                annotate("text", size = 4.5, x = max(pred1$fit)-0.03, y = -75, label = paste("n = ",nstations,sep="") ) + 
                                                annotate("text", size = 4.5, x = max(pred1$fit)-0.03, y = -85, label = paste(sign,sep="")) + 
                                                theme_classic()   
                                    
                                } else {
                                    
                                    plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = fit - se, xmax = fit + se), fill = "black", data = pred1, alpha = 0.25) +
                                                geom_path(aes(y = y, x = fit), data = pred1, colour = "black") +
                                                scale_y_continuous(position = "right", limits = c(-90,90), name = "Latitude", expand = c(0,0),
                                                    breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N")) + 
                                                xlab( expression("Fitted median ESD" ~ log(µm)) ) +
                                                annotate("text", size = 4.5, x = min(pred1$fit)+0.05, y = -65, label = as.expression(bquote(R^2~"="~.(r2.1))) ) + 
                                                annotate("text", size = 4.5, x = min(pred1$fit)+0.05, y = -75, label = paste("n = ",nstations,sep="") ) + 
                                                annotate("text", size = 4.5, x = min(pred1$fit)+0.05, y = -85, label = paste(sign,sep="")) + 
                                                theme_classic()   
                                    
                                }
                                
                                # Save plots
                                setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Latitudinal_gradients/", sep = ""))
                                #panel <- ggarrange(map1, plot1, ncol = 2, nrow = 1, widths = c(3,1))
                                ggsave(plot = map1, filename = paste("map_",var,"_",net,".pdf", sep = ""), dpi = 300, width = 7, height = 5)       
                                ggsave(plot = plot1, filename = paste("plot_zonal_",var,"_",net,".pdf", sep = ""), dpi = 300, width = 3, height = 4)       
                                
                            } else {
                                
                                message(paste("NOT ENOUGH STATIONS", sep = ""))
                                
                            }
         
                            
                        } else {
                            
                            message(paste("NOT ENOUGH STATIONS", sep = ""))
                            
                        }
                                               
            } # eo for loop - net types
        
} # eo fun
# Apply fun above in for loop
variables <- "ESD_Copepoda"
for(var in variables) {
    distrib.analyzer(var = var, data = esd)   
}


### 15/07/20: Compute rank corr coeff between abund and median ESD for total zooplankton
counts.bongo <- read.table("table_counts_groups_bongo.txt", h = T, sep = "\t")
counts.wp2 <- read.table("table_counts_groups_wp2.txt", h = T, sep = "\t")
counts.regent <- read.table("table_counts_groups_regent.txt", h = T, sep = "\t")

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
    
        message(paste("Zooplankton ab vs ESD", sep = ""))      
        gr <- "Zooplankton"      
        gr.size <- sizes[,c("Station",gr)]
        gr.abund <- abs[,c("Station",gr)]
        # na.omit
        gr.size <- na.omit(gr.size)
        gr.abund <- na.omit(gr.abund)
        # And order per station
        gr.size <- gr.size[order(gr.size$Station),]
        gr.abund <- gr.abund[order(gr.abund$Station),]
                  
        # Restrict co common stations to compute correlation
        common_stations <- intersect(unique(gr.size$Station), unique(gr.abund$Station))
        gr.size2 <- gr.size[which(gr.size$Station %in% common_stations),]
        gr.abund2 <- gr.abund[which(gr.abund$Station %in% common_stations),]
        # There can be several abund measurements for abund
        gr.abund3 <- data.frame(gr.abund2 %>% group_by(Station) %>% summarize(abund = mean(!! sym(gr)) ) ) # noice
        # If both esd and abund gave same dimensios and theer's enough stations (n = 10)
        dim(gr.size) ; dim(gr.abund3)
                        
        rho <- cor(gr.size2[,2], gr.abund3[,2], method = "spearman")
        pval <- cor.test(x = gr.size2[,2], y = gr.abund3[,2], method = "spearman")$p.value
        # Put in a ddf
        res.cor <- data.frame(group = gr, net = net, n = nrow(gr.size2), rho = rho, pval = pval)
        return(res.cor)
        
    } # eo 1st FUN
    
) # eo 1st LAPPLY
ddf <- bind_rows(res)
dim(ddf)
rm(res)

#         group    net   n       rho         pval
# 1 Zooplankton  bongo 117 0.1736911 6.108950e-02
# 2 Zooplankton    wp2 152 0.3915662 6.089625e-07
# 3 Zooplankton regent 120 0.3488553 9.425303e-05

# -----------------------------------------------------------------------------------

### Figure 4: Correlation heatmaps between Hydro and abundance and then hydro and median ESD

library("corrplot")
library("corrgram")
library("ggcorrplot")
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


### A) Abund vs environmental covariates
abund <- read.table("table_abund+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
esd <- read.table("table_ESD+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
dim(abund) ; dim(esd)
# colnames(abund) --> 5:45 + 47,48 + 50:54 + 56:59
# colnames(esd) --> 5:45 + 47,48 + 50:54 + 56:59
net <- "wp2"

### Draw the correlation coeff heatmaps for each net

for(net in unique(abund$net)) {
    
        message(paste("", sep = ""))
        message(paste("Examining covariance between abundance and size for ",net, sep = ""))
        message(paste("", sep = ""))
        names <- colnames(abund)[c(5:21,23:27,29:37,43:45,47:59)] #; names
        abs <- abund[abund$net == net,names]
        # summary(abs)
        # 2nd lapply based on the groups 
        groups <- colnames(abs)[c(1:31)]
        groups <- str_replace_all(groups,"Ab_","")
        colnames(abs)[c(1:31)] <- groups
        
        # Remove columns that ONLY have NAs (zoo groups)
        abs <- abs %>% select_if(~sum(!is.na(.)) > 0)
        mydata <- na.omit(abs)

        cormat <- round(cor(mydata, method = "spearman"),2)
        p.mat <- cor_pmat(mydata, method = "spearman", conf.level = 0.95)
        # dim(cormat) ; dim(p.mat)
        # Re-order?
        cor_tri <- get_upper_tri(cormat)
        cor_tri <- melt(cor_tri, na.rm = T)
        pval_tri <- get_upper_tri(p.mat)
        pval_tri <- melt(pval_tri, na.rm = T)
        colnames(cor_tri)[3] <- "rho"
        colnames(pval_tri)[3] <- "pval"
        cor_tri$pval <- pval_tri$pval
        rm(pval_tri,p.mat,cormat) ; gc()
        # unique(cor_tri$Var1)
        # unique(cor_tri$Var2)
        
        ### First heatmap: main zoo groups versus hydro
        zoo_groups <- c("Zooplankton","Copepoda","Chaetognatha","Gel_carn","Gel_FF","Crust","Pteropoda","Rhizaria","Small_grazers","Calanoida","Poecilostomatoida","Cyclopoida")
        hydro_vars <- c("Temperature","O2","Salinity","MLD","Chla","bbp470","Micro","Nano","Pico","PAR2","NO2NO3","PO4","SiO2")
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
         
        ### Add signif label
        cormat$signif <- NA
        cormat[cormat$pval <= 0.001,"signif"] <- "***"
        cormat[cormat$pval > 0.001 & cormat$pval <= 0.01,"signif"] <- "**"
        cormat[cormat$pval > 0.01 & cormat$pval < 0.05,"signif"] <- "*"
        cormat[cormat$pval >= 0.05,"signif"] <- "-"
        
        # Change Net name for the ggtitle
        if(net == "wp2") {
            titlemap <- "WP2 (200µm)"
        } else if(net == "bongo") {
            titlemap <- "Bongo (300µm)"
        } else {
            titlemap <- "Régent (680µm)"
        } 
        
        colnames(cormat)[c(1,2)] <- c("Group","Covariate")
        
        plot <- ggplot(cormat, aes(factor(Covariate), factor(Group), fill = rho)) + geom_tile(color = "white") +
                scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-0.65,0.65), 
                        name = paste("Spearman' rank\ncorrelation\ncoefficient", sep = ""),
                        breaks = c(-0.65,-0.45,-0.25,0,0.25,0.45,0.65), 
                        guide = guide_colourbar(ticks = T, ticks.colour = "black", frame.colour = "black")) + 
                xlab("") + ylab("") + 
                theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
                        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))+
                coord_fixed() + geom_text(aes(factor(Covariate), factor(Group), label = signif), color = "black", size = 4) +
                ggtitle(titlemap)
            
        
        # Save plots
        ggsave(plot = plot, filename = paste("heatmap_corr_abund_hydro_",net,"_V1.pdf", sep = ""), dpi = 300, width = 8, height = 8)
        
} # eo for loop
    
    

### B) median ESD vs environmental covariates

### In a lapply, per group and net, filter ESD and hydro data. Restrict to stations that have n > 25 for ESD and abund
nets <- c("bongo","wp2","regent")

net <- "wp2"
res <- lapply(nets, function(net) {
    
        message(paste("", sep = ""))
        message(paste("Examining covariance between abundance and size for ",net, sep = ""))
        message(paste("", sep = ""))
        subset <- esd[esd$net == net,]
        
        # 2nd lapply based on the groups 
        names <- colnames(subset)[c(5:42)]
        groups <- str_replace_all(names,"ESD_","")
        colnames(subset)[c(5:42)] <- groups
 
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
       
        common.groups <- c("Zooplankton",intersect(groups, unique(counts$group)))
        
        # gr <- "Copepoda"
        grps <- lapply(common.groups, function(gr) {
            
                    message(paste(gr, sep = ""))   
                    if( gr ==  "Zooplankton") {
                        
                        # Restrict to vars of intrest for corr heatmap
                        subset2 <- subset[,c("Station",gr,"Temperature","O2","Salinity","PAR2",
                                    "MLD","Chla","bbp470","Micro","Nano","Pico","NO2NO3","SiO2","PO4")]
                        # Remove rows with missing values
                        subset2 <- na.omit(subset2)
              
                    }   else {
                        
                        stations2keep <- counts[counts$group == gr & counts$n >= 20 ,"station"]
                        # Restrict to vars of intrest for corr heatmap
                        subset2 <- subset[subset$Station %in% stations2keep, c("Station",gr,"Temperature","O2","Salinity","PAR2",
                                    "MLD","Chla","bbp470","Micro","Nano","Pico","NO2NO3","SiO2","PO4")]
                        # Remove rows with missing values
                        subset2 <- na.omit(subset2)
                        
                    }      
                        
                    # dim(gr.size2) ; dim(gr.abund2)
                    if( nrow(subset2) > 10 ) {
                        
                        cormat <- round(cor(subset2[,c(2:length(subset2))], method = "spearman"),2)
                        p.mat <- cor_pmat(subset2[,c(2:length(subset2))], method = "spearman", conf.level = 0.95)
                        # dim(cormat) ; dim(p.mat)
                        # Re-order?
                        cor_tri <- get_upper_tri(cormat)
                        cor_tri <- melt(cor_tri, na.rm = T)
                        pval_tri <- get_upper_tri(p.mat)
                        pval_tri <- melt(pval_tri, na.rm = T)
                        colnames(cor_tri)[3] <- "rho"
                        colnames(pval_tri)[3] <- "pval"
                        cor_tri$pval <- pval_tri$pval
                        # Remove pairwise correlations including only env covariates
                        cor_tri <- cor_tri[cor_tri$Var1 == gr,]
                        cor_tri <- cor_tri[-which(cor_tri$Var2 == gr),]
                        colnames(cor_tri)[c(1,2)] <- c("Group","Covariate")
                        
                        # Add signif label
                        cor_tri$signif <- NA
                        cor_tri[cor_tri$pval <= 0.001,"signif"] <- "***"
                        cor_tri[cor_tri$pval > 0.001 & cor_tri$pval <= 0.01,"signif"] <- "**"
                        cor_tri[cor_tri$pval > 0.01 & cor_tri$pval < 0.05,"signif"] <- "*"
                        cor_tri[cor_tri$pval >= 0.05,"signif"] <- "-"
  
                        # Add net
                        cor_tri$net <- net
                        
                        # Return
                        return(cor_tri)
                        
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
head(ddf)

### Change some names 
ddf$Group <- factor(ddf$Group)
ddf$Covariate <- factor(ddf$Covariate)
levels(ddf$Group)[levels(ddf$Group) == "Gel_carn"] <- "Cnidaria"
levels(ddf$Group)[levels(ddf$Group) == "Gel_FF"] <- "Tunicata"
levels(ddf$Group)[levels(ddf$Group) == "Crust"] <- "Eumalacostraca"
levels(ddf$Group)[levels(ddf$Group) == "Small_grazers"] <- "Ostracoda+Cladocera"
# Same for hydrobio data
levels(ddf$Covariate)[levels(ddf$Covariate) == "O2"] <- "Oxygen"
levels(ddf$Covariate)[levels(ddf$Covariate) == "Chla"] <- "Chlorophylla"
levels(ddf$Covariate)[levels(ddf$Covariate) == "Micro"] <- "%Micro"
levels(ddf$Covariate)[levels(ddf$Covariate) == "Nano"] <- "%Nano"
levels(ddf$Covariate)[levels(ddf$Covariate) == "Pico"] <- "%Pico"
levels(ddf$Covariate)[levels(ddf$Covariate) == "PAR2"] <- "PAR"
levels(ddf$Covariate)[levels(ddf$Covariate) == "NO2NO3"] <- "Nitrate"
levels(ddf$Covariate)[levels(ddf$Covariate) == "PO4"] <- "Phosphate"
levels(ddf$Covariate)[levels(ddf$Covariate) == "SiO2"] <- "Silicate"

# Change order of factors along x a,d y axes
myzoogroups <- c("Zooplankton","Copepoda","Cnidaria","Eumalacostraca","Rhizaria","Calanoida","Cyclopoida","Poecilostomatoida")

mycovids <- c("Temperature","Oxygen","Salinity","MLD","Chlorophylla","bbp470","%Micro","%Nano","%Pico","PAR","Nitrate","Silicate","Phosphate")          


# Filter vars out: those 2 keep for main igures and those to keep for SI
for(net in nets) {
    
    subset <- ddf[ddf$net == net,]
    data2plot <- subset[subset$Group %in% myzoogroups,]
    # re-order for plot
    data2plot$Group <- factor(data2plot$Group, levels = myzoogroups)
    data2plot$Covariate <- factor(data2plot$Covariate, levels = mycovids)
    
    # Change Net name for the ggtitle
    if(net == "wp2") {
        titlemap <- "WP2 (200µm)"
    } else if(net == "bongo") {
        titlemap <- "Bongo (300µm)"
    } else {
        titlemap <- "Régent (680µm)"
    } 

    plot <- ggplot(data2plot, aes(factor(Covariate), factor(Group), fill = rho)) + geom_tile(color = "white") +
                scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-0.85,0.85), 
                        name = paste("Spearman' rank\ncorrelation\ncoefficient", sep = ""),
                        breaks = c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8), 
                        guide = guide_colourbar(ticks = T, ticks.colour = "black", frame.colour = "black")) + 
                xlab("") + ylab("") + 
                theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
                        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))+
                coord_fixed() + geom_text(aes(factor(Covariate), factor(Group), label = signif), color = "black", size = 4) +
                ggtitle(titlemap)
                
    ggsave(plot = plot, filename = paste("heatmap_corr_ESD_hydro_",net,"_V2.pdf", sep = ""), dpi = 300, width = 8, height = 8)            
    
}


### Check which groups display the main patterns
ddf[ddf$Group %in% myzoogroups & ddf$signif != "-",]
# Check N signif relationships
summary(factor(ddf[ddf$Group %in% myzoogroups & ddf$signif != "-","net"]))
#  bongo regent    wp2 
#    38     45     40 

   
# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------

### 31/07/20: For preparing supplementary materials
library("ggrepel")
### Figure S1: 
abund <- read.table("table_abund+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
dim(abund); str(abund)

abund$net <- factor(abund$net)
levels(abund$net)[levels(abund$net) == "regent"] <- "Régent"
levels(abund$net)[levels(abund$net) == "bongo"] <- "Bongo"
levels(abund$net)[levels(abund$net) == "wp2"] <- "WP2"

#quartz()
map <- ggplot() + geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey85", colour = "grey50", size = 0.3) +
  	geom_point(aes(x = Longitude, y = Latitude, shape = net, fill = net), data = abund) +
  	scale_shape_manual(name = "", labels = c("Bongo","Régent","WP2"), values = c(21:23)) + 
    scale_fill_manual(name = "", labels = c("Bongo","Régent","WP2"), values = c("#1f78b4","#33a02c","#e31a1c")) + 
    geom_text_repel(aes(x = Longitude, y = Latitude, label = Station), data = abund, size = 2.5) + 
    coord_quickmap() + scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,180),
           	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
    scale_y_continuous(name = "", labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N"), limits = c(-90,90), 
        	breaks = c(-90,-60,-30,0,30,60,90), expand = c(0,0)) +
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
  		panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
    facet_wrap(factor(net) ~. , nrow = 3, ncol = 1)

ggsave(plot = map, filename = "Fig.S1_map_stations.pdf", dpi = 300, width = 10, height = 14)

# -----------------------------------------------------------------------------------

### Figure S10: heatmaps of Spearman's rank correlation between environmental covariates per net

abund <- read.table("table_abund+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
dim(abund); str(abund)

net <- "wp2"

for(net in unique(abund$net)) {
    
        message(paste("", sep = ""))
        message(paste("Examining covariance between abundance and size for ",net, sep = ""))
        message(paste("", sep = ""))
        names <- colnames(abund)[c(43:45,47:53,55:57)] # names 2 keep (env covariates)       
        mydata <- na.omit(abund[abund$net == net,names])

        cormat <- round(cor(mydata, method = "spearman"),2)
        p.mat <- cor_pmat(mydata, method = "spearman", conf.level = 0.95)
        # dim(cormat) ; dim(p.mat)
        # Re-order?
        cor_tri <- get_upper_tri(cormat)
        cor_tri <- melt(cor_tri, na.rm = T)
        pval_tri <- get_upper_tri(p.mat)
        pval_tri <- melt(pval_tri, na.rm = T)
        colnames(cor_tri)[3] <- "rho"
        colnames(pval_tri)[3] <- "pval"
        cor_tri$pval <- pval_tri$pval
        rm(pval_tri,p.mat,cormat) ; gc()
        # unique(cor_tri$Var1)
        # unique(cor_tri$Var2)
        
        cormat <- cor_tri[!(cor_tri$Var1 == cor_tri$Var2),]
        
        # Same for hydrobio data
        levels(cormat$Var1)[levels(cormat$Var1) == "O2"] <- "Oxygen"
        levels(cormat$Var1)[levels(cormat$Var1) == "Chla"] <- "Chlorophyll a"
        levels(cormat$Var1)[levels(cormat$Var1) == "Micro"] <- "%Micro"
        levels(cormat$Var1)[levels(cormat$Var1) == "Nano"] <- "%Nano"
        levels(cormat$Var1)[levels(cormat$Var1) == "Pico"] <- "%Pico"
        levels(cormat$Var1)[levels(cormat$Var1) == "PAR2"] <- "PAR"
        levels(cormat$Var1)[levels(cormat$Var1) == "NO2NO3"] <- "Nitrate"
        levels(cormat$Var1)[levels(cormat$Var1) == "PO4"] <- "Phosphate"
        levels(cormat$Var1)[levels(cormat$Var1) == "SiO2"] <- "Silicate"
        
        levels(cormat$Var2)[levels(cormat$Var2) == "O2"] <- "Oxygen"
        levels(cormat$Var2)[levels(cormat$Var2) == "Chla"] <- "Chlorophyll a"
        levels(cormat$Var2)[levels(cormat$Var2) == "Micro"] <- "%Micro"
        levels(cormat$Var2)[levels(cormat$Var2) == "Nano"] <- "%Nano"
        levels(cormat$Var2)[levels(cormat$Var2) == "Pico"] <- "%Pico"
        levels(cormat$Var2)[levels(cormat$Var2) == "PAR2"] <- "PAR"
        levels(cormat$Var2)[levels(cormat$Var2) == "NO2NO3"] <- "Nitrate"
        levels(cormat$Var2)[levels(cormat$Var2) == "PO4"] <- "Phosphate"
        levels(cormat$Var2)[levels(cormat$Var2) == "SiO2"] <- "Silicate"
         
        ### Add signif label
        cormat$signif <- NA
        cormat[cormat$pval <= 0.001,"signif"] <- "***"
        cormat[cormat$pval > 0.001 & cormat$pval <= 0.01,"signif"] <- "**"
        cormat[cormat$pval > 0.01 & cormat$pval < 0.05,"signif"] <- "*"
        cormat[cormat$pval >= 0.05,"signif"] <- "-"
        
        # Change Net name for the ggtitle
        if(net == "wp2") {
            titlemap <- "WP2 (200µm)"
        } else if(net == "bongo") {
            titlemap <- "Bongo (300µm)"
        } else {
            titlemap <- "Régent (680µm)"
        } 
        
        colnames(cormat)[c(1,2)] <- c("Covariate1","Covariate2")
        
        plot <- ggplot(cormat, aes(factor(Covariate2), factor(Covariate1), fill = rho)) + geom_tile(color = "white") +
                scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), 
                        name = paste("Spearman' rank\ncorrelation\ncoefficient", sep = ""),
                        breaks = c(-1.0,-0.75,-0.50,-0.25,0,0.25,0.50,0.75,1.0), 
                        guide = guide_colourbar(ticks = T, ticks.colour = "black", frame.colour = "black")) + 
                xlab("") + ylab("") + 
                theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
                        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))+
                coord_fixed() + geom_text(aes(factor(Covariate2), factor(Covariate1), label = rho), color = "black", size = 3) +
                ggtitle(titlemap)            
        
        # Save plots
        ggsave(plot = plot, filename = paste("heatmap_corr_hydro_",net,".pdf", sep = ""), dpi = 300, width = 8, height = 8)
        
} # eo for loop
    

# -----------------------------------------------------------------------------------

### Figure S2: panel of cubic-transformed abundances

summary(abund[abund$net == "bongo","Ab_Poecilostomatoida"])
abund[abund$net == "bongo","Ab_Poecilostomatoida"] <- (abund[abund$net == "bongo","Ab_Corycaeidae"]) + (abund[abund$net == "bongo","Ab_Oncaeidae"])
summary(abund[abund$net == "bongo","Ab_Poecilostomatoida"])

### Same with other nets just to be sure
summary(abund[abund$net == "wp2","Ab_Sapphirinidae"])
summary(abund[abund$net == "regent","Ab_Sapphirinidae"])


abund[abund$net == "wp2","Ab_Poecilostomatoida"] <- (abund[abund$net == "wp2","Ab_Corycaeidae"]) + (abund[abund$net == "wp2","Ab_Oncaeidae"]) + (abund[abund$net == "wp2","Ab_Sapphirinidae"])
abund[abund$net == "regent","Ab_Poecilostomatoida"] <- (abund[abund$net == "regent","Ab_Corycaeidae"]) + (abund[abund$net == "regent","Ab_Oncaeidae"]) + (abund[abund$net == "regent","Ab_Sapphirinidae"])
summary(abund[abund$net == "wp2","Ab_Poecilostomatoida"])
summary(abund[abund$net == "regent","Ab_Poecilostomatoida"])

data <- abund

var <- "Ab_Poecilostomatoida"
net <- "wp2"

distrib.analyzer <- function(var, data) {
    
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
                        title2 <- paste(title," abundance", sep = "")
                   
                        nstations <- nrow(subset)        
                        
                        if(net == "wp2") {
                            net2 <- "WP2 (200µm)"
                        } else if(net == "bongo") {
                            net2 <- "Bongo (300µm)"
                        } else {
                            net2 <- "Régent (680µm)"
                        }
                      
                     
                        ### Do we even have more than 30 stations?
                        if( nstations >= 30 ) {
                            
                            #lab <- paste0("Abundance<br>(ind.m<sup>2</sup>)")                   
                            # Map 'em
                            map1 <- ggplot() + geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey85", colour = "grey50", size = 0.3) +
                                geom_point(aes(x = Longitude, y = Latitude, fill = get(var)), data = subset, colour = "black", pch = 21, size = 2) + 
                                scale_fill_distiller(name = paste0("Abundance"), palette = "Spectral") + 
                            	coord_quickmap() + scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,180),
                                       	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                            	scale_y_continuous(name = "", labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N"), limits = c(-90,90), 
                            	      	breaks = c(-90,-60,-30,0,30,60,90), expand = c(0,0)) +
                            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                            			panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
                                        ggtitle(paste(title," - ",net2, sep = ""))
                        
                            # Fit gams if enough stations without NAs
                            require("mgcv")
                            gam1 <- mgcv::gam(data = subset, get(var) ~ s(Latitude, bs="tp"), method = "REML")
                            # Extract deviance explained : str(summary(gam1))
                            r2.1 <- round(summary(gam1)$r.sq,2)
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
                            library("ggtext")
                            library("ggpubr")
                            maxfit <- max(pred1$fit)
                            
                                
                            plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = fit - se, xmax = fit + se), fill = "black", data = pred1, alpha = 0.25) +
                                            geom_path(aes(y = y, x = fit), data = pred1, colour = "black") +
                                            scale_y_continuous(position = "right", limits = c(-90,90), name = "Latitude", expand = c(0,0),
                                                breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N")) + 
                                            xlab( expression("Fitted abundance" ~ (ind.m^{3})) ) +
                                            annotate("text", size = 3, x = maxfit - (0.1*(maxfit)), y = -65, label = as.expression(bquote(R^2~"="~.(r2.1))) ) + 
                                            annotate("text", size = 3, x = maxfit - (0.1*(maxfit)), y = -75, label = paste("n = ",nstations,sep="") ) + 
                                            annotate("text", size = 3, x = maxfit - (0.1*(maxfit)), y = -85, label = paste(sign,sep="")) + 
                                            theme_classic()
                                
                            # Save plots
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Latitudinal_gradients/", sep = ""))
                            panel <- ggarrange(map1, plot1, ncol = 2, nrow = 1, widths = c(2.75,1))
                            #ggsave(plot = map1, filename = paste("map_",var,"_",net,".jpg", sep = ""), dpi = 300, width = 6, height = 4)       
                            #ggsave(plot = plot1, filename = paste("plot_zonal_",var,"_",net,".jpg", sep = ""), dpi = 300, width = 3, height = 4)   
                            ggsave(plot = panel, filename = paste("plot_zonal_",var,"_",net,".jpg", sep = ""), dpi = 300, width = 8.5, height = 3)           

                                
                            } else {
                                
                                message(paste("NOT ENOUGH STATIONS", sep = ""))
                                
                            }
                                                       
} # eo fun

### Apply fun above in for loop: whoch columns to keep?
# c(10,11,13:42)
colnames(abund)[c(10,13,14,22)] <- c("Ab_Eumalacostraca","Ab_Ostracoda+Cladocera","Ab_Copepoda (unidentified)","Ab_Calanoida (unidentified)")

variables <- colnames(abund)[c(10,11,13:27,29:37)] ; variables

net <- "bongo"

for(var in variables) {
    distrib.analyzer(var = var, data = abund)   
}

### Assemble manually on a panel (Fig. S2)


### Figure S4: panel of log-transformed median ESD
setwd(WD)
counts.bongo <- read.table("table_counts_groups_bongo.txt", h = T, sep = "\t")
counts.wp2 <- read.table("table_counts_groups_wp2.txt", h = T, sep = "\t")
counts.regent <- read.table("table_counts_groups_regent.txt", h = T, sep = "\t")

data <- esd
var <- "ESD_Rhizaria"
net <- "bongo"

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
                        levels(counts$group)[levels(counts$group) == "Calanoida_unid"] <- "Calanoida (unidentified)"
                        levels(counts$group)[levels(counts$group) == "Copepoda_unid"] <- "Copepoda (unidentified)"
                        levels(counts$group)[levels(counts$group) == "Crustaceans"] <- "Eumalacostraca"
                        levels(counts$group)[levels(counts$group) == "Grazers"] <- "Ostracoda+Cladocera"
                        levels(counts$group)[levels(counts$group) == "Gel_carn"] <- "Cnidaria"
                        levels(counts$group)[levels(counts$group) == "Gel_FF"] <- "Tunicata"
   
                        # Change varname for title/ axes caption
                        title <- str_replace(var,"ESD_","")
                        
                        #title1 <- paste(title,"\nmedian abundance", sep = "")
                        #title2 <- paste("Fitted ",title," median abundance", sep = "")
                        if(net == "wp2") {
                            titlemap <- "WP2 (200µm)"
                        } else if(net == "bongo") {
                            titlemap <- "Bongo (300µm)"
                        } else {
                            titlemap <- "Régent (680µm)"
                        }
                   
                        if(var == "ESD_Zooplankton") {
                            
                             stations2keep <- unique(subset$Station)
                             nstations <- length(stations2keep)
                             subset2 <- subset
                             
                        } else {
                            
                             stations2keep <- counts[counts$group == title & counts$n >= 20 ,"station"]
                             nstations <- length(stations2keep)
                             subset2 <- subset[which(subset$Station %in% stations2keep),]
                             
                        } #                     
                     
                        ### Do we even have more than 30 stations?
                        if(nstations >= 30) {
                            
                            # Map 'em
                            map1 <- ggplot() + geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey85", colour = "grey50", size = 0.3) +
                                geom_point(aes(x = Longitude, y = Latitude, fill = get(var)), data = subset2, colour = "black", pch = 21, size = 2) + 
                                scale_fill_distiller(name = "Median ESD", palette = "Spectral") + 
                            	coord_quickmap() + scale_x_continuous(name = "", breaks = c(-180,-120,-60,0,60,120,180),
                                       	labels = c("180°W","120°W","60°W","0°W","60°E","120°E","180°E"), expand = c(0,0)) +
                            	scale_y_continuous(name = "", labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N"), limits = c(-90,90), 
                            	      	breaks = c(-90,-60,-30,0,30,60,90), expand = c(0,0)) +
                            	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                            			panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
                                        ggtitle(paste(title," - ",titlemap, sep = ""))
                        
                            # Fit gams if enough stations without NAs
                            if( nrow(subset2) >= 30 ) {
                                
                                require("mgcv")
                                gam1 <- mgcv::gam(data = subset2, get(var) ~ s(Latitude, bs="tp"), method = "REML")
                                # Extract deviance explained : str(summary(gam1))
                                r2.1 <- round(summary(gam1)$r.sq,2)
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
                                library("ggtext")
                                library("ggpubr")
                                
                                if( net != "regent" ) {
                                    
                                    maxfit <- max(pred1$fit)
                                    
                                    plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = fit - se, xmax = fit + se), fill = "black", data = pred1, alpha = 0.25) +
                                                geom_path(aes(y = y, x = fit), data = pred1, colour = "black") +
                                                scale_y_continuous(position = "right", limits = c(-90,90), name = "Latitude", expand = c(0,0),
                                                    breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N")) + 
                                                xlab( expression("Fitted median ESD" ~ log(µm)) ) +
                                                annotate("text", size = 3, x = maxfit-(0.01*(maxfit)), y = -65, label = as.expression(bquote(R^2~"="~.(r2.1))) ) + 
                                                annotate("text", size = 3, x = maxfit-(0.01*(maxfit)), y = -75, label = paste("n = ",nstations,sep="") ) + 
                                                annotate("text", size = 3, x = maxfit-(0.01*(maxfit)), y = -85, label = paste(sign,sep="")) + 
                                                theme_classic()   
                                    
                                } else {
                                    
                                    minfit <- min(pred1$fit)
                                    
                                    plot1 <- ggplot() + geom_ribbon(aes(y = y, xmin = fit - se, xmax = fit + se), fill = "black", data = pred1, alpha = 0.25) +
                                                geom_path(aes(y = y, x = fit), data = pred1, colour = "black") +
                                                scale_y_continuous(position = "right", limits = c(-90,90), name = "Latitude", expand = c(0,0),
                                                    breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","0°N","30°N","60°N","90°N")) + 
                                                xlab( expression("Fitted median ESD" ~ log(µm)) ) +
                                                annotate("text", size = 3, x = minfit+(0.01*(minfit)), y = -65, label = as.expression(bquote(R^2~"="~.(r2.1))) ) + 
                                                annotate("text", size = 3, x = minfit+(0.01*(minfit)), y = -75, label = paste("n = ",nstations,sep="") ) + 
                                                annotate("text", size = 3, x = minfit+(0.01*(minfit)), y = -85, label = paste(sign,sep="")) + 
                                                theme_classic()   
                                    
                                }
                                
                                # Save plots
                                setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Latitudinal_gradients/", sep = ""))
                                panel <- ggarrange(map1, plot1, ncol = 2, nrow = 1, widths = c(2.75,1))
                                ggsave(plot = panel, filename = paste("plot_zonal_",var,"_",net,".jpg", sep = ""), dpi = 300, width = 8.5, height = 3)                
                                
                            } else {
                                
                                message(paste("NOT ENOUGH STATIONS", sep = ""))
                                
                            }
         
                        } else {
                            
                            message(paste("NOT ENOUGH STATIONS", sep = ""))
                            
                        }
                                               
            } # eo for loop - net types
        
} # eo fun

# Apply fun above in for loop
colnames(esd)[c(8:10,13:14,22)] <- c("ESD_Cnidaria","ESD_Tunicata","ESD_Eumalacostracat","ESD_Ostracoda+Cladocera","ESD_Copepoda (unidentified)","ESD_Calanoida (unidentified)")

variables <- colnames(esd)[c(5,7:27,29:37)] ; variables

for(var in variables) {
    distrib.analyzer(var = var, data = esd)   
}



