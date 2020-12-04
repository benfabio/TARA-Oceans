
### 01/09/2020: Script to explore an alternative representation of the DTW clusters (MDS/nMDS)
### Aims to: 
### - Perform DTW clustering to classify the resp variabels based on the shapes of their reponses
### - Extract the distance matrix resulting from the clustering (so euclid distance matrix between groups (nets)
### - Plot distances between zoo groups on a MDS 

### Last update: 02/09/2020

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
library("dtwclust")
library("parallel")

world <- map_data("world") # coastline for maps
WD <- getwd()

# --------------------------------------------------------------------------------------------------------------------------------

### Load abund and esd data
abund <- read.table("table_abund+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
esd <- read.table("table_ESD+hydro_allnets_18_06_20.txt", sep = "\t", h = T)

nets <- c("wp2","bongo","regent")
n <- "wp2"

res <- mclapply(nets, function(n) {
    
                setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",n, sep = ""))
                files <- dir()[grep("gam_",dir())]
                # Exclude those for Harpact
                f <- files[85]
                
                gams <- lapply(files, function(f) {
                            
                            message(paste("Retrieving GAM results for ",f, sep = ""))
                            model <- get(load(f))
                            
                            require("mgcv")
                            plot.data <- plot.gam(model, pages = 1)
                            # str(plot.data)
                            dev.off()
                            
                            # Extract the terms of the model from the file name
                            terms <- do.call(cbind, strsplit(as.character(f),"\\+"))
                            # 1st element contains response var and first predictor
                            terms[1,] <- str_replace_all(terms[1,], "gam_", "")
        
                            # Watchout, if Copepoda_unid is the group, then there's an extra underscore to consider
                            if( grepl("Copepoda_unid", terms[1,]) | grepl("Calanoida_unidentified", terms[1,]) | 
                                grepl("Gel_carn", terms[1,])| grepl("Gel_FF", terms[1,]) | grepl("Small_grazers", terms[1,]) ) {
                                
                                resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2],unlist(strsplit(terms[1,],"_"))[3],sep = "_")
                                pred1 <- unlist(strsplit(terms[1,],"_"))[4]
                                # And extract the last pred from the last element of terms
                                predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
                                
                            } else {
                                
                                resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2], sep = "_")
                                pred1 <- unlist(strsplit(terms[1,],"_"))[3]
                                # And extract the last pred from the last element of terms
                                predlast <- unlist(strsplit(terms[nrow(terms),]
                                ,"_"))[1]
                                
                            } # eo if else loop - grepl("other_Copepoda", terms[1,]) 
                                
                            preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],predlast)
                       
                            smooth.fits <- data.frame(resp = resp, net = n, 
                                    pred1 = plot.data[[1]]$fit, pred2 = plot.data[[2]]$fit, pred3 = plot.data[[3]]$fit,
                                    pred4 = plot.data[[4]]$fit, pred5 = plot.data[[5]]$fit, pred6 = plot.data[[6]]$fit,
                                    pred7 = plot.data[[7]]$fit, pred8 = plot.data[[8]]$fit, pred9 = plot.data[[9]]$fit
                            ) # eo ddf
                       
                            # Extract information about the GAM
                            smooth.fits$formula <- paste(preds, collapse = "+")
                            smooth.fits$Deviance <- summary(model)$dev.expl # can serve as weight for the ordination analysis
                            
                            rownames(smooth.fits) <- NULL
                            return(smooth.fits)
                            
                    } # eo FUN
                    
                ) # eo lapply
                
                # Rbind 
                table <- do.call(rbind, gams)
                # Check
                # summary(table)
                
                return(table)
                rm(gams);gc()
    
        }, mc.cores = 3
        
) # eo mclapply
# Rbind 
table <- bind_rows(res)
head(table); dim(table)
rm(res); gc()
summary(table)

# vector of resp vars you're interested in here
table$resp <- factor(table$resp)
levels(table$resp)[levels(table$resp) == "Ab_Calanoida_unidentified"] <- "Ab_Calanoida_unid"
levels(table$resp)[levels(table$resp) == "Ab_Gel_carn"] <- "Ab_Cnidaria"
levels(table$resp)[levels(table$resp) == "Ab_Gel_FF"] <- "Ab_Tunicata"
levels(table$resp)[levels(table$resp) == "Ab_Small_grazers"] <- "Ab_Cladocera+Ostracoda"
levels(table$resp)[levels(table$resp) == "Ab_Crust"] <- "Ab_Eumalacostraca"
levels(table$resp)[levels(table$resp) == "ESD_Calanoida_unidentified"] <- "ESD_Calanoida_unid"
levels(table$resp)[levels(table$resp) == "ESD_Gel_carn"] <- "ESD_Cnidaria"
levels(table$resp)[levels(table$resp) == "ESD_Gel_FF"] <- "ESD_Tunicata"
levels(table$resp)[levels(table$resp) == "ESD_Small_grazers"] <- "ESD_Cladocera+Ostracoda"
levels(table$resp)[levels(table$resp) == "ESD_Crust"] <- "ESD_Eumalacostraca"

# Plot distribution of models' R2 and explained deviance per nets
resp_esd <- unique(table$resp)[grep('ESD_',unique(table$resp))]
resp_ab <- unique(table$resp)[grep('Ab_',unique(table$resp))]


### Perform dtwclust (multivariate example) 

### A) ESD vars  ------------------------------------------------------------------------------------

vars <- resp_esd[c(1,3,4,7,8,10:12,15,17)] ; vars

### Separate in two:
oxy.models <- table[table$formula == "O2+Salinity+MLD+PAR2+NO2NO3+Chla+bbp470+Micro+Nano",]
temp.models <- table[table$formula == "Temperature+Salinity+MLD+PAR2+NO2NO3+Chla+bbp470+Micro+Nano",]
### Change their column names
colnames(oxy.models)[c(3:11)] <- c("Oxygen","Salinity","MLD","PAR","NO2NO3","Chlorophylla","bbp470","%Micro","%Nano")
colnames(temp.models)[c(3:11)] <- c("Temperature","Salinity","MLD","PAR","NO2NO3","Chlorophylla","bbp470","%Micro","%Nano")
oxy.models2 <- oxy.models[oxy.models$resp %in% vars,]
temp.models2 <- temp.models[temp.models$resp %in% vars,]    
oxy.models2$Temperature <- temp.models2$Temperature

### Add an ID: resp+net
esd.model <- oxy.models2
esd.model$ID <- factor(paste(esd.model$resp, esd.model$net, sep = "_"))

### Only account for those ESD models that show >50% explaiend deviance
mods2choose <- unique(esd.model[esd.model$Deviance >= 0.50,"ID"]) ; mods2choose
esd.model2 <- esd.model[esd.model$ID %in% mods2choose,]
# For each unique(esd.model$ID), extract the corresponding smooth terms series and concatenate in lapply
list.esd <- lapply(unique(esd.model2$ID), function(id) {
                data <- esd.model2[esd.model2$ID == id,c(14,3:11)]
                return(data)
        } # eo lapply
) # eo lapply
names(list.esd) <- unique(esd.model2$ID)

### Perform ts clustering
clusters.pam <- tsclust(list.esd, k = 4L, distance = "dtw_basic", centroid = "pam", type = "partitional")
ddf_k <- data.frame(Resp = unique(esd.model2$ID), k = clusters.pam@cluster)
ddf_k[order(ddf_k$k),]

### Try cmdscale (MDS)
mds <- cmdscale(as.dist(clusters.pam@distmat), k = 2, eig = T)
ddf <- data.frame(mds$point)
colnames(ddf) <- c("x","y")
ddf$group <- rownames(ddf)
# data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))
ddf$zoo.group <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,2]
ddf$net <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,3]
levels(ddf$net)[levels(ddf$net) == "regent"] <- "Régent"
levels(ddf$net)[levels(ddf$net) == "wp2"] <- "WP2"
levels(ddf$net)[levels(ddf$net) == "bongo"] <- "Bongo"
ddf$id <- factor(paste(ddf$zoo.group," ","(",ddf$net,")", sep = ""))

### Provide cluster from ddf_k
ddf$cluster <- NA
for(g in unique(ddf$group)) {
    k <- ddf_k[ddf_k$Resp == g,"k"]
    ddf[ddf$group == g,"cluster"] <- k
}
ddf

plot.mds <- ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_point(aes(x = x, y = y, fill = factor(cluster)), data = ddf, colour = "black", size = 3, pch = 21) + 
    geom_text_repel(aes(x = x, y = y, label = id), data = ddf) + scale_fill_brewer(name = "Cluster", palette = "Spectral") +
    xlab("MDS1") + ylab("MDS2") + theme_bw() + ggtitle("MDS - median ESD")
quartz()
plot.mds

### Try nMDS
# library("vegan")
# #?metaMDS
# nmds <- metaMDS(as.dist(clusters.pam@distmat), k = 2, try = 20, trymax = 20)
#
# ddf <- data.frame(nmds$points)
# colnames(ddf) <- c("x","y")
# ddf$group <- rownames(ddf)
# # data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))
# ddf$zoo.group <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,2]
# ddf$net <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,3]
#
# levels(ddf$net)[levels(ddf$net) == "regent"] <- "Régent"
# levels(ddf$net)[levels(ddf$net) == "wp2"] <- "WP2"
# levels(ddf$net)[levels(ddf$net) == "bongo"] <- "Bongo"
#
# ddf$cluster <- NA
# ddf[ddf$group %in% c("ESD_Calanoida_wp2","ESD_Chaetognatha_bongo","ESD_Copepoda_bongo",
#     "ESD_Eumalacostraca_bongo","ESD_Eumalacostraca_regent","ESD_Poecilostomatoida_wp2",
#     "ESD_Rhizaria_wp2","ESD_Zooplankton_bongo"),"cluster"] <- 1
# ddf[ddf$group == "ESD_Cnidaria_wp2","cluster"] <- 2
# ddf[ddf$group == "ESD_Chaetognatha_regent","cluster"] <- 3
# ddf[ddf$group %in% c("ESD_Calanoida_regent","ESD_Cnidaria_bongo","ESD_Copepoda_regent"),"cluster"] <- 4
#
# ddf$id <- factor(paste(ddf$zoo.group," ","(",ddf$net,")", sep = ""))
#
# #quartz()
# plot.nmds <- ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
#     geom_point(aes(x = x, y = y, fill = factor(cluster)), data = ddf, colour = "black", size = 3, pch = 21) +
#     geom_text_repel(aes(x = x, y = y, label = id), data = ddf) + scale_fill_brewer(name = "Cluster", palette = "Spectral") +
#     xlab("nMDS1") + ylab("nMDS2") + theme_bw() + ggtitle(paste("nMDS (stress = ",round(nmds$stress,3),") - median ESD", sep = ""))

# Save plots
ggsave(plot = plot.mds, "plot_MDS_clusters_esd.jpg", dpi = 300, width = 7, height = 7)
#ggsave(plot = plot.nmds, "plot_nMDS_clusters_esd.jpg", dpi = 300, width = 7, height = 7)


### 02/09/2020: Plot facet curves
head(esd.model2)
esd.model2$cluster_pam <- NA
for(id in unique(esd.model2$ID)) {
    k <- ddf_k[ddf_k$Resp == id,"k"]
    esd.model2[esd.model2$ID == id,"cluster_pam"] <- k
}
summary(esd.model2$cluster_pam)

### Forgot to keep the x axis ^^' (1:100)
esd.model2$x <- rep(x = c(1:100), times = 13)

# Melt to have covariates as long format
m.esd <- melt(esd.model2, id.vars = c("cluster_pam","resp","net","Deviance","formula","ID","x"))
head(m.esd)

### Change the ID labels for better plot 
m.esd$net <- factor(m.esd$net)
levels(m.esd$net)[levels(m.esd$net) == "wp2"] <- "(WP2)"
levels(m.esd$net)[levels(m.esd$net) == "regent"] <- "(Régent)"
levels(m.esd$net)[levels(m.esd$net) == "bongo"] <- "(Bongo)"

m.esd$ID2 <- factor(paste(str_replace_all(m.esd$resp,"_"," "), m.esd$net, sep = " "))
# Add cluster the resp belongs to 
m.esd$ID2 <- paste(m.esd$ID2," - ",m.esd$cluster_pam, sep = "")
unique(m.esd$ID2)

facet <- ggplot(data = m.esd) + geom_line(aes(x = x, y = value, colour = factor(ID2))) +
    xlab("") + ylab("") + theme_bw() + facet_grid(factor(variable) ~ factor(cluster_pam), scales = "free") + 
    scale_colour_discrete(name = "") + 
    geom_hline(yintercept = 0, linetype = "dashed")

quartz()
facet
    
ggsave(plot = facet, filename = "panel_smooth_clustersxcovariates_ESD.pdf", dpi = 300, height = 14, width = 10)


### 02/09/2020: Re-plot the facetting without clusters
head(m.esd)

quartz()
m.esd$ID3 <- factor(paste(str_replace_all(m.esd$resp,"ESD_",""), m.esd$net, sep = " "))

facet2 <- ggplot(data = m.esd) + geom_line(aes(x = x, y = value), colour = "black") +
    xlab("") + ylab("") + theme_bw() + facet_grid(factor(ID3) ~ factor(variable), scales = "free") + 
    scale_colour_discrete(name = "") + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey33")
    
ggsave(plot = facet2, filename = "panel_smooth_groupsxcovariates_ESD.pdf", dpi = 300, height = 25, width = 15)
    

### B) Abundances vars  ------------------------------------------------------------------------------------

vars2 <- resp_ab[c(4,8,9,12,13,16,17,24,26,28,32)] ; vars2
# Provide Temp smooths to the oxygen table 
oxy.models2 <- oxy.models[oxy.models$resp %in% vars2,]
temp.models2 <- temp.models[temp.models$resp %in% vars2,]    
oxy.models2$Temperature <- temp.models2$Temperature
ab.model <- oxy.models2
ab.model$ID <- factor(paste(ab.model$resp, ab.model$net, sep = "_"))
mods2choose <- unique(ab.model[ab.model$Deviance >= 0.40,"ID"]) ; mods2choose
ab.model2 <- ab.model[ab.model$ID %in% mods2choose,]
ab.model2 <- ab.model2 %>% distinct()
# For each unique(esd.model$ID), extract the corresponding smooth terms series and concatenate in lapply
list.ab <- lapply(unique(ab.model2$ID), function(id) {
                data <- ab.model2[ab.model2$ID == id,c(14,3:11)]
                return(data)
        } # eo lapply
) # eo lapply
names(list.ab) <- unique(ab.model2$ID)

### Perform PAM clustering
clusters.pam <- tsclust(list.ab, k = 5L, distance = "dtw_basic", centroid = "pam", type = "partitional")
ddf_k <- data.frame(Resp = unique(ab.model2$ID), k = clusters.pam@cluster)[order(data.frame(Resp = unique(ab.model2$ID), k = clusters.pam@cluster)$k),]

### Try cmdscale (MDS)
mds <- cmdscale(as.dist(clusters.pam@distmat), k = 2, eig = T)
ddf <- data.frame(mds$point)
colnames(ddf) <- c("x","y")
ddf$group <- rownames(ddf)
ddf$zoo.group <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,2]
ddf$net <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,3]
levels(ddf$net)[levels(ddf$net) == "regent"] <- "Régent"
levels(ddf$net)[levels(ddf$net) == "wp2"] <- "WP2"
levels(ddf$net)[levels(ddf$net) == "bongo"] <- "Bongo"
ddf$cluster <- NA
ddf[ddf$group %in% c("Ab_Calanoida_regent","Ab_Copepoda_regent","Ab_Zooplankton_regent"),"cluster"] <- 1
ddf[ddf$group %in% c("Ab_Chaetognatha_bongo","Ab_Pteropoda_bongo","Ab_Rhizaria_bongo","Ab_Rhizaria_wp2"),"cluster"] <- 2
ddf[ddf$group %in% c("Ab_Chaetognatha_regent","Ab_Cnidaria_wp2","Ab_Pteropoda_regent","Ab_Pteropoda_wp2","Ab_Tunicata_bongo"),"cluster"] <- 3
ddf[ddf$group %in% c("Ab_Tunicata_wp2","Ab_Cyclopoida_wp2"),"cluster"] <- 4
ddf[ddf$group %in% c("Ab_Poecilostomatoida_wp2"),"cluster"] <- 5
ddf$id <- factor(paste(ddf$zoo.group," ","(",ddf$net,")", sep = ""))

# quartz()
plot.mds <- ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_point(aes(x = x, y = y, fill = factor(cluster)), data = ddf, colour = "black", size = 3, pch = 21) + 
    geom_text_repel(aes(x = x, y = y, label = id), data = ddf) + scale_fill_brewer(name = "Cluster", palette = "Spectral") +
    xlab("MDS1") + ylab("MDS2") + theme_bw() + ggtitle("MDS - abundances")

quartz()
plot.mds

ggsave(plot = plot.mds, "plot_MDS_clusters_abv2.jpg", dpi = 300, width = 7, height = 7)


# ### Try nMDS
# nmds <- metaMDS(as.dist(clusters.pam@distmat), k = 2, try = 20, trymax = 20)
# ddf <- data.frame(nmds$points)
# colnames(ddf) <- c("x","y")
# ddf$group <- rownames(ddf)
# ddf$zoo.group <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,2]
# ddf$net <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,3]
# levels(ddf$net)[levels(ddf$net) == "regent"] <- "Régent"
# levels(ddf$net)[levels(ddf$net) == "wp2"] <- "WP2"
# levels(ddf$net)[levels(ddf$net) == "bongo"] <- "Bongo"
# ddf$cluster <- NA
# ddf[ddf$group %in% c("Ab_Calanoida_regent","Ab_Copepoda_regent","Ab_Zooplankton_regent"),"cluster"] <- 1
# ddf[ddf$group %in% c("Ab_Chaetognatha_bongo","Ab_Pteropoda_bongo","Ab_Rhizaria_bongo","Ab_Rhizaria_wp2"),"cluster"] <- 2
# ddf[ddf$group %in% c("Ab_Chaetognatha_regent","Ab_Cnidaria_wp2","Ab_Pteropoda_regent","Ab_Pteropoda_wp2","Ab_Tunicata_bongo"),"cluster"] <- 3
# ddf[ddf$group %in% c("Ab_Tunicata_wp2","Ab_Cyclopoida_wp2"),"cluster"] <- 4
# ddf[ddf$group %in% c("Ab_Poecilostomatoida_wp2"),"cluster"] <- 5
# ddf$id <- factor(paste(ddf$zoo.group," ","(",ddf$net,")", sep = ""))
#
# #quartz()
# plot.nmds <- ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
#     geom_point(aes(x = x, y = y, fill = factor(cluster)), data = ddf, colour = "black", size = 3, pch = 21) +
#     geom_text_repel(aes(x = x, y = y, label = id), data = ddf) + scale_fill_brewer(name = "Cluster", palette = "Spectral") +
#     xlab("nMDS1") + ylab("nMDS2") + theme_bw() + ggtitle(paste("nMDS (stress = ",round(nmds$stress,3),") - abundances", sep = ""))
# #
#ggsave(plot = plot.nmds, "plot_nMDS_clusters_ab.jpg", dpi = 300, width = 7, height = 7)

### Better panel plot, like with ESD
head(ab.model2)
ab.model2$cluster_pam <- NA
#intersect(unique(ab.model2$ID), unique(ddf$group))

for(id in unique(ab.model2$ID)) {
    message(paste(id, sep = ""))
    k <- ddf[ddf$group == id,"cluster"]
    ab.model2[ab.model2$ID == id,"cluster_pam"] <- k
}
summary(ab.model2$cluster_pam)

### Forgot to keep the x axis ^^' (1:100)
ab.model2$x <- rep(x = c(1:100), times = 15)

# Melt to have covariates as long format
m.ab <- melt(ab.model2, id.vars = c("cluster_pam","resp","net","Deviance","formula","ID","x"))
head(m.ab)

### Change the ID labels for better plot 
m.ab$net <- factor(m.ab$net)
levels(m.ab$net)[levels(m.ab$net) == "wp2"] <- "(WP2)"
levels(m.ab$net)[levels(m.ab$net) == "regent"] <- "(Régent)"
levels(m.ab$net)[levels(m.ab$net) == "bongo"] <- "(Bongo)"

m.ab$ID2 <- factor(paste(str_replace_all(m.ab$resp,"_"," "), m.ab$net, sep = " "))
m.ab$ID2 <- str_replace_all(m.ab$ID2,"Ab","Abundance")
m.ab$ID2 <- paste(m.ab$ID2," - ",m.ab$cluster_pam, sep = "")
unique(m.ab$ID2)

facet.ab <- ggplot(data = m.ab) + geom_line(aes(x = x, y = value, colour = factor(ID2))) +
    xlab("") + ylab("") + theme_bw() + facet_grid(factor(variable) ~ factor(cluster_pam), scales = "free") + 
    scale_colour_discrete(name = "") + 
    geom_hline(yintercept = 0, linetype = "dashed")
#
quartz()
facet.ab

ggsave(plot = facet.ab, filename = "panel_smooth_clustersxcovariates_abundV2.pdf", dpi = 300, height = 14, width = 10)

### 02/09/2020: Re-plot the facetting without clusters
head(m.ab)

m.ab$ID3 <- factor(paste(str_replace_all(m.ab$resp,"Ab_",""), m.ab$net, sep = " "))
unique(m.ab$ID3)

facet2 <- ggplot(data = m.ab) + geom_line(aes(x = x, y = value), colour = "black") +
    xlab("") + ylab("") + theme_bw() + facet_grid(factor(ID3) ~ factor(variable), scales = "free") + 
    scale_colour_discrete(name = "") + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey33")
    
ggsave(plot = facet2, filename = "panel_smooth_groupsxcovariates_Ab.pdf", dpi = 300, height = 25, width = 15)


### Need to find a way to supply actual x values 
