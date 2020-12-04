
### 02/09/2020: Script to explore an alternative representation of the DTW clusters (MDS/nMDS)
### Aims to: plot the groups' smooth curves as a function of actual predictors values (no DTW clustering)

### Last update: 03/09/2020

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
# n <- "wp2"

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
                            
                            ### Extract the env covariates data usd for training 
                            # str(model$model)
                            # head(model$model[,c(2:length(model$model))])
                            covariates <- model$model[,c(2:length(model$model))]
                            colnames(covariates) <- preds
                            # Create a new data.frame of env covariates with the same length as smooth.fits and cbind them
                            new.covariates <- data.frame( 
                                cov1 = seq(from = min(covariates[,1]), to = max(covariates[,1]), length.out = nrow(smooth.fits)),
                                cov2 = seq(from = min(covariates[,2]), to = max(covariates[,2]), length.out = nrow(smooth.fits)),
                                cov3 = seq(from = min(covariates[,3]), to = max(covariates[,3]), length.out = nrow(smooth.fits)),
                                cov4 = seq(from = min(covariates[,4]), to = max(covariates[,4]), length.out = nrow(smooth.fits)),
                                cov5 = seq(from = min(covariates[,5]), to = max(covariates[,5]), length.out = nrow(smooth.fits)),
                                cov6 = seq(from = min(covariates[,6]), to = max(covariates[,6]), length.out = nrow(smooth.fits)),
                                cov7 = seq(from = min(covariates[,7]), to = max(covariates[,7]), length.out = nrow(smooth.fits)),
                                cov8 = seq(from = min(covariates[,8]), to = max(covariates[,8]), length.out = nrow(smooth.fits)),
                                cov9 = seq(from = min(covariates[,9]), to = max(covariates[,9]), length.out = nrow(smooth.fits))
                            )
                                                        
                            smooth.fits <- cbind(smooth.fits, new.covariates)
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


### A) ESD vars  ------------------------------------------------------------------------------------

vars <- resp_esd[c(1,3,4,7,8,10:12,15,17)] ; vars

### Separate in two:
oxy.models <- table[table$formula == "O2+Salinity+MLD+PAR2+NO2NO3+Chla+bbp470+Micro+Nano",]
temp.models <- table[table$formula == "Temperature+Salinity+MLD+PAR2+NO2NO3+Chla+bbp470+Micro+Nano",]

### Change their column names
colnames(oxy.models)[c(3:11)] <- c("Oxygen","Salinity","MLD","PAR","NO2NO3","Chlorophylla","bbp470","%Micro","%Nano")
colnames(temp.models)[c(3:11)] <- c("Temperature","Salinity","MLD","PAR","NO2NO3","Chlorophylla","bbp470","%Micro","%Nano")

colnames(oxy.models)[c(14:22)] <- c("Oxygen_fit","Salinity_fit","MLD_fit","PAR_fit","NO2NO3_fit","Chlorophylla_fit","bbp470_fit","%Micro_fit","%Nano_fit")
colnames(temp.models)[c(14:22)] <- c("Temperature_fit","Salinity_fit","MLD_fit","PAR_fit","NO2NO3_fit","Chlorophylla_fit","bbp470_fit","%Micro_fit","%Nano_fit")

oxy.models2 <- oxy.models[oxy.models$resp %in% vars,]
temp.models2 <- temp.models[temp.models$resp %in% vars,]  
  
oxy.models2$Temperature <- temp.models2$Temperature  
oxy.models2$Temperature_fit <- temp.models2$Temperature_fit

### Add an ID: resp+net
esd.model <- oxy.models2
esd.model$ID <- factor(paste(esd.model$resp, esd.model$net, sep = "_"))
head(esd.model)

### Only account for those ESD models that show >50% explaiend deviance
mods2choose <- unique(esd.model[esd.model$Deviance >= 0.50,"ID"]) ; mods2choose
esd.model2 <- esd.model[esd.model$ID %in% mods2choose,]

# # For each unique(esd.model$ID), extract the corresponding smooth terms series and concatenate in lapply
# list.esd <- lapply(unique(esd.model2$ID), function(id) {
#                 data <- esd.model2[esd.model2$ID == id,c(14,3:11)]
#                 return(data)
#         } # eo lapply
# ) # eo lapply
# names(list.esd) <- unique(esd.model2$ID)

### Figure out a way to have the following table: per ID (groupxnet) and per covariate have y (smooth terms) and x (obs covariate range) columns
head(esd.model2)

m.esd <- melt(esd.model2, id.vars = c("resp","net","Deviance","formula","ID","Temperature_fit","Oxygen_fit","Salinity_fit","MLD_fit",
                                    "PAR_fit","NO2NO3_fit","Chlorophylla_fit","bbp470_fit","%Micro_fit","%Nano_fit"))
colnames(m.esd)[c(16,17)] <- c("Covariate","smooths")
#head(m.esd)
m.esd$obs <- NA
# i <- unique(m.esd$ID)[1]
for(i in unique(m.esd$ID)) {
    
    temp <- m.esd[m.esd$ID == i & m.esd$Covariate == "Temperature","Temperature_fit"]
    do2 <- m.esd[m.esd$ID == i & m.esd$Covariate == "Oxygen","Oxygen_fit"]
    sal <- m.esd[m.esd$ID == i & m.esd$Covariate == "Salinity","Salinity_fit"]
    mld <- m.esd[m.esd$ID == i & m.esd$Covariate == "MLD","MLD_fit"]
    par <- m.esd[m.esd$ID == i & m.esd$Covariate == "PAR","PAR_fit"]
    no2no3 <- m.esd[m.esd$ID == i & m.esd$Covariate == "NO2NO3","NO2NO3_fit"]
    chla <- m.esd[m.esd$ID == i & m.esd$Covariate == "Chlorophylla","Chlorophylla_fit"]
    bbp <- m.esd[m.esd$ID == i & m.esd$Covariate == "bbp470","bbp470_fit"]
    micro <- m.esd[m.esd$ID == i & m.esd$Covariate == "%Micro","%Micro_fit"]
    nano <- m.esd[m.esd$ID == i & m.esd$Covariate == "%Nano","%Nano_fit"]
    
    m.esd[m.esd$ID == i & m.esd$Covariate == "Temperature","obs"] <- temp
    m.esd[m.esd$ID == i & m.esd$Covariate == "Oxygen","obs"] <- do2
    m.esd[m.esd$ID == i & m.esd$Covariate == "Salinity","obs"] <- sal
    m.esd[m.esd$ID == i & m.esd$Covariate == "MLD","obs"] <- mld
    m.esd[m.esd$ID == i & m.esd$Covariate == "PAR","obs"] <- par
    m.esd[m.esd$ID == i & m.esd$Covariate == "NO2NO3","obs"] <- no2no3
    m.esd[m.esd$ID == i & m.esd$Covariate == "Chlorophylla","obs"] <- chla
    m.esd[m.esd$ID == i & m.esd$Covariate == "bbp470","obs"] <- bbp
    m.esd[m.esd$ID == i & m.esd$Covariate == "%Micro","obs"] <- micro
    m.esd[m.esd$ID == i & m.esd$Covariate == "%Nano","obs"] <- nano
    
} # eo for loop 

head(m.esd)

# Change the ID labels for better plot 
m.esd$net <- factor(m.esd$net)
levels(m.esd$net)[levels(m.esd$net) == "wp2"] <- "(WP2)"
levels(m.esd$net)[levels(m.esd$net) == "regent"] <- "(Régent)"
levels(m.esd$net)[levels(m.esd$net) == "bongo"] <- "(Bongo)"
m.esd$ID3 <- factor(paste(str_replace_all(m.esd$resp,"ESD_",""), m.esd$net, sep = " "))
#head(m.esd)

facet2 <- ggplot(data = m.esd) + geom_line(aes(x = obs, y = smooths), colour = "black") +
    xlab("") + ylab("") + theme_bw() + facet_grid(factor(ID3) ~ factor(Covariate), scales = "free") + 
    scale_colour_discrete(name = "") + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey33")
    
ggsave(plot = facet2, filename = "panel_smooth_groupsxcovariates_ESD.pdf", dpi = 300, height = 23, width = 17)
   
   


### B) Abundance vars  ------------------------------------------------------------------------------------

vars2 <- resp_ab[c(4,8,9,12,13,16,17,24,26,28,32)] ; vars2

oxy.models <- table[table$formula == "O2+Salinity+MLD+PAR2+NO2NO3+Chla+bbp470+Micro+Nano",]
temp.models <- table[table$formula == "Temperature+Salinity+MLD+PAR2+NO2NO3+Chla+bbp470+Micro+Nano",]
colnames(oxy.models)[c(3:11)] <- c("Oxygen","Salinity","MLD","PAR","NO2NO3","Chlorophylla","bbp470","%Micro","%Nano")
colnames(temp.models)[c(3:11)] <- c("Temperature","Salinity","MLD","PAR","NO2NO3","Chlorophylla","bbp470","%Micro","%Nano")
colnames(oxy.models)[c(14:22)] <- c("Oxygen_fit","Salinity_fit","MLD_fit","PAR_fit","NO2NO3_fit","Chlorophylla_fit","bbp470_fit","%Micro_fit","%Nano_fit")
colnames(temp.models)[c(14:22)] <- c("Temperature_fit","Salinity_fit","MLD_fit","PAR_fit","NO2NO3_fit","Chlorophylla_fit","bbp470_fit","%Micro_fit","%Nano_fit")
oxy.models2 <- oxy.models[oxy.models$resp %in% vars2,]
temp.models2 <- temp.models[temp.models$resp %in% vars2,]  
oxy.models2$Temperature <- temp.models2$Temperature  
oxy.models2$Temperature_fit <- temp.models2$Temperature_fit
ab.model <- oxy.models2
ab.model$ID <- factor(paste(ab.model$resp, ab.model$net, sep = "_"))
mods2choose <- unique(ab.model[ab.model$Deviance >= 0.40,"ID"]) ; mods2choose
ab.model2 <- ab.model[ab.model$ID %in% mods2choose,]
ab.model2 <- ab.model2 %>% distinct()

# Melt to have covariates as long format
m.ab <- melt(ab.model2, id.vars = c("resp","net","Deviance","formula","ID","Temperature_fit","Oxygen_fit","Salinity_fit","MLD_fit",
                                    "PAR_fit","NO2NO3_fit","Chlorophylla_fit","bbp470_fit","%Micro_fit","%Nano_fit"))
colnames(m.ab)[c(16,17)] <- c("Covariate","smooths")
#head(m.ab)
m.ab$obs <- NA

for(i in unique(m.ab$ID)) {
    
    temp <- m.ab[m.ab$ID == i & m.ab$Covariate == "Temperature","Temperature_fit"]
    do2 <- m.ab[m.ab$ID == i & m.ab$Covariate == "Oxygen","Oxygen_fit"]
    sal <- m.ab[m.ab$ID == i & m.ab$Covariate == "Salinity","Salinity_fit"]
    mld <- m.ab[m.ab$ID == i & m.ab$Covariate == "MLD","MLD_fit"]
    par <- m.ab[m.ab$ID == i & m.ab$Covariate == "PAR","PAR_fit"]
    no2no3 <- m.ab[m.ab$ID == i & m.ab$Covariate == "NO2NO3","NO2NO3_fit"]
    chla <- m.ab[m.ab$ID == i & m.ab$Covariate == "Chlorophylla","Chlorophylla_fit"]
    bbp <- m.ab[m.ab$ID == i & m.ab$Covariate == "bbp470","bbp470_fit"]
    micro <- m.ab[m.ab$ID == i & m.ab$Covariate == "%Micro","%Micro_fit"]
    nano <- m.ab[m.ab$ID == i & m.ab$Covariate == "%Nano","%Nano_fit"]
    
    m.ab[m.ab$ID == i & m.ab$Covariate == "Temperature","obs"] <- temp
    m.ab[m.ab$ID == i & m.ab$Covariate == "Oxygen","obs"] <- do2
    m.ab[m.ab$ID == i & m.ab$Covariate == "Salinity","obs"] <- sal
    m.ab[m.ab$ID == i & m.ab$Covariate == "MLD","obs"] <- mld
    m.ab[m.ab$ID == i & m.ab$Covariate == "PAR","obs"] <- par
    m.ab[m.ab$ID == i & m.ab$Covariate == "NO2NO3","obs"] <- no2no3
    m.ab[m.ab$ID == i & m.ab$Covariate == "Chlorophylla","obs"] <- chla
    m.ab[m.ab$ID == i & m.ab$Covariate == "bbp470","obs"] <- bbp
    m.ab[m.ab$ID == i & m.ab$Covariate == "%Micro","obs"] <- micro
    m.ab[m.ab$ID == i & m.ab$Covariate == "%Nano","obs"] <- nano
    
} # eo for loop 

# Change the ID labels for better plot 
m.ab$net <- factor(m.ab$net)
levels(m.ab$net)[levels(m.ab$net) == "wp2"] <- "(WP2)"
levels(m.ab$net)[levels(m.ab$net) == "regent"] <- "(Régent)"
levels(m.ab$net)[levels(m.ab$net) == "bongo"] <- "(Bongo)"
m.ab$ID3 <- factor(paste(str_replace_all(m.ab$resp,"Ab_",""), m.ab$net, sep = " "))
#head(m.ab)

facet2 <- ggplot(data = m.ab) + geom_line(aes(x = obs, y = smooths), colour = "black") +
    xlab("") + ylab("") + theme_bw() + facet_grid(factor(ID3) ~ factor(Covariate), scales = "free") + 
    scale_colour_discrete(name = "") + 
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey33")
    
ggsave(plot = facet2, filename = "panel_smooth_groupsxcovariates_Ab.pdf", dpi = 300, height = 25, width = 17)
   
   
