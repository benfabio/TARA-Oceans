
### 29/06/2020: Script to analyze the updated zooplankton imaging data from TARA Oceans (WP2, Bongo & Regent).
### Aims to: 
### - For each net, choose the set of env predictors (uncorrelated and transformed) to be used for training the GAMs
### - Train GAMs to predict ESD and abundacnes for all mesozooplankton groups based on the various sets of predictors
### - Explore the GAMs results (signif of smooth terms, r2, resp curves)
### Choose emerging results and summarize in a Table/ Figure
### - Perform DTW clustering to classify the resp variabels based on the shapes of their reponses

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

world <- map_data("world") # coastline for maps
WD <- getwd()

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) First, analyze the correlations heatmaps to define the pools of predictors
setwd(WD)
abund <- read.table("table_abund+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
esd <- read.table("table_ESD+hydro_allnets_18_06_20.txt", sep = "\t", h = T)

#summary(abund[abund$net == "wp2",])
# - Temperature or O2
# - Salinity
# - MLD
# - PAR2
# - Micro + Nano (exclude Pico)
# - NO2NO3 (exclude SiO2 and PO4)
# - Chl-a
# - bbp470

### So 2 sets of predictors: one with SST and the other with O2
colnames(abund)
list.preds <- list(
        c("Temperature","Salinity","MLD","PAR2","NO2NO3","Chla","bbp470","Micro","Nano"), 
        c("O2","Salinity","MLD","PAR2","NO2NO3","Chla","bbp470","Micro","Nano")
) # eo list.pred.bongo

### And define vector of resp vars
responses <- colnames(abund)[c(5:27,29:37)] ; responses
# net <- "wp2"
# variables <- list.preds[[1]]
# v <- "Ab_Rhizaria"
gam.fitter <- function(net, variables) {
    
                setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net, sep = ""))
                message(paste("\nTraining GAMs for net = ",net, sep = ""))
                subset <- abund[abund$net == net,]
                
                # Run GAMs for each element of this list
                require("parallel")
                # i <- 1
                mclapply(c(1:length(list.preds)), function(i) {
                    
                        preds <- list.preds[[i]]
                        formula <- paste(preds, collapse = "+")
                        npred <- 9 
                        
                        # Train GAM model for each variable
                        for(v in variables) {
                                
                                subset2 <- na.omit(subset[,c(v,preds)])
                                nobs <- nrow(subset2)
                                
                                # Control for a minimum of 15 obs, or if the resp var is a flat one (abundances always = 0)
                                if( nobs >= 25 & sum(subset2[,v]) > 0 ) {
                                    
                                    # The number of coeff in the model should lower than the nb of obs -> nobs*9
                                    message(paste("Training GAM with list ",i,"  | ",v, " ~ ",formula, sep = ""))
                                    k <- round(nobs/npred, digits = 0)
                                    model <- mgcv::gam(get(v) ~ s(get(preds[1]),bs="tp",k=k)+s(get(preds[2]),bs="tp",k=k)+s(get(preds[3]),bs="tp",k=k)+
                                                s(get(preds[4]),bs="tp",k=k)+s(get(preds[5]),bs="tp",k=k)+s(get(preds[6]),bs="tp",k=k)+
                                                s(get(preds[7]),bs="tp",k=k)+s(get(preds[8]),bs="tp",k=k)+s(get(preds[9]),bs="tp",k=k), 
                                                data = subset2, select = T, method = "REML")
                                            
                                    save(model, file = paste("gam_",v,"_",formula,"_",net,".Rdata", sep = "") )      
                                    
                                } else {
                                    
                                    message(paste("Stopped because not enough data || n = ",nobs, sep = ""))
                                    
                                } # eo if else loop - for controlling nobs     
                                
                        } # eo for v in variables
                    
                    }, mc.cores = 2
                
                ) # eo mclapply
    
} # eo gam.fitter FUN
# Apply to each 3 net
gam.fitter(net = "wp2", variables = responses)
gam.fitter(net = "bongo", variables = responses)
gam.fitter(net = "regent", variables = responses)



### Cool, now same as above but for ESD ! Same but with a twist: only use those stations with Nind > 20
setwd(WD)
responses <- colnames(esd)[c(5:27,29:37)] ; responses
# Get counts per station
counts.bongo <- read.table("table_counts_groups_bongo.txt", h = T, sep = "\t")
counts.wp2 <- read.table("table_counts_groups_wp2.txt", h = T, sep = "\t")
counts.regent <- read.table("table_counts_groups_regent.txt", h = T, sep = "\t")
# dim(counts.bongo) ; dim(counts.wp2) ; dim(counts.regent)
# head(counts.bongo)
# net <- "wp2"
# v <- "ESD_Rhizaria"
# i <- 1
gam.fitter2 <- function(net, variables) {
    
                setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net, sep = ""))
                message(paste("\nTraining GAMs for net = ",net, sep = ""))
                subset <- esd[esd$net == net,]
                
                # Select appropriate net counts data tabke
                if(net == "bongo") {
                    counts <- counts.bongo
                } else if(net == "wp2") {
                    counts <- counts.wp2
                } else {
                    counts <- counts.regent
                } # eo if else loop
                # Change some levels in counts so they match the colnames of abund/esd
                levels(counts$group)[levels(counts$group) == "Calanoida_unid"] <- "Calanoida_unidentified"
                levels(counts$group)[levels(counts$group) == "Crustaceans"] <- "Crust"
                levels(counts$group)[levels(counts$group) == "Grazers"] <- "Small_grazers"
                
                # Run GAMs for each element of this list
                require("parallel")
                # i <- 1
                mclapply(c(1:length(list.preds)), function(i) {
                    
                        preds <- list.preds[[i]]
                        formula <- paste(preds, collapse = "+")
                        npred <- 9 
                        
                        # Train GAM model for each variable
                        for(v in variables) {
                                
                                subset2 <- na.omit(subset[,c("Station",v,preds)])
                                ### And reduce based on counts for the group
                                group <- str_replace_all(v,"ESD_","")
                                
                                # Check if group is present in counts data
                                if( group %in% unique(counts$group) | group == "Zooplankton" ) {
                                    
                                    if( group == "Zooplankton" ) {
                                        
                                        subset3 <- subset2[,c(v,preds)]
                                        nobs <- nrow(subset3)
                                
                                        # Control for a minimum of 15 obs, or if the resp var is a flat one (abundances always = 0)
                                        if( nobs >= 20 & sum(subset3[,v]) > 0 ) {
                                    
                                            # The number of coeff in the model should lower than the nb of obs -> nobs*9
                                            message(paste("Training GAM with list ",i,"  | ",v, " ~ ",formula, sep = ""))
                                            k <- round(nobs/npred, digits = 0)
                                            model <- mgcv::gam(get(v) ~ s(get(preds[1]),bs="tp",k=k)+s(get(preds[2]),bs="tp",k=k)+s(get(preds[3]),bs="tp",k=k)+
                                                    s(get(preds[4]),bs="tp",k=k)+s(get(preds[5]),bs="tp",k=k)+s(get(preds[6]),bs="tp",k=k)+
                                                    s(get(preds[7]),bs="tp",k=k)+s(get(preds[8]),bs="tp",k=k)+s(get(preds[9]),bs="tp",k=k), 
                                                    data = subset3, select = T, method = "REML")
                                            
                                            save(model, file = paste("gam_",v,"_",formula,"_",net,".Rdata", sep = "") )      
                                            
                                        } else {
                                    
                                             message(paste("Stopped because ", group, " is not in counts table", sep = ""))
                                    
                                        } # eo if else loop - for controlling nobs     
                                        
                                    } else {
                                        
                                        stations2keep <- counts[counts$group == group & counts$n >= 25 ,"station"]
                                        subset3 <- subset2[subset2$Station %in% stations2keep, c(v,preds)]
                                        nobs <- nrow(subset3)
                                
                                        # Control for a minimum of 15 obs, or if the resp var is a flat one (abundances always = 0)
                                        if( nobs >= 20 & sum(subset3[,v]) > 0 ) {
                                    
                                            # The number of coeff in the model should lower than the nb of obs -> nobs*9
                                            message(paste("Training GAM with list ",i,"  | ",v, " ~ ",formula, sep = ""))
                                            k <- round(nobs/npred, digits = 0)
                                            model <- mgcv::gam(get(v) ~ s(get(preds[1]),bs="tp",k=k)+s(get(preds[2]),bs="tp",k=k)+s(get(preds[3]),bs="tp",k=k)+
                                                    s(get(preds[4]),bs="tp",k=k)+s(get(preds[5]),bs="tp",k=k)+s(get(preds[6]),bs="tp",k=k)+
                                                    s(get(preds[7]),bs="tp",k=k)+s(get(preds[8]),bs="tp",k=k)+s(get(preds[9]),bs="tp",k=k), 
                                                    data = subset3, select = T, method = "REML")
                                            
                                            save(model, file = paste("gam_",v,"_",formula,"_",net,".Rdata", sep = "") )      
                                    
                                        } else {
                                    
                                            message(paste("Stopped because not enough data || n = ",nobs, sep = ""))
                                    
                                        } # eo if else loop - for controlling nobs  
                
                                    } # eo if else loop - for when group == "Zooplankton" 
                                                        
                                } else {
                                    
                                     message(paste("Stopped because ", group, " is not in counts table", sep = ""))
                                    
                                }   # eo if else loop - for controlling if group is present in the counts data.table  
                                
                        } # eo for v in variables
                    
                    }, mc.cores = 2
                
                ) # eo mclapply
    
} # eo gam.fitter FUN

gam.fitter2(net = "wp2", variables = responses)
gam.fitter2(net = "bongo", variables = responses)
gam.fitter2(net = "regent", variables = responses)

# --------------------------------------------------------------------------------------------------------------------------------

### And now, examine:
# - for which reponse var, which is the best model combination
# - which are the best explaining variables (importance in RF, p-values in GAMs)
# - then, plot response curves

### For each net, retrieve the results of the GAMs
nets <- c("wp2","bongo","regent")
n <- "bongo"
require("parallel")
res <- mclapply(nets, function(n) {
    
                setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",n, sep = ""))
                files <- dir()[grep("gam_",dir())]
                # Exclude those for Harpact
                # f <- files[85]
                gams <- lapply(files, function(f) {
                            
                            message(paste("Retrieving GAM results for ",f, sep = ""))
                            model <- get(load(f))
                            # Extract the terms of the model from the file name
                            terms <- do.call(cbind, strsplit(as.character(f),"\\+"))
                            # 1st element contains response var and first predictor
                            terms[1,] <- str_replace_all(terms[1,], "gam_", "")
                            
                            # plot.data <- plot(model)
                            plot.data[[2]]$fit
                            
                            # unlist(strsplit(terms[1,],"_"))
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
                                predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
                                
                            } # eo if else loop - grepl("other_Copepoda", terms[1,]) 
                                
                            preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],predlast)
                       
                            # Extract r2 and AIC of the GAM model
                            aic <- round(model$aic, 3)
                            r2 <- round(summary(model)$r.sq, 3)
                            formula <- paste(preds, collapse = "+")
                            pvalues <- summary(model)$s.pv
                            Fstats <- summary(model)$s.table[,3]
                            dev <- round(summary(model)$dev.expl, 3)
                            nobs <- summary(model)$n
                            
                            table <- data.frame(resp = resp, preds = preds, net = n, formula = formula, Nobs = nobs, Dev = dev, R2 = r2, AIC = aic, pval = pvalues, F = Fstats)
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

table$resp <- factor(table$resp)
# levels(table$resp)
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

dim(table[which(table$resp %in% resp_esd & table$net == "wp2"),])
dim(table[which(table$resp %in% resp_esd & table$net == "bongo"),])
dim(table[which(table$resp %in% resp_esd & table$net == "regent"),])


# dim(table[which(table$resp %in% resp_ab),])

# quartz()
# ggplot(aes(x = factor(net),y = R2), data = table[which(table$resp %in% resp_esd),]) + geom_violin(aes(fill = factor(net))) +
#     geom_boxplot(fill = "white", colour = "black", width = 0.2) + xlab("") + ylab("GAMs Adjusted R2") +
#     scale_fill_discrete(name = "") + theme_classic()

median(table[which(table$resp %in% resp_esd),"Dev"])
IQR(table[which(table$resp %in% resp_esd),"Dev"])
quantile(table[which(table$resp %in% resp_esd),"Dev"])

quartz()
ggplot(aes(x = factor(net),y = Dev), data = table[which(table$resp %in% resp_esd),]) + geom_violin(aes(fill = factor(net))) +
    geom_boxplot(fill = "white", colour = "black", width = 0.2) + xlab("") + ylab("GAMS explained deviance (%)") +
    scale_fill_discrete(name = "") + theme_classic()
#    
quartz()
ggplot(aes(x = factor(net),y = Dev), data = table[which(table$resp %in% resp_ab),]) + geom_violin(aes(fill = factor(net))) +
    geom_boxplot(fill = "white", colour = "black", width = 0.2) + xlab("") + ylab("GAMS explained deviance (%)") +
    scale_fill_discrete(name = "") + theme_classic()

### --> Régent GAMs clearly better than WP2 and then Bongo to explain size structure
### --> no difference regarding abundance pattern ?
### Associated tests
kruskal.test(data = table[which(table$resp %in% resp_esd),], Dev ~ factor(net) )
# Post hoc tests
library("PMCMR")
posthoc.kruskal.dunn.test(Dev ~ factor(net), data = table[which(table$resp %in% resp_esd),], p.adjust = "bonf")

### And test %Dev between models including O2 vs Temperature
kruskal.test(data = table[which(table$resp %in% resp_esd & table$preds %in% c("O2","Temperature")),], Dev ~ factor(preds) )
# ANd per net? 
kruskal.test(data = table[which(table$resp %in% resp_esd & table$preds %in% c("O2","Temperature") & table$net == "wp2"),], Dev ~ factor(preds) )
kruskal.test(data = table[which(table$resp %in% resp_esd & table$preds %in% c("O2","Temperature") & table$net == "bongo"),], Dev ~ factor(preds) )
kruskal.test(data = table[which(table$resp %in% resp_esd & table$preds %in% c("O2","Temperature") & table$net == "regent"),], Dev ~ factor(preds) )

### Check distribution of explained deviance of ESD and Abund 
dev.esd <- data.frame(table[which(table$resp %in% resp_esd),] %>% group_by(net) %>% summarize(med.R2 = median(R2), IQR.R2 = IQR(R2), med.dev = median(Dev), IQR.dev = IQR(Dev)) )
dev.abund <- data.frame(table[which(table$resp %in% resp_ab),] %>% group_by(net) %>% summarize(med.R2 = median(R2), IQR.R2 = IQR(R2), med.dev = median(Dev), IQR.dev = IQR(Dev)) )
dev.esd
dev.abund


### Examine, for each net, the resp vars that are best modelled
best.models <- data.frame(table %>% group_by(net,resp,formula) %>% summarize(r2 = unique(R2), Dev = unique(Dev)) )
# Show best models according to net
best.models[best.models$n == "wp2" & best.models$resp %in% resp_esd,c("resp","Dev")][order(best.models[best.models$n == "wp2" & best.models$resp %in% resp_esd,"Dev"], decreasing = T),]
best.models[best.models$n == "bongo" & best.models$resp %in% resp_esd,c("resp","Dev")][order(best.models[best.models$n == "bongo" & best.models$resp %in% resp_esd,"Dev"], decreasing = T),]
best.models[best.models$n == "regent" & best.models$resp %in% resp_esd,c("resp","Dev")][order(best.models[best.models$n == "regent" & best.models$resp %in% resp_esd,"Dev"], decreasing = T),]

dev.esd <- data.frame(table[which(table$resp %in% resp_esd),] %>% group_by(formula) %>% summarize(med.R2 = median(R2), IQR.R2 = IQR(R2), med.dev = median(Dev), IQR.dev = IQR(Dev)) )
dev.esd

### Change some labels first
levels(table$preds)[levels(table$preds) == "Chla"] <- "Chlorophyll a"
levels(table$preds)[levels(table$preds) == "Micro"] <- "%Micro"
levels(table$preds)[levels(table$preds) == "Nano"] <- "%Nano"
levels(table$preds)[levels(table$preds) == "O2"] <- "Oxygen"
levels(table$preds)[levels(table$preds) == "PAR2"] <- "PAR"
# And net labels
table$net <- factor(table$net)
levels(table$net)[levels(table$net) == "regent"] <- "Régent"
levels(table$net)[levels(table$net) == "bongo"] <- "Bongo"
levels(table$net)[levels(table$net) == "wp2"] <- "WP2"


### Ok, now you want to examine the ranks of the covariates (not only those that are significant)
plot1 <- ggplot(aes(x = reorder(preds, rank, median), y = rank, fill = factor(preds)), data = table[which(table$resp %in% resp_esd),]) + 
    geom_boxplot(colour = "black", notch = F) + scale_fill_brewer(name = "", palette = "Spectral") + 
    xlab("Covariates") + ylab("Normalized rank") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1) )
    
plot2 <- ggplot(aes(x = reorder(preds, rank, median), y = rank, fill = factor(preds)), data = table[which(table$resp %in% resp_esd),]) + 
    geom_boxplot(colour = "black", notch = F) + scale_fill_brewer(name = "", palette = "Spectral") + 
    xlab("Covariates") + ylab("Normalized rank") + theme_bw() + facet_wrap(.~ factor(net)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1) )
# Save
ggsave(plot = plot1, filename = "plot_ranks_allnets.pdf", dpi = 300, height = 4, width = 9)
ggsave(plot = plot2, filename = "plot_ranks_per_net.pdf", dpi = 300, height = 4, width = 9)

### Display the ranks of the signif smooth terms per response var and net
table[which(table$resp == "ESD_Poecilostomatoida" & table$pval < 0.05 & table$net == "WP2"),c("formula","preds","pval","rank")]

### And when filtering out the groups not shown in Table 1
ggplot(aes(x = reorder(preds, rank, median), y = rank, fill = factor(preds)),
        data = table[which(table$resp %in% resp_esd[c(1,3,4,7,8,10:12,15,17)]),]) + 
    geom_boxplot(colour = "black", notch = F) + scale_fill_brewer(name = "", palette = "Spectral") + 
    xlab("Covariates") + ylab("Normalized rank") + theme_bw() + facet_wrap(.~ factor(net)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1) )


### Test if abundance GAMs display more explanatory power than the median ESD ones
table$type <- NA
table[table$resp %in% resp_esd,"type"] <- "Size"
table[table$resp %in% resp_ab,"type"] <- "Abundance"
table$type <- factor(table$type)

median(table[table$type == "Size","Dev"]) ; IQR(table[table$type == "Size","Dev"])
median(table[table$type == "Abundance","Dev"]) ; IQR(table[table$type == "Abundance","Dev"])
# --> ESD models show higher Dev
kruskal.test(data = table, Dev ~ factor(type) )

# Same across all nets?
median(table[table$type == "Size" & table$net == "WP2","Dev"]) ; median(table[table$type == "Abundance" & table$net == "WP2","Dev"])
median(table[table$type == "Size" & table$net == "Bongo","Dev"]) ; median(table[table$type == "Abundance" & table$net == "Bongo","Dev"])
median(table[table$type == "Size" & table$net == "Régent","Dev"]) ; median(table[table$type == "Abundance" & table$net == "Régent","Dev"])

kruskal.test(data = table[table$net == "WP2",], Dev ~ factor(type) )
kruskal.test(data = table[table$net == "Bongo",], Dev ~ factor(type) )
kruskal.test(data = table[table$net == "Régent",], Dev ~ factor(type) )
# All signif
posthoc.kruskal.dunn.test(Dev ~ factor(type), data = table, p.adjust = "bonf")


### 21/07/2020: Describe Abundance-based GAMs
# Display the ranks of the signif smooth terms per response var and net
tbl <- table[which(table$resp == "Ab_Cyclopoida" & table$net == "Régent" & table$pval < 0.05),c("formula","preds","pval","rank")]
na.omit(tbl[tbl$formula == unique(tbl$formula)[1],][order(tbl$rank, decreasing = TRUE),]) # for O2 models
na.omit(tbl[tbl$formula == unique(tbl$formula)[2],][order(tbl$rank, decreasing = TRUE),]) # for Temp models

# Test variations of %Dev across nets for abundance-based GAMs
library("PMCMR")
kruskal.test(data = table[which(table$resp %in% resp_ab),], Dev ~ factor(net) )
posthoc.kruskal.dunn.test(Dev ~ factor(net), data = table[which(table$resp %in% resp_ab),], p.adjust = "bonf")
dev.abund <- data.frame(table[which(table$resp %in% resp_ab),] %>% group_by(net) %>% summarize(med.dev = median(Dev), IQR.dev = IQR(Dev)) )
dev.abund

### And test %Dev between models including O2 vs Temperature
kruskal.test(data = table[which(table$resp %in% resp_ab & table$preds %in% c("Oxygen","Temperature")),], Dev ~ factor(preds) )

### Ok, now you want to examine the ranks of the covariates (not only those that are significant)
plot1 <- ggplot(aes(x = reorder(preds, rank, median), y = rank, fill = factor(preds)), data = table[which(table$resp %in% resp_ab),]) + 
    geom_boxplot(colour = "black", notch = F) + scale_fill_brewer(name = "", palette = "Spectral") + 
    xlab("Covariates") + ylab("Normalized rank") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1) )
    
plot2 <- ggplot(aes(x = reorder(preds, rank, median), y = rank, fill = factor(preds)), data = table[which(table$resp %in% resp_ab),]) + 
    geom_boxplot(colour = "black", notch = F) + scale_fill_brewer(name = "", palette = "Spectral") + 
    xlab("Covariates") + ylab("Normalized rank") + theme_bw() + facet_wrap(.~ factor(net)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1) )
# Save
ggsave(plot = plot1, filename = "plot_ranks_allnets_abund.pdf", dpi = 300, height = 4, width = 9)
ggsave(plot = plot2, filename = "plot_ranks_per_net_abund.pdf", dpi = 300, height = 4, width = 9)

### And when filtering out the groups not shown in Table 1/2
plot3 <- ggplot(aes(x = reorder(preds, rank, median), y = rank, fill = factor(preds)),
        data = table[which(table$resp %in% resp_esd[c(1,3,4,7,8,10:12,15,17)]),]) + 
    geom_boxplot(colour = "black", notch = F) + scale_fill_brewer(name = "", palette = "Spectral") + 
    xlab("Covariates") + ylab("Normalized rank") + theme_bw() + facet_wrap(.~ factor(net)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1) )

plot4 <- ggplot(aes(x = reorder(preds, rank, median), y = rank, fill = factor(preds)),
        data = table[which(table$resp %in% resp_ab[c(4,8,9,12,13,16,17,24,26,28,32)]),]) + 
    geom_boxplot(colour = "black", notch = F) + scale_fill_brewer(name = "", palette = "Spectral") + 
    xlab("Covariates") + ylab("Normalized rank") + theme_bw() + facet_wrap(.~ factor(net)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1) )
# Save
ggsave(plot = plot3, filename = "plot_ranks_allnets_esd_v2.pdf", dpi = 300, height = 4, width = 9)
ggsave(plot = plot4, filename = "plot_ranks_per_net_abund_v2.pdf", dpi = 300, height = 4, width = 9)


### In a for loop, display the outputs for each resp
# responses <- unique(best.models$resp)[grep("Ab_",unique(best.models$resp))] ; responses
# net <- "regent"
# r <- "Ab_Tunicata"
#
# for(r in responses[c(4,7,8,11,12,15,16,23:27)] ) {
#
#         message(paste("Displaying best models for ",r," for ",net," net -------------------------------------", sep = ""))
#         sub <- best.models[best.models$resp == r & best.models$net == net,]
#         message("   ")
#         # And then display the ranking of vars n table
#         print( sub )
#
#         if( nrow(sub) > 0 ) {
#
#                 for(i in c(1:nrow(sub)) ) {
#                     form <- sub[i,"formula"]
#                     message(paste("Displaying predictors ranking for ",r," modelled with ",form, sep = ""))
#                     sub.table <- table[which(table$formula == form & table$resp == r & table$net == net),]
#                     print( sub.table[order(sub.table$F, decreasing = T),c("preds","pval","F")] )
#                 } # eo for loop
#
#         } # eo if loop
#
#         message("   ")
#         message("   ")
#
# } # eo for loop - r in resp

### Summarize model's stats in a table for the SI:
dim(table); head(table)
unique(table$resp)
# perform a dcast to create a table summarizing the skillz of each model for each response variable:
# Resp/ Net/ Nobs/ Dev%/ Adjusted r2/ AIC/ p-val of each preds

casted <- dcast(data = table[,c(1:9)], formula = resp+formula+net+Nobs+Dev+R2+AIC ~ preds, value.var = "pval", fun.aggregate = mean)
dim(casted)

### Save 'casted' object and prepare it on excle as a proper file for SI
write.table(casted, file = "table_summary_GAMs.csv", sep = ";")

### OK, next step: plot response curves of the GAMs

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
library("marmap")
library("mgcv")
library("nlme")
library("MuMIn")
library("ggrepel")

world <- map_data("world") # coastline for maps
WD <- getwd()

abund <- read.table("table_abund+hydro_allnets_18_06_20.txt", sep = "\t", h = T)
esd <- read.table("table_ESD+hydro_allnets_18_06_20.txt", sep = "\t", h = T)

counts.bongo <- read.table("table_counts_groups_bongo.txt", h = T, sep = "\t")
counts.wp2 <- read.table("table_counts_groups_wp2.txt", h = T, sep = "\t")
counts.regent <- read.table("table_counts_groups_regent.txt", h = T, sep = "\t")

### For each net separately, go to the corresponding dir, list all ESD GAMs, extract, and plot rcurves
nets <- c("bongo","wp2","regent")
#net <- "wp2"

for(net in nets) {
    
        message(paste("Loading GAM models for net == ",net, sep = ""))
        message(paste("", sep = ""))
        message(paste("", sep = ""))
    
        setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net, sep = ""))
        files.gams <- dir()[grep("gam_",dir())]
        files.gams2 <- files.gams[grep("ESD_",files.gams)] # files.gams2
    
        # Select appropriate net counts data table
        if(net == "bongo") {
                counts <- counts.bongo
        } else if(net == "wp2") {
                counts <- counts.wp2
        } else {
                counts <- counts.regent
        } # eo if else loop
        # Change some levels in counts so they match the colnames of abund/esd
        levels(counts$group)[levels(counts$group) == "Calanoida_unid"] <- "Calanoida_unidentified"
        levels(counts$group)[levels(counts$group) == "Crustaceans"] <- "Crust"
        levels(counts$group)[levels(counts$group) == "Grazers"] <- "Small_grazers"    
    
        # And now, for each of these models, plot rcurves 
        # m <- files.gams2[1]
        for(m in files.gams2) {
        
             message(paste("Loading model : ",m, sep = ""))
             setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net, sep = ""))
             model <- get(load(m))
         
             # Extract the terms of the model from the filename
             terms <- do.call(cbind, strsplit(as.character(m),"\\+"))
             terms[1,] <- str_replace_all(terms[1,], "gam_", "")
        
             # Watchout, if Copepoda_unid is the group, then there's an extra underscore to consider
             if( grepl("Copepoda_unid", terms[1,]) | grepl("Calanoida_unidentified", terms[1,]) | 
                 grepl("Gel_carn", terms[1,])| grepl("Gel_FF", terms[1,]) | grepl("Small_grazers", terms[1,]) ) {
             
                     resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2],unlist(strsplit(terms[1,],"_"))[3],sep = "_")
                     pred1 <- unlist(strsplit(terms[1,],"_"))[4]
                     predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
                 
             } else {
             
                     resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2], sep = "_")
                     pred1 <- unlist(strsplit(terms[1,],"_"))[3]
                     predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
             
             } # eo if else loop - grepl("other_Copepoda", terms[1,])   
         
             preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],predlast)      
        
             # Subset 'esd' ddf to draw resp curves
             subset <- na.omit(esd[esd$net == net,c("Station",resp,preds)])
             # And reduce based on counts for the group
             group <- str_replace_all(resp,"ESD_","")
         
             # Check if group is present in counts data
             if( group %in% unique(counts$group) | group == "Zooplankton" ) {
             
                 if( group == "Zooplankton" ) {
                 
                     subset2 <- subset[,c(resp,preds)]
                     nobs <- nrow(subset2)
         
                     # Control for a minimum of 15 obs, or if the resp var is a flat one (abundances always = 0)
                     if( nobs >= 20 & sum(subset2[,resp]) > 0 ) {
             
                         k <- round(nrow(subset2)/10, digits = 0)
                     
                         if( colnames(subset2)[2] == "O2" ) {
                         
                             # Train the modle again but by changing the names of the terms to draw the resp curves
                             colnames(subset2)[c(1:10)] <- c("y","Oxygen","Salinity","MLD","PAR","NO2NO3","Chlorophylla",
                                                             "Backscattering","Microphytoplankton","Nanophytoplankton")
                             
                             model.new <- mgcv::gam(formula = y ~ s(Oxygen,bs="tp",k=k) + s(Salinity,bs="tp",k=k) + s(MLD,bs="tp",k=k) + 
                                        s(PAR,bs="tp",k=k) + s(NO2NO3,bs="tp",k=k) + s(Chlorophylla,bs="tp",k=k) + s(Backscattering,bs="tp",k=k) + 
                                        s(Microphytoplankton,bs="tp",k=k) + s(Nanophytoplankton,bs="tp",k=k), data = subset2, select = TRUE, method = "REML")           
                            
                            # Save the 10 plots on 1 page
                            message(paste("Plotting response curves for model : ",m, sep = ""))
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net,"/rcurves/", sep = ""))
                        
                            ### Plot smooth termsq
                            pdf(file = paste("plot_rcurves_",resp,"_Oxygen.pdf", sep = ""), width = 10, height = 10)
                                plot(model.new, shade = T, residuals = F, page = 1)
                            dev.off()
                                            
                        } else if ( colnames(subset2)[2] == "Temperature" ) {
                         
                            # Train the modle again but by changing the names of the terms to draw the resp curves
                            colnames(subset2)[c(1:10)] <- c("y","Temperature","Salinity","MLD","PAR","NO2NO3","Chlorophylla",
                                                            "Backscattering","Microphytoplankton","Nanophytoplankton")
                            
                            model.new <- mgcv::gam(formula = y ~ s(Temperature,bs="tp",k=k) + s(Salinity,bs="tp",k=k) + s(MLD,bs="tp",k=k) + 
                                       s(PAR,bs="tp",k=k) + s(NO2NO3,bs="tp",k=k) + s(Chlorophylla,bs="tp",k=k) + s(Backscattering,bs="tp",k=k) + 
                                       s(Microphytoplankton,bs="tp",k=k) + s(Nanophytoplankton,bs="tp",k=k), data = subset2, select = TRUE, method = "REML")           
                        
                            # Save the 10 plots on 1 page
                            message(paste("Plotting response curves for model : ",m, sep = ""))
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net,"/rcurves/", sep = ""))
                        
                            ### Plot smooth termsq
                            pdf(file = paste("plot_rcurves_",resp,"_Temperature.pdf", sep = ""), width = 10, height = 10)
                                plot(model.new, shade = T, residuals = F, page = 1 )
                            dev.off()                   
                    
                        } # eo if else loop
                        
                        setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net, sep = ""))
                                                                      
                 } else {
             
                        message(paste("Stopped because ", group, " is not in counts table", sep = ""))
             
                 } # eo if else loop - for controlling nobs     
                 
             } else {
                 
                 stations2keep <- counts[counts$group == group & counts$n >= 25 ,"station"]
                 subset2 <- subset[subset$Station %in% stations2keep, c(resp,preds)]
                 nobs <- nrow(subset2)
         
                 # Control for a minimum of 15 obs, or if the resp var is a flat one (abundances always = 0)
                 if( nobs >= 20 & sum(subset2[,resp]) > 0 ) {
             
                     k <- round(nrow(subset2)/10, digits = 0)
                     
                     if( colnames(subset2)[2] == "O2" ) {
                         
                            # Train the modle again but by changing the names of the terms to draw the resp curves
                            colnames(subset2)[c(1:10)] <- c("y","Oxygen","Salinity","MLD","PAR","NO2NO3","Chlorophylla",
                                                         "Backscattering","Microphytoplankton","Nanophytoplankton")
                         
                            model.new <- mgcv::gam(formula = y ~ s(Oxygen,bs="tp",k=k) + s(Salinity,bs="tp",k=k) + s(MLD,bs="tp",k=k) + 
                                        s(PAR,bs="tp",k=k) + s(NO2NO3,bs="tp",k=k) + s(Chlorophylla,bs="tp",k=k) + s(Backscattering,bs="tp",k=k) + 
                                        s(Microphytoplankton,bs="tp",k=k) + s(Nanophytoplankton,bs="tp",k=k), data = subset2, select = TRUE, method = "REML")           
                      
                            # Save the 10 plots on 1 page
                            message(paste("Plotting response curves for model : ",m, sep = ""))
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net,"/rcurves/", sep = ""))
                        
                            ### Plot smooth termsq
                            pdf(file = paste("plot_rcurves_",resp,"_Oxygen.pdf", sep = ""), width = 10, height = 10)
                                plot(model.new, shade = T, residuals = F, page = 1)
                            dev.off()
                                                      
                        } else if ( colnames(subset2)[2] == "Temperature" ) {
                         
                            # Train the modle again but by changing the names of the terms to draw the resp curves
                            colnames(subset2)[c(1:10)] <- c("y","Temperature","Salinity","MLD","PAR","NO2NO3","Chlorophylla",
                                                         "Backscattering","Microphytoplankton","Nanophytoplankton")
                         
                            model.new <- mgcv::gam(formula = y ~ s(Temperature,bs="tp",k=k) + s(Salinity,bs="tp",k=k) + s(MLD,bs="tp",k=k) + 
                                        s(PAR,bs="tp",k=k) + s(NO2NO3,bs="tp",k=k) + s(Chlorophylla,bs="tp",k=k) + s(Backscattering,bs="tp",k=k) + 
                                        s(Microphytoplankton,bs="tp",k=k) + s(Nanophytoplankton,bs="tp",k=k), data = subset2, select = TRUE, method = "REML")           
                      
                            # Save the 10 plots on 1 page
                            message(paste("Plotting response curves for model : ",m, sep = ""))
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net,"/rcurves/", sep = ""))
                        
                            ### Plot smooth termsq
                            pdf(file = paste("plot_rcurves_",resp,"_Temperature.pdf", sep = ""), width = 10, height = 10)
                                plot(model.new, shade = T, residuals = F, page = 1 )
                            dev.off()
                    
                        } # eo if else loop
                        
                        setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net, sep = ""))
                
                 } else {
             
                     message(paste("Stopped because not enough data || n = ",nobs, sep = ""))
             
                 } # eo if else loop - for controlling nobs  

             } # eo if else loop - for when group == "Zooplankton" 
                                 
         } else {
             
              message(paste("Stopped because ", group, " is not in counts table", sep = ""))
             
         }   # eo if else loop - for controlling if group is present in the counts data.table  
                  
    } # eo for loop     
    
} # eo for loop 


### Same as above but for abundances GAMs
for(net in nets) {
    
        message(paste("Loading GAM models for net == ",net, sep = ""))
        message(paste("", sep = ""))
        message(paste("", sep = ""))
    
        setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net, sep = ""))
        files.gams <- dir()[grep("gam_",dir())]
        files.gams2 <- files.gams[grep("Ab_",files.gams)] # files.gams2
    
        # And now, for each of these models, plot rcurves 
        # m <- files.gams2[1]
        for(m in files.gams2) {
        
             message(paste("Loading model : ",m, sep = ""))
             setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net, sep = ""))
             model <- get(load(m))
         
             # Extract the terms of the model from the filename
             terms <- do.call(cbind, strsplit(as.character(m),"\\+"))
             terms[1,] <- str_replace_all(terms[1,], "gam_", "")
        
             # Watchout, if Copepoda_unid is the group, then there's an extra underscore to consider
             if( grepl("Copepoda_unid", terms[1,]) | grepl("Calanoida_unidentified", terms[1,]) | 
                 grepl("Gel_carn", terms[1,])| grepl("Gel_FF", terms[1,]) | grepl("Small_grazers", terms[1,]) ) {
             
                     resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2],unlist(strsplit(terms[1,],"_"))[3],sep = "_")
                     pred1 <- unlist(strsplit(terms[1,],"_"))[4]
                     predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
                 
             } else {
             
                     resp <- paste(unlist(strsplit(terms[1,],"_"))[1], unlist(strsplit(terms[1,],"_"))[2], sep = "_")
                     pred1 <- unlist(strsplit(terms[1,],"_"))[3]
                     predlast <- unlist(strsplit(terms[nrow(terms),],"_"))[1]
             
             } # eo if else loop - grepl("other_Copepoda", terms[1,])   
         
             preds <- c(pred1,terms[2,],terms[3,],terms[4,],terms[5,],terms[6,],terms[7,],terms[8,],predlast)      
        
             group <- str_replace_all(resp,"Ab_","")

             subset2 <- na.omit(abund[abund$net == net,c(resp,preds)])
             nobs <- nrow(subset2)
         
             # Control for a minimum of 15 obs, or if the resp var is a flat one (abundances always = 0)
             if( nobs >= 20 & sum(subset2[,resp]) > 0 ) {
             
                         k <- round(nrow(subset2)/10, digits = 0)
                     
                         if( colnames(subset2)[2] == "O2" ) {
                         
                             # Train the modle again but by changing the names of the terms to draw the resp curves
                             colnames(subset2)[c(1:10)] <- c("y","Oxygen","Salinity","MLD","PAR","NO2NO3","Chlorophylla",
                                                             "Backscattering","Microphytoplankton","Nanophytoplankton")
                             
                             model.new <- mgcv::gam(formula = y ~ s(Oxygen,bs="tp",k=k) + s(Salinity,bs="tp",k=k) + s(MLD,bs="tp",k=k) + 
                                        s(PAR,bs="tp",k=k) + s(NO2NO3,bs="tp",k=k) + s(Chlorophylla,bs="tp",k=k) + s(Backscattering,bs="tp",k=k) + 
                                        s(Microphytoplankton,bs="tp",k=k) + s(Nanophytoplankton,bs="tp",k=k), data = subset2, select = TRUE, method = "REML")           
                            
                            # Save the 10 plots on 1 page
                            message(paste("Plotting response curves for model : ",m, sep = ""))
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net,"/rcurves/", sep = ""))
                        
                            ### Plot smooth termsq
                            pdf(file = paste("plot_rcurves_",resp,"_Oxygen.pdf", sep = ""), width = 10, height = 10)
                                plot(model.new, shade = T, residuals = F, page = 1)
                            dev.off()
                                            
                        } else if ( colnames(subset2)[2] == "Temperature" ) {
                         
                            # Train the modle again but by changing the names of the terms to draw the resp curves
                            colnames(subset2)[c(1:10)] <- c("y","Temperature","Salinity","MLD","PAR","NO2NO3","Chlorophylla",
                                                            "Backscattering","Microphytoplankton","Nanophytoplankton")
                            
                            model.new <- mgcv::gam(formula = y ~ s(Temperature,bs="tp",k=k) + s(Salinity,bs="tp",k=k) + s(MLD,bs="tp",k=k) + 
                                       s(PAR,bs="tp",k=k) + s(NO2NO3,bs="tp",k=k) + s(Chlorophylla,bs="tp",k=k) + s(Backscattering,bs="tp",k=k) + 
                                       s(Microphytoplankton,bs="tp",k=k) + s(Nanophytoplankton,bs="tp",k=k), data = subset2, select = TRUE, method = "REML")           
                        
                            # Save the 10 plots on 1 page
                            message(paste("Plotting response curves for model : ",m, sep = ""))
                            setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net,"/rcurves/", sep = ""))
                        
                            ### Plot smooth termsq
                            pdf(file = paste("plot_rcurves_",resp,"_Temperature.pdf", sep = ""), width = 10, height = 10)
                                plot(model.new, shade = T, residuals = F, page = 1)
                            dev.off()                   
                    
                        } # eo if else loop
                        
                        setwd(paste("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/Models/GAMs/",net, sep = ""))
                                                                                      
                 } else {
             
                     message(paste("Stopped because not enough data || n = ",nobs, sep = ""))
             
                 } # eo if else loop - for controlling nobs  
     
    } # eo for loop     
    
} # eo for loop 



# --------------------------------------------------------------------------------------------------------------------------------

### 21/07/2020: To summarize the information stored in the smooth terms curves --> extarct all the fitted data (rcurves) and 
### cbind them in a data.frame as follows: Resp/Net/Pred1/Pred2/Pred3 etc.
### --> then, perform multivariate analyses (PCA, weight by relative rank !) and clustering to identify mdoes of variances !

nets <- c("wp2","bongo","regent")
n <- "bongo"
require("parallel")

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


### Test for dtwclust (multivariate example):
library("dtwclust")

# # Test with your own data?
# resp1 <- temp.models[temp.models$resp == vars[1] & temp.models$net == "wp2",c(3:11)]
# resp2 <- temp.models[temp.models$resp == vars[2] & temp.models$net == "wp2",c(3:11)]
# resp3 <- temp.models[temp.models$resp == vars[3] & temp.models$net == "wp2",c(3:11)]
# resp4 <- temp.models[temp.models$resp == vars[4] & temp.models$net == "wp2",c(3:11)]
# resp5 <- temp.models[temp.models$resp == vars[5] & temp.models$net == "wp2",c(3:11)]
# resp6 <- temp.models[temp.models$resp == vars[6] & temp.models$net == "wp2",c(3:11)]
# resp7 <- temp.models[temp.models$resp == vars[7] & temp.models$net == "wp2",c(3:11)]
# resp8 <- temp.models[temp.models$resp == vars[8] & temp.models$net == "wp2",c(3:11)]
# resp9 <- temp.models[temp.models$resp == vars[9] & temp.models$net == "wp2",c(3:11)]
# # dim(resp1); dim(resp2); dim(resp3); dim(resp4)
# ts.list <- list(resp1 = resp1, resp2 = resp2, resp3 = resp3, resp4 = resp4, resp5 = resp5, resp6 = resp6, resp7 = resp7, resp8 = resp8, resp9 = resp9)
# names(ts.list) <- vars[c(1:9)]
# str(ts.list) #; head(ts.list)
# # Cluster time series based on dtw, try k = 3
# mvc <- tsclust(ts.list, k = 6, distance = "dtw", args = tsclust_args(dist = list(sigma = 100)), type = "partitional")
# plot(mvc)
# data.frame(Resp = vars, k = mvc@cluster)[order(data.frame(Resp = vars, k = mvc@cluster)$k),]

### On a good path :-) evaluate clustering outputs
# ?compare_clusterings #  Compare many different clustering algorithms with support for parallelization.
# Fuzzy preprocessing: calculate autocorrelation up to 50th lag
# acf_fun <- function(series, ...) {
#     lapply(series, function(x) {
#              as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf)
#     })
# } # eo FUN

# # Define overall configuration !!! BEWARE here you spficy the range of k --> should be lower than N vars (so in our case from 2 to 8)
# cfgs <- compare_clusterings_configs(
#          types = c("p", "h", "f", "t"),
#          k = 2L:8L,
#          controls = list(
#              partitional = partitional_control(
#                  iter.max = 30L,
#                  nrep = 1L
#              ),
#              hierarchical = hierarchical_control(
#                  method = "all"
#              ),
#              fuzzy = fuzzy_control(
#                  # notice the vector
#                  fuzziness = c(2,2.5),
#                  iter.max = 30L
#              ),
#              tadpole = tadpole_control(
#                  # notice the vectors
#                  dc = c(1.5,2),
#                  window.size = 2L:8L
#              )
#          ),
#          preprocs = pdc_configs(
#              type = "preproc",
#              # shared
#              none = list(),
#              zscore = list(center = c(FALSE)),
#              # only for fuzzy
#              fuzzy = list(
#                  acf_fun = list()
#              ),
#              # only for tadpole
#              tadpole = list(
#                  reinterpolate = list(new.length = 205L)
#              ),
#              # specify which should consider the shared ones
#              share.config = c("p", "h")
#          ),
#          distances = pdc_configs(
#              type = "distance",
#              sbd = list(),
#              fuzzy = list(
#                  L2 = list()
#              ),
#              share.config = c("p","h")
#          ),
#          centroids = pdc_configs(
#              type = "centroid",
#              partitional = list(
#                  pam = list()
#              ),
#              # special name 'default'
#              hierarchical = list(
#                  default = list()
#              ),
#              fuzzy = list(
#                  fcmdd = list()
#              ),
#              tadpole = list(
#                  default = list(),
#                  shape_extraction = list(znorm = TRUE)
#              )
#          )
# )
#
# num_configs <- sapply(cfgs, attr, which = "num.configs")
# cat("\nTotal number of configurations without considering optimizations:",sum(num_configs),"\n\n")
#
# # Define evaluation functions based on CVI: Variation of Information (only crisp partition)
# vi_evaluators <- cvi_evaluators("VI", ground.truth = vars)
# score_fun <- vi_evaluators$score
# pick_fun <- vi_evaluators$pick
#
# require("doParallel")
# registerDoParallel(cl <- makeCluster(detectCores()))
#
# # Compare the
# comparison_long <- compare_clusterings(ts.list, types = c("p", "h", "f", "t"),
#                         configs = cfgs, seed = 293L, trace = T,
#                         score.clus = score_fun, pick.clus = pick_fun,
#                         return.objects = T
# ) # eo compare_clusterings
#
# # Using all external CVIs and majority vote
# external_evaluators <- cvi_evaluators("external", ground.truth = CharTrajLabels)
# score_external <- external_evaluators$score
# pick_majority <- external_evaluators$pick
#
# comparison_majority <- compare_clusterings(CharTraj, types = c("p", "h", "f", "t"),
#                         configs = cfgs, seed = 84L, trace = T,
#                         score.clus = score_external, pick.clus = pick_majority,
#                         return.objects = T
# ) # eo compare_clusterings
#
#
# # Display best results and close parallel computing
# plot(comparison_majority$pick$object)
# print(comparison_majority$pick$config)
#
# stopCluster(cl); registerDoSEQ()
#
# ### Perform the same as above but based on Temperature only
# resp1 <- temp.models[temp.models$resp == vars[1] & temp.models$net == "wp2","Temperature"]
# resp2 <- temp.models[temp.models$resp == vars[2] & temp.models$net == "wp2","Temperature"]
# resp3 <- temp.models[temp.models$resp == vars[3] & temp.models$net == "wp2","Temperature"]
# resp4 <- temp.models[temp.models$resp == vars[4] & temp.models$net == "wp2","Temperature"]
# resp5 <- temp.models[temp.models$resp == vars[5] & temp.models$net == "wp2","Temperature"]
# resp6 <- temp.models[temp.models$resp == vars[6] & temp.models$net == "wp2","Temperature"]
# resp7 <- temp.models[temp.models$resp == vars[7] & temp.models$net == "wp2","Temperature"]
# resp8 <- temp.models[temp.models$resp == vars[8] & temp.models$net == "wp2","Temperature"]
# resp9 <- temp.models[temp.models$resp == vars[9] & temp.models$net == "wp2","Temperature"]
# # dim(resp1); dim(resp2); dim(resp3); dim(resp4)
# ts.list2 <- list(resp1 = resp1, resp2 = resp2, resp3 = resp3, resp4 = resp4, resp5 = resp5, resp6 = resp6, resp7 = resp7, resp8 = resp8, resp9 = resp9)
# names(ts.list2) <- vars[c(1:9)]
# str(ts.list2) #; head(ts.list)
#
# # Using all external CVIs and majority vote
# external_evaluators <- cvi_evaluators("external", ground.truth = vars)
# score_external <- external_evaluators$score
# pick_majority <- external_evaluators$pick
#
# require("doParallel")
# registerDoParallel(cl <- makeCluster(detectCores()))
#
# comparison_majority <- compare_clusterings(ts.list2, types = c("p","h"),
#                         configs = cfgs, seed = 84L, trace = T,
#                         score.clus = score_external, pick.clus = pick_majority,
#                         return.objects = T
# ) # eo compare_clusterings
#
#
# plot(comparison_majority$pick$object)
# print(comparison_majority$pick$config)
#
# stopCluster(cl); registerDoSEQ()

### Try this: 
# ?tsclust
# klus_1 <- tsclust(ts.list, k = 3L:8L, distance = "dtw_basic", centroid = "pam", type = "partitional")
# names(klus_1) <- paste0("k_", 3L:8L)
# cvi.table1 <- data.frame(sapply(klus_1, cvi, type = "internal"))
# cvi.table1$index <- factor(rownames(cvi.table1))
#
# klus_2 <- tsclust(ts.list, k = 3L:8L, distance = "dtw_basic", centroid = "pam", type = "hierarchical")
# names(klus_2) <- paste0("k_", 3L:8L)
# cvi.table2 <- data.frame(sapply(klus_2, cvi, type = "internal"))
# cvi.table2$index <- factor(rownames(cvi.table2))
#
# # Plot with facet per index
# m.cvi.table1 <- melt(cvi.table1, id.var = "index")
# m.cvi.table2 <- melt(cvi.table2, id.var = "index")
# colnames(m.cvi.table1)[2] <- "k"
# colnames(m.cvi.table2)[2] <- "k"
#
# # Plot
# ggplot(data = m.cvi.table1, aes(x = factor(k), y = value)) + geom_point() + geom_path() +
#         xlab("N clusters") + ylab("Index value") + theme_bw() +
#         facet_wrap(. ~ factor(index), scales = "free")
# #
# ggplot(data = m.cvi.table2, aes(x = factor(k), y = value)) + geom_point() + geom_path() +
#         xlab("N clusters") + ylab("Index value") + theme_bw() +
#         facet_wrap(. ~ factor(index), scales = "free")


### Choose vars to keep (large groups in Tables 1/2)
vars <- resp_esd[c(1,3,4,7,8,10:12,15,17)] ; vars

### Separate in two:
oxy.models <- table[table$formula == "O2+Salinity+MLD+PAR2+NO2NO3+Chla+bbp470+Micro+Nano",]
temp.models <- table[table$formula == "Temperature+Salinity+MLD+PAR2+NO2NO3+Chla+bbp470+Micro+Nano",]
### Change their column names
colnames(oxy.models)[c(3:11)] <- c("Oxygen","Salinity","MLD","PAR","NO2NO3","Chlorophylla","bbp470","%Micro","%Nano")
colnames(temp.models)[c(3:11)] <- c("Temperature","Salinity","MLD","PAR","NO2NO3","Chlorophylla","bbp470","%Micro","%Nano")
# Do they have same dimension?
dim(oxy.models[oxy.models$resp %in% vars,]); dim(temp.models[temp.models$resp %in% vars,]) # good
# summary(oxy.models[oxy.models$resp %in% vars,]) ; summary(temp.models[temp.models$resp %in% vars,])

# Provide Temp smooths to the oxygen table 
oxy.models2 <- oxy.models[oxy.models$resp %in% vars,]
temp.models2 <- temp.models[temp.models$resp %in% vars,]    
oxy.models2$Temperature <- temp.models2$Temperature
head(oxy.models2)

### Add an ID: resp+net
esd.model <- oxy.models2
esd.model$ID <- factor(paste(esd.model$resp, esd.model$net, sep = "_"))
# unique(esd.model$ID)

### 24/07/2020: Only account for those ESD models that show >50% explaiend deviance
mods2choose <- unique(esd.model[esd.model$Deviance >= 0.50,"ID"]) ; mods2choose
esd.model2 <- esd.model[esd.model$ID %in% mods2choose,]
# For each unique(esd.model$ID), extract the corresponding smooth terms series and concatenate in lapply
list.esd <- lapply(unique(esd.model2$ID), function(id) {
                data <- esd.model2[esd.model2$ID == id,c(14,3:11)]
                return(data)
        } # eo lapply
) # eo lapply
# Examine list
#str(list.esd)
# Attribute names
names(list.esd) <- unique(esd.model2$ID)
# cool ! test clusters like above

### PAM clustering (k = 2 to 10)
klus_1 <- tsclust(list.esd, k = 2L:10L, distance = "dtw_basic", centroid = "pam", type = "partitional")
names(klus_1) <- paste0("k_", 2L:10L)
cvi.table1 <- data.frame(sapply(klus_1, cvi, type = "internal"))
cvi.table1$index <- factor(rownames(cvi.table1))

### HAC clustering
klus_2 <- tsclust(list.esd, k = 2L:10L, distance = "dtw_basic", type = "hierarchical")
names(klus_2) <- paste0("k_", 2L:10L)
cvi.table2 <- data.frame(sapply(klus_2, cvi, type = "internal"))
cvi.table2$index <- factor(rownames(cvi.table2))

# The indices marked with an exclamation mark (!) calculate (or
# re-use if already available) the whole distance matrix between the
# series in the data. If you were trying to avoid this in the first
# place, then these CVIs might not be suitable for your application.
#
# The indices marked with a question mark (?) depend on the
# extracted centroids, so bear that in mind if a hierarchical
# procedure was used and/or the centroid function has associated
# randomness (such as ‘shape_extraction()’ with series of different length).
#
# The indices marked with a tilde (~) require the calculation of a
# global centroid. Since ‘DBA()’ and ‘shape_extraction()’ (for
# series of different length) have some randomness associated, these
# indices might not be appropriate for those centroids.
#
#         • Crisp partitions
#
#             • ‘"Sil"’ (!): Silhouette index (Rousseeuw (1987); to be
#               maximized).
#
#             • ‘"D"’ (!): Dunn index (Arbelaitz et al. (2013); to be
#               maximized).
#
#             • ‘"COP"’ (!): COP index (Arbelaitz et al. (2013); to be
#               minimized).
#
#             • ‘"DB"’ (?): Davies-Bouldin index (Arbelaitz et al.
#               (2013); to be minimized).
#
#             • ‘"DBstar"’ (?): Modified Davies-Bouldin index (DB*) (Kim
#               and Ramakrishna (2005); to be minimized).
#
#             • ‘"CH"’ (~): Calinski-Harabasz index (Arbelaitz et al.
#               (2013); to be maximized).
#
#             • ‘"SF"’ (~): Score Function (Saitta et al. (2007); to be
#               maximized; see notes).
              
# Plot with facet per index
m.cvi.table1 <- melt(cvi.table1, id.var = "index")
colnames(m.cvi.table1)[2] <- "k"
m.cvi.table1$k <- str_replace_all(m.cvi.table1$k,"k_","")

m.cvi.table2 <- melt(cvi.table2, id.var = "index")
colnames(m.cvi.table2)[2] <- "k"
m.cvi.table2$k <- str_replace_all(m.cvi.table2$k,"k_","")

### Change the names of the indices for better ones and remove COP and SF
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "CH"] <- "Calinski-Harabasz index (maximize)"
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "Sil"] <- "Silhouette index (maximize)"
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "D"] <- "Dunn index (maximize)"
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "DB"] <- "Davies-Bouldin index (minimize)"
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "DBstar"] <- "DB* (minimize)"

levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "CH"] <- "Calinski-Harabasz index (maximize)"
levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "Sil"] <- "Silhouette index (maximize)"
levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "D"] <- "Dunn index (maximize)"
levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "DB"] <- "Davies-Bouldin index (minimize)"
levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "DBstar"] <- "DB* (minimize)"

quartz()
plot <- ggplot(data = m.cvi.table1[!(m.cvi.table1$index %in% c("SF","COP")),],
                aes(x = as.numeric(k), y = value)) + geom_point() +
        xlab("N clusters (DTW+PAM)") + ylab("Index value") + theme_bw() + 
        scale_x_continuous(labels = c(2:10), breaks = c(2:10)) + 
        facet_wrap(. ~ factor(index), scales = "free")
        
ggsave(plot = plot, filename = "panel_cvi_dtw_pam_ESD.pdf", dpi = 300, height = 5, width = 9)
        
#
quartz()
ggplot(data = m.cvi.table2[!(m.cvi.table2$index %in% c("SF","COP")),],
            aes(x = as.numeric(k), y = value)) + geom_point() +
        xlab("N clusters (HAC)") + ylab("Index value") + theme_bw() + 
        scale_x_continuous(labels = c(2:15), breaks = c(2:15)) + 
        facet_wrap(. ~ factor(index), scales = "free")

### CONCLUSION
# - PAM should be: 3, 4 or 6
# - HAC should be: ... strong disagreement :s

### Try PAM
clusters.pam <- tsclust(list.esd, k = 4L, distance = "dtw_basic", centroid = "pam", type = "partitional")
plot(clusters.pam)
ddf_k <- data.frame(Resp = unique(esd.model2$ID), k = clusters.pam@cluster)[order(data.frame(Resp = unique(esd.model2$ID), k = clusters.pam@cluster)$k),]
#ddf_k

### Plot prototype of each cluster :-)
for(i in 1:length(clusters.pam@cluster)) { 
    quartz() ; plot(clusters.pam, type = "series", clus = i)
    quartz() ; plot(clusters.pam, type = "centroids", clus = i)
}

### --> 4 looks like the best

### Try HAC
# clusters.hac <- tsclust(list.esd, k = 4, distance = "dtw_basic", args = tsclust_args(dist = list(sigma = 100)), type = "hierarchical")
# # plot(clusters.hac, hang = 1)
# data.frame(Resp = unique(esd.model2$ID), k = clusters.hac@cluster)[order(data.frame(Resp = unique(esd.model2$ID), k = clusters.hac@cluster)$k),]
#
# ### Plot prototype of each cluster :-)
# for(i in 1:length(clusters.hac@cluster)) {
#         quartz() ; plot(clusters.hac, type = "series", clus = i)
# }
### forget HAC clustering

str(clusters.pam)
str(clusters.pam@distmat)
class(clusters.pam@distmat)
mat <- data.frame(as.dist(clusters.pam@distmat))
mat$group <- rownames(mat)
melt <- melt(mat, id.vars = "group") ; head(melt)

ggplot(data = melt) + geom_tile(aes(x = factor(group), y = factor(variable),
        fill = value/max(melt$value) ), colour = "black") + coord_fixed() +
    xlab("") + ylab("") + scale_fill_viridis(name = "Distance") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Or, cmdscale
cmd <- cmdscale(as.dist(clusters.pam@distmat), k = 2, eig = T)
class(cmd)
str(cmd)
ddf <- data.frame(cmd$point)
colnames(ddf) <- c("x","y")
ddf$group <- rownames(ddf)
# data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))
ddf$zoo.group <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,2]
ddf$net <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,3]

ddf$cluster <- NA
for(g in unique(ddf$zoo.group)) {
    k <- ddf_k[ddf_k$Resp == g,"k"]
    ddf[ddf$group == g,"cluster"] <- k
}


ddf$cluster <- NA
ddf[ddf$group %in% c("ESD_Calanoida_wp2","ESD_Chaetognatha_bongo","ESD_Copepoda_bongo",
    "ESD_Eumalacostraca_bongo","ESD_Eumalacostraca_regent","ESD_Poecilostomatoida_wp2","ESD_Rhizaria_wp2","ESD_Zooplankton_bongo"),"cluster"] <- 1
#
ddf[ddf$group == "ESD_Cnidaria_wp2","cluster"] <- 2
ddf[ddf$group == "ESD_Chaetognatha_regent","cluster"] <- 3
ddf[ddf$group %in% c("ESD_Calanoida_regent","ESD_Cnidaria_bongo","ESD_Copepoda_regent"),"cluster"] <- 4

ddf$id <- factor(paste(ddf$zoo.group," ","(",ddf$net,")", sep = ""))

ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_point(aes(x = x, y = y, fill = factor(cluster)), data = ddf, colour = "black", size = 3, pch = 21) + 
    geom_text_repel(aes(x = x, y = y, label = id), data = ddf) + 
    theme_bw()


### 24/07/2020: Replot the profiles per cluster (column) and covariate (rows)
# First, supply ddf_k$k to 'esd.model2'
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
#dim(m.esd[m.esd$ID == "ESD_Calanoida_wp2" & m.esd$variable == "Oxygen",])

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

ggsave(plot = facet, filename = "panel_smooth_clustersxcovariates_ESD.pdf", dpi = 300, height = 13, width = 9)

### Cool.



### ---------------------------------------------------------------

### Do the same for abundances
vars2 <- resp_ab[c(4,8,9,12,13,16,17,24,26,28,32)] ; vars2

# Provide Temp smooths to the oxygen table 
oxy.models2 <- oxy.models[oxy.models$resp %in% vars2,]
temp.models2 <- temp.models[temp.models$resp %in% vars2,]    
oxy.models2$Temperature <- temp.models2$Temperature
head(oxy.models2)

### Add an ID: resp+net
ab.model <- oxy.models2
ab.model$ID <- factor(paste(ab.model$resp, ab.model$net, sep = "_"))
unique(ab.model$ID)

### 24/07/2020: Only account for those ESD models that show >50% explaiend deviance
mods2choose <- unique(ab.model[ab.model$Deviance >= 0.40,"ID"]) ; mods2choose
ab.model2 <- ab.model[ab.model$ID %in% mods2choose,]
# Need to remove the duplicated rows though
ab.model2 <- ab.model2 %>% distinct()
# For each unique(esd.model$ID), extract the corresponding smooth terms series and concatenate in lapply
list.ab <- lapply(unique(ab.model2$ID), function(id) {
                data <- ab.model2[ab.model2$ID == id,c(14,3:11)]
                return(data)
        } # eo lapply
) # eo lapply
# Examine list
# str(list.ab)
# Attribute names
names(list.ab) <- unique(ab.model2$ID)
# cool ! test clusters like above

### PAM clustering
klus_1 <- tsclust(list.ab, k = 2L:10L, distance = "dtw_basic", centroid = "pam", type = "partitional")
names(klus_1) <- paste0("k_", 2L:10L)
cvi.table1 <- data.frame(sapply(klus_1, cvi, type = "internal"))
cvi.table1$index <- factor(rownames(cvi.table1))

# ### HAC clustering
# klus_2 <- tsclust(list.ab, k = 2L:10L, distance = "dtw_basic", type = "hierarchical")
# names(klus_2) <- paste0("k_", 2L:10L)
# cvi.table2 <- data.frame(sapply(klus_2, cvi, type = "internal"))
# cvi.table2$index <- factor(rownames(cvi.table2))

# Plot with facet per index
m.cvi.table1 <- melt(cvi.table1, id.var = "index")
colnames(m.cvi.table1)[2] <- "k"
m.cvi.table1$k <- str_replace_all(m.cvi.table1$k,"k_","")

# m.cvi.table2 <- melt(cvi.table2, id.var = "index")
# colnames(m.cvi.table2)[2] <- "k"
# m.cvi.table2$k <- str_replace_all(m.cvi.table2$k,"k_","")

### Change the names of the indices for better ones and remove COP and SF
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "CH"] <- "Calinski-Harabasz index (maximize)"
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "Sil"] <- "Silhouette index (maximize)"
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "D"] <- "Dunn index (maximize)"
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "DB"] <- "Davies-Bouldin index (minimize)"
levels(m.cvi.table1$index)[levels(m.cvi.table1$index) == "DBstar"] <- "DB* (minimize)"

# levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "CH"] <- "Calinski-Harabasz index (maximize)"
# levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "Sil"] <- "Silhouette index (maximize)"
# levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "D"] <- "Dunn index (maximize)"
# levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "DB"] <- "Davies-Bouldin index (minimize)"
# levels(m.cvi.table2$index)[levels(m.cvi.table2$index) == "DBstar"] <- "DB* (minimize)"

# Plot
quartz()
plot <- ggplot(data = m.cvi.table1[!(m.cvi.table1$index %in% c("SF","COP")),],
            aes(x = as.numeric(k), y = value)) + geom_point() +
        xlab("N clusters (PAM)") + ylab("Index value") + theme_bw() + 
        scale_x_continuous(labels = c(2:15), breaks = c(2:15)) + 
        facet_wrap(. ~ factor(index), scales = "free")
ggsave(plot = plot, filename = "panel_cvi_dtw_pam_abund.pdf", dpi = 300, height = 5, width = 9)

# quartz()
# ggplot(data = m.cvi.table2[!(m.cvi.table2$index %in% c("SF","COP")),],
#             aes(x = as.numeric(k), y = value)) + geom_point() +
#         xlab("N clusters (HAC)") + ylab("Index value") + theme_bw() +
#         scale_x_continuous(labels = c(2:15), breaks = c(2:15)) +
#         facet_wrap(. ~ factor(index), scales = "free")

### CONCLUSIONS:
# - PAM should be 6 clearly --> clearer patterns than HAC so stick with it
# - HAC should be ... 3?

### Try PAM
clusters.pam <- tsclust(list.ab, k = 5L, distance = "dtw_basic", centroid = "pam", type = "partitional")
plot(clusters.pam)
ddf_k <- data.frame(Resp = unique(ab.model2$ID), k = clusters.pam@cluster)[order(data.frame(Resp = unique(ab.model2$ID), k = clusters.pam@cluster)$k),]

### Plot prototype of each cluster :-)
# for(i in 1:length(clusters.pam@cluster)) {
#         #quartz() ; plot(clusters.pam, type = "centroid", clus = i)
#         quartz() ; plot(clusters.pam, type = "centroid", clus = i)
# }

#str(clusters.pam)
#str(clusters.pam@distmat)
class(clusters.pam@distmat)
mat <- data.frame(as.dist(clusters.pam@distmat))
mat$group <- rownames(mat)
melt <- melt(mat, id.vars = "group") ; head(melt)

ggplot(data = melt) + geom_tile(aes(x = factor(group), y = factor(variable),
        fill = value/max(melt$value) ), colour = "black") + coord_fixed() +
    xlab("") + ylab("") + scale_fill_viridis(name = "Distance") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Or, cmdscale
cmd <- cmdscale(as.dist(clusters.pam@distmat), k = 2, eig = T)
#class(cmd)
#str(cmd)
ddf <- data.frame(cmd$point)
colnames(ddf) <- c("x","y")
ddf$group <- rownames(ddf)
# data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))
ddf$zoo.group <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,2]
ddf$net <- data.frame(do.call(rbind, strsplit(ddf$group, split = "_")))[,3]

### Provide cluster to ddf from ddf_k
head(ddf) ; head(ddf_k)

ddf$cluster <- NA
for(g in unique(ddf$group)) {
    k <- ddf_k[ddf_k$Resp == g,"k"]
    ddf[ddf$group == g,"cluster"] <- k
}

# ddf[ddf$group %in% c("Ab_Calanoida_regent","Ab_Copepoda_regent","Ab_Zooplankton_regent"),"cluster"] <- 5
# ddf[ddf$group %in% c("Ab_Chaetognatha_bongo","Ab_Pteropoda_bongo","Ab_Rhizaria_bongo"),"cluster"] <- 4
# ddf[ddf$group %in% c("Ab_Chaetognatha_regent","Ab_Cnidaria_wp2","Ab_Pteropoda_regent","Ab_Tunicata_bongo"),"cluster"] <- 1
# ddf[ddf$group %in% c("Ab_Tunicata_wp2"),"cluster"] <-
# ddf[ddf$group %in% c("Ab_Poecilostomatoida_wp2"),"cluster"] <-

ddf$id <- factor(paste(ddf$zoo.group," ","(",ddf$net,")", sep = ""))

ggplot() + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_point(aes(x = x, y = y, fill = factor(cluster)), data = ddf, colour = "black", size = 3, pch = 21) + 
    geom_text_repel(aes(x = x, y = y, label = id), data = ddf) + 
    theme_bw()


### Better panel plot, like with ESD
head(ab.model2)
ab.model2$cluster_pam <- NA
for(id in unique(ab.model2$ID)) {
    k <- ddf_k[ddf_k$Resp == id,"k"]
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

### 02/09/2020: Move Ab_Cyclopoida (WP2)
to 
head(m.ab)

facet <- ggplot(data = m.ab) + geom_line(aes(x = x, y = value, colour = factor(ID2))) +
    xlab("") + ylab("") + theme_bw() + facet_grid(factor(variable) ~ factor(cluster_pam), scales = "free") + 
    scale_colour_discrete(name = "") + 
    geom_hline(yintercept = 0, linetype = "dashed")

ggsave(plot = facet, filename = "panel_smooth_clustersxcovariates_abund.pdf", dpi = 300, height = 13, width = 9)

### Cool.


