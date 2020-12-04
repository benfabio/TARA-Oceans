
##### 29/11/19: R script to examine the overlap between the hydrological parameters measured by the 
### To do list:
#   - Combine the hydro final with the hydro meadured by imagery (overall env conditions at stations versus those measured at sampling event)
#   - Compare their distrbutions (violins+boxplots) for each net and parameters
#   - Also compare the transformed 

### last update: 26/05/2020

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
library("ggrepel")

world <- map_data("world") # coastline for maps
WD <- getwd()

### ---------------------------------------------------------------------------------------------------------------------------

### Load the hydrological data from each net
imagery <- read.csv("TARA_OCEANS_table_imagery_data.csv", h = T, sep = ";")
colnames(imagery)[3] <- "net"
colnames(imagery)

### Load the mean hydrological data given by Fabien Lombard (mean conditions over a sampling station)
hydro <- read.csv("Hydro_final_TARA_08_03_20.csv", h = T, sep = ",")
dim(hydro)
str(hydro)
colnames(hydro)

unique(imagery$station) # length(unique(imagery$station)) --> 166
unique(hydro$Station)   # length(unique(hydro$Station))   --> 210

# Exclude stations in hydro which are not in imagery (should be zero...)
unique(imagery[imagery$station %in% unique(hydro$Station),"station"]) # 166
unique(hydro[hydro$Station %in% unique(imagery$station),"Station"]) # 166
# OK :-) 

### Variables to examine and their equivalent in "hydro", compare n NAs
Nbongo <- nrow(imagery[imagery$net == "bongo",]); Nbongo
stations.bongo <- unique(imagery[imagery$net == "bongo","station"])
Nwp2 <- nrow(imagery[imagery$net == "wp2",]); Nwp2
stations.wp2 <- unique(imagery[imagery$net == "wp2","station"])
Nregent <- nrow(imagery[imagery$net == "regent",]); Nregent
stations.regent <- unique(imagery[imagery$net == "regent","station"])

# "x_Micro"         --> micro     
sum(is.na(imagery[imagery$net == "bongo","x_Micro"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"micro"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","x_Micro"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"micro"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","x_Micro"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"micro"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

# "x_Nano"          --> nano
sum(is.na(imagery[imagery$net == "bongo","x_Nano"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"nano"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","x_Nano"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"nano"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","x_Nano"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"nano"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

# "x_Pico"          --> pico   
sum(is.na(imagery[imagery$net == "bongo","x_Pico"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"pico"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","x_Pico"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"pico"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","x_Pico"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"pico"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
                
# "min_o2_depth"    --> NA
# "DCM"             --> NA                          
# "z_eu"            --> Zeu
sum(is.na(imagery[imagery$net == "bongo","z_eu"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Zeu"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","z_eu"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Zeu"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","z_eu"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Zeu"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

# "T_10m"           --> Tzeu
sum(is.na(imagery[imagery$net == "bongo","T_10m"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Tzeu"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","T_10m"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Tzeu"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","T_10m"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Tzeu"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
                    
# "T_zeu"           --> Tzeu
sum(is.na(imagery[imagery$net == "bongo","T_zeu"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Tzeu"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","T_zeu"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Tzeu"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","T_zeu"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Tzeu"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

# "T_mld"           --> Tzeu
sum(is.na(imagery[imagery$net == "bongo","T_mld"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Tzeu"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","T_mld"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Tzeu"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","T_mld"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Tzeu"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
                  
# "S_10m"           --> Szeu
sum(is.na(imagery[imagery$net == "bongo","S_10m"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","S_10m"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","S_10m"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

# "S_zeu"           --> Szeu
sum(is.na(imagery[imagery$net == "bongo","S_zeu"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","S_zeu"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","S_zeu"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
                 
# "S_mld"           --> Szeu
sum(is.na(imagery[imagery$net == "bongo","S_mld"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","S_mld"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","S_mld"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
sum(is.na(imagery[imagery$net == "bongo","S_mld"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","S_mld"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","S_mld"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Szeu"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
                 
# "o2_10m"          --> O2_10m
sum(is.na(imagery[imagery$net == "bongo","o2_10m"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"O2_10m"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","o2_10m"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"O2_10m"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","o2_10m"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"O2_10m"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

# "o2_zeu"          --> O2_10m
sum(is.na(imagery[imagery$net == "bongo","o2_zeu"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"O2_10m"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","o2_zeu"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"O2_10m"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","o2_zeu"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"O2_10m"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
                
# "o2_mld"          --> O2_10m
sum(is.na(imagery[imagery$net == "bongo","o2_mld"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"O2_10m"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","o2_mld"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"O2_10m"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","o2_mld"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"O2_10m"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

# "Chla_10m"        --> Chl10m ? chla_ss ? chla_int ? chlatot ? Chla_hplc ?
sum(is.na(imagery[imagery$net == "bongo","Chla_10m"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Chl10m"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","Chla_10m"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Chl10m"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","Chla_10m"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Chl10m"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

# "Chla_zeu"        --> Chl10m ? chla_ss ? chla_int ? chlatot ? Chla_hplc ?

# "Chla_mld"        --> Chl10m ? chla_ss ? chla_int ? chlatot ? Chla_hplc ? 
  
# "Chla"            --> Chl10m ? chla_ss ? chla_int ? chlatot ? Chla_hplc ? 
sum(is.na(imagery[imagery$net == "bongo","Chla"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"chla_ss"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","Chla"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"chla_ss"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","Chla"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"chla_ss"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
  
# "PO4"             --> PO4 
sum(is.na(imagery[imagery$net == "bongo","PO4"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"PO4"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","PO4"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"PO4"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","PO4"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"PO4"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
                 
# "NO3"             --> NO3_10m
sum(is.na(imagery[imagery$net == "bongo","NO3"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"NO3_10m"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","NO3"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"NO3_10m"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","NO3"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"NO3_10m"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

# "Si"              --> Si
sum(is.na(imagery[imagery$net == "bongo","Si"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"Si"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","Si"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"Si"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","Si"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"Si"]))/ nrow(hydro[hydro$Station %in% stations.regent,])
                 
# "Fluo"            --> NA
# "bac660"          --> NA                 
# "PAR"             --> PAR
sum(is.na(imagery[imagery$net == "bongo","PAR"]))/ Nbongo 
sum(is.na(hydro[hydro$Station %in% stations.bongo,"PAR"]))/ nrow(hydro[hydro$Station %in% stations.bongo,])

sum(is.na(imagery[imagery$net == "wp2","PAR"]))/ Nwp2 
sum(is.na(hydro[hydro$Station %in% stations.wp2,"PAR"]))/ nrow(hydro[hydro$Station %in% stations.wp2,])

sum(is.na(imagery[imagery$net == "regent","PAR"]))/ Nregent 
sum(is.na(hydro[hydro$Station %in% stations.regent,"PAR"]))/ nrow(hydro[hydro$Station %in% stations.regent,])

### --> SO the only increase in benefit is actually for PAR...

### Add PAR to imagery and re-save 
imagery$PAR2 <- NA
for(s in unique(imagery$station) ) {
    
        par2add <- hydro[hydro$Station == s,"PAR"]
        
        if( exists("par2add") ) {
            message(paste("Adding PAR value for station ", s, sep = ""))
            imagery[imagery$station == s,"PAR2"] <- par2add
        } else {
            message(paste("No PAR value to add for station ", s, sep = ""))
            imagery[imagery$station == s,"PAR2"] <- NA
        } # eo else if loop
        
} # eo for loop - s in unique(imagery$station)
summary(imagery[,c("PAR","PAR2")])
write.table(imagery, file = "TARA_OCEANS_table_imagery_dataV2.txt", sep = ";")
# Save imagery :-)


### ----------------------------------------------------

### 17/03/2020: Examine sample size for each categories across stations (to determine categories to consider for ESD patterns)
raw <- read.table("TARA_Zooplankton_raw_data_V2.txt", h = T, sep = "\t")
dim(raw); colnames(raw)
head(raw)
unique(raw$Dataset) # wp2      bongo    regent   multinet
summary(raw[raw$Dataset == "wp2",])
raw[raw$Dataset == "wp2","Cyclopoida.Corycaeidae"]
### For WP2, OK to look at: Calanoida/Cyclopoida/Oithonidae/Corycaeidae/Oncaeidae

summary(raw[raw$Dataset == "bongo",])
raw[raw$Dataset == "bongo","Calanoida.Augaptilidae"]

### Check if zero were translated into NaN as they should
raw[raw$Dataset == "bongo",c("Station","Calanoida.Calanidae")]
imagery2[imagery2$net == "bongo",c("station","ESD_Calanidae")]
imagery2[imagery2$net == "bongo",c("station","Abund_Calanidae")]

### And the untransformed data?
imagery <- read.table("TARA_OCEANS_table_imagery_data_V3.txt", h = T, sep = "\t")
dim(imagery); head(imagery); colnames(imagery)

### Make sure N in 'raw' match abundances in imagery, for various nets
raw[raw$Dataset == "wp2",c("Station","Calanoida.Augaptilidae")]
imagery[imagery$net == "wp2",c("station","Abund_Augaptilidae")]

raw[raw$Dataset == "bongo",]

### Melt and plot for each three main nets
colnames(raw)
m_raw <- melt(raw, id.vars = c("Station","Dataset"))
head(m_raw)
colnames(m_raw) <- c("Station","Net","Group","n")

unique(m_raw$Group) # rename
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Acartiidae"] <- "Acartiidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Aetideidae"] <- "Aetideidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Augaptilidae"] <- "Augaptilidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Calanidae"] <- "Calanidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Candaciidae"] <- "Candaciidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Centropagidae"] <- "Centropagidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Eucalanidae"] <- "Eucalanidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Euchaetidae"] <- "Euchaetidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Paracalanidae"] <- "Paracalanidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Pontellidae"] <- "Pontellidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Rhincalanidae"] <- "Rhincalanidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida.Temoridae"] <- "Temoridae"
levels(m_raw$Group)[levels(m_raw$Group) == "Cyclopoida.Oithonidae"] <- "Oithonidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Cyclopoida.Corycaeidae"] <- "Corycaeidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Cyclopoida.Oncaeidae"] <- "Oncaeidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Cyclopoida.Sapphirinidae"] <- "Sapphirinidae"
levels(m_raw$Group)[levels(m_raw$Group) == "Calanoida"] <- "Calanoida (other)"
levels(m_raw$Group)[levels(m_raw$Group) == "Cyclopoida"] <- "Cyclopoida (other)"
levels(m_raw$Group)[levels(m_raw$Group) == "Other.Copepoda"] <- "Copepoda (other)"

unique(m_raw$Net) # rename
levels(m_raw$Net)[levels(m_raw$Net) == "wp2"] <- "WP2"
levels(m_raw$Net)[levels(m_raw$Net) == "bongo"] <- "Bongo"
levels(m_raw$Net)[levels(m_raw$Net) == "regent"] <- "Régent"

quartz()
ggplot(aes(x = factor(Group), y = log(n), fill = factor(Group)), data = m_raw[m_raw$Net %in% c("Bongo","WP2","Régent"),]) +
    geom_boxplot(colour = "black") + xlab("") + ylab("") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log(30)), linetype = "dashed") + 
    facet_wrap(~factor(Net), ncol = 2, scales = "fixed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()

require("forcats")
# fct_reorder(factor(Group), log(n),.desc = T) us this as x argument to reorer boxplots
quartz()
ggplot(aes(x = fct_reorder(factor(Group), log(n),.desc = T), y = log(n), fill = factor(Group)), data = m_raw[m_raw$Net %in% c("Bongo","WP2","Régent"),]) +
    geom_boxplot(colour = "black") + xlab("") + ylab("") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log(30)), linetype = "dashed") + 
    facet_wrap(~factor(Net), ncol = 2, scales = "fixed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
#
#quartz()
plotA <- ggplot(aes(x = fct_reorder(factor(Group), log(n),.desc = T), y = log(n), fill = factor(Group)), data = m_raw[m_raw$Net == "Bongo",]) +
    geom_boxplot(colour = "black") + xlab("") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log(30)), linetype = "dashed") + 
    scale_y_continuous(limits = c(0,12.5), name = "N individuals (logged)\nin Bongo samples") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))   
#
#quartz()
plotB <- ggplot(aes(x = fct_reorder(factor(Group), log(n),.desc = T), y = log(n), fill = factor(Group)), data = m_raw[m_raw$Net == "WP2",]) +
    geom_boxplot(colour = "black") + xlab("") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log(30)), linetype = "dashed") + 
    scale_y_continuous(limits = c(0,12.5), name = "N individuals (logged)\nin WP2 samples") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
#
#quartz()
plotC <- ggplot(aes(x = fct_reorder(factor(Group), log(n),.desc = T), y = log(n), fill = factor(Group)), data = m_raw[m_raw$Net == "Régent",]) +
    geom_boxplot(colour = "black") + xlab("") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log(30)), linetype = "dashed") + 
    scale_y_continuous(limits = c(0,12.5), name = "N individuals (logged)\nin Régent samples") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
#

require("ggpubr")
# quartz()
panel <- ggarrange(plotA, plotB, plotC, ncol = 1, nrow = 3, common.legend = T)
ggsave(plot = panel, filename = "plots_distrib_Ngroups_net.pdf", width = 7, height = 12, dpi = 300)

# Nice
summary(imagery[imagery$net == "bongo",c("station","ESD_Cyclopoida")])
summary(imagery[imagery$net == "bongo",c("station","ESD_Oncaeidae")])

summary(imagery[imagery$net == "wp2",c("station","ESD_Cyclopoida")])
summary(imagery[imagery$net == "wp2",c("station","ESD_Oithonidae")])
### If ESD_Cyclopoida was indeed based on the ESD measurements of all organisms belonging to the Cyclopoida, there would be as many data points as in ESD_Oithonidae
    
