
### 28/05/2020: Script to analyze the updated zooplankton imaging data from TARA Oceans (WP2, Bongo & Regent).
### Aims to: 
### - For each net, load the median size data and the folders containing the colnames (variables/ taxonomic group names)
### - Put them together and filter the variables you don't need (Chrlorophytes etc. )
### - Use the 'numbers' file to compute the N per group and stations


### Last update: 16/06/2020

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("RColorBrewer")
library("reshape2")
library("viridis")
library("scales")
library("maps")

world <- map_data("world")
WD <- getwd()

# --------------------------------------------------------------------------------------------------------------------------------

### Define vector of net names to perform all subsequent data manipulation in a lapply()
nets <- c("WP2","bongo","regent")
#n <- "WP2"
res <- lapply(nets, function(n) {
    
        message(paste("Preparing data for net ",n,sep = ""))
        
        # Check data structure for WP2 net, should be the same for all nets
        setwd(paste(WD,"/",n,"/", sep = "")) # ; dir()
        # Load
        pft_colnames <- read.csv(paste(n,"speciesnamespft.csv", sep = ""), h = F, sep = ",") #; dim(pft_colnames)
        families_colnames <- read.csv(paste(n,"speciesnamespftfamilies.csv", sep = ""), h = F, sep = ";") #; dim(families_colnames)
        order_colnames <- read.csv(paste(n,"speciesnamespftorder.csv", sep = ""), h = F, sep = ";") #; dim(order_colnames)
        # And load the data
        esd_pft <- read.csv(paste(n,"mediansizespft.csv", sep = ""), h = F, sep = ",") #; dim(esd_pft) ; nrow(pft_colnames)
        esd_ord <- read.csv(paste(n,"mediansizespftorder.csv", sep = ""), h = F, sep = ",") #; dim(esd_ord) ; nrow(order_colnames)
        esd_fam <- read.csv(paste(n,"mediansizespftfamilies.csv", sep = ""), h = F, sep = ",") #; dim(esd_fam) ; nrow(families_colnames)
        # summary(esd_pft)
        colnames(esd_pft) <- pft_colnames[,1] ; colnames(esd_ord) <- order_colnames[,1] ; colnames(esd_fam) <- families_colnames[,1]
        # Check
        # summary(esd_pft); summary(esd_ord); summary(esd_fam)

        ### Remove those columns that only have NAs
        not_all_na <- function(x) any(!is.na(x))
        esd_pft <- esd_pft %>% select_if(not_all_na)
        esd_ord <- esd_ord %>% select_if(not_all_na)
        esd_fam <- esd_fam %>% select_if(not_all_na)
        # summary(esd_pft$chaetognatha); summary(esd_ord$chaetognatha); summary(esd_fam$chaetognatha)
        # Gut, the valeus are the same for common groups, now choose the columns to choose for each table
        # colnames(esd_pft) 
        # all_copepoda, chaetognatha, gelatinous_carnivorous, gelatinous_filter_feeders, large_crustracean, pteropoda, rhizaria, small_grazers
        # colnames(esd_ord); so first 5
        # colnames(esd_fam) 
        ### so first 20 except 7
        sizes <- esd_pft[,c("all_copepoda","chaetognatha","gelatinous_carnivorous","gelatinous_filter_feeders",
                            "large_crustracean","pteropoda","rhizaria","small_grazers")]
        # summary(sizes)
        sizes <- cbind(sizes, esd_ord[,c(1:5)] )
        sizes <- cbind(sizes, esd_fam[,c(1:6,8:20)] )
        # summary(sizes)
        colnames(sizes)[c(1:13)] <- c("Copepoda","Chaetognatha","Gel_carn","Gel_FF","Crust","Pteropoda",
                                    "Rhizaria","Small_grazers","Copepoda_unid","Calanoida","Cyclopoida","Harpacticoida",
                                    "Poecilostomatoida")
        # colnames(sizes)
        if(n == "WP2") {
            n <- "wp2"
        }
        sizes$net <- n
        sizes$station <- c(1:nrow(sizes))
        return(sizes)
    
        } # eo n
    
) # eo lapply
# Rbind, check
table <- dplyr::bind_rows(res)
dim(table); str(table)
summary(table) # length(table) 37 columns, ONE of those is never empty because ot's the net name 
# Count NAs per rows, remove stations that only have NA (count = 36)
table$nNA <- apply(table, 1, function(x) sum(is.na(x)))
nrow(table[table$nNA == 36,]) # 241 ! 

table2 <- table[table$nNA != 36,]
dim(table2)
summary(table2)
length( unique(table2$station) ) # 168 stations !

### Check if there are stations with very few esd estimates
table2[table2$nNA > 30 & table2$nNA != 36,] # OK, a station with only large Calanoida (Calanidae+Metridinidae)

### Now, go back to the dataset you were using previously and check if dimensions fit
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/data/") ; dir()
data <- read.table("TARA_OCEANS_table_imagery2_transformed.txt", h = T, sep = "\t")
dim(data) # 469
length( unique(as.numeric(data$station)) ) # 166
# 163 stations, so 2 new stations were added ^^
colnames(data) # cbind columns 7, 8 and 67:105

### Check if the env measurements are net-dependent in 'data 
data[which(data$net == "wp2" & data$station == 92),c("T_10m","S_10m","Si")]
data[which(data$net == "bongo" & data$station == 92),c("T_10m","S_10m","Si")]
data[which(data$net == "regent" & data$station == 92),c("T_10m","S_10m","Si")]

### Create empty columns with names from colnames(data)[c(7,8,67:105)]
table3 <- table2
table3[,colnames(data)[c(7,8,67:105)]] <- NA
head(table3)
colnames(table3)
# Define common stations
commons <- intersect(unique(table3$station), unique(data$station))
commons ; length(commons)

### Fill with double for loop: per cells and per net
n <- "bongo"
c <- 58
for(n in c("wp2","bongo","regent")) {
    
    message(paste("                   ", sep = ""))
    message(paste("Filling for net == ", n, sep = ""))
    message(paste("                   ", sep = ""))
    
    for(c in commons) {
    
        message(paste(c, sep = ""))
        dat2replace <- data[which(data$station == c & data$net == n),colnames(data)[c(7,8,67:105)]]
        
        if(nrow(dat2replace) > 0) {
        
            table3[which(table3$station == c & table3$net == n),c(40:80)] <- data[which(data$station == c & data$net == n),colnames(data)[c(7,8,67:105)]]            
        
        } else {
            
            table3[which(table3$station == c & table3$net == n),c(40:80)] <- NA
            
        }
    
    } # eo 2nd for loop
    
} # eo 1st for loop
# Check
summary(table3)
# Check per net
summary(table3[table3$net == "bongo",c(40:80)])
summary(table3[table3$net == "wp2",c(40:80)])
summary(table3[table3$net == "regent",c(40:80)])

### Make a map to check if sizes make sense
quartz()
ggplot() + geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
    geom_point(aes(x = x, y = y, fill = sqrt(Calanoida)), data = table3[table3$net == "wp2",], colour = "black", pch = 21, size = 2.75) + 
    scale_fill_distiller(name = "Calanoida ESD", palette = "Spectral") + coord_quickmap() + 
    theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
        panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )


### --------------------------------------------------------------------------------------------------------------------------------

### 29/05/2020: Compute N per groups and for each net
setwd(paste(WD,"/","WP2","/", sep = "")) ; dir()
nb <- read.csv("WP2numbers.csv", h = F, sep = ",")
head(nb); dim(nb)

# Does it match the names in "WP2speciesnames.csv"? 
names <- read.csv("WP2speciesnames.csv", h = T, sep = ";")
head(names); dim(names) # seems ok
colnames(nb) <- names[,1]
# Add a vector of stations
nb$station <- c(1:nrow(nb))

### Check if the counts make sense 
nb$living_Eukaryota_Opisthokonta_Holozoa_Metazoa_Arthropoda_Crustacea_Maxillopoda_Copepoda_Calanoida_Temoridae
nb$living_Eukaryota_Opisthokonta_Holozoa_Metazoa_Arthropoda_Crustacea_Maxillopoda_Copepoda_Calanoida_Calanidae
# weiirddd
nb$living_Eukaryota_Opisthokonta_Holozoa_Metazoa_Arthropoda_Crustacea_Maxillopoda_Copepoda

# OK this works, melt to summarize counts per groups of interest
melt <- melt(nb, id.vars = "station")
colnames(melt)[c(2,3)] <- c('group','n')
head(melt) ; dim(melt)
#summary(melt) ; melt[1000:1500,]
### Inform their pft/ ord/ fam
melt$pft <- NA
melt$ord <- NA
melt$fam <- NA
for(g in unique(melt$group)) {
        
        # g <- unique(melt$group)[42] ; g
        message(paste(g, sep = ""))
        pft <- unique(names[names$original_label == g,"group_pft"])
        ord <- unique(names[names$original_label == g,"group_ord"])
        fam <- unique(names[names$original_label == g,"group_fam"])
        melt[melt$group == g,"pft"] <- as.character(pft)
        melt[melt$group == g,"ord"] <- as.character(ord)
        melt[melt$group == g,"fam"] <- as.character(fam)
        
} # eo for loop - g in groups
head(melt) ; dim(melt) 
melt[melt$fam == "Calanidae","n"]
melt[melt$fam == "Temoridae","n"]
melt[melt$fam == "Oncaeidae","n"]
# Check if pfts are consistent
melt[melt$fam == "rhizaria","n"]
melt[melt$ord == "rhizaria","n"]
melt[melt$pft == "rhizaria","n"]
dim(melt[melt$pft == "rhizaria",]) ; dim(melt[melt$fam == "Oithonidae",])

# Dcast to summarize sums per station/ group
cast_pft <- dcast(melt[,c("station","pft","n")], factor(station) ~ factor(pft), fun.aggregate = sum)
head(cast_pft); dim(cast_pft)
### Looks good, remove those stations that have no counts (only zeroes)
cast_pft$tot <- rowSums(as.matrix(cast_pft[,c(2:length(cast_pft))]))
cast_pft$tot
cast_pft2 <- cast_pft[which(cast_pft$tot > 0),]
dim(cast_pft2) # cool
colnames(cast_pft2)[1] <- "station"
# Remove: 4-6, 10-11 and 15 (tot)
cast_pft2 <- cast_pft2[,c(1:3,7:9,12:14)]
# summary(cast_pft2)
colnames(cast_pft2)
colnames(cast_pft2) <- c("station","Copepoda","Chaetognatha","Gel_carn","Gel_FF","Crustaceans","Pteropoda","Rhizaria","Grazers")
# Ok, finally, melt and plot
counts_pft_wp2 <- melt(cast_pft2, id.vars = "station")
head(counts_pft_wp2)
colnames(counts_pft_wp2) <- c("station","group","n")

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_pft_wp2) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(25)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()


### Do the same with ord instead of pft ------------------------------------------------------------
cast_ord <- dcast(melt[,c("station","ord","n")], factor(station) ~ factor(ord), fun.aggregate = sum)
head(cast_ord); dim(cast_ord)
### Looks good, remove those stations that have no counts (only zeroes)
cast_ord$tot <- rowSums(as.matrix(cast_ord[,c(2:length(cast_ord))]))
cast_ord$tot
cast_ord2 <- cast_ord[which(cast_ord$tot > 0),]
dim(cast_ord2) # cool
colnames(cast_ord2)[1] <- "station"
colnames(cast_ord2)
# Remove: 8-10, 14-15, 19
cast_ord2 <- cast_ord2[,c(1:7,11:13,16:18)]
# summary(cast_pft2)
colnames(cast_ord2)
colnames(cast_ord2) <- c("station","Calanoida","Cyclopoida","Harpacticoida","Poecilostomatoida","Chaetognatha",
                        "Copepoda_unid","Gel_carn","Gel_FF","Crustaceans","Pteropoda","Rhizaria","Grazers")
# Ok, finally, melt and plot
counts_ord_wp2 <- melt(cast_ord2, id.vars = "station")
head(counts_ord_wp2)
colnames(counts_ord_wp2) <- c("station","group","n")

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_ord_wp2) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(25)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()


### Do the same with fam instead of ord ------------------------------------------------------------
cast_fam <- dcast(melt[,c("station","fam","n")], factor(station) ~ factor(fam), fun.aggregate = sum)
head(cast_fam); dim(cast_fam)
### Looks good, remove those stations that have no counts (only zeroes)
cast_fam$tot <- rowSums(as.matrix(cast_fam[,c(2:length(cast_fam))]))
cast_fam$tot
cast_fam2 <- cast_fam[which(cast_fam$tot > 0),]
dim(cast_fam2) # cool
colnames(cast_fam2)[1] <- "station"
colnames(cast_fam2)
# Remove: 11-13, 24-25, 34
cast_fam2 <- cast_fam2[,c(1:10,14:23,26:33)]
# summary(cast_pft2)
colnames(cast_fam2)
colnames(cast_fam2)[c(5,8,9,13,14,17,23,25,27)] <- c("Calanoida_unid","Chaetognatha","Copepoda_unid",
                                                    "Gel_carn","Gel_FF","Crustaceans","Pteropoda","Rhizaria","Grazers")
# Ok, finally, melt and plot
counts_fam_wp2 <- melt(cast_fam2, id.vars = "station")
head(counts_fam_wp2)
colnames(counts_fam_wp2) <- c("station","group","n")

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_fam_wp2) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(25)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()
    
### Gather all 3 countsdata frame
counts_fam_wp2[counts_fam_wp2$group == "Gel_FF","n"]
counts_ord_wp2[counts_ord_wp2$group == "Gel_FF","n"]
counts_pft_wp2[counts_pft_wp2$group == "Gel_FF","n"]
# Ok, makes sense...combine in one all the non overlapping groups
counts_wp2 <- counts_pft_wp2
commons1 <- intersect(unique(counts_pft_wp2$group), unique(counts_ord_wp2$group))
diff1 <- unique(counts_ord_wp2$group)[!(unique(counts_ord_wp2$group) %in% commons1)]

commons2 <- intersect(unique(counts_pft_wp2$group), unique(counts_fam_wp2$group))
diff2 <- unique(counts_fam_wp2$group)[!(unique(counts_fam_wp2$group) %in% commons2)]

counts_wp2 <- rbind(counts_wp2, counts_ord_wp2[counts_ord_wp2$group %in% diff1,])
counts_wp2 <- rbind(counts_wp2, counts_fam_wp2[counts_fam_wp2$group %in% diff2,])
unique(counts_wp2$group)

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_wp2) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(20)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Distribution of number of individuals per group - WP2 stations")

counts_wp2[counts_wp2$group == "Gel_FF",] # dim(counts_wp2[counts_wp2$group == "Gel_FF",])
### Per group, compute % of stations with n > 30 --> summarize with conditional count
perc <- data.frame(counts_wp2 %>% group_by(group) %>% summarize(perc = (length(unique(station[n >= 20]))/ length(unique(station)))*100 ) )
perc

### Go save counts
setwd("/Users/fabiobenedetti/Desktop/~work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/data/")
write.table(counts_wp2, file = "table_counts_groups_wp2.txt", sep = "\t")


### -----------------------------------------------------

### Same as above but for Bongo
setwd(paste(WD,"/","bongo","/", sep = "")) ; dir()
nb <- read.csv("bongonumbers.csv", h = F, sep = ",")
head(nb); dim(nb)

# Does it match the names in "WP2speciesnames.csv"? 
names <- read.csv("bongospeciesnames.csv", h = T, sep = ";")
head(names); dim(names) # seems ok
colnames(nb) <- names[,1]
# Add a vector of stations
nb$station <- c(1:nrow(nb))

# OK this works, melt to summarize counts per groups of interest
melt <- melt(nb, id.vars = "station")
colnames(melt)[c(2,3)] <- c('group','n')
head(melt) ; dim(melt)
#summary(melt) ; melt[1000:1500,]
### Inform their pft/ ord/ fam
melt$pft <- NA
melt$ord <- NA
melt$fam <- NA
for(g in unique(melt$group)) {
        
        # g <- unique(melt$group)[42] ; g
        message(paste(g, sep = ""))
        pft <- unique(names[names$original_label == g,"group_pft"])
        ord <- unique(names[names$original_label == g,"group_ord"])
        fam <- unique(names[names$original_label == g,"group_fam"])
        melt[melt$group == g,"pft"] <- as.character(pft)
        melt[melt$group == g,"ord"] <- as.character(ord)
        melt[melt$group == g,"fam"] <- as.character(fam)
        
} # eo for loop - g in groups
head(melt) ; dim(melt) 
melt[melt$fam == "Calanidae","n"]
melt[melt$fam == "Temoridae","n"]
melt[melt$fam == "Oncaeidae","n"]
# Check if pfts are consistent
melt[melt$fam == "rhizaria","n"]
melt[melt$ord == "rhizaria","n"]
melt[melt$pft == "rhizaria","n"]
# dim(melt[melt$pft == "rhizaria",]) ; dim(melt[melt$fam == "Oithonidae",])

# Dcast to summarize sums per station/ group
cast_pft <- dcast(melt[,c("station","pft","n")], factor(station) ~ factor(pft), fun.aggregate = sum)
head(cast_pft); dim(cast_pft)
### Looks good, remove those stations that have no counts (only zeroes)
cast_pft$tot <- rowSums(as.matrix(cast_pft[,c(2:length(cast_pft))]))
cast_pft$tot
cast_pft2 <- cast_pft[which(cast_pft$tot > 0),]
dim(cast_pft2) # cool
colnames(cast_pft2)[1] <- "station"
colnames(cast_pft2)
# Remove: 4-7, 11-12, 16
cast_pft2 <- cast_pft2[,c(1:3,8:10,13:15)]
# summary(cast_pft2)
colnames(cast_pft2)
colnames(cast_pft2) <- c("station","Copepoda","Chaetognatha","Gel_carn","Gel_FF","Crustaceans","Pteropoda","Rhizaria","Grazers")
# Ok, finally, melt and plot
counts_pft_bg <- melt(cast_pft2, id.vars = "station")
head(counts_pft_bg)
colnames(counts_pft_bg) <- c("station","group","n")

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_pft_bg) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(25)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()


### Do the same with ord instead of pft ------------------------------------------------------------
cast_ord <- dcast(melt[,c("station","ord","n")], factor(station) ~ factor(ord), fun.aggregate = sum)
head(cast_ord); dim(cast_ord)
### Looks good, remove those stations that have no counts (only zeroes)
cast_ord$tot <- rowSums(as.matrix(cast_ord[,c(2:length(cast_ord))]))
cast_ord$tot
cast_ord2 <- cast_ord[which(cast_ord$tot > 0),]
dim(cast_ord2) # cool
colnames(cast_ord2)[1] <- "station"
colnames(cast_ord2)
# Remove: 9-12,16-17,21
cast_ord2 <- cast_ord2[,c(1:8,13:15,18:20)]
# summary(cast_pft2)
colnames(cast_ord2)
colnames(cast_ord2) <- c("station","Calanoida","Cyclopoida","Harpacticoida","Monstrilloida","Poecilostomatoida","Chaetognatha",
                        "Copepoda_unid","Gel_carn","Gel_FF","Crustaceans","Pteropoda","Rhizaria","Grazers")
# Ok, finally, melt and plot
counts_ord_bg <- melt(cast_ord2, id.vars = "station")
head(counts_ord_bg)
colnames(counts_ord_bg) <- c("station","group","n")

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_ord_bg) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(25)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()


### Do the same with fam instead of ord ------------------------------------------------------------
cast_fam <- dcast(melt[,c("station","fam","n")], factor(station) ~ factor(fam), fun.aggregate = sum)
head(cast_fam); dim(cast_fam)
### Looks good, remove those stations that have no counts (only zeroes)
cast_fam$tot <- rowSums(as.matrix(cast_fam[,c(2:length(cast_fam))]))
cast_fam$tot
cast_fam2 <- cast_fam[which(cast_fam$tot > 0),]
dim(cast_fam2) # cool
colnames(cast_fam2)[1] <- "station"
colnames(cast_fam2)
# Remove: 13-16,29-30,39
cast_fam2 <- cast_fam2[,c(1:12,17:28,31:38)]
# summary(cast_pft2)
colnames(cast_fam2)
colnames(cast_fam2)[c(7,10,11,15,16,19,27,29,31)] <- c("Calanoida_unid","Chaetognatha","Copepoda_unid",
                                                     "Gel_carn","Gel_FF","Crustaceans","Pteropoda","Rhizaria","Grazers")
# Ok, finally, melt and plot
counts_fam_bg <- melt(cast_fam2, id.vars = "station")
head(counts_fam_bg)
colnames(counts_fam_bg) <- c("station","group","n")

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_fam_bg) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(25)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()
    
### Gather all 3 countsdata frame
counts_fam_bg[counts_fam_bg$group == "Gel_FF","n"]
counts_ord_bg[counts_ord_bg$group == "Gel_FF","n"]
counts_pft_bg[counts_pft_bg$group == "Gel_FF","n"]
# Ok, makes sense...combine in one all the non overlapping groups
counts_bg <- counts_pft_bg
commons1 <- intersect(unique(counts_pft_bg$group), unique(counts_ord_bg$group))
diff1 <- unique(counts_ord_bg$group)[!(unique(counts_ord_bg$group) %in% commons1)]

commons2 <- intersect(unique(counts_pft_bg$group), unique(counts_fam_bg$group))
diff2 <- unique(counts_fam_bg$group)[!(unique(counts_fam_bg$group) %in% commons2)]

counts_bg <- rbind(counts_bg, counts_ord_bg[counts_ord_bg$group %in% diff1,])
counts_bg <- rbind(counts_bg, counts_fam_bg[counts_fam_bg$group %in% diff2,])
unique(counts_bg$group)

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_bg) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(20)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Distribution of number of individuals per group - Bongo stations")

### Per group, compute % of stations with n > 30 --> summarize with conditional count
perc <- data.frame(counts_bg %>% group_by(group) %>% summarize(perc = (length(unique(station[n >= 20]))/ length(unique(station)))*100 ) )
perc

### Go save counts
setwd("/Users/fabiobenedetti/Desktop/~work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/data/")
write.table(counts_bg, file = "table_counts_groups_bongo.txt", sep = "\t")


### -----------------------------------------------------

### Last: do it for Régent
setwd(paste(WD,"/","regent","/", sep = "")) ; dir()
nb <- read.csv("regentnumbers.csv", h = F, sep = ",")
head(nb); dim(nb)

# Does it match the names in "WP2speciesnames.csv"? 
names <- read.csv("regentspeciesnames.csv", h = T, sep = ";")
head(names); dim(names) # seems ok
colnames(nb) <- names[,1]
# Add a vector of stations
nb$station <- c(1:nrow(nb))

# OK this works, melt to summarize counts per groups of interest
melt <- melt(nb, id.vars = "station")
colnames(melt)[c(2,3)] <- c('group','n')
head(melt) ; dim(melt)
#summary(melt) ; melt[1000:1500,]
### Inform their pft/ ord/ fam
melt$pft <- NA
melt$ord <- NA
melt$fam <- NA
for(g in unique(melt$group)) {
        
        # g <- unique(melt$group)[42] ; g
        message(paste(g, sep = ""))
        pft <- unique(names[names$original_label == g,"group_pft"])
        ord <- unique(names[names$original_label == g,"group_ord"])
        fam <- unique(names[names$original_label == g,"group_fam"])
        melt[melt$group == g,"pft"] <- as.character(pft)
        melt[melt$group == g,"ord"] <- as.character(ord)
        melt[melt$group == g,"fam"] <- as.character(fam)
        
} # eo for loop - g in groups
head(melt) ; dim(melt) 
melt[melt$fam == "Calanidae","n"]
melt[melt$fam == "Temoridae","n"]
melt[melt$fam == "Oncaeidae","n"]
# Check if pfts are consistent
melt[melt$fam == "rhizaria","n"]
melt[melt$ord == "rhizaria","n"]
melt[melt$pft == "rhizaria","n"]
# dim(melt[melt$pft == "rhizaria",]) ; dim(melt[melt$fam == "Oithonidae",])

# Dcast to summarize sums per station/ group
cast_pft <- dcast(melt[,c("station","pft","n")], factor(station) ~ factor(pft), fun.aggregate = sum)
head(cast_pft); dim(cast_pft)
### Looks good, remove those stations that have no counts (only zeroes)
cast_pft$tot <- rowSums(as.matrix(cast_pft[,c(2:length(cast_pft))]))
cast_pft$tot
cast_pft2 <- cast_pft[which(cast_pft$tot > 0),]
dim(cast_pft2) # cool
colnames(cast_pft2)[1] <- "station"
colnames(cast_pft2)
# Remove: 4-5, 9-10, 14
cast_pft2 <- cast_pft2[,c(1:3,6:8,11:13)]
# summary(cast_pft2)
colnames(cast_pft2)
colnames(cast_pft2) <- c("station","Copepoda","Chaetognatha","Gel_carn","Gel_FF","Crustaceans","Pteropoda","Rhizaria","Grazers")
# Ok, finally, melt and plot
counts_pft_reg <- melt(cast_pft2, id.vars = "station")
head(counts_pft_reg)
colnames(counts_pft_reg) <- c("station","group","n")

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_pft_reg) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(25)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()


### Do the same with ord instead of pft ------------------------------------------------------------
cast_ord <- dcast(melt[,c("station","ord","n")], factor(station) ~ factor(ord), fun.aggregate = sum)
head(cast_ord); dim(cast_ord)
### Looks good, remove those stations that have no counts (only zeroes)
cast_ord$tot <- rowSums(as.matrix(cast_ord[,c(2:length(cast_ord))]))
cast_ord$tot
cast_ord2 <- cast_ord[which(cast_ord$tot > 0),]
dim(cast_ord2) # cool
colnames(cast_ord2)[1] <- "station"
colnames(cast_ord2)
# Remove: 9-10,14-15,19
cast_ord2 <- cast_ord2[,c(1:8,11:13,16:18)]
# summary(cast_pft2)
colnames(cast_ord2)
colnames(cast_ord2) <- c("station","Calanoida","Cyclopoida","Harpacticoida","Monstrilloida","Poecilostomatoida","Chaetognatha",
                        "Copepoda_unid","Gel_carn","Gel_FF","Crustaceans","Pteropoda","Rhizaria","Grazers")
# Ok, finally, melt and plot
counts_ord_reg <- melt(cast_ord2, id.vars = "station")
head(counts_ord_reg)
colnames(counts_ord_reg) <- c("station","group","n")

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_ord_reg) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(25)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()


### Do the same with fam instead of ord ------------------------------------------------------------
cast_fam <- dcast(melt[,c("station","fam","n")], factor(station) ~ factor(fam), fun.aggregate = sum)
head(cast_fam); dim(cast_fam)
### Looks good, remove those stations that have no counts (only zeroes)
cast_fam$tot <- rowSums(as.matrix(cast_fam[,c(2:length(cast_fam))]))
cast_fam$tot
cast_fam2 <- cast_fam[which(cast_fam$tot > 0),]
dim(cast_fam2) # cool
colnames(cast_fam2)[1] <- "station"
colnames(cast_fam2)
# Remove: 12-13,25-26,34
cast_fam2 <- cast_fam2[,c(1:11,14:24,27:33)]
# summary(cast_pft2)
colnames(cast_fam2)
colnames(cast_fam2)[c(6,9,10,14,15,18,24,26,28)] <- c("Calanoida_unid","Chaetognatha","Copepoda_unid",
                                                     "Gel_carn","Gel_FF","Crustaceans","Pteropoda","Rhizaria","Grazers")
# Ok, finally, melt and plot
counts_fam_reg <- melt(cast_fam2, id.vars = "station")
head(counts_fam_reg)
colnames(counts_fam_reg) <- c("station","group","n")

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_fam_reg) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(25)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ coord_flip()
    
### Gather all 3 countsdata frame
counts_fam_reg[counts_fam_reg$group == "Gel_FF","n"]
counts_ord_reg[counts_ord_reg$group == "Gel_FF","n"]
counts_pft_reg[counts_pft_reg$group == "Gel_FF","n"]
# Ok, makes sense...combine in one all the non overlapping groups
counts_reg <- counts_pft_reg
commons1 <- intersect(unique(counts_pft_reg$group), unique(counts_ord_reg$group))
diff1 <- unique(counts_ord_reg$group)[!(unique(counts_ord_reg$group) %in% commons1)]

commons2 <- intersect(unique(counts_pft_reg$group), unique(counts_fam_reg$group))
diff2 <- unique(counts_fam_reg$group)[!(unique(counts_fam_reg$group) %in% commons2)]

counts_reg <- rbind(counts_reg, counts_ord_reg[counts_ord_reg$group %in% diff1,])
counts_reg <- rbind(counts_reg, counts_fam_reg[counts_fam_reg$group %in% diff2,])
unique(counts_reg$group)

quartz()
ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n), fill = factor(group)), data = counts_reg) +
    geom_boxplot(colour = "black") + xlab("") + ylab("N individuals") + theme_classic() + 
    scale_fill_discrete(name = "") + geom_hline(aes(yintercept = log1p(20)), linetype = "dashed") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Distribution of number of individuals per group - Regent stations")

### Per group, compute % of stations with n > 30 --> summarize with conditional count
perc <- data.frame(counts_reg %>% group_by(group) %>% summarize(perc = (length(unique(station[n >= 20]))/ length(unique(station)))*100 ) )
perc

### Go save counts
setwd("/Users/fabiobenedetti/Desktop/~work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/data/")
write.table(counts_reg, file = "table_counts_groups_regent.txt", sep = "\t")

### Remake the 3 main plots, assemble in a panel with ggarrange and send to Manoel & Fabien
counts_reg <- read.table("table_counts_groups_regent.txt", h = T, sep = "\t")
counts_bg <- read.table("table_counts_groups_bongo.txt", h = T, sep = "\t")
counts_wp2 <- read.table("table_counts_groups_wp2.txt", h = T, sep = "\t")
dim(counts_reg) ; dim(counts_bg) ; dim(counts_wp2)

plotA <- ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n) ), data = counts_bg) +
    geom_boxplot(colour = "black", fill = "grey70") + xlab("") + ylab("N individuals") + theme_classic() + 
    geom_hline(aes(yintercept = log1p(25)), linetype = "dashed", colour = "#e31a1c") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ggtitle("Distribution of number of individuals per group across Bongo stations")

plotB <- ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n) ), data = counts_wp2) +
    geom_boxplot(colour = "black", fill = "grey70") + xlab("") + ylab("N individuals") + theme_classic() + 
    geom_hline(aes(yintercept = log1p(25)), linetype = "dashed", colour = "#e31a1c") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ggtitle("Distribution of number of individuals per group across WP2 stations")

plotC <- ggplot(aes(x = fct_reorder(factor(group), log1p(n),.desc = T), y = log1p(n) ), data = counts_reg) +
    geom_boxplot(colour = "black", fill = "grey70") + xlab("") + ylab("N individuals") + theme_classic() + 
    geom_hline(aes(yintercept = log1p(25)), linetype = "dashed", colour = "#e31a1c") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ggtitle("Distribution of number of individuals per group across Regent stations")

# Arrange in panel
library("ggpubr")
#quartz()
panel <- ggarrange(plotA, plotB, plotC, labels = LETTERS[1:3], nrow = 3, ncol = 1, align = "hv")
ggsave(plot = panel, filename = "panel_distrib_N_groups_allnets.pdf", dpi = 300, height = 12, width = 8)

### And summarize % stations with N > 30
perc_bg <- data.frame(counts_bg %>% group_by(group) %>% summarize(perc = (length(unique(station[n >= 25]))/ length(unique(station)))*100 ) )
perc_wp2 <- data.frame(counts_wp2 %>% group_by(group) %>% summarize(perc = (length(unique(station[n >= 25]))/ length(unique(station)))*100 ) )
perc_rg <- data.frame(counts_reg %>% group_by(group) %>% summarize(perc = (length(unique(station[n >= 25]))/ length(unique(station)))*100 ) )

# Check by ordering
perc_bg[order(perc_bg$perc, decreasing = T),]
perc_wp2[order(perc_wp2$perc, decreasing = T),]
perc_rg[order(perc_rg$perc, decreasing = T),]


# --------------------------------------------------------------------------------------------------------------------------------

### 02/06/2020: Combine with all median sizes zoo
nets <- c("WP2","bongo","regent")
n <- "bongo"

res <- lapply(nets, function(n) {
    
        message(paste("Preparing data for net ",n,sep = ""))
        
        # Check data structure for WP2 net, should be the same for all nets
        setwd(paste(WD,"/",n,"/", sep = "")) # ; dir()
        # Load
        pft_colnames <- read.csv(paste(n,"speciesnamespft.csv", sep = ""), h = F, sep = ",") #; dim(pft_colnames)
        families_colnames <- read.csv(paste(n,"speciesnamespftfamilies.csv", sep = ""), h = F, sep = ";") #; dim(families_colnames)
        order_colnames <- read.csv(paste(n,"speciesnamespftorder.csv", sep = ""), h = F, sep = ";") #; dim(order_colnames)
        # And load the data
        esd_pft <- read.csv(paste(n,"mediansizespft.csv", sep = ""), h = F, sep = ",") #; dim(esd_pft) ; nrow(pft_colnames)
        esd_ord <- read.csv(paste(n,"mediansizespftorder.csv", sep = ""), h = F, sep = ",") #; dim(esd_ord) ; nrow(order_colnames)
        esd_fam <- read.csv(paste(n,"mediansizespftfamilies.csv", sep = ""), h = F, sep = ",") #; dim(esd_fam) ; nrow(families_colnames)
        colnames(esd_pft) <- pft_colnames[,1] ; colnames(esd_ord) <- order_colnames[,1] ; colnames(esd_fam) <- families_colnames[,1]
        # summary(esd_pft); summary(esd_ord); summary(esd_fam)
        
        ### Add the all zoo median size
        all_zoo <- read.csv(paste(n,"mediansizesallzoo.csv", sep = ""), h = F, sep = ",")
        # dim(all_zoo); summary(all_zoo)
        colnames(all_zoo) <- "Zooplankton"

        ### Remove those columns that only have NAs
        not_all_na <- function(x) any(!is.na(x))
        esd_pft <- esd_pft %>% select_if(not_all_na)
        esd_ord <- esd_ord %>% select_if(not_all_na)
        esd_fam <- esd_fam %>% select_if(not_all_na)
        # summary(esd_pft$chaetognatha); summary(esd_ord$chaetognatha); summary(esd_fam$chaetognatha)
        sizes <- esd_pft[,c("all_copepoda","chaetognatha","gelatinous_carnivorous","gelatinous_filter_feeders",
                            "large_crustracean","pteropoda","rhizaria","small_grazers")]
        # summary(sizes)
        sizes <- cbind(all_zoo,sizes)
        sizes <- cbind(sizes, esd_ord[,c(1:5)] )
        sizes <- cbind(sizes, esd_fam[,c(1:6,8:20)] )
        colnames(sizes)[c(2:14)] <- c("Copepoda","Chaetognatha","Gel_carn","Gel_FF","Crust","Pteropoda",
                                    "Rhizaria","Small_grazers","Copepoda_unid","Calanoida","Cyclopoida","Harpacticoida",
                                    "Poecilostomatoida")
        # colnames(sizes)
        if(n == "WP2") {
            n <- "wp2"
        }
        sizes$net <- n
        sizes$Station <- c(1:nrow(sizes))
        
        ### Combine with hydrodata (in abundance folder)
        setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/") 
        hydro <- read.csv("Hydro_final_TARA_08_03_20.csv", h = T, sep = ",")
        # names <- data.frame(t(read.csv("headershydro.csv", h = F, sep = ";")[1,]))
        # colnames(names) <- "name" ; names <- gsub("\\'", "", as.character(names$name))
        names <- colnames(hydro)[c(2:8,39:53,56:71)] # env vars to provide
        # Match with ESD data
        #hydro$Station <- as.factor(hydro$Station)      
        # Identify common stations
        commons <- intersect(unique(sizes$Station), unique(hydro$Station))
        # colnames(hydro)
        sizes[,names] <- NA
        # Fill with for loop
        for(c in commons){
            message(paste("Station ", c, sep = "")) # c <- 25
            sizes[which(sizes$Station == c),names] <- hydro[which(hydro$Station == c),names]
        } # eo for loop
        
        return(sizes)
    
        } # eo n
    
) # eo lapply
# Rbind, check
table.sizes <- dplyr::bind_rows(res) ; rm(res)
dim(table.sizes); str(table.sizes)
summary(table.sizes)

### Save 
write.table(table.sizes, file = "table_ESD+hydro_allnets_16_06_20.txt", sep = "\t")

# Got to abund directory
setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/data/preanalyzed_data_FL_28_05_20/abundance")
dir()
# Like size data above, use a lapply to cbind all abund data per net
nets <- c("WP2","bongo","regent")
n <- "WP2"
res <- lapply(nets, function(n) {
    
        message(paste("Preparing abundance data for net ", n, sep = ""))
        
        # Load the colnames. Have to change dir for pft because Manoela hasn't save them somehow
        setwd(paste(WD,"/abundance/", sep = ""))
        pft_names <- data.frame(t(read.csv(paste(n,"pftnames.csv", sep = ""), h = F, sep = ";")))
        colnames(pft_names) <- "name"
        pft_names <- gsub("\\'", "", as.character(pft_names$name))
        
        fam_names <- data.frame(t(read.csv(paste(n,"familiesnames.csv", sep = ""), h = F, sep = ";")[1,]))
        ord_names <- data.frame(t(read.csv(paste(n,"ordersnames.csv", sep = ""), h = F, sep = ";")[1,]))
        colnames(fam_names) <- "name"  ;  colnames(ord_names) <- "name"
        # families_colnames ; order_colnames
        fam_names <- gsub("\\'", "", as.character(fam_names$name))
        ord_names <- gsub("\\'", "", as.character(ord_names$name))
        # ord_names; fam_names; pft_names
        # length(pft_names); length(ord_names); length(fam_names)
        # Add an underscore to colnames
        pft_names <- str_replace_all(as.character(pft_names), " ", "_")
        ord_names <- str_replace_all(as.character(ord_names), " ", "_")
        fam_names <- str_replace_all(as.character(fam_names), " ", "_")
        
        # And load the data
        abund_pft <- read.csv(paste(n,"abundancespft.csv", sep = ""), h = F, sep = ",")
        abund_fam <- read.csv(paste(n,"abundancespftfamilies.csv", sep = ""), h = F, sep = ",")
        abund_ord <- read.csv(paste(n,"abundancespftorder.csv", sep = ""), h = F, sep = ",") 
        # dim(abund_pft) ; dim(abund_ord); dim(abund_fam)
        ### BEWARE: for Régent data, the 3rd col of the abundances table is an extra empty one and it majes the dimensions not fit the colnames
        if(n == "regent") {
            abund_pft <- abund_pft[,c(1:2,4:length(abund_pft))]
        }
        
        if( length(pft_names) == length(abund_pft) ) {
            message(paste("Matching colnames for PFT || ", n, sep = ""))
            colnames(abund_pft) <- pft_names
        }
        if( length(ord_names) == length(abund_ord) ) {
            message(paste("Matching colnames for Orders || ", n, sep = ""))
            colnames(abund_ord) <- ord_names
        }
        if( length(fam_names) == length(abund_fam) ) {
            message(paste("Matching colnames for Families || ", n, sep = ""))
            colnames(abund_fam) <- fam_names
        }
        # summary(abund_pft) ; summary(abund_ord) ; summary(abund_fam)
        
        ### Remove those columns that only have NAs
        not_all_na <- function(x) any(!is.na(x))
        abund_pft <- abund_pft %>% select_if(not_all_na)
        abund_ord <- abund_ord %>% select_if(not_all_na)
        abund_fam <- abund_fam %>% select_if(not_all_na)
        
        ### Cbind all columns of interest
        #colnames(abund_pft)
        abund <- abund_pft[,c("all_copepoda","chaetognatha","gelatinous_carnivorous","filter_feeders",
                            "large_crustracean","pteropoda","rhizaria","small_grazers")]
        # summary(sizes)
        abund <- cbind(abund, abund_ord[,c(1:5)] )
        abund <- cbind(abund, abund_fam[,c(1:6,8:20)] )
        # dim(abund)
        # Rename
        colnames(abund)[c(1:13)] <- c("Copepoda","Chaetognatha","Gel_carn","Gel_FF","Crust","Pteropoda","Rhizaria",
                        "Small_grazers","Copepoda_unid","Calanoida","Cyclopoida","Harpacticoida","Poecilostomatoida")
        
        # Quickly check if the sum of Copepoda_unid+Calanoida+Cyclopoida+Harpacticoida+Poecilostomatoida = Copepoda
        #summary(rowSums(as.matrix(abund[,c("Copepoda_unid","Calanoida","Cyclopoida","Harpacticoida","Poecilostomatoida")])))                
        #summary(abund$Copepoda)
        #plot(abund$Copepoda,rowSums(as.matrix(abund[,c("Copepoda_unid","Calanoida","Cyclopoida","Harpacticoida","Poecilostomatoida")])) )
        ### GOOD :)
        
        # Derive total zoo abundance
        abund$Zooplankton <- rowSums(as.matrix(abund[,c(1:8)]))
        # summary(abund$Zooplankton)
        
        ### Provide net
        if(n == "WP2") {
            n <- "wp2"
        }
        abund$net <- n
        
        ### Load vector of stations that shoukd match dimensions of 'abund'
        stations <- read.csv(paste(n,"stations.csv", sep= "") , h = F, sep = ";")
        colnames(stations) <- "station"
        if( nrow(abund) == nrow(stations) ) {
            message(paste("Matching number of rows for stations || ", n, sep = ""))
            abund$Stations <- factor(stations$station)   
        }
        
        ### And now combine with the hydro data from March (08_03_20)
        setwd("/Users/fabiobenedetti/Desktop/work/PostDocs/TARA/Brandao&al._MS#1/REVISIONS/") 
        hydro <- read.csv("Hydro_final_TARA_08_03_20.csv", h = T, sep = ",")
        # names <- data.frame(t(read.csv("headershydro.csv", h = F, sep = ";")[1,]))
        # colnames(names) <- "name" ; names <- gsub("\\'", "", as.character(names$name))
        names <- colnames(hydro)[c(2:8,39:53,56:71)] # env vars to provide
        # Match with ESD data
        #hydro$Station <- as.factor(hydro$Station)      
        # Identify common stations
        commons <- intersect(unique(abund$Station), unique(hydro$Station))
        # colnames(hydro)
        abund[,names] <- NA
        # Fill with for loop
        for(c in commons){
            message(paste("Station ", c, sep = "")) # c <- 25
            abund[which(abund$Station == c),names] <- hydro[which(hydro$Station == c),names]
        } # eo for loop
        # summary(abund)
        
        return(abund)
    
        } # eo n
    
) # eo lapply
# Rbind
table.abund <- dplyr::bind_rows(res) ; rm(res)
dim(table.abund); str(table.abund) ; colnames(table.abund)
summary(table.abund)

### Save 
write.table(table.abund, file = "table_abund+hydro_allnets_16_06_20.txt", sep = "\t")


summary()

### But older hydro have less NAs than before...
### Just go back to the older version I guess...i don't have time for this


### OK, close and move to Script 4

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
