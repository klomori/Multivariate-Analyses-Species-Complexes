#######################################
#### Other Rockfish Combined Data  ####
#### Data prep & analyses          ####
#### October 2019: Kristen Omori   ####
#######################################

#### Load libraries ####
library(tm)
library(remotes)
library(dplyr)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(grid)
library(ggvegan)
library(vegan)
library(cowplot)
library(pvclust)

`%notin%` <- function(x,y) !(x %in% y) 
'%nin%' <- Negate('%in%')
pos.count.fn <- function(x) {sum(x >0, na.rm=T)}

#### Main directory ####
rf.dir.sum <- ("G:/My Drive/Rockfish Background & Data/Rockfish R files & summary/")

#### Read in Data Reshaping and MV functions ####

source("G:/My Drive/Rockfish Background & Data/Rockfish R files & summary/OR_MV_functions_new.R")

#################################
#### Trawl & Longline survey ####
#################################

#### Read in areas & species ####
FMP <- c("GOA")  # otherwise can use AI or BSAI
year.condition.survey <- 1 # 1= all, 2= year..., 3= years...
start.year.ll.survey <- 1995

#species.code <- c(30330, 30340, 30360, 30200)   # NEED TO CHANGE
spp.name <- c("rockfish, black", "rockfish, blackgill", "rockfish, blue", "rockfish, bocaccio",
              "rockfish, canary", "rockfish, chilipepper", "rockfish, china", "rockfish, copper",
              "rockfish, darkblotched", "rockfish, dark", "rockfish, dusky","rockfish, greenstripe", 
              "rockfish, harlequin", "rockfish, northern", "rockfish, pygmy", "rockfish, quillback", 
              "rockfish, redbanded", "rockfish, redstripe", "rockfish, rosethorn", 
              "rockfish, sharpchin", "rockfish, silvergray", "rockfish, splitnose", 
              "rockfish, stripetail", "rockfish, tiger", "rockfish, vermilion", "rockfish, widow", 
              "rockfish, yelloweye", "rockfish, yellowmouth", "rockfish, yellowtail")
code <- c(30330, 30340, 30360, 30400, 30410, 30260, 30370, 30120, 30170, 30151, 30152, 
          30200, 30535, 30420, 30550, 30320, 30475, 30430, 30270, 30560, 30100, 30190, 
          30490, 30380, 30350, 30220, 30470, 30600, 30240)
OR.grp <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
all.species <- data.frame("spp.name" = spp.name, "code" = code, "OR.grp"= OR.grp)

species.code <- all.species[all.species$OR.grp == 1, "code"]   ## can change for just species wanted
all.species <- all.species[all.species$OR.grp == 1, ]  ## select only species with OR
all.species$Lcode <- LETTERS[1:dim(all.species)[1]]

old.mngmt <- c("Slope", "Slope", "Demersal", "Slope", "Demersal", "Demersal", "Slope",
               "Slope", "Slope", "Slope", "Slope", "Demersal", "Slope", "Slope", "Demersal",
               "Slope", "Slope", "Slope", "Slope", "Demersal", "Slope", "Slope", "Demersal",
               "Slope", "Slope")

all.species$name <- removeWords(as.character(all.species$spp.name), "rockfish, ")

#### Read in Cleaned Data Files ####

source("G:/My Drive/Rockfish Background & Data/Rockfish R files & summary/OR_survey_datainput_cleaning_new.R")

head(trawlsurv.haul.expand.dat,3)
head(llsurv.haul.expand.dat,3)

#########################
#### Fisheries Data  ####
#########################

#### Read in data (cleaned and ready to be analyzed) ####

source("G:/My Drive/Rockfish Background & Data/Rockfish R files & summary/OR_fisheries_datainput_cleaning_new.R")
setwd(rf.dir.sum)


#####################################
#### Creating new single dataset ####
#####################################

## Adding Unique Unit Column ##

## Unit = spatial- temporal- gear: NMFS management unit - year - month - gear
head(trawlsurv.haul.expand.dat)
head(llsurv.haul.expand.dat)
head(fish.dat)

## Create unique 'Unit' coluumn ##

trawl.unit.dat <- factors.dat.fn(trawlsurv.haul.expand.dat) # starts 1984
ll.unit.dat <- factors.dat.fn(llsurv.haul.expand.dat) # 1995
fish.unit.dat <- factors.dat.fn(fish.dat) # 2010

## Create unique 'Unit' with season (not month)

## Create Unit dataframe

trawl.prop.dat <- PA_prop_datprep_fn(trawl.unit.dat)
ll.prop.dat <- PA_prop_datprep_fn(ll.unit.dat)
fish.prop.dat <- PA_prop_datprep_fn(fish.unit.dat)

## Combine all datasets into one ##

trawl.prop.dat <- rownames_to_column(trawl.prop.dat, "Unit")
ll.prop.dat <- rownames_to_column(ll.prop.dat, "Unit")
fish.prop.dat <- rownames_to_column(fish.prop.dat, "Unit")

all.prop.dat <- bind_rows(trawl.prop.dat, ll.prop.dat, fish.prop.dat) %>% 
  remove_rownames() %>%
  column_to_rownames(var="Unit") %>%
  replace(is.na(.),0)

## Subset for rare species

all.prop.sub.dat <- subset_norare_fn(all.prop.dat)

all.prop.sub.dat <- all.prop.sub.dat %>% rename(canary_d= canary, quillback_d= quillback, rosethorn_d= rosethorn,
                                               tiger_d =tiger, yelloweye_d= yelloweye, china_d= china)

names(all.prop.dat); length(names(all.prop.dat))
names(all.prop.sub.dat); length(names(all.prop.sub.dat))

## Season string
month.string <- data.frame("Month"= month.name, 
                           "Season"= c(rep("Winter", 3), rep("Spring", 3), rep("Summer", 3), rep("Fall", 3)))

#### Read in ENSO and PDO data files ####
ENSO.dat <- read.csv("ENSO_index_years_months.csv") %>% gather(key= "Month", value= "ENSO_index", -Year)
PDO.dat <- read.csv("PDO_index_years_months.csv") %>% gather(key= "Month", value= "PDO_index", -Year)

set.seed(51584)
##########################################
#### Datasets for analyses Q & R mode ####
##########################################

## Q mode => sites= rows, species= cols; R mode => species= rows, sites= cols

q.sub.dat <- all.prop.sub.dat ## x=sites, y=species
r.sub.dat <- t(all.prop.sub.dat) ## x=species, y=sites

q.root.dat <- rootroot_scale_fn(all.prop.sub.dat) # root-root transformed
r.root.dat <- t(rootroot_scale_fn(all.prop.sub.dat)) # root-root transformed

## Distance matrices
q.chord.mat <- dist.ldc(q.sub.dat, method= "chord")  # chord distance matrix for proportions
r.chord.mat <- dist.ldc(r.sub.dat, method= "chord")  # chord distance matrix for proportions
q.chisq.mat <- dist.ldc(q.root.dat, method= "chisquare") # chi_sq distance matrix
r.chisq.mat <- dist.ldc(r.root.dat, method= "chisquare") # chi_sq distance matrix

r.chord.mat2 <- dist.ldc(r.root.dat, method= "chord")
#q.pa.dat <- Pres_Abs_fn(q.sub.dat) # Makes presence-absence matrix
#q.sorensen.mat <- Sorensen_PA_mat_fn(q.pa.dat) # sorensen/ bray-curtis dissimilarity mat
#r.pa.dat <- Pres_Abs_fn(r.sub.dat) # Makes presence-absence matrix
#r.sorensen.mat <- Sorensen_PA_mat_fn(r.pa.dat) # sorensen/ bray-curtis dissimilarity mat


#######################################
#### Wards Hierarchical Clustering ####
#######################################

kclust.decision <- function(trans.dat, dist.dat, k.clust, max.n) {
  wards.results <- hclust(dist.dat, method= "ward.D2")
  print(fviz_nbclust(trans.dat, FUN= hcut, method="silhouette", k.max= max.n))
  tree.cut.temp <- cutree(wards.results, k= k.clust)
  table(tree.cut.temp)
}

wards.graph.fn <- function(trans.dat, dist.dat, k.clust, mode, x.temp, dat.set) {
  wards.results <- hclust(dist.dat, method= "ward.D2")
  tree.cut.temp <- cutree(wards.results, k= k.clust)

  plot(wards.results, hang= 0.25, xlab= x.temp, main= paste("Wards: ",mode,"- ", dat.set))
  rect.hclust(wards.results, k= k.clust, border= c(2:(k.clust+1)))
}

wards.unit.results.fn <- function(trans.dat, dist.dat, k.clust) {
  wards.results <- hclust(dist.dat, method= "ward.D2")
  tree.cut.temp <- cutree(wards.results, k= k.clust)
  
  dat.temp.new <- trans.dat %>% mutate(wards.cluster= tree.cut.temp)
  cluster.results.temp <- as.data.frame(tree.cut.temp) %>%
    rownames_to_column("Unit") %>%
    separate("Unit", c("Year", "Month", "Stratum", "Gear"), remove=FALSE) %>%
    rename(ClusterNum= tree.cut.temp)
}

wards.sp.results.fn <- function(trans.dat, dist.dat, k.clust) {
  wards.results <- hclust(dist.dat, method= "ward.D2")
  tree.cut.temp <- cutree(wards.results, k= k.clust) %>% as.data.frame() %>%
    rownames_to_column("names") %>% rename(ClustNum= ".")
}

### Q mode: 
## q.sub.dat & q.chord.mat; q.root.dat & q.chisq.mat
kclust.decision(trans.dat= q.sub.dat, dist.dat= q.chord.mat, k.clust=2, max.n=50)
wards.graph.fn(dist.dat= q.chord.mat, k.clust= 2, mode="Q-mode", dat.set="Chord", x.temp="Unit")
q.chord.results <- wards.unit.results.fn(trans.dat= q.sub.dat, dist.dat= q.chord.mat, k.clust=2)
head(q.chord.results)

kclust.decision(trans.dat= q.root.dat, dist.dat= q.chisq.mat, k.clust=7, max.n=50)
wards.graph.fn(dist.dat= q.chisq.mat, k.clust= 7, mode="Q-mode", dat.set="Chisq", x.temp="Unit")
q.chisq.results <- wards.unit.results.fn(trans.dat= q.root.dat, dist.dat= q.chisq.mat, k.clust=7)


#all.prop.sub.dat[which(row.names(all.prop.sub.dat) %in% c("2012_December_610_HAL")),]

### R mode:
## r.sub.dat & r.chord.mat; r.root.dat & r.chisq.mat
kclust.decision(trans.dat= r.sub.dat, dist.dat= r.chord.mat, k.clust=3, max.n=12)
wards.graph.fn(dist.dat= r.chord.mat, k.clust= 3, mode="R-mode-Sub", dat.set="Chord", x.temp="Species")
r.chord.results <- wards.sp.results.fn(trans.dat= r.sub.dat, dist.dat= r.chord.mat, k.clust=2)
r.chord.results

kclust.decision(trans.dat= r.root.dat, dist.dat= r.chisq.mat, k.clust=3, max.n=12)
wards.graph.fn(dist.dat= r.chisq.mat, k.clust= 3, mode="R-mode", dat.set="Chisq", x.temp="Species")
r.chisq.results <- wards.sp.results.fn(trans.dat= r.root.dat, dist.dat= r.chisq.mat, k.clust=2)
r.chisq.results

kclust.decision(trans.dat= r.root.dat, dist.dat= r.chord.mat2, k.clust=3, max.n=12)
wards.graph.fn(dist.dat= r.chord.mat2, k.clust= 3, mode="R-mode", dat.set="Chord-Root", x.temp="Species")
r.chord2.results <- wards.sp.results.fn(trans.dat= r.root.dat, dist.dat= r.chord.mat2, k.clust=2)
r.chord2.results

#source("G:/My Drive/Rockfish Background & Data/Rockfish R files & summary/Wards_bootstrapping.R")

#r.fit <- pvclust(t(r.sub.dat), method.hclust = "ward.D2", method.dist= "correlation")
#r.fit2 <- pvclust.ko(r.sub.dat, method.hclust = "ward.D2", method.dist= "chord", use.cor="pairwise.complete.obs")
#r.fit2 <- pvclust(r.sub.dat, method.hclust = "ward.D2", method.dist= "euclidean", use.cor="pairwise.complete.obs") #allobs

#plot(r.fit2)
#pvrect(r.fit, alpha= .95)
library(fpc)
r.fit.temp <- clusterboot(r.chord.mat, distances=TRUE, clustermethod=disthclustCBI, k=3, 
                          cut="number", method="ward.D2", B=1000, showplots=FALSE)
#print(r.fit.temp)
str(r.fit.temp)
plot(r.fit.temp$result$result)
k.clust.temp <- r.fit.temp$result$nc
rect.hclust(r.fit.temp$result$result, k= k.clust.temp, border= c(2:(2+1)))

# groups:
r.Jaccard.bootstrap <- data.frame(cluster= c(1:k.clust.temp) ,bootmean= round(r.fit.temp$bootmean,3))

r.fit.groups <- data.frame(cluster=r.fit.temp$partition) %>% rownames_to_column("species") %>%
  arrange(cluster) %>% full_join(r.Jaccard.bootstrap)
r.fit.groups

saveRDS(r.fit.temp, file= paste0("G:\\My Drive\\MV Form Sp Complex Project\\Results\\Graphs&Tables\\Wards_Rmode_units.rds"))

#savehistory(file="MV_analysis_allData_PVCLUST")
#q.pa.wards.results <- hclust(q.sorensen.mat, method= "average")
#fviz_nbclust(q.pa.dat, FUN= hcut, method="silhouette", k.max=50)
#num.clust <- 3
#sub.grp.q.pa <- cutree(q.pa.wards.results, k=num.clust) ; table(sub.grp.q.pa)
#plot(q.pa.wards.results, hang= 0.25, xlab= "Units", 
#     main= paste("Wards: Q-mode, P-A sub"))
#rect.hclust(q.pa.wards.results, k= num.clust, border= c(2:(num.clust+1)))
s
#############################
#### Kmediods Clustering ####
#############################

# Datasets: q.pa.dat, r.pa.dat, q.sub.dat, r.sub.dat
### Q mode: 
## q.sub.dat & q.chord.mat; q.root.dat & q.chisq.mat

### R mode:
## r.sub.dat & r.chord.mat; r.root.dat & r.chisq.mat
## Don't make a function... need to inspect the # suggested clusters to decide

## Q mode: q.sub.dat, q.root.dat, 

dat.temp <- q.sub.dat  # q.sub.dat, q.root.dat, r.sub.dat, r.root.dat  
dat.name1 <- "Q-mode, " # "Q-mode, " , "R-mode, "
dat.name2 <- "Prop sub" # "Prop sub", "Root sub"
dim(dat.temp)

#sil.results <- fviz_nbclust(dat.temp, kmeans, method='silhouette', k.max=45) # at 50 for Q
sil.results <- fviz_nbclust(dat.temp, cluster::pam, method='silhouette', k.max=15) # at 50 for Q
num.clust <- as.numeric(sil.results$data$clusters[which.max(sil.results$data$y)]) # either 3 or 18
num.clust <- 2
#write.csv(sil.results$data, "Kmediods_sil_width.csv")
kmeans.results <- kmeans(dat.temp, centers = num.clust, nstart= 10)
fig.print <- fviz_cluster(kmeans.results, data= dat.temp, ellipse.type= "convex", 
                          palette="jco", ggtheme = theme_minimal(),
                          main = paste("Kmeans:", dat.name1, dat.name2),
                          #geom=c("text"))
                          geom=c("point"))
### Kmediods

# Calculates Euclidean distance in function (don't need to standardize, already done)
kmediods.results <- pam(dat.temp, k= num.clust, metric= "euclidean", stand= FALSE)

kmethod.results <- kmediods.results$clustering  #kmeans.results$clusters
kmethod.results.all <- kmediods.results  #kmeans.results

temp.cluster.results <- as.data.frame(kmethod.results) %>% 
  rownames_to_column("Unit") %>%
  separate("Unit", c("Year", "Month", "Stratum", "Gear"), remove=FALSE) %>%
  rename(ClusterNum= "kmethod.results") %>% 
  full_join(month.string) %>% 
  as.data.frame()

fig.print <- fviz_cluster(kmethod.results.all, data= dat.temp, ellipse.type= "convex", 
                          palette="jco", ggtheme = theme_minimal(),
                          main = paste("Kmeans:", dat.name1, dat.name2),
                          #geom=c("text"))
                          geom=c("point"))

head(temp.cluster.results)
                
fig.print

kmeds.q.prop.dat2 <- temp.cluster.results
kmeds.q.prop.fig2 <- fig.print
kmeds.q.prop.dat6 <- temp.cluster.results
kmeds.q.prop.fig6 <- fig.print
kmeds.q.prop.dat8 <- temp.cluster.results
kmeds.q.prop.fig8 <- fig.print

kmeds.q.root.dat3 <- temp.cluster.results
kmeds.q.root.fig3 <- fig.print
kmeds.q.root.dat5 <- temp.cluster.results
kmeds.q.root.fig5 <- fig.print

write.csv(kmeans.q.prop.dat2, "Results_Q_kmeds_prop_dat2.csv", row.names= F)
write.csv(kmeans.q.prop.dat6, "Results_Q_kmeds_prop_dat6.csv", row.names= F)
write.csv(kmeans.q.prop.dat8, "Results_Q_kmeds_prop_dat8.csv", row.names= F)
write.csv(kmeans.q.root.dat3, "Results_Q_kmedss_root_dat3.csv", row.names= F)
write.csv(kmeans.q.root.dat5, "Results_Q_kmeds_root_dat5.csv", row.names= F)

kmeds.q.prop.dat2 <- temp.cluster.results
kmeds.q.prop.fig2 <- fig.print

kmeds.q.prop.dat5 <- temp.cluster.results
kmeds.q.prop.fig5 <- fig.print

write.csv(kmeds.q.prop.dat2, "Results_Q_kmeds_prop_dat2.csv", row.names= F)
write.csv(kmeds.q.prop.dat5, "Results_Q_kmeds_prop_dat5.csv", row.names= F)

## Saving Kmeans cluster results
# PA using k=3 clusters (could also be 18+ look at graph) (note did not replace kmeans with kmediods name, but uses kmediods)
kmeds.q.prop.dat <- kmeds.q.prop.dat5
kmeds.q.prop.fig <- kmeds.q.prop.fig5
num.temp.prop <- 5
temp1 <- kmeds.q.prop.dat %>% group_by(ClusterNum) %>% count(Year, name= "TotYear") %>% 
  spread(ClusterNum, TotYear) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp2 <- kmeds.q.prop.dat %>% group_by(ClusterNum) %>% count(Month, name= "TotMonth") %>% 
  spread(ClusterNum, TotMonth) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp3 <- kmeds.q.prop.dat %>% group_by(ClusterNum) %>% count(Stratum, name= "TotStratum") %>% 
  spread(ClusterNum, TotStratum) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp4 <- kmeds.q.prop.dat %>% group_by(ClusterNum) %>% count(Gear, name= "TotGear") %>% 
  spread(ClusterNum, TotGear) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp5 <- kmeds.q.prop.dat %>% group_by(ClusterNum) %>% count(Season, name= "TotSeason") %>% 
  spread(ClusterNum, TotSeason) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp6 <- kmeds.q.prop.dat %>% unite("Strat_Gear",c(Stratum, Gear),sep="_", remove=F) %>%
  group_by(ClusterNum) %>% count(Strat_Gear, name= "TotStratGear") %>% 
  spread(ClusterNum, TotStratGear) %>% replace(is.na(.), 0 ) %>% as.data.frame()
kmeds.q.prop.table <- bind_rows(temp1, temp2, temp3, temp4, temp5, temp6)

write.csv(kmeds.q.prop.table, paste0("Results_Q_kmediods_prop_table",num.temp.prop,".csv"), row.names = F)

## Root

kmeds.q.root.dat <- kmeds.q.root.dat5
kmeds.q.root.fig <- kmeds.q.root.fig5
num.temp.root <- 5
temp1 <- kmeds.q.root.dat %>% group_by(ClusterNum) %>% count(Year, name= "TotYear") %>% 
  spread(ClusterNum, TotYear) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp2 <- kmeds.q.root.dat %>% group_by(ClusterNum) %>% count(Month, name= "TotMonth") %>% 
  spread(ClusterNum, TotMonth) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp3 <- kmeds.q.root.dat %>% group_by(ClusterNum) %>% count(Stratum, name= "TotStratum") %>% 
  spread(ClusterNum, TotStratum) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp4 <- kmeds.q.root.dat %>% group_by(ClusterNum) %>% count(Gear, name= "TotGear") %>% 
  spread(ClusterNum, TotGear) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp5 <- kmeds.q.root.dat %>% group_by(ClusterNum) %>% count(Season, name= "TotSeason") %>% 
  spread(ClusterNum, TotSeason) %>% replace(is.na(.), 0 ) %>% as.data.frame()
temp6 <- kmeds.q.root.dat %>% unite("Strat_Gear",c(Stratum, Gear),sep="_", remove=F) %>%
  group_by(ClusterNum) %>% count(Strat_Gear, name= "TotStratGear") %>% 
  spread(ClusterNum, TotStratGear) %>% replace(is.na(.), 0 ) %>% as.data.frame()
kmeds.q.root.table <- bind_rows(temp1, temp2, temp3, temp4, temp5, temp6)
write.csv(kmeds.q.root.table, paste0("Results_Q_kmeds_root_table",num.temp.root,".csv"), row.names = F)

## Examine gear, month, year in clusters

kmediods.q.prop.dat2 <- read.csv("Results_Q_kmediods_prop_dat2.csv")
kmediods.q.prop.dat5 <- read.csv("Results_Q_kmediods_prop_dat5.csv")

kmediods.q.prop.table2 <- read.csv("Results_Q_kmediods_prop_table2.csv")
kmediods.q.prop.table5 <- read.csv("Results_Q_kmediods_prop_table5.csv")



kmeds.graph.table.fn <- function(dataset) {
  dat.temp <- dataset
  dat.temp.save <- list(NA)
  
  temp1 <- dat.temp %>% select(-c("Year", "Month", "Stratum", "Season", "Strat_Gear")) %>% 
    rename_at( vars(starts_with("X")), list(~str_replace(., "X", "Cluster_"))) %>% 
    drop_na() %>% gather(key=ClusterNum, value="value", factor_key=T, -Gear)
  temp2 <- dat.temp %>% select(-c("Year", "Month", "Gear", "Season", "Strat_Gear")) %>% 
    rename_at( vars(starts_with("X")), list(~str_replace(., "X", "Cluster_"))) %>% 
    drop_na() %>% gather(key=ClusterNum, value="value", factor_key=T, -Stratum)
  dat.temp.save[[1]] <- temp1
  dat.temp.save[[2]] <- temp2
  return(dat.temp.save)
  
}


kmediods.prop.temp2.list <- kmeds.graph.table.fn(dataset= kmediods.q.prop.table2)
kmediods.prop.temp5.list <- kmeds.graph.table.fn(dataset= kmediods.q.prop.table5)

windows(record=T)

graph.output.dir <- "G:\\My Drive\\MV Form Sp Complex Project\\Results\\Graphs&Tables\\"

kmeds.graph.fn <- function(data.list, title.temp, title.temp2, ncol.temp) {
  dat.gear <- data.list[[1]]
  dat.stratum <- data.list[[2]]
  
  graph.gear <- ggplot(dat.gear, aes(x= Gear, y = value)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Gear") + 
    ylab("# of Gear Types in Cluster") +
    theme_bw() +
    facet_wrap(.~ClusterNum, scales= "fixed", ncol= ncol.temp) +
    ggtitle(paste0("Q Kmediods Results- # Gears in cluster: ", title.temp)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle= 90, vjust = 0.5))
  
  graph.stratum <- ggplot(dat.stratum, aes(x= Stratum, y = value)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Stratum") + 
    ylab("# of Stratum in Cluster") +
    theme_bw() +
    facet_wrap(.~ClusterNum, scales= "fixed", ncol= ncol.temp) +
    ggtitle(paste0("Q Kmediods Results- # Stratum in cluster: ", title.temp)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle= 90, vjust = .5))
  
  plot.save <- plot_grid(graph.gear, graph.stratum, labels="AUTO", nrow=2)
  ggsave(paste0(graph.output.dir,"Q-mode_Kmediods_units_", title.temp2,".jpeg"), plot.save)
}

#kmeds.graph.fn(data.list= kmeds.graph.fn.root.temp3.list, title.temp= "Root", title.temp2= "Root_3", ncol.temp=3)
#kmeds.graph.fn(data.list= kmeds.graph.fn.root.temp5.list, title.temp= "Root", title.temp2= "Root_5", ncol.temp=5)
#kmeds.graph.fn(data.list= kmeds.graph.fn.prop.temp2.list, title.temp= "Prop", title.temp2= "Prop_2", ncol.temp=2)
#kmeds.graph.fn(data.list= kmeds.graph.fn.prop.temp6.list, title.temp= "Prop", title.temp2= "Prop_6", ncol.temp=3)
#kmeds.graph.fn(data.list= kmeds.graph.fn.prop.temp8.list, title.temp= "Prop", title.temp2= "Prop_8", ncol.temp=4)

kmeds.graph.fn(data.list= kmediods.prop.temp2.list, title.temp= "Prop", title.temp2= "Prop_2", ncol.temp=2)
kmeds.graph.fn(data.list= kmediods.prop.temp5.list, title.temp= "Prop", title.temp2= "Prop_5", ncol.temp=5)

## kmediods by % composition

d2 <- d %>% 
  group_by(groupchange,Symscore3) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

kmeds.perc.graph.fn <- function(data.list, title.temp, title.temp2, ncol.temp) {
  dat.gear <- data.list[[1]]
  dat.stratum <- data.list[[2]]
  
  dat.gear.temp <- dat.gear %>% group_by(ClusterNum) %>% summarise(tot_counts= sum(value))
  dat.gear.temp.merge <- dat.gear %>% full_join(dat.gear.temp) %>% mutate(perc= round(value/tot_counts,3))
  
  dat.stratum.temp <- dat.stratum %>% group_by(ClusterNum) %>% summarise(tot_counts= sum(value))
  dat.stratum.temp.merge <- dat.stratum %>% full_join(dat.stratum.temp) %>% mutate(perc= round(value/tot_counts,3))
  
  graph.gear <- ggplot(dat.gear.temp.merge, aes(x= Gear, y = perc)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Gear") + 
    ylab("% Gear Types in Cluster") +
    theme_bw() +
    facet_wrap(.~ClusterNum, scales= "fixed", ncol= ncol.temp) +
    ggtitle(paste0("Q Kmediods Results- % Gears in cluster: ", title.temp)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle= 90, vjust = 0.5))
  
  graph.stratum <- ggplot(dat.stratum.temp.merge, aes(x= Stratum, y = perc)) +
    geom_bar(stat = "identity", position = "dodge") +
    xlab("Stratum") + 
    ylab("% Stratum in Cluster") +
    theme_bw() +
    facet_wrap(.~ClusterNum, scales= "fixed", ncol= ncol.temp) +
    ggtitle(paste0("Q Kmediods Results- % Stratum in cluster: ", title.temp)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle= 90, vjust = .5))
  
  plot.save <- plot_grid(graph.gear, graph.stratum, labels="AUTO", nrow=2)
  ggsave(paste0(graph.output.dir,"Q-mode_Kmediods_Perc_", title.temp2,".jpeg"), plot.save)
}


kmeds.perc.graph.fn(data.list= kmediods.prop.temp2.list, title.temp= "Prop", title.temp2= "Prop_2", ncol.temp=2)
kmeds.perc.graph.fn(data.list= kmediods.prop.temp5.list, title.temp= "Prop", title.temp2= "Prop_5", ncol.temp=4)

##### Other attempts with season and adding ENSO and PDO
dat.enter.graph <-  kmeds.root.temp # kmeds.root.temp kmeds.prop.temp
dat.title <- "Root"

ggplot(dat.enter.graph, aes(x= Season, y = value)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Season") + 
  ylab("Number of Seasons in Cluster") +
  theme_bw() +
  facet_wrap(.~ClusterNum, scales= "free") +
  ggtitle(paste("Q Kmeans Results- # Seasons in each cluster: ", dat.title)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle= 45, vjust = .5))

#### ADDING ENSO AND PDO TEMP ####
index.dat <- kmeds.q.prop.dat5 %>% left_join(PDO.dat) %>% left_join(ENSO.dat)

dat.temp.kmeds <- index.dat %>% select(ClusterNum, PDO_index, ENSO_index) %>%
  mutate(ClusterNum= as.factor(ClusterNum))

mean.fn <- function(x) {mean(x, na.rm=T)}

ggplot(dat.temp.kmeds, aes(x= ClusterNum, y= PDO_index)) +
  geom_violin(trim=FALSE) +
  stat_summary(fun.y= mean.fn, geom= "point", colour="red") +
  xlab("Cluster") +
  ylab("ENSO_index") +
  theme_bw() +
  ggtitle("Q-mode: Kmeans") +
  theme(plot.title = element_text(hjust = 0.5))

###############################################################
######## Examine species composition in each cluster   ########

pos.count.fn <- function(x) {sum(x >0, na.rm=T)}
# Q-prop

kmeds.q.sp.prop.dat.fn <- function(dataset) {
  dat.temp <- dataset
  
  dat.temp.counts <- q.sub.dat %>% rownames_to_column(var= "Unit") %>%
    full_join(subset(dat.temp, select= c("Unit", "ClusterNum"))) %>%
    group_by(ClusterNum) %>% summarise_all(pos.count.fn) %>%
    column_to_rownames(var= "ClusterNum") %>% as.data.frame()
  dat.temp.percent <- dat.temp.counts/dat.temp.counts$Unit
  return(dat.temp.percent)
  
}

tot.units.function <- function(dataset) {
  dat.temp <- dataset
  
  dat.temp.counts <- q.sub.dat %>% rownames_to_column(var= "Unit") %>%
    full_join(subset(dat.temp, select= c("Unit", "ClusterNum"))) %>%
    group_by(ClusterNum) %>% summarise_all(pos.count.fn) %>%
    column_to_rownames(var= "ClusterNum") %>% as.data.frame()
  dat.temp.tot.units <- dat.temp.counts$Unit
  return(dat.temp.tot.units)
}



kmediods.sp.prop.2 <- kmeds.q.sp.prop.dat.fn(dataset= kmediods.q.prop.dat2)
kmediods.sp.prop.5 <- kmeds.q.sp.prop.dat.fn(dataset= kmediods.q.prop.dat5)

tot.units.sp.prop.2b <- tot.units.function(dataset= kmediods.q.prop.dat2)
tot.units.sp.prop.5b <- tot.units.function(dataset= kmediods.q.prop.dat5)


write.csv(kmediods.sp.prop.2, "Results_Q_kmediods_proppercent_table2.csv", row.names = T)
write.csv(kmediods.sp.prop.5, "Results_Q_kmediods_proppercent_table5.csv", row.names = T)


species.proportion.graph.fn <- function(dataset, tot.units, title.temp) {
  dat.enter.temp <- dataset
  dat.enter <- dat.enter.temp %>% select(-Unit) %>% 
    rownames_to_column("ClusterNum") %>% gather(key= "species", value="value", factor_key=T, -ClusterNum) %>%
    mutate(ClusterNum= as.factor(ClusterNum))
  
 var.names.temp <-  paste("Cluster: ",unique(dat.enter$ClusterNum), "_(",tot.units, ")", sep="")
 names(var.names.temp) <- c(unique(dat.enter$ClusterNum))
 
 plot.temp <- ggplot(dat.enter, aes(x= species, y = value)) +
   geom_bar(stat = "identity", position = "dodge") +
   geom_hline(yintercept= c(.25, 0.75), linetype= "dashed", color= "grey") +
   geom_hline(yintercept= .5, color= "grey") +
   coord_flip() +
   xlab("species") + 
   ylab("Proportion present in total units") +
   theme_bw() +
   facet_wrap(.~ClusterNum, scales= "fixed", nrow=1, labeller= labeller(ClusterNum= var.names.temp)) +
   ggtitle(paste("Presence of catch in each area by species: ", title.temp)) +
   theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle= 90, vjust= 0))
 
 #plot.save <- plot_grid(graph.gear, graph.stratum, labels="AUTO", nrow=2)
 ggsave(paste0(graph.output.dir,"Fig8A_SpComp_", title.temp,".jpeg"), plot.temp)
 
}

species.proportion.graph.fn(dataset= kmediods.sp.prop.2, tot.units= tot.units.sp.prop.2b, title.temp= "Prop_2_KMediods")
species.proportion.graph.fn(dataset= kmediods.sp.prop.5, tot.units= tot.units.sp.prop.8b, title.temp= "Prop_5_KMediods")
#species.proportion.graph.fn(dataset= kmeans.sp.root.3, tot.units= tot.units.sp.root.3, title.temp= "Root_3")
#species.proportion.graph.fn(dataset= kmeans.sp.root.5, tot.units= tot.units.sp.root.5, title.temp= "Root_5")


 
data.set.enter <- kmeans.q.root.grpperc # kmeans.q.sub.grpperc2 kmeans.q.sub.grpperc7 kmeans.q.root.grpperc

datprop.enter <- data.set.enter %>% select(-Unit) %>% 
  rownames_to_column("ClusterNum") %>% gather(key= "species", value="value", factor_key=T, -ClusterNum) %>%
  mutate(ClusterNum= as.factor(ClusterNum))


dat.enter <- datprop.enter  # datroot.enter  datprop.enter
title.enter <- "Root"  # Root Prop

tot.units.clust <-  tot.units.clust.root  # tot.units.clust.root tot.units.clust.sub2 tot.units.clust.sub7

var.names <- paste("Cluster: ",unique(dat.enter$ClusterNum), "_(",tot.units.clust, ")", sep="")
names(var.names) <- c(unique(dat.enter$ClusterNum))

ggplot(dat.enter, aes(x= species, y = value)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("species") + 
  ylab("Proportion present in total units") +
  theme_bw() +
  facet_wrap(.~ClusterNum, scales= "fixed", ncol=7, labeller= labeller(ClusterNum= var.names)) +
  ggtitle(paste("Presence of catch in each area by species: ", title.enter)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle= 90, vjust= 0))

ggplot(dat.enter, aes(x= species, y = value, fill= ClusterNum)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("species") + 
  ylab("Proportion present in total units") +
  theme_bw() +
  ggtitle(paste("Presence of catch in each area by species: ", title.enter)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle= 90, vjust= 0))

dat.enter2 <- datroot.enter2  # datroot.enter datprop.enter2
title.enter2 <- "Root"  # Root Prop

var.names2 <- paste("Cluster: ",unique(dat.enter2$ClusterNum), "_(",tot.units.clust, ")", sep="")
names(var.names2) <- c(unique(dat.enter2$ClusterNum))

ggplot(dat.enter2, aes(x= species, y = value)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("species") + 
  ylab("Present in total units") +
  theme_bw() +
  facet_wrap(.~ClusterNum, scales= "fixed", ncol=1, labeller= labeller(ClusterNum= var.names2)) +
  ggtitle(paste("Presence of catch in each area by species: ", title.enter2)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle= 90, vjust = 0))

ggplot(dat.enter, aes(x= species, y = ClusterNum, size=value)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(1, 10)) +
  coord_flip() +
  xlab("species") + 
  ylab("Cluster") +
  theme_bw() +
  #facet_wrap(.~ClusterNum, scales= "fixed", ncol=1, labeller= labeller(ClusterNum= var.names)) +
  ggtitle(paste("Presence of catch in each area by species: ", title.enter)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle= 90, vjust= 0))

ggplot(cpue.sp.dat2, aes(x=species.new, y= Gear, size=LogTotCPUE)) +
  geom_point(alpha=0.5) +
  coord_flip() +
  xlab('Species') +
  ylab("Log(Total CPUE)") +
  facet_wrap(.~NMFS.Area, nrow=1) +
  theme_bw()

## Graphics for number of units (or sum of values) for each species in each stratum and gear

q.sub.stratcounts <- q.sub.dat %>% rownames_to_column(var= "Unit") %>% 
  separate("Unit", c("Year", "Month", "Stratum", "Gear"), remove=FALSE) %>%
  select(-c(Unit,Year, Month, Gear)) %>%
  group_by(Stratum) %>% summarise_all(sum) %>%
  gather(key= Stratum, factor_key=T) %>% rename(species= Stratum) %>%
  mutate(Statum2= rep(c("610", "620", "630", "640", "650"), length(names(q.sub.dat)))) %>%
  as.data.frame() 
q.sub.geartcounts <- q.sub.dat %>% rownames_to_column(var= "Unit") %>% 
  separate("Unit", c("Year", "Month", "Stratum", "Gear"), remove=FALSE) %>%
  select(-c(Unit,Year, Month, Stratum)) %>%
  group_by(Gear) %>% summarise_all(sum)
  summarise_all(sum) %>% as.data.frame()





################################################
#### CCA (Canonical Correspondence Analysis ####
################################################

## Variable dataset

factors.q.sub.dat <- as.data.frame(q.sub.dat) %>% 
  rownames_to_column("Unit") %>%
  separate("Unit", c("Year", "Month", "Stratum", "Gear"), remove=FALSE)

factors.dat <- factors.q.sub.dat %>%
  select(Year, Month, Stratum, Gear)

### Run CCA change only following lines
dat.enter <- q.sub.dat # q.sub.dat, q.root.dat
dat.name <- "Prop sub" # "Prop sub", "Root sub"

## Run CCA 
cca.results <- cca(dat.enter ~ Stratum + Gear, data= factors.dat) ## removed Month
cca.results.temp <- cca(dat.enter ~ Stratum + Gear, data= factors.dat) ## removed Month

#saveRDS(cca.results, file= paste0("G:\\My Drive\\MV Form Sp Complex Project\\Results\\Graphs&Tables\\CCA_proportion_dat.rds"))

anova(cca.results.temp, step=1000)
#anova(cca.results.temp, step=1000, by="axis")
vif.cca(cca.results.temp)

info <- (cca.results$CCA$tot.chi / cca.results$tot.chi)*100

## Save CCA centroids results (CCA values for each factor component)

cca.factors.centroids <- cca.results$CCA$centroids
general.cca.results <- cca.results
axis.eig.val <- cca.results$CCA$eig[1:4]/cca.results$tot.chi

varpart()

## plot cca results w/ % explained
ford.temp <- fortify(model=cca.results, axes= c(1:12))
ford <- fortify(cca.results, axes = c(1,2))  # fortify the ordination          ## AXIS CHANGE
ford.sites <- subset(ford, Score == 'sites') %>% 
  separate("Label", c("Year", "Month", "Stratum", "Gear"), remove=FALSE)

take <- c('CCA1', 'CCA2') # which columns contain the scores we want          ## AXIS CHANGE
arrows <- subset(ford, Score == 'centroids')  # take only biplot arrow scores
## multiplier for arrows to scale them to the plot range
mul <- ggvegan:::arrowMul(arrows[, take], subset(ford, select = take, Score == 'sites'))
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows

cca.plot <- ggplot() +
  geom_point(data = ford.sites,  # plotting sites
             mapping = aes(x = CCA1, y = CCA2, colour= factor(Gear))) +        ## AXIS CHANGE & Factor
  geom_text(data= subset(ford, Score == 'species'),
            mapping= aes(label = Label, x = CCA1, y = CCA2), colour= 'blue') +       ## AXIS CHANGE
  geom_segment(data = arrows, mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2), ## AXIS CHANGE
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text(data = arrows, # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1*1.1 , y = CCA2 *1.1 )) +         ## AXIS CHANGE
  coord_fixed() + 
  theme_bw() +
  labs(title= paste("CCA- ", dat.name), subtitle= paste(round(info, 3), "% explained") )

cca.plot

## Saving  CCA results
dat.name
#cca.results.root <- cca.results
#cca.factors.centriod.root <- as.data.frame(cca.results$CCA$centroids)
#min.cca.root <- as.data.frame(t(as.matrix(apply(ford.temp[,-c(1,2)], 2, min))))
#min.cca.root$tempname <- "min"
#max.cca.root <- as.data.frame(t(as.matrix(apply(ford.temp[,-c(1,2)], 2, max))))
#max.cca.root$tempname <- "max"
#abs.cca.root <- as.data.frame(t(as.matrix(apply())))
#cca.factors.centriod.root <- as.data.frame(cca.factors.centriod.root) %>% 
#  rownames_to_column(var= "tempname" ) %>%
#  full_join(min.cca.root) %>%
#  full_join(max.cca.root) %>% column_to_rownames(var= "tempname")


cca.results.prop <- cca.results
cca.factors.centriod.prop <- as.data.frame(cca.results$CCA$centroids)
min.cca.prop <- as.data.frame(t(as.matrix(apply(ford.temp[,-c(1,2)], 2, min))))
min.cca.prop$tempname <- "min"
max.cca.prop <- as.data.frame(t(as.matrix(apply(ford.temp[,-c(1,2)], 2, max))))
max.cca.prop$tempname <- "max"
cca.factors.centriod.prop <- as.data.frame(cca.factors.centriod.prop) %>% 
  rownames_to_column(var= "tempname" ) %>%
  full_join(min.cca.prop) %>%
  full_join(max.cca.prop) %>% column_to_rownames(var= "tempname")

getwd()
write.csv(cca.factors.centriod.root, "CCA_centroidfactors_root.csv")
write.csv(cca.factors.centriod.prop, "CCA_centroidfactors_prop.csv")

cca.factors.root.read.dat <- read.csv("CCA_centroidfactors_root.csv")
cca.factors.root.dat <- cca.factors.root.read.dat %>% rename(Factor= X) %>%
  select(Factor, CCA1, CCA2, CCA3) %>% 
  mutate(CCA1= abs(CCA1), CCA2= abs(CCA2),CCA3= abs(CCA3)) %>%
  filter(Factor !=c("min", "max"))
cca.factors.root.dat1 <- cca.factors.root.dat %>% gather(key="Factor", factor_key=T) %>%
  rename(CCA_axis= Factor) %>% mutate(Factor= rep(cca.factors.root.dat[,1], 3))


ggplot(cca.factors.root.dat1, aes(x= Factor, y = value)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Factors") + 
  ylab("CCA values") +
  theme_bw() +
  facet_wrap(.~CCA_axis, scales= "free", ncol=3) +
  ggtitle("CCA values for Root") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle= 90, vjust = 0))

##############
#### NMDS ####
##############

dat.enter <- q.sub.dat # q.sub.dat, q.root.dat
dat.name <- "Prop sub" # "Prop sub", "Root sub"

envir.dat <- factors.dat <- factors.q.sub.dat %>%
  select(Unit, Year, Month, Stratum, Gear)

dist.nmds <-  "bray"
k.val <- 2

nmds.results <- metaMDS(dat.enter, dist= dist.nmds, k= k.val, trymax= 50)
ordiplot(nmds.results, type="n", main=paste("NMDS- ", dat.name,"- Stress= ", 
                                                   round(nmds.results$stress, 3)))
#orditorp(nmds.results, display= "sites", air= 0.01, col="black")
ordiplot(nmds.results, type="p", main=paste("NMDS- ", dat.name,"- Stress= ", 
                                            round(nmds.results$stress, 3)))
orditorp(nmds.results, display="species",col="red",air=0.01, cex=1.1)

## With environmental dataset
env.nmds.results <- envfit(nmds.results, envir.dat)
plot(env.nmds.results)
env.nmds.results

## k.val= 3 converged for the root data
## k.val=  converged for the prop data

##########

### General summaries for datasets
head(q.sub.dat)

pos.count.fn(q.sub.dat)

temp <- apply(q.sub.dat, 2, function(x) {sum(x>0, na.rm=T)}) %>% as.data.frame() %>%
  rename("TotPos" = ".") %>% rownames_to_column(var="species") %>% arrange(desc(TotPos))
temp

temp2 <- colSums(q.sub.dat) %>% as.data.frame() %>% rename("TotSum" = ".") %>%  
  rownames_to_column(var="species") %>% arrange(desc(TotSum))
temp2

## Summary of number of gear vs stratum graph

gearstrat.temp <- q.sub.dat %>% rownames_to_column("Unit") %>% select(Unit) %>%
  separate(Unit, c("Year", "Month", "Stratum", "Gear")) %>% group_by(Gear, Stratum) %>%
  summarise(Counts = length(Month)) %>% as.data.frame

gearstrat.graph <- ggplot(gearstrat.temp, aes(x= Stratum, y = Counts)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Stratum") + 
  ylab("Counts") +
  theme_bw() +
  facet_wrap(.~Gear, scales= "fixed", ncol=4) +
  ggtitle("Gears-Statum Counts") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle= 90, vjust = 0))
gearstrat.graph

#### Species composition based on dataset ####

head(all.prop.sub.dat)
# do histograms? or bargraphs?? 

bin.fn <- function(x) {cut(x, breaks=c(-1,0.2, 0.4, 0.6, 0.8, Inf), labels= c("1", "2", "3", "4", "5"))}

sub.sp.com.dat <- all.prop.sub.dat %>% rownames_to_column(var="Unit") %>% 
  gather(key="species", value= "Prop",  factor_key=T, -Unit) %>%
  separate(Unit, c("Year", "Month", "Stratum", "Gear")) %>%
  mutate(Bin= bin.fn(Prop))

sp.comp.gear.dat <- sub.sp.com.dat %>% group_by(Gear, species, Bin) %>% 
  summarise(BinCount= sum(as.numeric(Bin), 20), Log_BinCount= log(BinCount-10)) %>% 
  as.data.frame()
sp.comp.stratum.dat <- sub.sp.com.dat %>% group_by(Stratum, species, Bin) %>% 
  summarise(BinCount= sum(as.numeric(Bin), 10), Log_BinCount= log(BinCount)) %>% as.data.frame()
#temp <- sub.sp.com.dat %>% filter(Prop > 0.75)

sp.comp.gear.plot <- ggplot(sp.comp.gear.dat, aes(x=species, y= Bin, size=BinCount)) +
  geom_point(alpha=0.5) +
  scale_size(range= c(5, 100), name= "Bin Count") +
  scale_size_area(breaks= c(0, 50, 100, 200, 1325)) +
  #scale_size(range= c(0, 8), name= "Log(Bin Count)") +   ## for natural log
  #scale_size_area(breaks= c(0, 2, 4, 6, 8)) +  ## for natural log
  coord_flip() +
  xlab('Species') +
  ylab("Bin Group") +
  facet_wrap(.~Gear, nrow=1) +
  theme_bw()
sp.comp.gear.plot

sp.comp.stratum.plot <- ggplot(sp.comp.stratum.dat, aes(x=species, y= Bin, size=BinCount)) +
  geom_point(alpha=0.5) +
  scale_size(range= c(0, 925), name= "Bin Count") +
  scale_size_area(breaks= c(0, 25, 100, 500, 925)) +
  #scale_size(range= c(0, 7), name= "Log(Bin Count)") +   ## for natural log
  #scale_size_area(breaks= c(2, 3, 4, 5, 6)) +  ## for natural log
  coord_flip() +
  xlab('Species') +
  ylab("Bin Group") +
  facet_wrap(.~Stratum, nrow=1) +
  theme_bw()
sp.comp.stratum.plot

###########################################################################################
##
## Other random summaries for paper
##
###########################################################################################

head(fish.dat)

unique(fish.dat$Year)

## Avg Catch per year

## Fisheries

avg.catch.gear.area.fish <- fish.dat %>% 
  filter(!species.new== "northern") %>%
  group_by(Year, Gear, NMFS.Area) %>%
  summarise(Tot.Catch.kg= sum(Catch.kg, na.rm=T)) %>%
  replace_na(list(Tot.Catch.kg=0)) %>%
  group_by(Gear, NMFS.Area) %>%
  summarise(avg.annual.catch= round(mean(Tot.Catch.kg),2)) %>% as.data.frame()

avg.catch.gear.fish <- fish.dat %>% group_by(Year, Gear) %>%
  filter(!species.new== "northern") %>%
  summarise(Tot.Catch.kg= sum(Catch.kg, na.rm=T)) %>%
  replace_na(list(Tot.Catch.kg=0)) %>%
  group_by(Gear) %>%
  summarise(avg.annual.catch= round(mean(Tot.Catch.kg),2)) %>% as.data.frame()
length(unique(fish.dat$Year))

## Surveys

head(trawlsurv.haul.expand.dat,3)
head(llsurv.haul.expand.dat,3)

length(unique(trawlsurv.haul.expand.dat$Year))
length(unique(llsurv.haul.expand.dat$Year))

## trawl
avg.catch.gear.area.trawl <- trawlsurv.haul.expand.dat %>% group_by(Year, Stratum) %>%
  filter(!name== "northern") %>%
  summarise(Tot.Catch.kg= sum(Weight..kg., na.rm=T)) %>%
  replace_na(list(Tot.Catch.kg=0)) %>%
  group_by(Stratum) %>%
  summarise(avg.annual.catch= round(mean(Tot.Catch.kg),2)) %>% as.data.frame()
avg.catch.gear.trawl <- trawlsurv.haul.expand.dat %>% group_by(Year) %>%
  filter(!name== "northern") %>%
  summarise(Tot.Catch.kg= sum(Weight..kg., na.rm=T)) %>%
  replace_na(list(Tot.Catch.kg=0)) %>%
  summarise(avg.annual.catch= round(mean(Tot.Catch.kg),2)) %>% as.data.frame()
num.trawls.avg <- trawlsurv.haul.expand.dat %>% group_by(Year) %>%
  summarise(Num.hauls= length(unique(Unit_small))) %>%
  summarise(Avg.hauls= mean(Num.hauls)) ; num.trawls.avg
num.sp.avg <- trawlsurv.haul.expand.dat %>% group_by(Year) %>%
  summarise(Num.sp= length(unique(name))) %>%
  summarise(Avg.hsp= mean(Num.sp)) ; num.sp.avg

sp.catch.trawl.all <- trawlsurv.haul.expand.dat %>% 
  filter(!name== "northern") %>%
  group_by(name) %>%
  summarise(Tot.Catch= sum(Weight..kg., na.rm=T), Tot.CPUE= sum(CPUE, na.rm=T)) %>%
  arrange(desc(Tot.CPUE)) %>% ungroup() %>%
  mutate(Cum.CPUE= cumsum(Tot.CPUE), Cum.Perc= round(Cum.CPUE/ sum(Tot.CPUE, na.rm=T)*100, 2),
         num.rank= row_number()) %>%
  arrange(desc(Tot.Catch)) %>%
  mutate(Cum.Catch= cumsum(Tot.Catch), Cum.Perc.Catch= round(Cum.Catch/ sum(Tot.Catch, na.rm=T)*100, 2) ) %>%
  as.data.frame()

sp.catch.trawl.list <- list(NA)
year.temp <- unique(trawlsurv.haul.expand.dat$Year)
for(i in 1:length(year.temp)) {

  sp.catch.trawl.temp <- trawlsurv.haul.expand.dat %>% 
    filter(Year== year.temp[i]) %>%
    filter(!name== "northern") %>%
    group_by(name) %>%
    summarise(Tot.Catch= sum(Weight..kg., na.rm=T), Tot.CPUE= sum(CPUE, na.rm=T)) %>%
    arrange(desc(Tot.CPUE)) %>% ungroup() %>%
    mutate(Cum.CPUE= cumsum(Tot.CPUE), Cum.Perc= round(Cum.CPUE/ sum(Tot.CPUE, na.rm=T)*100, 2),
           num.rank= row_number()) %>%
    arrange(desc(Tot.Catch)) %>%
    mutate(Cum.Catch= cumsum(Tot.Catch), Cum.Perc.Catch= round(Cum.Catch/ sum(Tot.Catch, na.rm=T)*100, 2) ) %>%
    mutate(Year= paste0(year.temp[i])) %>%
    as.data.frame()
    
  sp.catch.trawl.list[[i]] <- sp.catch.trawl.temp
}
sp.catch.trawl.full <- reduce(sp.catch.trawl.list, full_join)
write.csv(sp.catch.trawl.full, "G:\\My Drive\\MV Form Sp Complex Project\\Results\\Sp_Cum_Catch_Trawl.csv")
sp.catch.trawl.full2 <- sp.catch.trawl.full %>% filter(Cum.Perc < 96) %>% arrange(Year, num.rank) %>% as.data.frame() 
write.csv(sp.catch.trawl.full2, "G:\\My Drive\\MV Form Sp Complex Project\\Results\\Sp_Cum_Catch_Trawl_v2.csv")


## longline
head(ll.temp2)
#llsurv.haul.expand.dat <- ll.temp2

num.ll.avg <- llsurv.haul.expand.dat %>%
  unite("Unit_temp", c(Year, Unit_small), remove=FALSE) %>%
  group_by(Year) %>%
  summarise(Num.sets= length(unique(Unit_temp))) %>%
  summarise(Avg.sets= mean(Num.sets)) ; num.ll.avg
num.sp.avg.ll <- llsurv.haul.expand.dat %>% group_by(Year) %>%
  summarise(Num.sp= length(unique(name))) %>%
  summarise(Avg.hsp= mean(Num.sp)) ; num.sp.avg.ll


avg.catch.gear.area.trawl <- llsurv.haul.expand.dat %>% group_by(Year, Stratum) %>%
  summarise(Tot.Catch.kg= sum(Weight..kg., na.rm=T)) %>%
  replace_na(list(Tot.Catch.kg=0)) %>%
  group_by(Stratum) %>%
  summarise(avg.annual.catch= round(mean(Tot.Catch.kg),2)) %>% as.data.frame()
avg.catch.gear.ll <- llsurv.haul.expand.dat %>% group_by(Year) %>%
  summarise(Tot.Catch= sum(Catch, na.rm=T)) %>%
  replace_na(list(Tot.Catch=0)) %>%
  summarise(avg.annual.catch= round(mean(Tot.Catch),2)) %>% as.data.frame()
