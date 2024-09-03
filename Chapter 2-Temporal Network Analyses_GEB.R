interval_dat=read.csv("Data/Time_Intervals_introduced.csv", check.names = FALSE)
climate_dat=read.csv("Data/Halls Cave Climate Data for 16 Time Intervals.csv", check.names = FALSE)
names(climate_dat) <- iconv(names(climate_dat), to = "ASCII", sub = "")
mat_raw=as.matrix(read.csv("Data/Matrix.csv", row.names=1))
spp_attrs=read.csv("Data/spp_attributes.csv")
plei_sites=read.csv("Data/Pleistocene_Sites.csv")
holo_sites=read.csv("Data/Holocene_Sites.csv")
mods_per_bin=read.csv("Data/mods_per_bin.csv")
library(igraph)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(lme4) #package for running mixed-effects models
library(lmerTest) #package for extracting p-values form lme4 models
library(forcats)
library(bipartite)
library(cowplot)
library(ppcor)


# Re-add primary resources ####
prim.res_names=rownames(mat_raw)[which(rownames(mat_raw)%in%colnames(mat_raw)==FALSE)]
primary_mat=matrix(0, ncol=length(prim.res_names), nrow=nrow(mat_raw), dimnames=list(rownames(mat_raw), prim.res_names))
mat=cbind(primary_mat, mat_raw)

prim.res_names
append(as.character(interval_dat[,i]),prim.res_names)

time.list=list()
for(i in 1:ncol(interval_dat)){
  time.list[[i]]=append(as.character(interval_dat[,i]),prim.res_names)
  time.list[[i]]=time.list[[i]][time.list[[i]]!=""]
}
time.list

# List of species for each time bin and find species in this list that don't match the species in matrices ####
time.list

lapply(time.list, function(x) x[which(is.na(match(x, rownames(mat))))])

# These are species that are in the matrix twice based on isotopic signatures in pleistocene vs. holocene. We need to know how to consider this when doing the time interval analysis. 

  ## TIME BIN 9 (11339-12656ybp) contains the start of the Holocene     (~11700 ybp)


for(i in 1:16){
  if(i<9){
    time.list[[i]]=gsub("Lynx_rufus", "Lynx_rufus_holocene",time.list[[i]])
    time.list[[i]]=gsub("Lepus_californicus", "Lepus_californicus_holocene",time.list[[i]])
    time.list[[i]]=gsub("Panthera_onca", "Panthera_onca_holocene",time.list[[i]])
    
    
    
  } else{
    time.list[[i]]=gsub("Lynx_rufus", "Lynx_rufus_pleistocene",time.list[[i]])
    time.list[[i]]=gsub("Lepus_californicus", "Lepus_californicus_pleistocene",time.list[[i]])
    time.list[[i]]=gsub("Panthera_onca", "Panthera_onca_pleistocene",time.list[[i]])
    
  }
}
time.list
lapply(time.list, function(x) x[which(is.na(match(x, rownames(mat))))])

# Now, use this list of species to make matrix for each time bin ####

time.mats=lapply(time.list, function(x) mat[rownames(mat) %in% x, rownames(mat) %in% x])

sapply(time.mats, dim)

#names(time.mats)=colnames(interval_dat)
time.mats
str(time.mats)

# Make a food web network of each interval ####
time_nets=list()
for(i in 1:length(time.mats)){
  time_nets[[i]]=graph_from_adjacency_matrix(time.mats[[i]], "directed")
  V(time_nets[[i]])$category=spp_attrs$Trophic_Category[match(V(time_nets[[i]])$name, spp_attrs$Binomial)]
}

# Plotting Time Bin Networks ####
#pdf("Quentin&Dai/New_Plots/time_nets_trophic_groups.pdf", width=6)
par(mfrow=c(1,1), mar=c(1,1,2,1), pty = "s")
for(i in 1:length(time_nets)){
  time_colors=vector(length=length(V(time_nets[[i]])))
  time_colors[V(time_nets[[i]])$category=="browser"]="yellow"
  time_colors[V(time_nets[[i]])$category=="grazer"]="blue"
  time_colors[V(time_nets[[i]])$category=="carnivore"]="red"
  time_colors[V(time_nets[[i]])$category=="mixed feeder"]="brown"
  time_colors[V(time_nets[[i]])$category=="omnivore"]="orange"
  time_colors[V(time_nets[[i]])$category=="frugivore_granivore"]="green"
  time_colors[V(time_nets[[i]])$category=="insectivore"]="purple"
  time_colors[match(V(time_nets[[i]])$name, prim.res_names)]="white"
  V(time_nets[[i]])$color=time_colors      
  
  colnames(interval_dat) <- iconv(names(interval_dat), to = "ASCII", sub = "")
  
  time_colors[time_colors=="FALSE"]=NA
  plot(time_nets[[i]], vertex.label="", edge.arrow.size=0.5, vertex.color=time_colors)
}
#main=colnames(interval_dat)[i],

colnames()


dev.off()

# Issue [FIXED]: If things are looking weird in the network (loner nodes) #### 
  ## Use degree() function to find isolate species and further generate a list of each interval and find the source of the issure

## degree(time_nets[[i]])

## time.isolates.list=list()
## for(i in 1:length(time.mats)){
##   time.isolates.list[[i]]=colSums(time.mats[[i]])
## }
## time.isolates.list 
#Notice: The issue is that the all herbivores/frugivores/granivores are isolates. The primary resources are missing.


# Calculating Modularity of networks per bin ####

time_mods=vector(length=length(time_nets))
for(i in 1:length(time_nets)){
  time_mods[i]=modularity(cluster_edge_betweenness(time_nets[[i]], directed=T))
}
time_mods
write.csv(time_mods, file = "Quentin&Dai/New_Data/time_mods.csv")

#Checking General Plot
plot(time_mods, type="o", pch=19, xlab="Time intervals (recent to old)", ylab="Food Web Modularity (Edge Betweenness)")

#Creating a Bipartite Network
interval_dat_list=list()
for(i in 1:length(time.list)){
  interval_dat_list[[i]]=data.frame(spp=time.list[[i]]) %>% left_join(., spp_attrs, by=c("spp"="Binomial"))
}
interval_dat_list[[1]]


bipart.mats=list()
spp_counts=vector()
for(k in 1:length(time.list)){
  consumer.rows=which(interval_dat_list[[k]]$Trophic_Category=="carnivore"|interval_dat_list[[k]]$Trophic_Category=="omnivore")
  
  resource.rows=which(interval_dat_list[[k]]$Trophic_Category=="browser"|interval_dat_list[[k]]$Trophic_Category=="grazer"|interval_dat_list[[k]]$Trophic_Category=="mixed feeder"|interval_dat_list[[k]]$Trophic_Category=="frugivore_granivore"|interval_dat_list[[k]]$Trophic_Category=="insectivore")
  
  consumer.spp=interval_dat_list[[k]]$spp[consumer.rows]
  resource.spp=interval_dat_list[[k]]$spp[resource.rows]
  
  bipart.mats[[k]]=time.mats[[k]][match(resource.spp, rownames(time.mats[[k]])), match(consumer.spp, rownames(time.mats[[k]]))]
  spp_counts[[k]]=length(consumer.spp)+length(resource.spp)
}

plotweb(bipart.mats[[1]], method="normal",ybig=0.1, arrow="no", col.interaction="black", col.high="tomato", col.low="lightblue",labsize=1,  adj.low=c(1,0), adj.high=c(0,0.5), text.rot = 90, y.lim=c(-1,2.5))


pdf("Quentin&Dai/New_Plots/bipartite.pdf", width=10)
par(mfrow=c(2,2), mar=c(1,1,1,1))
for(i in 1:length(time_nets)){
  plotweb(bipart.mats[[i]], method="normal",ybig=0.1, col.interaction="black", arrow="no",  col.high="tomato", col.low="lightblue",labsize=1,  adj.low=c(1,0), adj.high=c(0,0.5), text.rot = 90, y.lim=c(-1,2.5))
}
dev.off()


#Calculating NOS() Modularity
??NOS

N_bar=sapply(bipart.mats, function(x) NOS(x)$Nbar)
plot(N_bar, type="o", pch=19, col="purple", ylab="Average NOS index (undirected)")

#Pool together data
interval.cat=colnames(interval_dat)
interval.start=-as.numeric(str_sub(interval.cat, 1, str_locate(interval.cat, "_")[,1]-1))
interval.end=-as.numeric(str_sub(interval.cat, str_locate(interval.cat, "_")[,1]+1, str_length(interval.cat)-3))


sum.data=data.frame(interval.name=interval.cat, interval.start=interval.start, interval.end=interval.end, spp.richness=sapply(time.list, length), web.modularity=time_mods, NOS=N_bar) %>%
  mutate(interval.mid=(interval.start+interval.end)/2)


ggplot(data=sum.data, aes(x=interval.mid)) +
  #geom_vline(xintercept = interval.start, color="gray") +
  geom_rect(xmin=interval.start, xmax=interval.end, ymin=-0.15, ymax=Inf, fill="gray", alpha=0.5) +
  geom_rect(xmin=interval.start[9], xmax=interval.end[9], ymin=-0.15, ymax=Inf, fill="#ffeda0", alpha=0.5) +
  geom_point(aes(y=web.modularity), color="blue") +
  geom_line(aes(y=web.modularity, color="Modularity")) +
  geom_point(aes(y=NOS), color="red") +
  geom_line(aes(y=NOS, color="NOS Index")) +
  scale_color_manual(values=c("blue", "red")) +
  theme_cowplot() +
  xlab("Time Interval (year before present)") +
  ylab("Modularity or NOS Index")

#Null Model Testing
null.nos=lapply(bipart.mats, function(x){
  a=nullmodel(x, N=100, method=1)
  sapply(a, function(y) NOS(y)$Nbar)
})

?nullmodel
NOS_null_meanse=sapply(null.nos, function(x) mean_se(x))
sum.data$NOS.mean.null=unlist(NOS_null_meanse[1,])
sum.data$NOS.se.high.null=unlist(NOS_null_meanse[2,])
sum.data$NOS.se.low.null=unlist(NOS_null_meanse[3,])


ggplot(data=sum.data, aes(x=interval.mid)) +
  #geom_vline(xintercept = interval.start, color="gray") +
  geom_rect(xmin=interval.start, xmax=interval.end, ymin=-0.15, ymax=Inf, fill="gray", alpha=0.5) +
  geom_rect(xmin=interval.start[9], xmax=interval.end[9], ymin=-0.15, ymax=Inf, fill="#ffeda0", alpha=0.5) +
  geom_point(aes(y=NOS), color="red") +
  geom_line(aes(y=NOS, color="NOS Index")) +
  geom_point(aes(y=NOS.mean.null), color="black") +
  geom_errorbar(aes(ymin=NOS.se.low.null, ymax=NOS.se.high.null), color="black") +
  scale_color_manual(values=c("red")) +
  theme_cowplot() +
  xlab("Time Interval (year before present)") +
  ylab("NOS Index")


null.model = lapply(bipart.mats, function(x){
  a=nullmodel(x, N=100, method=1)
})

plotweb(null.model[[1]][[4]], method="normal",ybig=0.1, arrow="no", col.interaction="black", col.high="tomato", col.low="lightblue",labsize=1,  adj.low=c(1,0), adj.high=c(0,0.5), text.rot = 90, y.lim=c(-1,2.5))

plotweb(bipart.mats[[1]], method="normal",ybig=0.1, arrow="no", col.interaction="black", col.high="tomato", col.low="lightblue",labsize=1,  adj.low=c(1,0), adj.high=c(0,0.5), text.rot = 90, y.lim=c(-1,2.5))




# Testing P-Value Distributions ####

#Calculate P-value for the 12th time bin ##

length(null.nos[[12]])

li = list(1:100)
for(i in length(li)){
  if(i >= 50){
    li[[i]]/100
    li[[i]]=li[[i]]/100
    
  } else{
    1-(li[[i]]/100)
    li[[i]]=1-(li[[i]]/100)
  }
}
li
null.list = null.nos[[12]]

data.frame(NOS=sort(unlist(null.list), decreasing = TRUE),
           Pvalue=unlist(li))

test.p = data.frame(NOS=sort(unlist(null.list), decreasing = TRUE),
                    Pvalue=unlist(li))

write.csv(test.p, file = "Quentin&Dai/New_Data/othertest.csv")


# The column names are numbers and aren't read correctly, convert the names into "ASCII". This will fix the weird symbols in front of the names
names(interval_dat) <- iconv(colnames(interval_dat), to = "ASCII", sub = "")

# Create a dataframe of Time and Climate ####
time.dat = data.frame(Interval_ID=1:16,
           interval=names(interval_dat),
           interval.start=interval.start, 
           interval.end=interval.end,
           TM = time_mods,
           SpC = mods_per_bin$Species_Counts, # Making Species counts to add as a column
           Tavg = climate_dat$Calculated_Average_Temp,
           Tmax = climate_dat$CCSM_Max_Temp_C,
           Tmax_SD = climate_dat$CCSM_Max_Temp_SD,
           Tmin = climate_dat$CCSM_Min_Temp_C,
           Tmin_SD = climate_dat$CCSM_Min_Temp_SD,
           Pr = climate_dat$CCSM_Precip_mm,
           PrCv = climate_dat$CCSM_Precip_CV) %>%
  mutate(interval.mid=(interval.start+interval.end)/2)
write.csv(time.dat, file = "Quentin&Dai/New_Data/time_dat.csv")



# We must account for any temporal autocorrelation (looking at change in TM with the change in other variables)
  # Using the lead() function (due to interval 16 being the start), we can subtract to get the difference between intervals. Therefore, calculating the change that occurs with a change in modularity over time

?lag

time.dat.change = data.frame(Interval_ID=1:16, 
                      interval=names(interval_dat),
                      interval=names(interval_dat),
                      interval.start=interval.start, 
                      interval.end=interval.end,
                      TM.change = time_mods-lag(time_mods),
                      Tavg.change = climate_dat$Calculated_Average_Temp-lag(climate_dat$Calculated_Average_Temp),
                      SpC.change = mods_per_bin$Species_Counts-lag(mods_per_bin$Species_Counts),
                      Pr.change = climate_dat$CCSM_Precip_mm-lag(climate_dat$CCSM_Precip_mm),
                      Tmax.change = climate_dat$CCSM_Max_Temp_C-lag(climate_dat$CCSM_Max_Temp_C),
                      Tmax_SD.change = climate_dat$CCSM_Max_Temp_SD-lag(climate_dat$CCSM_Max_Temp_SD),
                      Tmin.change = climate_dat$CCSM_Min_Temp_C-lag(climate_dat$CCSM_Min_Temp_C),
                      Tmin_SD.change = climate_dat$CCSM_Min_Temp_SD-lag(climate_dat$CCSM_Min_Temp_SD),
                      PrCv.change = climate_dat$CCSM_Precip_CV-lag(climate_dat$CCSM_Precip_CV)) %>%
  mutate(interval.mid=(interval.start+interval.end)/2)
write.csv(time.dat.change, file = "Quentin&Dai/New_Data/time_dat_change.csv")

  




  
#run correlation
par(mfrow=c(1,1), mar=c(1,1,1,1))

cor.test(time.dat.change$Tavg.change, time.dat.change$TM.change, method = "spearman")
plot(time.dat.change$Tavg.change, time.dat.change$TM.change, type="p", col="red", pch=16, xlab = "Mean Temperature (°C)", ylab = "Modularity")
abline(lm(time.dat.change$TM.change~time.dat.change$Tavg.change))


cor.test(time.dat.change$Pr.change, time.dat.change$TM.change, method = "spearman")
plot(time.dat.change$Pr.change, time.dat.change$TM.change, type="p", col="blue", pch=16, xlab = "Precipitation (mm)", ylab = "Modularity")
abline(lm(time.dat.change$TM.change~time.dat.change$Pr.change))


cor.test(time.dat.change$SpC.change, time.dat.change$TM.change, method = "spearman")
plot(time.dat.change$SpC.change, time.dat.change$TM.change, type="p", col="green", pch=16, xlab = "Species Richness", ylab = "Modularity")
abline(lm(time.dat.change$TM.change~time.dat.change$SpC.change))

dev.off()

cor.test(time.dat$Tavg, time.dat$TM)
plot(time.dat$Tavg, time.dat$TM)
abline(lm(time.dat$TM~time.dat$Tavg))


# Plotting Modularity over Time Intervals ####
par(mfrow=c(1,1), mar=c(4,4,1,1)) #Resetting plot margins to default

fct_inorder(interval)

ggplot(time.dat.change, aes(x = fct_inorder(interval), y = TM.change)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.background = element_blank(), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Modularity") +
  coord_cartesian(ylim = c(0.0235, 0.5))

ggplot(time.dat, aes(x = fct_inorder(interval), y = TM, group = 1)) +
  geom_rect(xmin=interval.start, xmax=interval.end, ymin=-0.15, ymax=Inf, fill="gray", alpha=0.5) +
  geom_rect(xmin=interval.start[9], xmax=interval.end[9], ymin=-0.15, ymax=Inf, fill="#ffeda0", alpha=0.5) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.background = element_blank(), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Modularity")


#################################################################################

test.1 = ggplot(data=time.dat.change, aes(x=interval.mid)) +
  geom_rect(xmin=interval.start, xmax=interval.end, ymin=-0.25, ymax=Inf, fill="gray", alpha=0.5) +
  geom_rect(xmin=interval.start[9], xmax=interval.end[9], ymin=-0.25, ymax=Inf, fill="#ffeda0", alpha=0.5) +
  geom_point(aes(y=TM.change), color="black") +
  geom_line(aes(y=TM.change)) +
  
  scale_color_manual(values=c("black")) +
  theme_cowplot() +
  xlab("Time Interval (year before present)") +
  ylab("Modularity") +
  theme(axis.text.x=element_blank()) +
  xlab("")


test.2 = ggplot(data=time.dat.change, aes(x=interval.mid)) +
  geom_rect(xmin=interval.start, xmax=interval.end, ymin=-2.5, ymax=Inf, fill="gray", alpha=0.5) +
  geom_rect(xmin=interval.start[9], xmax=interval.end[9], ymin=-2.5, ymax=Inf, fill="#ffeda0", alpha=0.5) +
  geom_point(aes(y=Tavg.change), color="black") +
  geom_line(aes(y=Tavg.change)) +
  scale_color_manual(values=c("black")) +
  theme_cowplot() +
  xlab("Time Interval (year before present)") +
  ylab("Mean Temperature (°C)") +
  theme(axis.text.x=element_blank()) +
  xlab("")

test.3 = ggplot(data=time.dat.change, aes(x=interval.mid)) +
  geom_rect(xmin=interval.start, xmax=interval.end, ymin=-25, ymax=Inf, fill="gray", alpha=0.5) +
  geom_rect(xmin=interval.start[9], xmax=interval.end[9], ymin=-25, ymax=Inf, fill="#ffeda0", alpha=0.5) +
  geom_point(aes(y=Pr.change), color="black") +
  geom_line(aes(y=Pr.change)) +
  scale_color_manual(values=c("black")) +
  theme_cowplot() +
  xlab("Time Interval (year before present)") +
  ylab("Precipitation (mm)") 


test.4 = ggplot(data=time.dat.change, aes(x=interval.mid)) +
  geom_rect(xmin=interval.start, xmax=interval.end, ymin=-25, ymax=Inf, fill="gray", alpha=0.5) +
  geom_rect(xmin=interval.start[9], xmax=interval.end[9], ymin=-25, ymax=Inf, fill="#ffeda0", alpha=0.5) +
  geom_point(aes(y=SpC.change), color="black") +
  geom_line(aes(y=SpC.change)) +
  scale_color_manual(values=c("black")) +
  theme_cowplot() +
  xlab("Time Interval (year before present)") +
  ylab("Species Richness")


plot_grid(test.1, test.2, test.3,
          nrow = 3,
          labels = "",
          label_size = 12,
          align = "v",
          rel_heights=c(1,1,1))


plot_grid(test.1, test.4,
          nrow = 2,
          labels = "",
          label_size = 12,
          align = "v",
          rel_heights=c(1,1,1))

#Reverse the order of the time bins using the "rev()" function (to be more chronologically accurate)

fct_rev(fct_inorder(interval))

mod.time = ggplot(time.dat.change, aes(x = fct_rev(fct_inorder(interval)), y = TM.change)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.background = element_blank(), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Modularity") +
  coord_cartesian(ylim = c(0.0235, -0.5))

mod.time

#Make a line plot
mod.time.line = ggplot(time.dat.change, aes(x = fct_rev(fct_inorder(interval)), y = TM.change, group = 1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.background = element_blank(), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Modularity") +
  theme(axis.text.x=element_blank()) +
  xlab("")

mod.time.line

#Plot temperature and precipitation to modularity
temp.mod = ggplot(time.dat.change, aes(x = Tavg.change, y = TM.change)) +
  geom_line() +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black")) +
  labs(x="Average Temperature (°C)", y="Modularity")

temp.mod

precip.mod = ggplot(time.dat.change, aes(x = Pr.change, y = TM.change)) +
  geom_line() +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black")) +
  labs(x="Precipitation (mm)", y="Modularity")

precip.mod

#Plot temperature and precipitation over time bins
temp.time = ggplot(time.dat.change, aes(x = fct_rev(fct_inorder(interval)), y = Tavg.change, group = 1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.background = element_blank(), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Mean Temperature (°C)") +
  theme(axis.text.x=element_blank()) +
  xlab("")

temp.time

Precip.time = ggplot(time.dat.change, aes(x = fct_rev(fct_inorder(interval)), y = Pr.change, group = 1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.background = element_blank(), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Precipitation (mm)")

Precip.time


plot_grid(mod.time.line, temp.time, Precip.time, 
          nrow = 3,
          labels = "",
          label_size = 12,
          align = "v",
          rel_heights=c(1,1,1.7))

#Plotting total species per bin

ggplot(time.dat.change, aes(x = fct_rev(fct_inorder(interval)), y = SpC.change, fill = SpC.change)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.background = element_blank(), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Number of Species") +
  coord_cartesian(ylim = c(15, -25))

na.omit(Tavg.change)

# The data needs to be standardized so that the values are leveled.

# Let's first scale the data by z-transformation ####


# Making vectors of each variable ####
TM.change	<- na.exclude(time.dat.change$TM)
Tavg.change <- na.exclude(time.dat.change$Tavg.change)
SpC.change	<-  na.exclude(time.dat.change$SpC.change)
Pr.change	<-  na.exclude(time.dat.change$Pr.change)
Tmax.change	<-  na.exclude(time.dat.change$Tmax.change)
Tmax_SD.change	<-  na.exclude(time.dat.change$Tmax_SD.change)
Tmin.change	<-  na.exclude(time.dat.change$Tmin.change)
Tmin_SD.change	<-  na.exclude(time.dat.change$Tmin_SD.change)
PrCv.change	<-  na.exclude(time.dat.change$PrCv.change)

# Z-transforming Data ####
TM.z <- (TM.change - mean(TM.change)) / sd(TM.change)
Tavg.z <- (Tavg.change - mean(Tavg.change)) / sd(Tavg.change)
SpC.z	<- (SpC.change - mean(SpC.change)) / sd(SpC.change)
Pr.z	<- (Pr.change - mean(Pr.change)) / sd(Pr.change)
Tmax.z	<- (Tmax.change - mean(Tmax.change)) / sd(Tmax.change)
Tmax_SD.z	<- (Tmax_SD.change - mean(Tmax_SD.change)) / sd(Tmax_SD.change)
Tmin.z	<- (Tmin.change - mean(Tmin.change)) / sd(Tmin.change)
Tmin_SD.z	<- (Tmin_SD.change - mean(Tmin_SD.change)) / sd(Tmin_SD.change)
PrCv.z	<- (PrCv.change - mean(PrCv.change)) / sd(PrCv.change)

# Checking the z-transformation ####
mean(TM.z) #the mean of the variable should now be 0 (close to 0)
sd(TM.z) #the standard deviation of the variable should now be 1.


# Making a new data frame of z-transformed data ####

time.dat_z = data.frame(Interval=time.dat$interval, 
                        TM.z=TM.z,
                        Tavg.z=Tavg.z,
                        SpC.z=SpC.z,	
                        Pr.z=Pr.z,	
                        Tmax.z=Tmax.z,	
                        Tmax_SD.z=Tmax_SD.z, 
                        Tmin.z=Tmin.z,	
                        Tmin_SD.z=Tmin_SD.z,	
                        PrCv.z=PrCv.z)
write.csv(time.dat_z, file = "time_dat_z.csv")

# Forming a combined plot ####
test_1=ggplot(time.dat, aes(x = fct_rev(fct_inorder(interval)), y = TM.z, group = 1)) +
  geom_line(size = 1) +
  geom_point(aes(y=TM.z), color="black") +
  geom_point(aes(y=TM.z, color="Mod")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Modularity")

ggplot(time.dat, aes(x = fct_rev(fct_inorder(interval)), y = Tavg.z, group = 1)) +
  geom_line(size = 1) +
  geom_point(aes(y=Tavg.z), color="red") +
  geom_line(aes(y=Tavg.z, color="Avg Temp (°C)")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Modularity")


ggplot(time.dat, aes(x = fct_rev(fct_inorder(interval)), y = Pr.z, group = 1)) +
  geom_line(size = 1) +
  geom_point(aes(y=Pr.z), color="blue") +
  geom_line(aes(y=Pr.z, color="Precip (mm)")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.line = element_line(color = "black")) +
  labs(x="Time Bins", y="Z-Trasformed Variables")


plot_grid()


# Test Spearman-Rank Correlation ####

 #Spearman_Rank Correlation here is a statistical test being used to understand how correlated temperature, precipitation, and species richness is to the change we see in modularity

spearman_Temp_Mod <- cor.test(time.dat$TM,time.dat$Tavg,method="spearman");
screen_results4 <- c(paste("Temperature explain ",round(100*((spearman_Temp_Mod$estimate)^2),1),"% of the variance in ranked modularity-scored.",sep=""),
                     paste("p = ",round(spearman_Temp_Mod$p.value/(10^floor(log10(spearman_Temp_Mod$p.value))),2)*(10^floor(log10(spearman_Temp_Mod$p.value))),".",sep=""));
print(screen_results4);

spearman_Precip_Mod <- cor.test(time.dat$TM,time.dat$Pr,method="spearman");
screen_results4 <- c(paste("Precipitation explain ",round(100*((spearman_Precip_Mod$estimate)^2),1),"% of the variance in of the variance in ranked modularity-scored.",sep=""),
                     paste("p = ",round(spearman_Precip_Mod$p.value/(10^floor(log10(spearman_Precip_Mod$p.value))),2)*(10^floor(log10(spearman_Precip_Mod$p.value))),".",sep=""));
print(screen_results4);

spearman_rich_Mod <- cor.test(time.dat$TM,time.dat$SpC,method="spearman");
screen_results4 <- c(paste("Species Richness explain ",round(100*((spearman_rich_Mod$estimate)^2),1),"% of the variance in of the variance in ranked modularity-scored.",sep=""),
                     paste("p = ",round(spearman_rich_Mod$p.value/(10^floor(log10(spearman_rich_Mod$p.value))),2)*(10^floor(log10(spearman_rich_Mod$p.value))),".",sep=""));
print(screen_results4);

## Using the standardized data ##
spearman_Temp_Mod <- cor.test(time.dat.change$Tavg.change, time.dat.change$TM.change,method="spearman");
screen_results4 <- c(paste("Temperature explain ",round(100*((spearman_Temp_Mod$estimate)^2),1),"% of the variance in ranked modularity-scored.",sep=""),
                     paste("p = ",round(spearman_Temp_Mod$p.value/(10^floor(log10(spearman_Temp_Mod$p.value))),2)*(10^floor(log10(spearman_Temp_Mod$p.value))),".",sep=""));
print(screen_results4);

spearman_Precip_Mod <- cor.test(time.dat.change$Pr.change, time.dat.change$TM.change,method="spearman");
screen_results4 <- c(paste("Precipitation explain ",round(100*((spearman_Precip_Mod$estimate)^2),1),"% of the variance in of the variance in ranked modularity-scored.",sep=""),
                     paste("p = ",round(spearman_Precip_Mod$p.value/(10^floor(log10(spearman_Precip_Mod$p.value))),2)*(10^floor(log10(spearman_Precip_Mod$p.value))),".",sep=""));
print(screen_results4);

spearman_rich_Mod <- cor.test(time.dat.change$SpC.change, time.dat.change$TM.change,method="spearman");
screen_results4 <- c(paste("Species Richness explain ",round(100*((spearman_rich_Mod$estimate)^2),1),"% of the variance in of the variance in ranked modularity-scored.",sep=""),
                     paste("p = ",round(spearman_rich_Mod$p.value/(10^floor(log10(spearman_rich_Mod$p.value))),2)*(10^floor(log10(spearman_rich_Mod$p.value))),".",sep=""));
print(screen_results4);

# So, what statistic correlates best with modularity? ####
 # Here, all variables of climate and precipitation is included

key_stats <- c("Tavg.change",	"Pr.change", "SpC.change", "Tmax.change", "Tmax_SD.change", "Tmin.change", "Tmin_SD.change");
spearman_summaries <- data.frame(stat=as.character(key_stats),rho=as.numeric(rep(0,length(key_stats))),pvalue=as.numeric(rep(0,length(key_stats))));
for (ks in 1:length(key_stats)) {
  spearman_test <- cor.test(time.dat.change[,match(key_stats[ks],colnames(time.dat.change))],time.dat.change$TM.change,method="spearman");
spearman_summaries$rho[ks] <- round(spearman_test$estimate,3);
spearman_summaries$pvalue[ks] <- spearman_test$p.value;
if (spearman_summaries$pvalue[ks]>10^-2)  {
  spearman_summaries$pvalue[ks] <- round(spearman_summaries$pvalue[ks],3)
  } else  {
    spearman_summaries$pvalue[ks] <- round(spearman_summaries$pvalue[ks]/(10^floor(log10(spearman_summaries$pvalue[ks]))),2)*(10^floor(log10(spearman_summaries$pvalue[ks])))
  }
}
spearman_summaries$rho2 <- spearman_summaries$rho^2;
write.csv(spearman_summaries, file = "spearman.csv")

# Plot ###
# par(mfrow=c(4,2), mar = c(1, 5, 1, 10))
# plot(time.dat$TM,time.dat$Tavg,pch=21,bg="green", xlab ="TM", ylab ="Tavg");
# abline(lm(time.dat$Tavg~time.dat$TM))
# 
# plot(time.dat$TM,time.dat$Pr,pch=21,bg="yellow", xlab ="TM", ylab ="Pr");
# abline(lm(time.dat$Pr~time.dat$TM))
# 
# plot(time.dat$TM,time.dat$SpC,pch=21,bg="blue", xlab = "TM", ylab = "SpC");
# abline(lm(time.dat$SpC~time.dat$TM))
# 
# plot(time.dat$TM,time.dat$Tmax_SD,pch=21,bg="red", xlab ="TM", ylab ="Tmax_SD");
# abline(lm(time.dat$Tmax_SD~time.dat$TM))
# 
# plot(time.dat$TM,time.dat$Tmax,pch=21,bg="red", xlab ="TM", ylab ="Tmax");
# abline(lm(time.dat$Tmax~time.dat$TM))
# 
# plot(time.dat$TM,time.dat$Tmin_SD,pch=21,bg="purple", xlab ="TM", ylab ="Tmin_SD");
# abline(lm(time.dat$Tmin_SD~time.dat$TM))
# 
# plot(time.dat$TM,time.dat$Tmin,pch=21,bg="purple", xlab ="TM", ylab ="Tmin");
# abline(lm(time.dat$Tmin~time.dat$TM))
# 
# plot(time.dat$TM,time.dat$PrCv,pch=21,bg="black", xlab ="TM", ylab ="PrCv");
# abline(lm(time.dat$PrCv~time.dat$TM))

# Spearman Partial Correlations ####
pcor_1 = pcor.test(na.exclude(time.dat.change$TM.change),na.exclude(time.dat.change$Tavg.change),na.exclude(time.dat.change$Pr.change),method="spearman");

pcor_2 = pcor.test(na.exclude(time.dat.change$TM.change),na.exclude(time.dat.change$Pr.change),na.exclude(time.dat.change$Tavg.change),method="spearman");

pcor_3 = pcor.test(na.exclude(time.dat.change$TM.change),na.exclude(time.dat.change$Tavg.change),na.exclude(time.dat.change$SpC.change),method="spearman");

pcor_4 = pcor.test(na.exclude(time.dat.change$TM.change),na.exclude(time.dat.change$SpC.change),na.exclude(time.dat.change$Tavg.change),method="spearman");

pcor_5 = pcor.test(na.exclude(time.dat.change$TM.change),na.exclude(time.dat.change$Pr.change),na.exclude(time.dat.change$SpC.change),method="spearman");

pcor_6 = pcor.test(na.exclude(time.dat.change$TM.change),na.exclude(time.dat.change$SpC.change),na.exclude(time.dat.change$Pr.change),method="spearman");

pcor_all = data.frame(rbind(pcor_1,pcor_2,pcor_3,pcor_4,pcor_5,pcor_6))

pcor_all$correlation <- c("Mod + Temp - Precip", "Mod + Precip - Temp", "Mod + Temp - Rich", "Mod + Rich - Temp", "Mod + Precip- Rich", 'Mod + Rich - Precip')

pcor_all
write.csv(pcor_all, file = "partial_correlations.csv")


# multiple regression ####
mod_temp_lm <- summary(lm(time.dat$TM ~ time.dat$Tavg));
mod_temp_lm$r.squared;
mod_temp_lm$adj.r.squared;
mod_temp_lm$coefficients;
mod_temp_lm$cov.unscaled;

mod_temp_precip_lm <- summary(lm(time.dat$TM ~ (time.dat$Tavg+time.dat$Pr)));
mod_temp_precip_lm$r.squared;
mod_temp_precip_lm$adj.r.squared;
mod_temp_precip_lm$coefficients;
mod_temp_precip_lm$cov.unscaled;

mod_temp_precip_sp_lm <- summary(lm(time.dat$TM ~ (time.dat$Tavg+time.dat$Pr+time.dat$SpC)));
mod_temp_precip_sp_lm$r.squared;
mod_temp_precip_sp_lm$adj.r.squared;
mod_temp_precip_sp_lm$coefficients;
mod_temp_precip_sp_lm$cov.unscaled;

multiple_regressions_all <- summary(lm(time.dat$TM ~ (time.dat$Tavg+time.dat$Pr+time.dat$SpC+time.dat$Tmax+time.dat$Tmax_SD+time.dat$Tmin+time.dat$Tmin_SD+time.dat$Pr+time.dat$PrCv)));
multiple_regressions_all$r.squared;
multiple_regressions_all$adj.r.squared;
multiple_regressions_all$coefficients;
multiple_regressions_all$cov.unscaled;

# Tmin returns NA in regression due to singularities. This means the predictor Tmin seems to be perfectly correlated. A mixed-effects model could likely represent these to account for multiple observations affecting modularity. #


# Mixed-Effects Model: Partial Correlations ####

# The model has these components:
# Response variable = Modularity
# 6 predictor variables = "Precipitation", "Temperature (Avg)" "Richness", "Composition", "Mass Kurtosis" "Mass Skew"

# One random effect = " "?
  #Try pre- and post- ext............(WORKS!)
  #Try interval number 1-16......(Nope!)

## time_mods Species_Counts	climate.CCSM_Precip_mm	climate.CCSM_Precip_CV	climate.CCSM_Max_Temp_C	climate.CCSM_Max_Temp_SD	climate.CCSM_Min_Temp_C	climate.CCSM_Min_Temp_SD	climate.Calculated_Average_Temp




# Adding a status column to identify/group data as pre- and post- extinction 

time.dat$status <- c("Pre_Extinction", "Pre_Extinction", "Pre_Extinction", "Pre_Extinction", "Pre_Extinction", "Pre_Extinction", 'Pre_Extinction', "Pre_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction")

time.dat_z$status <- c("Pre_Extinction", "Pre_Extinction", "Pre_Extinction", "Pre_Extinction", "Pre_Extinction", "Pre_Extinction", 'Pre_Extinction', "Pre_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction", "Post_Extinction")




summary(lm(TM ~ Pr, data = time.dat))
summary(lm(TM.z ~ Pr.z, data = time.dat_z))


lmer(TM ~ Pr  + Tavg + SpC + (1 | status), data = time.dat)
summary(lmer(TM ~ Pr  + Tavg + SpC + (1 | status), data = time.dat))

lmer(TM.z ~ Pr.z + Tavg.z + SpC.z + (1 | status), data = time.dat_z)
summary(lmer(TM.z ~ Pr.z + Tavg.z + SpC.z + (1 | status), data = time.dat_z))



lmer(TM ~ scale(Pr)  + scale(Tavg) + scale(SpC) + (1 | status), data = time.dat)
summary(lmer(TM ~ scale(Pr)  + scale(Tavg) + scale(SpC) + (1 | status), data = time.dat))

summary(lmer(TM.z ~ Pr.z + Tavg.z + SpC.z + Tmax.z + Tmax_SD.z + Tmin.z + Tmin_SD.z + PrCv.z + (1 | status), data = time.dat_z))
