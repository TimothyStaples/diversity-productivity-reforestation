#################################################   ####
# Code goal: Evaluate all richness and trait models ####
# Author:    Tim Staples                            ####
# Date:      02/09/2015                             ####  
################################################### ####
# Libraries ####
rm(list=ls())
library(raster) # for environmental rasters
library(rgdal) # for environmental rasters
library(maptools) # manipulating environmental rasters and plotting
library(usdm) # for variance-inflation factor comparison
library(vegan) # for rarefied richness etc
library(nlme) # mixed-effect modelling
library(lme4) # mixed-effect modelling
library(MuMIn) # model selection and averaging
library(gamm4) # additive modelling
library(parallel) # for parallel computation of some apply and model selection functions
library(FD) # to calculate functional diversity indices
library(fields) # for k-means clustering
library(multcomp) # for glht
library(geoR) # for variogram
library(plotrix) # for gradient rectangles in plots

setwd("/home/timothy/Dropbox/Tim/PhD/Data/diversity-productivity-reforestation/")
# Global functions ####
sapply(list.files(path="./Functions", pattern=".R", full.names=TRUE),
       source)

# ####
# DATA PREP ####
# ####
# INITIAL DATA IMPORT ####
  
  species<-read.csv("./Data/species.list.csv", header=T)
  site<-read.csv("./Data/plotscale.csv", header=T)
  plant<-read.csv("./Data/plantscale.csv", header=T)
  
# read in existing environment data if not being recalculated

# Check whether any of our plots don't match up between plot info and plant info - 0s mean everything matches
unique(plant$pl.p[is.na(match(plant$pl.p, site$pl.p))])
unique(site$pl.p[is.na(match(site$pl.p, plant$pl.p))])
unique(plant$planting[is.na(match(plant$planting, site$planting))])
unique(site$planting[is.na(match(site$planting, plant$planting))])

#             EXTRACT ENVIRONMENTAL CONDITIONS AT PLOTS ####

shape.file.dir<-paste(ifelse(Sys.info()['sysname']=="Linux", "/home/Storage HDD/", "F:/"),
            "University files/Shape files", sep="")

data.coords<-as.data.frame(cbind(site$long, site$lat))
colnames(data.coords)<-c("Longitude","Latitude")

point.values<-sapply(list.files(path=shape.file.dir,
                                pattern=".tif"), function(x){
                                  extract(raster(paste0(shape.file.dir,"/",x)), data.coords, method="simple")
                                })

colnames(point.values)<-substr(colnames(point.values), 
                               1, regexpr(".tif", colnames(point.values))-1)

site<-cbind(site, point.values)

IBRA<-readShapeSpatial(paste0(shape.file.dir,
                              "/Shape files/",
                              "IBRA7_regions/ibra7_regions.shp"))
IBRAsub<-readShapeSpatial(paste0(shape.file.dir,
                                 "/Shape files/",
                                 "IBRA7_subregions_states/IBRA7_subregions_states.shp"))

data.coords.nona<-data.coords[rowSums(is.na(data.coords))==0,]
data.coords.plp<-site$pl.p[rowSums(is.na(data.coords))==0]
coordinates(data.coords.nona)<-c("Longitude","Latitude")

IBRA.points<-data.frame(IBRA=over(data.coords.nona,IBRA)$REG_CODE_7,
                        IBRAsub=over(data.coords.nona,IBRAsub)$SUB_CODE_7,
                        pl.p=site$pl.p[rowSums(is.na(data.coords))==0])

site1<-merge(site, IBRA.points, all.x=TRUE, all.y=FALSE, by.x="pl.p", by.y="pl.p")

site1<-cbind(pl.p=site1$pl.p,
             site1[,c("IBRA", "IBRAsub")],
             site1[,!colnames(site1) %in% c("pl.p","IBRA", "IBRAsub")])

head(site1)

site<-site1

#             CALCULATE PLANTING BIOMASS & PLANT DENSITY ####


# Pull out just the planting to start
mixed.site<-droplevels(subset(site, plot.type=="Planting"))
mixed.plant<-droplevels(plant[plant$pl.p %in% mixed.site$pl.p,])

# make a plot-level mass variable by summing plant-scale mass per plot
# we're back using the big mixed.plant dataframe as we want to include
# the mass of dead/unknown plants
mass<-aggregate(tot.ag.mass ~ pl.p, data=mixed.plant, sum)
hist(mass$tot.ag.mass, breaks=100)

# add plot-scale biomas to plot-scale dataframe
mixed.site<-merge(mixed.site, mass, by.x="pl.p", by.y="pl.p")

# now we want biomass per hectare, so divide biomass by plot area
mixed.site$plot.area<-as.numeric(as.character(mixed.site$plot.area))
mixed.site$biomass.area<-mixed.site$tot.ag.mass/mixed.site$plot.area
hist(mixed.site$biomass.area, breaks=100)

# remove all unknown/dead entries (so we are just calculating alpha on ID'd species)
mixed.plant.alive<-mixed.plant[mixed.plant$species %in% 
                               species$species[species$dead==0],]

# it might also be interesting to have plant density as a model predictor
density<-aggregate(species ~ pl.p, data=mixed.plant.alive, length)
colnames(density)<-c("pl.p","plant.count")
hist(density$plant.count, breaks=100)

# add this measure to plot-scale dataframe
mixed.site<-merge(mixed.site, density, by.x="pl.p", by.y="pl.p")
colnames(mixed.site)

# Some plantings only have a small sample size - what proportion do we lose by 
# setting a minimum value? What proportion of plots do we need to exclude in 
# order to achieve a minimum tree count per plot?
mixed.site.temp<-mixed.site
mixed.site.temp$size.prop<-NA
for (i in (1:length(mixed.site.temp$size.prop))){
mixed.site.temp$size.prop[i]<-length(mixed.site.temp$plant.count[
                                            mixed.site.temp$plant.count<
                                            mixed.site.temp$plant.count[i]])/
                                            length(mixed.site.temp$plant.count)}

plot(mixed.site.temp$size.prop ~ mixed.site.temp$plant.count, xlim=c(0,100))
# from this plot it looks like we don't want to have our minimum sample size as
# anything more than 30, and we probably don't want it anything less than 10 - 
# let's run with 30 and see how it goes

# we want to remove plantings with only a small number of plants (less than 30)
mixed.site.sub<-droplevels(mixed.site[mixed.site$plant.count>=30,])
summary(mixed.site.sub$plant.count)

# sub-set plant scale dataset to match
mixed.plant.sub<-droplevels(mixed.plant.alive[mixed.plant.alive$pl.p %in%
                                              mixed.site.sub$pl.p,])

# get a plant number per hectare density score
mixed.site.sub$density<-mixed.site.sub$plant.count/mixed.site.sub$plot.area
hist(mixed.site.sub$density, breaks=100)

#             CALCULATE SPECIES RICHNESS ####

# we need to calculate a site-species matrix without deads to get rarefied richness
# we need to do something with the unknowns, so I'm going to count them as a single 'species',
# as the data collector knew that they weren't whatever else was IDed so it seems reasonable and
# conservative to count them as 1 additional species
mixed.plant.sub$species[mixed.plant.sub$species %in% 
                        species$species[species$unknown==1]]="Unknown"

ssmat<-table(mixed.plant.sub$pl.p, mixed.plant.sub$species)

# See how number of stems per plot is distributed
hist(mixed.site.sub$plant.count[mixed.site.sub$plant.count<100], breaks=100)

# calculate alpha diversity
alpha<-rowSums(ifelse(ssmat>0,1,0))

# rarefy richness based on 30 plants per plot
rare.rich<-as.data.frame(t(rarefy(x=ssmat, sample=30, se=TRUE)))
rare.rich$alpha<-alpha
rare.rich$pl.p<-rownames(rare.rich)
colnames(rare.rich)<-c("rare.rich", "rare.rich.se", "alpha", "pl.p")

# make sure we're not missing any plots in either dataset
unique(rare.rich$pl.p[is.na(match(rare.rich$pl.p, mixed.site.sub$pl.p))])
unique(mixed.site$pl.p[is.na(match(mixed.site.sub$pl.p, rare.rich$pl.p))])

mixed.site.sub<-merge(mixed.site.sub, rare.rich)
colnames(mixed.site.sub)

#                               EXPLORE RICHNESS ####

# look at correlation between richness, evenness and plant count
with(mixed.site.sub, cor(cbind(alpha, rare.rich, plant.count)))

plot(alpha ~ log(mixed.site.sub$plant.count))
plot(rare.rich$rare.rich ~ log(mixed.site.sub$plant.count))

#             CALCULATE TRAITS ####
#                               CONTINUOUS TRAIT IMPORT #####

# import all lists of trait values


raw.trait<-list(height=read.csv("./Data/max.height.csv", header=T),
                      sla=read.csv("./Data/sla.csv", header=T),
                      seed.mass=read.csv("./Data/seed.mass.csv", header=T),
                      wood.density=read.csv("./Data/wood.density.csv", header=T))

trait.mean.list<-lapply(raw.trait, function(x){
  sapply(split(x[,2], f=x$species), function(y){mean(y, na.rm=TRUE)})
})

trait.means<-data.frame(species=species$species)

trait.mean.data<-sapply(trait.mean.list, function(x){
  temp.trait<-data.frame(species=names(x),
                         trait=x)
  temp<-merge(trait.means, temp.trait, all.x=TRUE)
  return(as.vector(temp[,2]))
  })
trait.means<-cbind(trait.means, trait.mean.data)
colnames(trait.means)<-c("species", "MH","SLA","SM","WD")

# now we want to get average values for the genus only IDs
# let's add genus to our trait list
trait.means<-merge(trait.means, species[,colnames(species) %in% 
                                       c("species", "genus", "genus.only")], 
                  all.x=TRUE)

# Now we can apply these values to the plantscale dataframe
# (giving each individual plant it's relevant trait score)
mixed.plant.sub<-merge(mixed.plant.sub, trait.means[,colnames(trait.means) %in%
                                                    c("species","SLA","WD","SM","MH")], 
                       all.x=TRUE)

# Now I want to set up an indicator for each trait telling me how many species it was calculated from
# 0 means it is a species-level value
trait.n<-as.data.frame(matrix(0, ncol=4, nrow=length(mixed.plant.sub[,1]), 
                        dimnames=list(NULL, paste(c("SLA","WD","SM","MH"), ".n", sep=""))))
# I also want the radius of the circle that was used to calculate genus means
trait.dist<-as.data.frame(matrix(NA, ncol=4, nrow=length(mixed.plant.sub[,1]), 
                              dimnames=list(NULL, paste(c("SLA","WD","SM","MH"), ".dist", sep=""))))

mixed.plant.sub<-cbind(mixed.plant.sub, trait.n, trait.dist)

#                               CALCULATE GENUS MEANS ####

# now here we want to set genus means, but universal means are a bad idea, especially in genera with a lot of
# variation (e.g., Eucalyptus). So what we'll do is a hierarchical process. First, I'm assuming if there's another species
# in the sample plot as a genus ID, the data collector recognised the genus ID was NOT the same species, so I don't want
# to apply the species trait values to the genus ID. Instead what I'll do is apply a mean from a discreet region.

# for starters we want all genus IDs, and species without trait values
genus<-unique(mixed.plant.sub$species[mixed.plant.sub$genus.only==1])

species.missing.trait<-
  unique(mixed.plant.sub$species[
    rowSums(is.na(mixed.plant.sub[,colnames(mixed.plant.sub) %in% 
                                   c("SLA","WD","SM","MH")]))>0 &
  mixed.plant.sub$genus.only==0 &
  mixed.plant.sub$unknown==0])


# next thing we need is a distance matrix for our plantings
grid.plp<-expand.grid(mixed.site.sub$pl.p, mixed.site.sub$pl.p, stringsAsFactors=FALSE)
grid.longlat<-data.frame(long1=mixed.site.sub$long[match(grid.plp[,1], mixed.site.sub$pl.p)],
                         lat1=mixed.site.sub$lat[match(grid.plp[,1], mixed.site.sub$pl.p)],
                         long2=mixed.site.sub$long[match(grid.plp[,2], mixed.site.sub$pl.p)],
                         lat2=mixed.site.sub$lat[match(grid.plp[,2], mixed.site.sub$pl.p)])
head(grid.longlat)

# custom function to turn degrees to radians
deg2rad <- function(deg) {(deg * pi) / (180)} 

R = 6371e3 # earth circumference in metres
lat1<-deg2rad(grid.longlat$lat1) # φ1 (lat of plot 1 in radians)
lat2<-deg2rad(grid.longlat$lat2) # φ2 (lat of plot 2 in radians)
dlat<-with(grid.longlat, deg2rad(lat2 - lat1)) # Δφ (lat difference in radians)
dlong<-with(grid.longlat, deg2rad(long2 - long1)) # Δλ (long difference in radians)

# funky trigonometry (converted from javascript function found at
# http://www.movable-type.co.uk/scripts/latlong.html)
a<-sin(dlat/2) * sin (dlat/2) +
  cos(lat2) * cos (lat2) *
  sin(dlong/2) * sin(dlong/2)
c<- 2 * atan2(sqrt(a), sqrt(1-a))
d <- (R * c) /1000 # distance in km

pl.p.dist<-matrix(d, ncol=length(mixed.site.sub$pl.p), 
                  nrow=length(mixed.site.sub$pl.p),
                  dimnames=list(mixed.site.sub$pl.p, mixed.site.sub$pl.p),
                  byrow=TRUE)

# now for each plot we need to find each species missing a trait value, then attempt to
# find other congeners in a certain radius
species.by.pl.p<-split(mixed.plant.sub, f=mixed.plant.sub$pl.p)

# identify which genus IDs / species are missing traits
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
setDefaultCluster(cl)
# export data and load libraries in cluster
clusterExport(varlist=c("mixed.plant.sub", "pl.p.dist", "species", "species.by.pl.p"))
genus.means.list<-parLapply(cl=cl, 1:length(species.by.pl.p), function(row.num){
  
  ### 1. GET PL.P ###
  pl.p<-species.by.pl.p[[row.num]]
  
  ### 2. GET PL.P CHARACTERISTICS ###
  
  # get species with missing traits
  species.traits<-pl.p[!duplicated(pl.p$species),
                       colnames(pl.p) %in% 
                       c("species", "SLA","WD","SM","MH")]
  missing.traits<-is.na(species.traits[,-1])
  missing.species<-species.traits[rowSums(is.na(species.traits[,-1]))>0,1]

  if(length(missing.species)==0){return(NULL)}
  
  # get distance scores for pl.p
  pl.p.distances<-pl.p.dist[row.num,]

  # identify genus of species
  missing.genus<-species$genus[species$species %in% missing.species]

  # now for each species, attempt to find other members of genus within a given radius
  
    # 3. IDENTIFY TRAIT VALUES FOR EACH GENUS WITHIN GIVEN RADIUS ###

  pl.p.mean.traits<-lapply(missing.genus, function(genus){
    
    candidate.species<-species$species[species$genus==genus]
  
    if(length(candidate.species)==1){return(rep(NA,12))}
    
    # for a number of radii, calculate mean traits and sample size
    cand.traits<-t(sapply(c(50,100,200,500,1000,10000), function(radius){
      
                  # subset pl.ps within radii
                  pl.p.radii<-mixed.plant.sub[
                                       mixed.plant.sub$pl.p %in%
                                       names(pl.p.distances)[pl.p.distances<radius] &
                                       mixed.plant.sub$species %in% candidate.species,]
                  if(dim(pl.p.radii)[1]==0){return(c(rep(NA,4),rep(0,4), rep(NA,4)))} # if no species from genus, return NAs
                  
                  if(dim(pl.p.radii)[1]>0){
                    
                    pl.p.radii<-pl.p.radii[!duplicated(pl.p.radii$species),]
                    
                    trait.means<-colMeans(pl.p.radii[,colnames(pl.p.radii) %in%
                                                      c("SLA","WD","SM","MH")], 
                                          na.rm=TRUE)
                  n<-apply(pl.p.radii[,colnames(pl.p.radii) %in% 
                                        c("SLA","WD","SM","MH")], 
                           MARGIN=2, function(x){length(x[!is.na(x)])})
                  
                    return(c(trait.means, n, rep(radius,4)))
                  }
              }))
    
    # now accept the first non-NA value
    smallest.distance<-apply(cand.traits[,1:4], MARGIN=2, function(x){match("TRUE", !is.na(x))})
    
    if(sum(!is.na(smallest.distance))>0){ # if there's a value for all traits
    
    genus.means<-cand.traits[cbind(smallest.distance, 1:4)]
    genus.n<-cand.traits[cbind(smallest.distance, 5:8)]
    genus.distance<-cand.traits[cbind(smallest.distance, 9:12)]
    return(c(genus.means, genus.n, genus.distance))
    }
    
    if(sum(!is.na(smallest.distance))==0){ # if there's no value for any traits (like a genus with no species)
      return(rep(NA,12))
    }

    })
   
  pl.p.mean.traits.df<-do.call("rbind", pl.p.mean.traits)
  pl.p.mean.traits.df<-as.data.frame(matrix(pl.p.mean.traits.df, ncol=12))
  colnames(pl.p.mean.traits.df)<-c("SLA","WD","SM","MH",
                                   "SLA.n","WD.n","SM.n","MH.n",
                                   "SLA.rad","WD.rad","SM.rad","MH.rad")
  rownames(pl.p.mean.traits.df)<-missing.species 
  return(pl.p.mean.traits.df)
  
    })
stopCluster(cl=cl)

# now we need to overwrite each pl.p's genus traits with new values
# except where traits already exist (for some species IDs)
pl.p.missing.traits<-mapply(pl.p=species.by.pl.p, 
                            trait.means=genus.means.list, 
                            function(pl.p, trait.means){
                              
                              if(is.null(trait.means)){return(NULL)}
                              
                              species<-split(pl.p, f=droplevels(pl.p$species))
                              
                              # first match up the rows of the pl.p to the missing traits
                              trait.match<-match(rownames(trait.means), names(species))
                              
                              # now we need to know whether any of the pl.p species 
                              # already have species level traits
                              trait.cols<-match(c("SLA","WD","SM","MH"), colnames(pl.p))
                              
                              # now we need to override trait values that are NA
                              species.missing.trait<-mapply(missing.species=species[trait.match], 
                                                            traits=split(trait.means, rownames(trait.means)), 
                                                            function(missing.species, traits){
                                                              
                                                              traits.df<-traits[rep(seq_len(nrow(traits)), 
                                                                                    each=length(missing.species[,1])),]
                                                              
                                                              #over-write trait values
                                                              missing.traits<-is.na(missing.species[1,
                                                                                                    grepl("SLA|WD|SM|MH", 
                                                                                                          colnames(missing.species))][,1:4])
                                                              
                                                              trait.cols<-colnames(missing.species) %in% colnames(missing.traits)
                                                              n.cols<-grepl("\\.n", colnames(missing.species))
                                                              dist.cols<-grepl("\\.dist", colnames(missing.species))
                                                              
                                                              missing.species[,trait.cols][,missing.traits]=traits.df[,1:4][,missing.traits]
                                                              missing.species[,n.cols][,missing.traits]=traits.df[,5:8][,missing.traits]
                                                              missing.species[,dist.cols][,missing.traits]=traits.df[,9:12][,missing.traits]
                                                              
                                                              return(missing.species)  
                                                            }, SIMPLIFY=FALSE)
                              
                              
                              pl.p.new<-rbind(do.call("rbind", species[-trait.match]), 
                                              do.call("rbind", species.missing.trait))
                              return(pl.p.new)
                            })

# now we need to combine our list (excluding the Nulls, with the pl.ps with no missing traits)
# first let's get the NULLs out
null.pl.p<-sapply(pl.p.missing.traits, 
                  function(x){if(is.null(x)){return(1)} else{return(0)}})

null.pl.p.names<-names(null.pl.p[null.pl.p==1])

mixed.plant.sub<-rbind(mixed.plant.sub[mixed.plant.sub$pl.p %in% 
                                       null.pl.p.names,],
                        do.call("rbind", pl.p.missing.traits[null.pl.p==0]))

#                               CATEGORICAL TRAIT IMPORT ####

# I've got a list of which species are nitrogen fixers, and categorised
# species as either shrubs (max height of <6m) and trees (max height >6m)
# so we can calculate a proportion of each of these for each plot and 
# analyse them using a logistic regression

cat.trait.species<-species[,colnames(species) %in% c("species", "n.fixer", "tree")]

mixed.plant.sub<-merge(mixed.plant.sub, cat.trait.species)

cat.traits.pl.p<-as.data.frame(t(sapply(split(mixed.plant.sub, mixed.plant.sub$pl.p), 
                        
                        function(x){
                          c(mean(x$n.fixer, na.rm=TRUE),
                            mean(x$tree, na.rm=TRUE))
                        }
                        )))
  
colnames(cat.traits.pl.p)<-c("nfixer.prop", "tree.prop")
cat.traits.pl.p$pl.p<-rownames(cat.traits.pl.p)

head(cat.traits.pl.p)

mixed.site.sub<-merge(mixed.site.sub, cat.traits.pl.p)

#                               TRANSFORM TRAITS ####

par(mfrow=c(2,2))
sapply(c("SLA","WD","SM","MH"), function(x){
  hist(unique(mixed.plant.sub[,colnames(mixed.plant.sub) %in% x]), 
       xlab="", main="")
  mtext(side=3, text=x)
})

# traits are pretty skewed, especially seed mass. We'll need to log-transform first
# otherwise our standardisation will be off by orders of magnitude.
log.traits<-sapply(c("SLA","WD","SM","MH"), function(x){
  log(mixed.plant.sub[,colnames(mixed.plant.sub) %in% x])
})

par(mfrow=c(2,2))
sapply(c("SLA","WD","SM","MH"), function(x){
  hist(unique(log.traits[,colnames(log.traits) %in% x]), 
       xlab="", main="")
  mtext(side=3, text=x)
})

# save our raw traits
mixed.plant.sub.raw<-mixed.plant.sub

# override raw traits with standardised values
mixed.plant.sub[,match(colnames(log.traits), colnames(mixed.plant.sub))]=
  log.traits

#             CALCULATE UNIVARIATE FUNCTIONAL DIVERSITY SCORES ####

# All of these come from "A user's guide to functional diversity indices"
# by D. Schleuter, M. Daufresne, F. Massol, C. Argillier. 2010. Ecological Monographs

#                               COMMUNITY-WEIGHTED MEANS ####
traits<-c("SLA","WD","SM","MH")

fr.pl.p.list<-split(mixed.plant.sub, mixed.plant.sub$pl.p)
fr.pl.p.traits<-lapply(fr.pl.p.list, function(x){x[,colnames(x) %in% traits]})

cwm<-t(sapply(fr.pl.p.traits, function(x){colMeans(x, na.rm=TRUE)}))
cwm<-as.data.frame(cwm)
colnames(cwm)<-c("SLA.cwm", "WD.cwm", "SM.cwm", "MH.cwm")
cwm$pl.p<-rownames(cwm)

mixed.site.sub<-merge(mixed.site.sub, cwm, all.x=TRUE)

#                               FUNCTIONAL RANGE (Mason et al 2005, FRr) ####

# "FR indices measure how much of the niche space is occupied by the species 
# present. They are usually interpreted by ecologists as an indicator for 
# potentially used/unused niche space and thus e.g. for productivity, buffering 
# against environmental fluctuations or vulnerability to invasion (Mason et al. 
# 2005). FR is naturally positively correlated to the number of species present
# (the more species there are, the larger the functional space occupied when 
# species traits are somewhat randomly distributed). However, two communities 
# with the same number of species may have different FR when functional traits 
# of species are more closely clustered in one community than in the other. FR 
# is not weighted by species abundance" (Schleuter et al, 2010).

# my first choice would be to use Schleuter's index, as it factors in gaps in 
# the trait space, but it requires intra-specific variation, which I don't have
# enough of to use. Instead I'll calculate FRr using Mason's protocol:
# calculated by subtracting the largest trait value in the community by the 
# smallest.

# I'm not dividing by the entire species list because my traits will be 
# standardised prior to modelling

# returns lists nestled within a dataframe, one for each plot
species.by.pl.p<-split(x=mixed.plant.sub[,colnames(mixed.plant.sub) %in% 
                                           c("species", traits)],         
                               f=mixed.plant.sub$pl.p)

# remove duplicate entries in each list to get unique species
unique.species.by.pl.p<-lapply(species.by.pl.p, 
                               function(x){x[!duplicated(x),]})

# make species the rownames of the dfs in the list
unique.species.by.pl.p<-lapply(unique.species.by.pl.p, 
                               function(x){rownames(x)<-x[,1]; x<-x[,-1]})

# calculate min and max values for each plot

FRr<-as.data.frame(t(sapply(unique.species.by.pl.p, function(x){
           min<-apply(x, 2, FUN=function(y){min(y, na.rm=TRUE)})
           max<-apply(x, 2, FUN=function(y){max(y, na.rm=TRUE)})
           
           return(max-min)
           })))

# Score of 0 indicate monocultures, as the max and min are the same (difference 
# of 0 divided by overall trait ranges)
colnames(FRr)<-paste0(colnames(FRr),".FRr")
FRr$pl.p<-rownames(FRr)

# merge these values into our site-level dataframe
mixed.site.sub<-merge(mixed.site.sub, FRr, all.x=TRUE)
head(mixed.site.sub)

#                               FUNCTIONAL EVENNESS (Mouillot, FEs) ####

# FE indices (AKA functional regularity) measure whether mean species traits are distributed 
# regularly within the occupied trait space, i.e. with equal distances between nearest neighbors 
# and equal abundances (a high FE index usually means a very regular distribution; a low FE index, 
# the existence of separate clouds of species and/or abundances).

# Actually called 

# modified from http://villeger.sebastien.free.fr/Rscripts.html
# proposed by Mouillot et al. 2005 (Oecologia 142: 353-359)                               
# and modified in Villeger et al., 2008, Marine Ecology Progress Series (364: 135-146) 
              
# essentially you capture the pair-wise trait differences, and divide the sum of these differences
# by 1/(richness-1)

# TRAIT LIST #

# First we need to get a list with all communities #
fr.pl.p.list<-split(mixed.plant.sub, mixed.plant.sub$pl.p)

# extract all the trait values, as well as the abundance of each species
fr.trait.list<-lapply(fr.pl.p.list, function(x){
  pl.p.trait<-x[!duplicated(x[colnames(x)=="species"]),colnames(x) %in% c("species", traits)]
  pl.p.count<-as.matrix(table(droplevels(x[colnames(x)=="species"])))
  pl.p.trait$count<-pl.p.count[match(pl.p.trait$species, rownames(pl.p.count))]
  return(pl.p.trait)
  })

# calculate the FEs score for each trait (function)
FEs.fun<-function(trait, count){
  
  pairwise.temp<-data.frame(trait=trait, count=count)[!is.na(trait),]
  pairwise.temp<-pairwise.temp[order(pairwise.temp$trait),]
    
  if(length(pairwise.temp$trait)>1){  
    
    tr<-pairwise.temp[,1]
    ab<-pairwise.temp[,2]
    
    # filter to keep only species present
    pres<-which(ab>0) ; abF<-ab[pres] ; trF<-tr[pres]  
    
    # number of species present
    s<-length(abF) ; os<-1/(s-1)
    
    # sorted traits values and relative abundances of species present 
    o<-order(trF,decreasing=F) ; to<-trF[o] ; abo<-abF[o]/sum(abF)
    
    # computation of index
    EW<- abs(to[-1]-to[-s]) / (abo[-1]+abo[-s])
    PEW<-EW/sum(EW)
    minPEW<-sapply(PEW, function(x) { min(x,os) } )
    
    if(length(minPEW)==1){return(NA)} else {return( (sum(minPEW)-os) / (1-os))}
     # in original formula FRo=sum(minPEW)
  }
  
  if(length(pairwise.temp$trait)<=1){
    return(NA)
  }
}

# run function on each plot for each trait
FEs<-as.data.frame(t(sapply(fr.trait.list, function(x){
  pair.temp<-as.data.frame(x)
  return(c(FEs.fun(pair.temp$SLA, pair.temp$count),
           FEs.fun(pair.temp$WD, pair.temp$count),
           FEs.fun(pair.temp$SM, pair.temp$count),
           FEs.fun(pair.temp$MH, pair.temp$count)))})))

colnames(FEs)<-c("SLA.FEs","WD.FEs","SM.FEs","MH.FEs")

FEs$pl.p<-rownames(FEs)

mixed.site.sub<-merge(mixed.site.sub, FEs, all.x=TRUE)

#                               FUNCTIONAL DIVERGENCE (Schleuter et al 2010, FDs) ####

# "FD indices, measure the variance of the species' functions and the position of their
# clusters in trait space (a high FD is caused by the clustering of species and/or abundances at
# the edges of the traits' space)" (Schleuter et al 2010).

# calculate these for each planting - returns dataframe with 5 columns
FDs.data<-lapply(species.by.pl.p, 
                   function(x){
                  
                      apply(x[,-1], MARGIN=2, 
                           FUN=function(y){
                           temp<-quantile(y, 
                                          probs=c(0,0.25,0.75,1), 
                                          na.rm=TRUE)
                           return((temp[3]-temp[2])/
                                  (temp[4]-temp[1]))
                           })
                      })
FDs<-as.data.frame(do.call("rbind", FDs.data))
FDs[is.na(FDs)]=0
                            
FDs$pl.p<-names(FDs.data)
colnames(FDs)[-5]<-paste(c("MH","SLA","SM","WD"), ".FDs", sep="")

mixed.site.sub<-merge(mixed.site.sub, FDs, all.x=TRUE)

#                               EXPLORE UNIVARIATE FD ####

# the most important thing to look for is correlation between trait measures
# and with richness, and plant number

trait.cormat<-cor(mixed.site.sub[,grepl("SLA|WD|SM|MH", colnames(mixed.site.sub))
                                  | colnames(mixed.site.sub) %in%
                                    c("rare.rich", "plant.count"),], 
    use="complete.obs")

ifelse(abs(trait.cormat)>0.4,trait.cormat,0)

#             CALCULATE MULTIVARIATE FUNCTIONAL DIVERSITY SCORES ####

# So far all of these come from "A user's guide to functional diversity indices"
# by D. Schleuter, M. Daufresne, F. Massol, C. Argillier. 2010. Ecological Monographs

# for these metrics we will need to work on the subset of plantings with 
# complete trait values
mixed.plant.sub$n.traits<-
  rowSums(!is.na(
    mixed.plant.sub[,colnames(mixed.plant.sub) %in% c("SLA","WD","SM","MH")]))

# get a logical vector where the plot has all traits for >90% of species
mixed.plant.sub.90<-sapply(split(mixed.plant.sub, f=mixed.plant.sub$pl.p),
                           function(x){
                             y<-x$n.traits
                             length(y[y==4])/length(y)})

# subset to sites with >90% trait coverage
pl.p90<-names(mixed.plant.sub.90[mixed.plant.sub.90>0.9])

# we used to subset out plots without complete trait coverage, but now I have
# assigned genus means we can use all plots, just subsetting out unknown or dead
# stems within each plot

# we now need to subset NAs out of the plant data in order to calculate multi-trait FD
mixed.plant.sub<-mixed.plant.sub[
  rowSums(is.na(mixed.plant.sub[match(c("SLA", "WD", "SM", "MH"),
                                         colnames(mixed.plant.sub))]))==0,]

# split up our plant-scale data frame by plot
traits.by.pl.p<-split(x=mixed.plant.sub[,!is.na(match(colnames(mixed.plant.sub), c("species", traits)))],
                      f=mixed.plant.sub$pl.p,
                      drop=TRUE)

abundance.by.pl.p<-lapply(traits.by.pl.p, 
                          function(x){
                            temp<-tapply(x$species,x$species, length)
                            temp.nona<-temp[!is.na(temp)]
                            return(matrix(temp.nona, byrow=TRUE, ncol=length(temp.nona),
                                         dimnames=list(NULL, names(temp.nona))))
                            })

# get rid of duplicates so we only have unique trait values
unique.pl.p.traits<-lapply(traits.by.pl.p, FUN=function(x){x[!duplicated(x),]})

# extract species names for use later
pl.p.traits.species<-lapply(unique.pl.p.traits, FUN=function(x){as.character(x[,1])})

# turn the species names into rownames (the FD package requires this)
unique.pl.p.traits<-lapply(unique.pl.p.traits, FUN=function(x){rownames(x)<-x[,1]; x<-x[,-1]})

# order both lists by species name so things match up
abundance.by.pl.p<-mapply(x=abundance.by.pl.p, y=pl.p.traits.species, function(x,y){
                                                                temp<-x[order(y)]
                                                                attr(temp, "names")=y[order(y)]
                                                                return(temp)})

unique.pl.p.traits<-lapply(unique.pl.p.traits, function(x){x[order(rownames(x)),]})

# we can use the FD package to calculate these, but it can't handle monoculture communities.
# I need to remember that the plot names are contained only in the abundance list, not the trait list
mixed.comms<-sapply(abundance.by.pl.p, function(x){length(x)>1})

abundance.by.pl.p.multisp<-abundance.by.pl.p[mixed.comms]
unique.pl.p.traits.multisp<-unique.pl.p.traits[mixed.comms]

#                               FUNCTIONAL RANGE (Villeger et al 2008, FRv) ####

# Schlueter et al functional range index needs intraspecific variation, which we don't have.
# we'll have to use the old convex-hull method of Villeger.

# this means we can only use plantings with more than 1 species
FRv<-mapply(x=unique.pl.p.traits.multisp, y=abundance.by.pl.p.multisp, 
            FUN=function(x,y){dbFD(x, a=y, w.abun=TRUE, calc.FRic=TRUE, stand.FRic=FALSE,
                                   calc.FGR=FALSE, calc.CWM=FALSE, calc.FDiv=FALSE, 
                                   messages=FALSE)$FRic})

FRv.df<-data.frame(pl.p=substr(names(FRv),1,regexpr("Community1", names(FRv))-2),
                    FRv=FRv)

mixed.site.sub<-merge(mixed.site.sub, FRv.df, all.x=TRUE)

# NAs represent monoculture plots which have a functional range of 0
mixed.site.sub$FRv[is.na(mixed.site.sub$FRv)]<-0

#                               FUNCTIONAL EVENNESS (Villeger et al 2008, FEm) #### 

FEm<-mapply(x=unique.pl.p.traits.multisp, y=abundance.by.pl.p.multisp, 
            FUN=function(x,y){dbFD(x, a=y, w.abun=TRUE, calc.FRic=FALSE, 
                                   calc.FGR=FALSE, calc.CWM=FALSE, calc.FDiv=FALSE, 
                                   messages=FALSE)$FEve})

FEm.df<-data.frame(pl.p=substr(names(FEm),1,regexpr("Community1", names(FEm))-2),
                FEm=FEm)

mixed.site.sub<-merge(mixed.site.sub, FEm.df, all.x=TRUE)

#                               FUNCTIONAL DIVERGENCE (Villeger et al 2008, FDm) #### 

# we can use the FD package to calculate this
?dbFD
FDm<-mapply(x=unique.pl.p.traits.multisp, y=abundance.by.pl.p.multisp, 
            FUN=function(x,y){dbFD(x, a=y, w.abun=TRUE, messages=FALSE)$FDiv})

FDm<-do.call("c", FDm)

FDm.df<-data.frame(pl.p=substr(names(FDm),1,regexpr("Community1", names(FDm))-2),
                    FDm=FDm)

mixed.site.sub<-merge(mixed.site.sub, FDm.df, all.x=TRUE)

# NAs represent monoculture plots which have a functional range of 0
mixed.site.sub$FDm[is.na(mixed.site.sub$FDm)]<-0
mixed.site.sub$FEm[is.na(mixed.site.sub$FEm)]<-0

#                               EXPLORE MULTIVARIATE FD ####

multi.fd.cor<-cor(mixed.site.sub[,colnames(mixed.site.sub) %in% 
                                           c("FRv", "FEm", "FDm", "rare.rich", "plant.count")], 
                  use="complete.obs")
ifelse(abs(multi.fd.cor)>0.10, multi.fd.cor, 0)

#             STANDARDISE ENVIRONMENTAL VARIABLES ####

mixed.site.sub$mass.per.year<-log(mixed.site.sub$biomass.area / mixed.site.sub$age)
# We need to look at regression slopes (and diagnostic plots) for each of our 
# predictors to work out if we need to transform them.

# We've already tranformed the traits, so we'll include functional diversity
# indices but we won't be able to transform them again.

# start by getting a vector of variables.
shape.files<-list.files(path=paste(ifelse(Sys.info()['sysname']=="Linux", 
                                          "/home/Storage HDD/", "F:/"),
                                   "University files/Shape files", sep=""), 
                        pattern="\\.tif")

predictors<-c(substr(shape.files,
                1, regexpr("\\.tif", shape.files)-1),
           "age","density","rare.rich")

names<-c("Moisture availability", "Carbon", "Clay", "Elevation relief",
         "Elevation", "Nitrogen", "Phosphorus", "pH", "Annual rainfall",
         "Rainfall seasonality", "Prescott index", "Sand", "Silt", "Slope",
         "Soil depth", "Solar radiation", "Max temp", "Min temp", "TWI", 
         "Age", "Density", "Richness")

predictors<-predictors[c(1,20,21)]
names<-names[c(1,20,21)]

pdf(paste0("./Plots/Transformation check ", Sys.Date(), ".pdf"), 
    height=6, width=4, useDingbats=FALSE)

par(mar=c(2,1,0,0), mfrow=c(3,2), ps=8, tck=-0.025, mgp=c(3,0,0))
#sapply(predictors, function(x){
x<-1
for(i in 1:length(predictors)){
  print(i)
  x<-predictors[i]
  model.data<-mixed.site.sub
  temp.pred<-model.data[,colnames(model.data) %in% x]
  
  # RAW MODEL #
  raw.lme<-lme(as.formula(paste("mass.per.year ~", x)),
               random=~1|IBRAsub/planting,
               data=model.data,
               na.action="na.omit")
  
  raw.gamm<-gamm(as.formula(paste0("mass.per.year ~ s(",x,", bs=\"cr\")")),
                 random=list(IBRAsub=~1, planting=~1),
                 data=model.data)
  

  # PLOTS 
  par(mar=c(3,3,1,0), ps=8, tck=-0.025, mgp=c(3,0.4,0))
  
  with(model.data, plot(mass.per.year ~ temp.pred,
       xlab=x, las=1, col="grey50"), axes=FALSE)
  
  axis(side=1)
  axis(side=2, mgp=c(3,0.4,0))
  
  new.dat<-data.frame(y<-seq(min(temp.pred, na.rm=TRUE), 
                             max(temp.pred, na.rm=TRUE), length.out=100))
  colnames(new.dat)<-x
  lme.pred<-predict(raw.lme, newdata=new.dat, level=0)
  gam.pred <- predict(raw.gamm, newdata = new.dat, type = "response")
  points(lme.pred ~ new.dat[,1], type="l", lwd=1.5)
  points(gam.pred ~ new.dat[,1], type = "l", col="red", lwd=1.5)
  mtext(side=3, text=paste0(names[i]))
  
  if(grepl("SLA|WD|SM|MH", x)){next}
  
   # LOG MODEL
  
  model.data[,colnames(model.data) %in% x]<-
    log(mixed.site.sub[,colnames(mixed.site.sub) %in% x])
  temp.pred<-model.data[,colnames(model.data) %in% x]
  
  log.lme<-lme(as.formula(paste("mass.per.year ~", x)),
                random=~1|IBRAsub/planting,
                data=model.data,
                na.action="na.omit")
  
  log.gamm<-gamm(as.formula(paste0("mass.per.year ~ s(",x,", bs=\"cr\")")),
                  random=list(IBRAsub=~1, planting=~1),
                  data=model.data)
  
  with(model.data, plot(mass.per.year ~ temp.pred,
                        xlab=x, las=1, col="grey50"), axes=FALSE)
  
  axis(side=1)
  axis(side=2, mgp=c(3,0.4,0))
  
  new.dat<-data.frame(y<-seq(min(temp.pred, na.rm=TRUE), 
                             max(temp.pred, na.rm=TRUE), length.out=100))
  colnames(new.dat)<-x
  lme.pred<-predict(log.lme, newdata=new.dat, level=0)
  gam.pred <- predict(log.gamm, newdata = new.dat, type = "response")
  
  points(lme.pred ~ new.dat[,1], type="l", lwd=1.5)
  points(gam.pred ~ new.dat[,1], type = "l", col="red", lwd=1.5)
  mtext(side=3, text=paste0("ln(",names[i],")"))
  
 }
dev.off()

log.predictor.names<-c("aridity.index", "prec.annual", "prescott", "age",
                       "density", "arid.man", "arid.monthmean")

log.predictors<-log(mixed.site.sub[,colnames(mixed.site.sub) %in% 
                                     log.predictor.names])

colnames(log.predictors)<-paste0("log.",colnames(log.predictors))

mixed.site.sub<-cbind(mixed.site.sub, log.predictors)

all.predictors<-c(colnames(log.predictors), predictors[!predictors %in% 
                                                          log.predictor.names])
# Remove diversity indices
covariates<-all.predictors[!grepl("SLA|WD|SM|MH|rare.rich", 
                                  all.predictors) &
                           !all.predictors %in% c("FRv","FEm","FDm")]

#                               EXPLORE ENVIRONMENTAL VARIABLES ####

#                            5.1.2 PAIRS PLOT #

# so let's pairs plot our environmental variables

# subset environmental covariates
covariate.names.sorted<-covariates[c(1,5, # structure
                                     2,3,4,13, # moisture
                                     6,10,11,12, # soil chemistry
                                     7,14,15, # soil texture
                                     18:20, # solar energy
                                     8,9,16,17,21)] # topography

head(mixed.site.sub.env)

mixed.site.sub.env<-mixed.site.sub[, c("log.age", "log.density", "log.aridity.index",
                                       "log.prec.annual", "log.prescott", "prec.season",
                                       "carbon", "nitrogen", "phosphorus", "ph",
                                       "clay", "sand", "silt", "solar.annual",
                                       "tmax.annual", "tmin.annual", 
                                       "elevation.relief", "elevation", "slope",
                                       "soil.depth", "TWI")]

head(mixed.site.sub.env)

colnames(mixed.site.sub.env)[colnames(mixed.site.sub.env)=="phosphorus"]<-"soil.p"
colnames(mixed.site.sub.env)[colnames(mixed.site.sub.env)=="elevation.relief"]<-"relief"

nrow<-length(c(1:length(colnames(mixed.site.sub.env))))
comb<-expand.grid(c(1:nrow),c(1:nrow))
upper.tri<-comb[,2]<comb[,1]
env.pairs<-apply(comb, MARGIN=1, FUN=function(x){mixed.site.sub.env[,x]})
env.names<-cbind(colnames(mixed.site.sub.env),
                 c("Age", "Density", 
                   "Moisture\navailability", "Rainfall", "Prescott", "Rainfall\nSeasonality",
                   "C", "N", "P", "pH",
                   "Clay", "Sand", "Silt", 
                   "Solar", "Mean max\n temperature", "Mean min\n temperature", 
                   "Elevation\nrelief", "Elevation", "Slope", "Soil depth", "TWI"))

pdf(paste("./Plots/Env pairs ", Sys.Date(), ".pdf"), height=10, width=14, useDingbats=FALSE)
par(mfrow=c(nrow,nrow), mar=c(0,0,0,0), oma=c(3,3,6,3))
# run pairs plot code

log<-mapply(x=env.pairs, y=upper.tri, FUN=function(x,y){
  
  # HISTOGRAM ON DIAGONAL
  
  if(grepl(colnames(x)[1], colnames(x)[2])) {
    hist(x[,1], col="grey80", axes=F, main=NULL, border="grey50")
    box()
    text(x=par("usr")[1]+0.5*(diff(par("usr")[c(1:2)])),y=par("usr")[3]+0.5*(diff(par("usr")[c(3:4)])), 
         labels=env.names[match(colnames(x)[1], env.names[,1]),2], font=2)
    
  }
  
  # SCATTER PLOT ON LOWER TRIANGULAR
  if(!grepl(colnames(x)[1], colnames(x)[2]) & !y){
    plot(x, axes=F, col=rgb(0.8,0.8,0.8,0.2), pch=16)
    box()
  }
  
  # CORRELATION COEFFICIENT ON UPPER TRIANGULAR
  if(y){
    plot(x=NULL,y=NULL,xlim=c(0,1),ylim=c(0,1), axes=F)
    
    text(x=0.5, y=0.5,labels=round(cor(x, use="complete.obs")[2,1],3), adj=0.5, 
         cex=(1.5^(abs(cor(x, use="complete.obs")[2,1])*3))/2)
    box()
    
  }})
dev.off()

#                            5.1.2 ENVIRONMENTAL CORRELATIONS #

env.cormat<-cor(mixed.site.sub.env, use="complete.obs")
ifelse(abs(env.cormat)>0.6, env.cormat, 0)

# lots of strong correlations, as expected - can we mitigate some of them
# by transforming groups of environmental variables onto PCA axes?

write.csv(mixed.site.sub, "./Data/mixed site sub.csv")
write.csv(mixed.plant.sub, "./Data/mixed plant sub.csv")
write.csv(mixed.plant.sub.raw, "./Data/mixed plant sub (raw traits).csv")

# ####
# IMPORT PROCESSED DATA ####
  mixed.site.sub<-read.csv("./Data/mixed site sub.csv")
  mixed.plant.sub<-read.csv("./Data/mixed plant sub.csv")
  mixed.plant.sub.raw<-read.csv("./Data/mixed plant sub (raw traits).csv")
# ####
#  MODEL PREP ####
#             VARIABLE NAME SORTING ####

log.predictors<-colnames(mixed.site.sub)[grepl("log", colnames(mixed.site.sub))]

covariates<-c(raw.predictors[!grepl(paste(substr(log.predictors, 5,100), collapse="|"), raw.predictors)],
              log.predictors)

colnames(mixed.site.sub)
covariates<-c("carbon", "clay", "elevation.relief", "elevation", "nitrogen", "phosphorus",
              "ph", "log.prec.annual", "prec.season", "sand", "silt", "slope", "soil.depth",
              "solar.annual","tmax.annual","tmin.annual","TWI","log.aridity.index",
              "log.age","log.density")
response.variable<-"biomass.area"
spatial.variable<-c("IBRA","IBRAsub","planting","plot")

# Fixed effects
primary.variable<-c("rare.rich", 
                    colnames(mixed.site.sub)[grepl("SLA|WD|SM|MH", colnames(mixed.site.sub))],
                    "FRic","FEm","RaoQ")

all.predictors<-c(covariates, primary.variable)

#             DETERMINE BEST MODEL TYPE ####
#             ESTABLISHING BEST COVARIATES ####
max.cor<-0.6
model.random.effects<-"IBRAsub/planting"

if(redo.best.covariates==TRUE){

#                           1.2.2 FULL SET WITH EXCLUSION MATRIX #

model.data<-mixed.site.sub[, colnames(mixed.site.sub) %in% 
                             c("biomass.area", "age", "pl.p","planting","IBRAsub","IBRA", covariates)]
model.data<-model.data[rowSums(is.na(model.data))==0,]

model.data$mass.per.year<-log(model.data$biomass.area / model.data$age)

center.data.output<-sapply(model.data[,covariates], 
                           function(x){
                             (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
                             })
summary(center.data.output)

model.data<-as.data.frame(cbind(center.data.output,
                                model.data[!colnames(model.data) %in% covariates]))

# prescott is tricky to explain. Let's see how things look without it.
model.data<-model.data[,!colnames(model.data) %in% c("log.prescott")]
dredge.covariates<-covariates[!covariates %in% c("log.prescott")]

# Now let's try a larger dredge with all the environmental variables, excluding 
# the ones with correlations over 0.6 (or -0.6)
# Let's also get a correlation matrix with 1s and 0s so we can exclude particular 
#combinations of covariates

cov.mat<-cor(model.data[,colnames(model.data) %in% dredge.covariates], 
             use="complete.obs")
cov.mat<-ifelse(upper.tri(cov.mat, diag=TRUE),NA,cov.mat)
cov.mat<-ifelse(abs(cov.mat)>max.cor,FALSE,TRUE)
colnames(cov.mat)<-colnames(model.data[,colnames(model.data) %in% dredge.covariates])
rownames(cov.mat)<-colnames(model.data[,colnames(model.data) %in% dredge.covariates])

env.covariates<-dredge.covariates[!dredge.covariates %in% c("log.age", "log.density")]

m1<-lme(as.formula(paste0("mass.per.year ~",
                          paste(dredge.covariates, collapse=" + "))),
        data=model.data,
        random=~1|IBRAsub/planting)
summary(m1)

 # set up clusters
  no_cores <- detectCores() - 4
  cl <- makeCluster(no_cores)
  setDefaultCluster(cl)
  # export data and load libraries in cluster
  clusterExport(varlist=c("model.data", "cov.mat"))
  clusterCall(cl, "library", "nlme", "MuMIn", character.only = TRUE)
  options(na.action = "na.fail")
  # run dredge
  start<-Sys.time()
  m1.full.pd<-pdredge(m1,cluster=cl, subset=cov.mat, m.lim=c(1,8),
                      fixed=c("log.age", "log.density"))
  stopCluster(cl=cl)
  saveRDS(m1.full.pd, "./Outputs/Covariate dredge.rds")
  finish<-Sys.time()
  print(finish-start)
}

if(redo.best.covariates==FALSE){
  
  m1.full.pd<-readRDS("./Outputs/Covariate dredge.rds")
  
}

#                           1.2.4 MODEL COMPARISON #

dredge.covariates<-covariates[!covariates %in% c("log.prescott")]


model.data<-mixed.site.sub[, colnames(mixed.site.sub) %in% 
                             c("biomass.area", "age", "pl.p","planting","IBRAsub","IBRA", covariates)]
model.data<-model.data[rowSums(is.na(model.data))==0,]

model.data$mass.per.year<-log(model.data$biomass.area / model.data$age)


center.data.output<-sapply(model.data[,covariates], 
                           function(x){
                             (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
                           })
summary(center.data.output)

model.data<-as.data.frame(cbind(center.data.output,
                                model.data[!colnames(model.data) %in% covariates]))


# so now we have a relative measure of the model fit, we need to evaluate this in an absolute way,
# using the R2 code loaded under section #0, for the first 10 models of each dredge
full.candidates<-get.models(m1.full.pd, subset=c(1:10))
# this returns a list 10 long with the best models, which we can evaluate using R2...
model.fits<-as.data.frame(t(sapply(full.candidates, r.squaredGLMM)))
model.fits$number<-1:10

model.fits<-model.fits[order(model.fits[,1], decreasing=TRUE),]
# the AIC columns in this combined data frame are useless as they're relative to each model type, but the Marginal and Conditional
# columns give an R^2 value for the fixed effects and fixed + random effects respectively.
best.model<-full.candidates[[model.fits$number[1]]]
best.covariates<-names(best.model$fixDF$terms)[-1]

#             MEAN CLUSTERING OF ENVIRONMENT ####
# set random seed so mean clustering is always the same (otherwise the same clusters get
# given different cluster IDs - makes it harder to plot)
seed<-12345

# Let's look at means clustering so we can have separate slopes for our diversity variables
# based on environmental conditions
colnames(mixed.site.sub)

best.covariates <- c("log.age", "log.density", "log.aridity.index", "tmin.annual", "sand")

best.env<-mixed.site.sub[, colnames(mixed.site.sub) %in% best.covariates & 
  !grepl("age|density", colnames(mixed.site.sub))]

best.env<-sapply(best.env,
                 function(x){
                   (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
                   })

# set up ss with 1 cluster, and the within group ss with various cluster numbers
wss.best<-c((nrow(best.env)-1)*sum(apply(best.env,2,var)), 
            sapply(2:30, function(x){sum(kmeans(best.env,centers=x)$withinss)}))

plot(1:30, wss.best, type="b", 
     xlab="Number of Clusters", 
     ylab="Within groups sum of squares")

# 10 clusters seems optimal
set.seed(seed)
env.clust.best<-kmeans(best.env, 8, iter.max=100, nstart=25)

mixed.site.sub$best.cluster<-env.clust.best$cluster

with(mixed.site.sub, plot(lat ~ long, col=best.cluster))
table(mixed.site.sub$best.cluster)

write.csv(mixed.site.sub[,c("pl.p","best.cluster")], "./Outputs/pl.p cluster.csv")

# ####
# MODELLING   ####
#             CORRECT BIOMASS FOR AGE ####

colnames(mixed.site.sub)

# raw biomass ha-1
biomass.area<-log(mixed.site.sub$biomass.area)

# correct for age then transform: ln(kg ha-1 year-1)
log.after<-log(mixed.site.sub$biomass.area / mixed.site.sub$age)

# transform biomass then correct for age: ln( ha-1) year-1
log.before<-log(mixed.site.sub$biomass.area) / mixed.site.sub$age

par(mfrow=c(2,3))
plot(log.after ~ biomass.area)
plot(log.before ~ biomass.area)
plot(log.after ~ log.before)

hist(log.after)
hist(log.before)
hist(biomass.area)

# transform biomass into t ha-1 year-1 rather than kg for ease of plotting
mixed.site.sub$mass.per.year<-log((mixed.site.sub$biomass.area/1000) / mixed.site.sub$age)

plot(mixed.site.sub$biomass.area ~ mixed.site.sub$age)

pdf("./Plots/mass per year biomass comparison.pdf",
    height=4, width=4)
par(mar=c(3,3,1,1), ps=8, mgp=c(3,0.5,0), tck=-0.015, las=1)
plot(mixed.site.sub$mass.per.year ~ log(mixed.site.sub$biomass.area),
     xlab="", ylab="")
mtext(side=1, text=expression("Stand biomass ln(kg ha"^-1*")"), line=1.5)
mtext(side=2, text=expression("Stand biomass accrual ln(kg ha"^-1*" year"^-1*")"),
    line=1, las=0)

dev.off()

#             REMOVE MONOCULTURES ####
# ### ### ### ###

# remove monocultures
mixed.site.sub.multi<-droplevels(mixed.site.sub[mixed.site.sub$alpha>1,])

# anonymise site names
mixed.site.sub.multi$planting <- as.numeric(mixed.site.sub.multi$planting)
mixed.site.sub.multi$pl.p <- as.numeric(mixed.site.sub.multi$pl.p)

write.csv(mixed.site.sub.multi, "./Outputs/plot.level.modelling.data.csv")

# 0. INTERCEPT ONLY MODEL - NO MONOCULTURES ####
model.data<-mixed.site.sub.multi[,c(best.covariates, "mass.per.year",
                                    "IBRAsub","planting", "best.cluster",
                                    "rare.rich", "SLA.cwm", "WD.cwm", "SM.cwm",
                                    "MH.cwm")]

model.data<-model.data[complete.cases(model.data),]
center.output<-sapply(model.data[, best.covariates], 
                      function(x){
                        (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
                      })

model.data<-as.data.frame(cbind(model.data[,!colnames(model.data) %in% 
                                             colnames(center.output)],
                                center.output))

null.model<-lme(mass.per.year ~ 1, 
                random=as.formula(paste("~1|", model.random.effects, sep="")), 
                data=model.data)
summary(null.model)
AIC(null.model)
r.squaredGLMM(null.model)
plot(null.model)

capture.output(summary(null.model),
               file=paste0("./Outputs/null model coefs", Sys.Date(),".txt"))

#             AGE ONLY MODEL ####

age.model<-lme(mass.per.year ~ log.age, 
               random=as.formula(paste("~1|", model.random.effects, sep="")), 
               data=model.data)
summary(age.model)
AIC(age.model)
r.squaredGLMM(age.model)

capture.output(summary(age.model),
               file=paste0("./Outputs/age model coefs", Sys.Date(),".txt"))

#             ENVIRONMENT ONLY MODELS ####

env.model<-lme.model.selection("mass.per.year",
                               best.covariates,
                               NULL,
                               two.way.interactions(best.covariates),
                               c("IBRAsub","planting"),
                               model.data,
                               2)

summary(env.model)
AIC(env.model)
r.squaredGLMM(env.model)

capture.output(summary(env.model),
               file=paste0("./Outputs/env model coefs", Sys.Date(),".txt"))

#             SPECIES RICHNESS MODELS ####

model.data<-model.data[complete.cases(model.data),]

center.output<-sapply(model.data[, c(best.covariates, "rare.rich")], 
                      function(x){
                        (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
                      })

model.data<-as.data.frame(cbind(model.data[,!colnames(model.data) %in% 
                                             colnames(center.output)],
                                center.output))

rich.model<-lme.model.selection("mass.per.year",
                                c(best.covariates, "rare.rich"),
                                NULL,
                                two.way.interactions(c(best.covariates, 
                                                       "rare.rich")),
                                c("IBRAsub","planting"),
                                model.data,
                                2)
summary(rich.model)
AIC(rich.model)
r.squaredGLMM(rich.model)

capture.output(summary(rich.model),
               file=paste0("./Outputs/species richness coefs", Sys.Date(),".txt"))

#                             SPLIT INTO CLUSTER REGIONS ####
env.model.update<-update(rich.model, .~. -rare.rich)

# extract out population level residuals to feed into richness model
model.data$residuals<-residuals(env.model.update, level=0)

model.data<-model.data[order(model.data$best.cluster),]
model.data$best.cluster.fact<-as.factor(model.data$best.cluster)

# test how clusters differ from eachother in biomass
cluster.model.nomono<-lme(residuals ~ best.cluster.fact, 
                   random=as.formula(paste("~1|", model.random.effects, sep="")), 
                   data=model.data,
                   method="ML")
summary(cluster.model.nomono)
r.squaredGLMM(cluster.model.nomono)

boxplot(model.data$residuals ~ model.data$best.cluster.fact)

cluster.hyptest.nomono<-glht(cluster.model.nomono, 
                      linfct= mcp(best.cluster.fact = "Tukey"))
summary(cluster.hyptest.nomono)

# Add in richness to see whether richness prod relationship differs in clusters
richness.only.cluster.model.nomono<-lme(residuals ~ rare.rich * best.cluster.fact, 
                                 random=as.formula(paste("~1|", model.random.effects, sep="")), 
                                 data=model.data, method="ML")
summary(richness.only.cluster.model.nomono)
r.squaredGLMM(richness.only.cluster.model.nomono)

richness.cluster.hyptest.nomono<-glht(richness.only.cluster.model.nomono, 
                               linfct = c("rare.rich = 0",
                                          paste("rare.rich + ", paste("rare.rich:best.cluster.fact", 2:8, sep=""), "= 0")))
summary(richness.cluster.hyptest.nomono)

capture.output(summary(richness.cluster.hyptest.nomono),
               file=paste0("./Outputs/species richness cluster coefs", Sys.Date(),".txt"))

#             MASS RATIO MODELS ####
mass.ratio.vars<-colnames(mixed.site.sub)[grepl("cwm", colnames(mixed.site.sub))]

model.data<-mixed.site.sub.multi[c("mass.per.year", mass.ratio.vars, best.covariates,
                                   "IBRAsub", "planting", "best.cluster")]

model.data<-model.data[complete.cases(model.data),]

center.output<-sapply(model.data[, c(best.covariates, mass.ratio.vars)], 
                      function(x){
                        (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
                      })

model.data<-as.data.frame(cbind(model.data[,!colnames(model.data) %in% 
                                             colnames(center.output)],
                                center.output))

vif(model.data[,c(best.covariates,mass.ratio.vars)])

model.data[!complete.cases(model.data),]

mass.ratio.env<-update(env.model, as.formula(paste0(".~. + ",
                                                    paste0(mass.ratio.vars,
                                                           collapse=" + "))))

summary(mass.ratio.env)
AIC(mass.ratio.env)
trait.coefs.nomono<-summary(mass.ratio.env)$tTable[
                               order(abs(summary(mass.ratio.env)$tTable[,1]),
                                     decreasing=TRUE),]
trait.coefs<-round(trait.coefs.nomono, 3)
write.csv(trait.coefs, "./Outputs/trait env model coefs nomono.csv")

capture.output(summary(mass.ratio.env),
               file=paste0("./Outputs/mass ratio coefs", Sys.Date(),".txt"))

# so now we've got two significant interaction effects between traits and structure

#                             SPLIT INTO CLUSTER REGIONS ####

# Because I found a significant relationship between richness and biomass in a single
# region, there are two things I can investigate:

# 1. Can I identify the fuctional components of that significant relationship?
# 2. Are other trait-productivity relationships hiding in non-significant global trends?

# subset our data frame so the residuals match up (taken from model selection function)
mass.ratio.env.only<-lme.model.selection(response="mass.per.year",
                                            main.fixed=best.covariates,
                                            covariate.fixed=NULL,
                                            interactions=c("log.age:tmin.annual", 
                                                           "log.density:log.aridity.index",
                                                           "log.density:log.age"),
                                            random.effects=c("IBRAsub","planting"),
                                            dataframe=model.data,
                                            AIC.threshold=2)

model.data.sub<-mass.ratio.env.only$data
model.data.sub$residuals<-residuals(mass.ratio.env.only, level=0)
head(model.data.sub)

# sort data by cluster and make it a factor (so model assesses different slopes per cluster)
model.data.sub<-model.data.sub[order(model.data.sub$best.cluster),]
model.data.sub$best.cluster.fact<-as.factor(model.data.sub$best.cluster)

# now model richness across clusters

mass.ratio.cluster<- lme(as.formula(paste("residuals ~ ", 
                                           "(", paste(mass.ratio.vars, 
                                                      collapse=" + "), 
                                           ")",
                                           "* best.cluster.fact")),
                          random=~1|IBRAsub/planting,
                          data=model.data.sub,
                          na.action="na.fail")

model.data<-mass.ratio.cluster$data

int.terms<-expand.grid(mass.ratio.vars,
                       "best.cluster.fact", 
                       sort(unique(model.data$best.cluster.fact))[-1])
int.terms<-paste0(int.terms[,1],":",int.terms[,2],int.terms[,3])

mass.ratio.cluster.hyptest<-glht(mass.ratio.cluster, 
                                 linfct = c(paste(mass.ratio.vars, "= 0"),
                                            paste(mass.ratio.vars, "+" ,
                                                  int.terms,
                                                  "= 0")))
summary(mass.ratio.cluster.hyptest)$model$fitted

mass.ratio.cluster.coefs<-summary(mass.ratio.cluster.hyptest)
str(summary(mass.ratio.cluster.hyptest))

capture.output(summary(mass.ratio.cluster.hyptest),
               file=paste0("./Outputs/mass ratio cluster coefs", Sys.Date(),".txt"))

#             FUNCTIONAL DIVERSITY MODELS ####
fun.div.vars<-c("FRv","FEm","FDm")

with(mixed.site.sub.multi, cor(cbind(FRv, FEm, FDm, log.age, log.density,
                                      log.aridity.index, tmin.annual, sand),
                                use="complete.obs"))

model.data<-mixed.site.sub.multi[,c("mass.per.year", fun.div.vars, 
                                     best.covariates, mass.ratio.vars,
                             "rare.rich", "IBRAsub", "planting", "best.cluster")]

center.output<-sapply(model.data[, c(best.covariates, fun.div.vars)], 
                      function(x){
                        (x-mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
                      })

model.data<-as.data.frame(cbind(model.data[,!colnames(model.data) %in% 
                                             colnames(center.output)],
                                center.output))

model.data<-model.data[complete.cases(model.data),]

vif(model.data[,c(best.covariates,fun.div.vars)])

fun.div.env<-lme.model.selection(response="mass.per.year",
                                    main.fixed=fun.div.vars,
                                    covariate.fixed=best.covariates,
                                    interactions=c(paste(expand.grid(fun.div.vars, 
                                                                     best.covariates)[,1],
                                                         expand.grid(fun.div.vars, 
                                                                     best.covariates)[,2], sep=":"),
                                                   "log.age:tmin.annual", 
                                                   "log.density:log.aridity.index",
                                                   "log.density:log.age"),
                                    random.effects=c("IBRAsub","planting"),
                                    dataframe=model.data,
                                    AIC.threshold=2)
summary(fun.div.env)
fun.div.coefs<-summary(fun.div.env)$tTable[
  order(abs(summary(fun.div.env)$tTable[,1]),
        decreasing=TRUE),]
trait.coefs<-round(fun.div.coefs, 3)
write.csv(trait.coefs, "./Outputs/fun div model coefs.csv")
r.squaredGLMM(fun.div.env)

capture.output(summary(fun.div.env),
               file=paste0("./Outputs/functional diversity coefs", Sys.Date(),".txt"))

# so now we've got two significant interaction effects between traits and structure
# and only the age:tmin env interactions stays in the model. This suggests that variance
# in traits can explain the other interactions.

#                             SPLIT INTO CLUSTER REGIONS ####

# Because I found a significant relationship between richness and biomass in a single
# region, there are two things I can investigate:

# 1. Can I identify the fuctional components of that significant relationship?
# 2. Are other trait-productivity relationships hiding in non-significant global trends?

# subset our data frame so the residuals match up (taken from model selection function)
fun.div.env.only<-lme.model.selection(response="mass.per.year",
                                         main.fixed=best.covariates,
                                         covariate.fixed=NULL,
                                         interactions=c("log.age:tmin.annual", 
                                                        "log.density:log.aridity.index",
                                                        "log.density:log.age"),
                                         random.effects=c("IBRAsub","planting"),
                                         dataframe=model.data,
                                         AIC.threshold=2)

model.data.sub<-fun.div.env.only$data
model.data.sub$residuals<-residuals(fun.div.env.only, level=0)
head(model.data.sub)

# sort data by cluster and make it a factor (so model assesses different slopes per cluster)
model.data.sub<-model.data.sub[order(model.data.sub$best.cluster),]
model.data.sub$best.cluster.fact<-as.factor(model.data.sub$best.cluster)

head(model.data.sub)
# now model richness across clusters

fun.div.cluster<- lme(as.formula(paste("residuals ~ ", 
                                          "(", paste(fun.div.vars, 
                                                     collapse=" + "), 
                                          ")",
                                          "* best.cluster.fact")),
                         random=~1|IBRAsub/planting,
                         data=model.data.sub,
                         na.action="na.fail")

model.data<-fun.div.cluster$data

int.terms<-expand.grid(fun.div.vars,
                       "best.cluster.fact", 
                       sort(unique(model.data$best.cluster.fact))[-1])
int.terms<-paste0(int.terms[,1],":",int.terms[,2],int.terms[,3])

fun.div.cluster.hyptest<-glht(fun.div.cluster, 
                                 linfct = c(paste(fun.div.vars, "= 0"),
                                            paste(fun.div.vars, "+" ,
                                                  int.terms,
                                                  "= 0")))
summary(fun.div.cluster.hyptest)

capture.output(summary(fun.div.cluster.hyptest),
               file=paste0("./Outputs/fun div cluster coefs", Sys.Date(),".txt"))

#             COMPARE FINAL MODELS ####

# mono lists

summary(mass.ratio.env)
final.model.list<-lapply(list(null.model,
                              age.model,
                              env.model,
                              rich.model,
                              mass.ratio.env,
                              fun.div.env),
                         function(x){
                           update(x, .~., data=model.data,
                                          method="ML")})
names(final.model.list)<-c("null","age","env","rich","mass","fdiv")

VarCorr(mass.ratio.env)

sapply(final.model.list, AIC)
sapply(final.model.list, logLik)

null.comparison<-lapply(final.model.list[-1], function(x){
  anova(final.model.list[[1]],
        x)
})
  
env.comparison<-lapply(final.model.list[-c(1,2,3)], function(x){
  anova(final.model.list[[3]],
        x)
})

lapply(final.model.list, VarCorr)


#             WRITE COEFFICIENTS ####

model.list<-list(null.model, age.model, env.model, 
                 rich.model, mass.ratio.env, fun.div.env)
names(model.list)<-c("null", "age", "env", "rich", "mass", "fundiv")

# model coefs
lapply(1:length(model.list),
       function(x){       
         temp<-round(summary(model.list[[x]])$tTable, 3)
         write.csv(temp, paste0("./Outputs/", names(model.list[x]), "model coefs.csv"))
                   })

model.list<-list(richness.cluster.hyptest.nomono, 
                 mass.ratio.cluster.hyptest, 
                 fun.div.cluster.hyptest)
names(model.list)<-c("rich", "mass", "fundiv")
# cluster coefs
lapply(1:length(model.list),
       function(x){       

         temp.stats<-do.call("cbind", sapply(c("coefficients","sigma","tstat","pvalues"), function(y){
           summary(model.list[[x]])$test[y]
         }))
         
         temp.stats<-round(temp.stats, 3)
         
         write.csv(temp.stats, paste0("./Outputs/", names(model.list[x]), "cluster coefs .csv"))
       })

# ####
# POST-MODEL PLOTTING ####
plot.data<-mixed.plant.sub.raw
plot.data<-merge(mixed.plant.sub.raw, mixed.site.sub[,c("pl.p","best.cluster")])


#             ORDINATION OF CLUSTERS ####
model.plant<-droplevels(mixed.plant.sub[mixed.plant.sub$pl.p %in% mixed.site.sub.multi$pl.p,])
ssmat<-with(model.plant[model.plant$genus.only==0 &
                          model.plant$dead==0 &
                          model.plant$unknown==0,], 
     table(pl.p, species))

ssmat<-ssmat[rowSums(ssmat)>0,colSums(ssmat)>0]
dim(ssmat)
ssmat.bin<-ifelse(ssmat>0,1,0)

# one plot is throwing the ordination off because it's the only plot with a 
# particular species (and it contains no other species)
dim(ssmat.bin)
ssmat.bin<-ssmat.bin[!rownames(ssmat.bin) %in% c("POAP7.B", "BWRP1.A", "Suttons.3"),]

ord<-metaMDS(ssmat.bin, distance="jaccard")

cluster.df<-read.csv("./Data/best cluster summary.csv", header=TRUE)

ssmat.cluster<-mixed.site.sub$best.cluster[match(rownames(ssmat.bin),
                                                 mixed.site.sub$pl.p)]
table(ssmat.cluster)


colour.mat<-cluster.df[ssmat.cluster, c("red","green","blue")]

pdf("./Plots/plot binary ordination.pdf", height=3.5, width=7, useDingbats=FALSE)
par(mfrow=c(1,2), las=1, tck=-0.015, mgp=c(3,0.5,0), ps=8,
    mar=c(2.5,3,1,0.5))

plot(ord$points,  
     pch=16, col=rgb(colour.mat[,1], colour.mat[,2], colour.mat[,3], 0.7),
     axes=FALSE, xlab="", ylab="")
axis(side=1, mgp=c(3,0.1,0))
axis(side=2)
box()
mtext(side=1, text="MDS1", line=1)
mtext(side=2, text="MDS2", line=1.75, las=0)
text(x=relative.axis.point(0.05,"x"),
     y=relative.axis.point(0.95, "y"),
      labels="(a)", font=2, cex=1.25)


plot(ord$points, type="n", axes=FALSE, xlab="", ylab="")

sapply(1:max(ssmat.cluster), function(x){
  
ordiellipse(ord, ssmat.cluster, draw="polygon", show.groups=x,
            col=rgb(cluster.df[x,"red"], cluster.df[x,"green"], cluster.df[x,"blue"]),
            alpha= 0.9*255)
})
axis(side=1, mgp=c(3,0.1,0))
axis(side=2)
box()
mtext(side=1, text="MDS1", line=1)
mtext(side=2, text="MDS2", line=1.75, las=0)
text(x=relative.axis.point(0.05,"x"),
     y=relative.axis.point(0.95, "y"),
     labels="(b)", font=2, cex=1.25)
dev.off()



#             COR MATRIX AND UNIVARIATE EFFECT SIZE PLOT ####
mixed.site.sub.multi<-mixed.site.sub.multi[complete.cases(mixed.site.sub.multi[,covariates]),]
mixed.site.centre <- mixed.site.sub.multi

# read in complete dredge on models
m1.full.pd<-readRDS("./Outputs/Covariate dredge.rds")

scores<-colSums(ifelse(!is.na(m1.full.pd[,colnames(m1.full.pd) %in% covariates]),
               TRUE, FALSE) * m1.full.pd$weight)
directions<-ifelse(colMeans(m1.full.pd[,colnames(m1.full.pd) %in% covariates], na.rm=TRUE)>0,
                   1, -1)

col.vect<-data.frame(value=c(seq(-1, 1, 0.01)),
                     col=c(colorRampPalette(c(rgb(1,0,0),"white"), bias=2)(100),
                           "white",
                           rev(colorRampPalette(c(rgb(0,0,1),"white"), bias=2)(100))))

#                             PLOT ####
cov.mat<-cor(mixed.site.centre[,names(sort(scores))], 
             use="complete.obs")

colnames(cov.mat)

cov.mat[lower.tri(cov.mat, diag=TRUE)]<-NA

codes<-data.frame(rcodes=names(sort(scores)),
                  str.codes=c("Nitrogen","Carbon","pH","Phosphorus", 
                              "Silt", "TWI", "Prec Seas",
                              "Clay", "Soil Depth", "Elevation",
                              "Slope", "Tmax", "Solar", "Elev Rel",
                              "Tmin", "Rainfall", "Moisture", "Sand",
                              "Age", "Density"),
                  colour=c(rep("grey50",14),"black","grey50",rep("black",4)),
                  font=c(rep(1, 14), 2, 1, rep(2,4)),
                  pch=c(rep(1,14),16,1, rep(16,4)))


pdf(paste0("./Plots/cov mat and dredge", Sys.Date(), ".pdf"), 
           width=7.086, height=5, useDingbats=FALSE)

framemat=rbind(c(0.11,0.675,0.1,0.9),
               c(0.79,0.895,0.1,0.9),
               c(0.905,0.925,0.3,0.7),
               c(0.11,0.885,0.1,0.9))

split.screen(framemat)

#                             MAIN PLOT ####
screen(1)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=8, tck=-0.01, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, xlim=c(-0.025,1.2), ylim=c(0,1.15), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
#axis(side=1)
#axis(side=2)
x.start<-0.2
x.finish<-1.2
y.start<-0
y.finish<-1

x.dist<-(x.finish-x.start)/(dim(cov.mat)[1])
y.dist<-(y.finish-y.start)/(dim(cov.mat)[1])

# get rect coordinates ready to draw
test<-cbind(expand.grid(seq(x.start, x.finish, length.out=20)[-20],
                       seq(y.start, y.finish, length.out=20)[-20]),
           expand.grid(seq(x.start, x.finish, length.out=20)[-1],
                       seq(y.start, y.finish, length.out=20)[-1]))
colnames(test)<-c("xleft","ybottom","xright","ytop")
test[which(is.na(t(cov.mat[-dim(cov.mat)[1],-1]))),]=NA

# get colours of cov mat
cov.mat.vect<-as.vector(t(cov.mat[-dim(cov.mat)[1], -1]))
cov.mat.vect<-cov.mat.vect[!is.na(cov.mat.vect)]

col.scores<-col.vect$col[match(round(cov.mat.vect, 2),
                               round(col.vect$value, 2))]

test.nona<-test[complete.cases(test),]

rect(xleft=par("usr")[1],
     xright=0.2,
     ybottom=par("usr")[3],
     ytop=y.finish + y.dist,
     border=NA, col="grey85")

rect(xleft=test.nona$xleft, 
     xright=test.nona$xright, 
     ybottom=test.nona$ybottom, 
     ytop=test.nona$ytop,
     col=as.character(col.scores))

axis(side=4, at=seq(y.start + 0.5*y.dist, 
                    y.finish + 0.5*y.dist, 
                    length.out=length(rownames(cov.mat))),
     labels=NA, las=1)

mtext(side=4, at=seq(y.start + 0.5*y.dist, 
                     y.finish + 0.5*y.dist, 
                     length.out=length(rownames(cov.mat))),
      text=codes$str.codes, adj=0.5, line=2.025,
      col=as.character(codes$colour), font=codes$font)

axis(side=1, at=seq(x.start +  0.5*x.dist, 
                    x.finish - 0.5*x.dist, 
                    length.out=length(colnames(cov.mat))-1),
     labels=NA, las=3)

par(xpd=NA)
text(x=seq(x.start +  0.5*x.dist, 
           x.finish - 0.5*x.dist, 
           length.out=length(colnames(cov.mat))-1),
     y=-0.025, labels= codes$str.codes[-1], srt=40,
     adj=1,
     col=as.character(codes$colour)[-1], font=codes$font[-1])
segments(x0=par("usr")[2], x1=par("usr")[2],
         y0=par("usr")[3], y1=y.finish+y.dist)

# custom axis lines

segments(x0=par("usr")[1], 
         x1=par("usr")[1],
         y0=par("usr")[3], 
         y1=y.finish + y.dist)

segments(x0=par("usr")[1], 
         x1=1,
         y0=y.finish + y.dist, 
         y1=y.finish + y.dist)

segments(x0=par("usr")[1], 
         x1=par("usr")[2],
         y0=par("usr")[3], 
         y1=par("usr")[3])

segments(x0=par("usr")[2], 
         x1=1.5,
         y0=y.finish + y.dist, 
         y1=y.finish + y.dist)

segments(x0=par("usr")[1], 
         x1=1.5,
         y0=par("usr")[3], 
         y1=par("usr")[3])

par(xpd=FALSE)

axis(side=3, at=c(0,0.2,0.4,0.6,0.8,1), pos=y.finish + y.dist,
     mgp=c(3,0.2,0))
axis(side=3, at=seq(0,1,0.1), pos=y.finish + y.dist,
     labels=NA, tck=-0.005)

mtext(side=3, text="Summed Akaike weights of models containing variable",
      line=-0.75)

axis(side=2, at=seq(y.start + 0.5*y.dist,
                    y.finish + 0.5*y.dist, length.out=20),
     labels=NA)

mtext(side=2, at=seq(y.start + 0.5*y.dist, 
                     y.finish + 0.5*y.dist, 
                     length.out=length(rownames(cov.mat))),
      text=codes$str.codes, adj=1, line=0.5,
      col=as.character(codes$colour), font=codes$font)

# ADD ON POINTS

sapply(1:length(scores), function(x){
  
  pos<-seq(y.start + 0.5*y.dist,
      y.finish + 0.5*y.dist, length.out=20)[x]
  
  rect(ybottom=pos-0.25*y.dist,
       ytop=pos+0.25*y.dist,
       xleft=par("usr")[1],
       xright=sort(scores)[x],
       col=ifelse(codes$pch[x]==16, "black", "white"))

})

#points(x=sort(scores), 
#       y= seq(y.start + 0.5*y.dist,
#              y.finish + 0.5*y.dist, length.out=20),
#       pch=codes$pch)

close.screen(1)

#                             DIVERSITY CORRELATIONS ####
screen(2)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=8, tck=-0.08, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, xlim=c(-0.025,1.2), ylim=c(0,1.15), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")

axis(side=2, at=seq(y.start + 0.5*y.dist, 
                    y.finish + 0.5*y.dist, 
                    length.out=length(rownames(cov.mat))),
     labels=NA, las=1)

# get rect coordinates ready to draw
test<-cbind(expand.grid(seq(par("usr")[1], par("usr")[2], length.out=5)[1:4],
                        seq(y.start, y.finish+y.dist, length.out=21)[1:20]),
            expand.grid(seq(par("usr")[1], par("usr")[2], length.out=5)[2:5],
                        seq(y.start, y.finish+y.dist, length.out=21)[2:21]))
colnames(test)<-c("xleft","ybottom","xright","ytop")
#test[which(is.na(t(cov.mat[-dim(cov.mat)[1],-1]))),]=NA

cov.mat<-cor(mixed.site.centre[,c("rare.rich","FRv", "FEm", "FDm",
                                           names(sort(scores)))])
cov.mat<-cov.mat[-(1:4),1:4]

cov.colours<-col.vect$col[match(round(t(cov.mat),2),
                                round(col.vect$value, 2))]

rect(xleft=test$xleft, 
     xright=test$xright, 
     ybottom=test$ybottom, 
     ytop=test$ytop,
     col=as.character(cov.colours))

axis(side=1, at=seq(par("usr")[1], par("usr")[2], length.out=9)[c(2,4,6,8)],
labels=NA)

par(xpd=NA)
text(x=seq(par("usr")[1], par("usr")[2], length.out=9)[c(2,4,6,8)],
     y=-0.025, labels= c("Richness",
                         "FRv", "FEm", "FDm"), srt=40,
     adj=1)
segments(x0=par("usr")[1], x1=par("usr")[1],
         y0=par("usr")[3], y1=y.finish+y.dist)
par(xpd=FALSE)

close.screen(2)

#                             LEGEND ####
screen(3)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=8, tck=-0.3, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(-1,1), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")

legend.list<-list(x=seq(0,1,length=1),
                  y=seq(-1,1,length=201),
                  z=matrix(rep(seq(-1,1,length=201), 1), nrow=1, ncol=201))

image(legend.list,
      col=as.character(col.vect$col), 
      axes=F, las=1, useRaster=TRUE, add=TRUE)

axis(side=4, labels=NA)
mtext(side=4, at=c(-1,-0.5,0,0.5,1), text=c(-1,-0.5,0,0.5,1),
      line=0.35)
mtext(side=4, line=1.45, las=0, text="Correlation coefficient")
box()
close.screen(3)

close.screen(all.screens=TRUE)
dev.off()

screen(4)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=8, tck=-0.3, mgp=c(3,0.5,0))
plot(x=NULL,y=NULL, xlim=c(0,1), ylim=c(0,1.15), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
par(xpd=NA)
rect(xleft=par("usr")[1],
     xright=par("usr")[2],
     ybottom=par("usr")[3],
     ytop=y.finish+y.dist, lwd=1.25)
par(xpd=FALSE)
close.screen(4)



#             NEW COMBINED RICHNESS PLOT ####
plant<-read.csv("./Data/plantscale.csv", header=TRUE)

ausmap<-readShapeSpatial("/home/timothy/University files - offline/PhD - offline/Shape files/Australia outline/aust_cd66states.shp")
proj4string(ausmap)<-"+proj=longlat +ellps=aust_SA +no_defs"
prec<-raster("/home/timothy/University files - offline/PhD - offline/Shape files/prec.annual.tif")
log.prec<-log(prec)
proj4string(log.prec)<-"+proj=longlat +ellps=aust_SA +no_defs"

mixed.plots<-mixed.site.sub[mixed.site.sub$alpha>1,]
mixed.planting.sub<-mixed.plots[!duplicated(mixed.plots$planting),]
cluster.mean.coords<-t(sapply(split(mixed.planting.sub, 
                                    f=mixed.planting.sub$best.cluster), 
                              function(x){c(mean(x$lat, na.rm=TRUE), 
                                            mean(x$long, na.rm=TRUE))}))

plot(ausmap)
with(mixed.planting.sub[mixed.planting.sub$best.cluster==8,], points(lat ~ long))

# mean environmental conditions at each cluster
mean.cluster.env<-as.data.frame(sapply(mixed.site.sub[,raw.predictors],
                        function(x){
                        tapply(x, mixed.site.sub$best.cluster, function(y){mean(y, na.rm=TRUE)}) 
                        }))

# 5 most common species in each cluster
plant.sub<-merge(plant,
                 mixed.site.sub[,colnames(mixed.site.sub) %in% 
                                 c("pl.p", "best.cluster")])

cluster.common.sp<-lapply(split(plant.sub, plant.sub$best.cluster),
                         function(x){
                            cbind(names(sort(table(droplevels(x$species)), 
                                       decreasing=TRUE))[1:10],
                                  sort(table(droplevels(x$species)), 
                                       decreasing=TRUE)[1:10])
                          })

# read in details of clusters
cluster.df<-read.csv("./Data/best cluster summary.csv", header=TRUE)
cluster.df<-as.data.frame(sapply(cluster.df, function(x){gsub("\\+AC0", "", x)}))
cluster.df[,c("red","green","blue")]<-sapply(cluster.df[,c("red","green","blue")],
                                             function(x){as.numeric(as.character(x))})

plot(ausmap, lwd=0.5, xlim=c(115,155), ylim=c(-45,-9), axes=FALSE, border="grey60")
map.usr<-par("usr")
# 13 long, 10 tall

#framemat<-rbind(c(0.075,0.99,0.35,0.95), # map
 #               c(0.075, 0.375,99, 0.05, 0.4), # coefficients
 # c(0.41, 0.54, 0.19, 0.29), # cluster 4 - Temperate forest
# c(0.28, 0.41, 0.32, 0.42), # cluster 2 - high altitude woodland
# c(0.25, 0.38, 0.12, 0.22), # cluster 1 - eastern mallee
# c(0.095, 0.225, 0.22, 0.32), # cluster 6 - southern mallee
# c(0.80, 0.93, 0.19, 0.29), # cluster 3 - Coastal plains

framemat<-rbind(c(0.075,0.99, 0.025, 0.625), # raster
                c(0.075,0.99, 0.025, 0.625), # map
                
                c(0.075,0.99, 0.62, 0.99), # box for top plots
                
                c(0.2, 0.475, 0.64, 0.97), # coefficients
                
                c(0.57, 0.725, 0.28, 0.39), # cluster 1 - Dry woodland
                c(0.53, 0.685, 0.0475, 0.1575), # cluster 2 - Wet woodland
                c(0.115, 0.265, 0.30, 0.41), # cluster 3 - Coastal plains
                c(0.80, 0.955, 0.05, 0.16), # cluster 4 - Temperate forest
                c(0.82, 0.975, 0.20, 0.31), # cluster 5 - Sub-tropical forest
                c(0.33, 0.485, 0.12, 0.23), # cluster 6 - Mediterranean woodland
                c(0.32, 0.475, 0.28, 0.39), # cluster 7 - Semi-arid woodland
                c(0.81,0.975,0.44,0.55),    # cluster 8 - Rainforest
                c(0.58, 0.71, 0.28, 0.38),  # cluster 9 
                c(0.76, 0.89, 0.07, 0.17),  # cluster 7
                
                c(0.56, 0.71, 0.84, 0.97), # Interaction plot (b)
                c(0.82, 0.97, 0.84, 0.97), # Interaction plot (c)
                c(0.56, 0.71, 0.64, 0.77), # Interaction plot (d)
                c(0.82, 0.845, 0.64, 0.77), # Interaction gradient plot (midpoints 0.82)
                
                c(0.14,0.315, 0.095,0.14), # rainfall density
                c(0.14,0.315,0.08,0.095), # rainfall legend
                c(0.075,0.99,0.025,0.59)) # final box

pdf(paste0("./Plots/Combined plot ", Sys.Date(), " - massperyear.pdf"), width=9, height=10, 
    useDingbats=FALSE)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=10, las=1, tck=-0.0075, mgp=c(3,0.4,0))
split.screen(framemat)

#                             MAP ####
screen(1)
par(ps=10)
plot(x=NULL,y=NULL, xlim=c(107.8123,162.1877), ylim=c(-46.44,-7.56), axes=FALSE)
plot(log.prec, col=colorRampPalette(c("white",rgb(0.9,0.9,1,1), rgb(0.15,0.15,1,1)))(200), 
     box=FALSE, add=TRUE, legend=FALSE, maxpixels=1000000, interpolate=TRUE)

mappar<-par("usr")
close.screen(1)

screen(2)
par(ps=10)
plot(x=NULL,y=NULL, xlim=c(107.8123,162.1877), ylim=c(-46.44,-7.56), axes=FALSE)
plot(ausmap, lwd=1, axes=FALSE, border="grey40", add=TRUE)

text(x=relative.axis.point(0.025, "x"),
     y=-9.5,
     labels="(e)", font=2, cex=1.25)
axis(side=1, at=seq(110,150,10), labels=parse(text=paste(seq(110,150,10), "*degree~E", sep="")))
axis(side=2, at=seq(-45,-10,5), labels=parse(text=paste(seq(-45,-10,5), "*degree~S", sep="")),
     mgp=c(3,0.7,0))

# plot the points for each cluster
lapply(cluster.df$cluster, function(x){
  cluster.sub<-mixed.planting.sub[mixed.planting.sub$best.cluster==x,]
  
  with(cluster.sub, points(lat ~ long, pch=21, lwd=0.5, 
                           bg=with(cluster.df, rgb(red[x], green[x], blue[x], 1)))) 
})

close.screen(2)

screen(3)
plot.new()
#rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
#     border=NA, col="white")
#box(lwd=1.5)
close.screen(3)

#                             COFFICIENT PLOT ####
screen(4)
rich.coefs<-coefTable(rich.model)[-1,-3]
position<-7-c(4,3,2,5.5,5,1.5,3.5,4.5,2.5)
#position<-7-c(3,1,2,5,4,6,2.5,3.5,1.5)

labels=c("Age (years)", expression("Density (plants ha"^"-1"*")"),
         "Moisture availability", expression("Min. temp ("*degree*"C)"),
         "Sand (%)", "Species richness",
         "Age:Density", "Age:Min. temp", "Density:Moist. avail.")
  
par(mar=c(0,0,0,0), las=1, tck=-0.02, mgp=c(3,0.2,0), ps=10)
plot(position ~ rich.coefs[,1], xlim=c(-0.25,0.55), 
     ylim=c(1.25,5.75), type="n", axes=FALSE,
     xlab="", ylab="")
rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
     border=NA, col="white")
points(position ~ rich.coefs[,1], pch=16,
     col=c(rep("black", 6), rep("grey50", 3)))
text(position[7:9]+0.2 ~ rich.coefs[7:9,1], labels=c("(c)","(d)","(b)"), col="grey50")
box()
segments(x0=0, x1=0, y0=0, y1=8, lty="dashed")
segments(y0=position, y1=position, lwd=1.5,
         x0=rich.coefs[,1]-1.98*rich.coefs[,2],
         x1=rich.coefs[,1]+1.98*rich.coefs[,2],
         col=c(rep("black", 6), rep("grey50", 3)))

axis(side=1, mgp=c(3,0.2,0))
axis(side=2, at=position[1:6], labels=labels[1:6], las=1, mgp=c(3,0.5,0))
axis(side=2, at=position[7:9], labels=labels[7:9], las=1, col.axis="grey50", 
     col.ticks="grey50", mgp=c(3,0.5,0))
mtext(side=1, text="Standardized slope estimate", line=1.25, cex=1.1)

par(xpd=NA)
text(x=relative.axis.point(-0.05,"x"), y=relative.axis.point(1.03, "y"),
    labels="(a)", font=2)
close.screen(4)

#                             RICHNESS PLOTS ####
for (i in 1:8){
  
  ### PREDICTIONS
  cluster.model.data<-richness.only.cluster.model.nomono$data
  
  cluster.data.sub<-cluster.model.data[cluster.model.data$best.cluster.fact==i,]
  cluster.random<-data.frame(planting=cluster.model.data$planting, 
                             IBRAsub=cluster.model.data$IBRAsub)
  
  cluster.random<-cluster.random[!duplicated(cluster.random$planting),]
  head(cluster.random)
  model.coefs<-c("rare.rich", "best.cluster.fact")
  
  prediction.df<-as.data.frame(matrix(0, nrow=100, ncol=length(model.coefs)))
  colnames(prediction.df)<-model.coefs
  prediction.df[,1]=seq(min(cluster.data.sub[, "rare.rich"]),
                        max(cluster.data.sub[, "rare.rich"]),
                        length.out=100)
  prediction.df[,2]=as.factor(i)
  
  predictions.temp<-predict(object=richness.only.cluster.model.nomono, 
                            newdata=prediction.df, level=0, se.fit=TRUE)
  predictions<-list(fit=predictions.temp$fit,
                    x=prediction.df[,1],
                    se=predictions.temp$se)
  
  screen(4+i)
  par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=10, las=1, tck=-0.025, mgp=c(3,0.5,0))
  
  # set up x-axis limits based on all richness data
  with(cluster.model.data, plot(residuals ~ rare.rich, las=1, type="n",
                              xlim=summary(cluster.model.data$rare.rich)[c(1,6)],
                              ylim=summary(cluster.model.data$residuals)[c(1,6)]+c(0,0.6), 
                              axes=FALSE))
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, col=rgb(1,1,1,0.85))
  
  axis(side=2, at=c(-2:2), labels=c(-2:2))
  axis(side=1, at=(c(1,5,10,15)-mean(mixed.site.sub$rare.rich))/sd(mixed.site.sub$rare.rich), 
       labels=c(1,5,10,15), mgp=c(3,0.1,0))
  
  with(cluster.data.sub, points(residuals ~ rare.rich, pch=16,
                                col=with(cluster.df, rgb(red[i], green[i], blue[i], 1))))
  #rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
   #    border=NA, col=rgb(1,1,1,0.5))
  
  polygon(y=c(predictions$fit+1.96*predictions$se,
              rev(predictions$fit-1.96*predictions$se),
              predictions$fit[1]+1.96*predictions$se[1]),
          x=c(predictions$x, rev(predictions$x), predictions$x[1]), 
          col=rgb(0.2,0.2,0.2,0.7), border=NA)
  points(predictions$fit ~ predictions$x, type="l")
  #segments(x0=-5, x1=25, y0=0, y1=0, lty="dashed")
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.91, "y"),
       labels=cluster.df$name[i], font=2, adj=0)
  
  if(i==8){
    mtext(side=1, text="Species richness", line=0.9)
    mtext(side=2, text="Residual productivity", las=3, line=2.35)
        mtext(side=2, text=expression("ln(t ha"^-1*" year"^-1*")"), 
          las=3, line=1.35)
    }
  
  box()
  
  close.screen(4+i)
}

#                             INTERACTION PLOTS ####
screen(15)
par(mar=c(0,0,0,0), ps=10, las=1, tck=-0.02, mgp=c(3,0.4,0))
temp.inter.list<-interaction.matrix.fun(rich.model,
                                        "log.aridity.index", "log.density",
                                        grid.size=150, quantile=0.01)

image(temp.inter.list, zlim=c(-0.6937553, 2.7072447),
      col=colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150), 
      axes=F, las=1, useRaster=TRUE)

contour(x=temp.inter.list,
        add=T, 
        levels=log(c(0.1,0.25,0.5,1,2.5,5,10,15)),
        labels=c(0.1,0.25,0.5,1,2.5,5,10,15),
        method="edge", lwd=0.5, labcex=0.9)

axis(side=1, mgp=c(3,0,0),
     at=(log(c(0.25,0.5,1,1.5))-mean(mixed.site.sub$log.aridity.index))/sd(mixed.site.sub$log.aridity.index),
     labels=c(0.25,0.5,1,1.5))
mtext(side=1, text="Moisture availability", line=0.85, las=0, cex=1.1)

axis(side=2, at=(log(c(250,500,1000,2500,5000))-mean(mixed.site.sub$log.density, na.rm=TRUE))/sd(mixed.site.sub$log.density, na.rm=TRUE),
     labels=c(250,500,1000,2500,5000))
mtext(side=2, text=expression("Density (plants ha"^"-1"*")"), line=2.15,
      las=0, cex=1.1)

par(xpd=NA)
text(x=relative.axis.point(-0.075,"x"), y=relative.axis.point(1.065, "y"),
     labels="(b)", font=2)
box()
close.screen(15)

screen(16)
par(mar=c(0,0,0,0), ps=10, las=1, tck=-0.02, mgp=c(3,0.4,0))
temp.inter.list<-interaction.matrix.fun(rich.model, 
                                        "log.age", "log.density", 
                                        grid.size=150, quantile=0.01)

image(temp.inter.list, zlim=c(-0.6937553, 2.7072447),
      col=colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150),
      axes=F, las=1, useRaster=TRUE)

contour(x=temp.inter.list,
        add=T, 
        levels=log(c(0.1,0.25,0.5,1,2.5,5,10,15)),
        labels=c(0.1,0.25,0.5,1,2.5,5,10,15),
        method="edge", lwd=0.5, labcex=0.9)

axis(side=1, mgp=c(3,0,0), at=(log(c(5,10,15,25))-mean(mixed.site.sub$log.age))/sd(mixed.site.sub$log.age),
     labels=c(5,10,15,25))
mtext(side=1, text="Age (years)", line=0.85, cex=1.1)
axis(side=2, at=(log(c(250,500,1000,2500,5000))-mean(mixed.site.sub$log.density, na.rm=TRUE))/sd(mixed.site.sub$log.density, na.rm=TRUE),
     labels=c(250,500,1000,2500,5000))
mtext(side=2, text=expression("Density (plants ha"^"-1"*")"), line=2.15,
      las=0, cex=1.1)

par(xpd=NA)
text(x=relative.axis.point(-0.075,"x"), y=relative.axis.point(1.065, "y"),
     labels="(c)", font=2)
box()

close.screen(16)


screen(17)
par(mar=c(0,0,0,0), ps=10, las=1, tck=-0.02, mgp=c(3,0.4,0))
temp.inter.list<-interaction.matrix.fun(rich.model, 
                                        "log.age", "tmin.annual", 
                                        grid.size=150, quantile=0.01)

image(temp.inter.list, zlim=c(-0.6937553, 2.7072447),
      col=colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150), 
      axes=F, las=1, useRaster=TRUE)

contour(x=temp.inter.list,
        add=T, 
        levels=log(c(0.1,0.25,0.5,1,2.5,5,10,15)),
        labels=c(0.1,0.25,0.5,1,2.5,5,10,15),
        method="edge", lwd=0.5, labcex=0.9)

axis(side=1, mgp=c(3,0,0), at=(log(c(5,10,15,25))-mean(mixed.site.sub$log.age))/sd(mixed.site.sub$log.age),
     labels=c(5,10,15,25))
mtext(side=1, text="Age (years)", line=0.85, cex=1.1)
axis(side=2, las=1, at=(c(5,7.5,10,12.5,15,17.5)-mean(mixed.site.sub$tmin.annual))/
                       sd(mixed.site.sub$tmin.annual),
                       labels=c(5,7.5,10,12.5,15,17.5))
mtext(side=2, text=expression("Min temp ("*degree*"C)"), line=2, las=0, cex=1.1)

par(xpd=NA)
text(x=relative.axis.point(-0.075,"x"), y=relative.axis.point(1.065, "y"),
     labels="(d)", font=2)
box()

close.screen(17)

#                             EXTRAS ####

screen(18)
library(plotrix)
par(mar=c(0,0,0,0), ps=10, las=1, tck=-0.25, mgp=c(3,0.4,0))
#plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(8.4,12.3), xaxs="i", yaxs="i", axes=FALSE)

colours<-colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150)
rect.coords<-seq(6.214,9.615, length.out=150)

image(z=t(as.matrix(rect.coords)), col=colours, axes=FALSE, useRaster=TRUE)

# lapply(2:length(rect.coords), function(x){
#   rect(xleft=0, xright=1, ytop=rect.coords[x], ybottom=rect.coords[x-1],
#        border=NA, col=colours[x])
# })
axis.ticks<-log(c(0.5,1,2.5,5,10)*1000)

axis.coords<-(axis.ticks-min(rect.coords))/(max(rect.coords)-min(rect.coords))
axis(side=4, at=axis.coords, 
     labels=c(0.5,1,2.5,5,10))
mtext(side=4, text="Plot productivity", las=3, line=1.8, cex=1.1)
mtext(side=4, text=expression("(t ha"^-1*" year"^-1*")"), las=3, 
      line=3, cex=1.1)
box()
close.screen(18)

screen(19)
par(mar=c(0,0,0,0), ps=10, las=1, tck=-0.02, mgp=c(3,0.4,0), xpd=TRUE)

limits<-c(cellStats(log.prec, stat="min"), cellStats(log.prec, stat="max"))
xlims<-par("usr")[1:2]
plot(x=NULL, y=NULL, ylim=c(0,1), xlim=limits, xaxs="i", yaxs="i", axes=FALSE)

polygon(density(log(mixed.site.sub$prec.annual)), col="grey90")
axis(side=2, at=mapply(proportion=c(0,1), axis=c("y","y"), relative.axis.point),
     tck=0, labels=NA)
axis(side=2, at=c(0,0.5,1), tck=-0.065, mgp=c(3,0.4,0))
mtext(side=2, line=1.65, text="Density", las=0, cex=1.1)
close.screen(19)

screen(20)

rain.cols<-colorRampPalette(c("white",rgb(0.9,0.9,1,1), rgb(0.15,0.15,1,1)))(200)

rect.coords<-seq(limits[1], limits[2], length.out=150)

image(z=as.matrix(rect.coords), col=rain.cols, axes=FALSE, useRaster=TRUE)
box()

axis.ticks<-log(c(250,500,1000,2000,4000))

axis.coords<-(axis.ticks-min(rect.coords))/(max(rect.coords)-min(rect.coords))
axis(side=1, at=axis.coords, 
     labels=NA, tck=-0.2)
par(xpd=NA)
text(x=axis.coords, y=relative.axis.point(-0.6, "y"),
     labels=c(250,500,1000,2000,4000), srt=30, adj=1)
par(xpd=TRUE)
mtext(side=1, line=1.5, text="Annual rainfall (mm)", cex=1.1)

close.screen(20)

screen(21)
plot.new()
box(lwd=1.5)
close.screen(21)
close.screen(all.screens=TRUE)
dev.off()

#             TRAIT MODEL COEFFICIENT PLOT ####
mass.df<-data.frame(estimate=mass.ratio.cluster.coefs$test$coefficients,
                    se=mass.ratio.cluster.coefs$test$sigma,
                    t=mass.ratio.cluster.coefs$test$tstat,
                    p=mass.ratio.cluster.coefs$test$pvalues,
                    lci=confint(mass.ratio.cluster.coefs)$confint[,2],
                    uci=confint(mass.ratio.cluster.coefs)$confint[,3])

mass.df$group<-1
mass.df$group <- as.numeric(substr(rownames(mass.df),
                                   nchar(rownames(mass.df)),
                                   nchar(rownames(mass.df))))
mass.df$group[is.na(mass.df$group)] = 1
mass.df$var<-as.factor(substr(rownames(mass.df),
                      1, regexpr("\\.", rownames(mass.df))-1))

div.df<-data.frame(estimate=summary(fun.div.cluster.hyptest)$test$coefficients,
                    se=summary(fun.div.cluster.hyptest)$test$sigma,
                    t=summary(fun.div.cluster.hyptest)$test$tstat,
                    p=summary(fun.div.cluster.hyptest)$test$pvalues,
                   lci=confint(fun.div.cluster.hyptest)$confint[,2],
                   uci=confint(fun.div.cluster.hyptest)$confint[,3])

div.df$group <- as.numeric(substr(rownames(div.df),
                                   nchar(rownames(div.df)),
                                   nchar(rownames(div.df))))
div.df$group[is.na(div.df$group)] = 1
div.df$var<-as.factor(substr(rownames(div.df),
                                1, 
                             ifelse(!grepl(" ", rownames(div.df)),
                                    nchar(rownames(div.df)),       
                                    regexpr(" ", rownames(div.df))-1)))


cluster.df<-read.csv("./Data/best cluster summary.csv", header=TRUE)

pdf("./Plots/ trait and fdiv coefficient plot.pdf", height=5, width=3.5)
par(mfrow=c(4,2), mar=c(0,0,0,0), oma=c(2.5,4,0.5,0.5), mgp=c(3,0.5,0),
    tck=-0.025, ps=8)

lapply(unique(mass.df$group), function(x){
  
temp.df<-rbind(mass.df[mass.df$group==x,],
               div.df[div.df$group==x,])
temp.df<-temp.df[order(as.character(temp.df$var)),]
temp.df$order<-c(3:1, 7.5:4.5)

temp.cols<-cluster.df[cluster.df$cluster == x,]
  
plot(x=NULL, y=NULL, xlim=c(-1.25, 1.5), ylim=c(0.75,8.5), axes=FALSE, ylab="", xlab="")

if(x %in% c(1, 3, 5, 7)){
  axis(side=2, at=temp.df$order, 
        labels=c("FDm",
                 "FEm",
                 "FRv",
                 "Max. height",
                 "SLA",
                 "Seed mass",
                 "Wood density"), las=1)
  
} else {axis(side=2, at=temp.df$order,
             labels=NA, las=1)
}

abline(v=0, col="grey60")
abline(h=3.75, lty="dashed")

segments(x0=temp.df$estimate + 1.96*temp.df$se,
         x1=temp.df$estimate - 1.96*temp.df$se,
         y0=temp.df$order, y1=temp.df$order,
         col=rgb(temp.cols$red, temp.cols$green, temp.cols$blue, 1))

segments(x0=temp.df$lci,
         x1=temp.df$uci,
         y0=temp.df$order, y1=temp.df$order,
         col=rgb(temp.cols$red, temp.cols$green, temp.cols$blue, 1))


points(x=temp.df$estimate, y=temp.df$order,
       pch=21, bg="white",
       col=rgb(temp.cols$red, temp.cols$green, temp.cols$blue, 1))

rect(xleft=par("usr")[1], xright=par("usr")[2],
     ybottom=relative.axis.point(0.875, "y"), ytop=par("usr")[4],
     col="white", border=NA)

text(x=relative.axis.point(0.02, "x"),
     y=relative.axis.point(0.925, "y"), 
     labels=paste0("(", letters[x], ") ", temp.cols$full.name), 
     col=rgb(temp.cols$red, temp.cols$green, temp.cols$blue, 1),
     adj=0, font=2)

if(x %in% 7:8){
  axis(side=1, at=c(-1,-0.5,0,0.5, 1), 
       labels=c(-1,-0.5,0,0.5, 1), mgp=c(3,0,0))
  
  if(x==7){mtext(side=1, line=1.1, at=par("usr")[2],
                 text="Standardized slope estimate", cex=0.8)}
  
}

if(x %in% 3:8){
  axis(side=3, at=c(-1,-0.5,0,0.5, 1), 
       labels=NA, tck=0.025, mgp=c(3,0,0))
}

box()

sig.df<-temp.df[temp.df$p <=0.05, ]

if(dim(sig.df)[1]==0){return(NULL)}

points(x=sig.df$estimate, y=sig.df$order,
       pch=16,
       col=rgb(temp.cols$red, temp.cols$green, temp.cols$blue, 1))
})
dev.off()

#             VARIANCE PARTITIONING OF RICHNESS ####
mixed.site.sub.multi <- droplevels(mixed.site.sub.multi)

table(mixed.site.sub.multi$collecting.agency)

rich.var.m <- lmer(rare.rich ~ 1 + (1|IBRAsub/planting) + (1|collecting.agency) + (1|year.planted),
                   data=mixed.site.sub.multi)

rich.var.m1 <- lmer(rare.rich ~ 1 + (1|IBRAsub/planting) + (1|year.planted),
                   data=mixed.site.sub.multi)

rich.var.m1 <- lmer(rare.rich ~ 1 + (1|best.cluster),
                    data=mixed.site.sub.multi)

summary(rich.var.m1)
summary(rich.var.m)

ranef(rich.var.m)

dev.off()
plot(mixed.site.sub.multi$lat ~ mixed.site.sub.multi$long,
     col=ifelse(mixed.site.sub.multi$collecting.agency=="Lachlan",
                "red",
                "black"))


# ####
summary(mixed.site.sub$alpha)
# SUMMARY STATS FOR PAPER ####

plant.sub<-plant[plant$pl.p %in% mixed.site.sub.multi$pl.p,]

plant.unique.species<-plant.sub[!plant.sub$genus.only & !plant.sub$unknown &
                                  !plant.sub$dead,]

sum(plant.sub$unknown) / dim(plant.sub)[1]


sum(mixed.plant.sub.raw$unknown)

tapply(mixed.site.sub$prec.annual, mixed.site.sub$best.cluster, mean)

summary(mixed.site.sub$prec.annual)
hist(mixed.site.sub$prec.annual)

#             CORRELATION MATRIX ####

# correlation matrix for our chosen variables
colnames(mixed.site.sub)
with(mixed.site.sub, cor(cbind(rare.rich, FRv, FEm, FDm,
                               log.age, log.density,
                               log.aridity.index, tmin.annual,
                               sand), use="complete.obs"))

#             OTHER STATS ####
summary(mixed.site.sub$year.planted)
mean(mixed.site.sub.multi$age)

hist(mixed.site.sub.multi$plot.area[mixed.site.sub.multi$plot.area<1])
mean(mixed.site.sub.multi$plot.area[mixed.site.sub.multi$plot.area<1], na.rm=TRUE)

mixed.site.sub.multi

dead.prop<-tapply(mixed.plant.sub$dead, mixed.plant.sub$pl.p, function(x){
  sum(x)/length(x)
})

dead.mass<-sapply(split(mixed.plant.sub, f=mixed.plant.sub$pl.p), function(x){
  mean(x$tot.ag.mass[x$dead==1])
})
hist(dead.prop)

sum(dead.prop<0.1) / length(dead.prop)

dead.prop<-tapply(mixed.plant.sub$dead, mixed.plant.sub$pl.p, function(x){
  sum(x)/length(x)
})




max(pl.p.dist[upper.tri(pl.p.dist)])
hist(pl.p.dist[upper.tri(pl.p.dist)])

max(mixed.site.sub$prec.annual)/min(mixed.site.sub$prec.annual)
max(mixed.site.sub$biomass.area, na.rm=TRUE)/min(mixed.site.sub$biomass.area, na.rm=TRUE)


mixed.temp<-mixed.site.sub[mixed.site.sub$alpha>1,]
mixed.plant.temp<-mixed.plant.sub[mixed.plant.sub$pl.p %in% mixed.temp$pl.p,]
length(unique(mixed.plant.temp$species))

length(unique(mixed.temp$planting))


est.traits<-mixed.plant.sub[rowSums(is.na(mixed.plant.sub[,c("SLA.dist","SM.dist",
                                                             "WD.dist","MH.dist")]))<4,]

sla.est<-est.traits[!is.na(est.traits$SLA.dist),]
summary(sla.est$SLA.dist)
length(sla.est$SLA.dist[sla.est$SLA.dist>50])/
  length(sla.est$SLA.dist)

summary(sla.est$SLA.n)
sd(sla.est$SLA.n)/sqrt(length(sla.est$SLA.n))

wd.est<-est.traits[!is.na(est.traits$WD.dist),]
summary(wd.est$WD.dist)
length(wd.est$WD.dist[wd.est$WD.dist>50])/
  length(wd.est$WD.dist)

summary(wd.est$WD.n)
sd(wd.est$WD.n)/sqrt(length(wd.est$WD.n))

sm.est<-est.traits[!is.na(est.traits$SM.dist),]
summary(sm.est$SM.dist)
length(sm.est$SM.dist[sm.est$SM.dist>50])/
  length(sm.est$SM.dist)

summary(sm.est$SM.n)
sd(sm.est$SM.n)/sqrt(length(sm.est$SM.n))

MH.est<-est.traits[!is.na(est.traits$MH.dist),]
summary(MH.est$MH.dist)
length(MH.est$MH.dist[MH.est$MH.dist>50])/
  length(MH.est$MH.dist)

summary(MH.est$MH.n)
sd(MH.est$MH.n)/sqrt(length(MH.est$MH.n))

# within 50km
# SLA: 80% within 50km on mean of 10 +- 0.1 species
# WD: 81.5% within 50km on mean of 9 +- 0.075 species
# SM: 80% within 50km on mean of 10 +- 0.1 species
# MH: 75% within 50km on mean of 10 +- 0.1 species

plot(density(sla.est$SLA.dist))

summary(mixed.site.sub$plot.area*10000)
sd(mixed.site.sub$plot.area*10000, na.rm=TRUE)/
  sqrt(length(mixed.site.sub$plot.area))

summary(mixed.site.sub$density/10000)


colnames(mixed.site.sub)
cor(mixed.site.sub[,c("rare.rich", "FRv", "FEm", "FDm")],
    use="complete.obs")

age.se<-sd(mixed.site.sub$log.age)/sqrt(length(mixed.site.sub$log.age))
exp(summary(mixed.site.sub$log.age))
exp(age.se)

sum(mixed.site.sub$plant.count)

exp(summary(mixed.site.sub$log.density))

table(mixed.site.sub$data.source)/length(mixed.site.sub$data.source)

length(unique(mixed.site.sub$planting))
summary(mixed.site.sub$year.planted)

year.table<-table(mixed.site.sub$year.planted)
sum(year.table[as.numeric(names(year.table))>=1997]) / sum(year.table)

hist(mixed.site.sub$measure.date)
table(mixed.site.sub$measure.date)
dev.off()

summary(mixed.site.sub$age)
hist(mixed.site.sub$age)
age.table<-table(mixed.site.sub$age)
sum(age.table[as.numeric(names(age.table))>=7 &
                as.numeric(names(age.table))<=15]) / sum(age.table)

summary(mixed.site.sub$plant.count)

with(mixed.plant.sub, table(genus.only, dead, unknown)) / length(mixed.plant.sub$species)

length(unique(mixed.plant.sub.not.genus$species))

summary(mixed.site.sub$plot.area)
hist(log(mixed.site.sub$plot.area))

# extract trait values just used in planting subset
height.planting<-height.raw[bin.match(height.raw$species, mixed.plant.sub$species),]
sla.planting<-sla.raw[bin.match(sla.raw$species, mixed.plant.sub$species),]
seed.mass.planting<-seed.mass.raw[bin.match(seed.mass.raw$species, mixed.plant.sub$species),]
wood.density.planting<-wood.density.raw[bin.match(wood.density.raw$species, mixed.plant.sub$species),]

length(mixed.plant.sub$species[grepl("Eucalyptus", mixed.plant.sub$species)])/length(mixed.plant.sub$species)
length(unique(mixed.plant.sub$species[grepl("Eucalyptus", mixed.plant.sub$species)]))

length(mixed.plant.sub$species[grepl("Acacia", mixed.plant.sub$species)])/length(mixed.plant.sub$species)
length(unique(mixed.plant.sub$species[grepl("Acacia", mixed.plant.sub$species)]))

write.csv(height.planting, "planting species heights.csv")
write.csv(sla.planting, "planting species sla.csv")
write.csv(seed.mass.planting, "planting species seed.mass.csv")
write.csv(wood.density.planting, "planting species wood.density.csv")

sum(kmeans(best.env,centers=8)$withinss)


trait.means

trait.means.sub<-trait.means[!grepl("sp\\.|Unknown|Dead", trait.means$species) &
                               trait.means$genus.only==0,]

trait.means.sub<-trait.means.sub[trait.means.sub$species %in% 
                                   unique(mixed.plant.sub$species),]

trait.means.mat<-as.matrix(trait.means.sub[,-1])

length(trait.means.mat[!is.na(trait.means.mat)])/length(trait.means.mat)

sum(rowSums(is.na(trait.means.mat))<1)
sum(rowSums(is.na(trait.means.mat))>0 & 
      rowSums(is.na(trait.means.mat))<4)

sum(rowSums(is.na(trait.means.mat))==4)



plot(mixed.site.sub$SM.cwm  ~ jitter(mixed.site.sub$best.cluster,
                                     amount=0.1))
colnames(mixed.site.sub)


# old trait plot ####
cluster.df<-read.csv("./Data/best cluster summary.csv", header=TRUE)
trait.coefs<-coefTable(trait.env.model.nomono)[-1,]
position<-17-c(4,5,6,7, # Mass ratio
               1,2,3, # FD
               14,12,13,15,16,
               8,10,9,11,
               13.5)
labels=rownames(trait.coefs)

cbind(trait.coefs, position)

labels=c("SLA CWM", "Wood density CWM", "Seed mass CWM", "Max height CWM",
         "Functional range", "Functional evenness", "Functional divergence",
         "Age (years)", expression("Density (stems ha"^"-1"*")"),
         expression("Min. temp ("*degree*"C)"), "Aridity index",
         "Sand (%)", 
         "SLA CWM:Age", "Seed mass CWM:Density",
         "SLA CWM:Min. temp", "Seed mass CWM:Age",
         "Age:Min. temp")
         
pdf(paste0("./Plots/Trait coefficients ", Sys.Date(),".pdf"), 
    height=7, width=10.5, useDingbats=FALSE)

framemat<-rbind(c(0.15,0.45,0.1,0.95), # Coef
                c(0.51,0.67,0.725,0.95), # Int 1
                c(0.73,0.89,0.725,0.95), # Int 2
                c(0.51,0.67,0.425,0.65), # Int 3
                c(0.73,0.89,0.425,0.65), # Int 4
                
                c(0.91,0.93,0.5,0.8), # Int gradient
                c(0.51,0.73,0.1,0.35), # Cluster 1
                c(0.73,0.95,0.1,0.35), # Cluster 2
                c(0.005,0.995,0.03,0.98)) # Box

split.screen(framemat)

#                             COEFFICIENT PLOT ####

screen(1)
par(mar=c(0,0,0,0), las=1, tck=-0.02, mgp=c(3,0.2,0), ps=8)
plot(position ~ trait.coefs[,1], xlim=c(-0.3,0.55), type="n", axes=FALSE,
     xlab="", ylab="")
rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
     border=NA, col="white")

text(x=0.4, y=15, label="Functional\ndiversity", 
     font=2, col="grey60")

rect(xleft=-1, xright=1,
     ytop=13.5, ybottom=9.5, col="grey90", 
     border=NA)
text(x=0.4, y=11.5, label="Mass ratio", 
     font=2, col="grey60")

text(x=0.4, y=7.5, label="Mass ratio\nx\nEnvironment\nand structure", 
     font=2, col="grey60")

rect(xleft=-1, xright=1,
     ytop=5.5, ybottom=0, col="grey90", 
     border=NA)
text(x=-0.175, y=2, label="Environment\nand structure", 
     font=2, col="grey60")

points(position ~ trait.coefs[,1], pch=16, cex=0.75,
       col=c(rep("black", 12), rep("grey50", 5)))

text(position[c(1,2)]+0.252 ~ trait.coefs[c(1,2),1], 
     labels=c("(f)","(g)"),
     col="black")

text(position[13:16]+0.225 ~ trait.coefs[13:16,1], 
     labels=c("(b)","(d)","(c)","(e)"), col="grey50")

box()
segments(x0=0, x1=0, y0=0, y1=18, lty="dashed")
segments(y0=position, y1=position,
         x0=trait.coefs[,1]-1.98*trait.coefs[,2],
         x1=trait.coefs[,1]+1.98*trait.coefs[,2],
         col=c(rep("black", 12), rep("grey50", 5)))

axis(side=1, mgp=c(3,0.2,0))
axis(side=2, at=position[1:12], labels=labels[1:12], las=1, mgp=c(3,0.5,0))
axis(side=2, at=position[13:17], labels=labels[13:17], las=1, 
     col.axis="grey50", col.ticks="grey50", mgp=c(3,0.5,0))
mtext(side=1, text="Standardised effect size", line=0.9)

text(x=relative.axis.point(0.05, "x"), 
     y=relative.axis.point(0.975, "y"),
     labels="(a)", font=2, cex=1.25)
close.screen(1)

#                             INTERACTION PLOTS ####

screen(2)
par(mar=c(0,0,0,0), ps=8, las=1, tck=-0.02, mgp=c(3,0.4,0))
temp.inter.list<-interaction.matrix.fun(trait.env.model.nomono, 
                                        x="SLA.cwm", y="log.age", 
                                        grid.size=150, quantile=0.01)
image(temp.inter.list, zlim=c(8.4,12.3),
      col=colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150), axes=F, las=1)

contour(x=temp.inter.list,add=T, 
        levels=log(c(25, 50, 75, 100, 125, 150)*1000),
        labels=c(25, 50, 75, 100, 125, 150),
        method="edge", lwd=0.5)

axis(side=1, mgp=c(3,0,0), at=(c(1.5,2,2.5,3,3.5,4)-mean(mixed.site.sub$SLA.cwm))/
                               sd(mixed.site.sub$SLA.cwm),
     labels=c(1.5,2,2.5,3,3.5,4))
mtext(side=1, text=expression("ln(SLA (mm"^"2"*" g"^"-1"*")) CWM"), line=0.85)

axis(side=2, at=(log(c(5,10,15,25))-mean(mixed.site.sub$log.age))/sd(mixed.site.sub$log.age),
     labels=c(5,10,15,25))
mtext(side=2, text="Age (years)", line=1.25, las=0)

mtext(side=3, at=relative.axis.point(0, "x"), text="(b)", font=2, 
      line=0, cex=1.25)
box()
close.screen(2)

screen(3)
par(mar=c(0,0,0,0), ps=8, las=1, tck=-0.02, mgp=c(3,0.4,0))
temp.inter.list<-interaction.matrix.fun(trait.env.model.nomono, 
                                        "SLA.cwm", "tmin.annual",
                                        grid.size=150, quantile=0.05)
image(temp.inter.list, zlim=c(8.4,12.3),
      col=colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150), axes=F, las=1)

#with(trait.env.model.nomono$data, points(tmin.annual ~ SLA.cwm))

contour(x=temp.inter.list,add=T, 
        levels=log(c(25, 50, 75, 100, 125, 150)*1000),
        labels=c(25, 50, 75, 100, 125, 150),
        method="edge", lwd=0.5)

axis(side=1, mgp=c(3,0,0), at=(c(1.5,2,2.5,3,3.5,4)-mean(mixed.site.sub$SLA.cwm))/
       sd(mixed.site.sub$SLA.cwm),
     labels=c(1.5,2,2.5,3,3.5,4))
mtext(side=1, text=expression("ln(SLA (mm"^"2"*" g"^"-1"*")) CWM"), line=0.85)

axis(side=2, las=1,
     at=(c(5,7.5,10,12.5,15,17.5)-mean(mixed.site.sub$tmin.annual))/
       sd(mixed.site.sub$tmin.annual),
     labels=c(5,7.5,10,12.5,15,17.5))
mtext(side=2, text=expression("Min temp ("*degree*"C)"), line=1.25, las=0)

mtext(side=3, at=relative.axis.point(0, "x"), text="(c)", font=2, line=0, cex=1.25)
box()
close.screen(3)

screen(4)
par(mar=c(0,0,0,0), ps=8, las=1, tck=-0.02, mgp=c(3,0.4,0))
temp.inter.list<-interaction.matrix.fun(trait.env.model.nomono, 
                                        "SM.cwm", "log.density",
                                        grid.size=150, quantile=0.025)
image(temp.inter.list, zlim=c(8.4,12.3),
      col=colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150), axes=F, las=1)

#with(trait.env.model.nomono$data, points(log.density ~ SM.cwm))

contour(x=temp.inter.list,add=T, 
        levels=log(c(25, 50, 75, 100, 125, 150)*1000),
        labels=c(25, 50, 75, 100, 125, 150),
        method="edge", lwd=0.5)

axis(side=1, mgp=c(3,0,0), at=(c(0,1,2,3)-mean(mixed.site.sub$SM.cwm))/
                               sd(mixed.site.sub$SM.cwm),
     labels=c(0,1,2,3))
mtext(side=1, text="ln(Seed mass (mg)) CWM", line=0.85)

axis(side=2, at=(log(c(250,1000,2500,5000))-mean(mixed.site.sub$log.density, na.rm=TRUE))/sd(mixed.site.sub$log.density, na.rm=TRUE),
     labels=c(0.25,1,2.5,5))
mtext(side=2, text=expression("Density (stems"^3*" ha"^"-1"*")"),
      line=1.25, las=0)

mtext(side=3, at=relative.axis.point(0, "x"), text="(d)", font=2, line=0, cex=1.25)
box()
close.screen(4)

screen(5)
par(mar=c(0,0,0,0), ps=8, las=1, tck=-0.02, mgp=c(3,0.4,0))
temp.inter.list<-interaction.matrix.fun(trait.env.model.nomono, 
                                        "SM.cwm", "log.age",
                                        grid.size=150, quantile=0.025)
image(temp.inter.list, zlim=c(8.4,12.3),
      col=colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150), axes=F, las=1)

#with(trait.env.model.nomono$data, points(log.age ~ SM.cwm))

contour(x=temp.inter.list,add=T, 
        levels=log(c(25, 50, 75, 100, 125, 150)*1000),
        labels=c(25, 50, 75, 100, 125, 150),
        method="edge", lwd=0.5)

axis(side=1, mgp=c(3,0,0), at=(c(0,1,2,3)-mean(mixed.site.sub$SM.cwm))/
       sd(mixed.site.sub$SM.cwm),
     labels=c(0,1,2,3))
mtext(side=1, text="ln(Seed mass (mg)) CWM", line=0.85)

axis(side=2, at=(log(c(5,10,15,25))-mean(mixed.site.sub$log.age))/sd(mixed.site.sub$log.age),
     labels=c(5,10,15,25))
mtext(side=2, text="Age (years)", line=1.25, las=0)
mtext(side=3, at=relative.axis.point(0, "x"), text="(e)", font=2, line=0, cex=1.25)
box()
close.screen(5)


screen(6)
library(plotrix)
par(mar=c(0,0,0,0), ps=8, las=1, tck=-0.25, mgp=c(3,0.4,0))
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(8.4,12.3), xaxs="i", yaxs="i", axes=FALSE)

colours<-colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150)
rect.coords<-seq(8.4, 12.3, length.out=150)

lapply(2:length(rect.coords), function(x){
  rect(xleft=0, xright=1, ytop=rect.coords[x], ybottom=rect.coords[x-1],
       border=NA, col=colours[x])
})
axis(side=4, at=log(c(5, 10, 25, 50, 75, 100, 150, 200)*1000), 
     labels=c(5, 10, 25, 50, 75, 100, 150, 200))
mtext(side=4, text=expression("Planting biomass (t ha"^-1*")"), las=3, line=2)
box()
close.screen(6)

#                             CLUSTER SLOPES ####

# 
plot.data<-trait.cluster.model.nomono$data

screen(7) # SLA.cwm
par(mar=c(0,0,0,0), ps=8, las=1, tck=-0.02, mgp=c(3,0.4,0))
with(plot.data, plot(residuals ~ SLA.cwm, pch=16,
                     col="grey90",
                     axes=FALSE, xlab="", ylab=""))
box()
axis(side=2)
     
axis(side=1, mgp=c(3,0.2,0), at=(c(1:4)-mean(mixed.site.sub$SLA.cwm))/
       sd(mixed.site.sub$SLA.cwm),
     labels=1:4)

mtext(side=1, text=expression("ln(SLA (mm"^"2"*" g"^"-1"*")) CWM"), 
      line=1.1)
mtext(side=2, text="Biomass residuals", line=1.5, las=0)

text(x=relative.axis.point(0.065, "x"), 
     y=relative.axis.point(0.93, "y"),
      labels="(f)", font=2, cex=1.25)

for (i in 1:length(unique(plot.data$best.cluster.fact))){
  
cluster.data.sub<-plot.data[plot.data$best.cluster.fact==i,]
cluster.random<-data.frame(planting=plot.data$planting, 
                           IBRAsub=plot.data$IBRAsub)

cluster.random<-cluster.random[!duplicated(cluster.random$planting),]
model.coefs<-colnames(attr(trait.cluster.model$terms, "factors"))
model.coefs<-model.coefs[!grepl(":", model.coefs)]

prediction.df<-as.data.frame(matrix(0, nrow=100, ncol=length(model.coefs)))
colnames(prediction.df)<-model.coefs
prediction.df[,"SLA.cwm"]=seq(min(cluster.data.sub[, "SLA.cwm"]),
                      max(cluster.data.sub[, "SLA.cwm"]),
                      length.out=100)
prediction.df[,"best.cluster.fact"]=as.factor(i)

predictions.temp<-predict(object=trait.cluster.model, 
                          newdata=prediction.df, level=0, se.fit=TRUE)
predictions<-list(fit=predictions.temp$fit,
                  x=prediction.df[,1],
                  se=predictions.temp$se)
if(i==8){
polygon(y=c(predictions$fit+1.96*predictions$se,
            rev(predictions$fit-1.96*predictions$se),
            predictions$fit[1]+1.96*predictions$se[1]),
        x=c(predictions$x, rev(predictions$x), predictions$x[1]), 
        col=with(cluster.df, rgb(0.47, 0.67, 0.46, 0.75)), 
        border=NA)


points(predictions$fit ~ predictions$x, type="l", lwd=2, 
       col=with(cluster.df, rgb(red[i], green[i], blue[i], 1)))
next
}

#polygon(y=c(predictions$fit+1.96*predictions$se,
#              rev(predictions$fit-1.96*predictions$se),
#              predictions$fit[1]+1.96*predictions$se[1]),
#          x=c(predictions$x, rev(predictions$x), predictions$x[1]), 
#          col=with(cluster.df, rgb(red[i], green[i], blue[i], 0.05)), 
#        border=NA)

  points(predictions$fit ~ predictions$x, type="l", lwd=2, lty="31",
         col=with(cluster.df, rgb(red[i], green[i], blue[i], 0.75)))
}
close.screen(7)

screen(8) # WD.cwm
par(mar=c(0,0,0,0), ps=8, las=1, tck=-0.02, mgp=c(3,0.4,0))
with(plot.data, plot(residuals ~ WD.cwm, pch=16,
                     col="grey90",
                     axes=FALSE, xlab="", ylab=""))
box()
axis(side=2, labels=NA)

axis(side=1, mgp=c(3,0.2,0), at=(c(1:4)-mean(mixed.site.sub$WD.cwm))/
       sd(mixed.site.sub$WD.cwm),
     labels=1:4)

mtext(side=1, text=expression("ln(Wood density (g cm"^"-3"*")) CWM"), 
      line=1.1)
text(x=relative.axis.point(0.065, "x"), 
     y=relative.axis.point(0.93, "y"),
     labels="(g)", font=2, cex=1.25)

for (i in c(1,3:length(unique(plot.data$best.cluster.fact)),2)){
  
  cluster.data.sub<-plot.data[plot.data$best.cluster.fact==i,]
  cluster.random<-data.frame(planting=plot.data$planting, 
                             IBRAsub=plot.data$IBRAsub)
  
  cluster.random<-cluster.random[!duplicated(cluster.random$planting),]
  model.coefs<-colnames(attr(trait.cluster.model$terms, "factors"))
  model.coefs<-model.coefs[!grepl(":", model.coefs)]
  
  prediction.df<-as.data.frame(matrix(0, nrow=100, ncol=length(model.coefs)))
  colnames(prediction.df)<-model.coefs
  prediction.df[,"WD.cwm"]=seq(min(cluster.data.sub[, "WD.cwm"]),
                                max(cluster.data.sub[, "WD.cwm"]),
                                length.out=100)
  prediction.df[,"best.cluster.fact"]=as.factor(i)
  
  predictions.temp<-predict(object=trait.cluster.model, 
                            newdata=prediction.df, level=0, se.fit=TRUE)
  predictions<-list(fit=predictions.temp$fit,
                    x=prediction.df[,"WD.cwm"],
                    se=predictions.temp$se)
  
  if(i==2){
    polygon(y=c(predictions$fit+1.96*predictions$se,
                rev(predictions$fit-1.96*predictions$se),
                predictions$fit[1]+1.96*predictions$se[1]),
            x=c(predictions$x, rev(predictions$x), predictions$x[1]), 
            col=with(cluster.df, rgb(0.74, 0.94, 0.53, 0.85)), 
            border=NA)

    points(predictions$fit ~ predictions$x, type="l", lwd=2,
           col=with(cluster.df, rgb(red[i], green[i], blue[i], 1)))
    next
  }
  
   points(predictions$fit ~ predictions$x, type="l", lwd=2, lty="31",
         col=with(cluster.df, rgb(red[i], green[i], blue[i], 0.75)))
}
close.screen(8)

screen(9)
par(mar=c(0,0,0,0))
box(lwd=1.5)
close.screen(9)
close.screen(all.screens=TRUE)

dev.off()

 ####

#             FUNCTIONAL DIVERSITY REGRESSION SLOPES ####
summary(fun.div.env)
plot.data<-fun.div.env$data

summary(plot.data$FEm)

framemat<-rbind(c(0.1,0.4,0.2,0.95),
                c(0.4,0.7,0.2,0.95),
                c(0.7,0.99,0.2,0.95))

pdf("./Plots/fun div regression slopes.pdf", height=2.25, width=6, useDingbats=FALSE)

split.screen(framemat)

# get residuals from environment model
env.plot.data<-env.model$data
env.plot.data$residuals<-residuals(env.model, level=0)

env.plot.data$pl.p<-mixed.site.sub[,"pl.p"][match(paste0(env.plot.data$biomass.area,
                                                    env.plot.data$planting),
                                             paste0(mixed.site.sub$biomass.area,
                                                    mixed.site.sub$planting))]

env.plot.data<-merge(env.plot.data, mixed.site.sub[,c("pl.p","FRv","FEm","FDm")])

fun.div.resid.model<-lme(residuals ~ FRv + FEm + FDm,
                random=~1|IBRAsub/planting, data=env.plot.data)

lapply(1:3, function(num){

  screen(num)
  par(mar=c(0,0,0,0), las=1, tck=-0.02, mgp=c(3,0.4,0), ps=10)
  
  x<-c("FRv","FEm","FDm")[num]

  new.dat<-as.data.frame(cbind(seq(min(env.plot.data[,x]),
                          c(16.5,
                            rep(max(env.plot.data[,x]),2))[num], 
                          length.out=100),
                          rep(mean(env.plot.data[,c("FRv","FEm","FDm")[-num][1]]), 100),
                          rep(mean(env.plot.data[,c("FRv","FEm","FDm")[-num][2]]), 100)))
                     
  colnames(new.dat)<-c(x, c("FRv","FEm","FDm")[-num])
  
  temp.pred<-predict(fun.div.resid.model, new.dat, se.fit=TRUE, level=0)
  
  if(num ==1){
  with(env.plot.data, plot(residuals ~ env.plot.data[,x],
                           xlim=c(0,17),
                           axes=FALSE, pch=16, col="grey90", xlab="", ylab=""))
  }
  
  if(num == 2){
    with(env.plot.data, plot(residuals ~ env.plot.data[,x],
                             xlim=c(-0.05,1.05),
                             axes=FALSE, pch=16, col="grey90", xlab="", ylab=""))
  }
  
  if(num == 3){
    with(env.plot.data, plot(residuals ~ env.plot.data[,x],
                             xlim=c(-0.05,1),
                             axes=FALSE, pch=16, col="grey90", xlab="", ylab=""))
  }
  

  polygon(x=c(new.dat[,1], rev(new.dat[,1])),
          y=c(temp.pred$fit + 1.96*temp.pred$se.fit,
            rev(temp.pred$fit - 1.96*temp.pred$se.fit)),
          col=rgb(0.2,0.2,0.2,0.5), border=NA)
  
  points(temp.pred$fit ~ new.dat[,x], type="l")

  if(num==1){
    axis(side=2, at=c(-2:2),
         labels=-2:2)
    mtext(side=2, text=expression("Residual biomass ln(kg ha"^-1*")"),
          line=1.5, las=0)
  } else{ 
    axis(side=2, at=c(-2:2), labels=NA)}
  
    axis(side=1, at=pretty(env.plot.data[,x], n=c(5,3,3)[num]),
         mgp=c(3,0,0))
  
    mtext(side=1, text=c("Functional range",
                          "Functional evenness",
                          "Functional divergence")[num],
         line=1)
   
   text(x=relative.axis.point(0.065, "x"),
        y=relative.axis.point(0.935, "y"),
        labels=c("(a)","(b)","(c)")[num], font=2)
   
   box()
   close.screen(num)
          
})
close.screen(all.screens=TRUE)
dev.off()


pdf("./Plots/fun div regression slopes no points.pdf", height=2.25, width=6)

split.screen(framemat)

# get residuals from environment model
env.plot.data<-env.model$data
env.plot.data$residuals<-residuals(env.model, level=0)

env.plot.data$pl.p<-mixed.site.sub[,"pl.p"][match(paste0(env.plot.data$biomass.area,
                                                         env.plot.data$planting),
                                                  paste0(mixed.site.sub$biomass.area,
                                                         mixed.site.sub$planting))]

env.plot.data<-merge(env.plot.data, mixed.site.sub[,c("pl.p","FRv","FEm","FDm")])

fun.div.resid.model<-lme(residuals ~ FRv + FEm + FDm,
                         random=~1|IBRAsub/planting, data=env.plot.data)

lapply(1:3, function(num){
  
  screen(num)
  par(mar=c(0,0,0,0), las=1, tck=-0.02, mgp=c(3,0.4,0), ps=8)
  
  x<-c("FRv","FEm","FDm")[num]
  
  new.dat<-as.data.frame(cbind(seq(min(env.plot.data[,x]),
                                   c(quantile(env.plot.data[,x],0.99),
                                     rep(max(env.plot.data[,x]),2))[num], 
                                   length.out=100),
                               rep(mean(env.plot.data[,c("FRv","FEm","FDm")[-num][1]]), 100),
                               rep(mean(env.plot.data[,c("FRv","FEm","FDm")[-num][2]]), 100)))
  
  colnames(new.dat)<-c(x, c("FRv","FEm","FDm")[-num])
  
  temp.pred<-predict(fun.div.resid.model, new.dat, se.fit=TRUE, level=0)
  
  with(env.plot.data, plot(NULL,
                           xlim=c(quantile(env.plot.data[,x],0),
                                  
                                  c(quantile(env.plot.data[,x],0.99),
                                    1.05,1)[num]),
                           ylim=c(-0.5,0.5),
                           axes=FALSE, pch=16, col="grey90", xlab="", ylab=""))
  
  abline(summary(fun.div.resid.model)$tTable[1,1], 0, lty="dashed")
  
  polygon(x=c(new.dat[,1], rev(new.dat[,1])),
          y=c(temp.pred$fit + 1.96*temp.pred$se.fit,
              rev(temp.pred$fit - 1.96*temp.pred$se.fit)),
          col=rgb(0.2,0.2,0.2,0.5), border=NA)
  
  points(temp.pred$fit ~ new.dat[,x], type="l")
  
  if(num==1){
    axis(side=2, at=c(-0.5,-0.25,0,0.25,0.5),
         labels=c(-0.5,-0.25,0,0.25,0.5))
    mtext(side=2, text=expression("Residual ln(biomass (t ha"^-1*"))"),
          line=1.5, las=0)
  } else{ 
    axis(side=2, at=c(-0.5,-0.25,0,0.25,0.5), labels=NA)}
  
  axis(side=1, at=pretty(env.plot.data[,x], n=c(5,3,3)[num]),
       mgp=c(3,0,0))
  
  mtext(side=1, text=c("Functional range",
                       "Functional evenness",
                       "Functional divergence")[num],
        line=1)
  
  text(x=relative.axis.point(0.05, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=c("a.","b.","c.")[num], font=2)
  
  box()
  close.screen(num)
  
})
close.screen(all.screens=TRUE)
dev.off()

summary((mixed.site.sub$FEm-mean(mixed.site.sub$FEm))/sd(mixed.site.sub$FEm))

# BOXPLOT OF CLUSTER AGES
cluster.df<-read.csv("./Data/best cluster summary.csv", header=TRUE)

pdf("./Plots/age replication by cluster.pdf", height=4, width=6, useDingbats=FALSE)
par(mar=c(3,8,1,1), las=1, tck=-0.02, mgp=c(3,0.2,0), ps=8)

plot(x=NULL, y=NULL, ylim=c(1,8), xlim=summary(mixed.site.sub$age)[c(1,6)],
     xlab="", ylab="", axes=FALSE)

segments(x0=17, x1=17, y0=0, y1=9, lty="dashed")
for (i in 1:8){
  
cluster.data.sub<-mixed.site.sub[mixed.site.sub$best.cluster==i,]
  
with(cluster.data.sub, points(x=age, y=jitter(rep(i, length(age)), amount=0.2), 
                              pch=16, cex=0.8,
                              col=with(cluster.df, rgb(red[i], 
                                                       green[i], 
                                                       blue[i], 
                                                       0.6))))
}

axis(side=1, mgp=c(3,0,0))
mtext(side=1, "Age (years)", line=1.25)

axis(side=2, at=1:8, labels=cluster.df$name, mgp=c(3,0.4,0))


box()

dev.off()

#             PLOTS CLUSTERS ON AUS MAP ####

ausmap<-readShapeSpatial(paste(ifelse(Sys.info()['sysname']=="Linux", "/home/Storage HDD", "F:/"),
                               "/University files/Shape files/Shape files/Australia outlines/aust_cd66states.shp", sep=""))

mixed.planting.sub<-mixed.site.sub[!duplicated(mixed.site.sub$pl.p),]

cluster.df<-read.csv("./Data/best cluster summary.csv", header=TRUE)

plot(ausmap, lwd=0.5, xlim=c(115,155), ylim=c(-45,-9), axes=FALSE, border="grey60")
map.usr<-par("usr")
# 13 long, 10 tall

pdf(paste0("./Plots/Clusters on map ", Sys.Date(), ".pdf"), width=6, height=5.5, 
    useDingbats=FALSE)
par(mar=c(3,3,3,3), oma=c(0,0,0,0), ps=8, las=1, tck=-0.0075, mgp=c(3,0.4,0))

plot(ausmap, lwd=1, axes=FALSE, xlim=c(113,155))

axis(side=1, at=seq(110,150,10), labels=parse(text=paste(seq(110,150,10), "*degree~E", sep="")))
axis(side=2, at=seq(-45,-10,5), labels=parse(text=paste(seq(-45,-10,5), "*degree~S", sep="")),
     mgp=c(3,0.7,0))
box(lwd=1)

# plot the points for each cluster and connecting lines
lapply(cluster.df$cluster, function(x){
  cluster.sub<-mixed.planting.sub[mixed.planting.sub$best.cluster==x,]
  
  with(cluster.sub, points(lat ~ long, pch=21, lwd=0.5, cex=0.8, 
                           bg=with(cluster.df, rgb(red[x], green[x], blue[x], 1)))) 
})
dev.off()

#             PLOT CLUSTERING PLOT ####
cluster.df<-read.csv("./Data/best cluster summary.csv", header=TRUE)

best.env<-mixed.site.sub[,best.covariates[!best.covariates %in% 
                                          c("log.age", "log.density")]]

best.pca<-eigen(cor(best.env))
best.pca$values/sum(best.pca$values)

best.env.pca<-(as.matrix(best.env) %*% best.pca$vectors)/10

axis1<-tapply(best.env.pca[,1], mixed.site.sub$best.cluster, mean)
axis2<-tapply(best.env.pca[,2], mixed.site.sub$best.cluster, mean)

pdf(paste0("./Plots/best cluster pca ", Sys.Date(), ".pdf"), height=3.5, width=5)
par(mar=c(2.25,2.25,1,1), ps=8, las=1, mgp=c(3,0.3,0), tck=-0.01)
plot(best.env.pca[,2] ~ best.env.pca[,1], type="n", xlab="", ylab= "", axes=FALSE)

arrows(x0=0, 
       y0=0, 
       x1=best.pca$vectors[,1]*3, 
       y1=best.pca$vectors[,2]*3, lwd=3, col="grey60")
text(x=best.pca$vectors[,1]*3+c(0.05,0,0), 
     y=best.pca$vectors[,2]*3+c(0.15,0,-0.25), 
     labels=c("Sand (%)", expression("Tmin ("*degree*"C)"), "Moisture\navailability"), 
     pos=c(3,1,2), col="grey60", font=2)
lapply(1:8, function(x){
  cluster.sub<-best.env.pca[mixed.site.sub$best.cluster==x,]
  points(cluster.sub[,2] ~ cluster.sub[,1], pch=21, lwd=0.5, cex=0.8, 
         bg=with(cluster.df, rgb(red[x], green[x], blue[x], 1),5)) 
})

axis(side=2)
axis(side=1, mgp=c(3,0,0))
mtext(side=1, text="PCA 1", line=1)
mtext(side=2, text="PCA 2", line=1, las=0)
box()
dev.off()

library(plot3D)
best.env<-as.data.frame(best.env)
with(best.env, plot(sand ~ log.aridity.index))

pdf("./Plots/clusters by env.pdf", height=4, width=5)
par(mar=c(0,3,3,0), mfrow=c(2,2), las=1, ps=8, tck=-0.025,
    mgp=c(0,0.4,0))

# sand vs tmin
with(best.env, plot(sand ~ tmin.annual, type="n", 
                    axes=FALSE, xlab="", ylab=""))
box()
axis(side=1, labels=NA)
axis(side=2, las=1)
mtext(side=2, text="Sand (%)", las=0, line=1.75)

axis(side=3, las=1)
mtext(side=3, text=expression("Min temp ("*degree*"C)"), las=0, line=1.5)

lapply(1:8, function(x){
  cluster.sub<-best.env[mixed.site.sub$best.cluster==x,]
  points(cluster.sub[,"sand"] ~ cluster.sub[,"tmin.annual"], 
         pch=21, lwd=0.5, cex=0.8, 
         bg=with(cluster.df, rgb(red[x], green[x], blue[x], 1),5)) 
})

par(mar=c(0,0,3,3))
# sand vs aridity
with(best.env, plot(sand ~ log.aridity.index, 
                    type="n", axes=FALSE))
box()
axis(side=2, labels=NA)
axis(side=1, las=1)
mtext(side=1, 
      text="Moisture availability", las=0, line=1.5)
axis(side=4, las=1)
mtext(side=4, text="Sand (%)", las=0, line=1.5)

lapply(1:8, function(x){
    cluster.sub<-best.env[mixed.site.sub$best.cluster==x,]
  points(cluster.sub[,"sand"] ~ cluster.sub[,"log.aridity.index"], 
         pch=21, lwd=0.5, cex=0.8, 
         bg=with(cluster.df, rgb(red[x], green[x], blue[x], 1),5)) 
})

par(mar=c(3,3,0,0))
# aridity vs tmin
with(best.env, plot(log.aridity.index ~ tmin.annual,
                    type="n", axes=FALSE, xlab="", ylab=""))
box()
axis(side=2)
mtext(side=2, 
      text="Moisture availability", las=0, line=1.75)

axis(side=1, las=1)
mtext(side=1, text=expression("Min temp ("*degree*"C)"), las=0, line=1.5)

lapply(1:8, function(x){
  cluster.sub<-best.env[mixed.site.sub$best.cluster==x,]
  points(cluster.sub[,"log.aridity.index"] ~ cluster.sub[,"tmin.annual"], 
         pch=21, lwd=0.5, cex=0.8, 
         bg=with(cluster.df, rgb(red[x], green[x], blue[x], 1),5)) 
})

par(mar=c(0,0,0,0))
plot.new()
par(xpd=NA)
legend(y=relative.axis.point(0.7, "y"),
       x=relative.axis.point(0.2, "x"),
       legend=paste0(cluster.df$name, " (", cluster.df$code, ")"),
       pch=21, pt.bg=rgb(cluster.df$red, cluster.df$green, cluster.df$blue),
       y.intersp=0.6, bty="n")
par(xpd=TRUE)

dev.off()

#             MISC PLOTS AND EXPLORATION ####

#                             SAND ####
cor(mixed.site.sub[,colnames(mixed.site.sub) %in% covariates], 
             use="complete.obs")

with(plot.data, plot(log(biomass.area) ~ sand, type="n"))

sapply(1:8, function(x){
  
  temp<-plot.data[plot.data$best.cluster==x,]
  
  with(temp, points(log(biomass.area) ~ sand, pch=21,
                    bg=rgb(cluster.df$red[x],
                           cluster.df$green[x],
                           cluster.df$blue[x],
                           1)))
})    

# looks like the Medit shrublands and coastal plains have high sand, and seem
# to have slightly higher residual biomass than the other clusters. It might be
# useful to look at a boxplot of average residuals for each cluster for their
# mean sand content

meansand<-sapply(1:8, function(x){
  mean(plot.data$sand[plot.data$best.cluster==x], na.rm=TRUE)
})

with(plot.data, plot(log(biomass.area) ~ sand, type="n"))
sapply(1:8, function(x){
  
  temp<-plot.data[plot.data$best.cluster==x,]
  
  with(temp, boxplot(log(biomass.area) ~ best.cluster, 
                    col=rgb(cluster.df$red[x],
                           cluster.df$green[x],
                           cluster.df$blue[x],
                           1),
                    at=meansand[x],
                    add=TRUE,
                    width=3,
                    na.action="na.omit"))
})  

# so there's 3 clusters of high sand content. Let's see how well biomass in these
# clusters is explained by our full model (minus sand)
model.data<-env.model$data
sandless.model<-update(env.model, .~. -sand)
plot.data<-sandless.model$data
plot.data$residuals<-residuals(sandless.model, level=0)

cluster.df<-read.csv("./Data/best cluster summary.csv", header=TRUE)
cluster.df<-as.data.frame(sapply(cluster.df, function(x){gsub("\\+AC0", "", x)}))
cluster.df[,c("red","green","blue")]<-sapply(cluster.df[,c("red","green","blue")],
                                             function(x){as.numeric(as.character(x))})

with(plot.data, plot(residuals ~ sand, type="n"))

sapply(1:8, function(x){
  
  temp<-plot.data[plot.data$best.cluster==x,]
  
  with(temp, points(residuals ~ sand, pch=21,
                    bg=rgb(cluster.df$red[x],
                           cluster.df$green[x],
                           cluster.df$blue[x],
                           1)))
})   

# looks like those high sand clusters have higher residuals than the other groups
# let's try our boxplot again
meansand<-sapply(1:8, function(x){
  mean(plot.data$sand[plot.data$best.cluster==x])
})

library(plotrix)

pdf(paste0("./Plots/sand clusters ", Sys.Date(), ".pdf"), width=5, height=4)
par(mar=c(3,3,1,1), ps=8, tck=-0.01, mgp=c(3,0.5,0), las=1)
with(plot.data, plot(sapply(split(plot.data$residuals, f=plot.data$best.cluster),
                            mean) ~ meansand, type="n", axes=FALSE, 
                     xlab="", ylab="", ylim=c(-0.55,0.55), xlim=c(-1.75, 1.55)))

segments(x0=-4, x1=3, y0=0, y1=0, lty="dashed", col="grey50")

draw.circle(x=mean(meansand[c(3,6,7)]),
            y=mean(plot.data$residuals[plot.data$best.cluster %in% c(3,6,7)]),
            radius=0.45, border="red", lwd=2)
text(x=mean(meansand[c(3,6,7)]), y=0.025, 
            labels="Mediterranean climates", col="red")

sapply(1:8, function(x){
  
  temp<-plot.data[plot.data$best.cluster==x,]
  
  resid.se<-sd(temp$residuals)/sqrt(length(temp$residuals))
  sand.se<-sd(temp$sand)/sqrt(length(temp$sand))
  
  # residuals CIs
  arrows(x0=meansand[x],
         x1=meansand[x],
         y0=mean(temp$residuals)+1.96*resid.se,
         y1=mean(temp$residuals)-1.96*resid.se,
         length=0.02, angle=90, code=3)
  
  # sand CIs
  arrows(x0=mean(temp$sand)+1.96*sand.se,
         x1=mean(temp$sand)-1.96*sand.se,
         y0=mean(temp$residuals),
         y1=mean(temp$residuals),
         length=0.02, angle=90, code=3)
  
  # points
  with(temp, points(mean(residuals) ~ meansand[x],
                     bg=rgb(cluster.df$red[x],
                             cluster.df$green[x],
                             cluster.df$blue[x],
                             1),
                     pch=21))
  
  # text
  if(x==3){
  with(temp, text(y=mean(residuals)+0.04, x=meansand[x]-0.05,
                  labels=cluster.df$code[x], adj=1))
    return(NULL)
  }
  
  with(temp, text(y=mean(residuals)+0.04, x=meansand[x]+0.05,
                  labels=cluster.df$code[x], adj=0))

})  

legend(y=relative.axis.point(1, "y"),
       x=relative.axis.point(0, "x"),
       legend=paste0(cluster.df$full.name, " (", cluster.df$code, ")"),
       pch=21, pt.bg=rgb(cluster.df$red, cluster.df$green, cluster.df$blue),
       bty="n", y.intersp=0.6)

axis(side=1, mgp=c(3,0,0), at=(c(30,40,50,60,70,80,90)-
                                 mean(mixed.site.sub$sand))/sd(mixed.site.sub$sand),
     labels=c(0.3,0.4,0.5,0.6,0.7,0.8,0.9))
axis(side=2)
mtext(side=2, text=expression("Residual productivity ln(t ha"^-1*" year"^-1*")"), 
      line=1.75, las=0)
mtext(side=1, text="Soil sand content (%)", line=0.9)

box()
dev.off()

#                             SEMI-ARID WOODLANDS RICHNESS CORRELATION ####

#                                             SUBSET AND EXPLORE DATA ####

plot.data<-richness.only.cluster.model$data
plot.data$pl.p<-mixed.site.sub$pl.p[match(rownames(plot.data),
                                          rownames(mixed.site.sub))]
plot.data$plot.area<-mixed.site.sub$plot.area[match(rownames(plot.data),
                                          rownames(mixed.site.sub))]

clus2.site<-plot.data[plot.data$best.cluster==7,]
clus2.plant<-droplevels(mixed.plant.sub[mixed.plant.sub$pl.p %in%
                             clus2.site$pl.p, ])

with(clus2.site, plot(residuals ~ rare.rich))

clus2.site$pl.p[order(clus2.site$residuals, decreasing=TRUE)]

clus2.plant[clus2.plant$pl.p=="CN921.1",
            c("species", "tot.ag.mass")]
mixed.site.sub[mixed.site.sub$pl.p=="CN921.1",]

library(ineq)
clus2.biomass.var<-tapply(clus2.plant$tot.ag.mass, clus2.plant$pl.p, Gini) 
                          
clus2.biomass.var<-data.frame(mass.var=clus2.biomass.var,
                              pl.p=names(clus2.biomass.var))

clus2.site<-merge(clus2.site, clus2.biomass.var, all.x=TRUE)
hist(clus2.site$mass.var)

clus2.site$rich.raw<-mixed.site.sub$rare.rich[match(clus2.site$pl.p,
                                                    mixed.site.sub$pl.p)]
hist(clus2.site$rich.raw, breaks=10)
clus2.site$rich.split<-cut(clus2.site$rich.raw, 
                           breaks=c(0,1,3,5,13))

cor(cbind(clus2.site[,c("residuals","mass.var", "rich.raw", "plot.area")]),
    use="complete.obs")

# Let's sort stems in each plot, look at what is contributing the most biomass
# (proportional) and try to extract particular species, genera or traits

clus2.plant.list<-split(clus2.plant, f=clus2.plant$pl.p)

clus2.list.sort<-lapply(clus2.plant.list, function(x){
  temp<-x[order(x[,"tot.ag.mass"], decreasing=TRUE),]
  
  equal.mass<-sum(temp$tot.ag.mass)/length(temp$species)
  
  temp$prop.mass<-temp$tot.ag.mass/equal.mass
  return(temp)
})

clus2.plant<-do.call("rbind", clus2.list.sort)

summary(clus2.plant$prop.mass)
hist(clus2.plant$prop.mass, breaks=20)

length(clus2.plant$prop.mass[clus2.plant$prop.mass>2])/
  length(clus2.plant$prop.mass)

sort(table(droplevels(clus2.plant$species[clus2.plant$prop.mass>3])),
     decreasing=TRUE)

clus2.plant<-merge(clus2.plant, mixed.site.sub[,c("pl.p","age")], all.x=TRUE)
clus2.plant<-merge(clus2.plant, species[,c("species","genus")], all.x=TRUE)


clus2.mixed<-droplevels(clus2.plant[clus2.plant$pl.p %in%
                           clus2.site$pl.p[clus2.site$rich.raw>1],])

table(droplevels(clus2.mixed$species[clus2.mixed$prop.mass>4]),
      clus2.mixed$age[clus2.mixed$prop.mass>4])

table(droplevels(clus2.mixed$genus[clus2.mixed$prop.mass>2]),
      clus2.mixed$age[clus2.mixed$prop.mass>2])

table(droplevels(sapply(split(clus2.mixed, f=clus2.mixed$pl.p), function(x){
  x$species[which.max(x$prop.mass)]})),
  sapply(split(clus2.mixed, f=clus2.mixed$pl.p), function(x){
    x$age[1]}))

table(droplevels(sapply(split(clus2.mixed, f=clus2.mixed$pl.p), function(x){
  x$genus[which.max(x$prop.mass)]})),
  sapply(split(clus2.mixed, f=clus2.mixed$pl.p), function(x){
    x$age[1]}))

plot(log(clus2.mixed$prop.mass) ~ clus2.mixed$MH)

#                                             PLOT ####

framemat<-rbind(c(0.1,0.55,0.2,0.9),
                c(0.625,0.95,0.2,0.9))

pdf(paste0("./Plots/cluster 7 breakdown ",Sys.Date(),".pdf"), height=3, width=8,
           useDingbats=FALSE)  

split.screen(framemat)

screen(1)
par(mar=c(0,0,0,5), las=1, tck=-0.02, mgp=c(3,0.4,0), ps=8)

with(clus2.site, plot(residuals ~ mass.var, pch=21,
                      xaxt="n",
                      bg=colorRampPalette(c("white",
                                            "burlywood1",
                                            "darkgoldenrod4"))(4)[rich.split]))
mtext(side=1, text="Gini coefficient", line=1.25)
mtext(side=2, text=expression("ln(biomass ha"^"-1"*") residuals"), line=1.5, las=0)

axis(side=1, mgp=c(3,0.2,0))

par(xpd=TRUE)
legend(x=0.9,
       y=relative.axis.point(0.5, "y"),
       legend=c("1","2-3","3-5","5-13"),
       fil=colorRampPalette(c("white",
                              "burlywood1",
                              "darkgoldenrod4"))(4),
       bty="n",
       yjust=0.5, title=expression(bold("Rarefied\nrichness")))

text(x=relative.axis.point(0.05, "x"), 
     y=relative.axis.point(0.95, "y"),
     labels="(a)", font=2, cex=1.25)

close.screen(1)       

clus2.model<-lme(residuals ~ (rich.raw + mass.var)^2,
          random=~1|IBRAsub/planting,
          data=clus2.site)
summary(clus2.model)

screen(2)
par(mar=c(0,0,0,0), las=1, tck=-0.02, mgp=c(3,0.4,0), ps=8)
temp.inter.list<-interaction.matrix.fun(clus2.model,
                                        "rich.raw", "mass.var",
                                        250, 0.01)
image(temp.inter.list,
      col=colorRampPalette(c(rgb(0.95,1,0.95,1),"darkseagreen1","darkgreen"))(150), 
      xaxt="n", las=1)
contour(x=temp.inter.list,add=T, nlevels=5)

with(clus2.site, points(mass.var ~ rich.raw))

axis(side=1, mgp=c(3,0.2,0))
mtext(side=1, text="Rarefied richness", line=1.25)

mtext(side=2, text="Gini coefficient", line=1.5, las=0)

text(x=relative.axis.point(0.05, "x"), 
     y=relative.axis.point(0.95, "y"),
     labels="(b)", font=2, cex=1.25)
box()
close.screen(2)
dev.off()

Lc(clus2.plant$tot.ag.mass[clus2.plant$pl.p=="ALWP16.a"],
   plot=TRUE)

r.squaredGLMM(clus2.model)

#                             TRAIT/FDIV BY ENV GROUP ####
head(mixed.site.sub.multi)

pdf("./Plots/richness trait etc by cluster.pdf", height=8, width=5)
par(mar=c(0,2,0,1), mfrow=c(4,2), oma=c(4.5,2,1,1), tcl=-0.25, 
    las=0, mgp=c(3,0.5,0), ps=8, las=1)
sapply(1:8, function(n){
  
  x <- c("SLA", "WD", "SM", "MH",
         "rare.rich", "FRv", "FEm", "FDm")[n]
           
  temp.data <- mixed.site.sub.multi
  
  temp.col <- rgb(cluster.df$red[temp.data$best.cluster],
                  cluster.df$green[temp.data$best.cluster],
                  cluster.df$blue[temp.data$best.cluster], 0.2)
  
          if(n %in% 1:4){
          
            site.means <- data.frame(tapply(mixed.plant.sub[,x],
                                            mixed.plant.sub$pl.p,
                                            mean))  
            site.means$pl.p <- rownames(site.means)
            colnames(site.means)[1] = 'var'
            
            site.means <- merge(site.means, mixed.site.sub.multi[,c("pl.p", "best.cluster")],
                                by.x="pl.p", by.y="pl.p", all.x=TRUE, all.y=FALSE, sort=FALSE)
            
            site.means <- site.means[!is.na(site.means$best.cluster),]
            
            plot(site.means$var ~ jitter(site.means$best.cluster, amount=0.25),
                 pch=16, col=temp.col, cex=0.5, axes=FALSE,
                 ylim=c(min(site.means$var),
                        max(site.means$var)*1.1))
            
            boxplot(site.means$var ~ site.means$best.cluster, 
                    axes=FALSE, add=TRUE, lwd=0.75)
          
            box()
            axis(side=1, at=1:8, labels=NA)  
            
            if(n==1){axis(side=2, at=log(c(2,5,10,20)), labels=c(2,5,10,20))}
            if(n==2){axis(side=2, at=log(c(0.5,0.6,0.7,0.8,0.9,1)), 
                          labels=c(0.5,0.6,0.7,0.8,0.9,1))}
            if(n==3){axis(side=2, at=log(c(1,10,100,400)), labels=c(1,10,100,400))}
            if(n==4){axis(side=2, at=log(c(2,5,10,20,50)), labels=c(2,5,10,20,50))}
            
            axis(side=2, at=log(c(seq(0.1,1,0.1),
                                  seq(1,10,1),
                                  seq(10,100,10),
                                  seq(100,1000,100))), labels=NA, tcl=-0.125)
            
          } else {
  
          temp.data$var = temp.data[,x]
          
          ylims <- c(min(temp.data$var),
                     max(temp.data$var)*1.1)
          
          if(n==6){ylims[2] = 20}
                     
          plot(temp.data$var ~ jitter(temp.data$best.cluster, amount=0.25),
                 pch=16, col=temp.col, cex=0.5, axes=FALSE,
                 ylim=ylims)
          
          boxplot(temp.data$var ~ temp.data$best.cluster, 
                  axes=FALSE, lwd=0.75, add=TRUE)
          
          axis(side=2)
          if(n %in% 7:8){
            axis(side=1, at=1:8, labels=NA)
            par(xpd=NA)
            text(x=1:8, y=relative.axis.point(-0.05, "y"),
                 labels=cluster.df$full.name, adj=1, srt=30)
            par(xpd=FALSE)
          } else {axis(side=1, at=1:8,labels=NA)}
          
          box()
          }
  
  mtext(side=2, line=1.5,
        text=list(expression("Specific leaf area (mm"^2*" g"^-2*")"),
                  expression("Wood density (g cm"^-3*")"),
                  "Seed mass (mg)",
                  "Maximum height (m)",
                  "Rarefied species richness (n = 30)",
                  "Functional range",
                  "Functional evenness",
                  "Functional divergence")[[n]], cex=0.8, las=0)
  
  text(x=relative.axis.point(0.02, "x"),
       y=relative.axis.point(0.95, "y"),
       labels=paste0("(", letters[n], ")"), font=2, adj=0)
          
         })
dev.off()



mixed.plant.sub<-merge(mixed.plant.sub, 
                        mixed.site.sub[,c("pl.p","best.cluster")],
                        all.x=TRUE)

SLA.by.group<-sapply(split(mixed.plant.sub,
                           mixed.plant.sub$best.cluster),
                     function(x){
                       unique(x$SLA)
                     })




#             RAINFALL AND DIVERSITY DENSITY PLOTS ####

mixed.site.sub.multi

pdf("./Plots/rainfall density.pdf", height=2, width=4)
par(mar=c(2.5,3,1,1), ps=8, las=1, tck=-0.02, mgp=c(3,0,0))

plot(x=NULL, y=NULL, ylim=c(0,1.2), 
     xlim=log(c(200,4000)), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
polygon(density(log(mixed.site.sub.multi$prec.annual)), col="grey90")

axis(side=1, at=log(c(200,300,400,500,1000,2000,3000,4000)),
     labels=c(200,300,400,500,1000,2000,3000,4000))
axis(side=1, at=log(seq(500,4000,100)), tck=-0.01, labels=NA)
axis(side=2, at=mapply(proportion=c(0,1), axis=c("y","y"), relative.axis.point),
     tck=0, labels=NA)
axis(side=2, at=c(0,0.5,1), tck=-0.02, mgp=c(3,0.4,0))
mtext(side=1, line=1, text="Annual rainfall (mm)")
mtext(side=2, line=1.5, at=0.5, text="Density", las=0)
box()
dev.off()

summary(mixed.site.sub.multi$alpha)
pdf("./Plots/alpha density.pdf", height=2, width=4)
par(mar=c(2.5,3,1,1), ps=8, las=1, tck=-0.02, mgp=c(3,0,0))

alpha.dens<-density(mixed.site.sub.multi$alpha, from=2)
plot(x=NULL, y=NULL, ylim=c(0,0.15),
     xlim=c(2,30), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
polygon(x=c(alpha.dens$x, alpha.dens$x[1]),
        y=c(alpha.dens$y, 0), col="grey90")

axis(side=1, at=seq(5,30,5),
     labels=seq(5,30,5))
axis(side=1, at=2:30, tck=-0.01, labels=NA)
axis(side=2, at=seq(0,0.15,0.05), mgp=c(3,0.4,0))

mtext(side=1, line=1, text="Plot richness")
box()
mtext(side=2, line=2, text="Density", las=0)
dev.off()

pdf("./Plots/age density.pdf", height=2, width=4)
par(mar=c(2.5,3,1,1), ps=8, las=1, tck=-0.02, mgp=c(3,0,0))

age.dens<-density(mixed.site.sub.multi$age, from=0)
plot(x=NULL, y=NULL, ylim=c(0,0.10),
     xlim=c(0,65), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
polygon(age.dens, col="grey90")
axis(side=1, at=seq(0,60,10))
axis(side=1, at=seq(5,55,5), tck=-0.01, labels=NA)
axis(side=2, at=seq(0,0.10,0.025), mgp=c(3,0.4,0))

mtext(side=1, line=1, text="Planting age (years)")
box()
mtext(side=2, line=2, text="Density", las=0)
dev.off()


# PRESENTATION PLOTS ####
#             COR MATRIX AND UNIVARIATE EFFECT SIZE PLOT ####
mixed.site.sub.multi<-mixed.site.sub.multi[complete.cases(mixed.site.sub.multi[,covariates]),]

sapply(as.data.frame(mixed.site.centre[,covariates]), mean)

mixed.site.centre<-sapply(mixed.site.sub.multi[,covariates], function(x){
  (x- mean(x)) / sd(x)
})

mixed.site.centre<-as.data.frame(cbind(mixed.site.sub.multi[,!colnames(mixed.site.sub.multi) %in%
                                                              covariates],
                                       mixed.site.centre))

# read in complete dredge on models
m1.full.pd<-readRDS("./Outputs/Covariate dredge.rds")

scores<-colSums(ifelse(!is.na(m1.full.pd[,colnames(m1.full.pd) %in% covariates]),
                       TRUE, FALSE) * m1.full.pd$weight)
directions<-ifelse(colMeans(m1.full.pd[,colnames(m1.full.pd) %in% covariates], na.rm=TRUE)>0,
                   1, -1)

col.vect<-data.frame(value=c(seq(-1, 1, 0.01)),
                     col=c(colorRampPalette(c(rgb(1,0,0),"white"), bias=2)(100),
                           "white",
                           rev(colorRampPalette(c(rgb(0,0,1),"white"), bias=2)(100))))

#                             PLOT ####
cov.mat<-cor(mixed.site.centre[,names(sort(scores))], 
             use="complete.obs")

cov.mat[lower.tri(cov.mat, diag=TRUE)]<-NA

codes<-data.frame(rcodes=names(sort(scores)),
                  str.codes=c("Nitrogen","Carbon","pH","Phosphorus","TWI","Soil Depth", 
                              "Elevation","Prec Seas",
                              "Silt","Clay","Slope","Solar","Elev Rel","Tmax",
                              "Tmin","Rainfall","Moisture","Sand","Age",
                              "Density"),
                  colour=c(rep("grey50",14),"black","grey50",rep("black",4)),
                  font=c(rep(1, 14), 2, 1, rep(2,4)),
                  pch=c(rep(1,14),16,1, rep(16,4)))

#cov.mat.90<-
#pdf(paste0("./Plots/cov mat and dredge", Sys.date(), ".pdf"), 
#   width=10, height=7, useDingbats=FALSE)

pdf(paste0("./Plots/cov mat and dredge pres", Sys.Date(), ".pdf"), 
    width=7.086, height=5, useDingbats=FALSE)

framemat=rbind(c(0.09,0.675,0.1,0.9),
               c(0.78,0.885,0.1,0.9),
               c(0.91,0.93,0.3,0.7),
               c(0.09,0.885,0.1,0.9))

split.screen(framemat)

#                             MAIN PLOT ####
screen(1)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=6, tck=-0.01, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, xlim=c(-0.025,1.2), ylim=c(0,1.15), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
#axis(side=1)
#axis(side=2)
x.start<-0.2
x.finish<-1.2
y.start<-0
y.finish<-1

x.dist<-(x.finish-x.start)/(dim(cov.mat)[1])
y.dist<-(y.finish-y.start)/(dim(cov.mat)[1])

# get rect coordinates ready to draw
test<-cbind(expand.grid(seq(x.start, x.finish, length.out=20)[-20],
                        seq(y.start, y.finish, length.out=20)[-20]),
            expand.grid(seq(x.start, x.finish, length.out=20)[-1],
                        seq(y.start, y.finish, length.out=20)[-1]))
colnames(test)<-c("xleft","ybottom","xright","ytop")
test[which(is.na(t(cov.mat[-dim(cov.mat)[1],-1]))),]=NA

# get colours of cov mat
cov.mat.vect<-as.vector(t(cov.mat[-dim(cov.mat)[1], -1]))
cov.mat.vect<-cov.mat.vect[!is.na(cov.mat.vect)]

col.scores<-col.vect$col[match(round(cov.mat.vect, 2),
                               round(col.vect$value, 2))]

test.nona<-test[complete.cases(test),]

rect(xleft=par("usr")[1],
     xright=0.1,
     ybottom=par("usr")[3],
     ytop=y.finish + y.dist,
     border=NA, col="grey85")

rect(xleft=test.nona$xleft, 
     xright=test.nona$xright, 
     ybottom=test.nona$ybottom, 
     ytop=test.nona$ytop,
     col=as.character(col.scores))


axis(side=4, at=seq(y.start + 0.5*y.dist, 
                    y.finish + 0.5*y.dist, 
                    length.out=length(rownames(cov.mat))),
     labels=NA, las=1)

mtext(side=4, at=seq(y.start + 0.5*y.dist, 
                     y.finish + 0.5*y.dist, 
                     length.out=length(rownames(cov.mat))),
      text=codes$str.codes, adj=0.5, line=1.825,
      col=as.character(codes$colour), font=codes$font)

axis(side=1, at=seq(x.start +  0.5*x.dist, 
                    x.finish - 0.5*x.dist, 
                    length.out=length(colnames(cov.mat))-1),
     labels=NA, las=3)

par(xpd=NA)
text(x=seq(x.start +  0.5*x.dist, 
           x.finish - 0.5*x.dist, 
           length.out=length(colnames(cov.mat))-1),
     y=-0.025, labels= codes$str.codes[-1], srt=40,
     adj=1,
     col=as.character(codes$colour)[-1], font=codes$font[-1])
segments(x0=par("usr")[2], x1=par("usr")[2],
         y0=par("usr")[3], y1=y.finish+y.dist)

# custom axis lines

segments(x0=par("usr")[1], 
         x1=par("usr")[1],
         y0=par("usr")[3], 
         y1=y.finish + y.dist)

segments(x0=par("usr")[1], 
         x1=1,
         y0=y.finish + y.dist, 
         y1=y.finish + y.dist)

segments(x0=par("usr")[1], 
         x1=par("usr")[2],
         y0=par("usr")[3], 
         y1=par("usr")[3])

segments(x0=par("usr")[2], 
         x1=1.5,
         y0=y.finish + y.dist, 
         y1=y.finish + y.dist)

segments(x0=par("usr")[1], 
         x1=1.5,
         y0=par("usr")[3], 
         y1=par("usr")[3])

par(xpd=FALSE)

axis(side=3, at=c(0,0.2,0.4,0.6,0.8,1), pos=y.finish + y.dist,
     mgp=c(3,0.2,0))
axis(side=3, at=seq(0,1,0.1), pos=y.finish + y.dist,
     labels=NA, tck=-0.005)

mtext(side=3, text="Variable importance",
      line=-0.75)

axis(side=2, at=seq(y.start + 0.5*y.dist,
                    y.finish + 0.5*y.dist, length.out=20),
     labels=NA)

mtext(side=2, at=seq(y.start + 0.5*y.dist, 
                     y.finish + 0.5*y.dist, 
                     length.out=length(rownames(cov.mat))),
      text=codes$str.codes, adj=1, line=0.5,
      col=as.character(codes$colour), font=codes$font)

# ADD ON POINTS

sapply(1:length(scores), function(x){
  
  pos<-seq(y.start + 0.5*y.dist,
           y.finish + 0.5*y.dist, length.out=20)[x]
  
  rect(ybottom=pos-0.25*y.dist,
       ytop=pos+0.25*y.dist,
       xleft=par("usr")[1],
       xright=sort(scores)[x],
       col=ifelse(codes$pch[x]==16, "black", "white"))
  
})

#points(x=sort(scores), 
#       y= seq(y.start + 0.5*y.dist,
#              y.finish + 0.5*y.dist, length.out=20),
#       pch=codes$pch)

close.screen(1)

#                             DIVERSITY CORRELATIONS ####
screen(2)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=6, tck=-0.08, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, xlim=c(-0.025,1.2), ylim=c(0,1.15), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")

axis(side=2, at=seq(y.start + 0.5*y.dist, 
                    y.finish + 0.5*y.dist, 
                    length.out=length(rownames(cov.mat))),
     labels=NA, las=1)

# get rect coordinates ready to draw
test<-cbind(expand.grid(seq(par("usr")[1], par("usr")[2], length.out=5)[1:4],
                        seq(y.start, y.finish+y.dist, length.out=21)[1:20]),
            expand.grid(seq(par("usr")[1], par("usr")[2], length.out=5)[2:5],
                        seq(y.start, y.finish+y.dist, length.out=21)[2:21]))
colnames(test)<-c("xleft","ybottom","xright","ytop")
#test[which(is.na(t(cov.mat[-dim(cov.mat)[1],-1]))),]=NA



cov.mat<-cor(mixed.site.centre[,c("rare.rich","FRv", "FEm", "FDm",
                                  names(sort(scores)))])


cov.mat<-cov.mat[-(1:4),1:4]

cov.colours<-col.vect$col[match(round(t(cov.mat),2),
                                round(col.vect$value, 2))]

rect(xleft=test$xleft, 
     xright=test$xright, 
     ybottom=test$ybottom, 
     ytop=test$ytop,
     col=as.character(cov.colours))

axis(side=1, at=seq(par("usr")[1], par("usr")[2], length.out=9)[c(2,4,6,8)],
     labels=NA)

par(xpd=NA)
text(x=seq(par("usr")[1], par("usr")[2], length.out=9)[c(2,4,6,8)],
     y=-0.025, labels= c("Richness",
                         "FRv", "FEm", "FDm"), srt=40,
     adj=1)
segments(x0=par("usr")[1], x1=par("usr")[1],
         y0=par("usr")[3], y1=y.finish+y.dist)
par(xpd=FALSE)

close.screen(2)

#                             LEGEND ####
screen(3)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=6, tck=-0.3, mgp=c(3,0.5,0))
plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(-1.02,1.02), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")

gradient.rect(xleft=0, xright=1, ybottom=-1.02, ytop=1.02,
              col=as.character(col.vect$col), nslices=length(col.vect$col)[1],
              gradient="y")
axis(side=4, labels=NA)
mtext(side=4, at=c(-1,-0.5,0,0.5,1), text=c(-1,-0.5,0,0.5,1),
      line=0.35)
mtext(side=4, line=1, las=0, text="Correlation coefficient")
close.screen(3)

close.screen(all.screens=TRUE)
dev.off()

screen(4)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=8, tck=-0.3, mgp=c(3,0.5,0))
plot(x=NULL,y=NULL, xlim=c(0,1), ylim=c(0,1.15), xaxs="i", yaxs="i", axes=FALSE,
     xlab="", ylab="")
par(xpd=NA)
rect(xleft=par("usr")[1],
     xright=par("usr")[2],
     ybottom=par("usr")[3],
     ytop=y.finish+y.dist, lwd=1.25)
par(xpd=FALSE)
close.screen(4)
dev.off()



#                             TOP CORRELATIONS IN MATRIX ####
colnames(mixed.site.sub.multi         )
write.csv(cor(mixed.site.sub.multi[,c("rare.rich", "FRv","FEm","FDm",
                            "log.density", "log.age", "sand",
                            "log.aridity.index", "tmin.annual")],
    use="complete.obs")[-c(1:4),1:4],
    "./Outputs/correlation matrix for presentation.csv")

# ####
# MODELLING - PNAS REVISIONS ####
#             LOOK FOR NEARBY FIRES ####
head(mixed.site.sub)
rain<-raster("/home/Storage HDD/University files/Shape files/prec.annual.tif")
# Import NASA MODIS data from 2000 to 2010
modis<-readShapeSpatial("./Data/fire_archive_M6_19363_edit")

# Each point at minimum represents 1km2 of fire, so let's conver the shape file
# into a 1km2 raster, using one of our 1km2 climate rasters.
modis.rast<-rasterize(modis, y=rain, field="BRIGHTNESS")

# now we extract the points at each of our site coordinates to see if a fire was
# registered between 2001 and 2010.
sub.data<-mixed.site.sub[rownames(mixed.site.sub) %in% rownames(model.data),
                         ]

site.coords<-mixed.site.sub[rownames(mixed.site.sub) %in% rownames(model.data),
                            c("long","lat")]
coordinates(site.coords)<-c("long","lat")

site.fire<-extract(x=modis.rast, y=site.coords)

temp<-as.matrix(table(is.na(site.fire), sub.data$best.cluster))
x<-temp[,1]
apply(temp, 2, function(x){
  x / sum(x)
})

pdf("./Plots/fire plot.pdf", height=3, width=5, useDingbats=FALSE)
split.screen(rbind(c(0.1,0.65,0.125,0.95),
                   c(0.775,0.99,0.125,0.95)))

screen(1)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=8, tck=-0.01, mgp=c(3,0.5,0))
with(mixed.site.sub[is.na(site.fire),], plot(log(biomass.area) ~ log.age,
                                             col=rgb(0.5,0.5,0.5,0.5), pch=16,
                                             xaxt="n"))
with(mixed.site.sub[!is.na(site.fire),], points(log(biomass.area) ~ log.age, 
                                                pch=16, col="red"))
mtext(side=2, text=expression("Plot biomass (ln(kg ha"^-1*"))"), las=0, line=1.25)

axis(side=1, at=log(c(2,5,10,20,50)), labels=c(2,5,10,20,50), mgp=c(3,0,0))
axis(side=1, at=log(c(seq(1,20,1), seq(20,100,10))), 
     labels=NA, tck=-0.005)
mtext(side=1, text="Planting age (years)", line=0.8)

text(x=relative.axis.point(0.05, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(a)", font=2)
close.screen(1)

screen(2)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), las=1, ps=8, tck=-0.025, mgp=c(3,0.5,0))
boxplot(resid(env.model) ~ as.factor(is.na(site.fire)), xaxt="n", cex=0.8)
axis(side=1, at=c(1,2), labels=c("Fire", 'No fire'), mgp=c(3,0.1,0))

mtext(side=2, text=expression("Residual plot biomass (ln(kg ha"^-1*"))"), las=0, line=1.75)
text(x=relative.axis.point(0.1, "x"),
     y=relative.axis.point(0.95, "y"),
     labels="(b)", font=2)
close.screen(all.screens=TRUE)
dev.off()



#             BIOMASS VARIATION MODELS ####

# Compare variation in residual biomass WITHIN plantings ####

# Let's get the variance in residual biomass from each planting (min of 2 plots),
# and plot that against richness, then we can do it in each cluster.

model.data<-droplevels(rich.model$data)
model.data$resids<-resid(rich.model)

length(levels(model.data$planting))
planting.stats<-do.call("rbind", 
                        lapply(split(model.data, f=model.data$planting), 
                               function(x){
                                 
                                 ifelse(dim(x)[1]<2, 1, 0)}))
table(planting.stats)

planting.stats<-do.call("rbind", 
                        lapply(split(model.data, f=model.data$planting), 
                                     function(x){
  
  if(dim(x)[1]<2){return(NULL)}
  
  return(data.frame(resid.var=var(x$resids),
                    resid.sd=sd(x$resids),
                    rich.mean=mean(x$rare.rich),
                    rich.sd=sd(x$rare.rich),
                    plot.n=dim(x)[1]))
                                       
}))

raw.mean.rich<-do.call("rbind", 
                       lapply(split(mixed.site.sub.multi, 
                                    f=mixed.site.sub.multi$planting), 
                              function(x){
                                
                                if(dim(x)[1]<2){return(NULL)}
                                
                                return(data.frame(raw.rich.mean=mean(x$rare.rich),
                                                  raw.rich.se=sd(x$rare.rich) / sqrt(dim(x)[1]),
                                                  planting=x$planting[1]))
                              }))

planting.stats$planting<-rownames(planting.stats)
planting.stats<-merge(planting.stats, model.data[!duplicated(model.data$planting),
                                                 c("planting", "best.cluster", "IBRAsub")],
                      all.x=TRUE, all.y=FALSE)
planting.stats<-merge(planting.stats, raw.mean.rich,
                      all.x=TRUE, all.y=FALSE)

dev.off()
plot(log(planting.stats$resid.var) ~ planting.stats$plot.n)
plot(log(planting.stats$resid.var) ~ planting.stats$raw.rich.mean)

cor(planting.stats[,c("plot.n", "rich.mean")])


temp<-lme(resid.var ~ raw.rich.mean, random=~1|IBRAsub, data=planting.stats)
summary(temp)
plot(temp)

temp<-lme(log(resid.var) ~ raw.rich.mean, random=~1|IBRAsub, data=planting.stats)
summary(temp)
plot(temp)

planting.stats$best.cluster<-as.factor(planting.stats$best.cluster)
temp<-lme(log(resid.var) ~ raw.rich.mean*best.cluster, random=~1|IBRAsub, data=planting.stats)
var.cluster.model<-glht(temp,
                        linfct = c("raw.rich.mean = 0",
                                   paste("raw.rich.mean + ", paste("raw.rich.mean:best.cluster", 2:8, sep=""), "= 0")))
write.csv(round(do.call("cbind", summary(var.cluster.model)$test[3:6]), 3),
          "./Outputs/within planting variance coefs.csv")



cluster.df<-read.csv("./Data/best cluster summary.csv")
pdf(paste0("./Plots/biomass variance cluster plot lme", Sys.Date(), ".pdf"),
    height=3.15, width=6, useDingbats=FALSE)

par(mfrow=c(2,4), mar=c(0,0,0,0), oma=c(3,4,2,2), ps=9, las=1, mgp=c(3,0.5,0),
    tck=-0.01)

cluster.var.models<-lapply(1:8, function(x){

  temp.data<-planting.stats[planting.stats$best.cluster==x,]

      rich.vals<-seq(min(temp.data$raw.rich.mean),
                 max(temp.data$raw.rich.mean),
                 length=200)
  temp.preds<-predict(temp, newdata=data.frame(raw.rich.mean=rich.vals,
                                               best.cluster=factor(rep(x, 200),
                                                                   levels(temp.data$best.cluster))),
          se.fit=TRUE, level=0)
  
  with(temp.data, plot(log(resid.var) ~ raw.rich.mean, pch=16,
                       col=rgb(cluster.df$red[x],
                               cluster.df$green[x],
                               cluster.df$blue[x]),
                       xlim=c(1,16.75), ylim=c(-12, 2),
                       axes=FALSE))
  
  if(x %in% c(5:8)){
  axis(side=1, at=c(1,5,10,15), labels=c(1,5,10,15), mgp=c(3,0,0))
  mtext(side=1, line=1.25, at=par("usr")[2], text=ifelse(x==6, "Mean rarefied planting richness", ""),
        cex=0.8)
  } else {axis(side=1, at=c(1,5,10,15), labels=NA)}
  
  if(x %in% c(1,5)){
    axis(side=2, at=log(c(1,0.1,0.01,0.001,0.0001,0.00001)),
         labels=c(1,0.1,0.01,0.001,0.0001,0.00001))
    mtext(side=2, line=2.75, at=par("usr")[4], las=0,
          text=ifelse(x==5, expression("Variaion in residual plot biomass (ln(kg ha"^-1*"))"), ""),
          cex=0.8)
  } else {axis(side=2, at=log(c(1,0.1,0.01,0.001,0.0001,0.00001)),
               labels=NA)}
  
  polygon(x=c(rich.vals, rev(rich.vals)),
          y=c(temp.preds$fit + 1.96*temp.preds$se.fit,
              rev(temp.preds$fit - 1.96*temp.preds$se.fit)),
          col=rgb(0.3,0.3,0.3,0.4), border=NA)
  points(temp.preds$fit ~ rich.vals, type="l", lwd=1.5)
  box()
  text(x=relative.axis.point(0.05, "x"),
       y=relative.axis.point(0.925, "y"),
       labels=cluster.df$name[x], adj=0)

  })
dev.off()
lapply(cluster.var.models, function(x){summary(x)$tTable})

# Compare mean biomass across plantings ####

# Some of our plots suggest that higher diversity might reduce the variability
# in productivity. There's a few ways we could test this: The simplest is to
# include a variance structure in our models and see if it soaks up variation in
# biomass.

# Because variance structures only work within groups, we can only get it to
# work for plots within plantings. Our alternative is to get mean biomass &
# richness for each planting and model it without random effects (gls)

model.data<-mixed.site.sub.multi

planting.stats<-do.call("rbind", 
                        lapply(split(model.data, f=model.data$planting), 
                               function(x){
                                 
                                 if(dim(x)[1]<2){return(NULL)}
                                 
                                 return(data.frame(mass.mean=mean(x$biomass.area),
                                                   mass.sd=sd(x$biomass.area),
                                                   density.mean=log(mean(x$density)),
                                                   rich.mean=mean(x$rare.rich),
                                                   rich.sd=sd(x$rare.rich),
                                                   plot.n=dim(x)[1]))
                                 
                               }))

# now we'll run our richness model using the new data, and start introducing 
# variance structures

# First let's merge across our environmental stats
planting.stats$planting<-rownames(planting.stats)
planting.stats<-merge(planting.stats, model.data[!duplicated(model.data$planting),
                                                 c("planting", "best.cluster", "IBRAsub",
                                                   "log.aridity.index", "sand", "tmin.annual",
                                                   "log.age")],
                      all.x=TRUE, all.y=FALSE, by.x="planting", by.y="planting")

# now run model without any variance structure
planting.stats<-planting.stats[complete.cases(planting.stats),]

base.model<-gls(log(mass.mean) ~ log.age + density.mean + log.aridity.index +
                  sand + tmin.annual + log.age:density.mean + density.mean:log.aridity.index +
                  log.age:tmin.annual + rich.mean,
                data=planting.stats)
summary(base.model)

# now add in different variance weightings for each cluster
cluster.model<-gls(log(mass.mean) ~ log.age + density.mean + log.aridity.index +
                  sand + tmin.annual + log.age:density.mean + density.mean:log.aridity.index +
                  log.age:tmin.annual + rich.mean, weights=varIdent(form=~1|best.cluster),
                  data=planting.stats)
summary(cluster.model)
anova(base.model, cluster.model)
# That does improve the model. now sub that for just variance changing over levels of richness

richvar.model<-gls(log(mass.mean) ~ log.age + density.mean + log.aridity.index +
                     sand + tmin.annual + log.age:density.mean + density.mean:log.aridity.index +
                     log.age:tmin.annual + rich.mean, weights=varFixed(~rich.mean),
                   data=planting.stats)
summary(richvar.model)
anova(base.model, richvar.model)

# that makes the model worse. Let's now try both together

richclust.model<-gls(log(mass.mean) ~ log.age + density.mean + log.aridity.index +
                     sand + tmin.annual + log.age:density.mean + density.mean:log.aridity.index +
                     log.age:tmin.annual + rich.mean, 
                   weights=varPower(form=~rich.mean|best.cluster),
                   data=planting.stats)
summary(richclust.model)
anova(base.model, richclust.model)
anova(cluster.model, richclust.model)

var.models<-list(base.model, cluster.model, richvar.model, richclust.model)

sapply(var.models, AIC)
anova(var.models[[3]], var.models[[1]])


write.csv(as.data.frame(anova(base.model, cluster.model, richvar.model, richclust.model)),
          "./Outputs/variance models.csv")


planting.stats$density.mean
# Let's use the same modelling process we used for the previous models, starting
# with the environment and richness model.

var.model<-update(rich.model, .~., weights=varFixed(~rare.rich))
summary(var.model)
anova(var.model, rich.model)
# variance structure doesn't improve model

# now let's try for each cluster.

# start with just different variance per cluster
var.clust.model<-update(richness.only.cluster.model.nomono,
                        .~., weights=varIdent(form=~1|best.cluster))
anova(richness.only.cluster.model.nomono, var.clust.model)

# that does improve model, suggest inherently different variance in residual biomass
# in different clusters - not surprising given residual plots in Fig 2.

# does this variance decrease with richness in different clusters

summary(richness.only.cluster.model.nomono)
var.clust.model2<-update(richness.only.cluster.model.nomono,
                        .~., weights=varExp(form=~rare.rich|best.cluster),
                        control=lmeControl(msMaxIter=100, maxIter=100, returnObject=TRUE))

var.clust.model3<-update(richness.only.cluster.model.nomono,
                         .~., weights=varPower(form=~rare.rich|best.cluster),
                         control=lmeControl(msMaxIter=100, maxIter=100, returnObject=TRUE))
anova(var.clust.model2, var.clust.model3)

anova(richness.only.cluster.model.nomono, var.clust.model)
# variance does change with richness differently across clusters.

# how about our functional diversity indices?
FRv.var.model<-update(fun.div.env, .~., weights=varFixed(~FRv))
FEm.var.model<-update(fun.div.env, .~., weights=varFixed(~FEm))
FDm.var.model<-update(fun.div.env, .~., weights=varFixed(~FDm))

anova(FRv.var.model, fun.div.env)
anova(FEm.var.model, fun.div.env)
anova(FDm.var.model, fun.div.env)

# not across Australia, but let's try the cluster models
FRv.var.clust.model<-update(fun.div.cluster,
                        .~., weights=varIdent(form=~FRv|best.cluster))
FEm.var.clust.model<-update(fun.div.cluster,
                            .~., weights=varIdent(form=~FEm|best.cluster))
FDm.var.clust.model<-update(fun.div.cluster,
                            .~., weights=varIdent(form=~FDm|best.cluster))

anova(FRv.var.clust.model, fun.div.cluster)
anova(FEm.var.clust.model, fun.div.cluster)
anova(FDm.var.clust.model, fun.div.cluster)


#             ANALYSIS OF DEAD STEMS ####
plot(x=NULL, y=NULL, ylim=log(c(0.01,3000)), xlim=c(0.75,2), type="n", axes=FALSE)
box()
abline(h=log(c(0.1,1,10,100,1000)), col="grey50")

boxplot(log(mixed.plant.sub$tot.ag.mass[mixed.plant.sub$dead==1]), 
        yaxt="n", add=TRUE)

points(log(mixed.plant.sub$tot.ag.mass[mixed.plant.sub$dead==1]) ~
         jitter(rep(1.5, sum(mixed.plant.sub$dead)), amount=0.2),
       pch=16, col=rgb(0.3,0.3,0.3,0.3))

limits<-cbind(c(exp(par("usr")[3]),0.1,1,10,100,1000),
              c(0.1,1,10,100,1000,exp(par("usr")[4])))
props<-apply(limits, 1, function(x){
  temp<-sum(mixed.plant.sub$tot.ag.mass[mixed.plant.sub$dead==1] > x[1] &
            mixed.plant.sub$tot.ag.mass[mixed.plant.sub$dead==1] <= x[2])
  temp / sum(mixed.plant.sub$dead==1)
})

text(x=1.85, y=log(limits[,1]) + 0.5*(log(limits[,2])-log(limits[1])),
     labels=round(props, 3))

axis(side=2, at=log(c(0.001,0.01,0.1,1,10,100,1000,10000)),
     labels=c(0.001,0.01,0.1,1,10,100,1000,10000), las=1)
axis(side=2, at=log(c(seq(0.01,0.1,0.01), seq(0.1,1,0.1), seq(1,10,1),
                      seq(10,100,10), seq(100,1000,100), seq(1000,10000,1000))),
     labels=NA, tck=-0.025)


# ####
