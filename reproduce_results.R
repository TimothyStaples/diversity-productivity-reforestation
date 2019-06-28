# #################################################################### ####
# Title: Productivity does not correlate with species and functional   ####
#        diversity in Australian reforestation plantings across a wide ####
#        climate gradient                                              #### 
# Author: Timothy L Staples                                            ####
# Details: Minimal code to reproduce results and figures               ####  
# #################################################################### ####
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

# IMPORT DATA TO REPRODUCE RESULTS ####

mixed.site.sub.multi <- read.csv("./Data/plot.level.modelling.data.csv")

model.random.effects <- "IBRAsub/planting"
best.covariates <- c("log.age", "log.density", "log.aridity.index", "tmin.annual",
                     "sand")
covariates <- read.csv("./Data/all.covariates.csv",
                       stringsAsFactors = FALSE)[,2]

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

# MODELLING ####

#             INTERCEPT ONLY MODEL ####

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
mass.ratio.vars<-colnames(mixed.site.sub.multi)[grepl("cwm", colnames(mixed.site.sub.multi))]

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
mass.ratio.cluster.hyptest

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
# PLOTTING ####
#             FIGURE 1 ####
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


pdf(paste0("./Plots/Figure 1", ".pdf"), 
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



#             FIGURE 2 ####

# Figure 2 cannot be reproduced in its entirety. We cannot provide the coordinates
# for the plantings, so the map in sub-panel (e) has been excluded from this version.

# read in details of clusters
cluster.df<-read.csv("./Data/best cluster summary.csv", header=TRUE)
cluster.df<-as.data.frame(sapply(cluster.df, function(x){gsub("\\+AC0", "", x)}))
cluster.df[,c("red","green","blue")]<-sapply(cluster.df[,c("red","green","blue")],
                                             function(x){as.numeric(as.character(x))})

plot(NULL, lwd=0.5, xlim=c(115,155), ylim=c(-45,-9), axes=FALSE)
map.usr<-par("usr")

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

pdf(paste0("./Plots/Figure 2", ".pdf"), width=9, height=10, 
    useDingbats=FALSE)
par(mar=c(0,0,0,0), oma=c(0,0,0,0), ps=10, las=1, tck=-0.0075, mgp=c(3,0.4,0))
split.screen(framemat)

#                             MAP ####
screen(1)
par(ps=10)
plot(x=NULL,y=NULL, xlim=c(107.8123,162.1877), ylim=c(-46.44,-7.56), axes=FALSE)
# Where raster layer would have been plotted
# plot(log.prec, col=colorRampPalette(c("white",rgb(0.9,0.9,1,1), rgb(0.15,0.15,1,1)))(200), 
#      box=FALSE, add=TRUE, legend=FALSE, maxpixels=1000000, interpolate=TRUE)
mappar<-par("usr")
close.screen(1)

screen(2)
par(ps=10)
plot(x=NULL,y=NULL, xlim=c(107.8123,162.1877), ylim=c(-46.44,-7.56), axes=FALSE)
# Where map shape file would have been plotted
# plot(ausmap, lwd=1, axes=FALSE, border="grey40", add=TRUE)

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
  axis(side=1, at=(c(1,5,10,15)-mean(mixed.site.sub.multi$rare.rich))/sd(mixed.site.sub.multi$rare.rich), 
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

#             FIGURE 3 ####
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

pdf("./Plots/Figure 3.pdf", height=5, width=3.5)
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
