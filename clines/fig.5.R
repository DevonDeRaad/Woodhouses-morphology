

install.packages("hzar")
## Load the package
library(hzar)

if(require(doMC)){
  registerDoMC();
}else if(require(doParallel)){
  registerDoParallel()
}else{
  registerDoSeq()
}

## A typical chain length
chainLength=1e5                       

## Make each model run off a separate seed
## mainSeed=
##   list(
##        A=c(978,544,99,596,528,124),
##        B=c(544,99,596,528,124,978),
##        C=c(99,596,528,124,978,544))

## Increment skip by 15 or so for every trait
rotateModelSeeds <- function(fiComp.1)
  hzar.multiFitRequest(fiComp.1,
                       rotateSeed=TRUE, skip=500,
                       baseChannel=NULL, each=1,
                       baseSeed=c(596,528,124,978,544,99))
## ## Save all plots in a pdf file?
#pdf(width=8, height=8, file="male.comp1.pdf")

#read.csv("~/Dropbox/scrub poster/manuscript/hzar.scrubbers/morph.csv")->oldmorph
#malemorph<-malescrub[,c(1,3,16)]
#oldmorph<-oldmorph[,1:3]
#newmorph<-merge(oldmorph,malemorph, by = "id")
#write.csv(newmorph, file = "~/Dropbox/scrub poster/manuscript/hzar.scrubbers/malemorph.csv")

#bring in data
SiteLoc<-read.csv("~/Dropbox/scrub poster/manuscript/hzar.scrubbers/siteloc.csv")
#for males
malemorph<-read.csv("~/Dropbox/scrub poster/manuscript/hzar.scrubbers/malemorph.csv")
#for females
femalemorph<-read.csv("~/Dropbox/scrub poster/manuscript/hzar.scrubbers/femalemorph.csv")
#blue plumage of the back for all individuals
S1blue<-read.csv("~/Dropbox/scrub poster/manuscript/hzar.scrubbers/s1.blue.csv")
head(S1blue)

## Trait Analysis of Comp.1

## Load example Quantitative trait data from the package
## Note that this is the individual data, with labels
## identifying the locality of each sample.

## Print a summary of the trait data

## Load example locality data, matching each localility to
## a site ID and a transect distance.

## Print the locality information

## Blank out space in memory to hold morphological analysis
if(length(apropos("^mkn$",ignore.case=FALSE)) == 0 ||
   !is.list(mkn) ) mkn <- list()
## We are doing just the one quantitative trait, but it is
## good to stay organized.
mkn$Comp.1 <- list();
## Space to hold the observed data
mkn$Comp.1$obs <- list();
## Space to hold the models to fit
mkn$Comp.1$models <- list();
## Space to hold the compiled fit requests
mkn$Comp.1$fitRs <- list();
## Space to hold the output data chains
mkn$Comp.1$runs <- list();
## Space to hold the analysed data
mkn$Comp.1$analysis <- list();


## Beard Length Trait from Brumfield et al 2001
mkn$Comp.1$obs <-
  hzar.doNormalData1DRaw(hzar.mapSiteDist(SiteLoc$Site,
                                          SiteLoc$Distance),
                         malemorph$Site,
                         malemorph$Comp.1)

## Look at a graph of the observed data
hzar.plot.obsData(mkn$Comp.1$obs);

## Make a helper function
mkn.loadComp.1model <- function(scaling,tails,
                             id=paste(scaling,tails,sep=".")){
  mkn$Comp.1$models[[id]] <<-
    hzar.makeCline1DNormal(mkn$Comp.1$obs, tails)
  ## As there is no quick option for "fixed" scaling, and the
  ## combined sites "A", "B" and "J","K" have a fair number of samples (> 20),
  ## fix the mean and variance of the left and right sides of
  ## the cline to the values observed by the combined sites.
  if (all(regexpr("fixed",scaling,ignore.case=TRUE) == 1 )){
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$muL <<- TRUE
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$muR <<- TRUE
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$varL <<- TRUE
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$varR <<- TRUE
  }
  ## Helper function to work around low sample size
  getCombo <- function(id1,id2,colName,frame=mkn$Comp.1$obs$frame)
    (frame[id1,colName]*frame[id1,"nEff"]+
       frame[id2,colName]*frame[id2,"nEff"])/
    (frame[id1,"nEff"]+frame[id2,"nEff"])
  
  ## Site A, B is the "left" side of the cline, so pull the
  ## fixed values from there.
  hzar:::meta.init(mkn$Comp.1$models[[id]])$muL <<-
    getCombo("A","B","mu")
  hzar:::meta.init(mkn$Comp.1$models[[id]])$varL <<-
    getCombo("A","B","var")
  ## Site J, K is the "right" side of the cline, so pull the
  ## fixed values from there.
  hzar:::meta.init(mkn$Comp.1$models[[id]])$muR <<-
    getCombo("I","J","mu")
  hzar:::meta.init(mkn$Comp.1$models[[id]])$varR <<-
    getCombo("I","J","var")
  ## Make a better estimate of varH using site D, E
  hzar:::meta.init(mkn$Comp.1$models[[id]])$varH <<-
    getCombo("E","F","var")-
    (getCombo("I","J","var")+getCombo("A","B","var"))/2
  
}
## mkn.loadComp.1model("fixed","none","modelI");
## mkn.loadComp.1model("free" ,"none","modelII");
## mkn.loadComp.1model("free" ,"both","modelIII");
mkn.loadComp.1model("fixed","none"  ,"fixN");
mkn.loadComp.1model("fixed","left"  ,"fixL");
mkn.loadComp.1model("fixed","right" ,"fixR");
mkn.loadComp.1model("fixed","mirror","fixM");
mkn.loadComp.1model("fixed","both"  ,"fixB");
mkn.loadComp.1model("free" ,"none"  ,"optN");
mkn.loadComp.1model("free" ,"left"  ,"opComp.1");
mkn.loadComp.1model("free" ,"right" ,"optR");
mkn.loadComp.1model("free" ,"mirror","optM");
mkn.loadComp.1model("free" ,"both"  ,"optB");


## Check the default settings
##print(mkn$Comp.1$models)

## Modify all models to focus on the region where the observed
## data were collected.
## Observations were between -10 and 130000 km.
mkn$Comp.1$models <- sapply(mkn$Comp.1$models,
                         hzar.model.addBoxReq,
                         -10 , 1100,
                         simplify=FALSE)

## Due to the large number of free variables, it is prudent to
## reduce the tune setting of optB from 1.5 to 1.1
hzar:::meta.tune(mkn$Comp.1$models$optB)<-1.1


## Check the updated settings
##print(mkn$Comp.1$models)

## Compile each of the models to prepare for fitting
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
mkn$Comp.1$fitRs$init <- sapply(mkn$Comp.1$models,
                             hzar.first.fitRequest.gC,
                             obsData=mkn$Comp.1$obs,
                             verbose=TRUE,
                             simplify=FALSE)

## Update the settings for the fitter if desired.
mkn$Comp.1$fitRs$init <- sapply(mkn$Comp.1$fitRs$init,
                             function(mdl) {
                               mdl$mcmcParam$chainLength <-
                                 chainLength; #1e5
                               mdl$mcmcParam$burnin <-
                                 chainLength %/% 10; #1e4
                               mdl },
                             simplify=FALSE)

mkn$Comp.1$fitRs$init <-
  rotateModelSeeds(mkn$Comp.1$fitRs$init)


## Check fit request settings
##print(mkn$Comp.1$fitRs$init)

## Do just one run of the models for an initial chain
mkn$Comp.1$runs$init <-
  hzar.doFit.multi(mkn$Comp.1$fitRs$init)
names(mkn$Comp.1$runs$init) <- names(mkn$Comp.1$fitRs$init)



## Compile a new set of fit requests using the initial chains 
mkn$Comp.1$fitRs$chains <-
  lapply(mkn$Comp.1$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
mkn$Comp.1$fitRs$chains <-
  hzar.multiFitRequest(mkn$Comp.1$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit


## Go ahead and run a chain of 3 runs for every fit request
mkn$Comp.1$runs$chains <-  hzar.doChain.multi(mkn$Comp.1$fitRs$chains,
                                           doPar=TRUE,
                                           inOrder=FALSE,
                                           count=3)

## Did fixN converge?
summary(do.call(mcmc.list,
                lapply(mkn$Comp.1$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did optN converge?
summary(do.call(mcmc.list,
                lapply(mkn$Comp.1$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did optB converge?
summary(do.call(mcmc.list,
                lapply(mkn$Comp.1$runs$chains[28:30],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Start aggregation of data for analysis

## Clear out a spot to collect the data for analysis (note that
## there is currenComp.1y no "null model" to compare against).
mkn$Comp.1$analysis$initDGs <- list(
)

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
mkn$Comp.1$analysis$initDGs <-
  c( mkn$Comp.1$analysis$initDGs,
     sapply(mkn$Comp.1$runs$init,
            hzar.dataGroup.add,
            simplify=FALSE))

## Create a hzar.obsDataGroup object from the four hzar.dataGroup
## just created, copying the naming scheme (modelI, modelII,
## modelIII).
mkn$Comp.1$analysis$oDG <-
  hzar.make.obsDataGroup(mkn$Comp.1$analysis$initDGs)
mkn$Comp.1$analysis$oDG <-
  hzar.copyModelLabels(mkn$Comp.1$analysis$initDGs,
                       mkn$Comp.1$analysis$oDG)

## Convert all 90 runs to hzar.dataGroup objects, adding them to
## the hzar.obsDataGroup object.
mkn$Comp.1$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(mkn$Comp.1$runs$chains,
                                hzar.dataGroup.add),
                         mkn$Comp.1$analysis$oDG);

## Check to make sure that there are only ten hzar.dataGroup
## objects in the hzar.obsDataGroup object.
print(summary(mkn$Comp.1$analysis$oDG$data.groups))


## Look at the variation in parameters for cline models
oDGkey <- which(!(names(mkn$Comp.1$analysis$oDG$data.groups) %in%
                    "nullModel"));
print(hzar.getLLCutParam(mkn$Comp.1$analysis$oDG$data.groups[oDGkey ],
                         c("center","width")));

## Compare the 3 cline models to the null model graphically
##hzar.plot.cline(mkn$Comp.1$analysis$oDG);

## Do model selection based on the AICc scores
print(mkn$Comp.1$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(mkn$Comp.1$analysis$oDG));

## Print out the model with the minimum AICc score
print(mkn$Comp.1$analysis$model.name <-
        rownames(mkn$Comp.1$analysis$AICcTable
        )[[ which.min(mkn$Comp.1$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
mkn$Comp.1$analysis$model.selected <-
  mkn$Comp.1$analysis$oDG$data.groups[[mkn$Comp.1$analysis$model.name]]

## Plot the maximum likelihood cline for the selected model
##hzar.plot.cline(mkn$Comp.1$analysis$model.selected);

print(hzar.get.ML.cline(mkn$Comp.1$analysis$model.selected))

## Plot the 95% credible cline region for the selected model
##hzar.plot.fzCline(mkn$Comp.1$analysis$model.selected);

par(mar=c(5,5,4,2)+.1)
par(ps=22, cex=1, cex.main=1)
male.pc1<-hzar.plot.fzCline(mkn$Comp.1$analysis$model.selected,
                  pch=19,
                  xlab="Distance (km)",
                  ylab="male PC1",
                  ylim=c(-3,3.5),xlim=c(0,1200))

#export as 8"x8" PDF
#dev.off()
## End Comp.1

## We are doing just the one quantitative trait, but it is
## good to stay organized.
mkn$Comp.1 <- list();
## Space to hold the observed data
mkn$Comp.1$obs <- list();
## Space to hold the models to fit
mkn$Comp.1$models <- list();
## Space to hold the compiled fit requests
mkn$Comp.1$fitRs <- list();
## Space to hold the output data chains
mkn$Comp.1$runs <- list();
## Space to hold the analysed data
mkn$Comp.1$analysis <- list();


## Beard Length Trait from Brumfield et al 2001
mkn$Comp.1$obs <-
  hzar.doNormalData1DRaw(hzar.mapSiteDist(SiteLoc$Site,
                                          SiteLoc$Distance),
                         femalemorph$Site,
                         femalemorph$Comp.1)

## Look at a graph of the observed data
hzar.plot.obsData(mkn$Comp.1$obs);

## Make a helper function
mkn.loadComp.1model <- function(scaling,tails,
                                id=paste(scaling,tails,sep=".")){
  mkn$Comp.1$models[[id]] <<-
    hzar.makeCline1DNormal(mkn$Comp.1$obs, tails)
  ## As there is no quick option for "fixed" scaling, and the
  ## combined sites "A", "B" and "J","K" have a fair number of samples (> 20),
  ## fix the mean and variance of the left and right sides of
  ## the cline to the values observed by the combined sites.
  if (all(regexpr("fixed",scaling,ignore.case=TRUE) == 1 )){
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$muL <<- TRUE
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$muR <<- TRUE
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$varL <<- TRUE
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$varR <<- TRUE
  }
  ## Helper function to work around low sample size
  getCombo <- function(id1,id2,colName,frame=mkn$Comp.1$obs$frame)
    (frame[id1,colName]*frame[id1,"nEff"]+
       frame[id2,colName]*frame[id2,"nEff"])/
    (frame[id1,"nEff"]+frame[id2,"nEff"])
  
  ## Site A, B is the "left" side of the cline, so pull the
  ## fixed values from there.
  hzar:::meta.init(mkn$Comp.1$models[[id]])$muL <<-
    getCombo("A","B","mu")
  hzar:::meta.init(mkn$Comp.1$models[[id]])$varL <<-
    getCombo("A","B","var")
  ## Site J, K is the "right" side of the cline, so pull the
  ## fixed values from there.
  hzar:::meta.init(mkn$Comp.1$models[[id]])$muR <<-
    getCombo("I","J","mu")
  hzar:::meta.init(mkn$Comp.1$models[[id]])$varR <<-
    getCombo("I","J","var")
  ## Make a better estimate of varH using site D, E
  hzar:::meta.init(mkn$Comp.1$models[[id]])$varH <<-
    getCombo("E","F","var")-
    (getCombo("I","J","var")+getCombo("A","B","var"))/2
  
}
## mkn.loadComp.1model("fixed","none","modelI");
## mkn.loadComp.1model("free" ,"none","modelII");
## mkn.loadComp.1model("free" ,"both","modelIII");
mkn.loadComp.1model("fixed","none"  ,"fixN");
mkn.loadComp.1model("fixed","left"  ,"fixL");
mkn.loadComp.1model("fixed","right" ,"fixR");
mkn.loadComp.1model("fixed","mirror","fixM");
mkn.loadComp.1model("fixed","both"  ,"fixB");
mkn.loadComp.1model("free" ,"none"  ,"optN");
mkn.loadComp.1model("free" ,"left"  ,"opComp.1");
mkn.loadComp.1model("free" ,"right" ,"optR");
mkn.loadComp.1model("free" ,"mirror","optM");
mkn.loadComp.1model("free" ,"both"  ,"optB");


## Check the default settings
##print(mkn$Comp.1$models)

## Modify all models to focus on the region where the observed
## data were collected.
## Observations were between -10 and 130000 km.
mkn$Comp.1$models <- sapply(mkn$Comp.1$models,
                            hzar.model.addBoxReq,
                            -10 , 1100,
                            simplify=FALSE)

## Due to the large number of free variables, it is prudent to
## reduce the tune setting of optB from 1.5 to 1.1
hzar:::meta.tune(mkn$Comp.1$models$optB)<-1.1


## Check the updated settings
##print(mkn$Comp.1$models)

## Compile each of the models to prepare for fitting
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
mkn$Comp.1$fitRs$init <- sapply(mkn$Comp.1$models,
                                hzar.first.fitRequest.gC,
                                obsData=mkn$Comp.1$obs,
                                verbose=TRUE,
                                simplify=FALSE)

## Update the settings for the fitter if desired.
mkn$Comp.1$fitRs$init <- sapply(mkn$Comp.1$fitRs$init,
                                function(mdl) {
                                  mdl$mcmcParam$chainLength <-
                                    chainLength; #1e5
                                  mdl$mcmcParam$burnin <-
                                    chainLength %/% 10; #1e4
                                  mdl },
                                simplify=FALSE)

mkn$Comp.1$fitRs$init <-
  rotateModelSeeds(mkn$Comp.1$fitRs$init)


## Check fit request settings
##print(mkn$Comp.1$fitRs$init)

## Do just one run of the models for an initial chain
mkn$Comp.1$runs$init <-
  hzar.doFit.multi(mkn$Comp.1$fitRs$init)
names(mkn$Comp.1$runs$init) <- names(mkn$Comp.1$fitRs$init)



## Compile a new set of fit requests using the initial chains 
mkn$Comp.1$fitRs$chains <-
  lapply(mkn$Comp.1$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
mkn$Comp.1$fitRs$chains <-
  hzar.multiFitRequest(mkn$Comp.1$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit


## Go ahead and run a chain of 3 runs for every fit request
mkn$Comp.1$runs$chains <-  hzar.doChain.multi(mkn$Comp.1$fitRs$chains,
                                              doPar=TRUE,
                                              inOrder=FALSE,
                                              count=3)

## Did fixN converge?
summary(do.call(mcmc.list,
                lapply(mkn$Comp.1$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did optN converge?
summary(do.call(mcmc.list,
                lapply(mkn$Comp.1$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did optB converge?
summary(do.call(mcmc.list,
                lapply(mkn$Comp.1$runs$chains[28:30],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Start aggregation of data for analysis

## Clear out a spot to collect the data for analysis (note that
## there is currenComp.1y no "null model" to compare against).
mkn$Comp.1$analysis$initDGs <- list(
)

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
mkn$Comp.1$analysis$initDGs <-
  c( mkn$Comp.1$analysis$initDGs,
     sapply(mkn$Comp.1$runs$init,
            hzar.dataGroup.add,
            simplify=FALSE))

## Create a hzar.obsDataGroup object from the four hzar.dataGroup
## just created, copying the naming scheme (modelI, modelII,
## modelIII).
mkn$Comp.1$analysis$oDG <-
  hzar.make.obsDataGroup(mkn$Comp.1$analysis$initDGs)
mkn$Comp.1$analysis$oDG <-
  hzar.copyModelLabels(mkn$Comp.1$analysis$initDGs,
                       mkn$Comp.1$analysis$oDG)

## Convert all 90 runs to hzar.dataGroup objects, adding them to
## the hzar.obsDataGroup object.
mkn$Comp.1$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(mkn$Comp.1$runs$chains,
                                hzar.dataGroup.add),
                         mkn$Comp.1$analysis$oDG);

## Check to make sure that there are only ten hzar.dataGroup
## objects in the hzar.obsDataGroup object.
print(summary(mkn$Comp.1$analysis$oDG$data.groups))


## Look at the variation in parameters for cline models
oDGkey <- which(!(names(mkn$Comp.1$analysis$oDG$data.groups) %in%
                    "nullModel"));
print(hzar.getLLCutParam(mkn$Comp.1$analysis$oDG$data.groups[oDGkey ],
                         c("center","width")));

## Compare the 3 cline models to the null model graphically
##hzar.plot.cline(mkn$Comp.1$analysis$oDG);

## Do model selection based on the AICc scores
print(mkn$Comp.1$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(mkn$Comp.1$analysis$oDG));

## Print out the model with the minimum AICc score
print(mkn$Comp.1$analysis$model.name <-
        rownames(mkn$Comp.1$analysis$AICcTable
        )[[ which.min(mkn$Comp.1$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
mkn$Comp.1$analysis$model.selected <-
  mkn$Comp.1$analysis$oDG$data.groups[[mkn$Comp.1$analysis$model.name]]

## Plot the maximum likelihood cline for the selected model
##hzar.plot.cline(mkn$Comp.1$analysis$model.selected);

print(hzar.get.ML.cline(mkn$Comp.1$analysis$model.selected))

## Plot the 95% credible cline region for the selected model
##hzar.plot.fzCline(mkn$Comp.1$analysis$model.selected);

par(mar=c(5,5,4,2)+.1)
par(ps=22, cex=1, cex.main=1)
male.pc1<-hzar.plot.fzCline(mkn$Comp.1$analysis$model.selected,
                            pch=19,
                            xlab="Distance (km)",
                            ylab="female PC1",
                            ylim=c(-3,3.5),xlim=c(0,1200))

#export as 8"x8" PDF
#dev.off()
## End Comp.1



## We are doing just the one quantitative trait, but it is
## good to stay organized.
mkn$Comp.1 <- list();
## Space to hold the observed data
mkn$Comp.1$obs <- list();
## Space to hold the models to fit
mkn$Comp.1$models <- list();
## Space to hold the compiled fit requests
mkn$Comp.1$fitRs <- list();
## Space to hold the output data chains
mkn$Comp.1$runs <- list();
## Space to hold the analysed data
mkn$Comp.1$analysis <- list()

## Beard Length Trait from Brumfield et al 2001
mkn$Comp.1$obs <-
  hzar.doNormalData1DRaw(hzar.mapSiteDist(SiteLoc$Site,
                                          SiteLoc$Distance),
                         S1blue$Site,
                         S1blue$S1.blue)

## Look at a graph of the observed data
hzar.plot.obsData(mkn$Comp.1$obs);

## Make a helper function
mkn.loadComp.1model <- function(scaling,tails,
                                id=paste(scaling,tails,sep=".")){
  mkn$Comp.1$models[[id]] <<-
    hzar.makeCline1DNormal(mkn$Comp.1$obs, tails)
  ## As there is no quick option for "fixed" scaling, and the
  ## combined sites "A", "B" and "J","K" have a fair number of samples (> 20),
  ## fix the mean and variance of the left and right sides of
  ## the cline to the values observed by the combined sites.
  if (all(regexpr("fixed",scaling,ignore.case=TRUE) == 1 )){
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$muL <<- TRUE
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$muR <<- TRUE
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$varL <<- TRUE
    hzar:::meta.fix(mkn$Comp.1$models[[id]])$varR <<- TRUE
  }
  ## Helper function to work around low sample size
  getCombo <- function(id1,id2,colName,frame=mkn$Comp.1$obs$frame)
    (frame[id1,colName]*frame[id1,"nEff"]+
       frame[id2,colName]*frame[id2,"nEff"])/
    (frame[id1,"nEff"]+frame[id2,"nEff"])
  
  ## Site A, B is the "left" side of the cline, so pull the
  ## fixed values from there.
  hzar:::meta.init(mkn$Comp.1$models[[id]])$muL <<-
    getCombo("A","B","mu")
  hzar:::meta.init(mkn$Comp.1$models[[id]])$varL <<-
    getCombo("A","B","var")
  ## Site J, K is the "right" side of the cline, so pull the
  ## fixed values from there.
  hzar:::meta.init(mkn$Comp.1$models[[id]])$muR <<-
    getCombo("I","J","mu")
  hzar:::meta.init(mkn$Comp.1$models[[id]])$varR <<-
    getCombo("I","J","var")
  ## Make a better estimate of varH using site D, E
  hzar:::meta.init(mkn$Comp.1$models[[id]])$varH <<-
    getCombo("E","F","var")-
    (getCombo("I","J","var")+getCombo("A","B","var"))/2
  
}
## mkn.loadComp.1model("fixed","none","modelI");
## mkn.loadComp.1model("free" ,"none","modelII");
## mkn.loadComp.1model("free" ,"both","modelIII");
mkn.loadComp.1model("fixed","none"  ,"fixN");
mkn.loadComp.1model("fixed","left"  ,"fixL");
mkn.loadComp.1model("fixed","right" ,"fixR");
mkn.loadComp.1model("fixed","mirror","fixM");
mkn.loadComp.1model("fixed","both"  ,"fixB");
mkn.loadComp.1model("free" ,"none"  ,"optN");
mkn.loadComp.1model("free" ,"left"  ,"opComp.1");
mkn.loadComp.1model("free" ,"right" ,"optR");
mkn.loadComp.1model("free" ,"mirror","optM");
mkn.loadComp.1model("free" ,"both"  ,"optB");


## Check the default settings
##print(mkn$Comp.1$models)

## Modify all models to focus on the region where the observed
## data were collected.
## Observations were between -10 and 130000 km.
mkn$Comp.1$models <- sapply(mkn$Comp.1$models,
                            hzar.model.addBoxReq,
                            -10 , 1100,
                            simplify=FALSE)

## Due to the large number of free variables, it is prudent to
## reduce the tune setting of optB from 1.5 to 1.1
hzar:::meta.tune(mkn$Comp.1$models$optB)<-1.1


## Check the updated settings
##print(mkn$Comp.1$models)

## Compile each of the models to prepare for fitting
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
mkn$Comp.1$fitRs$init <- sapply(mkn$Comp.1$models,
                                hzar.first.fitRequest.gC,
                                obsData=mkn$Comp.1$obs,
                                verbose=TRUE,
                                simplify=FALSE)

## Update the settings for the fitter if desired.
mkn$Comp.1$fitRs$init <- sapply(mkn$Comp.1$fitRs$init,
                                function(mdl) {
                                  mdl$mcmcParam$chainLength <-
                                    chainLength; #1e5
                                  mdl$mcmcParam$burnin <-
                                    chainLength %/% 10; #1e4
                                  mdl },
                                simplify=FALSE)

mkn$Comp.1$fitRs$init <-
  rotateModelSeeds(mkn$Comp.1$fitRs$init)


## Check fit request settings
##print(mkn$Comp.1$fitRs$init)

## Do just one run of the models for an initial chain
mkn$Comp.1$runs$init <-
  hzar.doFit.multi(mkn$Comp.1$fitRs$init)
names(mkn$Comp.1$runs$init) <- names(mkn$Comp.1$fitRs$init)



## Compile a new set of fit requests using the initial chains 
mkn$Comp.1$fitRs$chains <-
  lapply(mkn$Comp.1$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
mkn$Comp.1$fitRs$chains <-
  hzar.multiFitRequest(mkn$Comp.1$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit


## Go ahead and run a chain of 3 runs for every fit request
mkn$Comp.1$runs$chains <-  hzar.doChain.multi(mkn$Comp.1$fitRs$chains,
                                              doPar=TRUE,
                                              inOrder=FALSE,
                                              count=3)

## Did fixN converge?
summary(do.call(mcmc.list,
                lapply(mkn$Comp.1$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did optN converge?
summary(do.call(mcmc.list,
                lapply(mkn$Comp.1$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did optB converge?
summary(do.call(mcmc.list,
                lapply(mkn$Comp.1$runs$chains[28:30],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Start aggregation of data for analysis

## Clear out a spot to collect the data for analysis (note that
## there is currenComp.1y no "null model" to compare against).
mkn$Comp.1$analysis$initDGs <- list(
)

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
mkn$Comp.1$analysis$initDGs <-
  c( mkn$Comp.1$analysis$initDGs,
     sapply(mkn$Comp.1$runs$init,
            hzar.dataGroup.add,
            simplify=FALSE))

## Create a hzar.obsDataGroup object from the four hzar.dataGroup
## just created, copying the naming scheme (modelI, modelII,
## modelIII).
mkn$Comp.1$analysis$oDG <-
  hzar.make.obsDataGroup(mkn$Comp.1$analysis$initDGs)
mkn$Comp.1$analysis$oDG <-
  hzar.copyModelLabels(mkn$Comp.1$analysis$initDGs,
                       mkn$Comp.1$analysis$oDG)

## Convert all 90 runs to hzar.dataGroup objects, adding them to
## the hzar.obsDataGroup object.
mkn$Comp.1$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(mkn$Comp.1$runs$chains,
                                hzar.dataGroup.add),
                         mkn$Comp.1$analysis$oDG);

## Check to make sure that there are only ten hzar.dataGroup
## objects in the hzar.obsDataGroup object.
print(summary(mkn$Comp.1$analysis$oDG$data.groups))


## Look at the variation in parameters for cline models
oDGkey <- which(!(names(mkn$Comp.1$analysis$oDG$data.groups) %in%
                    "nullModel"));
print(hzar.getLLCutParam(mkn$Comp.1$analysis$oDG$data.groups[oDGkey ],
                         c("center","width")));

## Compare the 3 cline models to the null model graphically
##hzar.plot.cline(mkn$Comp.1$analysis$oDG);

## Do model selection based on the AICc scores
print(mkn$Comp.1$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(mkn$Comp.1$analysis$oDG));

## Print out the model with the minimum AICc score
print(mkn$Comp.1$analysis$model.name <-
        rownames(mkn$Comp.1$analysis$AICcTable
        )[[ which.min(mkn$Comp.1$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
mkn$Comp.1$analysis$model.selected <-
  mkn$Comp.1$analysis$oDG$data.groups[[mkn$Comp.1$analysis$model.name]]

## Plot the maximum likelihood cline for the selected model
##hzar.plot.cline(mkn$Comp.1$analysis$model.selected);

print(hzar.get.ML.cline(mkn$Comp.1$analysis$model.selected))

## Plot the 95% credible cline region for the selected model
##hzar.plot.fzCline(mkn$Comp.1$analysis$model.selected);

par(mar=c(5,5,4,2)+.1)
par(ps=22, cex=1, cex.main=1)
male.pc1<-hzar.plot.fzCline(mkn$Comp.1$analysis$model.selected,
                            pch=19,
                            xlab="Distance (km)",
                            ylab="S1.blue",
                            ylim=c(.22, .35),xlim=c(0,1200))

#export as 8"x8" PDF
#dev.off()
## End Comp.1


