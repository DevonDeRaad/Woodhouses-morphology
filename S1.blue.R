setwd("~/Desktop/hzar/hzar.scrubbers")

## Load the package
library(hzar);

if(require(doMC)){
  registerDoMC();
}else if(require(doParallel)){
  registerDoParallel()
}else{
  registerDoSeq()
}

## A typical chain length
chainLength=1e5;                       

## Make each model run off a separate seed
## mainSeed=
##   list(
##        A=c(978,544,99,596,528,124),
##        B=c(544,99,596,528,124,978),
##        C=c(99,596,528,124,978,544))

## Increment skip by 15 or so for every trait
rotateModelSeeds <- function(fiS1.blue)
  hzar.multiFitRequest(fiS1.blue,
                       rotateSeed=TRUE, skip=500,
                       baseChannel=NULL, each=1,
                       baseSeed=c(596,528,124,978,544,99))
## ## Save all plots in a pdf file?
## pdf(width=5, height=7, file="ManakinExamplePlots.pdf")

read.csv("morph.csv")->morph
read.csv("siteloc.csv")->SiteLoc

## Trait Analysis

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
mkn$S1.blue <- list();
## Space to hold the observed data
mkn$S1.blue$obs <- list();
## Space to hold the models to fit
mkn$S1.blue$models <- list();
## Space to hold the compiled fit requests
mkn$S1.blue$fitRs <- list();
## Space to hold the output data chains
mkn$S1.blue$runs <- list();
## Space to hold the analysed data
mkn$S1.blue$analysis <- list();


## Beard Length Trait from Brumfield et al 2001
mkn$S1.blue$obs <-
  hzar.doNormalData1DRaw(hzar.mapSiteDist(SiteLoc$Site,
                                          SiteLoc$Distance),
                         morph$Site,
                         morph$S1.blue)

## Look at a graph of the observed data
hzar.plot.obsData(mkn$S1.blue$obs);

## Make a helper function
mkn.loadS1.bluemodel <- function(scaling,tails,
                              id=paste(scaling,tails,sep=".")){
  mkn$S1.blue$models[[id]] <<-
    hzar.makeCline1DNormal(mkn$S1.blue$obs, tails)
  ## As there is no quick option for "fixed" scaling, and the
  ## combined sites "A", "B" and "J","K" have a fair number of samples (> 20),
  ## fix the mean and variance of the left and right sides of
  ## the cline to the values observed by the combined sites.
  if (all(regexpr("fixed",scaling,ignore.case=TRUE) == 1 )){
    hzar:::meta.fix(mkn$S1.blue$models[[id]])$muL <<- TRUE
    hzar:::meta.fix(mkn$S1.blue$models[[id]])$muR <<- TRUE
    hzar:::meta.fix(mkn$S1.blue$models[[id]])$varL <<- TRUE
    hzar:::meta.fix(mkn$S1.blue$models[[id]])$varR <<- TRUE
  }
  ## Helper function to work around low sample size
  getCombo <- function(id1,id2,colName,frame=mkn$S1.blue$obs$frame)
    (frame[id1,colName]*frame[id1,"nEff"]+
     frame[id2,colName]*frame[id2,"nEff"])/
       (frame[id1,"nEff"]+frame[id2,"nEff"])
  
  ## Site A, B is the "left" side of the cline, so pull the
    ## fixed values from there.
  hzar:::meta.init(mkn$S1.blue$models[[id]])$muL <<-
    getCombo("A","B","mu")
  hzar:::meta.init(mkn$S1.blue$models[[id]])$varL <<-
    getCombo("A","B","var")
  ## Site J, K is the "right" side of the cline, so pull the
  ## fixed values from there.
  hzar:::meta.init(mkn$S1.blue$models[[id]])$muR <<-
    getCombo("I","J","mu")
  hzar:::meta.init(mkn$S1.blue$models[[id]])$varR <<-
    getCombo("I","J","var")
  ## Make a better estimate of varH using site D, E
  hzar:::meta.init(mkn$S1.blue$models[[id]])$varH <<-
    getCombo("E","F","var")-
      (getCombo("I","J","var")+getCombo("A","B","var"))/2

}
## mkn.loadS1.bluemodel("fixed","none","modelI");
## mkn.loadS1.bluemodel("free" ,"none","modelII");
## mkn.loadS1.bluemodel("free" ,"both","modelIII");
mkn.loadS1.bluemodel("fixed","none"  ,"fixN");
mkn.loadS1.bluemodel("fixed","left"  ,"fixL");
mkn.loadS1.bluemodel("fixed","right" ,"fixR");
mkn.loadS1.bluemodel("fixed","mirror","fixM");
mkn.loadS1.bluemodel("fixed","both"  ,"fixB");
mkn.loadS1.bluemodel("free" ,"none"  ,"optN");
mkn.loadS1.bluemodel("free" ,"left"  ,"opS1.blue");
mkn.loadS1.bluemodel("free" ,"right" ,"optR");
mkn.loadS1.bluemodel("free" ,"mirror","optM");
mkn.loadS1.bluemodel("free" ,"both"  ,"optB");


## Check the default settings
##print(mkn$S1.blue$models)

## Modify all models to focus on the region where the observed
## data were collected.
## Observations were between -10 and 130000 km.
mkn$S1.blue$models <- sapply(mkn$S1.blue$models,
                          hzar.model.addBoxReq,
                          -10 , 1100,
                          simplify=FALSE)

## Due to the large number of free variables, it is prudent to
## reduce the tune setting of optB from 1.5 to 1.1
hzar:::meta.tune(mkn$S1.blue$models$optB)<-1.1


## Check the updated settings
##print(mkn$S1.blue$models)

## Compile each of the models to prepare for fitting
## Note that we are using hzar.first.fitRequest.gC for fitting
## guassian (aka "normal") clines.
mkn$S1.blue$fitRs$init <- sapply(mkn$S1.blue$models,
                         hzar.first.fitRequest.gC,
                         obsData=mkn$S1.blue$obs,
                         verbose=TRUE,
                         simplify=FALSE)

## Update the settings for the fitter if desired.
mkn$S1.blue$fitRs$init <- sapply(mkn$S1.blue$fitRs$init,
                             function(mdl) {
                               mdl$mcmcParam$chainLength <-
                                 chainLength; #1e5
                               mdl$mcmcParam$burnin <-
                                 chainLength %/% 10; #1e4
                               mdl },
                             simplify=FALSE)

mkn$S1.blue$fitRs$init <-
  rotateModelSeeds(mkn$S1.blue$fitRs$init)


## Check fit request settings
##print(mkn$S1.blue$fitRs$init)

## Do just one run of the models for an initial chain
mkn$S1.blue$runs$init <-
  hzar.doFit.multi(mkn$S1.blue$fitRs$init)
names(mkn$S1.blue$runs$init) <- names(mkn$S1.blue$fitRs$init)



## Compile a new set of fit requests using the initial chains 
mkn$S1.blue$fitRs$chains <-
  lapply(mkn$S1.blue$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
mkn$S1.blue$fitRs$chains <-
  hzar.multiFitRequest(mkn$S1.blue$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit


## Go ahead and run a chain of 3 runs for every fit request
mkn$S1.blue$runs$chains <-  hzar.doChain.multi(mkn$S1.blue$fitRs$chains,
                                            doPar=TRUE,
                                            inOrder=FALSE,
                                            count=3)

## Did fixN converge?
summary(do.call(mcmc.list,
                lapply(mkn$S1.blue$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did optN converge?
summary(do.call(mcmc.list,
                lapply(mkn$S1.blue$runs$chains[16:18],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Did optB converge?
summary(do.call(mcmc.list,
                lapply(mkn$S1.blue$runs$chains[28:30],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

## Start aggregation of data for analysis

## Clear out a spot to collect the data for analysis (note that
## there is currenS1.bluey no "null model" to compare against).
mkn$S1.blue$analysis$initDGs <- list(
  )

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
mkn$S1.blue$analysis$initDGs <-
  c( mkn$S1.blue$analysis$initDGs,
    sapply(mkn$S1.blue$runs$init,
           hzar.dataGroup.add,
           simplify=FALSE))

## Create a hzar.obsDataGroup object from the four hzar.dataGroup
## just created, copying the naming scheme (modelI, modelII,
## modelIII).
mkn$S1.blue$analysis$oDG <-
  hzar.make.obsDataGroup(mkn$S1.blue$analysis$initDGs)
mkn$S1.blue$analysis$oDG <-
    hzar.copyModelLabels(mkn$S1.blue$analysis$initDGs,
                         mkn$S1.blue$analysis$oDG)

## Convert all 90 runs to hzar.dataGroup objects, adding them to
## the hzar.obsDataGroup object.
mkn$S1.blue$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(mkn$S1.blue$runs$chains,
                                hzar.dataGroup.add),
                         mkn$S1.blue$analysis$oDG);

## Check to make sure that there are only ten hzar.dataGroup
## objects in the hzar.obsDataGroup object.
print(summary(mkn$S1.blue$analysis$oDG$data.groups))


## Look at the variation in parameters for cline models
oDGkey <- which(!(names(mkn$S1.blue$analysis$oDG$data.groups) %in%
                  "nullModel"));
print(hzar.getLLCutParam(mkn$S1.blue$analysis$oDG$data.groups[oDGkey ],
                         c("center","width")));

## Compare the 3 cline models to the null model graphically
##hzar.plot.cline(mkn$S1.blue$analysis$oDG);

## Do model selection based on the AICc scores
print(mkn$S1.blue$analysis$AICcTable <-
      hzar.AICc.hzar.obsDataGroup(mkn$S1.blue$analysis$oDG));

## Print out the model with the minimum AICc score
print(mkn$S1.blue$analysis$model.name <-
  rownames(mkn$S1.blue$analysis$AICcTable
           )[[ which.min(mkn$S1.blue$analysis$AICcTable$AICc )]])

## Extract the hzar.dataGroup object for the selected model
mkn$S1.blue$analysis$model.selected <-
  mkn$S1.blue$analysis$oDG$data.groups[[mkn$S1.blue$analysis$model.name]]

## Plot the maximum likelihood cline for the selected model
##hzar.plot.cline(mkn$S1.blue$analysis$model.selected);

print(hzar.get.ML.cline(mkn$S1.blue$analysis$model.selected))

## Plot the 95% credible cline region for the selected model
##hzar.plot.fzCline(mkn$S1.blue$analysis$model.selected);

par(mar=c(5,5,4,2)+.1)
par(ps=22, cex=1, cex.main=1)
hzar.plot.fzCline(mkn$S1.blue$analysis$model.selected,
                  pch=19,
                  xlab="Distance (km)",
                  ylab="Blue Saturation",
                  ylim=c(0.24,0.32),xlim=c(0,1200))

#export PDF 8"x8"
## End Quantitative Trait Analysis of S1B
dev.off()
#
