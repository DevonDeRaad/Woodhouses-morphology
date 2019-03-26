###last update: 16 March 2019
#scrub jay morphological analysis r script
# Create Figures 3 & 4 from the manuscript "Phenotypic clines across an unstudied hybrid zone in Woodhouse's Scrub-Jay (Aphelocoma woodhouseii)"
##libraries
library(ggplot2)
library(ggbiplot)
library(grid)
library(gridExtra)
library(rworldmap)
library(maps)
library(mapdata)
library(ggmap)
library(ggridges)

#read in data
scrubspecmorph <- read.csv("~/Dropbox/scrub poster/manuscript/scrubspecmorph.csv", header=T)

#######
###start analyzing morph
##normalitycheck
hist(scrubspecmorph$wing)
hist(scrubspecmorph$tail)
hist(scrubspecmorph$tarsus)
hist(scrubspecmorph$bl)
hist(scrubspecmorph$bw)
hist(scrubspecmorph$bd)
##all six characteristics look normalish

###compare wing length between subspecies
subspecwing<-ggplot(scrubspecmorph, aes(subspecies, wing)) + 
  geom_boxplot() +
  ylab("wing length (mm)")+
  xlab("subspec")
#plot
subspecwing

###compare lineage level differences in wing length
specwing<-ggplot(scrubspecmorph, aes(locality, wing)) + 
  geom_boxplot() +
  ylab("wing length (mm)")+
  xlab("group")
#plot
specwing

###
###investigate sexual dimorphism identified by pitelka (1951)
#compare tail length by sex
tailsex<-ggplot(scrubspecmorph, aes(sex, tail)) + 
  geom_boxplot() +
  ylab("tail length (mm)")+
  xlab("sex")
#plot
tailsex #noticeable difference

#wing length by sex
wingsex<-ggplot(scrubspecmorph, aes(sex, wing)) + 
  geom_boxplot() +
  ylab("wing length (mm)")+
  xlab("sex")
#plot
wingsex #males slightly larger once again


#Manova for male v female in all 6 characteristics
fem.male <- manova(cbind(wing, tail, tarsus, bl, bw, bd) ~ sex, data = scrubspecmorph)
summary(fem.male)

##significant differences in morphology based on sex
##males are generally larger bodied
##morph PCA based on sex
pca.morph<- prcomp(scrubspecmorph[, 8:13])
#PrintPCA
print(pca.morph)
#summarize PC's
summary(pca.morph)

#visualize the differences between sexes
ggbiplot(pca.morph, obs.scale = 1, var.scale = 1,
         groups = scrubspecmorph$sex, ellipse = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal')

##overall, there is a detectable difference in body size (especially wing and tail) between the sexes
#but the difference is variable, and the groups are largely overlapping.
#Because we have roughly equivalent sex ratios between lineages and subspecies, and the size differences
#are largely overlapping, we will keep all individuals together for manuscript analyses.

#####################################
##########################
#####################################
#plot winglength v latitude
plot(scrubspecmorph$lat, scrubspecmorph$wing)
#plot taillength v latitude
plot(scrubspecmorph$lat, scrubspecmorph$tail)

#negative relationship, meaning increasing wing/tail length moving north to south
#subset dataset by sex
malescrub<- subset(scrubspecmorph, scrubspecmorph$sex=="male")
femalescrub<-subset(scrubspecmorph, scrubspecmorph$sex=="female")
#plot males
plot(malescrub$lat, malescrub$wing)
plot(malescrub$lat, malescrub$tail)
#plot females
plot(femalescrub$lat, femalescrub$wing)
plot(femalescrub$lat, femalescrub$tail)
###relationships between size and geography hold whether divided by sex or not
###
###
#plot winglength v latitude color coded for sex
plot(scrubspecmorph$lat, scrubspecmorph$wing, col = scrubspecmorph$sex)
#plot taillength v latitude color coded for sex
plot(scrubspecmorph$lat, scrubspecmorph$tail, col = scrubspecmorph$sex)
#####size differences between the sexes noticeable but pattern identical
#####
#check out the spec variables based on subspecies
#boxplot based on the spec variable S1.blue (blue saturation of the mantle)
S1.blueboxgroup<-ggplot(scrubspecmorph, aes(locality, S1.blue)) + 
  geom_boxplot() +
  ylab("s1.blue")+
  xlab("group")
S1.blueboxgroup
#boxplot based on the spec variable S1.red (red saturation of the mantle)
S1.redboxgroup<-ggplot(scrubspecmorph, aes(locality, S1.red)) + 
  geom_boxplot() +
  ylab("s1.red")+
  xlab("group")
S1.redboxgroup
###
###double dipping on saturation on either end of the spectrum for a single plumage patch is pseudoreplication
###so we will drop S1.red


###PCA of morph data all specimens
pca.scrubmorph <- princomp(scrubspecmorph[,8:13])
#visualize
ggbiplot(pca.scrubmorph, obs.scale = 1, var.scale = 1,
         groups = scrubspecmorph$locality, ellipse = FALSE,
         var.axes = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'vertical')
#this is dope, super well broken out, with contact in the middle
#does seem to be basically a smooth cline from small to large
###PCA of morph data by sex
#subset
malescrub<- subset(scrubspecmorph, scrubspecmorph$sex=="male")
femalescrub<-subset(scrubspecmorph, scrubspecmorph$sex=="female")
#run
pca.malescrubmorph <- prcomp(malescrub[,8:13])
pca.femalescrubmorph <- prcomp(femalescrub[,8:13])
#visualize male
ggbiplot(pca.malescrubmorph, obs.scale = 1, var.scale = 1,
         groups = malescrub$locality, ellipse = FALSE,
         var.axes = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'vertical')
#visualize female
ggbiplot(pca.femalescrubmorph, obs.scale = 1, var.scale = 1,
         groups = femalescrub$locality, ellipse = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal')
#####
#####both show the same pattern of separation
#####
#Poster PCA including morph and spec variables
#####
#PCA of all variables morph + blue color using correlative matrix
#####PCA of all morphology only using blue not red
pca.scruball <- princomp(scrubspecmorph[,8:14], cor=TRUE)
#view loadings by variable
unclass(pca.scruball$loadings)
#create dataframe of PC scores
pcascoresbyindiv<-data.frame(pca.scruball$scores)
#insert id column from scrubspecmorph into pcascorebyindiv
pcascoresbyindiv$id <- scrubspecmorph$id
#merge datasets by the column ID to ensure that the correct values stay with correct individuals
PCAscrubspecmorph<-merge(scrubspecmorph, pcascoresbyindiv, by = "id")
View(PCAscrubspecmorph)
#now we can use the values of PC1 as a proxy for general body size in downstream analyses via 'PCAscrubspecmorph'

#visualize PCA of all 6 morph variables + blue saturation of mantle by lineage
ggbiplot(pca.scruball, obs.scale = 1, var.scale = 1,
         groups = scrubspecmorph$locality, ellipse = FALSE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal')
#PCA color coded by lineage with sex indicated by shape
ggbiplot(pca.scruball, obs.scale = 1, var.scale = 1,
         groups = scrubspecmorph$locality, ellipse = TRUE) +
  geom_point(aes(colour=scrubspecmorph$locality, shape=scrubspecmorph$sex), size = 3)+
  scale_color_discrete(name ="Locality", 
                       labels = c("Contact zone", "Sumichrasti group", "Woodhouse's group")) +
  theme(legend.direction = 'vertical') +
  scale_shape_discrete(name = "Sex", labels = c("Female", "Male"))+
  theme_bw()
##PCA color coded by subspecies
ggbiplot(pca.scruball, obs.scale = 1, var.scale = 1,
         groups = scrubspecmorph$subspecies, ellipse = FALSE) +
  geom_point(aes(colour=scrubspecmorph$subspecies, shape=scrubspecmorph$sex), size = 3)+
  scale_color_discrete(name ="Subspecies") +
  scale_color_manual(values=c("cornflowerblue", "blue3", "lightcoral", "tan4", "burlywood3"),
                     name="Subspecies",
                     limits=c("grisea", "cyanotis", "sumichrasti/cyanotis", "sumichrasti", "remota"))+
  scale_shape_discrete(name = "Sex", labels = c("Female", "Male"))+
  theme_bw()

##PCA of only individuals across contact zone transect
contactscrub<-scrubspecmorph[!(scrubspecmorph$subspecies %in% c("grisea", "remota")), ]
##PCA using only sumi-cy-hyb
pca.sumicyhyb <- princomp(contactscrub[,8:14], cor=TRUE)
#visualize
ggbiplot(pca.sumicyhyb, obs.scale = 1, var.scale = 1,
         groups = contactscrub$subspecies, ellipse = FALSE) +
  geom_point(aes(colour=contactscrub$subspecies, shape=contactscrub$sex), size = 3)+
  theme(legend.direction = 'vertical') +
  scale_shape_discrete(name = "Sex", labels = c("Female", "Male"))+
  theme_bw()

#
##
###Make manuscript Figure 3###
##
#
PCAallsubspecies<-ggbiplot(pca.scruball, obs.scale = 1, var.scale = 1,
                           groups = scrubspecmorph$subspecies, ellipse = FALSE) +
  geom_point(aes(colour=scrubspecmorph$subspecies, shape=scrubspecmorph$sex), size = 3)+
  scale_color_discrete(name ="Subspecies") +
  scale_color_manual(values=c("tomato2", "darkred", "green4", "navyblue", "cornflowerblue"),
                     name="Subspecies",
                     limits=c("grisea", "cyanotis", "sumichrasti/cyanotis", "sumichrasti", "remota"))+
  scale_shape_discrete(name = "Sex", labels = c("Female", "Male"))+
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2))+
  theme_bw()+
  theme(legend.text=element_text(size=14, face = "italic"),
        legend.position =c(.18,.18),
        legend.background = element_rect(fill="transparent"),
        legend.title=element_text(size=16))
#print
PCAallsubspecies


PCAaxestwoandthree<-ggbiplot(pca.scruball, choices = 2:3, obs.scale = 1, var.scale = 1,
                           groups = scrubspecmorph$subspecies, ellipse = FALSE) +
  geom_point(aes(colour=scrubspecmorph$subspecies, shape=scrubspecmorph$sex), size = 3)+
  scale_color_discrete(name ="Subspecies") +
  scale_color_manual(values=c("tomato2", "darkred", "green4", "navyblue", "cornflowerblue"),
                     name="Subspecies",
                     limits=c("grisea", "cyanotis", "sumichrasti/cyanotis", "sumichrasti", "remota"))+
  scale_shape_discrete(name = "Sex", labels = c("Female", "Male"))+
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2))+
  theme_bw()+
  theme(legend.text=element_text(size=14, face = "italic"),
        legend.position =c(.18,.18),
        legend.background = element_rect(fill="transparent"),
        legend.title=element_text(size=16))
PCAaxestwoandthree

#same PCA with slightly adjusted graphics  (used for 2017 AOS poster)
###
posterfigure<-ggbiplot(pca.scruball, obs.scale = 1, var.scale = 1,
                       groups = scrubspecmorph$subspecies, ellipse = FALSE) +
  geom_point(aes(colour=scrubspecmorph$subspecies, shape=scrubspecmorph$sex), size = 3)+
  scale_color_discrete(name ="Subspecies") +
  scale_color_manual(values=c("cornflowerblue", "blue3", "lightcoral", "tan4", "burlywood3"),
                     name="Subspecies",
                     limits=c("grisea", "cyanotis", "sumichrasti/cyanotis", "sumichrasti", "remota"))+
  scale_shape_discrete(name = "Sex", labels = c("Female", "Male"))+
  guides(color = guide_legend(order=1),
         shape = guide_legend(order=2))+
  theme_bw()+
  theme(legend.text=element_text(size=14, face = "italic"),
        legend.position =c(.18,.18),
        legend.background = element_rect(fill="transparent"),
        legend.title=element_text(size=16))

#print
posterfigure



#
##
###Make Figure 4

#separate sexes from the dataframe with PC values added as variables
malescrub<-subset(PCAscrubspecmorph, sex == "male")
femalescrub<-subset(PCAscrubspecmorph, sex == "female")

#plot body size for males
malebodsize<-ggplot(malescrub, aes(x = Comp.1, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("fill", "color")) +
  labs(x = "Male Body Size", y = "Group") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c("sumichrast's", "contact", "woodhouse's")) +
  theme_classic() +
  theme(legend.position = "none")

#plot body size for females
femalebodsize<-ggplot(femalescrub, aes(x = Comp.1, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("fill", "color")) +
  labs(x = "Female Body Size", y = "Group") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c("sumichrast's", "contact", "woodhouse's")) +
  theme_classic() +
  theme(legend.position = "none")

#plot back color for all
colorforall<-ggplot(PCAscrubspecmorph, aes(x = S1.blue, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("fill", "color")) +
  labs(x = "Blue Saturation of the Mantle", y = "Group") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c("sumichrast's", "contact", "woodhouse's")) +
  theme_classic() +
  theme(legend.position = "none")


grid.arrange(malebodsize, femalebodsize, colorforall, ncol = 1)







