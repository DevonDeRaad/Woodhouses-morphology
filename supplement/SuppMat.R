library(ggbiplot)
library(ggridges)
library(grid)
library(gridExtra)

#data analysis for supplementary material
#Phenotypic clines in Woodhouse's Scrub-Jay

#bring in data
scrubspecmorph<-read.csv("~/Dropbox/scrub poster/manuscript/scrubspecmorph.csv")
head(scrubspecmorph)

###start analyzing morph
##normalitycheck
hist(scrubspecmorph$wing)
hist(scrubspecmorph$tail)
hist(scrubspecmorph$tarsus)
hist(scrubspecmorph$bl)
hist(scrubspecmorph$bw)
hist(scrubspecmorph$bd)
hist(scrubspecmorph$S1.blue)
#visually verified normality in the distribution of each trait

#subset by sex
malescrub<- subset(scrubspecmorph, scrubspecmorph$sex=="male")
femalescrub<-subset(scrubspecmorph, scrubspecmorph$sex=="female")

#subset by group
woodhousesgroup<- subset(scrubspecmorph, scrubspecmorph$locality=="woodhouse")
sumichrastsgroup<-subset(scrubspecmorph, scrubspecmorph$locality=="sumichrasti")

#subset by subspecies
gris<- subset(scrubspecmorph, scrubspecmorph$subspecies=="grisea")
cyanotis<- subset(scrubspecmorph, scrubspecmorph$subspecies=="cyanotis")
hyb<- subset(scrubspecmorph, scrubspecmorph$subspecies=="sumichrasti/cyanotis")
sumi<- subset(scrubspecmorph, scrubspecmorph$subspecies=="sumichrasti")
remota<- subset(scrubspecmorph, scrubspecmorph$subspecies=="remota")


std.error(malescrub$wing)
std.error(malescrub$tail)
std.error(malescrub$tarsus)
std.error(malescrub$bl)
std.error(malescrub$bw)
std.error(malescrub$bd)
std.error(malescrub$S1.blue)

std.error(femalescrub$wing)
std.error(femalescrub$tail)
std.error(femalescrub$tarsus)
std.error(femalescrub$bl)
std.error(femalescrub$bw)
std.error(femalescrub$bd)
std.error(femalescrub$S1.blue)

std.error(woodhousesgroup$wing)
std.error(woodhousesgroup$tail)
std.error(woodhousesgroup$tarsus)
std.error(woodhousesgroup$bl)
std.error(woodhousesgroup$bw)
std.error(woodhousesgroup$bd)
std.error(woodhousesgroup$S1.blue)

std.error(sumichrastsgroup$wing)
std.error(sumichrastsgroup$tail)
std.error(sumichrastsgroup$tarsus)
std.error(sumichrastsgroup$bl)
std.error(sumichrastsgroup$bw)
std.error(sumichrastsgroup$bd)
std.error(sumichrastsgroup$S1.blue)

std.error(gris$wing)
std.error(gris$tail)
std.error(gris$tarsus)
std.error(gris$bl)
std.error(gris$bw)
std.error(gris$bd)
std.error(gris$S1.blue)

std.error(cyanotis$wing)
std.error(cyanotis$tail)
std.error(cyanotis$tarsus)
std.error(cyanotis$bl)
std.error(cyanotis$bw)
std.error(cyanotis$bd)
std.error(cyanotis$S1.blue)

std.error(hyb$wing)
std.error(hyb$tail)
std.error(hyb$tarsus)
std.error(hyb$bl)
std.error(hyb$bw)
std.error(hyb$bd)
std.error(hyb$S1.blue)

std.error(sumi$wing)
std.error(sumi$tail)
std.error(sumi$tarsus)
std.error(sumi$bl)
std.error(sumi$bw)
std.error(sumi$bd)
std.error(sumi$S1.blue)

std.error(remota$wing)
std.error(remota$tail)
std.error(remota$tarsus)
std.error(remota$bl)
std.error(remota$bw)
std.error(remota$bd)
std.error(remota$S1.blue)

#t.test for differences in each trait between sexes
t.test(malescrub$wing, femalescrub$wing)
t.test(malescrub$tail, femalescrub$tail)
t.test(malescrub$tarsus, femalescrub$tarsus)
t.test(malescrub$bl, femalescrub$bl)
t.test(malescrub$bw, femalescrub$bw)
t.test(malescrub$bd, femalescrub$bd)
t.test(malescrub$S1.blue, femalescrub$S1.blue)
#everything but back color is significantly different


#subset
sumgroup<- subset(scrubspecmorph, locality == "sumichrasti")
woodgroup<- subset(scrubspecmorph, locality == "woodhouse")
contactintermediateintrogressedhybridgroup<- subset(scrubspecmorph, locality == "contact")

#t.test for differences in each trait between groups
t.test(sumgroup$wing, woodgroup$wing)
t.test(sumgroup$tail, woodgroup$tail)
t.test(sumgroup$tarsus, woodgroup$tarsus)
t.test(sumgroup$bl, woodgroup$bl)
t.test(sumgroup$bw, woodgroup$bw)
t.test(sumgroup$bd, woodgroup$bd)
t.test(sumgroup$S1.blue, woodgroup$S1.blue)

#data vis to see if size and back color are corresponding
plot(contactintermediateintrogressedhybridgroup$wing ~contactintermediateintrogressedhybridgroup$S1.blue)
plot(contactintermediateintrogressedhybridgroup$tail ~contactintermediateintrogressedhybridgroup$S1.blue)
plot(contactintermediateintrogressedhybridgroup$tarsus ~contactintermediateintrogressedhybridgroup$S1.blue)
#looks pretty interesting, try for PC1/S1B

#build PCAscrubspecmorph
pca.scruball <- princomp(scrubspecmorph[,8:14], cor=TRUE)
#view loadings by variable
unclass(pca.scruball$loadings)
#create dataframe of PC scores
pcascoresbyindiv<-data.frame(pca.scruball$scores)
#insert id column from scrubspecmorph into pcascorebyindiv
pcascoresbyindiv$id <- scrubspecmorph$id
#merge datasets to add PC scores to scrubspecmorph
PCAscrubspecmorph<-merge(scrubspecmorph, pcascoresbyindiv, by = "id")
View(PCAscrubspecmorph)

#subset
pcahybs<- subset(PCAscrubspecmorph, locality == "contact")
malepcahybs<- subset(pcahybs, sex == "male")
femalepcahybs<- subset(pcahybs, sex == "female")

#plot
plot(pcahybs$S1.blue ~ pcahybs$Comp.1)
plot(malepcahybs$S1.blue ~ malepcahybs$Comp.1)
plot(femalepcahybs$S1.blue ~ femalepcahybs$Comp.1)

#simple regression
mod<-lm(Comp.1 ~ S1.blue, data = pcahybs)
summary(mod)

ggplot(pcahybs, aes(x = Comp.1, y = S1.blue)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")+
  labs(x= "Body Size (PC1)", y= "Blue Saturation of mantle(S1B)")+
  theme_bw()

#create PCA by sex for only body size
pca.morph<- princomp(scrubspecmorph[, 8:13], cor=TRUE)
#PrintPCA
print(pca.morph)
#summarize PC's
summary(pca.morph)
#visualize the differences between sexes
ggbiplot(pca.morph, obs.scale = 1, var.scale = 1,
         groups = scrubspecmorph$sex, ellipse = TRUE) +
  scale_color_discrete(name = '') +
  theme_bw()
#save as PDF 7"x7"

#use gg_ridges package to visualize differences in general body size and mantle color
##plot male wing length
malewing<-ggplot(malescrub, aes(x = wing, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "male wing length", y = "Group") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c("sumichrast's", "contact", "woodhouse's")) +
  theme_classic() +
  theme(legend.position = "none")

#plot female wing length
femwing<-ggplot(femalescrub, aes(x = wing, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "female wing length", y = "Group") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c("sumichrast's", "contact", "woodhouse's")) +
  theme_classic() +
  theme(legend.position = "none")

##plot male tarsus length
maletars<-ggplot(malescrub, aes(x = tarsus, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "male tarsus length", y = "") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c()) +
  theme_classic() +
  theme(legend.position = "none")

#plot female tarsus length
femtars<-ggplot(femalescrub, aes(x = tarsus, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "female tarsus length", y = "") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c()) +
  theme_classic() +
  theme(legend.position = "none")

##plot male tail length
maletail<-ggplot(malescrub, aes(x = tail, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "male tail length", y = "") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c()) +
  theme_classic() +
  theme(legend.position = "none")

#plot female tail length
femtail<-ggplot(femalescrub, aes(x = tail, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "female tail length", y = "") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c()) +
  theme_classic() +
  theme(legend.position = "none")

##plot male bill length
malebl<-ggplot(malescrub, aes(x = bl, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "male bill length", y = "Group") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c("sumichrast's", "contact", "woodhouse's")) +
  theme_classic() +
  theme(legend.position = "none")

#plot female bill length
fembl<-ggplot(femalescrub, aes(x = bl, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "female bill length", y = "Group") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c("sumichrast's", "contact", "woodhouse's")) +
  theme_classic() +
  theme(legend.position = "none")

##plot male bill width
malebw<-ggplot(malescrub, aes(x = bw, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "male bill width", y = "") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c()) +
  theme_classic() +
  theme(legend.position = "none")

#plot female bill width
fembw<-ggplot(femalescrub, aes(x = bw, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "female bill width", y = "") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c()) +
  theme_classic() +
  theme(legend.position = "none")

##plot male bill depth
malebd<-ggplot(malescrub, aes(x = bd, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "male bill depth", y = "") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c()) +
  theme_classic() +
  theme(legend.position = "none")

#plot female bill depth
fembd<-ggplot(femalescrub, aes(x = bd, y = locality, fill = locality, color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("color", "fill")) +
  labs(x = "female bill depth", y = "") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c()) +
  theme_classic() +
  theme(legend.position = "none")



grid.arrange(malewing, femwing, maletars, femtars, 
             maletail, femtail, malebl, fembl, malebw, 
             fembw, malebd, fembd, ncol = 2)

grid.arrange(malewing, maletars,  
             maletail, malebl, 
             malebw, malebd, ncol = 3)
#save as PDF 8"x5"
grid.arrange(femwing,femtars, 
             femtail, fembl, 
             fembw, fembd, ncol = 3)

#save as PDF 8"x5"





#plot PC1 for the three groups
ggplot(PCAscrubspecmorph, aes(x = Comp.1, y = locality, fill = locality, point_color = locality)) +
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = .4) +
  scale_color_manual(values = c("green4", "navyblue", "darkred"), aesthetics = c("fill", "point_color")) +
  labs(x = "Body Size", y = "Group") +
  scale_y_discrete(limits=c("sumichrasti", "contact", "woodhouse"), labels=c("sumichrast's", "contact", "woodhouse's")) +
  theme_classic() +
  theme(legend.position = "none")




#run MANOVA for Table 2 line 4 (sex)
#Manova for male v female in all 6 characteristics
fem.male <- manova(cbind(wing, tail, tarsus, bl, bw, bd, S1.blue) ~ sex, data = scrubspecmorph)
#get summary F statistic
summary(fem.male, test="Pillai")
#get F statistic and p value for each variable in order
summary.aov(fem.male)



#subset woodhouse's group, sumichrast's group, and bind the two dataframes together to drop contact zone individuals
woodscrub<- subset(scrubspecmorph, scrubspecmorph$locality=="woodhouse")
sumscrub<- subset(scrubspecmorph, scrubspecmorph$locality=="sumichrasti")
dropcontact<-rbind(woodscrub, sumscrub)
#drop unused levels (contact)
dropcontact<-droplevels(dropcontact)
#ttest to compare back color between sumi and woodhouse
t.test(dropcontact$S1.blue~dropcontact$locality)
#sigaf
#visualize
boxplot(S1.blue~locality, data = dropcontact)
#run manova
dropman <- cbind(dropcontact$wing,dropcontact$tail,dropcontact$tarsus,dropcontact$bl,dropcontact$bw,dropcontact$bd,dropcontact$S1.blue)
droppednova<- manova(dropman ~ dropcontact$locality)
#get summary F statistic
summary(droppednova, test="Pillai")
#get F statistic and p value for each variable in order
summary.aov(droppednova)



#subset dataset by sex
malescrub<- subset(scrubspecmorph, scrubspecmorph$sex=="male")
femalescrub<-subset(scrubspecmorph, scrubspecmorph$sex=="female")
#plot males
plot(malescrub$lat, malescrub$wing)
plot(malescrub$lat, malescrub$tail)
#plot females
plot(femalescrub$lat, femalescrub$wing)
plot(femalescrub$lat, femalescrub$tail)



#plot winglength v latitude
plot(scrubspecmorph$lat, scrubspecmorph$wing)
#plot taillength v latitude
plot(scrubspecmorph$lat, scrubspecmorph$tail)


###investigate sexual dimorphism identified in pitelka
#sex by tail
tailsex<-ggplot(scrubspecmorph, aes(sex, tail)) + 
  geom_boxplot() +
  ylab("tail length (mm)")+
  xlab("sex")
#plot
tailsex
#sex by wing
wingsex<-ggplot(scrubspecmorph, aes(sex, wing)) + 
  geom_boxplot() +
  ylab("wing length (mm)")+
  xlab("sex")
#plot
wingsex
##males larger than females
#Manova for male v female in all 6 characteristics
fem.male <- manova(cbind(wing, tail, tarsus, bl, bw, bd) ~ sex, data = scrubspecmorph)
summary(fem.male)
##significantly different in morphology based on sex
##males are generally larger bodied
##morph PCA based on sex
pca.morph<- princomp(scrubspecmorph[, 8:14], corr=TRUE)
#PrintPCA
print(pca.morph)
#summarize PC's
summary(pca.morph)
#visualize the differences between sexes
ggbiplot(pca.morph, obs.scale = 1, var.scale = 1,
         groups = scrubspecmorph$sex, ellipse = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal')
