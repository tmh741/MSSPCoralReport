library(tidyverse)
library(magrittr)
library(lme4)
library(arm)
library(vegan)
library(car)
library(rstanarm)
library(bayesplot)
library(loo)

## This takes a while to read in.
data <- read.csv("https://ndownloader.figshare.com/files/9335584")

#Clean Data
data2012 <- subset(data,OBS_YEAR >= 2012)
newdata <- subset(data2012,COUNT!=0)
newdata %<>% dplyr::select(REGION,ISLAND,SITE,REEF_ZONE,
                            OBS_YEAR,DEPTH_M,HARD_CORAL,MA,CCA,SAND,
                            OTHER,HABITAT_CODE,
                            MAX_DEPTH_M,MAX_HEIGHT,VISIBILITY_M,SPECIES,COUNT)
clean <- newdata %>% unite(col=Code, c("REGION","ISLAND","SITE",
                                            "REEF_ZONE","OBS_YEAR",
                                            "DEPTH_M","HARD_CORAL","MA",
                                            "CCA","SAND","OTHER",
                                            "HABITAT_CODE","MAX_DEPTH_M",
                                            "MAX_HEIGHT","VISIBILITY_M")
                                ,sep="_",remove=T)
clean %<>% group_by(Code) %>% summarize(RICHNESS=n(),TOTALFISH = sum(COUNT),
                                            SHANNON=diversity(COUNT)) %>%
  mutate(MENHINICK=RICHNESS/sqrt(TOTALFISH))

clean %<>% separate(col=Code, into = c("REGION","ISLAND","SITE",
                                           "REEF_ZONE","OBS_YEAR",
                                           "DEPTH_M","HARD_CORAL","MA",
                                           "CCA","SAND","OTHER",
                                           "HABITAT_CODE","MAX_DEPTH_M",
                                           "MAX_HEIGHT","VISIBILITY_M"), sep="_")
for (i in 1:ncol(clean)){
  clean[,i] <- type.convert(clean[,i])
}

clean$HARD.Z <- zscore(clean$HARD_CORAL)
clean$MA.Z <- zscore(clean$MA)
clean$CCA.Z <- zscore(clean$CCA)
clean$SAND.Z <- zscore(clean$SAND)
clean$OTHER.Z <- zscore(clean$OTHER)
clean$YEARC <- clean$OBS_YEAR - 2014
clean$DEPTH.Z <- zscore(clean$DEPTH_M)

clean$SHANNON.Z <- zscore(clean$SHANNON)


clean2 <- clean[clean$REGION!="PRIA",]
clean3 <- na.omit(dplyr::select(clean2,c(HARD_CORAL, MA, SAND, CCA, DEPTH.Z,YEARC,MAX_HEIGHT, ISLAND, REGION,MENHINICK)))


#Explore outcome variables and relationships.
#How does each parameter relate?
ggplot(clean) + geom_bar(aes(x=REGION))

ggplot(clean) + geom_histogram(aes(x=MENHINICK,fill=REGION)) + xlab("Menhinick Index") + ylab("Total Number")

ggplot(clean) + geom_histogram(aes(x=HARD_CORAL)) + xlab("Hard Coral cover (%)") + ylab("Count")
ggplot(clean) + geom_point(aes(x=DEPTH_M,y=MENHINICK))
ggplot(clean) + geom_histogram(aes(x=DEPTH_M)) + xlab("Hard Coral cover (%)") + ylab("Count")
ggplot(clean) + geom_histogram(aes(x=MAX_DEPTH_M))

ggplot(clean) + geom_point(aes(x=sqrt(TOTALFISH),y=RICHNESS,color=REGION)) + xlab("Counted Fish^0.5") + ylab("Number of Species")
ggplot(clean) + geom_bar(aes(x=REGION,y=TOTALFISH,fill=REGION),position="dodge",stat="identity")

ggplot(clean) + geom_point(aes(x=RICHNESS,y=exp(SHANNON),color=REGION))

ggplot(clean) + geom_point(aes(x=RICHNESS,y=MENHINICK,color=REGION))
ggplot(clean) + geom_point(aes(x=TOTALFISH,y=MENHINICK,color=REGION))

ggplot(clean) + geom_point(aes(x=MENHINICK,y=SHANNON,color=REGION))

ggplot(clean) + geom_bar(aes(x=REGION,y=TOTALFISH),stat="identity")
ggplot(clean) + geom_bar(aes(x=REGION,y=RICHNESS),stat="identity")

## Map data
locations <- data2012 %>% 
  dplyr::select(REGION, ISLAND, SITE,LATITUDE,LONGITUDE) %>% 
  distinct() %>%
  mutate(Lon = abs(LONGITUDE))

islands <- data2012 %>%
  dplyr::select(REGION, ISLAND) %>% 
  distinct() %>%
  group_by(REGION) %>%
  summarize(n=n()) %>%
  ungroup()

counts <- data2012 %>%
  dplyr::select(SPECIES, COUNT) %>% 
  group_by(COUNT) %>%
  summarize(n=n()) %>%
  ungroup()

sites <- data2012 %>%
  dplyr::select(ISLAND,SITE) %>%
  distinct() %>% group_by(ISLAND) %>% summarize(n=n())

longtransform <- function (x){
  if (x<0){
    x <- x + 360
  }
  return(x)
}

for (i in 1:nrow(locations)) {
  locations$test[i] <- longtransform(locations$LONGITUDE[i])
}

##pal1 <- colorFactor(palette=c('red','blue','purple','yellow','green'),domain=locations$REGION)
##leaflet(locations) %>% 
#  setView(lng = 180, lat = 10, zoom = 2.5)  %>% 
#  addTiles() %>% 
#  addCircles(data = locations, lat = ~ LATITUDE, lng = ~ test, weight = 1, 
#             radius = 30, popup = ~as.character(REGION),
#             label = ~as.character(REGION), 
#             color = ~ pal1(REGION), 
#             fillOpacity = 1)

ggplot(locations) + geom_point(aes(y=LATITUDE,x=test,color=REGION))

ggplot(clean3) + geom_histogram(aes(x=MENHINICK,fill=REGION),position="identity") + theme_bw() + facet_wrap(~REGION)



#Transform data
zscore <- function(x) {
  z <- (x - mean(x))/sd(x)
  return(z)
}

clean$HARD.Z <- zscore(clean$HARD_CORAL)
clean$MA.Z <- zscore(clean$MA)
clean$CCA.Z <- zscore(clean$CCA)
clean$SAND.Z <- zscore(clean$SAND)
clean$OTHER.Z <- zscore(clean$OTHER)
clean$YEARC <- clean$OBS_YEAR - 2014
clean$DEPTH.Z <- zscore(clean$DEPTH_M)

clean$SHANNON.Z <- zscore(clean$SHANNON)

clean2 <- clean[clean$REGION!="PRIA",]
clean3 <- na.omit(dplyr::select(clean2,c(HARD_CORAL, MA, SAND, CCA, DEPTH.Z,YEARC,MAX_HEIGHT, ISLAND, REGION,MENHINICK)))
clean3 <- clean3[order(clean3$REGION),]

ggplot(clean) + geom_histogram(aes(x=MENHINICK))

#Model pick. Which one?
men1 <- lmer(MENHINICK ~ HARD_CORAL + CCA + MA + SAND + DEPTH.Z + 
               YEARC + VISIBILITY_M + MAX_HEIGHT + REEF_ZONE +
                 (1|ISLAND),data=clean3)
men2 <- lmer(MENHINICK ~ HARD_CORAL + CCA + MA + SAND + DEPTH.Z + 
               YEARC + VISIBILITY_M + MAX_HEIGHT + REEF_ZONE +
              (1|ISLAND:REGION),data=clean3)
men3 <- lmer(MENHINICK ~ HARD_CORAL + CCA + MA + SAND + DEPTH.Z + 
               YEARC + VISIBILITY_M + MAX_HEIGHT +
               (1|ISLAND),data=clean3)
men4 <- lmer(MENHINICK ~ HARD_CORAL + CCA + MA + SAND + DEPTH.Z + 
               YEARC + MAX_HEIGHT +
               (1|ISLAND),data=clean3)
men5 <- lmer(MENHINICK ~ 
               I(HARD_CORAL/10) + I(MA/10) + I(SAND/10) + I(CCA/10) +
               DEPTH.Z + YEARC + MAX_HEIGHT +
               (1|ISLAND),data=clean3)
men6 <- lmer(MENHINICK ~ 
               RICHNESS + I(HARD_CORAL/10) + I(MA/10) + I(SAND/10) + I(CCA/10) +
               DEPTH.Z + YEARC + MAX_HEIGHT + 
               (1|ISLAND),data=clean3)


ggplot(clean) + geom_point(aes(y=RICHNESS,x=sqrt(TOTALFISH)))

AIC(men5,men6,men7)


###



men7 <- lmer(MENHINICK ~ 
               I(HARD_CORAL/10) + I(MA/10) + I(SAND/10) + I(CCA/10) +
               DEPTH.Z + YEARC + MAX_HEIGHT + 
               (1|ISLAND:REGION),data=clean3)

summary(men7) 
plot(fitted(men7),resid(men7))
ggplot() + geom_point(aes(x=fitted(men7),y=resid(men7),color=clean3$REGION))
qqnorm(resid(men7))
qqline(resid(men7))
binnedplot(fitted(men7),resid(men7,type="pearson"))
###

#Cook's Distance
ggplot(data.frame(rs=rstudent(men7),hat=hatvalues(men7),
                  cd=cooks.distance(men7)))+geom_point(shape=21)+
  aes(x=hat,y=rs,size=cd,color=clean2$REGION[-1])+xlab("hat value")+theme(legend.position="none")

#Plotting Menhinick Index at each region.
ggplot(clean2) + geom_violin(aes(x=REGION,y=MENHINICK,fill=REGION)) + theme_bw()

sims <- sim(men7,1000)
arm <- data.frame(
  fit=apply(fitted(sims,men7),1,function(x) quantile(x,0.5)),
  up=apply(fitted(sims,men7),1,function(x) quantile(x,0.975)),
  low=apply(fitted(sims,men7),1,function(x) quantile(x,0.025))
)


ggplot(data=arm,aes(x=as.numeric(rownames(arm)),y=fit,ymin=low,ymax=up)) + geom_point() + geom_linerange(alpha=0.5)


round(confint(men7),2)

ranef(men7)

## Confidence Intervals
hist(as.numeric(sims$sim_1),col="blue")
hist(clean2$MENHINICK,add=T,col="red")
plot(x=clean3$MENHINICK,y=sims$sim_1)

a <- fixef(men7)
c1 <- confint(men7)
b <- c1[-1:-2,]

d <- cbind(b,a)
dta <- data.frame(Value = a,
                 up= b[,1],
                 low=b[,2],
                 index=rownames(d),
                 Variables = Var)

Var = c("Intercept","Hard Coral", "Macroalgae", "Sand", "CCA", "Depth(z)", "Year", "Substrate Height")
dtn <- dta[-1,]

ggplot(dtn,aes(y=Variables,x=Value,xmin=low,xmax=up)) +
  geom_point() + geom_errorbarh()


## Model check graphics

plot(fitted(men7),resid(men7))

qqnorm(resid(men7))
qqline(resid(men7))
#binnedplot(fitted(men7),resid(men7,type="pearson"))


hist(x=as.numeric(simulate(men7)$sim_1))

ggplot() + geom_histogram(aes(x=simulate(men7)$sim_1),fill="blue",alpha=0.6) + geom_histogram(aes(x=clean3$MENHINICK),fill="red",alpha=0.6)


### Prediction
s_y <- sigma.hat(men7)$sigma$data
s_p <- sigma.hat(men7)$sigma$ISLAND
mu_i <- ranef(men7)$ISLAND

men7 <- lmer(MENHINICK ~I(HARD_CORAL/10) + I(MA/10) + I(SAND/10) + I(CCA/10) +
                        DEPTH.Z + YEARC + MAX_HEIGHT + 
                        (1|ISLAND),data=clean3)
stan_fit <- stan_lmer(MENHINICK ~I(HARD_CORAL/10) + I(MA/10) + I(SAND/10) + I(CCA/10) +
                   DEPTH.Z + YEARC + MAX_HEIGHT + 
                   (1|ISLAND),data=clean3)
stan_fit <- readRDS("stan_fit.rds")

ggplot(reshape2::melt(as.array(stan_fit)[,1,2:8])) +
  geom_density() + aes(x=value,color=parameters)

ggplot(reshape2::melt(as.array(stan_fit)[,1,9:15])) +
  geom_density() + aes(x=value,color=parameters,ncol=1)

# Graphics for data display
ggplot(clean3) + geom_histogram(aes(x=MENHINICK,fill=REGION),position="identity") + theme_bw() + facet_wrap(~REGION)

ggplot(clean3) + geom_histogram(aes(x=MENHINICK,fill=REGION),position="identity") + theme_bw() + facet_wrap(~ISLAND)



men7 <- lmer(MENHINICK ~ 
               I(HARD_CORAL/10) + I(MA/10) + I(SAND/10) + I(CCA/10) +
               DEPTH.Z + YEARC + MAX_HEIGHT + 
               (1|ISLAND),data=clean3)


length(unique(clean3$ISLAND))

hist(simulate(men7,use.u=T,newdata=clean3)$sim_1)
hist(clean3$MENHINICK,add=T)

sims <- simulate(men7,use.u=T,newdata=clean3)
simdata <- cbind(sims,clean3$REGION,clean3$ISLAND,clean3$MENHINICK)
colnames(simdata) <- c("sims", "REGION", "ISLAND","MENHINICK")

ggplot(simdata) + geom_density(aes(x=sims,color=REGION)) + geom_density(aes(x=MENHINICK,color="Measured Value"),linetype="dotted",color="darkblue") + theme_bw() + 
  facet_wrap(~factor(ISLAND,levels=names))

unique(factor(clean3$ISLAND,levels=names))
names <- unique(as.character(simdata$ISLAND))

ggplot(simdata) + geom_density(aes(x=sims,color=REGION)) + geom_density(aes(x=MENHINICK,color="Measured Value"),linetype="dotted",color="darkblue") + theme_bw() + facet_wrap(~ISLAND)
ggplot(simdata) + geom_density(aes(x=sims,color=REGION)) + geom_density(aes(x=MENHINICK,color="Measured Value"),linetype="dotted",color="darkblue") + theme_bw() + facet_wrap(~REGION)
ggplot(simdata) + geom_density(aes(x=sims),color="red") + geom_density(aes(x=MENHINICK),linetype="dashed",color="darkblue") + theme_bw()



ggplot(simdata) + geom_histogram(aes(x=MENHINICK,fill=REGION),alpha=0.5) + geom_histogram(aes(x=sims,fill=REGION),alpha=0.5) + theme_bw() + facet_wrap(~REGION)

ggplot(simdata) + geom_histogram(aes(x=sims,fill=REGION)) + geom_density(aes(x=sims,color=REGION),adjust=2) + theme_bw() + facet_wrap(~REGION)

ggplot(simdata) + geom_smooth(aes(x=sims,y=MENHENICK, fill=REGION))

ggplot(fortify(men7), aes(YEARC, MENHINICK,color=REGION)) +
  stat_summary(fun.data=mean_se,geom="pointrange") +
  stat_summary(aes(y=.fitted),fun.y=mean,geom="line")
  
years <- clean3 %>% group_by(OBS_YEAR,REGION,ISLAND) %>% summarize(n=n(),MENHINICK=mean(MENHINICK),RICHNESS=mean(RICHNESS),TOTALFISH=sum(TOTALFISH))

ggplot(years) + geom_bar(aes(x=OBS_YEAR,y=n,fill=REGION),stat="identity",width=0.8) + facet_wrap(~REGION)
ggplot(years) + geom_bar(aes(x=OBS_YEAR,y=MENHINICK,fill=REGION),stat="identity",width=0.8) + facet_wrap(~REGION)
ggplot(years) + geom_bar(aes(x=OBS_YEAR,y=RICHNESS,fill=REGION),stat="identity",width=0.8) + facet_wrap(~REGION)
ggplot(years) + geom_bar(aes(x=OBS_YEAR,y=TOTALFISH,fill=REGION),stat="identity",width=0.8) + facet_wrap(~REGION)

ggplot(clean3) + geom_histogram(aes(x=HARD_CORAL,fill=REGION)) + facet_wrap(~ISLAND)

hist(clean3$HARD_CORAL)

cleanpria <- na.omit(clean[clean$REGION=="PRIA",])
simp <- simulate(men7,newdata=cleanpria)
simpdataa <- cbind(simp,cleanpria$REGION,cleanpria$ISLAND,cleanpria$MENHINICK)
colnames(simpdataa) <- c("sims", "REGION", "ISLAND","MENHINICK")

ggplot(simpdataa) + geom_density(aes(x=sims,color=REGION)) + geom_density(aes(x=MENHINICK,color="Measured Value"),linetype="dashed",color="darkblue") + theme_bw() + facet_wrap(~ISLAND)

samples <- clean3 %>% group_by(ISLAND) %>% summarize(Sites=n())

ggplot(ranef(men7)) + geom_pointrange(aes())


posterior <- posterior_predict(stan_fit)

ppc_dens(y=stan_fit$y, yrep=posterior)

pp_check(stan_fit$y,posterior,fun="ppc_dens_overlay")

ppc_boxplot(y=stan_fit$y, yrep=posterior)

ggplot(fortify(ranef(men7)$'ISLAND:REGION')) + geom_pointrange(aes('Intercept'))

rr1 <- ranef(men7)
dd <- as.data.frame(rr1)


ggplot(dd,aes(y=grp,x=condval)) +
  geom_point() +
  geom_errorbarh(aes(xmin=condval -2*condsd,xmax=condval+2*condsd),height=0) + ylab("Island") +
  xlab("Random Effect")

ppc_dens_overlay(stan_fit$y,posterior[1:1000,],group=factor(clean3$ISLAND,levels=names))

loo(men7)


