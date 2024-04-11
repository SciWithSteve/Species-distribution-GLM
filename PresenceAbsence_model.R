setwd("C:/Users/steph/Documents/PhD/Experiments/Temperature development and survival/Hybrid SDM")

library(dismo) #load packages
library(car)
library(effects)
library(raster)
library(rgdal)
library(leaflet)
library(terra)
library(ozmaps)
library(spotoroo)
library(ggplot2)
library(ggspatial)



###Load in presence/absence data
PA_data_original<-read.csv("Trunc_PA.csv") #read PA data file
PA_data <- PA_data_original #copy data, so we can keep original dataframe (for plotting later, etc)
PA_data
coordinates(PA_data)= ~ lon+ lat #define coordinates
crs(PA_data) <- '+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs' #define coordinates reference system
PA_data$PresentAbsent <- factor(PA_data$PresentAbsent, levels=unique(PA_data$PresentAbsent)) #make sure PresentAbsent variable is a factor

###Prepare stack with predictor variables
#Get a vector map of Vic
#ozmap_states
#plot(ozmap_states)
#vic_map <- ozmap_states %>% filter(NAME %in% c("Victoria"))
plot(vic_map)
#OR
ozmap()
oz_states <- ozmaps::ozmap_states
oz_states
vic_state<-oz_states[2,]
vic_state
plot(vic_state)

#set working directory for predictor layer data
setwd("C:/Users/steph/Documents/PhD/Experiments/Temperature development and survival/wc2.1_30s_bio")

r1<-c("bio1.tif","bio2.tif","bio3.tif","bio4.tif","bio5.tif","bio6.tif",
      "bio7.tif","bio8.tif","bio9.tif","bio10.tif","bio11.tif","bio12.tif",
      "bio13.tif","bio14.tif","bio15.tif","bio16.tif","bio17.tif","bio18.tif",
      "bio19.tif") #prepare dataframe for predictor variables
r1<- stack(paste(getwd(),"/subset/", r1, sep="")) #upload predictors into a raster stack
r1<- crop(x = r1, y = extent(vic_state)) #crop variables to vic only
plot(r1) #visualise variables
r1low <- aggregate(r1, 10) #decrease resolution, better fit for resolution of PA data.
plot(r1low) #plot again, resolution has decreased
predictors<-r1low #put into stack named predictors

names(predictors) #list names of layers in stack

###option for more predictor variables###
r2<-c("aus_for18.tif") #make dataframe
r2<- stack(paste(getwd(),"/Veg/", r2, sep="")) #load variables into stack
plot(r2)
r2<- crop(x = r2, y = extent(vic_state)) #crop to vic
r2 <- terra::resample(r2, r1) #resample r2 to the the values (extent, etc) of r1,
#our first suite of variables using terra package. Takes a while.
s <- stack(r1, r2) #put both stacks together
predictors<-s #put into stack named predictors




###Extract predictor variable values at trap locations
#We need the values of the predictors at each trap location to build our SDM
presvals <- raster::extract(predictors, PA_data) #extract the predictor variable values at the presence/absence locations for C trunc.
#Code uses raster package, but should update to terra package
presvals #check that we have values for each predictor at each of the trap locations

presvals_scaled<-presvals #duplicate so we can create a dataframe for the predictor variables rescaled (so all are equal/standardised)
presvals_scaled<-raster::scale(presvals_scaled) #scale predictors (means all =0, sds all =1)
presvals_scaled

presvals_scaled.means<-apply(presvals_scaled,2,mean)  #save means of the transformed predictors
presvals_scaled.sds<-apply(presvals_scaled,2,sd)  #save means of the transformed predictors
print(round(presvals_scaled.means,2)) #check all means are now 0
print(round(presvals_scaled.sds,2)) #check all sds are now 1

presvals.means<-apply(presvals,2,mean)  #save mean of the original predictors
presvals.sds<-apply(presvals,2,sd) #save sds of the original predictors
predictors_scaled<-predictors #create new raster stack for transformed predictor rasters
for (ii in 1:19){
  predictors_scaled[[ii]]<-(predictors_scaled[[ii]]-presvals.means[ii])/presvals.sds[ii]
} #a loop to replace all raster stack values with transformed values, so we can make predictions with the same transformation we used for the model
predictors
predictors_scaled #check to see if scaling worked--looks about right

Trunc_PA_pred <- cbind(PA_data, presvals_scaled) #add scaled predictor values to dataframe with lat long of traps
Trunc_PA_pred
write.table(Trunc_PA_pred,file= "Trunc_PA_pred.csv", append=FALSE, sep= ",", row.names = FALSE, col.names=TRUE) #write as .csv file
summary(Trunc_PA_pred)
table(Trunc_PA_pred$PresentAbsent) #makes a cross-tabulation i.e. we have 32 1s (presences), 48 0s (absences) for our trap data



###Narrow down predictors
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...) #set up table for correlation figure
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

#check correlations: e.g. narrow down so we don't have very correlated predictor variables
library(corrplot)
corrplot(cor(presvals), type="upper")
subset_presvals<-presvals #lets iteratively remove predictor variables
subset_presvals<-subset(presvals, select=-c(bio1, bio19))
subset_presvals
corrplot(cor(subset_presvals), type="upper")
subset_presvals<-subset(subset_presvals, select=-c(bio5, bio18))
subset_presvals
corrplot(cor(subset_presvals), type="upper")
subset_presvals<-subset(subset_presvals, select=-c(bio2, bio17))
subset_presvals
corrplot(cor(subset_presvals), type="upper")
subset_presvals<-subset(subset_presvals, select=-c(bio3, bio16))
subset_presvals
corrplot(cor(subset_presvals), type="upper")
subset_presvals<-subset(subset_presvals, select=-c(bio10, bio14))
subset_presvals
corrplot(cor(subset_presvals), type="upper")
subset_presvals<-subset(subset_presvals, select=-c(bio11, bio13))
subset_presvals
corrplot(cor(subset_presvals), type="upper")
subset_presvals<-subset(subset_presvals, select=-c(bio9, bio12))
subset_presvals
corrplot(cor(subset_presvals), type="upper") #down to 5 variables, much less overlap
subset_presvals<-subset(subset_presvals, select=-c(bio6))
subset_presvals
corrplot(cor(subset_presvals), type="upper")
subset_presvals<-subset(subset_presvals, select=-c(bio7))
subset_presvals
corrplot(cor(subset_presvals), type="upper") #3 predictors, perhaps too few.

#Check actual values of correlations with correlation plot from psych package
library("psych")
corPlot(subset_presvals, cex=1.2)



par(mfrow=c(2,3))  #to do 2x3 plots in one figure, look at how C. truncatus presence/absence shifts with these variables
plot(Trunc_PA_pred$bio4,Trunc_PA_pred$PresentAbsent,pch=20,ylab="Presence/absence",xlab="bio4")
plot(Trunc_PA_pred$bio8,Trunc_PA_pred$PresentAbsent,pch=20,ylab="Presence/absence",xlab="bio8")
plot(Trunc_PA_pred$bio15,Trunc_PA_pred$PresentAbsent,pch=20,ylab="Presence/absence",xlab="bio15")
plot(Trunc_PA_pred$bio6,Trunc_PA_pred$PresentAbsent,pch=20,ylab="Presence/absence",xlab="bio6")
plot(Trunc_PA_pred$bio7,Trunc_PA_pred$PresentAbsent,pch=20,ylab="Presence/absence",xlab="bio7")

par(mfrow=c(2,3))  #to do 1 plot in the figure. Want to visualise how PA data plots against predictor
plot(predictors,'bio4') #plot a predictor
points(PA_data_original[PA_data_original[,3]==0,1:2], col="black", pch=3,  cex=0.8)  #add in points for presence locations
points(PA_data_original[PA_data_original[,3]==1,1:2], col="red",   pch=19, cex=0.8) #add in points for absence locations
legend('topleft', c('presence', 'absence'), bty='n', col=c('red', 'black'), pch=c(19, 3)) #add legend
plot(predictors,'bio8')
points(PA_data_original[PA_data_original[,3]==0,1:2], col="black", pch=3,  cex=0.8)  #add in points for presence locations
points(PA_data_original[PA_data_original[,3]==1,1:2], col="red",   pch=19, cex=0.8) #add in points for absence locations
legend('topleft', c('presence', 'absence'), bty='n', col=c('red', 'black'), pch=c(19, 3)) #add legend
plot(predictors,'bio15') #plot a predictor
points(PA_data_original[PA_data_original[,3]==0,1:2], col="black", pch=3,  cex=0.8)  #add in points for presence locations
points(PA_data_original[PA_data_original[,3]==1,1:2], col="red",   pch=19, cex=0.8) #add in points for absence locations
legend('topleft', c('presence', 'absence'), bty='n', col=c('red', 'black'), pch=c(19, 3)) #add legend
plot(predictors,'bio6')
points(PA_data_original[PA_data_original[,3]==0,1:2], col="black", pch=3,  cex=0.8)  #add in points for presence locations
points(PA_data_original[PA_data_original[,3]==1,1:2], col="red",   pch=19, cex=0.8) #add in points for absence locations
legend('topleft', c('presence', 'absence'), bty='n', col=c('red', 'black'), pch=c(19, 3)) #add legend
plot(predictors,'bio7') #plot a predictor
points(PA_data_original[PA_data_original[,3]==0,1:2], col="black", pch=3,  cex=0.8)  #add in points for presence locations
points(PA_data_original[PA_data_original[,3]==1,1:2], col="red",   pch=19, cex=0.8) #add in points for absence locations
legend('topleft', c('presence', 'absence'), bty='n', col=c('red', 'black'), pch=c(19, 3)) #add legend





###Explore models
#fit some toy models first with just 1 predictor
mod0 <- glm(PresentAbsent ~ 1, family = binomial(link="logit"), data = Trunc_PA_pred ) #null model, includes only intercept
mod1 <- glm(PresentAbsent ~ bio4, family = binomial(link="logit"), data = Trunc_PA_pred, trace=T) #single predictor
mod2 <- glm(PresentAbsent ~ bio8, family = binomial(link="logit"), data = Trunc_PA_pred, trace=T)
mod3 <- glm(PresentAbsent ~ bio15, family = binomial(link="logit"), data = Trunc_PA_pred, trace=T)
mod4 <- glm(PresentAbsent ~ bio6, family = binomial(link="logit"), data = Trunc_PA_pred, trace=T)
mod5 <- glm(PresentAbsent ~ bio7, family = binomial(link="logit"), data = Trunc_PA_pred, trace=T)
mod6 <- glm(PresentAbsent ~ bio4 + I(bio4^2), family = binomial(link="logit"), data = Trunc_PA_pred, trace=T) #try a quadratic predictor
mod7 <- glm(PresentAbsent ~ bio4+bio8+bio15+bio6+bio7, family = binomial(link="logit"), data = Trunc_PA_pred, trace=T) #multiple predictors
mod8 <- glm(PresentAbsent ~ bio4+bio8+I(bio8^2)+bio15+bio6+bio7, family = binomial(link="logit"), data = Trunc_PA_pred, trace=T) #multiple predictors

AIC(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8) #Summary of AIC values--lower is better
#mod 7 is lowest, as expected. Toyed with quatratic predictor in mod8, but seemed to raise AIC
#(we might expect temp variables to have an optimum, before dropping again, so good to check quadratic)

summary(mod7) #summary, is it better than null? Also prints AIC, etc.
library("performance")
r2(mod7)
mod7$fitted.values #check what the fitted values are at the trap sites
install.packages("relimp")
library("tornado")
imp <- importance(mod7, mod0)
plot(imp)
help(plot.importance_plot)


head(predict(mod7,type="response")) #response returns the fitted values on a probability scale (which is what we want to visualise)

mod7_pred <- predict(predictors_scaled, mod7, type="response")
par(mfrow=c(1,1))  #to do 1 plot in the figure
plot(mod7_pred, main="Predictions from mod7",zlim=c(0,1)) #plot prediction of a model
points(PA_data_original[PA_data_original[,3]==0,1:2], col="black", pch=3,  cex=0.8)  #add in points for presence locations
points(PA_data_original[PA_data_original[,3]==1,1:2], col="red",   pch=19, cex=0.8) #add in points for absence locations
legend('topleft', c('presence', 'absence'), bty='n', col=c('red', 'black'), pch=c(19, 3)) #add legend

#plot multiple models to compare visually
mod2_pred <- predict(predictors_scaled, mod2, type="response") #mod2 (predictor bio8, significant)
mod4_pred <- predict(predictors_scaled, mod4, type="response") #mod4 (predictor bio6, significant)
par(mfrow=c(1,3))  #to do 1 plot in the figure
plot(mod2_pred, main="Predictions from mod1",zlim=c(0,1))
plot(mod4_pred, main="Predictions from mod9",zlim=c(0,1))
plot(mod7_pred, main="Predictions from mod9",zlim=c(0,1))
AIC(mod2, mod4, mod7)

mod7_output<- rast(mod7_pred) #create a raster of the prediction output for mod7
mod7_output #check what it looks like. Scale about 0-1, as expected
par(mfrow=c(1,1))
plot(mod7_output) #check visually


###Initial model evaluation
library(car)
vif(mod7) #vif above 5-10 mean there are correlations between variables
#However, doesn't negatively affect model predictionss, we can see it looks to be performing, no ostensible overfitting

mod7.se<-summary(mod7)$coefficients[, 2]  #extract SEs from model object
mod7.se

Pres <- 1-predict(mod7, newdata=Trunc_PA_pred[Trunc_PA_pred$PresentAbsent==1,2:19],  type="response")
Abs  <- 1-predict(mod7, newdata=Trunc_PA_pred[Trunc_PA_pred$PresentAbsent==0,2:19],  type="response")

eval.mod7 <- evaluate(p=Pres, a=Abs)
eval.mod7 

str(eval.mod7)
eval.mod7@auc #equals 0.123, which means the classifier is predicting much worse than random (AUC=0.5)

threshold(eval.mod7)
par(mfrow=c(1, 2))

density(eval.mod7)
boxplot(eval.mod7)

cor.test(c(Pres, Abs), c(rep(1,length(Pres)),rep(0,length(Abs))))$estimate
eval.mod7@cor  # this should give the same value as above

dev.off()
par(mfrow=c(1, 2))
plot(eval.mod7, "ROC")
boxplot(eval.mod7)


###Map###

setwd("C:/Users/steph/Documents/PhD/Experiments/Temperature development and survival/Hybrid SDM")
library("tidyterra") #package to be able to use Spatrasters with tidyverse (ggplot)
library("ggforce")

basemap_vect <- ozmap_states %>% 
  vect()
basemap_vect #turn oz_states to vector

Trunc<-read.csv("Trunc_numbers.csv")
TruncIndividual<-read.csv("Presence_data.csv")

TruncPresent<- subset(Trunc, truncatus >= 1)
TruncPresent

TruncAbsent<- subset(Trunc, truncatus <1)
TruncAbsent


Output_map <-
  ggplot()+
  geom_spatraster(data = mod7_output) + #our model output
  geom_sf(data = Contour5, col="white", linetype=2) +
  geom_sf(data = basemap_vect, # vic state lines on top of raster
          fill = NA,
          col = "white",
          lwd = 0.75)+
  coord_sf(xlim = c(141.2, 149.8), ylim = c(-39.2, -34.12))  +
  
  #Gray border, titles etc
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "lightgrey",
                                    colour = "lightgrey",
                                    linewidth = 0.5, linetype = "solid"),
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                    colour = "white"), 
    panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                    colour = "white")) +

  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Carpophilus truncatus GLM projection", subtitle = "Victoria") +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  
  #C. truncatus presence dots
  geom_point(aes(x = Lon, y = Lat,
                 size = 1, stroke=0),
             colour = "#bd2f2f",  alpha = 0.5, show.legend=FALSE,
             data = TruncPresent) +
  geom_label(label="Present", colour = "#bd2f2f", x=152, y=-38.5)+
  annotate("rect",xmin=148,xmax=150, ymin=-39, ymax=-38.3, fill="white", alpha=0.5)+ #background for legend label
  annotate("text", colour = "#bd2f2f", x=149, y=-38.5, label= " • Present", fontface="bold") + #legend
  annotate("text", colour = "white", x=142.9, y=-34.4, label= "Mildura") + 
  annotate("text", colour = "white", x=145.3, y=-36.3, label= "Shepparton") + 
  annotate("text", colour = "white", x=144.8, y=-37.6, label= "Melbourne") + 
  
  #C. truncatus absent crosses
  geom_point(aes(x = Lon, y = Lat,
                 stroke=1.5, shape=4, alpha= 0.5),
             colour= "black", show.legend=FALSE,
             data = TruncAbsent) + scale_shape_identity() +
  annotate("text", colour = "black", x=149, y=-38.8, label= "x Absent", fontface="bold") +
  
  geom_circle(aes(x0 = 142.88, y0 = -34.7491, r = 1), colour="black", linetype=4, show.legend = FALSE, inherit.aes = FALSE)+

  labs(data=mod7_output, fill = "Probability\nof\npresence")+
  scale_fill_continuous(type="viridis", breaks = c(min(0.06),max(0.94)),
                        labels =c("High", "Low"),
                        guide=guide_colorbar(reverse = TRUE, ticks=FALSE))+
  labs(data=TruncPresent)
  

Output_map

