# Codes for the paper : "Population genomics in Sebastes spp."

[![https://upload.wikimedia.org/wikipedia/en/d/d6/DFO_Logo.png](https://upload.wikimedia.org/wikipedia/en/d/d6/DFO_Logo.png)](DFO)

_______________________________________________________________________________



#### Laura Benestan, Caroline Senay, Geneviève Parent

Rimouski, 2017-2020

Submited to Evolutionary Applications, 2020


_______________________________________________________________________________

Using multivariate analyses is not an easy task and finding the proper way to perform such analysis can be sometime trivial.
Here, we indicated how we used a db-RDA analysis following a step by step tutorial.
You will see during thsi tutorial how the tool of distance-based redundancy analysis can nicely demonstrate the influence of depth on teh genomic variation of Sebastes spp.

### Softwares

- [R Version 3.6.0](https://cran.r-project.org/)
	* R packages: codep, adespatial, adegraphics, vegan, ape, car, adegenet, dplyr

### Data import

Data need to be in genepop format. 
You can easly find functions that can transform your dataset into a genpop format.
For instance the `genomic_converter` function available in the elegant package [assigner](https://github.com/thierrygosselin/assigner)

Import the Genepop file into the R environment:
```{r}
genpop_sebastes <- read.genepop("24603snps_860ind_sebastes.gen",ncode = 3L) # for both species
```

### Calculate euclidian distances on a genepop object 

First, estimating individual genetic distances is a crucial step before performing the db-RDA.
The individual genetic distances will be considered as the **response variables**.
```{r}
distgenEUCL <- dist(genpop_sebastes, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
```

Second, import the envrionmental dataset, which will contained the **explanatory variables**
```{r}
Env <- read.table("Env_860ind.txt", header=TRUE)
```
[![https://fr.freepik.com/icones-gratuites/attention-signe_770843.htm](Attention)

Check that both datasets - response and explanatory variables - are in the same order.
Order <- read.table("Order.txt", header=FALSE)
Order$IND <- substr(Order$V1, 7,12)
Order$SPECIES <- substr(Order$V1,1,3)
Order$SPECIES_NUMBER <- ifelse(Order$SPECIES=="fas",1,2)

# Merge Order and Env to be in the right Order
Env_order <- merge(Order, Env, by.x="IND", by.y="label_bioinfo",sort=F)

#Recreate a Env file
Env <-select(Env_order,-V1)

######################
#Create MEM- Spatial variable
######################
#Look ast sites in space
#Keep just lat and long in Coor
Coor=Env[,5:6]
#look at sites spatial distribution 
plot(Coor, asp=1) 

###Compute spatial distances among sites accounting for the earth curvature
DistSpatial=gcd.hf(Env[,5:6]) 

#Compute MEM
# Keeping default setting for truncation (length of the longest edge of the minimum spanning tree will be used as the threshold) and just positive MEM
dbmem = dbmem(DistSpatial) 
#Look at general output
summary(dbmem)
#With this treashold, we get 16 MEM

#Visualising the links longer than the treashold
s.label(Coor, nb = attr(dbmem, "listw"))

#Visualising the 16 mem. Can create more spatial structure type when MEM are selected together in the same analysis
#The 1rst MEMs are large spatial scales, the last MEMs are small spatial scales
s.value(Coor, dbmem[,1:16])
#MEM can be used in stats like any other env variables

#db-Rda Response variables

##Pcoa on genetic distance matrix. Genetic distances are then in the multivariate space format, adapted to db-RDA
Pcoa=pcoa(distgenEUCL)
Pcoa

#"There were no negative eigenvalues. No correction was applied"   
#If negative or null eigenvalues are produced, they need to be expluded
#The first axe explains 499% of the variation. We will keep them all to analyse all genetic variation

#Extract Pcoa principal components, which will be the response variable in the db-RDA
X=Pcoa$vectors
#Look at genotypes distribution in relation to the first 2 Pcoa axes
plot(X[,1], X[,2])
#Clusters are clearly appearent

#explanatory variables
#Create a matrix with all expanatory variables, 16 MEMs, depth, and n-1 years coded in dummy variables (2001, 2002, 2008, 2013, 2014)  
Y=cbind(dbmem, Env[,7:12])

#Correlation among explanatory variables
##Look at correlation among explanatory variables
cor(Y)
#Nothing problematic. Not surprising cause MEM are orthogonal to each other

#db-RDA global model with all explanatory variables
rda1=rda(X, Y)

#Looking at VIf for multicolinearity within the model. Greater than around 10 is problematic
vif(rda1)
#We don't need to exclude variables because of correlation

#Explained variance by the global db-RDA model. Adjusted R2 accounts for the number of variables 
RsquareAdj(rda1)
#37% explained variance

#db-RDA Global model probability
anova(rda1, perm=999)
#p=0.001 

#Variable selection with OrdiR2Step. 
#This function allows to add and remove variable to maximise the explained variance. 
#To avoid overfitting, selected variables should not explained more than the global model (36%) 

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0<-rda(X ~ 1, Y)
#OrdiR2step will move towards the global model with all explanatory variables
rdaG<- rda(X ~ ., Y)

#Variables selection
#Variables will be selected until the 36% (rda global model) is reached
Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="both") 
#  X ~ MEM2 + MEM1 + Y2008 + Depth_m + Y2013 + Y2014 + Y2002 + Y2001
#summary table with selected variables               
Sel$anova

#Spatial scale of the MEM2 the 1rst variables that we get
s.value(Coor, dbmem[,2])
#Spatial scale of the MEM3
s.value(Coor, dbmem[,3])
#All MEMS
s.value(Coor, dbmem[,1:16])

#Model with only selected variables
Ysel=cbind(Y$MEM3,Y$MEM2,Y$Y2014, Y$Y2013,Y$MEM4,Y$MEM10,Y$MEM6,Y$MEM9,Y$Depth_m,Y$MEM15,Y$MEM1,Y$MEM11,Y$MEM8,Y$Y2002,Y$MEM5,Y$Y2001,Y$MEM12,Y$MEM14,Y$MEM7,Y$MEM13)
rdaS<- rda(X ,Ysel)
summary(rdaS, scaling=1)    

# Check the RDA summary
RsquareAdj(rdaS)

### Load RDA on 860 ind
load("RDA_all_860ind.Rdata")

################# VISUALISE DB-RDA results for 860ind ##################
#db-RDA biplot
site=cbind(scores(rdaS, display="sites", choices=c(1,2), scaling=1), Env[,13])
species <- site[,3]
#Splitting the database by species to give them different colors in biplot
Fasciatus=site[site[,3]==1,]
Mentella=site[site[,3]==2,]

#Scaling 1 is used, where  distances among objects approximate their Euclidean distances in the space of response variables 
plot(rdaS, scaling=2, main="", type="none", xlab=c("db-RDA-1"), ylab=c("db-RDA-2"), xlim=c(-20, 20), ylim=c(-5, 5))
col2 <- c("blue","red")

#Add point with a different color for each species, some point are getting lost under the other ones, but not a big deal I think
#Can multiply both axis by a constant to rescale it nicely.
points(Fasciatus[,1]*3, Fasciatus[,2]*3, col="black",bg="red", pch=21, cex=1.2) 
points(Mentella[,1]*3, Mentella[,2]*3, col="black",bg="blue", pch=21, cex=1.2) 

#Add arrows showing selected variables pointing to their dbRDA scores
#Again can rescale to minimise overlapping of both data source. got to multiply both axis
arrows(0,0, scores(rdaS, display="bp", choices=1, scaling=1)*70, scores(rdaS, display="bp", choices=2, scaling=1)*70, col="black", length=0.1)
#Add variables names
text(scores(rdaS, display="bp", choices=1, scaling=1)*75, scores(rdaS, display="bp", choices=2, scaling=1)*75, labels=c("MEM3","MEM2","2014","2013","MEM4","MEM10","MEM6","MEM9","Depth","MEM15","MEM1","MEM11","MEM8","2002","MEM5","2001","MEM12","MEM7","MEM13"), col="black", cex=0.8, pos=3)

### Save image of RDA
save.image("RDA_all_860ind.Rdata")

####################################### FASCIATUS ############################################## FASCIATUS
#Clean up workspace
rm(list=ls())

#Import Genetic distance matrix
DistGenfas <- read.table('distgenEUCL-sebastes_444ind.txt')

### Extract the order of individuals in fasciatus
Order <- read.table("Orderfasciatus.txt", header=FALSE)
Order$IND <- substr(Order$V1, 7,12)

#Import explanatory variables
Env <- read.table("Env_860ind.txt", header=TRUE)

# Merge Order and Env to be in the right Order
Env_order <- merge(Order, Env, by.x="IND",by.y="label_bioinfo",sort=F)

#Recreate a Env file
Envfas <-select(Env_order,-V1)

#Look at sites in space
#Keep just lat and long in Coor
Coorfas=Envfas[,3:4]
#look at sites spatial distribution 
plot(Coorfas, asp=1) 

###Compute spatial distances among sites accounting for the earth curvature
DistSpatialfas=gcd.hf(Envfas[,3:4]) 

#Compute MEM
# Keeping default setting for truncation (length of the longest edge of the minimum spanning tree will be used as the threshold) and just positive MEM
dbmemfas = dbmem(DistSpatialfas) 
#Look at general output
summary(dbmemfas)
#With this threashold, we get 16 MEM

#Visualising the links longer than the treashold
s.label(Coorfas, nb = attr(dbmemfas, "listw"))

#Visualising the 16 mem. Can create more spatial structure type when MEM are selected together in the same analysis
#The 1rst MEMs are large spatial scales, the last MEMs are small spatial scales
s.value(Coorfas, dbmemfas[,1:9]) #MEM can be used in stats like any other env variables

##Pcoa on genetic distance matrix. Genetic distances are then in the multivariate space format, adapted to db-RDA
Pcoafas=pcoa(DistGenfas)
Pcoafas

#"There were no negative eigenvalues. No correction was applied"   
#If negative or null eigenvalues are produced, they need to be excluded

#Extract Pcoa principal components, which will be the response variable in the db-RDA
Xfas=Pcoafas$vectors
#Look at genotypes distribution in relation to the first 2 Pcoa axes
plot(Xfas[,1], Xfas[,2])
#Clusters are clearly appearent

#explanatory variables
#Create a matrix with all expanatory variables, here for fasciatus the 2008 variable was excluded as 0 was the value for all individual (no variation)
Yfas=cbind(dbmemfas, Envfas%>% select(Y2001,Y2002,Y2013,Y2014))

#Correlation among explanatory variables
##Look at correlation among explanatory variables
cor(Yfas)
#Nothing problematic. Not surprising cause MEM are orthogonal to each other

#db-RDA global model with all explanatory variables
rda1=rda(Xfas, Yfas)

#Looking at VIf for multicolinearity within the model. Greater than around 10 is problematic
vif(rda1)
#We don't need to exclude variables because of correlation

#Explained variance by the global db-RDA model. Adjusted R2 accounts for the number of variables 
RsquareAdj(rda1)
#3.71% explained variance

#db-RDA Global model probability
anova(rda1, perm=999)
#p=0.001 

#Variable selection with OrdiR2Step. 
#This function allows to add and remove variable to maximise the explained variance. 
#To avoid overfitting, selected variables should not explained more than the global model (36%) 

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0fas<-rda(Xfas ~ 1, Yfas)
#OrdiR2step will move towards the global model with all explanatory variables
rdaGfas<- rda(Xfas ~ ., Yfas)

#Variables selection
#Variables will be selected until the 36% (rda global model) is reached
Selfas <- ordiR2step(rda0fas, scope = formula(rdaGfas), direction="both")                
#summary table with selected variables               
Selfas$anova

#Spatial scale of the MEM4 the 2nd variables that we get
s.value(Coorfas, dbmemfas[,4])
#Spatial scale of the MEM2
s.value(Coorfas, dbmemfas[,2])
#All MEMS
s.value(Coorfas, dbmemfas[,1:9])

#Model with only selected variables
Yselfas=cbind(Yfas$Y2002,Yfas$MEM4,Yfas$Y2014,Yfas$Y2001,Yfas$Y2013,Yfas$MEM2,Yfas$MEM8,Yfas$MEM5,Yfas$MEM1,Yfas$YMEM3,Yfas$MEM7,Yfas$MEM6)
rdaSfas<- rda(Xfas ,Yselfas)
summary(rdaSfas, scaling=1)        

### Load RDA for fasciatus
load("RDA_all_444ind.Rdata")

#Put together the site scores and species identity
sitefas=cbind(scores(rdaSfas, display="sites", choices=c(1,2), scaling=1), Envfas[,11])

### Add cluster info
cluster <- read.table("Admixture_results_fas_K4_444ind.txt",header=TRUE)
cluster_vector <- cluster$CLUSTER
col5 <- c("purple","orange1","red1","chartreuse4","deeppink2")

#Add point with a different color for each cluster, some point are getting lost under the other ones, but not a big deal I think
plot(rdaSfas, type="n", xlim=c(-10, 10))
points(rdaSfas, display="sites", pch=21, cex=1.3, col="black", scaling=2, bg=col5[cluster_vector]) # the wolves
#Again can rescale to minimise overlapping of both data source. got to multiply both axis
arrows(0,0, scores(rdaSfas, display="bp", choices=1, scaling=1)*70, scores(rdaSfas, display="bp", choices=2, scaling=1)*70, col="black", length=0.1)
#Add variables names
text(scores(rdaSfas, display="bp", choices=1, scaling=1)*50, scores(rdaSfas, display="bp", choices=2, scaling=1)*50, labels=c("2002","MEM4","2014","2001","2013","MEM2","MEM8","MEM5","MEM1","MEM3","MEM7","MEM6"), col="black", cex=0.8, pos=3)
#legend("bottomright", legend=c("K1","K2","K3","K4","K5"), bty="n", col="gray32", pch=21, cex=1, pt.bg=col5[cluster_vector])

### Save image of RDA
save.image("RDA_all_444ind.Rdata")

####################################### MENTELLA ####################################
#Import Genetic distance matrix
DistGenmen <- read.table('distgenEUCL-sebastes_416ind.txt')

### Extract the order of individuals in fasciatus
Order <- read.table("Ordermentella.txt", header=FALSE)
Order$IND <- substr(Order$V1, 7,12)

#Import explanatory variables
Env <- read.table("Env_860ind.txt", header=TRUE)

# Merge Order and Env to be in the right Order
Env_order <- merge(Order, Env, by.x="IND",by.y="label_bioinfo",sort=F)

#Recreate a Env file
Envmen<-select(Env_order,-V1)

#Look at sites in space
#Keep just lat and long in Coor
Coormen=Envmen[,3:4]
#look at sites spatial distribution 
plot(Coormen, asp=1) 

###Compute spatial distances among sites accounting for the earth curvature
DistSpatialmen=gcd.hf(Envmen[,3:4]) 

#Compute MEM
# Keeping default setting for truncation (length of the longest edge of the minimum spanning tree will be used as the threshold) and just positive MEM
dbmemmen = dbmem(DistSpatialmen) 
#Look at general output
summary(dbmemmen)
#With this threashold, we get 16 MEM

#Visualising the links longer than the treashold
s.label(Coormen, nb = attr(dbmemmen, "listw"))

#Visualising the 16 mem. Can create more spatial structure type when MEM are selected together in the same analysis
#The 1rst MEMs are large spatial scales, the last MEMs are small spatial scales
s.value(Coormen, dbmemmen[,1:9]) #MEM can be used in stats like any other env variables

##Pcoa on genetic distance matrix. Genetic distances are then in the multivariate space format, adapted to db-RDA
Pcoamen=pcoa(DistGenmen)
Pcoamen

#"There were no negative eigenvalues. No correction was applied"   
#If negative or null eigenvalues are produced, they need to be excluded

#Extract Pcoa principal components, which will be the response variable in the db-RDA
Xmen=Pcoamen$vectors
#Look at genotypes distribution in relation to the first 2 Pcoa axes
plot(Xmen[,1], Xmen[,2])
#Clusters are clearly appearent

#explanatory variables
#Create a matrix with all expanatory variables, here for fasciatus the 2008 variable was excluded as 0 was the value for all individual (no variation)
Ymen=cbind(dbmemmen, Envmen%>% select(Depth_m,Y2008,Y2013,Y2014))

#Correlation among explanatory variables
##Look at correlation among explanatory variables
cor(Ymen)
#Nothing problematic. Not surprising cause MEM are orthogonal to each other

#db-RDA global model with all explanatory variables
rda2=rda(Xmen, Ymen)

#Looking at VIf for multicolinearity within the model. Greater than around 10 is problematic
vif(rda2)
#We don't need to exclude variables because of correlation

#Explained variance by the global db-RDA model. Adjusted R2 accounts for the number of variables 
RsquareAdj(rda2)
#8.7% explained variance

#db-RDA Global model probability
anova(rda2, perm=999)
#p=0.001 

#Variable selection with OrdiR2Step. 
#This function allows to add and remove variable to maximise the explained variance. 
#To avoid overfitting, selected variables should not explained more than the global model (36%) 

#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0men<-rda(Xmen ~ 1, Ymen)
#OrdiR2step will move towards the global model with all explanatory variables
rdaGmen<- rda(Xmen ~ ., Ymen)

#Variables selection
#Variables will be selected until the 36% (rda global model) is reached
Selmen <- ordiR2step(rda0men, scope = formula(rdaGmen), direction="both")                
#summary table with selected variables               
Selmen$anova

#Spatial scale of the MEM2 the 2nd variables that we get
s.value(Coormen, dbmemmen[,2])
#Spatial scale of the MEM1
s.value(Coormen, dbmemmen[,1])
#All MEMS
s.value(Coormen, dbmemmen[,1:9])

#Model with only selected variables
Yselmen=cbind(Ymen$Depth_m, Ymen$MEM2, Ymen$MEM1, Ymen$MEM3, Ymen$Y2008, Ymen$MEM7,Ymen$Y2013, Ymen$Y2014, Ymen$MEM10,  Ymen$MEM5, Ymen$MEM4,  Ymen$MEM8)
rdaSmen<- rda(Xmen ,Yselmen)
summary(rdaSmen, scaling=1)        

### Add ecotype information
cluster <- read.table("Admiwture_mentella_K4.txt",header=TRUE)
cluster_vector <- cluster$CLUSTER
col3 <- c('darkblue',"gray","blue",'cyan2')

### Load RDA for mentella
load("RDA_all_416ind.Rdata")

#Put together the site scores and species identity
sitemen=cbind(scores(rdaSmen, display="sites", choices=c(1,2), scaling=1), Envmen[,11])

#Add point with a different color for each species, some point are getting lost under the other ones, but not a big deal I think
plot(rdaSmen, type="n", xlim=c(-15, 20))
plot(rdaSmen, scaling=2, main="", type="none", xlab=c("db-RDA-1"), ylab=c("db-RDA-2"), xlim=c(-5, 5), ylim=c(-10, 10))
points(rdaSmen, display="sites", pch=21, cex=1.3, col="black", scaling=2, bg=col3[cluster_vector]) 

#Again can rescale to minimise overlapping of both data source. got to multiply both axis
arrows(0,0, scores(rdaSmen, display="bp", choices=1, scaling=1)*70, scores(rdaSmen, display="bp", choices=2, scaling=1)*70, col="black", length=0.1)
#Add variables names
text(scores(rdaSmen, display="bp", choices=1, scaling=1)*75, scores(rdaSmen, display="bp", choices=2, scaling=1)*75, labels=c("Depth","MEM2","MEM1","MEM3","2008","MEM7","2013","2014","MEM10","MEM5","MEM4","MEM8"), col="black", cex=0.8, pos=3)
legend("bottomright", legend="Mentella", bty="n", col="gray32", pch=21, cex=1, pt.bg="blue")

### Save image of RDA
save.image("RDA_all_416ind.Rdata")

