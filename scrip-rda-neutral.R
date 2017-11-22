### Remove last features
rm(list=ls())

### Download libraries
library(vegan)
library(reshape2)
#library(rdaTest)
library(FactoMineR)

### Download data
data <- read.table("9770snps-neutral.frq.strat", header=T, sep="\t", dec=".")
#colnames(data)=c("SNP","POP","A1","A2","MAF","MAC","NCHR")
data2 <- data[,2:3]
data3 <- cbind(data2, data$MAF)
colnames(data3)=c("SNP","POP","MAF")

### Convert from long to wide
data4 <- dcast(data3, POP~SNP)
write.table(data4, "9770snps-all.txt", quote=F, sep="\t")

### Download geno, env and geo
# Genotypes object
geno <- data4

# Remove first colum (pop)
Geno=geno[,2:ncol(geno)]

# Environmental object
env=read.table("Environmental.txt", header=TRUE)

# Remove first colum (ENV)
Env=env[,2:8]

# geo object
geo=read.table("Geographical_data.txt", header=TRUE)
GEO=geo[,2:ncol(geo)]

# Do environmental matrix (PCA accept missing value)
Env_pca=PCA(Env, scale.unit=FALSE, ncp=10)
Env_pca$eig

### Keep the eigenvalue > 1 and indicate the number of eigenvalue to keep
PCeigen_1=3
ind=Env_pca$ind
load=Env_pca$var # check the loading of each environmental parameter
Env_load=load$cor[,1:PCeigen_1]
ENV=ind$coord[,1:PCeigen_1]

### Do geno matrix 
Geno_pca=prcomp(Geno)
summary(Geno_pca)

### Keep number of axes which explained percent of variance that explained more than 5% (or the percent of cumulative variance explained more than 90%).
# here we kept the percent of variance that explained more than 5%
Geno_rda=Geno_pca$x[,1:9]

### Do a RDA with PCA factors
rda_homard=rda(Geno_rda~ENV[,1]+ENV[,2]+ENV[,3], scale=T)
anova(rda_homard,by="margin",step=1000)

### Obtain the Pvalue for the relationship
anova(rda_homard, step=1000)

### Look the Rsquared
RsquareAdj(rda_homard)
# adj(r.squared) = 0.05
sommaire = summary(rda_homard)

### Make the orbistep to select the variables
ordistep(rda_homard,perm.max=200)
mod1 <- rda(Geno_rda ~ ., Env)
ordistep(mod1,perm.max=200)

### Do a RDA with environmental variables
rda_homard=rda(Geno_rda~Env$MEAN_WINTER+Env$MAX_YEAR, scale=T)
anova(rda_homard,by="margin",step=1000)
anova(rda_homard,step=1000)

### Look the Rsquared
RsquareAdj(rda_homard)
# adj(r.squared) = 0.07
sommaire = summary(rda_homard)

### Look the proportion explained in the "Importance of Component"
sommaire

# GRAPHIQUE RDA1*RDA2;
# -------------------;
axes=c(1,2);
R2 = c("9.4%","8.3%")
marqueur = sommaire$species[,axes]
envt = sommaire$biplot[,axes]
objet = sommaire$constraints[,axes]
row.names(objet) = c("ANT","BON","BOO","BRA","BRO","BUZ","CAN","MAG","CAR","GAS","LOB","MAL","MAR","OFF","RHO","SEA","SID","SJI","TRI")

#calcule la qualite graphique des especes (representativite)
qualite12 = marqueur[,1]^2 + marqueur[,2]^2
marqueur2 = marqueur[qualite12>0.05,]

#graphe des axes et des marqueurs 
par(new=F)
yrange1 = round(max(c(abs(marqueur[,2]),abs(objet[,2]))),digits=2)
xrange1 = round(max(c(abs(marqueur[,1]),abs(objet[,1]))),digits=2)
plot(marqueur2,type="p",xlab=paste("RDA",axes[1]," (",R2[1],")",sep=""),
     ylab=paste("RDA",axes[2]," (",R2[2],")",sep=""),ylim=c(-yrange1,yrange1),
     xlim=c(-xrange1,xrange1),pch="+",cex=0)
abline(h=0,lty=3)	

# cross-lines
abline(v=0,lty=3)
points(objet, pch=21,bg=c("black","black","white","black","white","white","black","black","black","black","white","black","white","white","white","white","black","black","black"))
text(objet, labels=row.names(objet), cex=0.8, pos=4)

### Make a graph 
par(new=T)
yrange2 = round(max(abs(envt[,2])),digits=2)
plot(envt,type="n",axes=F,ylab="",xlab="",ylim=c(-yrange2,yrange2),xlim=c(-1,1),pch=0,cex=0.75)
axis(4, col="black",col.ticks="black",col.axis="black"); axis(3,col="black",col.ticks="black",col.axis="black")

#b) fleches;
x1 <- envt[,1]*0.95;
y1 <- envt[,2]*0.95;
points(x1,y1,pch="");
arrows(x0=rep(0,12),y0=rep(0,12),x1=x1, y1=y1,length=0.05, col="black")

### Look which axe explains more the genetic variation
x1

#c) Identification des facteurs;
text(x1[1]+0.5,y1[1]-0.05,"MEAN WINTER (+0.56)", cex=0.75,pos=2,font=2)
text(x1[1]-0.7,y1[1]-0.15,"MAX YEAR (-0.46)", cex=0.75,pos=2,font=2)

#ajouter p-value et R2
R=expression(paste("Adj.R"^"2"))
text(0.4,0.8,R,pos=4)
text(0.7,0.8,"= 0.07",pos=4)
text(0.4,0.7,"P-value = 0.004",pos=4)

#################### PARTIAL RDA #######################

### Do a partial RDA with PC factors
rda_homard=rda(Geno_rda, ENV, GEO)

### Obtain the Pvalue for the relationship
anova(rda_homard, step=1000)
# pvalue = 0.05

### Look the Rsquared
RsquareAdj(rda_homard)

# adj(r.squared) = 0.02
sommaire = summary(rda_homard)

### Look the proportion explained in the "Importance of Component"
sommaire

### Do a RDA with environmental variables
rda_partial_homard=rda(Geno_rda~Env$MEAN_WINTER+Env$MAX_YEAR, GEO)
anova(rda_partial_homard,by="margin",step=1000)
anova(rda_partial_homard,step=1000)

# GRAPHIQUE RDA1*RDA2;
# -------------------;
axes=c(1,2);
R2 = c("14.3%","9.5%")
marqueur = sommaire$species[,axes]
envt = sommaire$biplot[,axes]
objet = sommaire$constraints[,axes]
row.names(objet) = c("ANT","BON","BOO","BRA","BRO","BUZ","CAN","MAG","CAR","GAS","LOB","MAL","MAR","OFF","RHO","SEA","SID","SJI","TRI")

#calcule la qualite graphique des especes (representativite)
qualite12 = marqueur[,1]^2 + marqueur[,2]^2
marqueur2 = marqueur[qualite12>0.05,]

# Make the graph with the axes and markers
par(new=F)
yrange1 = round(max(c(abs(marqueur[,2]),abs(objet[,2]))),digits=2)
xrange1 = round(max(c(abs(marqueur[,1]),abs(objet[,1]))),digits=2)
plot(marqueur2,type="p",xlab=paste("RDA",axes[1]," (",R2[1],")",sep=""),
     ylab=paste("RDA",axes[2]," (",R2[2],")",sep=""),ylim=c(-yrange1,yrange1),
     xlim=c(-xrange1,xrange1),pch="+",cex=0)
abline(h=0,lty=3) 

# Cross-lines
abline(v=0,lty=3)
points(objet, pch=21,bg=c("black","black","white","black","white","white","black","black","black","black","white","black","white","white","white","white","black","black","black"))
text(objet, labels=row.names(objet), cex=0.8, pos=4)


# Make a graph with the variables
#a) points;
par(new=T)
yrange2 = round(max(abs(envt[,2])),digits=2)
plot(envt,type="n",axes=F,ylab="",xlab="",ylim=c(-yrange2,yrange2),xlim=c(-1,1),pch=0,cex=0.75)
axis(4, col="black",col.ticks="black",col.axis="black"); axis(3,col="black",col.ticks="black",col.axis="black")

#b) fleches;
x1 <- envt[,1]*0.95;
y1 <- envt[,2]*0.95;
points(x1,y1,pch="");
arrows(x0=rep(0,12),y0=rep(0,12),x1=x1, y1=y1,length=0.05, col="black")

### Look which axe explains more the genetic variation
x1
### Then look which variables is related to each axe
Env_load

#c) Identification des facteurs;
text(x1[1],y1[1],"MIN YEAR (0.84)", cex=0.75,pos=2,font=2)
text(x1[1],y1[1]-0.1,"MIN SUMMER (0.80)", cex=0.75,pos=2,font=2)

#ajouter p-value et R2
R=expression(paste("adj.R"^"2"))
text(-0.5,1.0,R,pos=4)
text(-0.3,1.0,"= 0.09",pos=4)
text(-1.1,0.63,"P-value = 0.045",pos=4)

