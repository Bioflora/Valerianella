setwd("C:/Users/iarnelas/Documents/VALERIANELLA_R")
library(ade4)
library(vegan)
library (gower)
library(ape)
library(FD)
#calculo matriz dissimilarity Gower v. morfologicas
MorphSect <-read.table ("morphoSect.txt", header = TRUE)
MorphSect
MorphSect.dis <-gowdis(MorphSect)
MorphSect.dis
#PcoA morfología con dist Gower
pcoaSect<-pcoa(MorphSect.dis, correction="lingoes")
eigenvaluesSect<-pcoaSect$values
eigenvaluesSect
eigenvectorSect<-pcoaSect$vectors
eigenvectorSect
biplot(pcoaSect)
#distancias morfologicas usando eigenvector PCoA
eigenvM.dist<-dist(eigenvectorSect, method="euclidean", diag=FALSE, upper=FALSE)
eigenvM.dist
#PCA con distancias genéticas Fst
genfst<-read.table("distSectpb.txt", header=TRUE, fill=TRUE)
genfst.dist<-as.dist(genfst)
genfst.dist
pcagen<-prcomp(genfst.dist, scale=TRUE)
summary(pcagen)
pcagen
#extraer PCA eigenvector values 3 primeros ejes
scores(pcagen)
scores(pcagen, choices = 1:3)

##Variable morfológica con matriz distancias eigenvectors PcoA
#Marginal test solo genetica PCA, genetica significativo buen modelo
Fstfull<-capscale(eigenvM.dist~scores(pcagen, choices = 1:3))
Fstfull
anova.cca(Fstfull)
anova.cca(Fstfull, by = "axis", permu= 9999)
anova.cca(Fstfull, by = "terms", permu= 9999)
RsquareAdj(capscale(eigenvM.dist~scores(pcagen, choices = 1:3)))
plot(Fstfull)
scores = scores(Fstfull)
scores

##Variable morfológica con matriz distancias Gower
#Marginal test solo genetica PCA
Fstfull<-capscale(MorphSect.dis~scores(pcagen, choices = 1:3))
Fstfull
anova.cca(Fstfull)
anova.cca(Fstfull, by = "axis", permu= 9999)
anova.cca(Fstfull, by = "terms", permu= 9999)
RsquareAdj(capscale(MorphSect.dis~scores(pcagen, choices = 1:3)))
plot(Fstfull)
scores = scores(Fstfull)
scores

