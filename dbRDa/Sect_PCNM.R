setwd("C:/Users/iarnelas/Documents/VALERIANELLA_R")
coordSect<-read.table("coordSect.txt", header=TRUE, fill=TRUE)
library (vegan)
pcnm1 <- pcnm(dist(coordSect))
pcnm1
scores(pcnm1)
genfst<-read.table("distSectpb.txt", header=TRUE, fill=TRUE)
genfst.dist<-as.dist(genfst)
genfst.dist
dataSectM <-read.table ("morphoSect.txt", header = TRUE)
dataSectM
#Marginal test todas variables
Fstfull<-capscale(genfst.dist~scores(pcnm1) + dataSectM$Fdiv + dataSectM$Actyp + dataSectM$AcCo + dataSectM$AcSpT + dataSectM$ShAc + dataSectM$CxTy, data = dataSectM, add =TRUE)
Fstfull
anova.cca(Fstfull)
anova.cca(Fstfull, by = "axis", permu= 9999)
anova.cca(Fstfull, by = "terms", permu= 9999)
RsquareAdj(capscale(genfst.dist~dataSectM$AcCo, data =dataSectM, add=TRUE))#se escoge significativa y mayor valor F
plot(Fstfull)
scores = scores(Fstfull)
geoscores = cbind(scores$sites, scores(pcnm1) + dataSectM$Fdiv + dataSectM$Actyp + dataSectM$AcCo + dataSectM$AcSpT + dataSectM$ShAc + dataSectM$CxTy)
cor(geoscores)
#Marginal solo geo
Fstfull<-capscale(genfst.dist~scores(pcnm1))
Fstfull
anova.cca(Fstfull)
anova.cca(Fstfull, by = "axis", permu= 9999)
anova.cca(Fstfull, by = "terms", permu= 9999)
RsquareAdj(capscale(genfst.dist~scores(pcnm1)))
plot(Fstfull)
scores = scores(Fstfull)
geoscores = cbind(scores$sites, scores(pcnm1))
cor(geoscores)
