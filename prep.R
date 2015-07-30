rm(list=ls())
setwd("/Users/minjay/Downloads/PEF-data")
library(MBA)
library(fields)
library(rgdal)
library(ggplot2)
library(spBayes)
library(geoR)

PEF.shp <- readOGR("shapefiles","PEF-bounds")

# remove NAs
PEF.plots <- na.omit(read.csv("PEF-plots.csv"))
n.ref <- nrow(PEF.plots)
n.dist <- 151

PEF.LVIS <- read.csv("PEF-LVIS.csv")

LVIS <- rbind(as.matrix(PEF.plots[, 37:158]), as.matrix(PEF.LVIS[, 32:153]))

##PCA using prcomp
pca.test<-prcomp(LVIS, scale = T)
summary(pca.test)
biplot(pca.test)
plot(pca.test)
pc.score <- predict(pca.test)[1:451,1:2]
PEF.data <- cbind(PEF.plots[,1:7],pc.score)
PEF.predict <- cbind(PEF.LVIS[,1:2], predict(pca.test)[452:12865,1:2])
PEF.coords <- PEF.data[,3:4]
Lid.coords <- as.matrix(PEF.LVIS[,1:2])

lm1 <- lm(biomass.mg.ha ~ PC1 + PC2, data = PEF.data)
vg1 <- variog(coords = PEF.coords, data = resid(lm1))
plot(vg1)
spline.obj <- spline(vg1$u, vg1$v, xmin=0, xmax=max(vg1$u))
lines(spline.obj$x, spline.obj$y, col="red")

max.dist <- max(iDist(PEF.coords))

starting <- list("phi" = 3/1000, "sigma.sq" = 1, "tau.sq" = 1)
tuning <- list("phi" = 0.1, "sigma.sq" = 0.1, "tau.sq" = 0.1)
priors.1 <- list("beta.Norm"=list(rep(0,3), diag(1000,3)),
                 "phi.Unif"=c(3/3000, 3/1), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))

cov.model <- "exponential"

n.samples = 2000
n.report = 500
verbose = TRUE

m.1 <- spLM(y3~T.ref[, 1]+T.ref[, 2], coords=PEF.coords, starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)

burn.in <- 0.5*n.samples
m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)
coords <- cbind(PEF.LVIS[, 1], PEF.LVIS[, 2])
m.1.pred <- spPredict(m.1, pred.covars=cbind(rep(1, dim(PEF.LVIS)[1]), T[-c(1:n.ref), 1:2]), pred.coords=coords,
                      start=0.5*n.samples)
y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, mean)

quilt.plot(coords, y.hat)







ggplot(data = PEF.plots, aes(x = x.coords, y = y.coords, label = MU)) +
  geom_point() + geom_text(aes(label = MU), hjust=0, vjust=0)

PEF.plots.long <- reshape(PEF.plots, varying = colnames(PEF.plots)[8:158], v.names = "energy", 
                          timevar = "height", times = 1:151, direction = "long")

PEF.plots.long.sort <- PEF.plots.long[order(PEF.plots.long$id),]

head(PEF.plots.long.sort)

selected <- subset(PEF.plots.long.sort, PEF.plots.long.sort$MU %in% c("10", " 7B"))
p <- ggplot(data = selected, aes(x = height, y = energy, group = id, colour = MU))
p + geom_line()

##Plot of all curves
p <- ggplot(data = PEF.plots.long.sort, aes(x = height, y = energy, group = id, colour = MU))
p + geom_line()

coords <- PEF.plots[,c("x.coords", "y.coords")]
stems.ha <- PEF.plots$stems.ha

##get bounding box of PEF
b.box <- as.vector(t(bbox(PEF.shp)))

surf <- mba.surf(cbind(coords, stems.ha), no.X=200, no.Y=200, extend=TRUE, sp=TRUE, b.box=b.box)$xyz.est

##set the projection for the resulting sp object
proj4string(surf) <- proj4string(PEF.shp)

##make pretty figure
surf <- as.image.SpatialGridDataFrame(surf[PEF.shp,])
image.plot(surf, xaxs = "r", yaxs = "r", xlab="Easting (m)", ylab="Northing (m)", main="Trees per ha")
plot(PEF.shp, add=TRUE)


