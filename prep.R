

rm(list=ls())
setwd("/Users/minjay/Downloads/PEF-data")
library(MBA)
library(fields)
library(rgdal)
library(ggplot2)
library(spBayes)
library(geoR)
library(RColorBrewer)

PEF.shp <- readOGR("shapefiles","PEF-bounds")
MU.shp <- readOGR("shapefiles","MU-bounds")

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


PEF.coords <- as.matrix(PEF.data[,3:4])
Lid.coords <- as.matrix(PEF.LVIS[,1:2])

lm1 <- lm(biomass.mg.ha ~ PC1 + PC2, data = PEF.data)
vg1 <- variog(coords = coords, data = resid(lm1))
plot(vg1)
spline.obj <- spline(vg1$u, vg1$v, xmin=0, xmax=max(vg1$u))
lines(spline.obj$x, spline.obj$y)

D.max <- max(iDist(coords))


starting <- list("phi" = 3/(0.5*D.max), "sigma.sq" = 700, "tau.sq" = 500)
tuning <- list("phi" = 0.03, "sigma.sq" = 0.03, "tau.sq" = 0.03)
priors.1 <- list("beta.Norm"=list(rep(0,3), diag(1000,3)),
                 "phi.Unif"=c(3/D.max, 3/1), "sigma.sq.IG"=c(2, 700),
                 "tau.sq.IG"=c(2, 500))

cov.model <- "exponential"

n.samples = 2000
n.report = 500
verbose = TRUE

m.1 <- spLM(biomass.mg.ha ~ PC1 + PC2, data = PEF.data, coords = coords, starting=starting,
            tuning = tuning, priors = priors.1, cov.model = cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)

plot(m.1$p.theta.samples)

burn.in <- 0.5*n.samples
m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)



round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],4)

round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2) #regression coeffs

quants <- function(x){quantile(x, prob=c(0.5,0.025,0.975))}

w.summary <- apply(m.1$p.w.recover.samples, 1, quants)

##posterior median surface
par(mfrow=c(1,2))
quilt.plot(PEF.coords[,1], PEF.coords[,2], w.summary[1,])

w.mu.surf <- mba.surf(cbind(PEF.coords, w.summary[1,]), no.X=200, no.Y=200, extend=TRUE)$xyz.est
image.plot(w.mu.surf)
plot(PEF.shp, add = TRUE)
plot(MU.shp, border = "black", add = TRUE)

m.1.pred <- spPredict(m.1, pred.covars = as.matrix(cbind(rep(1, nrow(PEF.predict)), PEF.predict[,3:4])), pred.coords = Lid.coords,
                      start=0.9*n.samples)

y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, quants)


y.mu.surf <- mba.surf(cbind(Lid.coords, y.hat[1,]), no.X=200, no.Y=200, extend=TRUE)$xyz.est
y.ci.surf <- mba.surf(cbind(Lid.coords, (y.hat[3,] - y.hat[2,])), no.X=200, no.Y=200, extend=TRUE)$xyz.est

jpeg(file = "Biomass-Maps.jpg")
par(mfrow=c(1,2))
image.plot(y.mu.surf, xlab = "Easting (m)", ylab = "Northing (m)",
           main = "Predicted Biomass (Mg/ha)")
plot(PEF.shp, add = TRUE)
plot(MU.shp, border = "black", add = TRUE)
points(PEF.coords)
image.plot(y.ci.surf, xlab = "Easting (m)", ylab = "Northing (m)",
           main = "95 percent CI width")
plot(PEF.shp, add = TRUE)
plot(MU.shp, border = "black", add = TRUE)
points(PEF.coords)
dev.off()





