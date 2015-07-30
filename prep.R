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
<<<<<<< Updated upstream
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

plot(m.1$p.theta.samples)

burn.in <- 0.5*n.samples
m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)


m.1.pred <- spPredict(m.1, pred.covars = as.matrix(cbind(rep(1, nrow(PEF.predict)), PEF.predict[,3:4])), pred.coords = Lid.coords,
                      start=0.5*n.samples)
y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, mean)

quilt.plot(coords, y.hat)


round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)

round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2) #regression coeffs

quants <- function(x){quantile(x, prob=c(0.5,0.025,0.975))}

w.summary <- apply(m.1$p.w.recover.samples, 1, quants)

##posterior median surface
par(mfrow=c(1,2))
quilt.plot(x[,1], x[,2], w.summary[1,])

w.mu.surf <- mba.surf(cbind(x, w.summary[1,]), no.X=200, no.Y=200, extend=TRUE)$xyz.est
image.plot(w.mu.surf)





=======
  
coords <- as.matrix(PEF.data[,3:4])
>>>>>>> Stashed changes

lm1 <- lm(biomass.mg.ha ~ PC1 + PC2, data = PEF.data)
vg1 <- variog(coords = coords, data = resid(lm1))
plot(vg1, title = )

D.max <- max(iDist(coords))
starting <- list("phi" = 3/0.5*D.max, "sigma.sq" = 700, "tau.sq" = 500)
tuning <- list("phi" = 0.04, "sigma.sq" = 0.04, "tau.sq" = 0.04)

priors.1 <- list("beta.Norm"=list(rep(0,3), diag(1000,3)),
                 "phi.Unif"=c(3/D.max, 3/1), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))

cov.model <- "exponential"

n.samples = 2000
n.report = 500
verbose = TRUE

<<<<<<< Updated upstream
##Plot of all curves
p <- ggplot(data = PEF.plots.long.sort, aes(x = height, y = energy, group = id, colour = MU))
p + geom_line()

coords <- PEF.plots[,c("x.coords", "y.coords")]
stems.ha <- PEF.plots$stems.ha
=======
m.1 <- spLM(basal.area.m2.ha ~ PC1 + PC2, data = PEF.data, coords=coords, starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)
plot(m.1)
>>>>>>> Stashed changes

burn.in <- 0.5*n.samples
m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)
coords <- cbind(PEF.LVIS[, 1], PEF.LVIS[, 2])
spPredict(m.1, pred.covars=PEF.predict, pred.coords=Lid.coords, start=burn.in, thin=50, n.report=1)

y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, mean)

<<<<<<< Updated upstream

=======
quilt.plot(coords, y.hat)
>>>>>>> Stashed changes
