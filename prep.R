rm(list=ls())
setwd("/Users/minjay/Downloads/PEF-data")
library(MBA)
library(fields)
library(rgdal)
library(ggplot2)
library(spBayes)

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


lm1 <- lm(biomass.mg.ha ~ PC1 + PC2, data = PEF.data)
vg1 <- variog(coords = PEF.data[,1:2], data = resid(lm1))
plot(vg1)









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


