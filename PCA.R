library(lattice)

SDMdat = read.csv("SDMdat.csv", header = TRUE)

VAR.ir <- SDMdat[, 2:15]
SDM.species <- SDMdat[, 1]

#PCA via fit
fit <- princomp(VAR.ir, cor=TRUE)
fit$sdev ^2
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit.sc <- fit$scores # the principal components
write.csv(fit.sc, file = "fit_pca.csv")
biplot(fit)

library(psych)
fit.ro <- principal(VAR.ir, nfactors=14, rotate="varimax", scores = T)
fit.ro

plot(fit.ro$values, type="b", ylab="Eigenvalues",
     xlab="Component", lab=c(5,5,5))



#Run PCA via prcomp
VAR.pca <- prcomp(VAR.ir,
                 center = TRUE,
                 scale = TRUE,
                 retx=TRUE) 

plot(VAR.pca, type = "l")
print(VAR.pca)
summary(VAR.pca)
varimax7 <- varimax(VAR.pca$rotation)
print(varimax7)

#Visualise variable loadings
load    <- VAR.pca$rotation
sorted.loadings <- load[order(load[, 1]), 1]
myTitle <- "Loadings Plot for PC1" 
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")

sorted.loadings <- load[order(load[, 2]), 2]
myTitle <- "Loadings Plot for PC2"
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")

sorted.loadings <- load[order(load[, 3]), 3]
myTitle <- "Loadings Plot for PC3"
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")

sorted.loadings <- load[order(load[, 4]), 4]
myTitle <- "Loadings Plot for PC3"
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")

# Now draw the BiPlot
biplot(VAR.pca, choices = c(1,2), cex=c(1, 0.7))

#PCA plot with ellipse groupings for each species
library(devtools)
install_github("vqv/ggbiplot")

library(ggbiplot)
g <- ggbiplot(VAR.pca, choices = c(1,3), obs.scale = 1, var.scale = 1, alpha = 0,
              groups = SDM.species, ellipse = TRUE, labels = NULL,
              labels.size = 0, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

#Possible help with PCA interpretation
library (ggplot2)

theta <- seq(0,2*pi,length.out = 100)
circle <- data.frame(x = cos(theta), y = sin(theta))
p <- ggplot(circle,aes(x,y)) + geom_path()

loadings <- data.frame(VAR.pca$rotation, 
                       .names = row.names(VAR.pca$rotation))
p + geom_text(data=loadings, 
              mapping=aes(x = PC1, y = PC2, label = .names, colour = .names)) +
  coord_fixed(ratio=1) +
  labs(x = "PC1", y = "PC2")


#FactoMineR approach
library(FactoMineR)

pca3 = PCA(VAR.ir, graph = FALSE)
pca3$eig #eigenvalues
pca3$var$coord # correlations between variables and PCs
head(pca3$ind$coord) # PCs (scores)

#### PCA histograms ####

MWdat = read.csv("MWdat.csv", header = TRUE)
tim.pca2 <- as.numeric(MWdat$PC2)

d <- density(tim.pca2, na.rm=T)
plot(d, col="black", lwd=2)

h <- hist(tim.pca2, freq=FALSE, col="grey", xlim=c(-4,6), breaks=20)
lines(d, col="black", lwd=3)

#tiff(filename = "[filename].tiff", res = 600, pointsize = 2)
#plot(1:100)
#dev.off()