library(foreign)
library(reshape)
library(ggplot2)
library(fpc)
library(protoclust)
library(cluster)
library(mclust)

# ============================
# Read in Data Files
# ============================

# Contains Test Data
test <- read.arff(unz(description="Data.zip",filename="AstronomyTestData.txt"))

# Contains row numbers of the low noise observations in the test data
lownoise <- read.csv(unz(description="Data.zip",filename="LowNoiseData.txt"),header=F)

# subset
lntest=test[lownoise[,1],]


# ============================
# CLEANING DATA
# ============================

# remove ID and class information
# and select variables which are not constant
vars <- apply(na.omit(lntest[,c(-1,-87)]),2,var)

# See variables which are constant
# which(vars==0)

lndat <- na.omit(lntest[,c(-1,-87)])[,-which(vars==0)]


# ============================
# PCA
# ============================
lnpcs <- princomp(x=lndat,scores=TRUE,cor=TRUE)

lneigvec <- lnpcs$loading
lneigval <- (lnpcs$sdev)^2
#lnscores <- lnpcs$scores

# dim(lneigvec)
colnames(lneigvec) <- names(lndat)
# lneigvec[1:10,1:5]

rm(lnpcs)

# for some reason, cannot coerse "loadings" type object to a data.frame
write.csv(lneigvec, "lnPCAloadings.csv", row.names=TRUE)
lneigvec <- read.csv("lnPCAloadings.csv")
names(lneigvec)[1] <- "Variable"


# Using Variable as id variables
lneigvec.m <- melt(lneigvec, id.vars = "Variable", variable_name = "PC")
names(lneigvec.m)[3] <- "Loading"
lneigvec.m$Eigenvalue <- lneigval[lneigvec.m$PC]


# =========================
# Eigenvector Heat Plot - Width of tiles "proportional" to Eigenvalue
# =========================

ggplot(lneigvec.m, aes(PC, Variable)) +
  geom_tile(aes(fill = Loading, width = Eigenvalue)) +
  scale_fill_gradient2(low="darkblue", high="red", guide="colorbar") +
  xlab("Principal Component Eigenvector (Width Proportional to Eigenvalue)") +
  theme(axis.text.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank()) +
  ggtitle("Low Noise Data Set Heat Plot")

ggsave(filename = "LowNoiseHeatPlot.eps", width = 13.52, height = 9.39)

# =========================
# Scree Plot
# =========================
qplot(y = lneigval/sum(lneigval), x = seq(1,length(lneigval), by = 1), geom = "line") +
  geom_point() + xlab("Number of Eigenvalues") +
  ylab("Proportion of Variance Explained") +
  ggtitle("Low Noise Data Set Scree Plot")

ggsave(filename = "LowNoiseScreePlot.eps", width = 13.52, height = 9.39)

# =========================
# Cumulative Explained Plot
# =========================
qplot(y = cumsum(lneigval)/sum(lneigval), x = seq(1,length(lneigval), by = 1), geom = "line") +
  geom_point() + xlab("Number of Eigenvalues") +
  ylab("Proportion of Variance Explained") +
  ggtitle("Low Noise Data Set Cumulative Variance Explained")

ggsave(filename = "LowNoiseCumulativePlot.eps", width = 13.52, height = 9.39)


# ============================
# Clustering
# ============================

# We will cluster on 28 variables selected from PCA heat plot.
# Selecting all variables from "freq3_harmonic_rel_phase_1" to the end
# and variables "freq(1,2)_harmonic_rel_phase_(1,2,3)"

# reorder variable names
lnnewdat <- lndat[,sort(colnames(lndat))]
lnnewdat <- lnnewdat[,c(grepl(pattern = "harmonics_rel_phase",x = colnames(lnnewdat)),
                    c(64:82))]


# ============================
# kmeans
# ============================

# clustering is dumb
clusts <- kmeansruns(lnnewdat, krange = 2:25, criterion = 'asw', iter.max = 100, runs = 5, critout = TRUE)

clustsk <- data.frame(crit = clusts$crit[-1], k = 2:25)

qplot(y = crit, x = k, data = clustsk, geom = "line") +
  geom_point() + xlab("Number of Clusters") + ylab("Average Silhouette Width") +
  ggtitle("Average Silhouette Width vs Number of Clusters from k-means in Low Noise Data")

ggsave(filename = "LowNoisekmeans.eps", width = 13.52, height = 9.39)


# ============================
# hierarchical clustering
# ============================

# Hierarchical Clustering with Centroid Linkage
lnhclust = hclust(dist(lnnewdat), method="centroid")
# lnhclust25=cutree(lnhclust, k = 25) - for just 25 clusters


# MINIMAX Linkage
plnhclust <- protoclust(dist(lnnewdat))
#plnhclust25 <- protocut(plnhclust, k = 25)


# average silhouette values for centroid link and MINIMAX linkage
hsils <- matrix(rep(0, times = 2*length(2:25)),ncol=2)

for(i in 2:25){
  # centroid
  hsils[i-1,1] <- summary(silhouette(x=cutree(lnhclust, k = i),dist=dist(lnnewdat)))$avg
  # MINIMAX
  hsils[i-1,2] <- summary(silhouette(x=protocut(plnhclust, k = i)$cl,dist=dist(lnnewdat)))$avg
  print(i)
}

hsilsk <- as.data.frame(cbind(centroid = hsils[,1], minimax = hsils[,2], k = 2:25))

hsils.m <- melt(hsilsk, measure.vars = c("centroid","minimax"))

names(hsils.m) <- c("k", "linkage", "aws")

qplot(y = aws, x = k, data = hsils.m, geom = "line", group = linkage, color = linkage) +
  geom_point() + ylab("Average Silhouette Width") + xlab("Number of Clusters") +
  ggtitle("Average Silhouette Width vs Number of Clusters from Low Noise Hierarchical Clustering, by Linkage")

ggsave(filename = "LowNoiseHierarchical.jpg", width = 13.52, height = 9.39)





######Visualize Low Noise Clusters with PCA#########

# qplot(Comp.1,Comp.2,data=as.data.frame(lnscores),color=lnhclust25) +
#   ggtitle("Heirarchical Clustering with Centroid Linkage")
# 
# library(car)
# 
# clusts3 <- kmeansruns(lnnewdat, krange = 4, criterion = 'asw', iter.max = 100, runs = 5, critout = TRUE)
# 
# 
# qplot(Comp.1,Comp.2,data=as.data.frame(lnscores),color=as.factor(clusts3$cluster))
