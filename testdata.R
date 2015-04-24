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

# ============================
# CLEANING DATA
# ============================

# remove ID and class information
# and select variables which are not constant
vars <- apply(na.omit(test[,c(-1,-87)]),2,var)

# See variables which are constant
# which(vars==0)

dat <- na.omit(test[,c(-1,-87)])[,-which(vars==0)] #6631 removed

# ============================
# PCA
# ============================
pcs <- princomp(x=dat,scores=TRUE,cor=TRUE)

eigvec <- pcs$loading
eigval <- (pcs$sdev)^2
#scores <- pcs$scores

# dim(eigvec)
colnames(eigvec) <- names(dat)
# eigvec[1:10,1:5]

rm(pcs)

# for some reason, cannot coerse "loadings" type object to a data.frame
write.csv(eigvec, "PCAloadings.csv", row.names=TRUE)
eigvec <- read.csv("PCAloadings.csv")
names(eigvec)[1] <- "Variable"


# Using Variable as id variables
eigvec.m <- melt(eigvec, id.vars = "Variable", variable_name = "PC")
names(eigvec.m)[3] <- "Loading"
eigvec.m$Eigenvalue <- eigval[eigvec.m$PC]

# =========================
# Eigenvector Heat Plot - Width of tiles "proportional" to Eigenvalue
# =========================

ggplot(eigvec.m, aes(PC, Variable)) +
  geom_tile(aes(fill = Loading, width = Eigenvalue)) +
  scale_fill_gradient2(low="darkblue", high="red", guide="colorbar") +
  xlab("Principal Component Eigenvector (Width Proportional to Eigenvalue)") +
  theme(axis.text.x = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank())+
  ggtitle("Test Data Set Heat Plot")

ggsave(filename = "TestDataHeatPlot.jpg", width = 13.52, height = 9.39)

# =========================
# Scree Plot
# =========================
qplot(y = eigval/sum(eigval), x = seq(1,length(eigval), by = 1), geom = "line") +
  geom_point() + xlab("Number of Eigenvalues") +
  ylab("Proportion of Variance Explained") +
  ggtitle("Test Data Set Scree Plot")

ggsave(filename = "TestDataScreePlot.jpg", width = 13.52, height = 9.39)

# =========================
# Cumulative Explained Plot
# =========================
qplot(y = cumsum(eigval)/sum(eigval), x = seq(1,length(eigval), by = 1), geom = "line") +
  geom_point() + xlab("Number of Eigenvalues") +
  ylab("Proportion of Variance Explained") +
  ggtitle("Test Data Set Cumulative Variance Explained")

ggsave(filename = "TestDataCumulativePlot.jpg", width = 13.52, height = 9.39)

# ============================
# Clustering
# ============================

#########Partitioning dataset into five subsets###########
### Creating Five Sets#####
N=1:dim(dat)[1]
n1=sample(N,8699)
n2=sample(N[-n1],8699)
n3=sample(N[-c(n1,n2)],8699)
n4=sample(N[-c(n1,n2,n3)],8698)
n5=N[-c(n1,n2,n3,n4)]

# We will cluster on 28 variables selected from PCA heat plot.
# Selecting all variables from "freq3_harmonic_rel_phase_1" to the end
# and variables "freq(1,2)_harmonic_rel_phase_(1,2,3)"

# reorder variable names
newdat <- dat[,sort(colnames(dat))]
newdat <- newdat[,c(grepl(pattern = "harmonics_rel_phase",x = colnames(newdat)),
                    c(64:82))]

# creating data chunks
new1=newdat[n1,]
new2=newdat[n2,]
new3=newdat[n3,]
new4=newdat[n4,]
new5=newdat[n5,]

# ===============================
# Running K Means on subsets
# ===============================

kc1 <- kmeansruns(new1, krange = 2:25, criterion = 'asw', iter.max = 100, runs = 5, critout = TRUE)
kc2 <- kmeansruns(new2, krange = 2:25, criterion = 'asw', iter.max = 100, runs = 5, critout = TRUE)
kc3 <- kmeansruns(new3, krange = 2:25, criterion = 'asw', iter.max = 100, runs = 5, critout = TRUE)
kc4 <- kmeansruns(new4, krange = 2:25, criterion = 'asw', iter.max = 100, runs = 5, critout = TRUE)
kc5 <- kmeansruns(new5, krange = 2:25, criterion = 'asw', iter.max = 100, runs = 5, critout = TRUE)

silhouettes <- cbind(Subset1 = kc1$crit[-1], Subset2 = kc2$crit[-1],
                     Subset3 = kc3$crit[-1], Subset4 = kc4$crit[-1],
                     Subset5 = kc5$crit[-1])

silhouettes <- melt(silhouettes)
names(silhouettes) <- c("k", "Subset", "AverageSilhouette")

qplot(data = silhouettes, y = AverageSilhouette, x = k+1,
      color = Subset, lty = Subset, group = Subset, geom="line") +
  geom_point() + xlab("Number of Clusters") + ylab("Average Silhouette Value") +
  ggtitle("Average Silhouette vs Number of Clusters by Subset of Test Data")
  
ggsave(filename = "TestDatakmeans.jpg", width = 13.52, height = 9.39)

# ============================
# hierarchical clustering
# ============================

# Hierarchical Clustering with Centroid Linkage

hclust1 <- hclust(dist(new1), method="centroid")
hclust2 <- hclust(dist(new2), method="centroid")
hclust3 <- hclust(dist(new3), method="centroid")
hclust4 <- hclust(dist(new4), method="centroid")
hclust5 <- hclust(dist(new5), method="centroid")


# MINIMAX Linkage

phclust1 <- protoclust(dist(new1))
phclust2 <- protoclust(dist(new2))
phclust3 <- protoclust(dist(new3))
phclust4 <- protoclust(dist(new4))
phclust5 <- protoclust(dist(new5))

# average silhouette values for centroid link and MINIMAX linkage on all subsets

# This list contains all 10 cluster objects created above (5 from centroid, 5 from minimax)
thelist <- list(hclust1, phclust1, hclust2, phclust2, hclust3, phclust3, hclust4, phclust4, hclust5, phclust5)

# This list contains all of the data chunks
thedata <- list(new1, new2, new3, new4, new5)

# Allocating space for entire data set.
hsils <- matrix(rep(0, times = 5*3*length(2:25)),ncol=3)

for (j in 1:5){
  
  # Allocating space for the jth subset's data
  hsilss <- matrix(rep(0, times = 3*length(2:25)),ncol=3)
  
  for(i in 2:25){ # where i indexes the number of clusters k
    # centroid
    hsilss[i-1,1] <- summary(silhouette(x=cutree(thelist[[(2*j)-1]], k = i),dist=dist(thedata[[j]])))$avg
    # MINIMAX
    hsilss[i-1,2] <- summary(silhouette(x=protocut(thelist[[2*j]], k = i)$cl,dist=dist(thedata[[j]])))$avg
    # save which subset
    hsilss[i-1,3] <- j
    print(i)
  }
  
  # Putting the jth subset's results into the entire matrix
  hsils[24*(j-1) + (seq(2, 25 ,by = 1) - 1), ] <- hsilss
  print(j)
}

# create data frame from the data matrix, naming columns
hsilsk <- as.data.frame(cbind(centroid = hsils[,1], minimax = hsils[,2], subset = hsils[,3], k = rep(2:25,times=5)))

hsils.m <- melt(hsilsk, id.vars = c("subset","k"), measure.vars = c("centroid","minimax"))

names(hsils.m)[3:4] <- c("linkage","aws")

hsils.m$subset <- factor(hsils.m$subset)

qplot(y = aws, x = k, data = hsils.m, geom = "line", group = subset, color = subset) +
  facet_wrap(~linkage) + geom_point() + ylab("Average Silhouette Width") + xlab("Number of Clusters") +
  ggtitle("Average Silhouette Width vs Number of Clusters from Hierarchical Clustering, by Linkage")

ggsave(filename = "TestDataHierarchical.eps", width = 13.52, height = 9.39)

