#### CBF missing data sim files;

library(dtw)
#library(SimilarityMeasures)


getwd()
setwd('..')

load('CBF.rda')
class(CBF)
names(CBF)
str(CBF)

#cbf <- matrix(data = CBF$data_test, nrow = 900, ncol = 128, byrow = FALSE,  dimnames = list(CBF$labels_test))
cbf <- matrix(data = CBF$data_train, nrow = 30, ncol = 128, byrow = FALSE,  dimnames = list(CBF$labels_train))
dim(cbf)

summary(factor(CBF$labels_train))
cbf[1:10, 1:5]

###visualizing the 30 profiles:
par(mfrow=c(5,6))
for(i in 1:30){
  i = i
  plot(1:128, cbf[i,], ylab='', xlab=i, main= CBF$labels_train[i], type = 'l', ylim=c(-3,3))
}

par(mfrow=c(1,1))

###setting up sparse version of CBF:
sparsify <- seq(from=0,to=128, by=4 ); 
sparsify <- sparsify[sparsify != 0]
cbf_sparse <- cbf[, sparsify]
dim(cbf_sparse)
head(cbf_sparse)

par(mfrow=c(5,6))
for(i in 1:30){
  i = i
  plot(1:length(sparsify), cbf_sparse[i,], ylab='', xlab=i, main= dimnames(cbf_sparse)[[1]][i], type = 'l', ylim=c(-3,3))
}
par(mfrow=c(1,1))

###true labels for the groups:
true_labels <- CBF$labels_train
summary(factor(true_labels))

###examples
eucl_clust_sparse <- hclust(dist(cbf_sparse))
clusttemp <- cutree(eucl_clust_sparse, k=4)
adjustedRand(clusttemp, true_labels, randMethod = c('HA', 'FM', 'Jaccard'))  ###from clues package
ClusterPurity(clusttemp, true_labels)
ClusterEntropy(clusttemp, true_labels)

# tempfunc <- function(clusters, classes) {
#   sum(apply(table(classes, clusters), 2, max)) / length(clusters)
# }
# tempfunc(clusttemp, true_labels)


#########################
#### Frechet distance function: from kml package

##incorporating the Frech distance script used in kml pkg:
##do NOT load KML pkg; distFrechet function there produces unnecess output;
##load the fxn directly here:
distFrechet <- function(Px,Py,Qx,Qy,timeScale=0.0001,FrechetSumOrMax="max"){
  missingsP <- is.na(Px)|is.na(Py);Px <- Px[!missingsP];Py <- Py[!missingsP]
  missingsQ <- is.na(Qx)|is.na(Qy);Qx <- Qx[!missingsQ];Qy <- Qy[!missingsQ]
  Px <- Px*timeScale
  Qx <- Qx*timeScale
  
  maxP <- length(Px)
  maxQ <- length(Qx)
  Mdist <- Mfret <- matrix(0,maxP,maxQ,dimnames=c(list(paste("P",1:maxP,sep=""),paste("Q",1:maxQ,sep=""))))
  for(i in 1:maxP){
    for (j in 1:maxQ){
      Mdist[i,j] <- dist(rbind(c(Px[i],Py[i]),c(Qx[j],Qy[j])))
      if(i == 1 && j == 1){Mfret[1,1] = Mdist[1,1]}
      if(i > 1 && j == 1){Mfret[i,1] = do.call(FrechetSumOrMax , list( Mfret[i-1,1] , Mdist[i,1] ) )}
      if(i == 1 && j > 1){Mfret[1,j] = do.call(FrechetSumOrMax , list( Mfret[1,j-1] , Mdist[1,j] ) )}
      if(i > 1  && j > 1){Mfret[i,j] = do.call(FrechetSumOrMax , list( min(Mfret[i-1,j],Mfret[i-1,j-1],Mfret[i,j-1]) , Mdist[i,j] ) )}
    }
  }
  # print(Mdist)
  # print(Mfret)
  # cat("\n\n Result =",Mfret[maxP,maxQ],"\n");
  return(Mfret[maxP,maxQ])
}


## distFrechet(Px,Py,Qx, Qy, timeScale=0.1, FrechetSumOrMax = "max")

##testing on cbf:
##for sparse cbf, the time stamps are in vector 'sparsify'
distFrechet(sparsify, cbf_sparse[1,], sparsify, cbf_sparse[2,])

####need to generate the full cbf data matrix with random missingness
###there are 30 x 128 =  3840 positions in the data matrix
### 10% missingness removes 384 values
### 25% missingness removes  960 values
### 50% missingness removes 1920 values

####need to generate the sparse cbf data matrix with random missingness
###there are 30 x 32 =  960 positions in the data matrix
### 10% missingness removes 96 values
### 25% missingness removes  240 values
### 50% missingness removes 480 values

dim(cbf_sparse)
miss_index50pct <- sample(1:960, 480, replace = FALSE);
miss_index50pct <- sort(miss_index50pct); head(miss_index50pct)

cbf_sparsevec <- as.vector(t(cbf_sparse))
length(cbf_sparsevec)
cbf_sparsevec[1:20]; cbf_sparse[1:2, 1:20]
cbf_sparsevecNA <- cbf_sparsevec


cbf_sparsevecNA[miss_index50pct] <- NA
head(miss_index50pct); tail(miss_index50pct)
head(cbf_sparsevecNA); tail(cbf_sparsevecNA)
cbf_sparseNA <- matrix(cbf_sparsevecNA, nrow=dim(cbf_sparse)[1], byrow=TRUE)    ##this is the desired data matrix with the NAs
dim(cbf_sparseNA); cbf_sparseNA[1:5, 1:10]

val1 <- na.omit(cbf_sparseNA[1,])
pos1 <- which(is.na(cbf_sparseNA[1, ]) == FALSE)
pos1
cbf_sparseNA[2,]
distFrechet(pos1, val1, pos2, val2)
DTW(matrix(c(pos1, val1), ncol=2), matrix(c(pos2, val2), ncol=2))
rm(pos1, val1, pos2, val2)

###setting up sim for sparsified data, 50 pct missingness:
nsims <- 1000

###storage for the evals:
##evals: 'HA', 'FM', 'Jaccard', purity, entropy 
sparse_50pctMiss_Frech_evals <- matrix(0, nrow=nsims, ncol = 5)
sparse_50pctMiss_DTW_evals <- matrix(0, nrow=nsims, ncol = 5)

starttime = Sys.time()
for(k in 1:nsims){
 ##generate missing data from cbfsparse
  miss_index50pct <- sort(sample(1:960, 480, replace = FALSE))
  cbf_sparsevecNA <- as.vector(t(cbf_sparse))
  cbf_sparsevecNA[miss_index50pct] <- NA
  cbf_sparseNA <- matrix(cbf_sparsevecNA, nrow=dim(cbf_sparse)[1], byrow=TRUE)
  
  ##get dist matr for ea distance for this simulated dataset
  distmatFrech <- matrix(0, nrow=dim(cbf_sparse)[1], ncol=dim(cbf_sparse)[1])
  distmatDTW <- matrix(0, nrow=dim(cbf_sparse)[1], ncol=dim(cbf_sparse)[1])
  
  for(i in 1:(dim(distmatFrech)[1]-1)){
    vali <- na.omit(cbf_sparseNA[i,])
    posi <- which(is.na(cbf_sparseNA[i, ]) == FALSE)
    for(j in (i+1):dim(distmatFrech)[1]){
      valj <- na.omit(cbf_sparseNA[j,])
      posj <- which(is.na(cbf_sparseNA[j, ]) == FALSE)
      distmatFrech[i,j] <- distFrechet(posi, vali, posj, valj)
      distmatDTW[i,j] <- DTW(matrix(c(posi, vali), ncol=2), matrix(c(posj, valj), ncol=2))
  }
  }
clustFrech <- hclust(as.dist(t(distmatFrech)))
cut3Frech <- cutree(clustFrech, k=3)
#cut4Frech <- cutree(clustFrech, k=4)
clustDTW <- hclust(as.dist(t(distmatDTW)))
cut3DTW <- cutree(clustDTW, k=3)
#cut4DTW <- cutree(clustDTW, k=4)

sparse_50pctMiss_Frech_evals[k, ] <- c(adjustedRand(cut3Frech, true_labels, randMethod = c('HA', 'FM', 'Jaccard')),
 ClusterPurity(cut3Frech, true_labels), ClusterEntropy(cut3Frech, true_labels))
sparse_50pctMiss_DTW_evals[k, ] <- c(adjustedRand(cut3DTW, true_labels, randMethod = c('HA', 'FM', 'Jaccard')),
                                       ClusterPurity(cut3DTW, true_labels), ClusterEntropy(cut3DTW, true_labels))
rm(clustFrech, clustDTW)
}

Sys.time() - starttime

length(which(sparse_50pctMiss_Frech_evals[,1] > sparse_50pctMiss_DTW_evals[,1]))

save.image(file = 'SimOutputWorkspaces/CBF_Sim_Sparse50PctMiss.RData')
