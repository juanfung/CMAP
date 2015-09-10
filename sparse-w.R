## CMAP Land Use - code for analysis
## see prep.R for data prep

## --------------------------------------------------------------------------- #

## Outline
## 1. Subset relevant data
## 2. Generate descriptive statistics
## 3. Estimate spatial probits (using McSpatial or spatialprobit)

## --------------------------------------------------------------------------- #

## load packages, functions, etc.
## REMARK: "landuse.R" loads everything, this is unnecessary
source("load.R")

library(RcppAnnoy)
library(Matrix)

## --------------------------------------------------------------------------- #

## Load data
lu.10.dir = paste0(cmap.dir,"LU10/")
load(file=paste0(lu.10.dir,"lu10.Rda"),verbose=TRUE)

## --------------------------------------------------------------------------- #

## Parameters
vintage <- "2010" # Census {"2000","2010"}
n.cores <- 6 # number of cores for parallelization, where applicable

n <- nrow(lu.10) # number of observations

## N.B.: ANNOY indexing begins at 0!
idx <- c(0:(n-1))
lu.10$Index <- idx

save(lu.10,file=paste0(lu.10.dir,"lu10.Rda"))

## --------------------------------------------------------------------------- #

## select relevant columns
columns <-c("Index",            
            "Longitude", # Polygon centroid longitude
            "Latitude" # Polygon centroid latitude           
            )

##set.seed(2010)
##lu.10.sub <- lu.10[sample(1:nrow(lu.10.sub),500),columns]

lu.10.sub <- lu.10[,columns]
## current object S4 class; use only dataframe
lu.10.sub <- lu.10.sub@data
rm(lu.10)


## drop missing tract/block FE's
rownames(lu.10.sub) <- NULL
##lu.10.sub <- lu.10.sub[!is.na(lu.10.sub$TRACTCE10),]
##lu.10.sub <- lu.10.sub[!is.na(lu.10.sub$BLOCKCE10),]

## put data in mlogit "long" format
##lu.10.sub <- mlogit.data(lu.10.sub,shape="wide",choice="Use",alt.var="LANDUSE")

## --------------------------------------------------------------------------- #

## 2. Build sparse weight matrix based on (approximate) k-nearest neighbors
## Based on Spotify's ANNOY C++ library

## --------------------------------------------------------------------------- #


## using RcppAnnoy and Matrix:
## 1. create index file
## 2. build trees
## 3. find "knn"-nearest neighbors, for any "knn"
## 4. build (column-oriented) sparse weight matrix

## example:

d <- 2 # data dimension (lon/lat)
nt <- 50 # number of trees
knn <- 6 # number of nearest neighbors

a <- new(AnnoyEuclidean, d) # init using Euclidean distance (alt: cosine)

## create index file
for (i in seq(n)) {
    v <- as( lu.10.sub[lu.10.sub[["Index"]]==i-1, c("Longitude","Latitude")], "numeric" )
    a$addItem(i-1, v)
}

## build trees
a$build(nt)

## save
a$save("./Data/Misc/annoy-out.tree")

## build contiguity matrix
##Cmat <- Matrix(0,nrow=n,ncol=n) # unnecessary?

## build weight matrix
Wmat <- Matrix(0,nrow=n,ncol=n)
for (i in seq(n)) {
    ## get indices of knn nearest neighbors
    idx <- a$getNNsByItem(i, knn+1) # includes self
    idx <- idx[idx!=i] # drop self
    Wmat[i,idx] <- 1/length(idx) # row-normalized weights
}

## convert to CsparseMatrix?
Wmat <- as(Wmat, "CsparseMatrix")

## --------------------------------------------------------------------------- #

## wrap the steps above in a function or two

## load functions
## knnIndex: build index
## sparseWeight: compute knn-weight matrix
source("sparseWeight.R")

df <- lu.10[,c("Longitude","Latitude")]

## index takes time to build...
a <- knnIndex(df, nt=nt) 

## ...but once built, Wmat is fast
Wmat1 <- sparseWeight(df, knn=knn, a=a)

Wmat2 <- sparseWeight(df, knn=2*knn, a=a)

all.equal(Wmat,Wmat1)

## --------------------------------------------------------------------------- #

## save Wmat
save(Wmat, file=paste0(cmap.dir,"Misc/Wmat.Rda"))
