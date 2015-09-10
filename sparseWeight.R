## Build sparse weight matrix based on (approximate) k-nearest neighbors
## Based on Spotify's ANNOY C++ library

## using RcppAnnoy and Matrix:
## 1. create index file
## 2. build trees
## 3. find "knn"-nearest neighbors, for any "knn"
## 4. build (column-oriented) sparse weight matrix

## --------------------------------------------------------------------------- #

## load packages

require(RcppAnnoy)
require(Matrix)

## --------------------------------------------------------------------------- #

knnIndex <- function(df, nt=50, save.index=NULL)
{
    ## Create (and save) index file
    ## Inputs:
    ## df: (nxd) data.frame    
    ## nt: integer number of trees for ANNOY (default: nt=50)
    ## save.index: character path to save index file (default: don't save)
    ## Output:
    ## a: index file, can be queried at any time
    
    d <- ncol(df)
    n <- nrow(df)
    
    ## init using Euclidean distance index file
    ## TODO: choose {AnnoyEuclidean, AnnoyCosine}
    a <- new(AnnoyEuclidean, d) 
    
    ## create index file
    ## NB: indexing from 0
    for (i in seq(n)) {
        v <- as( df[ i-1, ], "numeric" )
        a$addItem( i-1, v )
    }
    
    ## build trees
    a$build(nt)
    
    ## save index in path provided as "save.index"
    if ( !is.null(save.index) ) a$save(paste0(save.index,"/annoy-out.tree"))

    return(a)
    
}

sparseWeight <- function(df, nt=50, knn=6, a=NULL, no.index=TRUE)
{
    ## Create weight matrix from "knnIndex" object
    ## Inputs:
    ## df: (nxd) data.frame    
    ## nt: integer number of trees for ANNOY (default: nt=50)
    ## a: Annoy index file, or path to load one    
    ## knn: integer number of nearest neighbors for W (default: knn=6)
    ## no.index: logical return weight matrix ONLY
    ## Outputs:
    ## Wmat: (nxn) CsparseMatrix, the weight matrix
    ## b: (optional) index file
    
    d <- ncol(df)
    n <- nrow(df)
    
    ## init using Euclidean distance
    ## TODO: choose {AnnoyEuclidean, AnnoyCosine}
    if (is.character(a)) {
        b <- new(AnnoyEuclidean, d)
        b$load(a)
    } else if (is.null(a)) {
        b <- knnIndex(df, nt)
    } else {
        b <- a
    }
    
    ## TODO: build contiguity matrix
    
    ## build (row-normalized) weight matrix
    ## TODO: allow other normalizations
    Wmat <- Matrix(0, nrow=n, ncol=n)
    
    for (i in seq(n)) {
        ## get indices of knn nearest neighbors
        idx <- b$getNNsByItem(i, knn+1) # includes self
        idx <- idx[idx!=i] # drop self
        Wmat[i,idx] <- 1/length(idx) # row-normalized weights
    }
    
    ## convert to CsparseMatrix
    Wmat <- as(Wmat, "CsparseMatrix")
    
    if (no.index) { # return Wmat only
        out <- Wmat 
    } else { # return Wmat and index
        out <- list(
            "W"=Wmat,
            "index"=b
        )
    }
    
    return(out)
    
}
