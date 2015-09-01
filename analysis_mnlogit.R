## CMAP Land Use - code for analysis
## see prep.R for data prep

## --------------------------------------------------------------------------- #

## Outline
## 1. Subset relevant data
## 2. Generate descriptive statistics
## 3. Estimate FE/ME logits (using mnlogit, mclogit), make tables and compare

## --------------------------------------------------------------------------- #

## load packages, functions, etc.
## REMARK: "landuse.R" loads everything, this is unnecessary
source("load.R")

## --------------------------------------------------------------------------- #

## Load data
lu.10.dir = paste0(cmap.dir,"LU10/")
load(file=paste0(lu.10.dir,"lu10.Rda"),verbose=TRUE)

## --------------------------------------------------------------------------- #

## Parameters
vintage <- "2010" # Census {"2000","2010"}
n.cores <- 6 # number of cores for parallelization, where applicable

## --------------------------------------------------------------------------- #

## select relevant columns
columns <-c("Use", # recoded response
            "LANDUSE", # original response label
            ##"dboundary", # Distance to nearest county line, version 1
            "boundary.distance", # Distance to nearest county line, version 2
            "boundary.name", # FIPS code for nearest county line            
            "Longitude", # Polygon centroid longitude
            "Latitude", # Polygon centroid latitude
            ##"LUI2010_ID", # Polygon unique ID (useless?)
            "County", # County
            "BLOCKCE10", # Census block
            "TRACTCE10", # Census tract
            "Chicago", # Chicago
            "Shape_Leng", # Polygon length
            "Shape_Area", # Polygon area
            "DENSCLASS", # residential density (class) per acre
            "place.distance", # distance to census place
            "cousub.distance" # distance to county subdivision
            )

set.seed(2010)
##lu.10.sub <- lu.10[sample(1:nrow(lu.10.sub),500),columns]
lu.10.sub <- lu.10[,columns]
## current object S4 class; use only dataframe
lu.10.sub <- lu.10.sub@data
rm(lu.10)

## Set Use base category to residential
## ignore for mnlogit
## lu.10.sub$Use <- relevel(lu.10.sub$Use, ref = "Residential")

## Hacky way for mnlogit
lu.10.sub$Use <- gsub("^(Commercial)$","3.\\1",lu.10.sub$Use)
lu.10.sub$Use <- gsub("^(Industrial)$","2.\\1",lu.10.sub$Use)
lu.10.sub$Use <- gsub("^(Residential)$","1.\\1",lu.10.sub$Use)

## --------------------------------------------------------------------------- #

## 2. Generate descriptive statistics


## --------------------------------------------------------------------------- #

## 3. Generate results

## --------------------------------------------------------------------------- #
## 3a. package mnlogit: fast logit with mixed effects
## based on package mlogit

## drop missing tract/block FE's
rownames(lu.10.sub) <- NULL
lu.10.sub <- lu.10.sub[!is.na(lu.10.sub$TRACTCE10),]
lu.10.sub <- lu.10.sub[!is.na(lu.10.sub$BLOCKCE10),]

## put data in mlogit "long" format
lu.10.sub <- mlogit.data(lu.10.sub,shape="wide",choice="Use",alt.var="LANDUSE")


## model formulas

## pure distance to county boundary effect
f1 <- formula(Use~1|boundary.distance)
## + County fixed effects
f2 <- formula(Use~1|boundary.distance+County)
## + County, Chicago fixed effects
f3 <- formula(Use~1|boundary.distance+County+Chicago)
## + which boundary?
f4 <- formula(Use~1|boundary.distance+County+Chicago+boundary.name)
## + poygon area
f5 <- formula(Use~1|boundary.distance+County+Chicago+boundary.name+Shape_Area)

## county boundary + "place" boundary (townships, etc.)
f6 <- formula(Use~1|boundary.distance+place.distance)
## + 
f7 <- formula(Use~1|boundary.distance+place.distance+County)
##
f8 <- formula(Use~1|boundary.distance+place.distance+County+Chicago)
##
f9 <- formula(Use~1|boundary.distance+place.distance+County+Chicago+boundary.name)
##
f10 <- formula(Use~1|boundary.distance+place.distance+County+Chicago+boundary.name+Shape_Area)

## county boundary + "county subdivision" boundary (inc cities)
f11 <- formula(Use~1|boundary.distance+cousub.distance)
## + 
f12 <- formula(Use~1|boundary.distance+cousub.distance+County)
##
f13 <- formula(Use~1|boundary.distance+cousub.distance+County+Chicago)
##
f14 <- formula(Use~1|boundary.distance+cousub.distance+County+Chicago+boundary.name)
##
f15 <- formula(Use~1|boundary.distance+cousub.distance+County+Chicago+boundary.name+Shape_Area)


## collect formulas into a list
forms <- lapply(ls(pattern="^f\\d+$"), function(x) get(x))

## name the list
labels <- as.character(1:(length(forms)))
names(forms) <- labels
## initialize list to store results
## out <- list()
out <- vector(mode="list", length=length(forms))

## iterate estimation of each formula in list
for (i in labels) {
    f <- forms[[i]]
    print(f)
    ## estimate
    fit <- mnlogit(f1,data=lu.10.sub,ncores=n.cores)
    ## save results
    out[[i]] <- fit
}

## clean up
rm(f,fit)

## ALTERNATIVE: Try sqldf

## ALTERNATIVE: update
## f0 <- formula(Use ~ 1 | boundary.distance)
## fit <- mnlogit(f0,data=lu.10.sub)

## f0.0 <- update(f0,~.+County-1)


## Save!
save(forms,labels,out,file=paste0(lu.10.dir,"out_mnlogit.Rda"))

## --------------------------------------------------------------------------- #
