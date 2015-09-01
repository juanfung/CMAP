## Process CMAP land use inventory data, 2010
## source: https://datahub.cmap.illinois.gov/group/land-use-inventories

## Outline
## 1. Create County borders shapefile
## 2. Prep land use data
## 3. Compute distances from parcels to nearest county border
## 4. Assign Census tracts and blocks; Chicago dummy
## 5. Aggregate Land Use indicators

## --------------------------------------------------------------------------- #

## load packages, functions, etc.
## REMARK: "landuse.R" loads everything, this is unnecessary
source("load.R")

## --------------------------------------------------------------------------- #

## FIPS County Codes
## State: IL=17
## Counties: Cook=031; DuPage=043; Kane=089; Kendall=093; Lake=097; McHenry=111; Will=197

counties <- c("Cook","DuPage","Kane","Kendall","Lake","McHenry","Will")
fips <- c("031","043","089","093","097","111","197")
ilfips <- paste0("17",fips)
names(fips) <- counties
names(ilfips) <- counties

## --------------------------------------------------------------------------- #

## Path to Census shapefiles
census.dir <- paste0(cmap.dir,"/Census/")
##tiger.path <- paste0(census.dir,"TIGER Shapefiles/")

## Benchmark:
## bench <- 2010

## Vintage:
vintage <- "2010" # "2000" or "2010"

tiger <- paste0(census.dir,vintage,"/") # eg, "/2000 Tiger Shapefiles"

## --------------------------------------------------------------------------- #

## 1. Create County Lines shapefile for computing distances from polygon centroids
## Use Census Tiger Line shapefiles at block group level for high-res County boundaries
## (benchmark Census 2010, vintage Census 2000/2010)

## select group of shapefiles {tabblock, tract, county,...}
coshp <- paste0("tabblock",substr(vintage,3,4))

## make County Lines object (and County Polygons)
colines <- countyLines(fips.list=ilfips,tl.res=coshp,tl.dir=tiger)
row.names(colines[["County.Lines"]]@lines) <- ilfips

## save
save(colines,file=paste0(cmap.dir,"Misc/colines",vintage,".Rda"))

## Make shapefile of County Lines

## create data.frame
df <- as.data.frame(cbind("COUNTYFP"=ilfips,
                          "Vintage"=as.factor(rep(vintage,length(fips)))))

## make sure row names match for next step
row.names(colines[["County.Lines"]]) <- names(ilfips)
## rownames(df) <- row.names(colines[["County.Lines"]])
## convert to SpatialLinesDataFrame
the.lines <- SpatialLinesDataFrame(colines[["County.Lines"]],data=df)
## write to shapefile
writeOGR(the.lines,dsn=cmap.dir,layer="the_lines",driver="ESRI Shapefile")

## repeat for County Polygons
rownames(df) <- row.names(colines[["County.Polygons"]])
the.polys <- SpatialLinesDataFrame(colines[["County.Polys"]],data=df)
writeOGR(the.polys,dsn=cmap.dir,layer="the_polys",driver="ESRI Shapefile")


## N.B.: Census not in same projection as CMAP

## --------------------------------------------------------------------------- #

## 1. Prep Land Use 2010 data
## N.B.:
## Shapefile created in ArcMap; original data in .gdb format

lu.10.dir = paste0(cmap.dir,"LU10/shapefile")

shp.files = list.files(lu.10.dir)
layer = shp.files[grep("\\.shp$",shp.files)]
layer = substr(layer,1,nchar(layer)-4)
lu.10 = readOGR(lu.10.dir,layer)

## Subset by {residential, commercial, industrial}

## LU codes:
## 1000 - Urbanized
## - 1100 - residential (1110 single-fam)
## - 1200 - commercial
## - 1300 - institutional
## - 1400 - industrial
## - 1500 - transportation, utility, waste
## 2000 - agricultural
## 3000 - open space
## 4000 - vacant
## 5000 - water
## 6000 - non-parcel areas
## 9999 - not classifiable

lu.10 = lu.10[grep("^1(1|2|4)",lu.10$LANDUSE),]

## label land use categories
lu.10[grep("^11",lu.10$LANDUSE),"Use"] = "Residential"
lu.10[grep("^12",lu.10$LANDUSE),"Use"] = "Commercial"
lu.10[grep("^14",lu.10$LANDUSE),"Use"] = "Industrial"
lu.10$Use = as.factor(lu.10$Use)

## drop unused levels
for (i in names(lu.10)) {
    if (is.factor(lu.10[[i]]))
        ## levels(droplevels(lu.10.sub[[i]]))
        lu.10[[i]] <- factor(lu.10[[i]])
}

## Assign County

## Note:
## LUI2010ID = 13-character code CCCTTSSLL1234
## CCC = Three-digit county FIPS code
## TT = Township number (county-designated, usually first two digits of the PIN number)
## SS = Section number (county-designated, usually third and fourth digits of the PIN)
## LL = General land use category (first two digits of land use code)
## 1234 = Sequential numbers, beginning with 0001, for each unique combination of county/township/section/land use class.

for (i in 1:length(fips)) {
    lu.10[grep(paste0("^",fips[i]),lu.10$LUI2010_ID),"County"] = names(ilfips[i])
}

lu.10$County = as.factor(lu.10$County)

## Compute centroid lon/lat

## project to Census CRS
lu.10 <- spTransform(lu.10, CRSobj = CRS(proj4string(colines[["County.Lines"]])))

## coordinates(polygon) returns centroid lon/lat; compare with ArcMap
centroids <- coordinates(lu.10) 
lu.10$Longitude <- centroids[,1]
lu.10$Latitude <- centroids[,2]
rm(centroids)

## --------------------------------------------------------------------------- #

## 3. Compute distance to nearest county line


## ALTERNATIVE:
## geosphere::dist2Line



lu.10$dboundary <- geoshape(lu.10$Longitude,lu.10$Latitude,linefile=colines[["County.Lines"]])

## ALTERNATIVE: Compute by county

a <- NULL
n <- 1
for (f in ilfips) {
    message("Computing distance to ",f)
    d <- geoshape(lu.10$Longitude,
                  lu.10$Latitude,
                  linefile=colines[["County.Lines"]][f])
    if (n==1) {
        a <- data.frame(d)
    } else {
        a <- cbind(a,d)
    }
    colnames(a)[n] <- paste0("d.",f)
    n <- n+1
}


minRows <- apply(a, MARGIN=1, FUN=which.min)
dRows <- apply(a, MARGIN=1, FUN=min)

a$boundary.distance <- dRows

a$boundary.name <- as.factor(ilfips[minRows])

row.names(a) <- row.names(lu.10)
lu.10 <- spCbind(lu.10,a)


## --------------------------------------------------------------------------- #

## 4. Assign Census Tract and Block

tl.shp <- list.files(tiger,pattern=paste(coshp,"shp$",sep="."))
vin <- substr(vintage,3,4)

lu.10.plus <- NULL
for (f in 1:length(ilfips)) {
    fp <- ilfips[f]
    ## Get shapefile name
    layer <- tl.shp[grep(fp,tl.shp)]
    layer <- substr(layer,1,nchar(layer)-4)

    ## Import shapefile
    tl <- readOGR(tiger,layer)

    ## Re-project shapefile to LU.10 projection
    tl <- spTransform(tl, CRSobj = CRS(proj4string(lu.10)))

    sub <- lu.10[lu.10$County==names(fp),]
    lmat <- SpatialPoints(sub[,c("Longitude","Latitude")],proj4string=CRS(proj4string(lu.10)))
    over.pts <- over(lmat,tl)

    sub <- spCbind(sub,over.pts[,paste0(c("TRACTCE","BLOCKCE"),vin)])

    ## first pass flag
    if (f==1) {
        lu.10.plus <- sub
    } else {
        lu.10.plus <- spRbind(lu.10.plus,sub)
        }
    ## Clean up
    rm(layer,tl,sub,lmat,over.pts)
}

lu.10 <- lu.10.plus
rm(lu.10.plus)


## Assign Chicago dummy

chi.dir <- paste0(census.dir,"Chicago/",vintage)
chi <- list.files(chi.dir,pattern="\\.shp$")
chi <- readOGR(chi.dir,substr(chi,1,nchar(chi)-4))

chi.tract <- unique(as.character(eval(parse(text=paste0("chi$TRACTCE",vin)))))

lu.10.chi <- lu.10[eval(parse(text=paste0("lu.10$TRACTCE",vin))) %in% chi.tract,]
lu.10.else <- lu.10[!(eval(parse(text=paste0("lu.10$TRACTCE",vin))) %in% chi.tract),]

lu.10.chi$Chicago <- as.factor("Yes")
lu.10.else$Chicago <- as.factor("No")

lu.10 <- spRbind(lu.10.chi,lu.10.else)

rm(chi,chi.tract)


## --------------------------------------------------------------------------- #

## 5. SAVE

lu.10.dir <- gsub("(.*)shapefile","\\1",lu.10.dir)

save(lu.10,file=paste0(lu.10.dir,"lu10.Rda"))

## --------------------------------------------------------------------------- #
