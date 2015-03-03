# Process CMAP land use inventory data
# source: https://datahub.cmap.illinois.gov/group/land-use-inventories

# Outline
# 1. Add County boundaries to CMAP data
# 2. Compute distances from parcels to nearest county (or state?) border
# 3. Assign Census tracts and blocks
# 4. Aggregate Land Use indicators

# --------------------------------------------------------------------------- #

# Load Packages
library(rgdal)
library(rgeos)
library(maptools)
library(McSpatial)

# --------------------------------------------------------------------------- #

# Path to data; change as necessary
path.to.data <- "/home/juanfung/Dropbox/Data"
path.to.code <- "/home/juanfung/Dropbox/Projects/McMillen/Code"

# --------------------------------------------------------------------------- #

# load functions #
source("countyLines.R")
source("countyBoundaries.R")

# --------------------------------------------------------------------------- #

# 1. Boundaries
# Path to CMAP shapefiles
cmap.dir <- paste0(path.to.data,"/CMAP/")

lu.dir <- c("lu01"="LU01","lu05"="LU05")

# Load CMAP shapefiles
land.use <- list()
for (i in lu.dir){
    shp.dir <- paste0(cmap.dir,i,"/shapefile")
    shp.files <- list.files(shp.dir)
    # shp.files <- shp.files[grep("(\\.zip$|\\.pdf|\\.doc|\\.txt|\\.htm)$",shp.files,invert=TRUE)]

    layer <- shp.files[grep("\\.shp$",shp.files)]
    layer <- substr(layer,1,nchar(layer)-4)
    lu <- readOGR(shp.dir,layer)
    land.use[i] <- lu
    rm(lu)
}

# --------------------------------------------------------------------------- #
# Overlay CMAP parcels with Census Tiger Line shapefiles for high-res County boundaries
# (benchmark Census 2010, vintage Census 2000/2010)

# FIPS Codes
# State: IL=17
# Counties: Cook=031; DuPage=043; Kane=089; Kendall=093; Lake=097; McHenry=111; Will=197
ilfips <- paste0("17",c("031","043","089","093","097","111","197"))

# Path to shapefiles
census.dir <- paste0(path.to.data,"/Census/")
tiger.path <- paste0(census.dir,"TIGER Shapefiles/")

# Benchmark:
# bench <- 2010

# Vintage:
vintage <- "2000" # "2010"

tiger <- paste0(tiger.path,vintage," Tiger Shapefiles") # eg, "/2000 Tiger Shapefiles"

# select group of shapefiles {tabblock, tract, county,...}
coshp <- paste0("tabblock",substr(vintage,3,4))

# Create County Lines shapefile
colines <- countyLines(fips.list=ilfips,tl.res=coshp,tl.dir=tiger)

save(colines,file=paste0(cmap.dir,"colines",vintage,".Rda"))

# N.B.: Census not in same projection as CMAP
# Function countyBoundaries takes care of projection

cmap.01 <- countyBoundaries(fips.list=ilfips,
                            lu.file=land.use[["LU01"]],
                            tl.poly=colines[["County.Polygons"]],
                            save.path=paste0(cmap.dir,lu.dir[1],"/"),
                            save.it=TRUE)

cmap.05 <- countyBoundaries(fips.list=ilfips,
                            lu.file=land.use[["LU05"]],
                            tl.poly=colines[["County.Polygons"]],
                            save.path=paste0(cmap.dir,lu.dir[2],"/"),
                            save.it=TRUE)

# COMBINE #

# create unique row names
row.names(cmap.01) <- paste(row.names(cmap.01),"LU01",sep=".")
row.names(cmap.05) <- paste(row.names(cmap.05),"LU05",sep=".")

# create year indicator
cmap.01$Year <- as.factor(rep("2001",nrow(cmap.01)))
cmap.05$Year <- as.factor(rep("2005",nrow(cmap.05)))

# standardize names
names(cmap.01) <- names(cmap.05)

cmap <- spRbind(cmap.01,cmap.05)
cmap$ID <- row.names(cmap)

# SAVE
save(cmap,file=paste0(cmap.dir,"cmap",vintage,".Rda"))

# --------------------------------------------------------------------------- #

# 2. Distance to boundary

# Compute centroids
centroids <- gCentroid(cmap,byid=TRUE)
cmap$Longitude <- coordinates(centroids)[,1]
cmap$Latitude <- coordinates(centroids)[,2]
rm(centroids)

# Compute distance from CMAP parcels (polygons) to nearest county boundary using McSpatial::geoshape
cmap$dboundary <- geoshape(cmap$Longitude,cmap$Latitude,linefile=colines[["County.Lines"]])

# SAVE
save(cmap,file=paste0(cmap.dir,"cmap",vintage,".Rda"))

# --------------------------------------------------------------------------- #

# 3. Get census tracts and blocks

tl.shp <- list.files(tiger,pattern=paste(coshp,"shp$",sep="."))

cmap.plus <- NULL
for (fips in ilfips) {
    # Get shapefile name
    layer <- tl.shp[grep(fips,tl.shp)]
    layer <- substr(layer,1,nchar(layer)-4)

    # Import shapefile
    tl <- readOGR(tiger,layer)

    # Re-project shapefile to CMAP projection
    tl <- spTransform(tl, CRSobj = CRS(proj4string(cmap)))

    sub <- cmap[cmap$COUNTYFP==fips,]
    lmat <- SpatialPoints(sub[,c("Longitude","Latitude")],proj4string=CRS(proj4string(cmap)))
    over.pts <- over(lmat,tl)

    sub <- spCbind(sub,over.pts[,c("TRACTCE00","BLOCKCE00")])

    # first pass flag
    if (fips==ilfips[1]) {
        cmap.plus <- sub
    } else {
        cmap.plus <- spRbind(cmap.plus,sub)
        }
    # Clean up
    rm(layer,tl,sub,lmat,over.pts)
}

cmap <- cmap.plus
rm(cmap.plus)

# SAVE 
save(cmap,file=paste0(cmap.dir,"cmap",vintage,".Rda"))

# --------------------------------------------------------------------------- #

# 4. Aggregate CMAP Land use indicators

# RESIDENTIAL (1100 Series)
# COMMERCIAL AND SERVICES (1200 Series)
# INSTITUTIONAL (1300 Series)
# INDUSTRIAL, WAREHOUSING AND WHOLESALE TRADE (1400 Series)
# TRANSPORTATION, COMMUNICATION, AND UTILITIES (1500 Series)
# AGRICULTURAL LAND (Series 2000) (excluding farm houses, 1120)
# OPEN SPACE (Series 3000)
# VACANT, WETLANDS, OR UNDER CONSTRUCTION (4000 Series)
# WATER (5000 Series)
# 9999.  OUTSIDE THE 7-COUNTY REGION

# split-apply-combine
use <- as.character(c(11:15,2,3,4,5))
names(use) <- c("Residential","Commercial","Institutional","Industrial","Transportation","Agricultural","Open","Vacant","Water")

df <- NULL
for (i in 1:length(use)) {
    sub <- cmap[grep(paste0("^",use[i]),cmap$LANDUSE),]
    sub$Use <- as.factor(rep(names(use[i]),nrow(sub)))
    if (i==1) {
        df <- sub
    } else {
        df <- spRbind(df,sub)
    }
    rm(sub)
}

cmap <- df
rm(df)

# SAVE
save(cmap,file=paste0(cmap.dir,"cmap",vintage,".Rda"))

# --------------------------------------------------------------------------- #
