# Overlay CMAP parcels with Census Tiger Line shapefiles (Census 2010) for high-res County boundaries
countyBoundaries <- function(fips.list,lu.file,tl.poly,save.path=getwd(),save.it=FALSE) { 
    # fips.list = list of 5-digit FIPS codes
    # lu.file = land use file
    # tl.poly = polygon object
    # save.path = "path/to/save/"
    # if save.it=TRUE, save luXX.Rda files in save.path

    require(rgdal)
    require(rgeos)
    require(maptools)
    
    # Create empty objects
    cmap <- NULL

    # Project
    tl.poly <- spTransform(tl.poly, CRSobj = CRS(proj4string(lu.file)))
   
    for (fips in fips.list) {
        # Intersect County Boundaries with CMAP and label by FIPS
        poly <- tl.poly[grep(fips,row.names(tl.poly)),]
        
        lu.sub <- lu.file[!is.na(over(lu.file,poly)),]
        lu.sub.df <- lu.sub@data
        lu.sub.df <- data.frame(lu.sub.df,"COUNTYFP"=rep(fips,nrow(lu.sub.df)))
        # lu.sub$COUNTYFP10 <- rep(fips,nrow(lu.sub))
        lu.sub <- SpatialPolygonsDataFrame(lu.sub,lu.sub.df)

        # Create unique polygon IDs
        row.names(lu.sub) <- paste(row.names(lu.sub),fips,sep=".")
        # names.polygons <- sapply(lu.sub@polygons, function(x) paste(slot(x,"ID"),fips,sep="."))

        # first pass flag
        if (fips==fips.list[1]) { 
            cmap <- lu.sub
        } else {
            cmap <- spRbind(cmap,lu.sub)
        }
        if (save.it) save(lu.sub,file=paste0(save.path,"lu",fips,".Rda"))
        
        rm(lu.sub,lu.sub.df,poly) 
    }

    # Create new "Object ID" from row.names
    cmap$ID <- row.names(cmap)
    return(cmap)
}
