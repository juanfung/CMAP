# Aggregate shapefiles into County Lines 
# Use Census block shapefiles for high-res County boundaries

countyLines <- function(fips.list,tl.res,tl.dir=getwd()) {  
    # fips.list = 5-digit FIPS code
    # tl.res = resolution (e.g., tabblock00, tract10, etc.
    # tl.dir = "/path/to/shapefiles/", default = getwd()

    require(rgdal)
    require(rgeos)
    require(maptools)

    # List shapefiles
    tl.shp <- list.files(tl.dir,pattern=paste(coshp,"shp$",sep="."))

    # Create empty objects
    copolygons <- NULL
    copolyset <- NULL
    colines <- NULL
    
    for (fips in fips.list) {
        # Get shapefile name
        layer <- tl.shp[grep(fips,tl.shp)]
        layer <- substr(layer,1,nchar(layer)-4)
    
        # Import shapefile
        tl <- readOGR(tl.dir,layer)
            
        # Subset EXCLUDE Lake Michigan
        tract <- grep("^TRACT",names(tl),value=T)
        tl <- tl[!(eval(parse(text=paste0("tl$",tract))) %in% c("000000","990000")),]
    
        # Dissolve
        pts <- coordinates(tl)
        IDOneBin <- cut(pts[,1], range(pts[,1]), include.lowest=TRUE)
        tl <- unionSpatialPolygons(tl,IDOneBin)

        # Create unique IDs
        row.names(tl) <- paste(row.names(tl),fips,sep=".")    #spChFIDs(tl,fips)
        tl.dissolve <- SpatialPolygons2PolySet(tl)    

        # Create SpatialLines object from dissolved polygon
        tl.lines <- as(tl,"SpatialLines")

        # Creat unique IDs
        row.names(tl.lines) <- paste(row.names(tl.lines),fips,sep=".")
        # spChFIDs(tl.lines,fips)

        # Merge boundary objects        
        # first pass flag!
        if (fips==fips.list[1]) { 
            copolygons <- tl
            copolyset <- tl.dissolve
            colines <- tl.lines
        } else {
            copolygons <- spRbind(copolygons,tl)
            copolyset <- rbind(copolyset,tl.dissolve)
            colines <- spRbind(colines,tl.lines)
        }
        message("County ",fips," done.")
        # Clean up!
        rm(layer,tl,tract,pts,IDOneBin,tl.dissolve,tl.lines)
    }

    return(list("County.Lines"=colines,"County.Polygons"=copolygons,"County.Polyset"=copolyset))
}
