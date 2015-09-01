## CMAP
## Load packages, functions, etc.

## --------------------------------------------------------------------------- #

## Load Packages

## GIS
library(rgdal)
library(rgeos)
library(maptools)
library(McSpatial)

## Analysis
library(dplyr)
library(sqldf)
library(mlogit)
library(mnlogit)
library(mclogit)
library(glmnet)
library(nnet)

## Plots
library(ggplot2)
library(ggmap)
library(lattice)
library(Cairo)

## --------------------------------------------------------------------------- #

## Load custom functions

source("countyLines.R")
source("countyBoundaries.R")

## --------------------------------------------------------------------------- #

## Paths to load data
## REMARK: move data to CMAP project directory for relative paths

## path.to.data <- "./Data/"

cmap.dir <- "./Data/"

## --------------------------------------------------------------------------- #
