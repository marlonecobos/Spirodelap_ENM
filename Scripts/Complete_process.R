# ------------------------------------------------------------------------------
# Project: Macroecological factors shaping the geographic distribution of 
#          Spirodela polyrhiza, with reflections on predictor selection in 
#          ecological niche modeling
# Author: Marlon E. Cobos and A. Townsend Peterson
# Date: 28/09/2022
# ------------------------------------------------------------------------------

# Description ------------------------------------------------------------------
# This script contains code to reproduce analyses and figures for the project:
# 
# Most of the data required can be downloaded using the code in this script, if 
# other data is required, instructions on how to get them are provided in the 
# relevant section.
#
# Note: this script will produce a large amount of results which will be written
# in your working directory.
# ------------------------------------------------------------------------------

# R packages needed ------------------------------------------------------------
## ellipsenm and kuenm are packages from GitHub. Instructions to install them 
## can be found in the following links:
## https://github.com/marlonecobos/ellipsenm#installing-the-package
## https://github.com/marlonecobos/kuenm#installing-the-package
##
## Other packages can be installed as: install.packages("package_name")

library(BIEN)
library(biosurvey)
library(ellipsenm)
library(gdalUtilities)
library(maps)
library(raster)
library(rgbif)
library(rgdal)
library(rgeos)
library(scales)
library(spocc)
library(viridis)
# ------------------------------------------------------------------------------


# Preparations -----------------------------------------------------------------
# set working directory
setwd("Your/Directory")

# initial subdirectory for data
dir.create("Data")

# loading functions prepared for running some analyses
source("GITHUB/REPO/00_Functions.R")
# ------------------------------------------------------------------------------


# Obtaining occurrence data ----------------------------------------------------
# data from BIEN
sp_bien <- BIEN_occurrence_species(species = "Spirodela polyrhiza", 
                                   cultivated = FALSE, native.status = TRUE, 
                                   observation.type = TRUE, all.taxonomy = TRUE)
## references from BIEN
spbcit <- BIEN_metadata_citation(sp_bien, bibtex_file = "Data/bien_ref.bib")


# data from GBIF
sp_spocc <- occ(query = "Spirodela polyrhiza", from = "gbif", limit = 100000, 
                has_coords = TRUE)

## getting only data frame with records
sp_gbif <- sp_spocc$gbif$data$Spirodela_polyrhiza

## references from GBIF
sp_search <- occ_search(taxonKey = sp_spocc$gbif$data$Spirodela_polyrhiza$taxonKey[1])
cit <- gbif_citation(sp_search)
sink("Data/gbif_ref.txt")
sapply(cit, print)
sink()
# ------------------------------------------------------------------------------


# Obtaining raster variables ---------------------------------------------------
# Worldclim layers (v 2)
cvars <- get_NWC_bio(period = "historical", res = "10m", 
                     output_dir = "Data/Current")
names(cvars) <- gsub("wc2.1_10m_", "", names(cvars))

# Solar radiation
download.file("http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_srad.zip",
              destfile = "Data/Current/wc2.1_10m_srad.zip", method = "auto", 
              quiet = FALSE, mode = "wb", cacheOK = TRUE)

## unzip
dfol <- "Data/Current/wc2.1_10m_srad"
dir.create(dfol)
unzip(zipfile = "Data/Current/wc2.1_10m_srad.zip", exdir = dfol)

## stack of SR layers
srad <- stack(list.files(dfol, pattern = ".tif$", full.names = TRUE))


# SoilGrids   
dir.create("Data/Soils")
## CEC
download.file(url = "https://files.isric.org/soilgrids/data/recent/CECSOL_M_sl1_250m_ll.tif",
              destfile = "Data/Soils/CECSOL_M_sl1_250m_ll.tif", 
              method = "auto", quiet = FALSE, mode = "wb", cacheOK = TRUE)

## SOM
download.file(url = "https://files.isric.org/soilgrids/data/recent/ORCDRC_M_sl1_250m_ll.tif",
              destfile = "Data/Soils/ORCDRC_M_sl1_250m_ll.tif", 
              method = "auto", quiet = FALSE, mode = "wb", cacheOK = TRUE)

## pHW
download.file(url = "https://files.isric.org/soilgrids/data/recent/PHIHOX_M_sl1_250m_ll.tif",
              destfile = "Data/Soils/PHIHOX_M_sl1_250m_ll.tif", 
              method = "auto", quiet = FALSE, mode = "wb", cacheOK = TRUE)

soils <- c("Data/Soils/CECSOL_M_sl1_250m_ll.tif", "Data/Soils/ORCDRC_M_sl1_250m_ll.tif",
           "Data/Soils/PHIHOX_M_sl1_250m_ll.tif")


# Phosphorus layers
## these data were downloaded manually in netCDF format from 
## https://daac.ornl.gov/SOILS/guides/Global_Phosphorus_Dist_Map.html 
## before running the following lines, download the zip folder, place it in 
## the folder "Data/Soils/", and unzip it if needed

layers <- c("tot", "lab", "org", "occ", "sec", "apa")
nc <- "Data/Soils/GLOBAL_PHOSPHORUS_DIST_MAP_122/data/pforms_den.nc"
phos <- lapply(layers, function(x) {
  pho <- brick(nc, varname = x)
  writeRaster(pho[[1]], filename = paste0("Data/Soils/global_phos_", x, ".tif"), 
              format = "GTiff")
  pho
})

## keeping only TP, LIP, and OP
phos <- do.call(stack, phos)[[1:3]]
names(phos) <- c("TP", "LIP", "OP")
# ------------------------------------------------------------------------------


# Occurrence data cleaning -----------------------------------------------------
# cleaning records from both data bases
## excluding records with no coordinates 
sp_bienc <- sp_bien[!is.na(sp_bien$longitude) | !is.na(sp_bien$latitude), ]

sp_gbifc <- sp_gbif[!is.na(sp_gbif$longitude) | !is.na(sp_gbif$latitude), ]

## excluding records older than 1970
sp_bienc <- sp_bienc[!is.na(sp_bienc$date_collected), ]
sp_bienc <- sp_bienc[sp_bienc$date_collected >= "1970-01-01", ]

sp_gbifc <- sp_gbifc[!is.na(sp_gbifc$year), ]
sp_gbifc <- sp_gbifc[sp_gbifc$year >= 1970, ]

## excluding records marked as absent or uncommon (GBIF)
sp_gbifc <- sp_gbifc[!sp_gbifc$occurrenceStatus %in% c("absent", "Uncommon"), ]

## excluding records with 0,0 coordinates
sp_bienc <- sp_bienc[sp_bienc$longitude != 0 & sp_bienc$latitude != 0, ]

sp_gbifc <- sp_gbifc[sp_gbifc$longitude != 0 & sp_gbifc$latitude != 0, ]

## excluding exact duplicates
sp_bienc <- sp_bienc[!duplicated(paste(sp_bienc$longitude, sp_bienc$latitude)), ]

sp_gbifc <- sp_gbifc[!duplicated(paste(sp_gbifc$longitude, sp_gbifc$latitude)), ]

## keeping only relevant fields
sp_bienr <- sp_bienc[, c("scrubbed_species_binomial", "longitude", "latitude")]

sp_gbifr <- sp_gbifc[, c("species", "longitude", "latitude")]

colnames(sp_bienr) <- c("species", "longitude", "latitude")


# merging databases
sp_all <- rbind(sp_bienr, sp_gbifr)

## excluding exact duplicates
sp_all <- sp_all[!duplicated(paste(sp_all$longitude, sp_all$latitude)), ]

## detecting records out of variables and moving if distance to layer < 10 km
sp_alltc <- to_closest(sp_all, "longitude", "latitude", cvars$bio_1, 10)

sp_alltc <- sp_alltc[, 1:3]


# spatial thinning
## excluding based on distance considering pixel size of variables
sp_allthin <- thin_data(sp_alltc, "longitude", "latitude", thin_distance = 50)
# ------------------------------------------------------------------------------


# Delimitation of areas for model calibration (M) ------------------------------
# loading Ecoregions 
## ecoregions where downloaded from:
## https://www.sciencebase.gov/catalog/item/508fece8e4b0a1b43c29ca22
## before running the following lines, download the zip folder, place it in 
## the folder "Data/", and unzip it if needed
ecoregions <- readOGR(dsn = "Data/TerrestrialEcos", layer = "wwf_terr_ecos")

# M hypothesis (buffers 500 km)
## new folder for Ms
dir.create("Shapefiles")

## buffer
buff <- buffer_area(data = sp_allthin, longitude = "longitude", 
                    latitude = "latitude", buffer_distance = 500, 
                    save = TRUE, name = "Shapefiles/buffer_area")

## concave polygons
conc <- concave_area(data = sp_allthin, longitude = "longitude", 
                     latitude = "latitude", length_threshold = 100, split = TRUE, 
                     cluster_method = "hierarchical", split_distance = 2500, 
                     buffer_distance = 500)
conc <- suppressWarnings(disaggregate(gBuffer(conc, width = 0)))
conc <- SpatialPolygonsDataFrame(conc, data = data.frame(RD = rep(1, length(conc))),
                                 match.ID = FALSE)

save_areas(name = "Shapefiles/concave_area", area_polygon = conc, 
           area_type = "concave_area")

## ecoregions (only 100 km to consider areas explored but not colonized)
ecor <- polygon_selection(data = sp_allthin, longitude = "longitude", 
                          latitude = "latitude", polygons = ecoregions)

un <- ifelse(is.na(ecor@data$REALM), "LA", as.character(ecor@data$REALM))
ecor <- gUnaryUnion(ecor, id = un, checkValidity = 2L)
ecor <- suppressWarnings(disaggregate(gBuffer(ecor, width = 100/111.32)))
ecor <- SpatialPolygonsDataFrame(ecor, data = data.frame(RD = rep(1, length(ecor))),
                                 match.ID = FALSE)

save_areas(name = "Shapefiles/ecoregions", area_polygon = ecor, 
           area_type = "polygon_selection")

## intersection of all previous hypothesis
inter <- gIntersection(buff, conc, byid = TRUE, drop_lower_td = TRUE)
inter <- gIntersection(inter, ecor, byid = TRUE, drop_lower_td = TRUE)
inter <- SpatialPolygonsDataFrame(inter, data = data.frame(RD = rep(1, length(inter))),
                                  match.ID = FALSE)

save_areas(name = "Shapefiles/intersection", area_polygon = inter, 
           area_type = "M_intersection")

## union of four M areas for a final projection (not an M hypothesis)
bcei <- gUnion(gUnion(buff, conc), gUnion(ecor, inter))
# ------------------------------------------------------------------------------


# Figure with all calibration areas and occurrences ----------------------------
marg <- rep(0, 4)

jpeg(file = paste0("Figures/Calibration_areas.jpg"), 
     width = 120, height = 200, units = "mm", res = 600)

par(mfrow = c(4, 1))

map(col = "gray85", mar = marg)
plot(buff, add = TRUE)
points(sp_allthin[, 2:3], col = alpha("gray15", 0.4), pch = 16, cex = 0.4)
legend("bottomleft", legend = "Buffer", bg = "white", box.col = "white")
box()

map(col = "gray85", mar = marg)
plot(conc, add = TRUE)
legend("bottomleft", legend = "Concave", bg = "white", box.col = "white")
box()

map(col = "gray85", mar = marg)
plot(ecor, add = TRUE)
legend("bottomleft", legend = "Ecoregions", bg = "white", box.col = "white")
box()

map(col = "gray85", mar = marg)
plot(inter, add = TRUE)
legend("bottomleft", legend = "Intersection", bg = "white", box.col = "white")
box()

dev.off()
# ------------------------------------------------------------------------------


# Variable preparation (I) -----------------------------------------------------
# solar radiation layer preparation
## getting values
rdvals <- values(srad)

## creating annual mean solar radiation
srmean <- srad[[1]]
srmean[] <- apply(rdvals, 1, mean)

## creating minimum solar radiation
srmin <- srad[[1]]
srmin[] <- apply(rdvals, 1, min)

## creating maximum solar radiation
srmax <- srad[[1]]
srmax[] <- apply(rdvals, 1, max)

## creating annual range solar radiation
srran <- srad[[1]]
srran <- srmax - srmin

## creating values to obtain metrics by quarters
srqts <- t(apply(rdvals, 1, window)) / 3

## creating mean for quarter with highest values
srqh <- srad[[1]]
srqh[] <- apply(srqts, 1, max)

## creating mean for quarter with lowest values
srql <- srad[[1]]
srql[] <- apply(srqts, 1, min)

## all solar variables
srall <- stack(srmean, srmax, srmin, srran, srqh, srql)
names(srall) <- c("AMSR", "SRMax", "SRMin", "RSR", "ASRQH", "ASRQL")

## writing layers
srdir <- "Data/Current/WC_solar_rad_sum"
dir.create(srdir)

wr <- lapply(1:nlayers(srall), function(x) {
  writeRaster(srall[[x]], filename = paste0(srdir, "/", names(srall)[x], ".tif"), 
              format = "GTiff")
})
# ------------------------------------------------------------------------------


# Variable preparation (II) ----------------------------------------------------
# variable aggregation to coarser resolutions
## bioclim to phosphorus resolution
### existent files
files <- list.files(path = "Data/Current/wc2.1_10m_bio/", pattern = ".tif$", 
                    full.names = TRUE)

### projection
dir.create("Data/WC_pr")
gfiles <- gsub("/Current/wc2.1_10m_bio/wc2.1_10m_", "/WC_pr/", files)

repro <- sapply(1:length(files), function(x) {
  gdal_translate(files[x], gfiles[x], a_srs = "EPSG:4326")
})

### new files
ndir <- "Data/WC_agg/"
dir.create(ndir)
nfiles <- gsub("/WC_pr/", "/WC_agg/", gfiles)

### aggregation
aggre <- sapply(1:length(files), function(x) {
  gdalwarp(gfiles[x], nfiles[x],r = "average", multi = TRUE, tr = res(phos$OC))
})


## solar radiation to phosphorus resolution
### existent files
files <- list.files(path = "Data/Current/WC_solar_rad_sum/", pattern = ".tif$", 
                    full.names = TRUE)

### projection
dir.create("Data/Current/WC_sr_pr")
gfiles <- gsub("/Current/WC_solar_rad_sum/", "/WC_sr_pr/", files)

repro <- sapply(1:length(files), function(x) {
  gdal_translate(files[x], gfiles[x], a_srs = "EPSG:4326")
})

### new files
ndir <- "Data/WC_sr_agg/"
dir.create(ndir)
nfiles <- gsub("/WC_sr_pr/", "/WC_sr_agg/", gfiles)

### aggregation
aggre <- sapply(1:length(files), function(x) {
  gdalwarp(gfiles[x], nfiles[x],r = "average", multi = TRUE, tr = res(phos$OC))
})


## SoilGrids to bioclim resolution
files <- list.files(path = "Data/Soils/", pattern = "ll.tif$", 
                    full.names = TRUE)

### projection
dir.create("Data/Soils_pr")
gfiles <- gsub("/Soils/", "/Soils_pr/", files)

repro <- sapply(1:length(files), function(x) {
  gdal_translate(files[x], gfiles[x], a_srs = "EPSG:4326")
})

### new files
ndir <- "Data/Soils_agg/"
dir.create(ndir)
nfiles <- gsub("/Soils/", "/Soils_agg/", files)

### aggregation
aggres <- sapply(1:length(files), function(x) {
  gdalwarp(gfiles[x], nfiles[x],r = "average", multi = TRUE, 
           tr = res(cvars$bio_1))
})
# ------------------------------------------------------------------------------


# Variable preparation (III) ---------------------------------------------------
# masking layers
## re-read aggregated layers
wc05 <- stack(list.files("Data/WC_agg/", pattern = ".tif$", full.names = TRUE))
sr05 <- stack(list.files("Data/WC_sr_agg/", pattern = ".tif$", full.names = TRUE))
sg10 <- stack(list.files("Data/Soils_agg/", pattern = ".tif$", full.names = TRUE))

## fixing names
names(sg10) <- c("CEC", "OC", "pH")
names(wc05) <- toupper(names(wc05))
names(cvars) <- toupper(names(cvars))

## making to consensus for all variables at each resolution
### 0.5 degrees
mskp05 <- !is.na(phos$LIP)
mskw05 <- !is.na(wc05$BIO_1)
msksr05 <- !is.na(sr05$AMSR)
msk05 <- mskp05 * mskw05 * msksr05
msk05[msk05[] == 0] <- NA

wc05m <- mask(wc05, msk05)
sr05m <- mask(sr05, msk05)
phosm <- mask(phos, msk05)

var05 <- stack(wc05m, sr05m, phosm)

### 10 min
msks10 <- !is.na(sg10$CEC)
mskw10 <- !is.na(cvars$BIO_1)
msksr10 <- !is.na(srall$AMSR)
msk10 <- msks10 * mskw10 * msksr10 # the warning is not a problem
msk10[msk10[] == 0] <- NA 

wc10 <- cvars * msk10 # the warning is not a problem
wc10m <- mask(wc10, msk10)
names(wc10m) <- names(cvars)

sr10 <- srall * msk10 # the warning is not a problem
sr10m <- mask(sr10, msk10)
names(sr10m) <- names(srall)

sg10m <- mask(sg10, msk10)

var10 <- stack(wc10m, sr10m, sg10m)

## writing variables masked to consensus
dir05 <- "Data/Layers_05_consensus"
dir10 <- "Data/Layers_10_consensus"

dir.create(dir05)
dir.create(dir10)

nms05 <- paste0(dir05, "/", names(var05), ".tif")
nms10 <- paste0(dir10, "/", names(var10), ".tif")

for (i in 1:nlayers(var05)) {
  writeRaster(var05[[i]], filename = nms05[i], format = "GTiff")
  writeRaster(var10[[i]], filename = nms10[i], format = "GTiff")
}


## list all M options
all_ms <- list(buffer = buff, concave = conc, ecorregions = ecor, 
               intersection = inter)

## masking to selected M (warnings from projections don't matter)
vmask <- lapply(all_ms, function(x) {
  list(variables_05 = mask(crop(var05, x), x), 
       variables_10 = mask(crop(var10, x), x))
})

## masking to union of Ms
v05_bcei <- mask(crop(var05, bcei), bcei)
v10_bcei <- mask(crop(var10, bcei), bcei)

## separating variables depending on resolution and M
variables05 <- lapply(vmask, function(x) x$variables_05)
variables10 <- lapply(vmask, function(x) x$variables_10)

## sorting variables
sortv05 <- c(1, 12:19, 2:11, 20, 24:25, 23, 21:22, 26:28)
sortv10 <- c(1, 12:19, 2:11, 20:28)
variables05 <- lapply(variables05, function(x) x[[sortv05]])
variables10 <- lapply(variables10, function(x) x[[sortv10]])

## writing raster layers 
### masked to four options of M
dir05 <- "Data/Layers_05_M"
dir10 <- "Data/Layers_10_M"

dir.create(dir05)
dir.create(dir10)

sfols05 <- paste0(dir05, "/", names(variables05))
sfols10 <- paste0(dir10, "/", names(variables10))

for (i in 1:length(variables05)) {
  dir.create(sfols05[i])
  names05 <- names(variables05[[i]])
  for (j in 1:nlayers(variables05[[i]])) {
    writeRaster(variables05[[i]][[j]], 
                filename = paste0(sfols05[i], "/", names05[j]), 
                format = "GTiff")
  }
  
  dir.create(sfols10[i])
  names10 <- names(variables10[[i]])
  for (j in 1:nlayers(variables10[[i]])) {
    writeRaster(variables10[[i]][[j]], 
                filename = paste0(sfols10[i], "/", names10[j]), 
                format = "GTiff")
  }
}

### masked to union of Ms
dir.create("Data/Layers_05_unionM")
dir.create("Data/Layers_10_unionM")

nms05um <- gsub("consensus", "unionM", nms05[i])
nms10um <- gsub("consensus", "unionM", nms10[i])

for (i in 1:nlayers(v05_bcei)) {
  writeRaster(v05_bcei[[i]], filename = nms05um, format = "GTiff")
  writeRaster(v10_bcei[[i]], filename = nms10um, format = "GTiff")
}
# ------------------------------------------------------------------------------


# Last step of data cleaning considering variables -----------------------------
# testing NAs in occurrence data after aggregation and masking
v <- lapply(vmask, function(x) {
  lapply(x, extract, sp_allthin[, 2:3], cellnumbers = T)
})
vt <- lapply(v, function(x) {
  lapply(x, function(y) {lapply(1:ncol(y), function(z) which(is.na(y[, z])))})
})

out05 <- unlist(lapply(vt, function(x) unique(unlist(x$variables_05))))
out05u <- unique(out05)

out10 <- unlist(lapply(vt, function(x) unique(unlist(x$variables_10))))
out10u <- unique(out10)

# excluding data out of layers (NAs) depending on resolution
sp_ready05 <- sp_allthin[-out05u, ]
sp_ready10 <- sp_allthin[-out10u, ]
# ------------------------------------------------------------------------------


# Exploring variables (I) ------------------------------------------------------
# variable values in M and occurrences
mvalues05 <- lapply(variables05, function(x) {
  m <- data.frame(background = "background", rasterToPoints(x))
  set.seed(1)
  m[sample(nrow(m), 10000), ]
})
mvalues10 <- lapply(variables10, function(x) {
  m <- data.frame(background = "background", rasterToPoints(x))
  set.seed(1)
  m[sample(nrow(m), 10000), ]
})

ovalues05 <- lapply(variables05, function(x) {
  data.frame(sp_ready05, extract(x, sp_ready05[, 2:3]))
})
ovalues10 <- lapply(variables10, function(x) {
  data.frame(sp_ready10, extract(x, sp_ready10[, 2:3]))
})

# exploring distribution of variable values in M and in occurrences
## calculating values
hist05 <- lapply(1:length(mvalues05), function(x) {
  vars <- lapply(4:ncol(mvalues05[[x]]), function(y) {
    M <- mvalues05[[x]][, y]; O <- ovalues05[[x]][, y]; mo <- c(M, O)
    ax <- pretty((min(mo) - 0.1):max(mo), n = 12) 
    df <- diff(ax[1:2])
    ax <- c(ax[1] - df, ax, ax[length(ax)] + df)
    list(hgM = hist(M, breaks = ax, plot = FALSE),
         hgO = hist(c(O, O, O), breaks = ax, plot = FALSE),
         meanM = mean(M), meanO = mean(O),
         clM = quantile(M, c(0.05, 0.95)),
         clO = quantile(O, c(0.05, 0.95))) 
  })
  names(vars) <- colnames(mvalues05[[x]])[4:ncol(mvalues05[[x]])]
  return(vars)
})
names(hist05) <- names(mvalues05)

hist10 <- lapply(1:length(mvalues10), function(x) {
  vars <- lapply(4:ncol(mvalues10[[x]]), function(y) {
    M <- mvalues10[[x]][, y]; O <- ovalues10[[x]][, y]; mo <- c(M, O)
    ax <- pretty((min(mo) - 0.1):max(mo), n = 12) 
    df <- diff(ax[1:2])
    ax <- c(ax[1] - df, ax, ax[length(ax)] + df)
    list(hgM = hist(M, breaks = ax, plot = FALSE),
         hgO = hist(c(O, O, O), breaks = ax, plot = FALSE),
         meanM = mean(M), meanO = mean(O),
         clM = quantile(M, c(0.05, 0.95)),
         clO = quantile(O, c(0.05, 0.95))) 
  })
  names(vars) <- colnames(mvalues10[[x]])[4:ncol(mvalues10[[x]])]
  return(vars)
})
names(hist10) <- names(mvalues10)


## plotting results
yn <- c(1, 6, 11, 16, 21, 26); pmar <- c(4.5, 4.3, 0.8, 0.6) 
pmf <- c(6, 5); pcex <- 0.4
mmcol <- "gray65"; mocol <- "gray25"
mcol <- alpha(mmcol, 0.4); ocol <- alpha(mocol, 0.4)
hnames <- names(hist10)

ps05 <- lapply(1:length(hnames), function(x) {
  #pdf(file = paste0("Figures/", hnames[x] , "_histograms_V05d.pdf"))
  jpeg(file = paste0("Figures/", hnames[x] , "_histograms_V05d.jpg"), 
       width = 166, height = 166, units = "mm", res = 600)
  plot_list_histMO(lits_histMO = hist05[[x]], par_mfrow = pmf, 
                   where_ylab = yn, par_cex = pcex, par_mar = pmar)
  par(mar = rep(0.5, 4))
  plot.new()
  legend("center", legend = c("Variable values in M (N = 10,000)", 
                              "Values in records (Frequency x 3)"),
         fill = c(mcol, ocol), border = c(mcol, ocol), bty = "n", cex = 1.1)
  plot.new()
  legend("center", legend = c("Mean of values in M", "Mean of values in records",
                              "95% CL for M", "95% CL for records"), bty = "n", 
         col = rep(c(mmcol, mocol), 2), lty = c(1, 1, 2, 2), cex = 1.1)
  dev.off()
})

ps10 <- lapply(1:length(hnames), function(x) {
  #pdf(file = paste0("Figures/", hnames[x] , "_histograms_V10m.pdf"))
  jpeg(file = paste0("Figures/", hnames[x] , "_histograms_V10m.jpg"), 
       width = 166, height = 166, units = "mm", res = 600)
  plot_list_histMO(lits_histMO = hist10[[x]], par_mfrow = pmf, 
                   where_ylab = yn, par_cex = pcex, par_mar = pmar)
  par(mar = rep(0.5, 4))
  plot.new()
  legend("center", legend = c("Variable values in M (N = 10,000)", 
                              "Values in records (Frequency x 3)"),
         fill = c(mcol, ocol), border = c(mcol, ocol), bty = "n", cex = 1.1)
  plot.new()
  legend("center", legend = c("Mean of values in M", "Mean of values in records",
                              "95% CL for M", "95% CL for records"), bty = "n", 
         col = rep(c(mmcol, mocol), 2), lty = c(1, 1, 2, 2), cex = 1.1)
  dev.off()
})
# ------------------------------------------------------------------------------


# Exploring variables (II) -----------------------------------------------------
# testing correlation of variables (correlation limit |0.8|)
dir.create("Results")
indir <- "Results/Correlation/"
dir.create(indir)
nc <- ncol(mvalues05$buffer)

jpeg(file = paste0("Figures/Variable_correlarion_V05d.jpg"), 
     width = 166, height = 166, units = "mm", res = 600)
par(mfrow = c(2, 2))
par(cex = 0.5)
cors05 <- lapply(1:length(hnames), function(x) {
  vals <- as.matrix(mvalues05[[x]][, 4:nc])
  colnames(vals)[20:28] <- c("AMSR", "SRMax", "SRMin", "RSR", "ASRQH", "ASRQL",
                             "TP", "LIP", "OP")
  cr05 <- variable_correlation(variables = vals, 
                               correlation_limit = 0.8, corrplot = TRUE, 
                               magnify_to = 3, tl.cex = 0.65, save = TRUE, 
                               name = paste0(indir, hnames[1], "_correlation_V05d"))
  legend(0, 4, legend = paste0(hnames[x], " M"), bty = "n")
  cr05
})
dev.off()


jpeg(file = paste0("Figures/Variable_correlarion_V10m.jpg"), 
     width = 166, height = 166, units = "mm", res = 600)
par(mfrow = c(2, 2))
par(cex = 0.5)
cors10 <- lapply(1:length(hnames), function(x) {
  vals <- as.matrix(mvalues10[[x]][, 4:nc])
  colnames(vals)[20:28] <- c("AMSR", "SRMax", "SRMin", "RSR", "ASRQH", "ASRQL",
                             "CEC", "OC", "pH")
  cr10 <- variable_correlation(variables = vals, 
                               correlation_limit = 0.8, corrplot = TRUE, 
                               magnify_to = 3, tl.cex = 0.65, save = TRUE, 
                               name = paste0(indir, hnames[1], "_correlation_V10m"))
  legend(0, 4, legend = paste0(hnames[x], " M"), bty = "n")
  cr10
})
dev.off()

names(cors05) <- hnames
names(cors10) <- hnames
# ------------------------------------------------------------------------------


# Initial selection of variables -----------------------------------------------
# excluding BIOs 8, 9, 18, 19 to reduce dimensions and avoid artifacts, as they 
# do not look particularly interesting in histograms
nams05 <- names(variables05$buffer)
nams10 <- names(variables10$buffer)
exbios <- paste0("BIO_", c(8, 9, 18, 19))

variables05 <- lapply(variables05, function(x) x[[which(!nams05 %in% exbios)]])
variables10 <- lapply(variables10, function(x) x[[which(!nams10 %in% exbios)]])

# excluding other variables
## keeping Labile inorganic P and organic P at 0.5 deg
## CEC, pH are and OC are not correlated
## given the distribution and active period of the species we want to include 
## quarterly variables to identify warmer periods, but also variables that help
## to represent environmental limits like minimum or maximum, this last were 
## prioritized
## we excluded variables so all M options work with the same initial
## non correlated variables 
nams05 <- names(variables05$buffer)
nams10 <- names(variables10$buffer)
ex05 <- c(paste0("BIO_", c(1, 3, 4, 7, 10, 11, 13, 16, 17)), "AMSR", "SRMax",
          "SRMin", "ASRQL", "TP")
ex10 <- c(paste0("BIO_", c(1, 3, 4, 7, 10, 11, 13, 16, 17)), "AMSR", "SRMax",
          "SRMin", "ASRQL")

variables05 <- lapply(variables05, function(x) x[[which(!nams05 %in% ex05)]])
variables10 <- lapply(variables10, function(x) x[[which(!nams10 %in% ex10)]])
# ------------------------------------------------------------------------------


# Preparing data for model calibration -----------------------------------------
# writing data ready for calibrations
detach("package:ellipsenm") # to avoid confusion in random splitting of records

dir05 <- "Var_selection05"
dir.create(dir05)

dir10 <- "Var_selection10"
dir.create(dir10)

sfols05 <- paste0(dir05, "/", names(variables05))
sfols10 <- paste0(dir10, "/", names(variables10))

for (i in 1:length(variables05)) {
  all_data05 <- prepare_swd(occ = sp_ready05, species = "species", 
                            longitude = "longitude", latitude = "latitude", 
                            data.split.method = "random", train.proportion = 0.5, 
                            raster.layers = variables05[[i]], sample.size = 20000, 
                            var.sets = "all_comb", min.number = 2, save = TRUE, 
                            name.occ = "Var_selection05/occ", 
                            back.folder = sfols05[i])
  
  all_data10 <- prepare_swd(occ = sp_ready10, species = "species", 
                            longitude = "longitude", latitude = "latitude", 
                            data.split.method = "random", train.proportion = 0.5, 
                            raster.layers = variables10[[i]], sample.size = 20000, 
                            var.sets = "all_comb", min.number = 2, save = TRUE, 
                            name.occ = "Var_selection10/occ", 
                            back.folder = sfols10[i])
}
# ------------------------------------------------------------------------------


# Model calibration using GLMs -------------------------------------------------
all_mnames <- c("buffer", "concave", "ecorregions", "intersection")

# data
## 0.5 degrees
adata05 <- read.csv("Var_selection05/occ_joint.csv")[, -1]
trdata05 <- read.csv("Var_selection05/occ_train.csv")[, -1]
tedata05 <- read.csv("Var_selection05/occ_test.csv")[, -1]

back05 <- paste0("Var_selection05/", all_mnames, "/Set_1013.csv")

## 10 min
adata10 <- read.csv("Var_selection10/occ_joint.csv")[, -1]
trdata10 <- read.csv("Var_selection10/occ_train.csv")[, -1]
tedata10 <- read.csv("Var_selection10/occ_test.csv")[, -1]

back10 <- paste0("Var_selection10/", all_mnames, "/Set_2036.csv")


# variable selection according to model performance (GLM)
## model calibration
rt <- c("lq", "lp", "q", "qp", "lqp")

### 0.5 degrees
od05 <- paste0("Var_selection05/cal_glm_", all_mnames)

glm05 <- list()

for (i in 1:length(back05)) {
  back <- read.csv(back05[i])[, -1]
  glm05[[i]] <- glm_calibration(all_data = adata05, train_data = trdata05, 
                                test_data = tedata05, background = back, 
                                variable_sets = "all_comb", min_number = 2,
                                response_types = rt, output_directory = od05[i])
  
}

names(glm05) <- all_mnames

### 10 min
od10 <- paste0("Var_selection10/cal_glm_", all_mnames)

glm10 <- list()

for (i in 3) {
  back <- read.csv(back10[i])[, -1]
  glm10[[i]] <- glm_calibration(all_data = adata10, train_data = trdata10, 
                                test_data = tedata10, background = back, 
                                variable_sets = "all_comb", min_number = 2,
                                response_types = rt, output_directory = od10[i])
  
}

names(glm10) <- all_mnames
# ------------------------------------------------------------------------------


# Model calibration using Maxent -----------------------------------------------
# variable selection according to model performance (Maxent, kuenm)
## model calibration
### 0.5 degrees
oj <- "Var_selection05/occ_joint.csv"
otr <- "Var_selection05/occ_train.csv"
ote <- "Var_selection05/occ_test.csv"
backs <- paste0("Var_selection05/", all_mnames)
bt <-  paste0("Var_selection05/batchc", all_mnames)
odc <-  paste0("Var_selection05/cand", all_mnames)
oe <- paste0("Var_selection05/c_res_", all_mnames)

rm <- c(0.1, 0.3, 0.6, 1, 2.5, 5)
fc <- c("lq", "lp", "q", "qp", "lqp")
mx <- "C:/Maxent/3.4.1"

for (i in 1:length(backs)) {
  kuenm_cal_swd(occ.joint = oj, occ.tra = otr, occ.test = ote, 
                back.dir = backs[i], batch = bt[i], out.dir.models = odc[i], 
                reg.mult = rm, f.clas = fc, maxent.path = mx, kept = T, 
                out.dir.eval = oe[i])
}

### 10 min
oj1 <- "Var_selection10/occ_joint.csv"
otr1 <- "Var_selection10/occ_train.csv"
ote1 <- "Var_selection10/occ_test.csv"
backs1 <- paste0("Var_selection10/", all_mnames)
bt1 <- paste0("Var_selection10/batchc", all_mnames)
odc1 <- paste0("Var_selection10/cand", all_mnames)
oe1 <- paste0("Var_selection10/c_res_", all_mnames)

for (i in 1:length(backs1)) {
  kuenm_cal_swd(occ.joint = oj1, occ.tra = otr1, occ.test = ote1, 
                back.dir = backs1[i], batch = bt1[i], out.dir.models = odc1[i], 
                reg.mult = rm, f.clas = fc, maxent.path = mx, kept = T, 
                out.dir.eval = oe1[i])
}
# ------------------------------------------------------------------------------


# Evaluating effects of variables in selected GLMs -----------------------------
### 0.5 degs
formulas_selected05 <- lapply(glm05, function(x) {x$selected_models$Formula})

effects05 <- lapply(1:length(formulas_selected05), function(x) {
  back <- read.csv(back05[x])[, -1]
  glm_enm(formula = formulas_selected05[[x]], all_data = adata05, 
          background = back, return_prediction = FALSE)
})
names(effects05) <- names(glm05)

#### summarizing results
names_e05 <- paste0("Var_selection05/predictor_effects05_", 
                    names(effects05), ".csv")

weff05 <- lapply(1:length(effects05), function(x) {
  eff <- as.data.frame(effects05[[x]]$effects)
  eff <- eff[order(eff$Deviance, na.last = FALSE, decreasing = TRUE), ]
  
  write.csv(eff, file = names_e05[x], row.names = TRUE)
})

save(glm05, formulas_selected05, effects05, 
     file = "Var_selection05/GLM05_results.RData")

### 10 min 
formulas_selected10 <- lapply(glm10, function(x) {x$selected_models$Formula})

effects10 <- lapply(1:length(formulas_selected10), function(x) {
  back <- read.csv(back10[x])[, -1]
  glm_enm(formula = formulas_selected10[[x]], all_data = adata10, 
          background = back, return_prediction = FALSE)
})
names(effects10) <- names(glm10)

#### summarizing results
names_e10 <- paste0("Var_selection10/predictor_effects10_", names(effects10), 
                    ".csv")

weff10 <- lapply(1:length(effects10), function(x) {
  eff <- as.data.frame(effects10[[x]]$effects)
  eff <- eff[order(eff$Deviance, na.last = FALSE, decreasing = TRUE), ]
  
  write.csv(eff, file = names_e10[x], row.names = TRUE)
})

save(glm10, formulas_selected10, effects10, 
     file = "Var_selection10/GLM10_results.RData")
# ------------------------------------------------------------------------------


# Evaluating effects of variables in selected Maxent models --------------------
### 0.5 degrees
occs05 <- read.csv(oj)
var_eff05 <- list()

for (i in 1:length(backs)) {
  selected <- read.csv(paste0(oe[i], "/selected_models.csv"))[1, 1]
  sel <- strsplit(selected, "_")[[1]]
  rm <- sel[2]
  fc <- sel[4]
  
  back <- read.csv(paste0(backs[i], "/", sel[5], "_", sel[6], ".csv"))
  vars <- colnames(back)[-(1:3)]
  whichb <- c(1:3, which(colnames(occs05) %in% vars))
  
  var_eff05[[i]] <- explore_var_contrib(occ = occs05[, whichb], 
                                        M_variables = back, maxent.path = mx, 
                                        reg.mult = rm, f.clas = fc, plot = F)
}

names(var_eff05) <-  all_mnames


### 10
occs10 <- read.csv(oj1)
var_eff10 <- list()

for (i in 1:length(backs1)) {
  selected <- read.csv(paste0(oe1[i], "/selected_models.csv"))[1, 1]
  sel <- strsplit(selected, "_")[[1]]
  rm <- sel[2]
  fc <- sel[4]
  
  back <- read.csv(paste0(backs1[i], "/", sel[5], "_", sel[6], ".csv"))
  vars <- colnames(back)[-(1:3)]
  whichb <- c(1:3, which(colnames(occs10) %in% vars))
  
  var_eff10[[i]] <- explore_var_contrib(occ = occs10[, whichb],
                                        M_variables = back, maxent.path = mx, 
                                        reg.mult = rm, f.clas = fc, plot = F)
}

names(var_eff10) <-  all_mnames


## summarizing results
### 0.5 degrees
pdf("Figures/Maxent_var_effects_05.pdf", width = 7, height = 3)
plot_contribution(var_eff05$buffer)
plot_contribution(var_eff05$concave)
plot_contribution(var_eff05$ecorregions)
plot_contribution(var_eff05$intersection)
dev.off()

### 10
pdf("Figures/Maxent_var_effects_10.pdf", width = 7, height = 2)
plot_contribution(var_eff10$buffer)
plot_contribution(var_eff10$concave)
plot_contribution(var_eff10$ecorregions)
plot_contribution(var_eff10$intersection)
dev.off()
# ------------------------------------------------------------------------------


# Figure summarizing all variables used in distinct processes ------------------
# colors for algorithms
## #1b9e77 = 1 = GLM; #7570b3 = 2 = Maxent; #d95f02 = 3 = Both

# codes for calibration areas
## 1 = Buffer; 2 = Concave; 3 = Ecoregions; 4 = Intersection

# codes for variables
## BIO2 = 1, BIO5 = 2, BIO6 = 3; BIO12 = 4; BIO14 = 5; BIO15 = 6; RSR = 7; 
## ASRQH = 8

## 0.5 degrees
## LIP = 9; OP = 10

## 10 minutes
## CEC = 9; OC = 10; pH = 11

# preparing data for plots
## 0.5 degrees
xval <- c(rep(1, 3), rep(2, 3), rep(3, 4), rep(4, 4), rep(5, 4), rep(6, 3), 
          rep(7, 4), rep(8, 4), 9, rep(10, 4), 11)
yval <- c(1:3, 1:2, 3, 1:4, 1:4, 1:4, 1, 3:4, 1:4, 1:4, 3, 1:4, 0) # 0 to add space
xy05 <- cbind(Variables = xval, Area = yval)
af05 <- as.factor(c(rep(1, 3), 1, 2, 1, rep(3, 10), 2, 3, 1, 1, 1, rep(3, 4),
                    1, 3, 3, 3, rep(2, 5)))

## 10 minutes
xval <- c(rep(1, 4), rep(2, 4), rep(3, 3), rep(4, 3), rep(5, 4), 6, 7, 
          rep(8, 4), rep(9, 4), rep(10, 3), rep(11, 4))
yval <- c(1:4, 1:4, 2:4, 1:2, 4, 1:4, 3, 1, 1:4, 1:4, 1, 3:4, 1:4)
xy10 <- cbind(Variables = xval, Area = yval)
af10 <- as.factor(c(rep(1, 4), 3, 3, 1, 3, rep(1, 7), 3, rep(1, 4), 3, 3, 2, 3, 
                    3, 1, 1, 3, 2, 1, rep(2, 5)))

# plotting
varnames05 <- c(paste("BIO", c(2, 5, 6, 12, 14, 15)), "RSR", "ASRQH", "LIP", 
                "OP")
varnames10 <- c(paste("BIO", c(2, 5, 6, 12, 14, 15)), "RSR", "ASRQH", "CEC", 
                "OC", "pH")
areanames <- c("Buffer", "Concave", "Ecoregions", "Intersection")
cols <- c("#b2df8a", "#a6cee3", "#1f78b4")
pchs <- 21
ptcex <- 2.5


jpeg(file = paste0("Figures/Summary_variables_selected2.jpg"), 
     width = 166, height = 150, units = "mm", res = 600)

layout(matrix(1:3, nrow = 3), heights = c(1, 10, 10))
par(cex = 0.73, mar = rep(0, 4))

plot.new()
legend(0.4, 0.79, legend = c("GLM", "Maxent", "Both"), pt.bg = cols, horiz = T, 
       bty = "n", pch = pchs, pt.cex = ptcex, inset = -0.04)

par(mar = c(4.5, 7.5, 1, 0.5))

plot(xy10, ylim = c(0.6, 4.4), xlim = c(0.85, 11.15), bty = "l", axes = F, xlab = "", 
     ylab = "", pch = pchs, bg = cols[af10], cex = ptcex)
abline(h = seq(1.5, 4.5, 1), col = "gray85")
abline(v = seq(1.5, 11.5, 1), col = "gray85")
axis(side = 1, at = 1:11, labels = varnames10)
axis(side = 2, at = 1:4, areanames, las = 2)
par(cex = 0.85); title(ylab = "Calibration areas (10 min)", line = 6.2)
par(cex = 0.73); box(bty = "l")

plot(xy05, ylim = c(0.6, 4.4), xlim = c(0.85, 11.15), bty = "l", axes = F, 
     xlab = "", ylab = "", pch = pchs, bg = cols[af05], cex = ptcex)
abline(h = seq(1.5, 4.5, 1), col = "gray85")
abline(v = seq(1.5, 11.5, 1), col = "gray85")
axis(side = 1, at = 1:10, labels = varnames05)
axis(side = 2, at = 1:4, areanames, las = 2)
par(cex = 0.85); title(xlab = "Variables")
title(ylab = "Calibration areas (30 min)", line = 6.2)
box(bty = "l")

dev.off()
# ------------------------------------------------------------------------------


# Geographic projection and comparison of best GLM models ----------------------
## all calibration exercises resulted in only one best model

## G variables (variables for geographic projection)

## 0.5 degrees
glm_gpred05 <- lapply(1:length(effects05), function (x) {
  pr <- predict_glm_enm(glm_enm = effects05[[x]], new_data = v05_bcei)
  
  indir <- paste0("Var_selection05/F_GLMmodel_", all_mnames[x])
  dir.create(indir)
  writeRaster(pr, filename = paste0(indir, "/pred_world.tif"), format = "GTiff")
  pr
})
names(glm_gpred05) <- all_mnames



## 10 min
glm_gpred10 <- lapply(1:length(effects10), function (x) {
  pr <- predict_glm_enm(glm_enm = effects10[[x]], new_data = v10_bcei)
  
  indir <- paste0("Var_selection10/F_GLMmodel_", all_mnames[x])
  dir.create(indir)
  writeRaster(pr, filename = paste0(indir, "/pred_world.tif"), format = "GTiff")
  pr
}) 
names(glm_gpred10) <- all_mnames
# ------------------------------------------------------------------------------


# Geographic projection and comparison of best Maxent models -------------------
## buffer areas resulted in 3 best models, only manually erased the second and 
## third and kept the best-best one

## preparing G variables 
dir.create("Var_selection05/G_variables")
dir.create("Var_selection10/G_variables")

### 0.5 degrees
for (i in 1:length(backs)) {
  sel <- read.csv(paste0(oe[i], "/selected_models.csv"), stringsAsFactors = F)
  sets <- strsplit(sel[, 1], split = "_")[[1]]
  
  set <- paste(sets[5:length(sets)], collapse = "_")
  bset <- read.csv(paste0(backs[i], "/", set, ".csv")) 
  gvs <- colnames(bset)[-(1:3)]
  
  dir.create(paste0("Var_selection05/G_variables/", set))
  inin <- paste0("Var_selection05/G_variables/", set, "/world")
  dir.create(inin)
  
  v05g <- v05_bcei[[gvs]] 
  nms05um <- paste0(inin, "/", names(v05g), ".asc")
  
  for (j in 1:nlayers(v05g)) {
    writeRaster(v05g[[j]], filename = nms05um[j], format = "ascii")
  }
}

### 10 min
for (i in 1:length(backs1)) {
  sel <- read.csv(paste0(oe1[i], "/selected_models.csv"), stringsAsFactors = F)
  sets <- strsplit(sel[, 1], split = "_")[[1]]
  
  set <- paste(sets[5:length(sets)], collapse = "_")
  bset <- read.csv(paste0(backs1[i], "/", set, ".csv")) 
  gvs <- colnames(bset)[-(1:3)]; print(gvs); cat("\n")
  
  dir.create(paste0("Var_selection10/G_variables/", set))
  inin <- paste0("Var_selection10/G_variables/", set, "/world")
  dir.create(inin)
  
  v10g <- v10_bcei[[gvs]] 
  nms10um <- paste0(inin, "/", names(v10g), ".asc")
  
  for (j in 1:nlayers(v10g)) {
    writeRaster(v10g[[j]], filename = nms10um[j], format = "ascii")
  }
}


## model 0.5 degrees
btf <- paste0("Var_selection05/batchf_", all_mnames)
fmod <- paste0("Var_selection05/F_model_", all_mnames)
gvar <- "Var_selection05/G_variables"

for (i in 1:length(backs)) {
  kuenm_mod_swd(occ.joint = oj, back.dir = backs[i], out.eval = oe[i], 
                batch = btf[i], rep.n = 1, rep.type = "Crossvalidate", 
                jackknife = F, out.format = "cloglog", project = T, 
                G.var.dir = gvar, ext.type = "ext", maxent.path = mx, 
                out.dir = fmod[i], wait = T)
}

## model 10 min 
btf1 <- paste0("Var_selection10/batchf_", all_mnames)
fmod1 <- paste0("Var_selection10/F_model_", all_mnames)
gvar1 <- "Var_selection10/G_variables"

for (i in 1:length(backs)) {
  kuenm_mod_swd(occ.joint = oj1, back.dir = backs1[i], out.eval = oe1[i], 
                batch = btf1[i], rep.n = 1, rep.type = "Crossvalidate", 
                jackknife = F, out.format = "cloglog", project = T, 
                G.var.dir = gvar1, ext.type = "ext", maxent.path = mx, 
                out.dir = fmod1[i], wait = T)
}
# ------------------------------------------------------------------------------


# Exploring ENMs in E space to understand variable roles -----------------------
## variables of interest
vin05 <- c("BIO_6", "BIO_12", "BIO_14", "RSR", "ASRQH", all_mnames)
vin10 <- c("BIO_5", "BIO_14", "ASRQH", "CEC", all_mnames)

## variables and predictions to data.frames (GLM)
glm_epred05 <- data.frame(na.omit(v05_bcei[]), na.omit(stack(glm_gpred05)[]))
glm_epred10 <- data.frame(na.omit(v10_bcei[]), na.omit(stack(glm_gpred10)[]))

## maxent predictions
max05 <- list.files(path = "Var_selection05", 
                    pattern = "Spirodela_polyrhiza_world.asc$", full.names = TRUE,
                    recursive = TRUE)
max05 <- stack(max05)
names(max05) <- all_mnames

max10 <- list.files(path = "Var_selection10", 
                    pattern = "Spirodela_polyrhiza_world.asc$", full.names = TRUE,
                    recursive = TRUE)
max10 <- stack(max10)
names(max10) <- all_mnames

## variables and predictions to data.frames (Maxent)
max_epred05 <- data.frame(na.omit(v05_bcei[]), na.omit(max05[]))
max_epred10 <- data.frame(na.omit(v10_bcei[]), na.omit(max10[]))


## keeping only relevant variables
glm_epred05 <- glm_epred05[, vin05]
glm_epred10 <- glm_epred10[, vin10]

max_epred05 <- max_epred05[, vin05]
max_epred10 <- max_epred10[, vin10]

## saving values 
save(glm_epred05, glm_epred10, max_epred05, max_epred10, 
     file = "Data/variable_values.RData")

## plotting 3d representations of enms in E space
## the code presented below allows to replicate the figures in the paper 
## as each figure can be rotated to find appropriate visualizations we present 
## the code for one figure (figures for distinct calibration areas can be 
## generated changing the value in the object cal_area)
## the part of the code used to save figures is commented

### Maxent results
cal_area <- "intersection"

### this is for results at 30 min
colf <- as.factor(max_epred05[, cal_area])
colpal <- magma(length(levels(colf)))

varcmbs <- combn(colnames(max_epred05)[1:5], 3)

open3d(windowRect = c(5, 5, 895, 1195))
layout3d(matrix(1:12, 4, byrow = TRUE))

pls <- lapply(1:ncol(varcmbs), function(x) {
  vs <- varcmbs[, x]
  vnam <- gsub("_", " ", vs)
  
  plot3d(max_epred05[, vs], col = alpha(colpal[colf], 0.5), xlab = "", ylab = "",
         zlab = "", pch = 16, size = 3, axes = FALSE)
  par3d(cex = 0.5); rgl::axes3d()
  rgl::title3d(xlab = vnam[1], ylab = vnam[2], zlab = vnam[3])
})

next3d()
bgplot3d({
  par(mar = rep(0, 4)); plot.new()
  legend_bar("left", col = magma(255), width_prop = 0.05, heigh_prop = 0.4, 
             title = "Suitability")
})

#rgl.snapshot(paste0("Figures/Sup_max05_", cal_area,"_e_space.png"))


### this is for results at 10 min
colf <- as.factor(max_epred10[, cal_area])
colpal <- magma(length(levels(colf)))

varcmbs <- combn(colnames(max_epred10)[1:4], 3)

open3d(windowRect = c(5, 5, 600, 600))
layout3d(matrix(c(1:5, 5), 2), widths = c(10, 10, 2))

pls <- lapply(1:ncol(varcmbs), function(x) {
  vs <- varcmbs[, x]
  vnam <- gsub("_", " ", vs)
  
  plot3d(max_epred10[, vs], col = alpha(colpal[colf], 0.5), xlab = "", ylab = "",
         zlab = "", pch = 16, size = 3, axes = FALSE)
  par3d(cex = 0.5); rgl::axes3d()
  rgl::title3d(xlab = vnam[1], ylab = vnam[2], zlab = vnam[3])
})

next3d()
bgplot3d({
  par(mar = rep(0, 4)); plot.new()
  legend_bar("left", col = magma(255), width_prop = 0.35, heigh_prop = 0.4, 
             title = "Suitability")
})

#rgl.snapshot(paste0("Figures/Sup_max10_", cal_area,"_e_space.png"))



### GLM results
cal_area <- "intersection"

### this is for results at 30 min
colf <- as.factor(glm_epred05[, cal_area])
colpal <- magma(length(levels(colf)))

varcmbs <- combn(colnames(glm_epred05)[1:5], 3)

open3d(windowRect = c(5, 5, 895, 1195))
layout3d(matrix(1:12, 4, byrow = TRUE))

pls <- lapply(1:ncol(varcmbs), function(x) {
  vs <- varcmbs[, x]
  vnam <- gsub("_", " ", vs)
  
  plot3d(glm_epred05[, vs], col = alpha(colpal[colf], 0.5), xlab = "", ylab = "",
         zlab = "", pch = 16, size = 3, axes = FALSE)
  par3d(cex = 0.5); rgl::axes3d()
  rgl::title3d(xlab = vnam[1], ylab = vnam[2], zlab = vnam[3])
})

next3d()
bgplot3d({
  par(mar = rep(0, 4)); plot.new()
  legend_bar("left", col = magma(255), width_prop = 0.35, heigh_prop = 0.4, 
             title = "Suitability")
})

#rgl.snapshot(paste0("Figures/Sup_glm05_", cal_area,"_e_space.png"))

### this is for results at 10 min
colf <- as.factor(glm_epred10[, cal_area])
colpal <- magma(length(levels(colf)))

varcmbs <- combn(colnames(glm_epred10)[1:4], 3)

open3d(windowRect = c(5, 5, 600, 600))
layout3d(matrix(c(1:5, 5), 2), widths = c(10, 10, 2))

pls <- lapply(1:ncol(varcmbs), function(x) {
  vs <- varcmbs[, x]
  vnam <- gsub("_", " ", vs)
  
  plot3d(glm_epred10[, vs], col = alpha(colpal[colf], 0.5), xlab = "", ylab = "",
         zlab = "", pch = 16, size = 3, axes = FALSE)
  par3d(cex = 0.5); rgl::axes3d()
  rgl::title3d(xlab = vnam[1], ylab = vnam[2], zlab = vnam[3])
})

next3d()
bgplot3d({
  par(mar = rep(0, 4)); plot.new()
  legend_bar("left", col = magma(255), width_prop = 0.35, heigh_prop = 0.4, 
             title = "Suitability")
})

#rgl.snapshot(paste0("Figures/Sup_glm10_", cal_area,"_e_space.png"))

# ------------------------------------------------------------------------------
