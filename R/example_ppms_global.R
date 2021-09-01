## Example global PPMs

## Set up work environment ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")

x <- c('data.table', 'sp', 'raster', 
       'sense', 'tools', 'bitops', 'RCurl', 
       'rgdal', 'gdalUtils', 'usethis', 'rgeos',
       'spatstat.data', 'nlme', 'rpart', 'spatstat', 'ppmlasso')
lapply(x, require, character.only = TRUE)
rm(x)

## Folder paths
data_dir <- "./data"
output_dir <- file.path(data_dir, "outputs")
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

## Functions
source("./gdal_calc.R")



## Prepare data ####
## You can ignore this part as it is mainly data processing and skip to line 121
## This is probably what takes most time for global analyses :) 
## so I've left some useful gdal processingn scripts in 
## (we can chat about this when relevant)

## >> Mask ~100km2 WGS84 ####
## step one_create mask from WorldClim layer
mask_file <- file.path(output_dir, "bio1_10km.bil")
global.mask <- raster(mask_file) 
global.mask[which(!is.na(global.mask[]))] <- 1
unique(values(global.mask))
wgs_crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
crs(global.mask) <- wgs_crs
writeRaster(global.mask, filename = file.path(output_dir, "globalmask_wgs.tif"), 
            format = "GTiff")

rm(mask_file, wgs_crs, global.mask)


## >> Covariate data - processing ####
bio1 <- file.path(output_dir, "bio1_10km.bil")
srtm <- file.path(data_dir, "srtm/mn30_grd/srtm.adf")
bio1_rcp85 <- file.path(data_dir, "processed/bc85bi701_treated.tif")

m <- raster(file.path(output_dir, "globalmask_wgs.tif"))
new_res <- res(m)
new_crs <- crs(m)
new_extent <- extent(m)

infile <- bio1_rcp85
outfile = file.path(output_dir, paste0(tools::file_path_sans_ext(basename(infile)), "_wgs.", tools::file_ext(infile)))
outfile = gsub("\\..*", ".tif", outfile)

system(paste0("gdalwarp -overwrite -ot Byte -tr ", 
              paste(new_res, collapse = " "), " -te ", 
              paste(new_extent[1], new_extent[3], 
                    new_extent[2], new_extent[4]), 
              " -t_srs '", new_crs, "' ",
              infile, " ", outfile))

infile <- outfile
outfile <- sub("_wgs", "_wgs_mask", outfile)
mask_file <- file.path(output_dir, "globalmask_wgs.tif")
system(paste0("gdal_calc.py -A ", infile, " -B ", mask_file,
              " --calc='((B==1)*A)+(-9999*(B!=1))' --NoDataValue=-9999",
              " --outfile=", outfile))

file.remove(infile)


## >> Find min non-NA set values across mask and covariates and sync NAs ####
## Covariate files
bio1 = file.path(output_dir, "bio1_10km_wgs_mask.tif")
srtm = file.path(output_dir, "srtm_wgs_mask.tif")
bio1_rcp85 = file.path(output_dir, "bc85bi701_treated_wgs_mask.tif")
infile = c(bio1, srtm, bio1_rcp85)
file.exists(infile)

## Create copy of mask file
mask_file2 <- file.path(output_dir, "globalmask_wgs_minNA.tif")
file.copy(mask_file, mask_file2)

## Update mask based on covariate layers for syncing NAs
for(j in 1:length(infile)){
  
  print(paste0("Updating mask using file = ", j, " of ", length(infile) ,":", infile[j]))
  
  # gdalcalc(calc="((A>-9999) & (B==1))", infile = list(A=infile[j], B=mask_file2),outfile = mask_file2,
  #          NoDataValue=-9999, overwrite=TRUE)
  
  system(paste0("gdal_calc.py -A ", infile[j], " -B ", mask_file2,
                " --calc='((A>-9999) & (B==1))' --NoDataValue=-9999",
                " --outfile=", mask_file2))
}

## Update covariate layers based onn updated mask for syncing NAs
for(j in 1:length(infile)){
  
  print(paste0("Masking file = ", j, " of ", length(infile) ,":", infile[j], " using updated mask"))
  
  gdalmask(infile = infile[j], mask = mask_file2, outfile = infile[j], output_Raster = FALSE, overwrite=TRUE, verbose=TRUE)
  
}

freq(raster(mask_file))
freq(raster(mask_file2))
raster(mask_file2)
unique(values(raster(mask_file2)))
rm(mask_file, mask_file2)



## PPMs ####
## >> Biodiversity data ####
gbif <- fread(file.path(data_dir, "2019-05-14_gbif_iucnsp.csv"))
gbif[, .N, species]
gbif <- gbif[,c(4,3,2,5,1)]

## >> Covariate data for modelling
bio1 = file.path(output_dir, "bio1_10km_wgs_mask.tif")
srtm = file.path(output_dir, "srtm_wgs_mask.tif")
bio1_rcp85 = file.path(output_dir, "bc85bi701_treated_wgs_mask.tif")
infile = c(bio1, srtm, bio1_rcp85)
file.exists(infile)
mask_file <- file.path(output_dir, "globalmask_wgs_minNA.tif")


## >> Quadrature (background) points ####
cov.mod <- stack(raster(infile[1]), raster(infile[2]))
names(cov.mod) <- c("bio1", "srtm")

global.mask <- raster(mask_file)
global.mask[which(is.na(global.mask[]))] <- 0
rpts <- rasterToPoints(global.mask, spatial=TRUE)
backxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
backxyz <- backxyz[,-1]
backxyz <- na.omit(cbind(backxyz, as.matrix(cov.mod)))
backxyz200k <- backxyz[sample(nrow(backxyz), 200000), ]
backxyz200k$Pres <- rep(0, dim(backxyz200k)[1])

## Checks
summary(backxyz200k)
plot(global.mask, legend = FALSE)
plot(rasterFromXYZ(backxyz200k[,1:3]), col = "black", add = TRUE, legend=FALSE)


## >> Prediction points ####
cov.pred <- stack(raster(infile[3]), raster(infile[2]))
names(cov.pred) <- c("bio1_2070_rcp85", "srtm")
predxyz <- data.frame(rpts@data, X=coordinates(rpts)[,1], Y=coordinates(rpts)[,2])     
predxyz <- predxyz[,-1]
predxyz <- na.omit(cbind(predxyz, as.matrix(cov.pred)))


## >> Define model parameters ####
## Define global ppm variables
ppm_terms <- names(backxyz200k) [1:(length(names(backxyz200k))-1)]

## Specify covariates with interactions
interaction_terms <- c("X", "Y")

## Fix n.fits = 20 and max.lambda = 100
lambdaseq <- round(sort(exp(seq(-10, log(100 + 1e-05), length.out = 20)), 
                        decreasing = TRUE),5)

## Estimate weights for background points (only works in WGS)
## calculate the area of the globe and then work out the weights based on the total area divided by number of points
global.mask <- raster(mask_file)
ar <- raster::area(global.mask)
ar <- mask(ar,global.mask)
totarea <- cellStats(ar,'sum')*1000 ## in meters
area_offset <- extract(ar, backxyz200k[,c('X','Y')], small = TRUE, fun = mean, na.rm = TRUE)*1000 ## in meters
bkgrd_wts <- c(totarea/area_offset)

  ## Estimate species weights - under dev, please ignore for now
  # spdat <- gbif
  # species_names <- levels(factor(spdat$species))
  # spwts <- list()
  # for(i in seq_along(species_names)){
  #       print(i)
  #       spxy <- spdat[spdat$species %in% species_names[i], c(4,3)]
  #       names(spxy) <- c("X", "Y")
  #       cellNo <- cellFromXY(ar,spxy)
  #       cellNoCounts <- table(cellNo)
  #       tmp_cell_area <- extract(ar, spxy, fun = mean, na.rm = TRUE)*1000
  #       tmp_dat <- data.frame(area=tmp_cell_area,cell_number=cellNo)
  #       tmp_wts <- ((tmp_dat$area*1000)/cellNoCounts[as.character(cellNo)])/totarea
  #       spwts[[i]] <- tmp_wts
  # }


## >> Fit PPMs ####
## this can be written as a function to apply to many species at once in parallel, 
## see: line 104 onwards https://github.com/payalbal/ppm_example/blob/master/R/2_run_ppms.R
## here I have set it out stepwise so you can see what is happening

spdat = as.data.frame(gbif) ## species data
species_names = levels(factor(spdat$species))
bkdat = backxyz200k ## background points and data
bkwts = bkgrd_wts ## background points weights
n.fits = 10
min.obs = 50

i = 1 

## Initialise log file
## Not necessary but this is handy if you're doing this for a lot of species
## If you don't want to use this, comment out lines 289-294
cat('Fitting a ppm to', species_names[i],'\nThis is the', i,'^th model of',length(species_names),'\n')
logfile <- paste0("./ppm_log_",gsub(" ","_",species_names[i]),"_",gsub("-", "", Sys.Date()), ".txt")
writeLines(c(""), logfile)

## Definne species specific data & remove NA (if any)
spxy <- spdat[spdat$species %in% species_names[i], c(1,2)]
names(spxy) <- c("X", "Y")


## Check if points fall outside the mask & decide what to do with these
## you could discard, to choose to move these to nearest 'land' cells (with caution!)
plot(global.mask)
points(spxy)


## Extract covariates for presence points using buffer
## If using land use as a covariate, 
## see line 140 in https://github.com/payalbal/ppm_example/blob/master/R/2_run_ppms.R (under dev)
## we can discuss this if this becomes relevant for your work
spxyz <- extract(cov.mod, spxy, buffer = 1000000, small = TRUE,
                 fun = mean, na.rm = TRUE)
spxyz <- cbind(spxy, spxyz)
spxyz <- na.omit(spxyz)


## Specify model based on number of observatios
## If number of observations < 20, do not fit a model.
## If number of observations >= 20, fit a model accordig to the followig rule:
##  1. If number of observations <= min.obs (default as 50), 
##    use only the spatial covariates i.e. lat and long. This gives 5 
##    terms: X, Y, X^2, Y^2, X:Y (linear, quadratic with interaction)
##  2. .... 50 - 90
##  2. If number of observatons > min.obs, 
##    use the one-in-ten rule (https://en.wikipedia.org/wiki/One_in_ten_rule)
##    where, an additional covariate is added for every 10 additonal obeservations.
##    Because we fit poly with degree 2, adding a covariate will add two 
##    terms x and x^2 to the model. Therefore, to operationalize this rule, 
##    we add 1 raw covariate for every 20 additional observations. All 10
##    covariates (in addition to lat long) are fit onlt at > 230 observations.
##    There are 12 covariates in all includign lat and log. 

nk <- nrow(spxyz)
min.obs <- 50
if (!(nk < 20)) {
  
  ## Add presence column to species occurrence table
  spxyz$Pres <- rep(1, dim(spxyz)[1])
  
  ## Specify weights for PPM data
  ppmxyz <- rbind(spxyz, bkdat)
  ppmxyz$wt <- NA
  ppmxyz[ppmxyz$Pres == 0,]$wt <- bkwts
  ppmxyz[ppmxyz$Pres == 1,]$wt <- 1e-6 #spwts[[i]]
  
  ## Specify PPM formula based on number of observations
  if(nk <= min.obs){
    ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                              ", degree = 2, raw = FALSE)",collapse =""))
  } else {
    if(nk > min.obs & nk <= min.obs + 40){
      ppmform <- formula(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                ", degree = 2, raw = FALSE) + factor(landuse)",collapse ="")) 
    } else {
      extra_covar <- ceiling((nk - min.obs)/20)
      if(extra_covar > 10) extra_covar <- 10 ## Fit all 10 covariates
      ppmform <- formula(paste0(paste0(" ~ poly(", paste0(interaction_terms, collapse = ", "),
                                       ", degree = 2, raw = FALSE) + factor(landuse)"), 
                                paste0(" + poly(", ppm_terms[!(ppm_terms %in% c(interaction_terms, "landuse"))][1:extra_covar],", degree = 2, raw = FALSE)", collapse = ""), collapse =""))
    }
  }
  
  ## Fit ppm & save output
  cat(paste("\nFitting ppm model for",i , " @ ", Sys.time(), "\n"), 
      file = logfile, append = T)
  cat(paste("   # original records for",i , " = ", dim(spxy)[1], "\n"), 
      file = logfile, append = T)
  cat(paste("   # extracted records for",i , " = ", dim(spxyz)[1], "\n"), 
      file = logfile, append = T)
  
  mod <- try(ppmlasso(formula = ppmform, data = ppmxyz, n.fits = n.fits, 
                      criterion = "bic", standardise = FALSE), silent=TRUE)
  
  gc()
  return(mod)
} else {
  return(NULL) ## i.e. if number of observations < 20, return NULL 
}

## Check & save output
mod
saveRDS(predmu, file = paste0("./MOD_", gsub(" ", "_", tolower(species_names[i]))))



## >> Prediction ####
newdata <- predxyz
bkdat <- backxyz200k
predmu <- list()

## Current predictions
preddat <- newdata[,c(1:2)]
predmu$current <- predict.ppmlasso(mod, newdata = preddat)

## Future predictions
preddat <- newdata
names(preddat) <- names(bkdat)[1:(length(names(bkdat))-1)]
predmu$rcp26 <- predict.ppmlasso(mod, newdata = preddat)
# predmu <- 1-exp(-predmu) ## gives relative probabilities

## Check & save output
predmu
saveRDS(predmu, file = paste0("./PRED_", gsub(" ", "_", tolower(species_names[i]))))