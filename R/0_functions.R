
## Author: Payal Bal
align.maskNA <- function(raster_stack, region_mask) {
  # update mask based on NAs in covariate stack
  
  print(paste0("# NAs in input mask: ", summary(region_mask)[6]))
  
  for (i in names(raster_stack)){
    if (!sum(is.na(region_mask@data@values)) == sum(is.na(raster_stack[[i]]@data@values))) {
      nona <- which(!is.na(values(region_mask))) # non-na values in mas
      nas <- which(is.na(values(raster_stack[[i]])[nona])) # na values in covariate
      xys <- xyFromCell(region_mask, nona)
      xys <- xys[nas,]
      values(region_mask)[cellFromXY(region_mask, xys)] <- NA
    }
  }
  
  new_mask <- region_mask
  
  print(paste0("# NAs in output mask: ", summary(new_mask)[6]))
  return(new_mask)
}


## Author: Nick Goulding
nearestLand <- function (points, raster, max_distance) {
  # get nearest non_na cells (within a maximum distance) to a set of points
  # points can be anything extract accepts as the y argument
  # max_distance is in the map units if raster is projected
  # or metres otherwise
  
  # function to find nearest of a set of neighbours or return NA
  nearest <- function (lis, raster) {
    neighbours <- matrix(lis[[1]], ncol = 2)
    point <- lis[[2]]
    # neighbours is a two column matrix giving cell numbers and values
    land <- !is.na(neighbours[, 2])
    if (!any(land)) {
      # if there is no land, give up and return NA
      return (c(NA, NA))
    } else{
      # otherwise get the land cell coordinates
      coords <- xyFromCell(raster, neighbours[land, 1])
      
      if (nrow(coords) == 1) {
        # if there's only one, return it
        return (coords[1, ])
      }
      
      # otherwise calculate distances
      dists <- sqrt((coords[, 1] - point[1]) ^ 2 +
                      (coords[, 2] - point[2]) ^ 2)
      
      # and return the coordinates of the closest
      return (coords[which.min(dists), ])
    }
  }
  
  # extract cell values within max_distance of the points
  neighbour_list <- extract(raster, points,
                            buffer = max_distance,
                            cellnumbers = TRUE)
  
  # add the original point in there too
  neighbour_list <- lapply(1:nrow(points),
                           function(i) {
                             list(neighbours = neighbour_list[[i]],
                                  point = as.numeric(points[i, ]))
                           })
  
  return (t(sapply(neighbour_list, nearest, raster)))
}


## Author: Payal Bal
catch_errors <- function(i, ppm_models, species_names, errorfile) {
  # catch errors in model outputs
  # -> list of outputs for models without errors, 
  # -> list of outputs for models with errors, 
  # -> text file with list of models with errors (species indices and species names)
  
  cat('Checking model for ', species_names[i],'for errors\n')
  for(i in 1:length(ppm_models)){
    if(!class(ppm_models[[i]])[1] == "try-error") {
      model_list[[n]] <- ppm_models[[i]]
      n <- n+1
    }else{
      print(paste0("Model ",i, " for '", spp[i], "' has errors"))
      cat(paste(i, ",", spp[i], "\n"),
          file = errorfile, append = T)
      error_list[[m]] <- ppm_models[[i]]
      m <- m+1
    }
  }
}

