# ppm_example
scripts for running point process models for many species

# Biodiversity data
Citation: GBIF.org (21 May 2018) GBIF Occurrence Download https://doi.org/10.15468/dl.nt0ftl
Downloaded from GBIF as a .csv and filtered as per requirements. We've cleaned data for incomplete sceintific names, goegraphical infomration, years, etc. Finally, we've retained species with >20 records.We have the scripts for data cleaning for global datasets. 

# Covariate data [under dev for global datasets]
We have scirpts for direct downloads and processing of global covariate data


# ------------------------
# TO BE RESOLVED:
1. Use tile 29 data for model building and vn data for prediction...
Probably not when the analysis become global because in that case we will have to use country + tile data for every region..
Maybe an option for species with less records when fittign models to regional subsets?


2 Quadrature (background) points
LATER: For global analyses, generate points based on species-regions as per convex hull+buffer aroubnd occurrence points.


3. Estimate species weights
LATER: not quite right yet


4.Move species occurrence points falling off the mask to nearest 'land' cells
We lose data again i.e. number of unique locations is reduced. This can be problematic for ppms...
  # nrow(unique(outside_pts))
  # nrow(unique(land))
  
  
5. Extract covariates for presence points
For landuse: Take the raster value with lowest distance to point AND non-NA value in the raster
This takes very logn to run for global layers; find alternative...


6. Catch errors in models
Think of how to reduce these. Possible quasiseperation issues due to spatially restricted...could resolved when biuilding models on reduced area?
At the moment, the script finds species with errors and reruns the models. If errors persist, species with errors are discarded.
