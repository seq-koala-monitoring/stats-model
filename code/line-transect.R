## Line Transect functions

library(sf)

# Create spatial object from eastings and northings information
fcn_line_transect_sf <- function(lineTransect) {
  lineTransectSpatial <- purrr::pmap(
    list(Start_Eastings = lineTransect$Start_Eastings, 
         Start_Northings = lineTransect$Start_Northings, 
         End_Eastings= lineTransect$End_Eastings, 
         End_Northings =lineTransect$End_Northings
    ), function(Start_Eastings, Start_Northings, End_Eastings, End_Northings) {
      st_linestring(matrix(c(Start_Eastings, Start_Northings, End_Eastings, End_Northings), 2, 2, byrow = T))
    })
  
  lineTransectSf <- st_sfc(lineTransectSpatial)
  st_crs(lineTransectSf) <- 7856 ## Set to GDA2020
  df <- st_sf(cbind(lineTransect,lineTransectSf))
  return(df)
}


