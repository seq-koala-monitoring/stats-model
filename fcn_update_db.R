#' Update the integrated database which contains koala survey data from 1996 onwards
#'
#' This function updates the integrated database with the latest koala survey data, combining both historical and recent information. It also updates the shapefile to reflect the data from the most current database.
#'
#' @param db character. File name (with extension) of the koala survey database to be added to the integrated database. DEFAULT: KoalaSurveyData2020_cur.accdb
#' @param db.integrated character. File name (with extension) of the koala survey database to be updated. DEFAULT: SEQKoalaDatabase_Integrated.accdb
#' @param survey_sites character. File name (with extension) of the latest spatial layer with koala survey sites. DEFAULT: KoalaSurveySites_231108.shp
#' @param update.spatial logical. If TRUE, updates the shapefile to reflect the data from the database to be added to the integrated database. DEFAULT: TRUE
#' @param spatial_transects character. File name (with extension) of the latest spatial layer with koala transects. DEFAULT: transects_v240517.shp
#' @param line_transect_width numeric. Measure of estimated line transect width. DEFAULT: 57.38
#' @param site.info character. File name of the most recent CSV file containing habitat data for koala survey sites. DEFAULT: site_info.csv
#' @param monitoring.units character. File name (with extension) of the latest spatial layer with the boundaries of monitoring units (i.e., previous genetic populations). DEFAULT: Population_boundaries_v2.shp
#'
#' @return An updated Microsoft Access database that integrates historical and recent koala survey data. An updated shapefile representing historical and recent transects as polygons.
#' @export
#'
#' @examples
#' # Run with default values
#' fcn_update_db()
fcn_update_db <- function(db = "KoalaSurveyData2020_cur.accdb",
                          db.integrated = "Integrated_SEQKoalaDatabase.accdb",
                          survey_sites = "KoalaSurveySites_231108.shp",
                          update.spatial = T,
                          spatial_transects = "transects_v240517.shp",
                          line_transect_width = 57.38,
                          site.info = "site_info.csv",
                          monitoring.units = "Population_boundaries_v2.shp"){
  
  
    # load function to create line transects as line feature
    source("code/line-transect.r")
  
    cat(
      "------------------------------------------------------------","\n",
      "IT MAY TAKE AROUND 10 MINS TO RUN DEPENDING ON YOUR COMPUTER","\n",
      "------------------------------------------------------------", "\n\n")
  
    # look for the database (db) within the working directory
    path <- list.files(pattern = db, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
    if (length(path) > 1) {
      stop(sprintf("There are multiple files named %s. Please, make sure to keep only the most up-to-date koala survey database in this working directory. Databases found in: %s", db))} # close {} of error message
    
    # create the db path for the odbc connection
    db_path <- paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};",
                      "DBQ=",
                      getwd(), "/", path)

    # establish odbc connection
    channel <- RODBC::odbcDriverConnect(db_path)
    
    # retrieve survey data table only
    dat <- RODBC::sqlFetch(channel, "tblKoalaSurveyData2020_cur")
    
    # close odbc connection
    RODBC::odbcClose(channel)
    
    # include number of sightings
    
    
    # adding the 2018 Daisy Hill Line Transects
    # create a temporary dataframe for merging
    tmp <- dat |> 
      dplyr::mutate(XID = as.character(ID),
             Start_Time = format(Start_Time, "%H:%M"),
             End_Time = format(End_Time, "%H:%M"),
             Time = format(Time, "%H:%M"),
             Year = as.numeric(format(Date, "%Y")))
    
    # standardize column names
    SOL_2018 <- read.csv("input/line_transects_2018/V2_MASTER_Daisy Hill_Line_tidy.csv") |>  
      dplyr::mutate(Date = as.POSIXct(Date),
             Weather = as.character(Weather),
             DBH = as.character(DBH),
             Wind = as.character(Wind),
             Cloud_Cover = as.character(Cloud_Cover),
             Canopy_Cover = as.character(Canopy_Cover)) |> 
      dplyr::select(-Concern, -Fatal_problem)
    
    # standardize the 2018 line transects dataframe 
    dat2 <- dplyr::bind_rows(tmp, SOL_2018) |> 
      dplyr::relocate(Year, .before = Date) |> 
      dplyr::relocate(Number_Sightings, .before = Sighting_Number) |>
      dplyr::rename(NumSightings = Number_Sightings) |> 
      dplyr::mutate(ID = dplyr::row_number())
    
    # include number of sightings per transect
    trans1 <- dat2 |>
      dplyr::mutate(Sighting_Number = ifelse(is.na(Sighting_Number),
                                      0, Sighting_Number)) |> 
      dplyr::group_by(Transect_ID) |>
      dplyr::mutate(NumSightings =
               if (dplyr::n() > 1) {
                 NumSightings = dplyr::n()
               }
             else if (dplyr::n() == 1 & Sighting_Number != 0) {
               NumSightings = 1
             }
             else {
               NumSightings = 0
             }) |> 
      dplyr::ungroup() |> 
      dplyr::mutate(Number_Sightings = ifelse(is.na(NumSightings),
                                       NumSightings,
                                       NumSightings)) |> 
      dplyr::relocate(Number_Sightings, .before = Sighting_Number) |> 
      dplyr::select(-NumSightings)
    
    
    # add length and area if not recorded
    trans.sf <- trans1 |> 
      dplyr::filter(!Method %in% c("IS", "VT")) |> 
      dplyr::group_by(Transect_ID) |> 
      dplyr::reframe(
        Start_Eastings = unique(Start_Eastings), 
        Start_Northings = unique(Start_Northings),
        End_Eastings = unique(End_Eastings),
        End_Northings = unique(End_Northings),
        Year = unique(Year)
      ) |> 
      fcn_line_transect_sf() |> 
      dplyr::mutate(TLength2 = as.numeric(st_length(geometry))) |> 
      dplyr::mutate(index = dplyr::row_number()) |> 
      sf::st_buffer(line_transect_width/2, endCapStyle = "FLAT")
    
    trans2 <- trans1 |> 
      dplyr::left_join(sf::st_drop_geometry(trans.sf) |> dplyr::select(Transect_ID, TLength2),
                by = "Transect_ID") |> 
      dplyr::mutate(T_Length = ifelse(is.na(T_Length),
                               TLength2, T_Length)) |>
      dplyr::mutate(T_Area = ifelse(is.na(T_Area),
                             T_Length * line_transect_width/1e4, T_Area),
             T_Width = ifelse(is.na(T_Width),
                              line_transect_width, T_Width)) |> 
      dplyr::select(-TLength2)
    
    
    #-----------------------------------------------------------
    # standardise Site ID
    # look for a file with the spatial representation of survey sites within the working directory
    path <- list.files(pattern = survey_sites, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
    if (length(path) > 1) {
      stop(sprintf("There are multiple files named %s. Please, make sure to keep only the most up-to-date spatial representation of survey sites in this working directory.", survey_sites))} # close {} of error message
    
    surveys <- sf::st_read(path) 
    sf::st_crs(surveys) <- "EPSG:7856"
    
    siteIds <-
      suppressWarnings({
        lapply(split(trans.sf, trans.sf$Transect_ID), function(d) {
          cat(
            paste0(
              "Matching Site_ID (",
              round(d$index/max(trans.sf$index)*100, 0),
              "% processed)","\n", "Transect_ID: ", d$Transect_ID, "\n\n"))
          
          sub <- sf::st_intersection(surveys, d) |>
            dplyr::select(Survey_Yr, Site, NAME, Subsite)
          
          if (nrow(sub) == 0) {
            Site_ID <- NA
            Subsite <- NA
            Site_Name <- NA
            
          } else if (nrow(sub) == 1) {
            Site_ID <- unique(sub$Site)
            Subsite <- unique(sub$Subsite)
            Site_Name <- unique(sub$NAME)
            
          } else {
            if (length(unique(sub$Site)) > 1) {
              sub <- sub |> 
                dplyr::mutate(inters.area = as.numeric(sf::st_area(geometry)) /
                                             1e4) |>
                dplyr::select(Survey_Yr, Site, NAME, Subsite, inters.area)
              
              max <- max(sub$inters.area, na.rm = T)
              sub <- sub[sub$inters.area == max,]
              
              year <-
                unique(unique(d$Year))
              survey_yr <- sub$Survey_Yr
              
              if (any(year %in% survey_yr)) {
                filter <- survey_yr == year
              } else {
                filter <-
                  survey_yr == survey_yr[which.max(survey_yr[which(survey_yr <= year)])]
              }
            } else {
              year <-
                unique(unique(d$Year))
              survey_yr <- sort(unique(sub$Survey_Yr))
              
              if (any(year %in% survey_yr)) {
                filter <- survey_yr == year
              } else {
                filter <-
                  survey_yr == survey_yr[which.max(survey_yr[which(survey_yr <= year)])]
              }
            }
            
            
            if (purrr::is_empty(filter)) {
              diff <- abs(year - survey_yr)
              filter <-
                survey_yr == survey_yr[which.min(diff)]
            }
            
            if (length(filter) != nrow(sub)) {
              sub <- sf::st_intersection(surveys, d) |>
                dplyr::mutate(inters.area = as.numeric(sf::st_area(geometry)) /
                                             1e4) |>
                dplyr::select(Survey_Yr, Site, NAME, Subsite, inters.area)
              sub <- sub[-which.min(sub$inters.area),]
              
              survey_yr <- sort(unique(sub$Survey_Yr))
              
              if (any(year %in% survey_yr)) {
                filter <- survey_yr == year
              } else {
                filter <-
                  survey_yr == survey_yr[which.max(survey_yr[which(survey_yr <= year)])]
              }
            }
            
            sub <- sub |>
              dplyr::distinct(Survey_Yr, Site, .keep_all = T)
            
            Site_ID <- sub[filter, "Site"]$Site
            Subsite <- sub[filter, "Subsite"]$Subsite
            Site_Name <- sub[filter, "NAME"]$NAME
          }
          
          df <- data.frame(
            Transect_ID = unique(d$Transect_ID),
            Site_ID = Site_ID,
            Subsite = Subsite,
            Site_Name = Site_Name
          )
          
          return(df)
        }) |> dplyr::bind_rows()
      })
    
    # Check if there are any duplicated Transect_IDs resulting from different Site_ID values for the same transect
    siteIds <- siteIds |>
      dplyr::mutate(dups = duplicated(Transect_ID))
    
    # Identify transects without a Site_ID
    siteNAs <- siteIds |>
      dplyr::filter(is.na(Site_ID)) |> 
      dplyr::pull(unique(Transect_ID))
    
    siteNAs <- dat |> 
      dplyr::filter(!Method %in% c("IS", "VT") & Transect_ID %in% siteNAs)
    
    siteNAs.n <- nrow(siteNAs)
    
    if (siteNAs.n > 0) {
      assign("siteNAs", siteNAs, envir = .GlobalEnv)
      
      warning(sprintf("There are %s transects without a matching site_ID in %s. Please ensure you are using the most up-to-date spatial representation of survey sites, as outdated information could make it hard to get accurate data for the model. To see which transects are missing a matching site_ID, check the siteNAs dataframe in the top right panel (i.e., your global environment)", siteNAs.n, survey_sites))
      } # close if (siteNAs.n > 0)
    
    # Assign the spatially retrieved Site_ID value to each transect of the data frame
    trans3 <- trans2 |> 
      dplyr::ungroup() |> 
      dplyr::left_join(siteIds, by = "Transect_ID") |> 
      
      # Resolve problems with duplicated column names
      dplyr::mutate(Site_Name = ifelse(is.na(Site_Name.y), Site_Name.x, Site_Name.y),
             Site_ID = ifelse(is.na(Site_ID.y), Site_ID.x, Site_ID.y),
             SiteID_Master = Site_ID.x,
             Subsite = ifelse(is.na(Subsite.y), NA, Subsite.y),
             Subsite_Master = Subsite.x) |> 
      dplyr::select(-Site_Name.x, -Site_Name.y, -Site_ID.x, 
             -Site_ID.y, -Subsite.x, -Subsite.y) |> 
      dplyr::relocate(Site_ID, .after = LGA) |> 
      dplyr::relocate(Subsite, .after = Site_ID) |> 
      dplyr::relocate(Site_Name, .after = Subsite) |> 
      dplyr::relocate(Number_Observers, .before = Observer_Obstacles) |> 
      
      # Remove unwanted columns
      dplyr::select(-Number_Observers_Master, -dups)
    
    # Define the source of the 2018 transects"V2_MASTER_Daisy Hill_Line"
    trans3 <- trans3 |> 
      dplyr::mutate(Source = ifelse(Year == 2018, "V2_MASTER_Daisy Hill_Line", "KoalaSurveyData2020_cur"))
    
    # Manually assign Site_ID and Subsite of a problematic transect
    trans3[trans3$Transect_ID == "221.0_10_SOL.20220503",
           "Site_ID"] <- 221
    trans3[trans3$Transect_ID == "221.0_10_SOL.20220503",
           "Subsite"] <- 3
    
    
    #-----------------------------------------------------------
    # UPDATE SHAPEFILE
    # look for the most up-to-date spatial representation of transects within the working directory
    if(update.spatial){
      path <- list.files(pattern = spatial_transects, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
      if (length(path) > 1) {
        stop(sprintf("There are multiple files named %s. Please, make sure to keep only the most up-to-date spatial representation of transects in this working directory", spatial_transects))} # close {} of error message
      
      transects.sf <- sf::st_read(path) |> 
        dplyr::filter(!Source %in% gsub(".accdb", "", db)) 
      sf::st_crs(transects.sf) <- "EPSG:7856"

      lineTransect <- trans3 |> 
        
        # Remove incidental sightings and vehicle transects
        dplyr::filter(Method %in% c("DOL", "SOL") & Year > 2018) |> 
        
        # Reframe to retain only columns relevant for modelling
        dplyr::group_by(Transect_ID) |> 
        dplyr::reframe(
          Date = as.character(unique(Date)),
          Start_Eastings = unique(Start_Eastings), 
          Start_Northings = unique(Start_Northings),
          End_Eastings = unique(End_Eastings),
          End_Northings = unique(End_Northings),
          Method = "SOL", 
          Source = gsub(".accdb", "", db)
        ) 
      # Create a simple feature object to be exported as a shapefile
      suppressWarnings({
      lineTransectSf <- fcn_line_transect_sf(lineTransect) |> 
        sf::st_buffer(line_transect_width/2, endCapStyle = "FLAT") |>
        dplyr::mutate(Date = format(as.POSIXct(Date,"%Y-%m-%d"), "%d/%m/%Y")) |> 
        dplyr::select(TrnscID = Transect_ID, Date, Method, Source)
      })
      
      # Append to the source spatial representation
      trans.final <- dplyr::bind_rows(transects.sf, lineTransectSf) |> 
        dplyr::mutate(Year = as.integer(substring(Date, nchar(Date) - 4 + 1)))
      
      # Export as shapefile
      if(!file.exists("output/transects_spatial_representation")){
        dir.create("output/transects_spatial_representation")}
      suppressWarnings({
      sf::st_write(trans.final, 
                   paste0("output/transects_spatial_representation/transects_v",
                          format(Sys.Date(), "%y%m%d"),
                          ".shp"), append = F)
      })
      
      files.from <- list.files(pattern = paste0("transects_v", format(Sys.Date(), "%y%m%d")),
                               path = "output/transects_spatial_representation",
                               full.names = T)
      files.to <- sub("/([^/]*)\\.", "/Integrated_SEQKoalaDatabase_Spatial.", files.from)
      files.to <- sub("output", "input", files.to)
      
      file.copy(from = files.from, 
                to = files.to,
                overwrite = T)
      
    } # close conditional to update spatial representation of transects
  
    # Fill any missing value in number of observers (i.e., SOL = 1 and DOL = 2)
    trans3[is.na(trans3$Number_Observers) & trans3$Method == "DOL", "Number_Observers"] <- 2
    trans3[is.na(trans3$Number_Observers) & trans3$Method == "SOL", "Number_Observers"] <- 1
    
    #-----------------------------------------------------------
    # Survey site information and monitoring unit
      path <- list.files(pattern = site.info, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
      if (length(path) > 1) {
        stop(sprintf("There are multiple files named %s. Please, make sure to keep only the most up-to-date information about survey site habitat in this working directory", site.info))} # close {} of error message
      
      # Load file with habitat information
      info <- read.csv(path)
      
      trans4 <- trans3 |>
        dplyr::left_join(
          info |> dplyr::select(LGA, KPA, Site.Number,
                          Land.Type, Habitat.Type, Site.Name),
          dplyr::join_by(Site_ID == Site.Number)) |>
        dplyr:: mutate(
          LGA = ifelse(is.na(LGA.x), LGA.y, LGA.x),
          KPA = ifelse(is.na(KPA.x),LGA.y, KPA.x),
          Hab_Type = ifelse(is.na(Hab_Type), Habitat.Type, Hab_Type),
          Land_Type = ifelse(is.na(Land_Type), Land.Type, Land_Type),
          Site_Name_Master = ifelse(is.na(Site_Name_Master), Site.Name, Site_Name_Master),
          Site_Name = ifelse(is.na(Site_Name), Site.Name, Site_Name)) |>
        dplyr::select(-LGA.y,-LGA.x,-KPA.x, -KPA.y, -Land.Type, -Habitat.Type, -Site.Name)
      
      # Assign monitoring unit from shapefile
      # This code follows the same logic applied to extract Site ID
      path <- list.files(pattern = monitoring.units, recursive = T, ignore.case = T, full.names = F, include.dirs = T) 
      path <- path[(substring(path, nchar(path) - 3 + 1) == "shp")]
      
      if (length(path) > 1) {
        stop(sprintf("There are multiple files named %s. Please, make sure to keep only the most up-to-date information about survey site habitat in this working directory", monitoring.units))} # close {} of error message
      
      units <- sf::st_read(path)
      sf::st_crs(units) <- 7856
      
      tmp <- trans4 |> 
        # Remove incidental sightings and vehicle transects
        dplyr::filter(Method %in% c("DOL", "SOL")) |> 
        # Reframe to retain only columns relevant for modelling
        dplyr::group_by(Transect_ID) |> 
        dplyr::reframe(
          Date = as.character(unique(Date)),
          Start_Eastings = unique(Start_Eastings), 
          Start_Northings = unique(Start_Northings),
          End_Eastings = unique(End_Eastings),
          End_Northings = unique(End_Northings),
          Method = "SOL", 
          Source = gsub(".accdb", "", db)
        ) |> 
        fcn_line_transect_sf() |>  
        dplyr::mutate(index = dplyr::row_number()) |> 
        sf::st_buffer(line_transect_width/2, endCapStyle = "FLAT")
      
      M.unit <- suppressWarnings({
        lapply(split(tmp, tmp$Transect_ID), function(d){
          cat(
            paste0(
              "Matching monitoring unit (",
              round(d$index/max(tmp$index)*100, 0),
              "% processed)","\n", "Transect_ID: ", d$Transect_ID, "\n\n"))
          
          sub <- sf::st_intersection(units, d) |> 
            dplyr::select(GENPOP_NAM, Transect_ID)
          
          if(nrow(sub) == 0) {
            Monitoring_Unit <- NA
            Transect_ID <- d$Transect_ID
            
          } else if (nrow(sub) == 1) {
            Monitoring_Unit <- unique(sub$GENPOP_NAM)
            Transect_ID <- unique(sub$Transect_ID)
            
          } else if(length(unique(sub$GENPOP_NAM)) > 1){
            sub <- sf::st_intersection(units, d) |>
              dplyr::mutate(inters.area = as.numeric(st_area(geometry)) / 1e4) |>
              dplyr::select(GENPOP_NAM, Transect_ID, inters.area) 
            
            max <- max(sub$inters.area, na.rm = T)
            sub <- sub[sub$inters.area == max, ]
            
            Monitoring_Unit <- sub$GENPOP_NAM
            Transect_ID <- unique(sub$Transect_ID)
          }
          
          df <- data.frame(Transect_ID = Transect_ID,
                           Monitoring_Unit = Monitoring_Unit)
          
          return(df)
        })
      }) |> dplyr::bind_rows()
      
      # Check any missing data
      M.unit.NAs <- M.unit |> 
        dplyr::filter(is.na(Monitoring_Unit))
      if (nrow(M.unit.NAs) > 0) {
        assign("M.unit.NAs", M.unit.NAs, envir = .GlobalEnv)
        stop(sprintf("There are %s transects without a matching monitoring unit. Please, check transects in the dataframe named M.unit.NAs before continuing.", nrow(M.unit.NAs)))} # close {} of error message
      
      # Assign a corresponding monitoring unit to the dataframe
      trans5 <- trans4 |> 
        dplyr::left_join(M.unit, by = "Transect_ID") |> 
        dplyr::rename(Monitoring_Unit = Monitoring_Unit.y) |> 
        dplyr::select(-Monitoring_Unit.x)
      
    #-----------------------------------------------------------
    LTPerpDist <- trans5 |>
      dplyr::filter(Method %in% c("SOL", "DOL") & Number_Sightings > 0) |> 
      dplyr::mutate(Perp_Dist = ifelse(is.na(Perp_Dist), 
                                (Distance_Observer_Koala * 
                                   sin(Angle_Koala * pi/180) / 1), # convert from degrees to radians
                                Perp_Dist)) |> 
      dplyr::relocate(Sighting_ID, .before = Transect_ID)
    
    # Check for missing data
    LTPerpDist.na <- LTPerpDist |> 
      dplyr::filter(is.na(Perp_Dist)) |> 
      
      # Remove one potentially duplicate entry (ID = 1816)
      dplyr::filter(!is.na(Koala_Eastings)) |> 
      sf::st_as_sf(coords = c("Koala_Eastings", "Koala_Northings"),
               crs = 7856) # In this case is impossible to calculate the perpendicular distance because the Angle_Koala is also missing data
    
    
    # Calculate perpendicular distances using spatial representations of the transects and of koala sightings
    # Create an spatial object with transects
    LTPerpDist.sf <- LTPerpDist |> 
      dplyr::group_by(Transect_ID) |> 
      dplyr::reframe(
        Site_ID = unique(Transect_ID), 
        Date = unique(Date),
        T_Length = unique(T_Length),
        Number_Sightings = unique(Number_Sightings),
        Number_Observers = unique(Number_Observers),
        Start_Eastings = unique(Start_Eastings), 
        Start_Northings = unique(Start_Northings),
        End_Eastings = unique(End_Eastings),
        End_Northings = unique(End_Northings),
        Year = as.integer(unique(Year)),
        Sighting_ID = unique(Sighting_ID)
      ) |> 
      fcn_line_transect_sf() |> 
      dplyr::filter(Sighting_ID %in% LTPerpDist.na$Sighting_ID)
    
    
    # Calculate distance between the koala sighting and the transect line based on the spatial representation
    dist <- lapply(split(LTPerpDist.na, LTPerpDist.na$Sighting_ID), function(x){
      
      dist <- x |>
        sf::st_distance(
          dplyr::filter(LTPerpDist.sf,
                 LTPerpDist.sf$Sighting_ID == x$Sighting_ID)
        ) |> 
        as.numeric()
      
      z <- data.frame(Sighting_ID = x$Sighting_ID,
                      Perp_Dist = dist)
      
      return(z)
    }) |> 
      dplyr::bind_rows()
    
    
    # Update the data frame to fill the gap
    trans6 <- trans5 |> 
      dplyr::left_join(dist, by = "Sighting_ID") |> 
      dplyr::mutate(Perp_Dist = ifelse(is.na(Perp_Dist.x),
                                Perp_Dist.y,
                                Perp_Dist.x)) |> 
      dplyr::select(-Perp_Dist.x, -Perp_Dist.y)
    
    if(sum(duplicated(LTPerpDist$Sighting_ID)) > 0) {
      warning(sprintf("%i duplicate Sighting_ID", sum(duplicated(lineTransect$Sighting_ID))))
    } else {cat("There is a unique Sighting_ID per koala sighting")}
    
    # Final formatting
        # Create a template
    templ <- trans1[1, ] |> 
      dplyr::mutate(across(everything(), ~ NA)) |> 
      # Rename original Gen_Pop as Monitoring_Unit_Master
      dplyr::rename(Monitoring_Unit_Master = Gen_Pop) |> 
      dplyr::select(-ID, -XID) |> 
      # Change column type for combining
      dplyr::mutate(Date = as.character(Date))
    
    
    # Keep incidental sightings and vehicle transects separate
    IS_VT <- trans6 |> 
      dplyr::filter(!Method %in% c("SOL", "DOL")) |> 
      dplyr::mutate(Date = as.POSIXct(Date, format = "%d/%m/%Y"),
                    Tran_Plot = ifelse(is.na(Tran_Plot), Tran_plot, Tran_Plot)) |> 
      dplyr::rename("\"Date\"" = Date, "\"Time\"" = Time) 
    
    
    trans7 <- trans6 |>
      # Tidy up date format
      dplyr::mutate(Date = format(Date, format = "%d/%m/%Y"))  |> 
      # Assign empty values ("") instead of missing data (NAs)
      dplyr::mutate(
        dplyr::across(
          dplyr::everything(),
          ~ ifelse(is.na(.), "", .))) |> 
      # Rename original Gen_Pop as Monitoring_Unit_Master
      dplyr::rename(Monitoring_Unit_Master = Gen_Pop) |> 
      # Tidying 
      dplyr::relocate(Perp_Dist, .after = Angle_Koala) |> 
      dplyr::select(-Site_Name_Master, -ID, -Distance_bt_Observers) |> 
      dplyr::mutate(Site_Area = as.numeric(Site_Area),
             Subsite = as.numeric(Subsite))
    
    suppressWarnings({
    final <- dplyr::bind_rows(templ, trans7) |> 
      dplyr::mutate(Date = as.POSIXct(Date, format = "%d/%m/%Y"),
             Tran_Plot = ifelse(is.na(Tran_Plot), Tran_plot, Tran_Plot),
             Site_ID = as.integer(Site_ID),
             Total_Field_Staff = as.integer(Total_Field_Staff),
             Start_Eastings = as.integer(Start_Eastings),
             End_Eastings = as.integer(End_Eastings),
             Start_Northings = as.integer(Start_Northings),
             End_Northings = as.integer(End_Northings),
             Number_Observers = as.integer(Number_Observers),
             Count = as.integer(Count),
             Observer_Eastings = as.integer(Observer_Eastings),
             Observer_Northings = as.integer(Observer_Northings),
             Koala_Eastings = as.integer(Koala_Eastings),
             Koala_Northings = as.integer(Koala_Northings),
             Koala_Elevation = as.integer(Koala_Elevation),
             Koala_Bearing = as.integer(Koala_Bearing),
             Angle_Koala = as.integer(Angle_Koala),
             Field_Bearing = as.integer(Field_Bearing),
             T_Length = as.double(T_Length),
             T_Area = as.double(T_Area),
             Perp_Dist = as.double(Perp_Dist),
             Tree_Height = as.double(Tree_Height),
             Koala_Height = as.double(Koala_Height)) |> 
      dplyr::rename("\"Date\"" = Date, "\"Time\"" = Time) |> 
      dplyr::relocate(Monitoring_Unit_Master, .after = Study_Area) |> 
      dplyr::relocate(Monitoring_Unit, .before = LGA) |> 
      dplyr::select(-Tran_plot, -XID)
    final <- final[-1, ]
    })
    
    # update db
    path <- list.files(pattern = db.integrated, recursive = T, ignore.case = T, full.names = F, include.dirs = T)
    if (length(path) > 1) {
      stop(sprintf("There are multiple files named %s. Please, make sure to keep only the most up-to-date koala survey database in this working directory. Databases found in: %s", db.integrated))} # close {} of error message
    
    # create a copy the database to avoid modifying the original
    if(!file.exists("output/integrated_database")){
      dir.create("output/integrated_database")}
    
    file.copy(from = path, 
              to = "output/integrated_database/SEQKoalaDatabase_Integrated.accdb")
    
    # update path
    path.upd <- "output/integrated_database/SEQKoalaDatabase_Integrated.accdb"
    
    # create the db path for the odbc connection
    db_path <- paste0("Driver={Microsoft Access Driver (*.mdb, *.accdb)};",
                      "DBQ=",
                      getwd(), "/", path.upd)
    
    # establish odbc connection
    channel <- RODBC::odbcDriverConnect(db_path)
    
    # remove old SOL table survey data table only
    RODBC::sqlDrop(channel, "SOL_compiled")
    
    RODBC::sqlSave(
        channel,
        dat = final,
        tablename = "test",
        fast = F, 
        safer = F, 
        rownames = F, 
        colnames = F,
        varTypes = c("\"Date\"" = "date")
      )
    
        # close odbc connection
    RODBC::odbcClose(channel)
    
    file.copy(from = "output/integrated_database/SEQKoalaDatabase_Integrated.accdb", 
              to = "input/databases/Integrated_SEQKoalaDatabase.accdb")
    
} # finish function

    