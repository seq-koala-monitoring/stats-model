SELECT UAoA_Compiled.Area_ID AS TransectID, First(UAoA_Compiled.Site_ID) AS SiteID, First(CStr(UAoA_Compiled.Date)) AS [Date], Val(Left(UAoA_Compiled.Area_ID,(InStr(1,UAoA_Compiled.Area_ID,'.')-1))) AS TrSiteID, First(UAoA_Compiled.Site_Area) AS TArea, First(UAoA_Compiled.Number_Sightings) AS Number_Sightings, Max(UAoA_Compiled.Number_Observers) AS Number_Observers, Max(UAoA_Compiled.Weather) AS Weather, Max(UAoA_Compiled.Cloud_Cover) AS Cloud_Cover, Max(UAoA_Compiled.Wind) AS Wind, Max(UAoA_Compiled.Canopy_Cover) AS Canopy_Cover, Max(UAoA_Compiled.Subcanopy_Cover) AS Subcanopy_Cover, First(UAoA_Compiled.Start_Time) AS Start_Time
FROM UAoA_Compiled
WHERE UAoA_Compiled.Method='UAoA'
GROUP BY UAoA_Compiled.Area_ID;

