SELECT ST_Compiled.Transect_ID AS TransectID, First(ST_Compiled.Site_ID) AS SiteID, First(CStr(ST_Compiled.Date)) AS [Date], Val(Left([Transect_ID],(InStr(1,[Transect_ID],'.')-1))) AS TrSiteID, First(ST_Compiled.T_Area) AS TArea, First(ST_Compiled.Number_Sightings) AS Number_Sightings, IIf(Max(ST_Compiled.Number_Observers) Is Null,1, Max(ST_Compiled.Number_Observers)) AS Number_Observers, Max(ST_Compiled.Weather) AS Weather, Max(ST_Compiled.Cloud_Cover) AS Cloud_Cover, Max(ST_Compiled.Wind) AS Wind, Max(ST_Compiled.Canopy_Cover) AS Canopy_Cover, Max(ST_Compiled.Subcanopy_Cover) AS Subcanopy_Cover, First(ST_Compiled.Start_Time) AS Start_Time, First(ST_Compiled.T_Area) AS Tlength
FROM ST_Compiled
WHERE (((ST_Compiled.Method)='ST') AND Fatal_problem='FALSE')
GROUP BY ST_Compiled.Transect_ID;
