SELECT SOL_Compiled.Transect_ID AS TransectID, SOL_Compiled.Sighting_ID AS SightingID, SOL_Compiled.Perp_Dist AS Perp_Dist, CStr(Date)
FROM SOL_Compiled
WHERE (VAL(SOL_Compiled.Sighting_Number)>0) AND (SOL_Compiled.Method='SOL' OR SOL_Compiled.Method='DOL');
