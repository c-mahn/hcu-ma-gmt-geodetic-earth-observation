@echo off
REM gmt begin geoid_height png
REM gmt begin grav_anom_surface png
gmt begin grav_anom_satellite png

echo Map is created ...
REM gmt basemap -JN12c -R-180/180/-90/90 -B+t"Geoid Heights" -B
REM gmt basemap -JN12c -R-180/180/-90/90 -B+t"Gravity Anomalies Surface" -B
gmt basemap -JN12c -R-180/180/-90/90 -B+t"Gravity Anomalies Satellite" -B

REM gmt grd2cpt geoid_height.nc -Ccyclic.cpt -Sh -E100 -Z -V
REM gmt grdimage geoid_height.nc -I+a315+ne0.2

REM gmt grd2cpt grav_anom_surface.nc -Ccyclic.cpt -Sh -E100 -Z -V
REM gmt grdimage grav_anom_surface.nc -I+a315+ne0.2

gmt grd2cpt grav_anom_satellite.nc -Ccyclic.cpt -Sh -E100 -Z -V
gmt grdimage grav_anom_satellite.nc -I+a315+ne0.2
	
gmt coast -W.1,120/120/120 -B

echo Legend is created...
gmt text -F+cBL+t"Editors: Christopher Mahn, Silas Teske, Joshua Wolf" -N -D2.3c/-1c
REM gmt colorbar -Dx0c/-1.5c+w12c/0.25c+h -B20+l"[m]" -V
gmt colorbar -Dx0c/-1.5c+w12c/0.25c+h -B0.05+l"[mm]" -V

gmt end show
@echo on
