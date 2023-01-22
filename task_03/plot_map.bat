@echo off
gmt gmtset FORMAT_GEO_MAP ddd

REM --- MODIFY VARIABLES BELOW FOR YOUR REGION/FILES!! ----
REM --- furthermore, you have to adjust the spacing of the grid lines
REM --- and label spacing of the colorbar	

REM --- filename of ascii file (without ending) that you want to plot
set file=DTU15MDT_2min.mdt
REM --- grid resolution in geographical degrees
set grid_res=30m
REM --- map projection
set map_proj=B-130/65/45/65/18c
REM --- region
set region=-145/-110/45/65
REM --- color palette
set color_palette=haxby
REM --- map title
set title=British Columbia, Canada

REM =============================================================
REM --- other variables ----
set grid_file=./data/%file%.nc

REM --- create map ---
echo Map is created:
gmt begin %file% png

REM --- create color palette table ---
gmt grdsample %grid_file% -G./data/DTU15MDT_30min.mdt.nc -I%grid_res% -R%region%
gmt makecpt -C%color_palette% -T-50/50/0.001 -Z
gmt grd2cpt %grid_file% -C%color_palette% -Z 

REM --- plot grid ---
gmt grdimage -J%map_proj% -R%region% %grid_file% -Q
gmt coast -Bxa5g5 -Bya5g5 -BWESN+t"%title%" -W0.25p,80/80/80 -Df -N1/1.25p,black -V 

REM --- plot legend including color scale ---
echo Legend is created...
gmt colorbar -Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V
gmt text -F+cBL+t"DTU15 MDT - mean dynamic topography" -N -D6.65c/-1c
gmt text -F+cBL+t"Editors: Christopher Mahn, Silas Teske, Joshua Wolf" -N -D5.15c/-1.5c

gmt end show
@echo on