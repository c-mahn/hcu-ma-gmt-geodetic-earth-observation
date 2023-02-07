@echo off
gmt gmtset FORMAT_GEO_MAP ddd

REM --- MODIFY VARIABLES BELOW FOR YOUR REGION/FILES!! ----
REM --- furthermore, you have to adjust the spacing of the grid lines
REM --- and label spacing of the colorbar
REM --- filename of nc file (without ending) that you want to plot
set file=test_vel
REM --- grid resolution in geographical degrees
set grid_res=0.5
REM --- map projection
set map_proj=M300/18c
REM --- region
set region=260/360/12/75
REM --- color palette
set color_palette=haxby
REM --- map title
set title=Gulf stream geostrophic velocity

REM =============================================================

gmt begin test_vel png

REM --- create map ---
echo Map is created:

REM --- create color palette table ---
gmt makecpt -C%color_palette%  -T0/0.8/0.001 -Z

REM --- plot grid ---
gmt grdimage -J%map_proj% -R%region% %file%.nc -Q
gmt coast -Bxa10g10 -Bya20g20 -G200/200/200 -BWESN+t"%title%" -W1p,80/80/80 -Di -V 

REM --- plot color scale ---
gmt colorbar -Dx0c/-2c+w17c/0.35c+h -B0.2+l"Velocity [m/s]" -V 

gmt end show
@echo on