@echo off
gmt gmtset FORMAT_GEO_MAP ddd

REM --- MODIFY VARIABLES BELOW FOR YOUR REGION/FILES!! ----
REM --- furthermore, you have to adjust the spacing of the grid lines
REM --- and label spacing of the colorbar
REM --- filename of nc file (without ending) that you want to plot
set file=test_mdt
REM --- filename of u and v velocities (without ending) for plotting the vectors
set file_u=test_u
set file_v=test_v
REM --- grid resolution in geographical degrees
set grid_res=0.5
REM --- map projection
set map_proj=M300/18c
REM --- region
set region=260/360/12/75
REM --- color palette
set color_palette=haxby
REM --- map title
set title=Gulf stream MDT

REM =============================================================

gmt begin test_mdt png

REM --- create map ---
echo Map is created:

REM --- create color palette table ---
gmt makecpt -C%color_palette% -T-1/1/0.001 -Z

REM --- plot grid ---
gmt grdimage -J%map_proj% -R%region% %file%.nc -Q
gmt grdvector %file_u%.nc %file_v%.nc -W0/0/0 -I2 -S2 -Q0.2+e
gmt coast -Bxa10g10 -Bya20g20 -G200/200/200 -BWESN+t"%title%" -W1p,80/80/80 -Di -V 

REM --- plot color scale ---
gmt colorbar -Dx0c/-2c+w17c/0.35c+h -B0.5+l"MDT [m]" -V 

gmt end show
@echo on