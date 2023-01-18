@echo off
gmt gmtset FORMAT_GEO_MAP ddd

REM --- MODIFY VARIABLES BELOW FOR YOUR REGION/FILES!! ----
REM --- furthermore, you have to adjust the spacing of the grid lines
REM --- and label spacing of the colorbar

REM --- filename of ascii file (without ending) that you want to plot
set file=test_noGF_BC_01
REM --- grid resolution in geographical degrees
set grid_res=0.5
REM --- filename of polygon file that outlines your averaging region
set file_poly=../data/region_polygon.txt
REM --- map projection (second line is example for Greenland)
set map_proj=M-120/18c
REM set map_proj=B-45/75/60/75/18c
REM --- region
set region=-175/-110/41/62
REM --- color palette
set color_palette=haxby
REM --- map title
set title=British Columbia, Canada

REM =============================================================
REM --- other variables ----
set ascii_file=%file%.txt
set grid_file=%file%.grd


gmt begin %file% png

echo Map is created ...
REM --- create grid from given ascii file ---
gmt xyz2grd %ascii_file% -R%region% -r -I%grid_res% -G%grid_file% -V

REM --- create map ---
echo Map is created:

REM --- create color palette table ---
REM gmt makecpt -C%color_palette% -T-42/-28/0.001 -Z
gmt grd2cpt %grid_file% -C%color_palette% -Z 

REM --- plot grid ---
gmt grdimage -J%map_proj% -R%region% %grid_file% -Q
gmt psxy %file_poly% -W3,red
gmt coast -Bxa5g5 -Bya5g5 -BWESN+t"%title%" -W0.25p,80/80/80 -Df -V 

REM --- plot legend including color scale ---
echo Legend is created...
REM gmt colorbar -Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V
gmt text -F+cBL+t"Equivalent water heights (EWH)" -N -D6.65c/-1c
gmt text -F+cBL+t"Editors: Christopher Mahn, Silas Teske, Joshua Wolf" -N -D5.15c/-1.5c

gmt end show
@echo on