gmt gmtset FORMAT_GEO_MAP ddd
set ascii_file=/output/ewh_2008-04_filtered_300km.csv
set grid_file=/output/ewh_2008-04_filtered_300km.grd
gmt begin ewh_2008-04_filtered_300km png
gmt xyz2grd %ascii_file% -R-145/-110/45/65 -r -I1 -G%grid_file% -V
gmt grd2cpt %grid_file% -Chaxby -Z
gmt grdimage -JB-130/65/45/65/18c -R-145/-110/45/65 %grid_file% -Q
gmt psxy data/region_polygon.txt -R-145/-110/45/65 -W3,red
gmt coast -Bxa5g5 -Bya5g5 -BWESN+t"British Columbia, Canada" -W0.25p,80/80/80 -Df -N1/1.25p,black -V
gmt text -F+cBL+t"Equivalent water heights (EWH)" -N -D6.65c/-1c
gmt text -F+cBL+t"Editors: Christopher Mahn, Silas Teske, Joshua Wolf" -N -D5.15c/-1.5c
gmt colorbar -Dx0c/-2c+w17c/0.35c+h -B0.5+l"EWH [m]" -V
gmt end