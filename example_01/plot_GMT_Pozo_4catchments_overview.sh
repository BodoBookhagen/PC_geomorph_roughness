#!/bin/bash
gmt gmtset MAP_FRAME_PEN    3
gmt gmtset MAP_FRAME_WIDTH    0.1
gmt gmtset MAP_FRAME_TYPE     plain
gmt gmtset FONT_TITLE    Helvetica-Bold
gmt gmtset FONT_LABEL    Helvetica-Bold 14p
gmt gmtset PS_PAGE_ORIENTATION    landscape
gmt gmtset PS_MEDIA    A4
gmt gmtset FORMAT_GEO_MAP    D
gmt gmtset MAP_DEGREE_SYMBOL degree
gmt gmtset PROJ_LENGTH_UNIT cm

#gmt grdconvert dtm_interp/Pozo_USGS_UTM11_NAD83_g_1m.tif=gd/1/0/-9999 dtm_interp/Pozo_USGS_UTM11_NAD83_g_1m.nc
POZO_DEM=dtm_interp/Pozo_USGS_UTM11_NAD83_g_1m.nc
POZO_DEM_HS=${POZO_DEM::-3}_HS.nc
if [ ! -e $POZO_DEM_HS ]
then
    echo "generate hillshade $DEM_GRID_HS"
    gmt grdgradient $POZO_DEM -Nt0.8 -Es75/55 -G$POZO_DEM_HS
fi

#Boundary (polygon) of SCI: /home/bodo/Dropbox/California/SCI/SCI_boundary_clip_UTM11N_NAD83.shp
#convert to GMT format
#ogr2ogr -f GMT SCI_boundary_clip_UTM11N_NAD83.gmt /home/bodo/Dropbox/California/SCI/SCI_boundary_clip_UTM11N_NAD83.shp
SCI_BOUNDARY=/raid-cachi/bodo/Dropbox/California/SCI/SCI_boundary_clip_UTM11N_NAD83.gmt

#Pozo catchment
#ogr2ogr -f GMT SCI_Pozo_catchment_UTM11N_NAD83.gmt /home/bodo/Dropbox/California/SCI/SCI_Pozo_catchment_UTM11N_NAD83.shp
POZO_BOUNDARY=/raid-cachi/bodo/Dropbox/California/SCI/SCI_Pozo_catchment_UTM11N_NAD83.gmt

#ogr2ogr -f GMT Pozo_DTM_noveg_UTM11_NAD83_cat16.gmt /raid-cachi/bodo/Dropbox/soft/github/PC_geomorph_roughness/example_01/Pozo_DTM_noveg_UTM11_NAD83_cat16.shp
EXAMPLE01=/raid-cachi/bodo/Dropbox/California/SCI/Pozo/Pozo_DTM_noveg_UTM11_NAD83_cat16.gmt

#ogr2ogr -f GMT SC12.gmt /raid-cachi/bodo/Dropbox/soft/github/PC_geomorph_roughness/example_02/SC12.shp
EXAMPLE02=/raid-cachi/bodo/Dropbox/California/SCI/Pozo/SC12.gmt

#ogr2ogr -f GMT Pozo_DTM_noveg_UTM11_NAD83_cat17.gmt /raid-cachi/bodo/Dropbox/soft/github/PC_geomorph_roughness/example_03/Pozo_DTM_noveg_UTM11_NAD83_cat17.shp
EXAMPLE03=/raid-cachi/bodo/Dropbox/California/SCI/Pozo/Pozo_DTM_noveg_UTM11_NAD83_cat17.gmt

#ogr2ogr -f GMT WestCanada.gmt /raid-cachi/bodo/Dropbox/soft/github/PC_geomorph_roughness/example_04/WestCanada.shp
EXAMPLE04=/raid-cachi/bodo/Dropbox/California/SCI/Pozo/WestCanada.gmt

echo "`tail -n +12 $EXAMPLE01 | awk '{ total += $1 } END { print total/NR }'` `tail -n +12 $EXAMPLE01 | awk '{ total += $2 } END { print total/NR }'` BR Example01" > catchment1.txt
echo "`tail -n +12 $EXAMPLE02 | awk '{ total += $1 } END { print total/NR }'` `tail -n +12 $EXAMPLE02 | awk '{ total += $2 } END { print total/NR }'` BR Example02" > catchment2.txt
echo "`tail -n +12 $EXAMPLE03 | awk '{ total += $1 } END { print total/NR }'` `tail -n +12 $EXAMPLE03 | awk '{ total += $2 } END { print total/NR }'` BR Example03" > catchment3.txt
echo "`tail -n +12 $EXAMPLE04 | awk '{ total += $1 } END { print total/NR }'` `tail -n +12 $EXAMPLE04 | awk '{ total += $2 } END { print total/NR }'` TR Example04" > catchment4.txt

#Preparing stream network:
#extracted stream from Matlab scripts (Neely et al., 2017) stored in SCI_1m_noveg_DTM_UTM11_NAD83_shapefiles.zip
#unzip  SCI_1m_noveg_DTM_UTM11_NAD83_shapefiles.zip
#SCI_FAC=shapefiles/SCI_1m_noveg_DTM_UTM11_NAD83_all_MS_proj.shp

#width of map in cm:
OVERVIEW_WIDTH=10
OVERVIEW_SCALE=1:25000
OVERVIEW_REGION=$POZO_DEM
OVERVIEW_XSTEPS=0.04
OVERVIEW_YSTEPS=0.04

echo "Creating map for Pozo"
POSTSCRIPT1=Pozo_DEM_examples_map.ps
TITLE="Pozo catchment, Santa Cruz Island, California, 1-m Lidar DEM"
CPT="dem2_color.cpt"
gmt grd2cpt $POZO_DEM -E25 -Cdem2 > $CPT
#additional color tables are: -Cdem1, -Cdem3, -Cdem4

gmt grdimage -Q -R$OVERVIEW_REGION $POZO_DEM -I$POZO_DEM_HS -C$CPT -Jx$OVERVIEW_SCALE -V -K --COLOR_BACKGROUND=white > $POSTSCRIPT1 
gmt psxy -Wthin,black -R -J $POZO_BOUNDARY -O -K >> $POSTSCRIPT1
gmt psxy -Wfat,blue -R -J $EXAMPLE01 -O -K >> $POSTSCRIPT1
gmt psxy -Wfat,red -R -J $EXAMPLE02 -O -K >> $POSTSCRIPT1
gmt psxy -Wfat,green -R -J $EXAMPLE03 -O -K >> $POSTSCRIPT1
gmt psxy -Wfat,black -R -J $EXAMPLE04 -O -K >> $POSTSCRIPT1
gmt pstext catchment1.txt -R$OVERVIEW_REGION -Jx$OVERVIEW_SCALE -Gwhite -Wthin,black -D-2c/+0.5cv -F+f16p,Helvetica,blue+j -V -O -K  >> $POSTSCRIPT1
gmt pstext catchment2.txt -R -Jx$OVERVIEW_SCALE -Gwhite -Wthin,black -D-1c/+0.5cv -F+f16p,Helvetica,red+j -V -O -K  >> $POSTSCRIPT1
gmt pstext catchment3.txt -R -Jx$OVERVIEW_SCALE -Gwhite -Wthin,black -D4c/-1cv -F+f16p,Helvetica,green+j -V -O -K  >> $POSTSCRIPT1
gmt pstext catchment4.txt -R -Jx$OVERVIEW_SCALE -Gwhite -Wthin,black -D-1c/+2cv -F+f16p,Helvetica,black+j -V -O -K  >> $POSTSCRIPT1

gmt pscoast -R$OVERVIEW_REGION -Ju11N/$OVERVIEW_SCALE -V -N1 -K -O -Df -Bx1m -By1m --FONT_ANNOT_PRIMARY=10p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT1
#add length scale
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=9p -LjRB+c32:34N+f+w1k+l1:25,000+u+o0.2i --FONT_LABEL=10p >> $POSTSCRIPT1
gmt psscale -R$OVERVIEW_REGION -V -J -DjTRC+o1.55c/0.3c/+w6c/0.3c+h -C$CPT -I -F+gwhite+r1p+pthin,black -Baf -By+lMeter --FONT=10p --FONT_ANNOT_PRIMARY=10p -O -K >> $POSTSCRIPT1

#creating map insert
RANGE=-125/-113/30/45
ticks=a5f5/a5f5
PROJECTION=M2
gmt psbasemap -Bp2dSEnw -Ju11N/1:3000000 -R$RANGE -P -K -V -X170c -Y10c --FONT_ANNOT_PRIMARY=9p --FONT_LABEL=10p >> $POSTSCRIPT1
gmt pscoast -J -R$RANGE -Df -A10 -Ggray -Na -P -K -O -V >> $POSTSCRIPT1
#adding red star around SCI:
echo "-119.9 33.9" >SCI_location.txt
gmt psxy -R$RANGE -J -Sa5c -Gred -O -K -V -P SCI_location.txt >> $POSTSCRIPT1

#convert to JPG
convert -rotate 90 -quality 100 -density 150 $POSTSCRIPT1 ${POSTSCRIPT1::-3}.jpg
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT1 ${POSTSCRIPT1::-3}.png


