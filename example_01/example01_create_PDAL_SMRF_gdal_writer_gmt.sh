#!/bin/bash
gmt gmtset MAP_FRAME_PEN    2
gmt gmtset MAP_FRAME_WIDTH    0.1
gmt gmtset MAP_FRAME_TYPE     plain
gmt gmtset FONT_TITLE    Helvetica-Bold
gmt gmtset FONT_LABEL    Helvetica-Bold 14p
gmt gmtset PS_PAGE_ORIENTATION    landscape
gmt gmtset PS_MEDIA    A4
gmt gmtset FORMAT_GEO_MAP    D
gmt gmtset MAP_DEGREE_SYMBOL degree
gmt gmtset PROJ_LENGTH_UNIT cm

SHAPEFILE=Pozo_DTM_noveg_UTM11_NAD83_cat16.shp
SHAPEFILE_GMT=${SHAPEFILE::-4}.gmt
if [ ! -e $SHAPEFILE_GMT ]
then
    echo "generate GMT file $SHAPEFILE_GMT"
    ogr2ogr -f GMT $SHAPEFILE_GMT $SHAPEFILE
fi

DEM1=dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_1m.tif
gdalwarp $DEM1 foo.tif -tap -tr 1 1 -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DEM1
DEM1_NC=${DEM1::4}.nc
gmt grdconvert $DEM1=gd/1/0/-9999 -G${DEM1_NC}
DEM1_GRD_HS=${DEM1::-4}_HS.nc
if [ ! -e $DEM1_GRD_HS ]
then
    echo "generate hillshade $DEM1_GRD_HS"
    gmt grdgradient $DEM1_NC -Nt0.8 -Es75/55 -G$DEM1_GRD_HS
fi


DEM6=dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_smrf_cl2_Zmean.tif
DEM6_NC=${DEM6::-4}.nc
gdalwarp $DEM6 foo.tif -tap -tr 1 1 -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DEM6
gmt grdconvert $DEM6=gd/1/0/-9999 -G${DEM6_NC}
DEM6_GRD_HS=${DEM6::-4}_HS.nc
if [ ! -e $DEM6_GRD_HS ]
then
    echo "generate hillshade $DEM6_GRD_HS"
    gmt grdgradient $DEM6_NC -Nt0.8 -Es75/55 -G$DEM6_GRD_HS
fi

DEM7=dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_smrf_cl2_intensity_mean.tif
DEM7_NC=${DEM7::-4}.nc
gdalwarp $DEM7 foo.tif -tap -tr 1 1 -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DEM7
gmt grdconvert $DEM7=gd/1/0/-9999 -G${DEM7_NC}


DEM8=dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_lasground_cl2_Zmean.tif
DEM8_NC=${DEM8::-4}.nc
gdalwarp $DEM8 foo.tif -tap -tr 1 1 -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DEM8
gmt grdconvert $DEM8=gd/1/0/-9999 -G${DEM8_NC}
DEM8_GRD_HS=${DEM8::-4}_HS.nc
if [ ! -e $DEM8_GRD_HS ]
then
    echo "generate hillshade $DEM8_GRD_HS"
    gmt grdgradient $DEM8_NC -Nt0.8 -Es75/55 -G$DEM8_GRD_HS
fi

gmt grdmath $DEM1 $DEM6 SUB = dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_1m_minus_Pozo_USGS_UTM11_NAD83_cat16_smrf_cl2_Zmean_1m_clip.nc

gmt grdmath $DEM1 $DEM8 SUB = dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_1m_minus_Pozo_USGS_UTM11_NAD83_cat16_lasground_cl2_Zmean_clip.nc

#YOU WILL NEED TO ADJUST THE FOLLOWING PARAMETERS
UTM_ZONE="11N"
J_SCALE="1:3000"
DD_SPACING="0.1m"
LENGTHScale_length_km=0.1

POSTSCRIPT_BASENAME="Pozo_cat16_PDAL"
POSTSCRIPT7=${POSTSCRIPT_BASENAME}_I.ps
DEM_CPT=dem2_color.cpt
gmt grd2cpt ${DEM7_NC} -E35 -D -N -Cplasma >$DEM_CPT
gmt grdimage ${DEM7_NC} -I$DEM6_GRD_HS -Jx${J_SCALE} -C$DEM_CPT -R$DEM7_NC -Q -B+t"DEM SMRF writer Intensity" -Xc -Yc -E300 -K >$POSTSCRIPT7
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT7
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT7
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DEM_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+lIntensity --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT7
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT7
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT7 ${POSTSCRIPT7::-3}.png

POSTSCRIPT_BASENAME="Pozo_cat16_PDAL"
POSTSCRIPT8=${POSTSCRIPT_BASENAME}_lasground_minus_smrf_p2g.ps
DEM_DIFF_CPT=dem_diff.cpt
gmt makecpt -T-1/1/0.1 -D  -Cseis >$DEM_DIFF_CPT
gmt grdimage dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_1m_minus_Pozo_USGS_UTM11_NAD83_cat16_smrf_cl2_Zmean_1m_clip.nc -I$DEM6_GRD_HS -Jx${J_SCALE} -C$DEM_DIFF_CPT -R$DEM6_NC -Q -B+t"lasground minus SMRF writer" -Xc -Yc -E300 -K >$POSTSCRIPT8
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT8
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT8
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DEM_DIFF_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+lMeter --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT8
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT8
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT8 ${POSTSCRIPT8::-3}.png

POSTSCRIPT_BASENAME="Pozo_cat16_PDAL"
POSTSCRIPT9=${POSTSCRIPT_BASENAME}_lasground_blast2dem_minus_p2g.ps
DEM_DIFF_CPT=dem_diff.cpt
gmt makecpt -T-1/1/0.1 -D  -Cseis >$DEM_DIFF_CPT
gmt grdimage dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_1m_minus_Pozo_USGS_UTM11_NAD83_cat16_lasground_cl2_Zmean_clip.nc -I$DEM8_GRD_HS -Jx${J_SCALE} -C$DEM_DIFF_CPT -R$DEM8_NC -Q -B+t"lasground: blast2dem minus gdal.writer" -Xc -Yc -E300 -K >$POSTSCRIPT9
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT9
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT9
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DEM_DIFF_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+lMeter --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT9
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT9
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT9 ${POSTSCRIPT9::-3}.png

convert ${POSTSCRIPT7::-3}.png ${POSTSCRIPT8::-3}.png ${POSTSCRIPT9::-3}.png -fuzz 1% -trim +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Cat16_smrf_PDAL_intensity_diff.png
