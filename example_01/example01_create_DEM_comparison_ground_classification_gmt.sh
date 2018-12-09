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


DEM2=dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_normal_1m.tif
gdalwarp $DEM2 foo.tif -tap -tr 1 1 -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DEM2
DEM2_NC=${DEM2::4}.nc
gmt grdconvert $DEM2=gd/1/0/-9999 -G${DEM2_NC}
DEM2_GRD_HS=${DEM2::-4}_HS.nc
if [ ! -e $DEM2_GRD_HS ]
then
    echo "generate hillshade $DEM2_GRD_HS"
    gmt grdgradient $DEM2_NC -Nt0.8 -Es75/55 -G$DEM2_GRD_HS
fi

DEM3=dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_lasground_new_1m.tif
gdalwarp $DEM3 foo.tif -tap -tr 1 1 -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DEM3
DEM3_NC=${DEM3::4}.nc
gmt grdconvert $DEM3=gd/1/0/-9999 -G${DEM3_NC}
DEM3_GRD_HS=${DEM3::-4}_HS.nc
if [ ! -e $DEM3_GRD_HS ]
then
    echo "generate hillshade $DEM3_GRD_HS"
    gmt grdgradient $DEM3_NC -Nt0.8 -Es75/55 -G$DEM3_GRD_HS
fi

gmt grdmath $DEM1 $DEM2 SUB = dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_1m_minus_Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_normal_1m.nc
gmt grdmath $DEM1 $DEM3 SUB = dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_1m_minus_Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_lasground_new_1m.nc

#YOU WILL NEED TO ADJUST THE FOLLOWING PARAMETERS
UTM_ZONE="11N"
J_SCALE="1:3000"
DD_SPACING="0.1m"
LENGTHScale_length_km=0.1

POSTSCRIPT_BASENAME="Pozo_cat16"
POSTSCRIPT1=${POSTSCRIPT_BASENAME}_lasground_adjusted.ps
DEM_CPT=dem2_gray.cpt
gmt grd2cpt ${DEM1_NC} -E25 -N -Cgray >$DEM_CPT
gmt grdimage $DEM1_GRD_HS -I$DEM1_GRD_HS -Jx${J_SCALE} -C$DEM_CPT -R$DEM1 -Q -B+t"Cat16 lasground adjusted" -Xc -Yc -E300 -K >$POSTSCRIPT1
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT1
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT1
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT1
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT1 ${POSTSCRIPT1::-3}.png

POSTSCRIPT_BASENAME="Pozo_cat16"
POSTSCRIPT2=${POSTSCRIPT_BASENAME}_lasground_normal.ps
DEM_CPT=dem2_gray.cpt
gmt grd2cpt ${DEM2_NC} -E25 -N -Cgray >$DEM_CPT
gmt grdimage $DEM2_GRD_HS -I$DEM2_GRD_HS -Jx${J_SCALE} -C$DEM_CPT -R$DEM1 -Q -B+t"Cat16 lasground normal" -Xc -Yc -E300 -K >$POSTSCRIPT2
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT2
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT2
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT2
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT2 ${POSTSCRIPT2::-3}.png

POSTSCRIPT_BASENAME="Pozo_cat16"
POSTSCRIPT3=${POSTSCRIPT_BASENAME}_lasground_new.ps
DEM_CPT=dem2_gray.cpt
gmt grd2cpt ${DEM3_NC} -E25 -N -Cgray >$DEM_CPT
gmt grdimage $DEM3_GRD_HS -I$DEM3_GRD_HS -Jx${J_SCALE} -C$DEM_CPT -R$DEM1 -Q -B+t"Cat16 lasground new" -Xc -Yc -E300 -K >$POSTSCRIPT3
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT3
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT3
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT3
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT3 ${POSTSCRIPT3::-3}.png

convert ${POSTSCRIPT1::-3}.png ${POSTSCRIPT2::-3}.png ${POSTSCRIPT3::-3}.png -fuzz 1% -trim +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Cat16_lasground_comparison.png

POSTSCRIPT_BASENAME="Pozo_cat16"
POSTSCRIPT4=${POSTSCRIPT_BASENAME}_lasground_minus_lasground_normal.ps
DEM_DIFF_CPT=dem_diff.cpt
gmt makecpt -T-1/1/0.1 -D  -Cseis >$DEM_DIFF_CPT
gmt grdimage dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_1m_minus_Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_normal_1m.nc -I$DEM1_GRD_HS -Jx${J_SCALE} -C$DEM_DIFF_CPT -R$DEM1 -Q -B+t"lasground minus lasground_normal" -Xc -Yc -E300 -K >$POSTSCRIPT4
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT4
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT4
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DEM_DIFF_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+lMeter --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT4
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT4
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT4 ${POSTSCRIPT4::-3}.png

POSTSCRIPT_BASENAME="Pozo_cat16"
POSTSCRIPT5=${POSTSCRIPT_BASENAME}_lasground_minus_lasground_new.ps
DEM_DIFF_CPT=dem_diff.cpt
gmt grdimage dtm_interp/Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_1m_minus_Pozo_USGS_UTM11_NAD83_cat16_clg_cl2_lasground_new_1m.nc -I$DEM1_GRD_HS -Jx${J_SCALE} -C$DEM_DIFF_CPT -R$DEM1 -Q -B+t"lasground minus lasground_new" -Xc -Yc -E300 -K >$POSTSCRIPT5
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT5
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT5
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DEM_DIFF_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+lMeter --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT5
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT5
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT5 ${POSTSCRIPT5::-3}.png

POSTSCRIPT_BASENAME="Pozo_cat16"
POSTSCRIPT6=${POSTSCRIPT_BASENAME}_DEM.ps
DEM_CPT=dem2_color.cpt
gmt grd2cpt ${DEM3_NC} -E35 -D -N -Cdem3 >$DEM_CPT
gmt grdimage ${DEM3_NC} -I$DEM1_GRD_HS -Jx${J_SCALE} -C$DEM_CPT -R$DEM1 -Q -B+t"DEM" -Xc -Yc -E300 -K >$POSTSCRIPT6
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT6
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mmF >> $POSTSCRIPT6
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DEM_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+lMeter --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT6
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT6
convert -rotate 90 -quality 100 -density 300 $POSTSCRIPT6 ${POSTSCRIPT6::-3}.png

convert ${POSTSCRIPT6::-3}.png ${POSTSCRIPT4::-3}.png ${POSTSCRIPT5::-3}.png -fuzz 1% -trim +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Cat16_lasground_diff.png
