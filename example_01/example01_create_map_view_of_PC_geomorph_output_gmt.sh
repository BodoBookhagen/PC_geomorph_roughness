#!/bin/bash
gmt gmtset MAP_FRAME_PEN    2
gmt gmtset MAP_FRAME_WIDTH    0.1
gmt gmtset MAP_FRAME_TYPE     plain
gmt gmtset FONT_TITLE    Helvetica-Bold 18p
gmt gmtset FONT_LABEL    Helvetica-Bold 14p
gmt gmtset PS_PAGE_ORIENTATION    landscape
gmt gmtset PS_MEDIA    A4
gmt gmtset FORMAT_GEO_MAP    D
gmt gmtset MAP_DEGREE_SYMBOL degree
gmt gmtset PROJ_LENGTH_UNIT cm

#The following parameters are called by pc_geomorph
TITLE="$1"
POSTSCRIPT_BASENAME="$2"
SHAPEFILE="$3"
DEM_MN_GRD="$4"
DEM_INTERP_GRD="$5"
SLP_LSTSQ_P1_GRD="$6"
SLP_LSTSQ_P2_GRD="$7"
SLP_LSTSQ_P1_RSME_GRD="$8"
SLP_LSTSQ_P2_RSME_GRD="$9"
CONTOUR_CURV_P2_GRD="${10}"
TANGENTIAL_CURV_P2_GRD="${11}"
PROFILE_CURV_P2_GRD="${12}"
NRLIDARPTS_GRD="${13}"
DZ_STDDEV_GRD="${14}"
DZ_IQR_GRD="${15}"
DZ_R9010_GRD="${16}"
DELTA_DEM="${17}"
DEM_RESOLUTION="${18}"

#YOU WILL NEED TO ADJUST THE FOLLOWING PARAMETERS
UTM_ZONE="11N"
J_SCALE="1:3000"
DD_SPACING="0.1m"
LENGTHScale_length_km=0.1

## Prepare DEMs and other GRIDS
DEM_MN_GRD_NC=${DEM_MN_GRD::-4}.nc
gdalwarp $DEM_MN_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DEM_MN_GRD
gmt grdconvert $DEM_MN_GRD=gd/1/0/-9999 -G${DEM_MN_GRD_NC}
DEM_MN_GRD_HS=${DEM_MN_GRD::-4}_HS.nc
if [ ! -e $DEM_MN_GRD_HS ]
then
    echo "generate hillshade $DEM_MN_GRD_HS"
    gmt grdgradient $DEM_MN_GRD_NC -Nt0.8 -Es75/55 -G$DEM_MN_GRD_HS
fi

DEM_INTERP_GRD_NC=${DEM_INTERP_GRD::-4}_${DEM_RESOLUTION}m.nc
gdalwarp $DEM_INTERP_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DEM_INTERP_GRD
gmt grdconvert $DEM_INTERP_GRD=gd/1/0/-9999 -G${DEM_INTERP_GRD_NC}
DEM_INTERP_GRD_HS=${DEM_INTERP_GRD_NC::-3}_HS.nc
if [ ! -e $DEM_INTERP_GRD_HS ]
then
    echo "generate hillshade $DEM_INTERP_GRD_HS"
    gmt grdgradient $DEM_INTERP_GRD_NC -Nt0.8 -Es75/55 -G$DEM_INTERP_GRD_HS
fi

gmt grdmath $DEM_INTERP_GRD_NC $DEM_MN_GRD_NC SUB = $DELTA_DEM

## Prepare SLOPE, Curvature, etc. GRIDSDEM_INTERP_GRD_NC=${DEM_INTERP_GRD::4}.nc
SLP_LSTSQ_P1_GRD_NC=${SLP_LSTSQ_P1_GRD::-4}.nc
gdalwarp $SLP_LSTSQ_P1_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $SLP_LSTSQ_P1_GRD
gmt grdconvert ${SLP_LSTSQ_P1_GRD}=gd/1/0/-9999 -G${SLP_LSTSQ_P1_GRD_NC}

SLP_LSTSQ_P1_RSME_GRD_NC=${SLP_LSTSQ_P1_RSME_GRD::-4}.nc
gdalwarp $SLP_LSTSQ_P1_RSME_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $SLP_LSTSQ_P1_RSME_GRD
gmt grdconvert ${SLP_LSTSQ_P1_RSME_GRD}=gd/1/0/-9999 -G${SLP_LSTSQ_P1_RSME_GRD_NC}

SLP_LSTSQ_P2_GRD_NC=${SLP_LSTSQ_P2_GRD::-4}.nc
gdalwarp $SLP_LSTSQ_P2_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $SLP_LSTSQ_P2_GRD
gmt grdconvert ${SLP_LSTSQ_P2_GRD}=gd/1/0/-9999 -G${SLP_LSTSQ_P2_GRD_NC}

SLP_LSTSQ_P2_RSME_GRD_NC=${SLP_LSTSQ_P2_RSME_GRD::-4}.nc
gdalwarp $SLP_LSTSQ_P2_RSME_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $SLP_LSTSQ_P2_RSME_GRD
gmt grdconvert ${SLP_LSTSQ_P2_RSME_GRD}=gd/1/0/-9999 -G${SLP_LSTSQ_P2_RSME_GRD_NC}

CONTOUR_CURV_P2_GRD_NC=${CONTOUR_CURV_P2_GRD::-4}.nc
gdalwarp $CONTOUR_CURV_P2_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $CONTOUR_CURV_P2_GRD
gmt grdconvert ${CONTOUR_CURV_P2_GRD}=gd/1/0/-9999 -G${CONTOUR_CURV_P2_GRD_NC}

TANGENTIAL_CURV_P2_GRD_NC=${TANGENTIAL_CURV_P2_GRD::-4}.nc
gdalwarp $TANGENTIAL_CURV_P2_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $TANGENTIAL_CURV_P2_GRD
gmt grdconvert ${TANGENTIAL_CURV_P2_GRD}=gd/1/0/-9999 -G${TANGENTIAL_CURV_P2_GRD_NC}

PROFILE_CURV_P2_GRD_NC=${PROFILE_CURV_P2_GRD::-4}.nc
gdalwarp $PROFILE_CURV_P2_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $PROFILE_CURV_P2_GRD
gmt grdconvert ${PROFILE_CURV_P2_GRD}=gd/1/0/-9999 -G${PROFILE_CURV_P2_GRD_NC}

NRLIDARPTS_GRD_NC=${NRLIDARPTS_GRD::-4}.nc
gdalwarp $NRLIDARPTS_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=2
mv foo.tif $NRLIDARPTS_GRD
gmt grdconvert ${NRLIDARPTS_GRD}=gd/1/0/-9999 -G${NRLIDARPTS_GRD_NC}

DZ_STDDEV_GRD_NC=${DZ_STDDEV_GRD::-4}.nc
gdalwarp $DZ_STDDEV_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DZ_STDDEV_GRD
gmt grdconvert ${DZ_STDDEV_GRD}=gd/1/0/-9999 -G${DZ_STDDEV_GRD_NC}

DZ_IQR_GRD_NC=${DZ_IQR_GRD::-4}.nc
gdalwarp $DZ_IQR_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DZ_IQR_GRD
gmt grdconvert ${DZ_IQR_GRD}=gd/1/0/-9999 -G${DZ_IQR_GRD_NC}

DZ_R9010_GRD_NC=${DZ_R9010_GRD::-4}.nc
gdalwarp $DZ_R9010_GRD foo.tif -tap -tr ${DEM_RESOLUTION} ${DEM_RESOLUTION} -r bilinear -cutline $SHAPEFILE -crop_to_cutline -co COMPRESS=DEFLATE -co ZLEVEL=7 -co predictor=3
mv foo.tif $DZ_R9010_GRD
gmt grdconvert ${DZ_R9010_GRD}=gd/1/0/-9999 -G${DZ_R9010_GRD_NC}

## Prepare Shapefile
SHAPEFILE_GMT=${SHAPEFILE::-4}.gmt
if [ ! -e $SHAPEFILE_GMT ]
then
    echo "generate GMT file $SHAPEFILE_GMT"
    ogr2ogr -f GMT $SHAPEFILE_GMT $SHAPEFILE
fi

## Interpolated DEM (best resolution)
POSTSCRIPT1=${POSTSCRIPT_BASENAME}_DEM.ps
echo " "
echo "Creating file $POSTSCRIPT1"
echo " "
DEM_CPT=dem2_color.cpt
gmt grd2cpt ${DEM_INTERP_GRD_NC} -E25 -N -Cdem4 >$DEM_CPT
gmt grdimage $DEM_INTERP_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$DEM_CPT -R$DEM_INTERP_GRD_NC -Q -B+t"${TITLE} interpolated DEM" -Xc -Yc -E300 -K >$POSTSCRIPT1
gmt psxy $SHAPEFILE_GMT -Wthick,black -L -J -R -O -K >>$POSTSCRIPT1
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF >> $POSTSCRIPT1
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT1
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DEM_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+lMeter --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT1
#creating map insert (only for first map)
RANGE=-122/-118/32/36
PROJECTION=M2
gmt psbasemap -Bp2dSEnw  -Ju${UTM_ZONE}/1:20000000 -R$RANGE -O -K -V -X1c -Y2.5c --FONT_ANNOT_PRIMARY=9p >> $POSTSCRIPT1
gmt pscoast -J -R -Df -Slightblue -I1 -Gdarkgray -Na -K -O -V >> $POSTSCRIPT1
#adding red star at SCI:
echo "-119.7 33.9" >SCI_location.txt
gmt psxy -R -J -Sa0.5c -Gred -O -K -V -P SCI_location.txt >> $POSTSCRIPT1
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT1 ${POSTSCRIPT1::-3}.png

## Delta DEM
POSTSCRIPT2=${POSTSCRIPT_BASENAME}_Delta_DEM.ps
echo " "
echo "Creating file $POSTSCRIPT2"
echo " "
DEM_DIFF_CPT=dem_diff_color.cpt
gmt makecpt -T-2/2/0.1 -D  -Cseis >$DEM_DIFF_CPT
gmt grdimage $DELTA_DEM -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$DEM_DIFF_CPT -R$DELTA_DEM -Q -B+t"${TITLE} Delta DEM" -Xc -Yc -E300 -K >$POSTSCRIPT2
gmt psxy $SHAPEFILE_GMT -Wthick,black -L -J -R -O -K >>$POSTSCRIPT2
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF >> $POSTSCRIPT2
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT2
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DEM_DIFF_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+lMeter --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O >> $POSTSCRIPT2
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT2 ${POSTSCRIPT2::-3}.png

## Slope - P1 - Least SQUARED
POSTSCRIPT3=${POSTSCRIPT_BASENAME}_SlopeLSTSQ_P1.ps
echo " "
echo "Creating file $POSTSCRIPT3"
echo " "
SLOPE_CPT=slope_color.cpt
#gmt makecpt -T0/1/0.01 -D  -Cplasma >$SLOPE_CPT
gmt grd2cpt -E25 $SLP_LSTSQ_P2_GRD_NC -D -Cplasma >$SLOPE_CPT
gmt grdimage $SLP_LSTSQ_P1_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$SLOPE_CPT -R$SLP_LSTSQ_P1_GRD_NC -Q -B+t"${TITLE} Slope-P1 from PC" -Xc -Yc -E300 -K >$POSTSCRIPT3
gmt psxy $SHAPEFILE_GMT -Wthick,black -L -J -R -O -K >>$POSTSCRIPT3
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF >> $POSTSCRIPT3
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT3
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$SLOPE_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"Slope-P1" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O >> $POSTSCRIPT3
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT3 ${POSTSCRIPT3::-3}.png

## Slope - P2 - Least SQUARED
POSTSCRIPT3b=${POSTSCRIPT_BASENAME}_SlopeLSTSQ_P2.ps
echo " "
echo "Creating file $POSTSCRIPT3b"
echo " "
gmt grdimage $SLP_LSTSQ_P2_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$SLOPE_CPT -R$SLP_LSTSQ_P2_GRD_NC -Q -B+t"${TITLE} Slope-P2 from PC" -Xc -Yc -E300 -K >$POSTSCRIPT3b
gmt psxy $SHAPEFILE_GMT -Wthick,black -L -J -R -O -K >>$POSTSCRIPT3b
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF >> $POSTSCRIPT3b
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT3b
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$SLOPE_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"Slope-P2" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O >> $POSTSCRIPT3b
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT3b ${POSTSCRIPT3b::-3}.png

## Slope - P1 - Least SQUARED - RMSE
POSTSCRIPT3rmse=${POSTSCRIPT_BASENAME}_P1_rmse.ps
echo " "
echo "Creating file $POSTSCRIPT3rmse"
echo " "
RMSE_CPT=p1_rmse_color.cpt
#gmt makecpt -T0/0.25/0.01 -D  -Cmagma >$RMSE_CPT
gmt grd2cpt -E25 $SLP_LSTSQ_P2_RSME_GRD_NC -D -Cmagma >$RMSE_CPT
gmt grdimage $SLP_LSTSQ_P1_RSME_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$RMSE_CPT -R$SLP_LSTSQ_P1_RSME_GRD_NC -Q -B+t"${TITLE} P1 RMSE from PC" -Xc -Yc -E300 -K >$POSTSCRIPT3rmse
gmt psxy $SHAPEFILE_GMT -Wthick,black -L -J -R -O -K >>$POSTSCRIPT3rmse
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF >> $POSTSCRIPT3rmse
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT3rmse
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$RMSE_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"P1 RMSE (m)" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O >> $POSTSCRIPT3rmse
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT3rmse ${POSTSCRIPT3rmse::-3}.png

## Slope - P2 - Least SQUARED - RMSE
POSTSCRIPT3brmse=${POSTSCRIPT_BASENAME}_P2_rmse.ps
echo " "
echo "Creating file $POSTSCRIPT3brmse"
echo " "
gmt grdimage $SLP_LSTSQ_P2_RSME_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$RMSE_CPT -R$SLP_LSTSQ_P2_RSME_GRD_NC -Q -B+t"${TITLE} P2 RMSE from PC" -Xc -Yc -E300 -K >$POSTSCRIPT3brmse
gmt psxy $SHAPEFILE_GMT -Wthick,black -L -J -R -O -K >>$POSTSCRIPT3brmse
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF >> $POSTSCRIPT3brmse
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT3brmse
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$RMSE_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"P2 RMSE (m)" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O >> $POSTSCRIPT3brmse
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT3brmse ${POSTSCRIPT3brmse::-3}.png

POSTSCRIPT4=${POSTSCRIPT_BASENAME}_curvature_contour.ps
echo " "
echo "Creating file $POSTSCRIPT4"
echo " "
CURVATURE_CPT=curvature_color.cpt
gmt makecpt -T-1/1/0.05 -D -Cpolar >$CURVATURE_CPT
#gmt grd2cpt $CONTOUR_CURV_P2_GRD_NC -E25 -D -T- -Cpolar >$CURVATURE_CPT
gmt grdimage $CONTOUR_CURV_P2_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$CURVATURE_CPT -R$NRLIDARPTS_GRD_NC -Q -B+t"${TITLE} Contour Curvature" -Xc -Yc -E300 -K >$POSTSCRIPT4
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT4
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF  >> $POSTSCRIPT4
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT4
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$CURVATURE_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"Curvature (1/m)" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT4
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT4 ${POSTSCRIPT4::-3}.png

POSTSCRIPT4b=${POSTSCRIPT_BASENAME}_curvature_tangential.ps
echo " "
echo "Creating file $POSTSCRIPT4b"
echo " "
gmt grdimage $TANGENTIAL_CURV_P2_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$CURVATURE_CPT -R$NRLIDARPTS_GRD_NC -Q -B+t"${TITLE} Tangential Curvature" -Xc -Yc -E300 -K >$POSTSCRIPT4b
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT4b
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF  >> $POSTSCRIPT4b
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT4b
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$CURVATURE_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"Curvature (1/m)" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT4b
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT4b ${POSTSCRIPT4b::-3}.png

POSTSCRIPT4c=${POSTSCRIPT_BASENAME}_curvature_profile.ps
echo " "
echo "Creating file $POSTSCRIPT4c"
echo " "
gmt grdimage $PROFILE_CURV_P2_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$CURVATURE_CPT -R$NRLIDARPTS_GRD_NC -Q -B+t"${TITLE} Profile Curvature" -Xc -Yc -E300 -K >$POSTSCRIPT4c
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT4c
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF  >> $POSTSCRIPT4c
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT4c
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$CURVATURE_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"Curvature (1/m)" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT4c
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT4c ${POSTSCRIPT4c::-3}.png

## Nr. of Lidar Points
POSTSCRIPT5=${POSTSCRIPT_BASENAME}_nrlidar_pts.ps
echo " "
echo "Creating file $POSTSCRIPT5"
echo " "
NRLIDARPTS_CPT=nr_lidar_points_color.cpt
gmt grd2cpt ${NRLIDARPTS_GRD_NC} -E10 -N -Crainbow >$NRLIDARPTS_CPT
gmt grdimage $NRLIDARPTS_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$NRLIDARPTS_CPT -R$NRLIDARPTS_GRD_NC -Q -B+t"${TITLE} Nr. Lidar Points/area" -Xc -Yc -E300 -K >$POSTSCRIPT5
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT5
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF  >> $POSTSCRIPT5
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT5
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$NRLIDARPTS_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"Nr. of points/area" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT5
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT5 ${POSTSCRIPT5::-3}.png

## Dz Std. Deviation
POSTSCRIPT6=${POSTSCRIPT_BASENAME}_dzstddev.ps
echo " "
echo "Creating file $POSTSCRIPT6"
echo " "
DZ_STDDEV_GRD_CPT=dz_std_dev_color.cpt
#gmt grd2cpt ${DZ_STDDEV_GRD} -E20 -N -Cplasma >$DZ_STDDEV_GRD_CPT
#gmt makecpt -T0/0.5/0.01 -D  -Cviridis >$DZ_STDDEV_GRD_CPT
gmt grd2cpt $DZ_IQR_GRD_NC -E25 -D -Cviridis >$DZ_STDDEV_GRD_CPT
gmt grdimage $DZ_STDDEV_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$DZ_STDDEV_GRD_CPT -R$DZ_STDDEV_GRD_NC -Q -B+t"${TITLE} Dz standard deviation" -Xc -Yc -E300 -K >$POSTSCRIPT6
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT6
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF >> $POSTSCRIPT6
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT6
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DZ_STDDEV_GRD_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"Dz. std. dev. (m)" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT6
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT6 ${POSTSCRIPT6::-3}.png

## Dz 90-10th percentile
POSTSCRIPT7=${POSTSCRIPT_BASENAME}_r9010p.ps
echo " "
echo "Creating file $POSTSCRIPT7"
echo " "
DZ_R9010p_GRD_CPT=$DZ_STDDEV_GRD_CPT
#gmt grd2cpt ${DZ_R9010_GRD_NC} -E20 -N -Cviridis >$DZ_R9010p_GRD_CPT
gmt grdimage $DZ_R9010_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$DZ_R9010p_GRD_CPT -R$DZ_R9010_GRD_NC -Q -B+t"${TITLE} Dz Range 10th-90th perc." -Xc -Yc -E300 -K >$POSTSCRIPT7
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT7
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF  >> $POSTSCRIPT7
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT7
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DZ_R9010p_GRD_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"Dz. 10-90p (m)" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT7
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT7 ${POSTSCRIPT7::-3}.png

## Dz IQR
POSTSCRIPT8=${POSTSCRIPT_BASENAME}_dz_IQR.ps
echo " "
echo "Creating file $POSTSCRIPT8"
echo " "
DZ_IQR_GRD_CPT=$DZ_STDDEV_GRD_CPT
#gmt grd2cpt ${DZ_IQR_GRD_NC} -E20 -N -Clava >$DZ_IQR_GRD_CPT
gmt grdimage $DZ_IQR_GRD_NC -I$DEM_INTERP_GRD_HS -Jx${J_SCALE} -C$DZ_IQR_GRD_CPT -R$DZ_IQR_GRD_NC -Q -B+t"${TITLE} Dz IQR" -Xc -Yc -E300 -K >$POSTSCRIPT8
gmt psxy $SHAPEFILE_GMT -Wthick,black -J -R -O -K >>$POSTSCRIPT8
gmt pscoast -R -Ju${UTM_ZONE}/${J_SCALE} -N1 -K -O -Df -B${DD_SPACING}SWne --FONT_ANNOT_PRIMARY=12p --FORMAT_GEO_MAP=ddd:mm:ssF  >> $POSTSCRIPT8
gmt psbasemap -R -J -O -K --FONT_ANNOT_PRIMARY=10p -LjLB+c34:36N+f+w${LENGTHScale_length_km}k+l${J_SCALE}+u+o0.5c --FONT_LABEL=10p >> $POSTSCRIPT8
gmt psscale -R -J -DjTC+o0.0c/0.3c/+w6c/0.3c+h -C$DZ_IQR_GRD_CPT -I -F+gwhite+r1p+pthin,black -Baf -By+l"Dz. IQR (m)" --FONT=12p --FONT_ANNOT_PRIMARY=12p --MAP_FRAME_PEN=1 --MAP_FRAME_WIDTH=0.1 -O -K >> $POSTSCRIPT8
convert -alpha off -rotate 90 -quality 100 -density 300 $POSTSCRIPT8 ${POSTSCRIPT8::-3}.png

## Combine Curvature
convert ${POSTSCRIPT4::-3}.png ${POSTSCRIPT4b::-3}.png ${POSTSCRIPT4c::-3}.png -fuzz 1% -trim -bordercolor white -border 0x25 +repage -append ${POSTSCRIPT_BASENAME}_3panel_curvature.png
convert ${POSTSCRIPT3rmse::-3}.png ${POSTSCRIPT3brmse::-3}.png -fuzz 1% -trim -bordercolor white -border 0x25 +repage -append ${POSTSCRIPT_BASENAME}_2panel_RMSE.png

##
convert ${POSTSCRIPT1::-3}.png ${POSTSCRIPT2::-3}.png -fuzz 1% -trim -background white -bordercolor white -border 0x25 +repage -append ${POSTSCRIPT_BASENAME}_2panel_DEMs.png
convert ${POSTSCRIPT3::-3}.png ${POSTSCRIPT3b::-3}.png -fuzz 1% -trim -bordercolor white -border 0x25 +repage -append ${POSTSCRIPT_BASENAME}_2panel_SLPs.png

#Combine into 4-panel figure:
convert ${POSTSCRIPT_BASENAME}_2panel_DEMs.png ${POSTSCRIPT_BASENAME}_2panel_SLPs.png -fuzz 1% -trim -bordercolor white -border 0x25 -background white +repage +append ${POSTSCRIPT_BASENAME}_4panel_DEMS_SLP_CURV.png

convert ${POSTSCRIPT5::-3}.png ${POSTSCRIPT6::-3}.png -fuzz 1% -trim -bordercolor white -border 0x25 +repage -append ${POSTSCRIPT_BASENAME}_2panel_NRLIDARPTS_DZ_STDDEV.png
convert ${POSTSCRIPT7::-3}.png ${POSTSCRIPT8::-3}.png -fuzz 1% -trim -bordercolor white -border 0x25 +repage -append ${POSTSCRIPT_BASENAME}_2panel_DZ9010P_IQR.png

#Combine into 4-panel figure:
convert ${POSTSCRIPT_BASENAME}_2panel_NRLIDARPTS_DZ_STDDEV.png ${POSTSCRIPT_BASENAME}_2panel_DZ9010P_IQR.png -alpha off -fuzz 1% -trim -border 0x25 +repage +append ${POSTSCRIPT_BASENAME}_4panel_NRLIDARPTS_DZ_STDDEV_DZ9010P_IQR.png

