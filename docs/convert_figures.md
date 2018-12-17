# prepare figures for PC_geomorph_roughness_manual
convert
convert -crop 600x212+1990+1990 +repage Pozo_Blanca_cl2_1.0m_DEM.png Pozo_Blanca_cl2_1.0m_pointdensity.png -border 0x25 +append org_DEM_PC_density_1m.png

convert Blanca_cl2_1.0m_DEM.png Blanca_cl2_1.0m_pointdensity.png -fuzz 1% -trim -border 0x25 +repage +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/org_DEM_PC_density_1m.png
convert Blanca_cl2_k8_1.0m_DEM.png Blanca_cl2_k8_1.0m_pointdensity.png -fuzz 1% -trim -border 0x25 +repage +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/k8_DEM_PC_density_1m.png
convert Blanca_cl2_p05_1.0m_DEM.png Blanca_cl2_p05_1.0m_pointdensity.png -fuzz 1% -trim -border 0x25 +repage +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/p05_DEM_PC_density_1m.png
convert Blanca_cl2_p08_1.0m_DEM.png Blanca_cl2_p08_1.0m_pointdensity.png -fuzz 1% -trim -border 0x25 +repage +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/p08_DEM_PC_density_1m.png

convert Blanca_cl2_k8_1.0m_pointdensity.png Blanca_cl2_k8_2.0m_pointdensity.png -fuzz 1% -trim -border 0x25 +repage +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/k8_DEM_PC_density_1m_2m.png
convert Blanca_cl2_k8_3.0m_pointdensity.png Blanca_cl2_k8_4.0m_pointdensity.png -fuzz 1% -trim -border 0x25 +repage +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/k8_DEM_PC_density_3m_4m.png
convert ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/k8_DEM_PC_density_1m_2m.png ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/k8_DEM_PC_density_3m_4m.png -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/k8_DEM_PC_density_1m_2m_3m_4m.png

convert -crop 2363x1697+180+330 PC_Pozo_color.jpg PC_Pozo_intensity.jpg +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/PC_Pozo_color_intensity.png

convert -crop 1422x1728+830+324 catchment16_intensity.jpg catchment16_classification.jpg +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/PC_cat16_color_intensity.png

convert -crop 1524x1549+732+479 cat16_color.png cat16_intensity.png cat16_classification.png +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/PC_cat16_lower_color_intensity_classification.png

convert -crop 1128x1835+1058+198 cat16_intensity.png cat16_smrf.png cat16_pmf.png +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/PC_cat16_intensity_smrf_pmf_classification.png

convert -crop 3015x1806+188+204 PC_zoom_org.png PC_zoom_radius10cm.png PC_zoom_radius50cm.png -border 0x25 +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Cat16_PC_zoom_poisson_sampling.png

convert -crop 1212x1891+1049+176 zoom_Z.png zoom_slope_LSTsquared.png +repage +append zoom_Z_slope_LSTsquared.png
convert -crop 1212x1891+1049+176 zoom_dZ_IQR.png zoom_curvature.png +repage +append zoom_dZ_IQ_curvature.png
convert zoom_Z_slope_LSTsquared.png zoom_dZ_IQ_curvature.png -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Cat16_PC_zoom_Z_slope_IQR_results.png

convert -crop 1326x2025+945+54 mean_elevation.png iqr.png +repage +append zoom_Z_IQR.png
convert -crop 1326x2025+945+54 p1_slope.png p2_slope.png +repage +append zoom_slopes.png
convert zoom_Z_IQR.png zoom_slopes.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Cat16_PC_zoom_Z_slope_IQR_results_2nd.png

convert Blanca_cl2_k8_1.0m_pointdensity.png Blanca_cl2_k8_2.0m_pointdensity.png -fuzz 1% -trim -border 0x25 +repage +append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/k8_DEM_PC_density_1m_2m.png

convert Ex01_cl2_1.0m_nrlidar_pts.png Ex01_cl2_2.0m_nrlidar_pts.png -fuzz 1% -trim +repage +append NrLidarPots_1_2m.png
convert Ex01_cl2_3.0m_nrlidar_pts.png Ex01_cl2_5.0m_nrlidar_pts.png -fuzz 1% -trim +repage +append NrLidarPots_3_5m.png
convert NrLidarPots_1_2m.png NrLidarPots_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/NrLidarPots_1_2_3_5m.png

convert Ex01_cl2_1.0m_dz_IQR.png Ex01_cl2_2.0m_dz_IQR.png -fuzz 1% -trim +repage +append dz_IQR_1_2m.png
convert Ex01_cl2_3.0m_dz_IQR.png Ex01_cl2_5.0m_dz_IQR.png -fuzz 1% -trim +repage +append dz_IQR_3_5m.png
convert dz_IQR_1_2m.png dz_IQR_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/dz_IQR_1_2_3_5m.png

convert Ex01_cl2_1.0m_SlopeLSTSQ_P2.png Ex01_cl2_2.0m_SlopeLSTSQ_P2.png -fuzz 1% -trim +repage +append SlopeP2_1_2m.png
convert Ex01_cl2_3.0m_SlopeLSTSQ_P2.png Ex01_cl2_5.0m_SlopeLSTSQ_P2.png -fuzz 1% -trim +repage +append SlopeP2_3_5m.png
convert SlopeP2_1_2m.png SlopeP2_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/SMRF_slopeP2_1_2_3_5m.png

convert Ex01_cl2_1.0m_SlopeLSTSQ_P1.png Ex01_cl2_2.0m_SlopeLSTSQ_P1.png -fuzz 1% -trim +repage +append SlopeP1_1_2m.png
convert Ex01_cl2_3.0m_SlopeLSTSQ_P1.png Ex01_cl2_5.0m_SlopeLSTSQ_P1.png -fuzz 1% -trim +repage +append SlopeP1_3_5m.png
convert SlopeP1_1_2m.png SlopeP1_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/SMRF_slopeP1_1_2_3_5m.png

convert Ex01_cl2_1.0m_P1_rmse.png Ex01_cl2_2.0m_P1_rmse.png -fuzz 1% -trim +repage +append RMSE_P1_1_2m.png
convert Ex01_cl2_3.0m_P1_rmse.png Ex01_cl2_5.0m_P1_rmse.png -fuzz 1% -trim +repage +append RMSE_P1_3_5m.png
convert RMSE_P1_1_2m.png RMSE_P1_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/SMRF_RMSE_P1_1_2_3_5m.png

convert Ex01_cl2_1.0m_P2_rmse.png Ex01_cl2_2.0m_P2_rmse.png -fuzz 1% -trim +repage +append RMSE_P2_1_2m.png
convert Ex01_cl2_3.0m_P2_rmse.png Ex01_cl2_5.0m_P2_rmse.png -fuzz 1% -trim +repage +append RMSE_P2_3_5m.png
convert RMSE_P2_1_2m.png RMSE_P2_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/SMRF_RMSE_P2_1_2_3_5m.png

#LASGROUND
convert Ex01_cl2_1.0m_SlopeLSTSQ_P2.png Ex01_cl2_2.0m_SlopeLSTSQ_P2.png -fuzz 1% -trim +repage +append SlopeP2_1_2m.png
convert Ex01_cl2_3.0m_SlopeLSTSQ_P2.png Ex01_cl2_5.0m_SlopeLSTSQ_P2.png -fuzz 1% -trim +repage +append SlopeP2_3_5m.png
convert SlopeP2_1_2m.png SlopeP2_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/lasground_slopeP2_1_2_3_5m.png

convert Ex01_cl2_1.0m_SlopeLSTSQ_P1.png Ex01_cl2_2.0m_SlopeLSTSQ_P1.png -fuzz 1% -trim +repage +append SlopeP1_1_2m.png
convert Ex01_cl2_3.0m_SlopeLSTSQ_P1.png Ex01_cl2_5.0m_SlopeLSTSQ_P1.png -fuzz 1% -trim +repage +append SlopeP1_3_5m.png
convert SlopeP1_1_2m.png SlopeP1_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/lasground_slopeP1_1_2_3_5m.png

convert Ex01_cl2_1.0m_P1_rmse.png Ex01_cl2_2.0m_P1_rmse.png -fuzz 1% -trim +repage +append RMSE_P1_1_2m.png
convert Ex01_cl2_3.0m_P1_rmse.png Ex01_cl2_5.0m_P1_rmse.png -fuzz 1% -trim +repage +append RMSE_P1_3_5m.png
convert RMSE_P1_1_2m.png RMSE_P1_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/lasground_RMSE_P1_1_2_3_5m.png

convert Ex01_cl2_1.0m_P2_rmse.png Ex01_cl2_2.0m_P2_rmse.png -fuzz 1% -trim +repage +append RMSE_P2_1_2m.png
convert Ex01_cl2_3.0m_P2_rmse.png Ex01_cl2_5.0m_P2_rmse.png -fuzz 1% -trim +repage +append RMSE_P2_3_5m.png
convert RMSE_P2_1_2m.png RMSE_P2_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/lasground_RMSE_P2_1_2_3_5m.png

convert Ex01_cl2_1.0m_dz_IQR.png Ex01_cl2_2.0m_dz_IQR.png -fuzz 1% -trim +repage +append dz_IQR_1_2m.png
convert Ex01_cl2_3.0m_dz_IQR.png Ex01_cl2_5.0m_dz_IQR.png -fuzz 1% -trim +repage +append dz_IQR_3_5m.png
convert dz_IQR_1_2m.png dz_IQR_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/lasground_dz_IQR_1_2_3_5m.png


# Example 03
#LASGROUND
convert Ex03_cl2_1.0m_SlopeLSTSQ_P2.png Ex03_cl2_2.0m_SlopeLSTSQ_P2.png -fuzz 1% -trim +repage +append SlopeP2_1_2m.png
convert Ex03_cl2_3.0m_SlopeLSTSQ_P2.png Ex03_cl2_5.0m_SlopeLSTSQ_P2.png -fuzz 1% -trim +repage +append SlopeP2_3_5m.png
convert SlopeP2_1_2m.png SlopeP2_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Ex03_lasground_slopeP2_1_2_3_5m.png

convert Ex03_cl2_1.0m_SlopeLSTSQ_P1.png Ex03_cl2_2.0m_SlopeLSTSQ_P1.png -fuzz 1% -trim +repage +append SlopeP1_1_2m.png
convert Ex03_cl2_3.0m_SlopeLSTSQ_P1.png Ex03_cl2_5.0m_SlopeLSTSQ_P1.png -fuzz 1% -trim +repage +append SlopeP1_3_5m.png
convert SlopeP1_1_2m.png SlopeP1_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Ex03_lasground_slopeP1_1_2_3_5m.png

convert Ex03_cl2_1.0m_P1_rmse.png Ex03_cl2_2.0m_P1_rmse.png -fuzz 1% -trim +repage +append RMSE_P1_1_2m.png
convert Ex03_cl2_3.0m_P1_rmse.png Ex03_cl2_5.0m_P1_rmse.png -fuzz 1% -trim +repage +append RMSE_P1_3_5m.png
convert RMSE_P1_1_2m.png RMSE_P1_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Ex03_lasground_RMSE_P1_1_2_3_5m.png

convert Ex03_cl2_1.0m_P2_rmse.png Ex03_cl2_2.0m_P2_rmse.png -fuzz 1% -trim +repage +append RMSE_P2_1_2m.png
convert Ex03_cl2_3.0m_P2_rmse.png Ex03_cl2_5.0m_P2_rmse.png -fuzz 1% -trim +repage +append RMSE_P2_3_5m.png
convert RMSE_P2_1_2m.png RMSE_P2_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Ex03_lasground_RMSE_P2_1_2_3_5m.png

convert Ex03_cl2_1.0m_dz_IQR.png Ex03_cl2_2.0m_dz_IQR.png -fuzz 1% -trim +repage +append dz_IQR_1_2m.png
convert Ex03_cl2_3.0m_dz_IQR.png Ex03_cl2_5.0m_dz_IQR.png -fuzz 1% -trim +repage +append dz_IQR_3_5m.png
convert dz_IQR_1_2m.png dz_IQR_3_5m.png +repage -append ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/Ex03_lasground_dz_IQR_1_2_3_5m.png

convert -crop 2082x1847+404+186 pozo.png +repage ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/pozotitle.png

#Example 04
cp -rv Ex04_cl2_1.0m_2panel_DEMs.png Ex04_cl2_1.0m_2panel_SLPs.png Ex04_cl2_1.0m_2panel_RMSE.png ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/
cp -rv Ex04_cl2_3.0m_2panel_RMSE.png  Ex04_cl2_3.0m_2panel_SLPs.png ~/Dropbox/soft/github/PC_geomorph_roughness/docs/figs/
