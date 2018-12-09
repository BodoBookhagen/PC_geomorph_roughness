#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 17:22:03 2017

@author: Bodo Bookhagen,    V0.1 Oct-Nov 2017
                            V0.2 Dec 2018


example call for k-nearest neighbor k=5 subsampling:
python -W ignore ~/github/PC_geomorph_roughness/pc_geomorph_roughness.py \
    --inlas Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC12.laz \
    --raster_m_range "1 10 1" \
    --shapefile_clip SC12.shp \
    --epsg_code 26911 --nr_of_cores 0 --create_geotiff 1 --create_gmt 1  \
    --create_shapefiles 0 --create_las 1 \
    --subsample_1m_pc_k 5 \
    2>&1 | tee Pozo_USGS_UTM11_NAD83_all_color_cl2_cat1_pc_geomorph_roughness_subsample_k5_1_10_1.log

example call for homogeneous subsampling based on probability to a fraction of 0.5 points:
python -W ignore ~/github/PC_geomorph_roughness/pc_geomorph_roughness.py \
    --inlas Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC12.laz \
    --raster_m_range "1 10 1" \
    --shapefile_clip SC12.shp \
    --epsg_code 26911 --nr_of_cores 0 --create_geotiff 1 --create_gmt 1 \
    --create_shapefiles 0 --create_las 1 \
    --subsample_1m_pc_p 0.5 \
    2>&1 | tee Pozo_USGS_UTM11_NAD83_all_color_cl2_cat1_pc_geomorph_roughness_subsample_p05_1_10_1.log
 
"""

from laspy.file import File
import copy, glob, time, sys
import numpy as np, os, argparse, pickle, h5py, subprocess, gdal, osr, datetime
from numpy.linalg import svd
from pykdtree.kdtree import KDTree
from scipy import interpolate
from scipy import linalg
import matplotlib as mpl
import matplotlib.cm as cm
import multiprocessing
from multiprocessing import Pool
from skimage import exposure
from itertools import combinations

### Function definitions
def cmdLineParser():
    # Command Line Parsing
    parser = argparse.ArgumentParser(description='PointCloud (PC) processing for DEM statistics. Deriving gridded ground data (elevation and slope) using centroid coordinates. B. Bookhagen (bodo.bookhagen@uni-potsdam.de), V0.2 Dec 2018.')
    # Important and required:
    parser.add_argument('-i', '--inlas', required=True, type=str, default='',  help='LAS/LAZ file with point-cloud data. This file must only contain ground points (class == 2)')
    parser.add_argument('--raster_m', type=float, default=1.0,  help='Raster spacing or diameter for subsampling seed points on LAS/LAZ PC. Usually 0.5 to 10 m, default = 1. Seed points are selected from radii half this diameter. ')
    parser.add_argument('--raster_m_range', type=str, default='',  help='Raster spacing for subsampling seed points on LAS/LAZ PC. Uses a list of ranges with spacing, e.g., --raster_m_range "1 10 1" will create raster files with spatial resolutions of 1 to 10 m in 1 m steps.')
    parser.add_argument('--subsample_1m_pc_k', type=int, default=0,  help='Number of points in radius=0.5m that are randomly subsampled from the full point cloud. This is useful if point-cloud density greatly varies, because statistics calculated for seed points with different point numbers may be biased. If subsample_pc_k > 0 then the point cloud will be homogenized by selecting k=n neighbors for each 1-m seed point. For example, if subsample_pc_k 10, then each 1m seed point will have only 10 neighbors. Subsampled point cloud is written to LAS file.')
    parser.add_argument('--subsample_1m_pc_p', type=float, default=0,  help='Factor to subsample point cloud based on probability. If --subsample_1m_pc_p 0.8, a point cloud with 80%% of the input points is generated and sampling of point cloud is based on density probability. That is, neighboring points for a seed point with a high number of neighbors a sampled less often, than a seed point with fewer neighbors. Will use original points, but creates a reduced point cloud. Calculates probability from 1m seed-point spacing. Subsampled point cloud is written to LAS file.')
    parser.add_argument('--redo_subsample_1m_pc_p', type=int, default=0,  help='Flag to redo random subsampling based on probability. By default, an existing file with a probability is loaded, if you set "--redo_subsample_1m_pc_p true", the random subsampling based on probability will be rerun and stored in a separate LAS file.')
    parser.add_argument('--k_nr_of_neighbors', type=int, default=100,  help='Number of neighbors for dynamic density estimation (k_nr_of_neighbors = 50 by default). Change to lower number for lower-density point clouds to increase procesing speed. For SfM PCs this should be set to 100 or higher.')
    parser.add_argument('--dem_fname', type=str, default='',  help='Filename of DEM to extract point spacing. Used to identify seed-point coordinates. Useful if a DEM exists and one wants to create point-cloud statistics aligned to the DEM grid.')
    parser.add_argument('--shapefile_clip', type=str, default='',  help='Name of shapefile to be used to clip interpolated surfaces. Make sure to give full pathname. This is likely the shapefile that has been previously generated to subset/clip the point-cloud data.')
    parser.add_argument('--epsg_code', type=int, default=26911,  help='EPSG code (integer) to define projection information. This should be the same EPSG code as the input LAS data (no re-projection is included in this code yet) and can be taken from LAS/LAZ input file. Add this to ensure that output shapefile and GeoTIFFs are properly geocoded.')
    parser.add_argument('--create_shapefiles', type=int, default=0,  help='Create point shapefiles in UTM (see --epsg_code) and Geographic-DD projection. These contain all attributes calculated during the processing (default no shapefiles are created: --create_shapefiles 0, set to --create_shapefiles 1 to generate shapefiles).')
    parser.add_argument('--create_geotiff', type=int, default=0,  help='Create interpolated geotif files from PC data (default no: --create_geotiff 0, set to --create_geotiff 1 to generate geotiff files). Note that creating geotiff files may increase processing time')
    parser.add_argument('--create_gmt', type=int, default=0,  help='Create gmt point or vector files for plotting with GMT shapefiles in UTM (see --epsg_code) and Geographic-DD projection. These contain all attributes calculated during the processing (default no: --create_gmt 0, set to --create_gmt 1 to generate GMT files).')
    parser.add_argument('--create_las', type=int, default=0,  help='Create LAS point file from seed points. The color shows mean elevation of the seed points (default no: --create_las 0, set to --create_las 1 to generate LAS files).')
    parser.add_argument('--mean_z_only', type=int, default=0,  help='Calculate mean elevation for grid cell size and no other parameters.')
    parser.add_argument('--nr_of_cores', type=int, default=0,  help='Max. number of cores to use for multi-core processing. Default is to use all cores (--nr_of_cores 0), set to --nr_of_cores 6 to use 6 cores. For some memory-intensive applications, it may be useful to reduce the number of cores.')
    parser.add_argument('--max_nr_of_neighbors_kdtree', type=int, default=100,  help='Setting the max. number of neighbors for KDTree search. This can remain at 100 points for airborne lidar data. For example, if you have a point density of 5 pts/m2, 100 pts are 20 m2. You may want to consider increasing this when using terrestrial lidar data, SfM data, or airborne data with high point densities.')
    parser.add_argument('--pt_lower_threshold', type=int, default=3,  help='Minimum number of points for performing plane fitting and slope normalization (lower point threshold). If there are less than pt_lower_threshold in the seed-point neighborhood (default --pt_lower_threshold 3), a point fitting is not performed and values are set to NaN.')
    parser.add_argument('--create_gmt_maps', type=str, default='',  help='BASH File with GMT commands for plotting maps. Full path and filename is required. Will need to be fine tuned (see example).')
    parser.add_argument('--gmt_title', type=str, default='',  help='GMT title to appear in output map.')
    parser.add_argument('--gmt_basename', type=str, default='',  help='GMT basename for Postscript filenames. ')
    parser.add_argument('--plot_plane_fits', type=int, default=0,  help='Create plots of plane fits for individual seed points with more than plot_plane_fits_nr_points (default = 10) neighborhood points in subdirectory "maps". Mostly for testing purposes, default is off (--plot_plane_fits 0).')
    parser.add_argument('--plot_plane_fits_nr_points', type=int, default=10,  help='Set number of neighborhood points to create plot for seed point. Default is --plot_plane_fits_nr_points 10. You will need to adjust this for larger neighborhood radii.')
    return parser.parse_args()

def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    source: https://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
    NOTE: This does not work well for nosiy point clouds such as lidar point clouds.
    This method is very sensitive to outliers!
    """
    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    try:
        assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    except AssertionError:
        return np.nan, np.nan, np.nan
    
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    plane_normal = svd(M)[0][:,-1]
    d = -ctr.dot(plane_normal)
    z = (-plane_normal[0] * points[0,:] - plane_normal[1] * points[1,:] - d) * 1. /plane_normal[2]
    errors = z - points[2,:]
    residual = np.linalg.norm(errors)

    return ctr, plane_normal, residual

def curvFit_lstsq_polygon(points, order=2):
    """
    Fitting a second order polynom to a point cloud and deriving the curvature in a simplified form.
    We follow Evans, I. S. (1980), An integrated system of terrain analysis and slope mapping, Z. Geomorphol., 36, 274–295.
    More details: https://gis.stackexchange.com/questions/37066/how-to-calculate-terrain-curvature
    """
    
    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    try:
        assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    except AssertionError:
        return np.nan, np.nan, np.nan

    #points_org = np.copy(points)
    #points[:,0] = points[:,0] - points[:,0].mean()
    #points[:,1] = points[:,1] - points[:,1].mean()
    #points[:,2] = points[:,2] - points[:,2].mean()
    points = points.T
    
#    # evaluate it on grid
#    X,Y = np.meshgrid(np.arange(np.nanmin(points[:,0]),np.nanmax(points[:,0]), current_rstep_size/10), np.arange(np.nanmin(points[:,1]),np.nanmax(points[:,1]), current_rstep_size/10))
#    XX = X.flatten()
#    YY = Y.flatten()
    if order == 1:
        # best-fit linear plane
        A = np.c_[points[:,0], points[:,1], np.ones(points.shape[0])]
        C,_,_,_ = linalg.lstsq(A, points[:,2])    # coefficients
        # evaluate it on grid
        #Z_pts = C[0]*X + C[1]*Y + C[2]
        # or expressed using matrix/vector product
        #Z_order1 = np.dot(np.c_[XX, YY, np.ones(XX.shape)], C).reshape(X.shape)
        #slope = np.mean(C[0:2])
        slope = np.sqrt( C[0]**2. + C[1]**2. )
        curv_meanc = np.nan
        curv_tan = np.nan
        curv_profc = np.nan
        Z_pts = C[0]*points[:,0] + C[1]*points[:,1] + C[2]
        errors = points[:,2] - Z_pts
        dZ_residuals = np.linalg.norm(errors)
    elif order == 2:
        # best-fit quadratic curve
        #Z = Dx² + Ey² + Fxy + Gx + Hy + I
        A = np.c_[points[:,0]**2, points[:,1]**2, points[:,0]*points[:,1], points[:,0], points[:,1], np.ones(points.shape[0])]
        C,_,_,_ = linalg.lstsq(A, points[:,2])    # coefficients
        Z_pts = C[0]*points[:,0]**2. + C[1]*points[:,1]**2. + C[2]*points[:,0]*points[:,1] + C[3]*points[:,0] + C[4]*points[:,1] + C[5]

        # evaluate it on a grid
        #Z_order2 = C[0]*X**2. + C[1]*Y**2. + C[2]*X*Y + C[3]*X + C[4]*Y + C[5]
        #Z = Dx² + Ey² + Fxy + Gx + Hy + I
        #Curvature = -2(D + E)
        #Slope = sqrt(G^2 + H ^2)
        curvature = -2 * (C[0] + C[1])
        curv_meanc = curvature
        curv_tan = np.nan
        curv_profc = np.nan
        slope = np.sqrt( C[3]**2 + C[4]**2 )

        #Z_pts = np.dot(np.c_[np.ones(points.shape[0]), points[:,0], points[:,1], points[:,0]*points[:,1], points[:,0]**2, points[:,1]**2], C)
        errors = points[:,2] - Z_pts
        dZ_residuals = np.linalg.norm(errors)
    elif order == 4:
        # best-fit fourth-order polynomial
        #Z = Ax²y² + Bx²y + Cxy² + Dx² + Ey² + Fxy + Gx + Hy + I
        #A = [(Z1 + Z3 + Z7 + Z9) / 4 - (Z2 + Z4 + Z6 + Z8) / 2 + Z5] / L4
        #B = [(Z1 + Z3 - Z7 - Z9) /4 - (Z2 - Z8) /2] / L3
        #C = [(-Z1 + Z3 - Z7 + Z9) /4 + (Z4 - Z6)] /2] / L3
        #D = [(Z4 + Z6) /2 - Z5] / L2
        #E = [(Z2 + Z8) /2 - Z5] / L2
        #F = (-Z1 + Z3 + Z7 - Z9) / 4L2
        #G = (-Z4 + Z6) / 2L
        #H = (Z2 - Z8) / 2L
        #I = Z5
        #Curvature = -2(D + E) * 100        
        A = np.c_[points[:,0]**2. * points[:,1]**2., points[:,0]**2. * points[:,1], points[:,0] * points[:,1]**2., \
                  points[:,0]**2., points[:,1]**2., points[:,0]*points[:,1], \
                  points[:,0], points[:,1], np.ones(points.shape[0])]
        C,_,_,_ = linalg.lstsq(A, points[:,2])    # coefficients
        Z_pts = C[0]*points[:,0]**2. * points[:,1]**2. + C[1]*points[:,0]**2.*points[:,1] + C[2]*points[:,0]*points[:,1]**2. \
            + C[3]*points[:,0]**2. + C[4]*points[:,1]**2. + C[5]*points[:,0]*points[:,1] \
            + C[6]*points[:,0] + C[7]*points[:,1] + C[8]
        fx=C[6]
        fy=C[7]
        fxx=C[3]
        fxy=C[5]
        fyy=C[4]        
        curv_profc = - (fx**2. * fxx + 2*fx * fy * fxy + fy**2. * fyy) / ( (fx**2. + fy**2.) * (1 + fx**2. + fy**2.)**(3/2) );
        #curv_planc = - (fy**2. * fxx + 2*fx * fy * fxy + fx**2. * fyy) / ( (fx**2. + fy**2.)**(3/2) )
        curv_meanc = - ((1 + fy**2) * fxx - 2 * fxy * fx * fy + (1 + fx**2) * fyy) / (2 * (fx**2 + fy**2 + 1)**(3/2) )
        curv_tan = - (fy**2. * fxx + 2*fx * fy * fxy + fx**2 * fyy) / ( (fx**2 + fy**2) * (1 + fx**2 + fy**2)**(1/2) )
 
        #curvature_total = C[4]**2. + 2*C[6]**2. + C[5]**2.
        curvature = -2 * (C[3] + C[4]) * 100
        slope = np.sqrt( fx**2. + fy**2. )
        errors = points[:,2] - Z_pts
        dZ_residuals = np.linalg.norm(errors)
    del A, Z_pts
    return slope, curv_meanc, curv_tan, curv_profc, dZ_residuals, errors, C

def gdal_grid_interpolate(x_coords, y_coords, ncols, nrows,resolution_m, layer_in, zfield, input_vrt, output_grid,radius1=0, radius2=0, grid_datatype='Float32', cliplayer='', interpolation_method='lin'):
    #Use gdal_grid to interpolate from point data to grid
    output_grid2 = output_grid[:-4] + 'b.tif'
    
    if interpolation_method == 'invdist':
        #inverse distance
        interpolation_method_string='invdist:radius1='+str(radius1)+':radius2='+str(radius1)+':angle=0.0:nodata=-9999'

    if interpolation_method == 'invdistnn':
        #inverse distance with nearest neighbor interpolation
        interpolation_method_string='invdistnn:radius='+str(radius1)+':max_points=100:nodata=-9999'

    if interpolation_method == 'lin':
        #linear interpolation:
        interpolation_method_string='linear:radius='+str(radius1)+':nodata=-9999'

    epsg_string='epsg:'+str(inps.epsg_code)
    if grid_datatype == 'Int16' or grid_datatype == 'Int8':
        predictor=2
    else:
        predictor=3

    if cliplayer != '':
        #consider adding option for multi-core processing: '--config', 'GDAL_NUM_THREADS ALL_CPUS', 
        cmd = ['gdal_grid', '-of', 'GTiff', '-co', 'PREDICTOR=%s'%(str(predictor)), \
               '-co', 'COMPRESS=DEFLATE', '-co', 'ZLEVEL=7', '-ot', grid_datatype, \
               '-txe', str(np.min(x_coords)), str(np.max(x_coords)), '-tye', str(np.max(y_coords)), str(np.min(y_coords)),  \
               '-zfield', zfield, '-a', interpolation_method_string, \
               '-outsize', str(ncols), str(nrows), '-a_srs', epsg_string, \
               input_vrt, output_grid2]
        cmd2 = ['gdalwarp', '-tr', str(resolution_m), str(resolution_m), \
                '-co', 'PREDICTOR=%s'%(str(predictor)), '-co', 'COMPRESS=DEFLATE', '-co', 'ZLEVEL=7', '-ot', grid_datatype, \
                '-srcnodata', '-9999', '-dstnodata', '-9999', \
                '-r', 'bilinear', '-cutline', cliplayer, \
                output_grid2, output_grid]                     
    else:
        cmd = ['gdal_grid', '-of', 'GTiff', '-co', 'PREDICTOR=%s'%(str(predictor)), \
               '-co', 'COMPRESS=DEFLATE', '-co', 'ZLEVEL=7', '-ot', grid_datatype, \
               '-txe', str(np.min(x_coords)), str(np.max(x_coords)), '-tye', str(np.max(y_coords)), str(np.min(y_coords)),  \
               '-zfield', zfield, '-a', interpolation_method_string, \
               '-outsize',str(ncols), str(nrows), '-a_srs', epsg_string,  \
               input_vrt, output_grid]
#    print('\n')
#    print(' '.join(cmd))
#    print('\n')
#    print(' '.join(cmd2))
    logfile_fname = os.path.join(inps.basedir, 'log') + '/gdal_grid_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '.txt'
    logfile_error_fname = os.path.join(inps.basedir, 'log') + '/gdal_grid_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '_err.txt'
    with open(logfile_fname, 'w') as out, open(logfile_error_fname, 'w') as err:
        subprocess_p = subprocess.Popen(cmd, stdout=out, stderr=err)
        subprocess_p.wait()
    logfile_fname = os.path.join(inps.basedir, 'log') + '/gdalwarp_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '.txt'
    logfile_error_fname = os.path.join(inps.basedir, 'log') + '/gdalwarp_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '_err.txt'
    with open(logfile_fname, 'w') as out, open(logfile_error_fname, 'w') as err:
        subprocess_p = subprocess.Popen(cmd2, stdout=out, stderr=err)
        subprocess_p.wait()
    if os.path.exists(output_grid2):
        os.remove(output_grid2)
    if os.path.exists(output_grid):
        ds = gdal.Open(output_grid)
        data_tif = np.array(ds.GetRasterBand(1).ReadAsArray()).astype(float)
        data_tif[np.where(data_tif == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
        ds = None
    else:
        print('gdal_grid and gdalwarp could not create %s.'%output_grid)
        data_tif=np.NaN
    return data_tif

def calc_stats_for_seed_points_wrapper(j):
    #print('starting {}/{}'.format(j+1, len(pos_array)))
    from_pos = pos_array[j] #Get start/end from position array
    to_pos = pos_array[j+1]
    #Setup array for seed point results:
    subarr = np.arange(from_pos,to_pos) #Slice the data into the selected part...
    pts_seed_stats_result = np.empty((subarr.shape[0], nr_of_datasets))
    
    #Setup array for PC results (X, Y, Z, Dz)
    dxyzn_subarr_result = np.empty((subarr.shape[0], dxyzn_max_nre, 4))
 
    for ii in range(subarr.shape[0]):
        pts_seed_stats_result[ii,:], dxyzn_subarr_result[ii,:,:] = calc_stats_for_seed_points(subarr[ii]) #Run point cloud processing for this inddex

    pickle_fn = os.path.join(pickle_dir, 'PC_seed_points_{}.pickle'.format(str(j).zfill(4)))
    pickle.dump((pts_seed_stats_result, dxyzn_subarr_result), open(pickle_fn,'wb'))
    if np.mod(j,10) == 0:
        print('{}, '.format(str(j).zfill(2)), end='', flush=True)
    pts_seed_stats_result = None
    dxyzn_subarr_result = None
        
def calc_stats_for_seed_points(k):
    ids2use = np.where(pc_xyz_distance[k] != np.inf)[0]
    pts_xyz = pc_xyz[pc_xyz_distance_id[k][ids2use]]
    #Make sure there are no points selected twice (remove them)
    pts_xyz = np.unique(pts_xyz,axis=0)
        
    nr_pts_xyz = pts_xyz.shape[0]
    if inps.mean_z_only == 1:
        #only calculate mean elevation for grid-cell size
        if pts_xyz.shape[0] < inps.pt_lower_threshold:
            pts_xyz_meanpt = np.nan
            #pts_xyz_normal = np.nan
            pts_seed_stats = np.array([pc_xyz_rstep_seed[k,0], pc_xyz_rstep_seed[k,1], pc_xyz_rstep_seed[k,2], 
                       np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 
                       np.nan, np.nan, np.nan, np.nan, np.nan, nr_pts_xyz, np.nan, np.nan, np.nan])
            dxyzn = np.empty((dxyzn_max_nre, 4))
            dxyzn.fill(np.nan)
        else:
            pts_seed_stats = np.array([pc_xyz_rstep_seed[k,0], pc_xyz_rstep_seed[k,1], pc_xyz_rstep_seed[k,2], 
                           np.mean(pts_xyz[:,0]), np.mean(pts_xyz[:,1]), np.mean(pts_xyz[:,2]), 
                           np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 
                           np.nan, np.nan, np.nan, np.nan, np.nan, nr_pts_xyz, np.nan, np.nan, np.nan])
            dxyzn = np.empty((dxyzn_max_nre, 4))
            dxyzn.fill(np.nan)
        return pts_seed_stats, dxyzn
        
    if pts_xyz.shape[0] < inps.pt_lower_threshold:
        #print('Less than %d points, plane fitting not meaningful for k = %s'%(inps.pt_lower_threshold,"{:,}".format(k)))
        pts_xyz_meanpt = np.nan
        #pts_xyz_normal = np.nan
        pts_seed_stats = np.array([pc_xyz_rstep_seed[k,0], pc_xyz_rstep_seed[k,1], pc_xyz_rstep_seed[k,2], 
                   np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 
                   np.nan, np.nan, np.nan, np.nan, np.nan, nr_pts_xyz, np.nan, np.nan, np.nan])
        dxyzn = np.empty((dxyzn_max_nre, 4))
        dxyzn.fill(np.nan)
    else:
        #not using planeFit, because it is very sensitive to outliers - only works well on
        #well-behaved terrain
        #pts_xyz_meanpt, pts_xyz_normal, plane_residual = planeFit(pts_xyz.T)

        #Instead use least square approach of fitting plane
        p1_slope_lstsq,_,_,_,p1_slope_lstsq_residual,p1_dz,C_p1 = curvFit_lstsq_polygon(pts_xyz.T, order=1)

        #calculate curvature
        p2_slope_lstsq, curvature_lstsq, _, _, curv_residuals,p2_dz,C_p2 = curvFit_lstsq_polygon(pts_xyz.T, order=2)

        pts_xyz_meanpt = np.mean(pts_xyz, axis=0)

        #normalize /detrend points with planeFit results - DO NOT USE FOR LIDAR PC
        #d = -pts_xyz_meanpt.dot(pts_xyz_normal)
        #z = (-pts_xyz_normal[0] * pts_xyz[:,0] - pts_xyz_normal[1] * pts_xyz[:,1] - d) * 1. /pts_xyz_normal[2]
        #plane_slope = pts_xyz_normal[2]
        #calculate offset for each point from plane
        #dz = pts_xyz[:,2] - z
    
        #Plotting for testing purposes
        if inps.plot_plane_fits == 1 and nr_pts_xyz >= 10:
            plot_seed_pts_neighborhood_3d_fn=os.path.join(map_dir,'PlaneFit_seed%08d.png'%k)
            plot_seed_pts_neighborhood_3d(pts_xyz, C_p1, C_p2, p1_dz, p2_dz, plot_seed_pts_neighborhood_3d_fn)
        
        #Using least squared linear (p1) fitting results for normalizing plane
        dz=p1_dz

        #stack points into X, Y, Z, delta-Z for each point
        dxyzn = np.empty((dxyzn_max_nre, 4))
        dxyzn.fill(np.nan)
        dxyzn[range(pts_xyz.shape[0]),:] = np.vstack([np.vstack((pts_xyz[:,0], pts_xyz[:,1], pts_xyz[:,2], dz)).T])
    
        #for each seed point, store relevant point statistics. Columns are:
        #0: Seed-X, 1: Seed-Y, 2: Seed-Z, 3: Mean-X, 4: Mean-Y, 
        #5: Mean-Z, 6: Z-min, 7: Z-max, 8: Dz-max, 9: Dz-min, 
        #10: Dz-std.dev, 11: Dz-range, 12: Dz-90-10th percentile range, 13: Dz-75-25th percentile range, 14: variance dz, 
        #15: Slp linear, 16: Slp-linear residual, 17: Slp_p2, 18: Slp_p2 residual, 19: nr. of lidar points, 
        #20: curvature_lstsq, 21: std. dev. of Z
        pts_seed_stats = np.array([pc_xyz_rstep_seed[k,0], pc_xyz_rstep_seed[k,1], pc_xyz_rstep_seed[k,2], 
                       pts_xyz_meanpt[0], pts_xyz_meanpt[1], pts_xyz_meanpt[2], 
                       np.min(pts_xyz, axis=0)[2], np.max(pts_xyz, axis=0)[2], dz.max(), dz.min(), \
                       np.std(dz), dz.max()-dz.min(), np.percentile(dz, 90)-np.percentile(dz,10), \
                       np.percentile(dz, 75)-np.percentile(dz,25), np.var(dz), \
                       p1_slope_lstsq, p1_slope_lstsq_residual, p2_slope_lstsq, curv_residuals, \
                       nr_pts_xyz, curvature_lstsq, np.nanstd(pts_xyz[:,2])])
    return pts_seed_stats, dxyzn


def plot_seed_pts_neighborhood_3d(pts_xyz, C_p1, C_p2, p1_dz, p2_dz, plot_seed_pts_neighborhood_3d_fn):
#For testing purposes: plot pointcloud
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.spatial.distance import pdist, squareform

# Alternative approach to plane fitting - but least squared approach is more robust
#    pts_x = pts_xyz[:,0]
#    pts_y = pts_xyz[:,1]
#    pts_z = pts_xyz[:,2]
#    tmp_A = []
#    tmp_b = []
#    for i in range(len(pts_x)):
#        tmp_A.append([pts_x[i], pts_y[i], 1])
#        tmp_b.append(pts_z[i])
#    b = np.matrix(tmp_b).T
#    A = np.matrix(tmp_A)
#    fit = (A.T * A).I * A.T * b
#    errors = b - A * fit
#    residual = np.linalg.norm(errors)
#    print("solution: %f x + %f y + %f = z" % (fit[0], fit[1], fit[2]))
#
#    d = -pts_xyz_meanpt.dot(pts_xyz_normal)
#    z = (-pts_xyz_normal[0] * pts_xyz[:,0] - pts_xyz_normal[1] * pts_xyz[:,1] - d) * 1. /pts_xyz_normal[2]
#    plane_slope = pts_xyz_normal[2]

    pts_xyz_meanpt = np.mean(pts_xyz, axis=0)
    # calculate plane parameters
    xlim = np.arange(pts_xyz[:,0].min(), pts_xyz[:,0].max(), 0.05)
    ylim = np.arange(pts_xyz[:,1].min(), pts_xyz[:,1].max(), 0.05)
    X,Y = np.meshgrid(xlim, ylim)
    Z_p1 = C_p1[0]*X + C_p1[1]*Y + C_p1[2]
    Z_p2 = C_p2[0]*X**2. + C_p2[1]*Y**2. + C_p2[2]*X*Y + C_p2[3]*X + C_p2[4]*Y + C_p2[5]

    fig = plt.figure(figsize=(16.53*1.5,11.69*1.5), dpi=150)
    ax = fig.add_subplot(121, projection='3d')
    ax.scatter(pts_xyz[:,0], pts_xyz[:,1], pts_xyz[:,2], c='b', marker='.', s=72*2)
    ax.scatter(pts_xyz_meanpt[0], pts_xyz_meanpt[1], pts_xyz_meanpt[2], c='r', marker='x',s=96*2)
    ax.plot_wireframe(X,Y,Z_p1, color='k')
    ax.plot_wireframe(X,Y,Z_p2, color='b')
    ax.set_xlabel('UTM-X (m)', fontsize=18)
    ax.set_ylabel('UTM-Y (m)', fontsize=18)
    ax.set_zlabel('UTM-Z (m)', fontsize=18)
    ax.set_title('Points and plane (Z_p1: k, Z_p2:blue)', fontsize=24)    

    ax2 = fig.add_subplot(122)
    A = np.array([pts_xyz[:,0], pts_xyz[:,1]]).T
    # explicitly calculate the whole n x n distance matrix
    dist_mat = squareform(pdist(A, metric="euclidean"))
    # mask the diagonal
    np.fill_diagonal(dist_mat, np.nan)
    # and calculate the minimum of each row (or column)
    distances = np.nanmin(dist_mat, axis=1)
    ax2.scatter(distances, p1_dz, c='k', marker='.', s=72*2)
    ax2.scatter(distances, p2_dz, c='b', marker='+', s=72*2)
    ax2.set_xlabel('Eucledian Distance between X-Y points (m)', fontsize=18)
    ax2.set_ylabel('normalized height by plane, Dz (m)', fontsize=18)
    ax2.set_title('Normalized Heights', fontsize=24)    
    ax2.axis('equal')
    ax2.grid()
    fig.savefig(plot_seed_pts_neighborhood_3d_fn, bbox_inches='tight')
    plt.close()


def pc_seed_random_k_subsampling(pc_xyz, pc_xyz_rstep_seed, pc_xyz_seed_pts_id, k):
    '''
    Sub-samples pointcloud (pc_xyz) for seed point list (pc_xyz_rstep_seed) and neighbors (pc_xyz_seed_pts_id).
    Subsampling is performed for k random neighbors. If k=nr. of points in neighborhood, all points in that
    neighborhood are chosen.
    
    pc_xyzg_k_random = pc_seed_random_k_subsampling(pc_xyz, pc_xyz_rstep_seed, pc_xyz_seed_pts_id, k)
    '''
    
    #iterate through n number of points (length of seed points)
    n = pc_xyz_rstep_seed.shape[0]
    pc_xyzg_k_random = np.empty((n*k,3))
    pc_xyzg_k_random.fill(np.nan)
    counter = 0
    for i in range(n):
        current_indices = np.array(pc_xyz_seed_pts_id[i][pc_xyz_seed_pts_id[i].mask==False])
        random_numbers = np.random.randint(0,len(current_indices),size=(k))
        if len(current_indices) >= k:
            pc_xyzg_k_random[counter:counter+k,0] = pc_xyz[ current_indices[ random_numbers.astype('int8') ], 0]
            pc_xyzg_k_random[counter:counter+k,1] = pc_xyz[ current_indices[ random_numbers.astype('int8') ], 1]
            pc_xyzg_k_random[counter:counter+k,2] = pc_xyz[ current_indices[ random_numbers.astype('int8') ], 2]
            counter = counter + k
        elif len(current_indices) < k:
            #print('\tpc_seed_random_k_subsampling: Did not find k=%d number of neighbors (found=%d).\n\tReduce nr. of neighborhood points or increase radius.'%(k, len(current_indices)))
            pc_xyzg_k_random[counter:counter+len(current_indices),0] = pc_xyz[ current_indices, 0]
            pc_xyzg_k_random[counter:counter+len(current_indices),1] = pc_xyz[ current_indices, 1]
            pc_xyzg_k_random[counter:counter+len(current_indices),2] = pc_xyz[ current_indices, 2]
            counter = counter + len(current_indices)
        else:
            pc_xyzg_k_random[counter:counter+len(current_indices),0] = pc_xyz[ current_indices, 0]
            pc_xyzg_k_random[counter:counter+len(current_indices),1] = pc_xyz[ current_indices, 1]
            pc_xyzg_k_random[counter:counter+len(current_indices),2] = pc_xyz[ current_indices, 2]
            counter = counter + len(current_indices)
        
    pc_xyzg_k_random = pc_xyzg_k_random[~np.isnan(pc_xyzg_k_random).any(axis=1)]
    return pc_xyzg_k_random


def pc_random_p_subsampling(pc_xyz, p, nr_of_out_points):
    '''
    Sub-samples indices of PC pc_xyzg_radius with probability weight p based 
    on point density of each point. Will result in greatly reduced point cloud. 
    Give nr_of_out_points for subsampled point cloud, usually len(p)/2
    
    call with a probability
    #pc_xyzg_radius_equal_nr_random = pc_random_p_subsampling(pc_xyzg_radius, pts_xyz, nr_of_out_points)
    '''
    
    #iterate through n number of points (length of seed  points)
    n = len(p)
    pc_xyz_p_random = np.empty((int(nr_of_out_points),3))
    pc_xyz_p_random.fill(np.nan)
    i = np.random.choice(n, size = int(nr_of_out_points), replace = False, p = p)
    pc_xyz_p_random[:,0] = pc_xyz[i,0]
    pc_xyz_p_random[:,1] = pc_xyz[i,1]
    pc_xyz_p_random[:,2] = pc_xyz[i,2]
    return pc_xyz_p_random

def pc_DynamicDensity(pc_xyz_kdtree, pts_seed_xyz, k = 10):
    radi, _ = pc_xyz_kdtree.query(pts_seed_xyz, k = k)
    dens =  k / np.pi / radi[:, -1]**2
    dens.shape = pts_seed_xyz.shape[0]
    disk = np.pi * radi[:, -1]**2
    probability = dens / disk    
    return dens, probability

#Start of the main program
if __name__ == '__main__': 

    inps = cmdLineParser()

## Testing
#    inps = argparse.ArgumentParser(description='PointCloud (PC) processing for DEM statistics. Deriving gridded ground data (elevation and slope) using centroid coordinates. B. Bookhagen (bodo.bookhagen@uni-potsdam.de), V0.1 Oct 2018.')
#    inps.nr_of_cores = 0
#    inps.inlas = 'Pozo_USGS_UTM11_NAD83_cat16_SMRF_cl2.las'
#    inps.raster_m_range='1 10 1'
#    inps.shapefile_clip = 'Pozo_DTM_noveg_UTM11_NAD83_cat16.shp'
#    inps.epsg_code=26911
#    inps.max_nr_of_neighbors_kdtree = 100
#    inps.pt_lower_threshold = 3
#    inps.create_geotiff = 1
#    inps.create_las = 1
#    inps.create_gmt = 0
#    inps.create_shapefiles = 0
#    inps.mean_z_only=0
#    inps.subsample_1m_pc_k = 10
#    inps.subsample_1m_pc_p = 0
#    inps.k_nr_of_neighbors = 10
#    inps.plot_plane_fits = 0
    

    if inps.subsample_1m_pc_k > 0 and inps.subsample_1m_pc_p > 0:
        print('subsample_1m_pc_k and subsample_1m_pc_p has been chosen - only one is possible. Exiting.')
        exit()
        
    if not inps.raster_m_range == '':
        #not empty, extracting start, end, and step sizes
        range_m_start = float(inps.raster_m_range.split(' ')[0])
        range_m_stop = float(inps.raster_m_range.split(' ')[1])
        range_m_step = float(inps.raster_m_range.split(' ')[2])
        
    inps.basedir = os.path.dirname(inps.inlas)
    if inps.basedir == '':
        inps.basedir = os.getcwd()
        
    pickle_dir = os.path.join(inps.basedir, 'pickle')
    if os.path.exists(pickle_dir) == False:
        os.mkdir(pickle_dir)

    map_dir = os.path.join(inps.basedir, 'maps')
    if os.path.exists(map_dir) == False:
        os.mkdir(map_dir)
    
    geotif_dir = os.path.join(inps.basedir, 'geotiff')
    if os.path.exists(geotif_dir) == False:
        os.mkdir(geotif_dir)
    
    vrt_dir = os.path.join(inps.basedir, 'vrt')
    if os.path.exists(vrt_dir) == False:
        os.mkdir(vrt_dir)

    hdf_dir = os.path.join(inps.basedir, 'hdf')
    if os.path.exists(hdf_dir) == False:
        os.mkdir(hdf_dir)

    las_dir = os.path.join(inps.basedir, 'LAS')
    if os.path.exists(las_dir) == False:
        os.mkdir(las_dir)

    log_dir = os.path.join(inps.basedir, 'log')
    if os.path.exists(log_dir) == False:
        os.mkdir(log_dir)
       
    ### Loading data and filtering
    if os.path.exists(inps.inlas) == False:
        print('\n%s does not exist. Exiting.'%inps.inlas)
        sys.exit
    else:    
        print('\nLoading input file: %s'%inps.inlas)
    inFile = File(inps.inlas, mode='r')
    pc_xyz = np.vstack((inFile.get_x()*inFile.header.scale[0]+inFile.header.offset[0], inFile.get_y()*inFile.header.scale[1]+inFile.header.offset[1], inFile.get_z()*inFile.header.scale[2]+inFile.header.offset[2])).transpose()
    pc_xy = np.vstack((inFile.get_x()*inFile.header.scale[0]+inFile.header.offset[0], inFile.get_y()*inFile.header.scale[1]+inFile.header.offset[1])).transpose()
    #pc_xyz is now a point cloud with x, y, z
    print('\tLoaded %s points'%"{:,}".format(pc_xyz.shape[0]))
    
    ### Random subsampling (if applied)
    if inps.subsample_1m_pc_k > 0:
        #randomly subsampled point cloud for 0.5m radius
        print('\tRandomly subsample pointcloud to k=%d points for r=0.5m or step size=1m ... '%(inps.subsample_1m_pc_k), end='', flush=True)
        output_las_fn = inps.inlas[:-4] + '_randomsubsampled_k%02d.las'%(inps.subsample_1m_pc_k)
        if os.path.exists(output_las_fn):
            print('\n\tLoading existing point cloud %s ... '%(output_las_fn), end='', flush=True)
            inFile = File(output_las_fn, mode='r')
            pc_xyz = np.vstack((inFile.get_x()*inFile.header.scale[0]+inFile.header.offset[0], inFile.get_y()*inFile.header.scale[1]+inFile.header.offset[1], inFile.get_z()*inFile.header.scale[2]+inFile.header.offset[2])).transpose()
            pc_xy = np.vstack((inFile.get_x()*inFile.header.scale[0]+inFile.header.offset[0], inFile.get_y()*inFile.header.scale[1]+inFile.header.offset[1])).transpose()
            print('done. Loaded %s points. '%"{:,}".format(pc_xyz.shape[0]), flush=True)
        else:
            [x_min, x_max] = np.min(pc_xyz[:,0]), np.max(pc_xyz[:,0])
            [y_min, y_max] = np.min(pc_xyz[:,1]), np.max(pc_xyz[:,1])
            [z_min, z_max] = np.min(pc_xyz[:,2]), np.max(pc_xyz[:,2])
            current_rstep_size = 1
            current_search_radius = current_rstep_size/2
            rstep_size = current_rstep_size
            x_elements = len(np.arange(x_min.round(), x_max.round(), rstep_size))
            y_elements = len(np.arange(y_min.round(), y_max.round(), rstep_size))        
            x_coords = np.arange(x_min.round(), x_max.round(), rstep_size) + rstep_size / 2
            y_coords = np.arange(y_min.round(), y_max.round(), rstep_size) + rstep_size / 2        
            xy_coordinates = np.array([(x,y) for x in x_coords for y in y_coords])            
            pc_xy_pykdtree = KDTree(pc_xy)
            [pc_xy_pykdtree_distance, pc_xy_pykdtree_id] = pc_xy_pykdtree.query(xy_coordinates, k=1)        
            pc_distances_lt_rstep_size = np.where(pc_xy_pykdtree_distance <= rstep_size)[0]        
            pc_xy_pykdtree_id = pc_xy_pykdtree_id[pc_distances_lt_rstep_size]                
            pc_xyz_rstep_seed = pc_xyz[pc_xy_pykdtree_id]        
    
            pc_xyz_kdtree = KDTree(pc_xyz)
            pc_xyz_distance, pc_xyz_seed_pts_id = pc_xyz_kdtree.query(pc_xyz_rstep_seed, 
                                                                        k=inps.max_nr_of_neighbors_kdtree, 
                                                                        distance_upper_bound=current_search_radius)
            pc_xyz_distance = pc_xyz_distance.astype('float32')
            pc_xyz_data_indices = [np.where(x==True)[0] for x in pc_xyz_distance != np.inf]       
            dxyzn_max_nre = np.max([x.shape[0] for x in pc_xyz_data_indices])
            pc_xyz_distance = pc_xyz_distance[:,0:dxyzn_max_nre]
            pc_xyz_seed_pts_id = pc_xyz_seed_pts_id[:,0:dxyzn_max_nre]
            #convert to masked array for         
            pc_xyz_seed_pts_id = np.ma.masked_equal(pc_xyz_seed_pts_id, pc_xyz_seed_pts_id[0,dxyzn_max_nre-1], copy=False)
            
            #pc_xyz_seed_pts_id = np.array([pc_xyz_distance_id[x] for x in pc_xyz_data_indices][0])
            #now call random subsampling. Pass seed point list and ids of each seed point neighbor:
            pc_xyz = pc_seed_random_k_subsampling(pc_xyz, pc_xyz_rstep_seed, 
                                                  pc_xyz_seed_pts_id, inps.subsample_1m_pc_k)
            pc_xy = pc_xyz[:,0:2]
            print('done.')
            output_las_fn = inps.inlas[:-4] + '_randomsubsampled_k%02d.las'%(inps.subsample_1m_pc_k)
            print('\tWriting subsampled point cloud to %s ... '%(output_las_fn), end='', flush=True)
            #normalize input and generate colors for height using colormap
            v = pc_xyz[:,2]
            #stretch to 10-90th percentile
            v_1090p = np.percentile(v, [10, 90])
            v_rescale = exposure.rescale_intensity(v, in_range=(v_1090p[0], v_1090p[1]))
            colormap_PuOr = mpl.cm.PuOr
            rgb = colormap_PuOr(v_rescale)
            #remove last column - alpha value
            rgb = (rgb[:, :3] * (np.power(2,16)-1)).astype('uint16')    
            outFile = File(output_las_fn, mode='w', header=inFile.header)
            new_header = copy.copy(outFile.header)
            #setting some variables
            new_header.created_year = datetime.datetime.now().year
            new_header.created_day = datetime.datetime.now().timetuple().tm_yday
            new_header.x_max = pc_xyz[:,0].max()
            new_header.x_min = pc_xyz[:,0].min()
            new_header.y_max = pc_xyz[:,1].max()
            new_header.y_min = pc_xyz[:,1].min()
            new_header.z_max = pc_xyz[:,2].max()
            new_header.z_min = pc_xyz[:,2].min()
            new_header.point_records_count = pc_xyz.shape[0]
            new_header.point_return_count = 0
            outFile.header.count = v.shape[0]
            new_header.scale=inFile.header.scale
            new_header.offset=inFile.header.offset
            outFile.X = (pc_xyz[:,0]-inFile.header.offset[0])/inFile.header.scale[0]
            outFile.Y = (pc_xyz[:,1]-inFile.header.offset[1])/inFile.header.scale[1]
            outFile.Z = (pc_xyz[:,2]-inFile.header.offset[2])/inFile.header.scale[2]
            outFile.Red = rgb[:,0]
            outFile.Green = rgb[:,1]
            outFile.Blue = rgb[:,2]    
            outFile.close()    
            print('done.')

    if inps.subsample_1m_pc_p > 0:
        #randomly subsampled point cloud for 0.5m radius with probability p
        nr_of_output_points = int(pc_xyz.shape[0] * inps.subsample_1m_pc_p)
        print('\tHomogeneous subsampling of pointcloud with probability to create an equal point density PC from %s points to %s points (factor: %s) ... '%("{:,}".format(pc_xyz.shape[0]), "{:,}".format(int(nr_of_output_points)), "{:,}".format(inps.subsample_1m_pc_p)), end='', flush=True )
        output_las_fn = inps.inlas[:-4] + '_psubsampled_n%d.las'%(nr_of_output_points)
        if os.path.exists(output_las_fn) and inps.redo_subsample_1m_pc_p == False:
            print('Loading existing point cloud %s ... '%(output_las_fn), end='', flush=True)
            inFile = File(output_las_fn, mode='r')
            pc_xyz = np.vstack((inFile.get_x()*inFile.header.scale[0]+inFile.header.offset[0], inFile.get_y()*inFile.header.scale[1]+inFile.header.offset[1], inFile.get_z()*inFile.header.scale[2]+inFile.header.offset[2])).transpose()
            pc_xy = np.vstack((inFile.get_x()*inFile.header.scale[0]+inFile.header.offset[0], inFile.get_y()*inFile.header.scale[1]+inFile.header.offset[1])).transpose()
            print('done.')
        else:
            pc_xyz_kdtree = KDTree(pc_xyz)
            pc_xyz_density_k, pc_xyz_density_k_p = pc_DynamicDensity(pc_xyz_kdtree, pc_xyz, k = inps.k_nr_of_neighbors)        
            p_norm = pc_xyz_density_k_p / pc_xyz_density_k_p.sum()
            pc_xyz = pc_random_p_subsampling(pc_xyz, p_norm, nr_of_output_points)
            pc_xy = pc_xyz[:,0:2]
            print('done.')
            output_las_fn = inps.inlas[:-4] + '_psubsampled_n%d.las'%(nr_of_output_points)
            print('\tWriting homogeneously subsampled point cloud to %s ... '%(output_las_fn), end='', flush=True)
            #normalize input and generate colors for height using colormap
            v = pc_xyz[:,2]
            #stretch to 10-90th percentile
            v_1090p = np.percentile(v, [10, 90])
            v_rescale = exposure.rescale_intensity(v, in_range=(v_1090p[0], v_1090p[1]))
            colormap_PuOr = mpl.cm.PuOr
            rgb = colormap_PuOr(v_rescale)
            #remove last column - alpha value
            rgb = (rgb[:, :3] * (np.power(2,16)-1)).astype('uint16')    
            outFile = File(output_las_fn, mode='w', header=inFile.header)
            new_header = copy.copy(outFile.header)
            #setting some variables
            new_header.created_year = datetime.datetime.now().year
            new_header.created_day = datetime.datetime.now().timetuple().tm_yday
            new_header.x_max = pc_xyz[:,0].max()
            new_header.x_min = pc_xyz[:,0].min()
            new_header.y_max = pc_xyz[:,1].max()
            new_header.y_min = pc_xyz[:,1].min()
            new_header.z_max = pc_xyz[:,2].max()
            new_header.z_min = pc_xyz[:,2].min()
            new_header.point_records_count = pc_xyz.shape[0]
            new_header.point_return_count = 0
            outFile.header.count = v.shape[0]
            outFile.X = pc_xyz[:,0]
            outFile.Y = pc_xyz[:,1]
            outFile.Z = pc_xyz[:,2]
            outFile.Red = rgb[:,0]
            outFile.Green = rgb[:,1]
            outFile.Blue = rgb[:,2]    
            outFile.close()    
            print('done.')
        
    ### pyKDTree setup and calculation
    #Generate KDTree for fast searching
    #cKDTree is faster than KDTree, but pyKDTree is fast then cKDTree
    print('\tGenerating XY-pyKDTree (2D) ... ',end='', flush=True)
    pc_xy_pykdtree = KDTree(pc_xy)
    print('done.')
    
    print('\tGenerating XYZ-pyKDTree (3D) ... ',end='', flush=True)
    pc_xyz_kdtree = KDTree(pc_xyz)
    print('done.')
    print('')

    #Setup coordinates
    [x_min, x_max] = np.min(pc_xyz[:,0]), np.max(pc_xyz[:,0])
    [y_min, y_max] = np.min(pc_xyz[:,1]), np.max(pc_xyz[:,1])
    [z_min, z_max] = np.min(pc_xyz[:,2]), np.max(pc_xyz[:,2])
    
    ### Search KDTree with points on a regularly-spaced raster
    #generating equally-spaced raster overlay from input coordinates with stepsize rstep_size
    #This will be used to query the point cloud. Step_size should be small enough and likely 1/2 of the output file resolution. 
    #Note that this uses a 2D raster overlay to slice a 3D point cloud.
    range_m_list = np.arange(range_m_start, range_m_stop+range_m_step, range_m_step)
    ts0 = time.time()
    for i in range(len(range_m_list)):
        current_rstep_size = range_m_list[i]
        print('At iteration %d of %d with grid size %0.2f m\n'%(i+1, len(range_m_list), current_rstep_size), end='', flush=True)
        current_search_radius = current_rstep_size/2
        rstep_size = current_rstep_size
        x_elements = len(np.arange(x_min.round(), x_max.round(), rstep_size))
        y_elements = len(np.arange(y_min.round(), y_max.round(), rstep_size))
        
        #get coordinate range and shift coordinates by half of the step size to make sure rater overlay is centered. 
        #This is not really necessary and only matters for very small point clouds with edge effects or for very large steps sizes:
        x_coords = np.arange(x_min.round(), x_max.round(), rstep_size) + rstep_size / 2
        y_coords = np.arange(y_min.round(), y_max.round(), rstep_size) + rstep_size / 2
        
        #create combination of all coordinates (this is using lists and could be optimized)
        xy_coordinates = np.array([(x,y) for x in x_coords for y in y_coords])

        #test if data were already processed and saved as H5:
        hdf5_fn = str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_seed_pts_stats_raster_%0.2fm.h5'%current_rstep_size
        pc_results_fn = os.path.join(hdf_dir, hdf5_fn)
        if os.path.exists(pc_results_fn) == True:
            print('\tLoading data from existing file:  %s '%hdf5_fn)
            #Still need to add
            #hdf_in = h5py.File(pc_results_fn,'r')
            
        #using the 2D KDTree to find the points that are closest to the defined 2D raster overlay
        [pc_xy_pykdtree_distance, pc_xy_pykdtree_id] = pc_xy_pykdtree.query(xy_coordinates, k=1)
        
        #take only points that are within the search radius (may omit some points at the border regions
        pc_distances_lt_rstep_size = np.where(pc_xy_pykdtree_distance <= rstep_size)[0]
        
        #the following list contains all IDs to the actual lidar points that are closest to the center of raster overlay. 
        #We will use these points as seed points for determining slopes, ground elevation and perform plane-fitting
        pc_xy_pykdtree_id = pc_xy_pykdtree_id[pc_distances_lt_rstep_size]
                
        #now select these points from the 3D pointcloud with X, Y, Z coordinates.
        #We refer to these as seed points from the rstep part of this script
        pc_xyz_rstep_seed = pc_xyz[pc_xy_pykdtree_id]

        #remove 2D pyKDTree from memory
        pc_xy_pykdtree_distance = None
        pc_xy_pykdtree_id = None
        
        ### Query points - use KDTree
        #find points from 3D seed / query points  / raster overlay with radius = inps.sphere_radius_m
        print('\tQuerying pyKDTree with radius %0.2f m... '%(current_search_radius), end='', flush=True)
        pc_xyz_distance, pc_xyz_distance_id  = pc_xyz_kdtree.query(pc_xyz_rstep_seed, k=inps.max_nr_of_neighbors_kdtree, distance_upper_bound=current_search_radius)
        print('done.')
        
        # Two steps to reduce memory use of pc_xyz_distance
        #First: pc_xyz_distance is of type float64 which is important for numerical purposes and very small distances. 
        #For lidar point clouds, this precision is not needed. We convert to float32 to easen to memory load.
        pc_xyz_distance = pc_xyz_distance.astype('float32')

        #Second: Find max. number of neighbors for each point - to cut down array size:
        #Because every point has a different number of neighbors (they are all stored in pc_xyz_distance), we create an array with size
        #of the max. number of points that are not inf
        dxyzn_max_nre = np.max([np.where(x==True)[0].shape[0] for x in pc_xyz_distance!=np.inf])
        pc_xyz_distance = pc_xyz_distance[:,0:dxyzn_max_nre]
            
        ## Setup variables
        pc_xyz_distance_nr = len(pc_xyz_distance)
        nr_of_seed_points = len(pc_xyz_rstep_seed)
        nr_of_datasets = 22 #nr of columns to save
        nr_of_processes = 100 #splitting the for loop into 100 processes and dividing array into 100 steps in pos_array
        dxyzn_nre = np.sum([len(x) for x in pc_xyz_distance])
        print('\tNumber of total points in neighborhood array to be processed (dxyzn_nre): %s'%"{:,}".format(dxyzn_nre))
        dxyzn_nre_pos_array = np.array(np.linspace(0, dxyzn_nre, nr_of_processes), dtype=int)
        pos_array = np.array(np.linspace(0, nr_of_seed_points, nr_of_processes), dtype=int) #This creates a position array so you can select from:to in each loop
        if inps.nr_of_cores == 0:
            inps.nr_of_cores = multiprocessing.cpu_count()
    
        
        print('\tCalculating PC statistics for radius %0.2f m using %d/%d cores for %s seed points [%%]: '%(current_search_radius, np.round(inps.nr_of_cores).astype(int), multiprocessing.cpu_count(), "{:,}".format(nr_of_seed_points)), end='', flush=True )
        #for each seed point, store relevant point statistics. Columns are:
        #0: Seed-X, 1: Seed-Y, 2: Seed-Z, 3: Mean-X, 4: Mean-Y, 
        #5: Mean-Z, 6: Z-min, 7: Z-max, 8: Dz-max, 9: Dz-min, 
        #10: Dz-std.dev, 11: Dz-range, 12: Dz-90-10th percentile range, 13: Dz-75-25th percentile range, 14: variance dz, 
        #15: slope_p1_leastsq, 16: slope_p1_residuals, 17: nr. of lidar points, 18: slope_p2_lstsq, 19: curvature_lstsq, 
        #20: curvature residuals, 21: std. dev. of Z
    
        ts = time.time()
        p = Pool(processes=np.round(inps.nr_of_cores).astype(int))
        for _ in p.imap_unordered(calc_stats_for_seed_points_wrapper, np.arange(0,len(pos_array)-1)):
            pass    
        print('\n',end='', flush=True)
        
        #combine pickle files
        pkls = glob.glob(os.path.join(pickle_dir, 'PC_seed_points_*')) #Now get all the pickle files we made
        pkls.sort() #make sure they're sorted
        dxyzn = np.empty((dxyzn_nre, 4)) #output for every lidar point (dz value)
        pts_seed_stats = np.empty((pc_xyz_distance_nr,nr_of_datasets)) #output for seed points
        count = 0
        dxyzn_counter = 0
        for fid in pkls:
            seed_res, dxyzn_res = pickle.load(open(fid,'rb')) #Loop through and put each pickle into the right place in the output array
            if seed_res.shape[0] != pos_array[count+1] - pos_array[count]:
                print('File %s, length of records do not match. file: %d vs pos_array: %d'%(fid, seed_res.shape[0], pos_array[count+1] - pos_array[count]))
                if seed_res.shape[0] < pos_array[count+1] - pos_array[count]:
                    pts_seed_stats[range(pos_array[count],pos_array[count+1]-1),:] = seed_res
                elif seed_res.shape[0] > pos_array[count+1] - pos_array[count]:
                    pts_seed_stats[range(pos_array[count],pos_array[count+1]),:] = seed_res[:-1]
            else:
                pts_seed_stats[range(pos_array[count],pos_array[count+1]),:] = seed_res
            #re-arrange dxyzn and remove nans
            dxyzn_reshape = dxyzn_res.reshape((dxyzn_res.shape[0]*dxyzn_res.shape[1], dxyzn_res.shape[2]))
            idx_nonan = np.where(np.isnan(dxyzn_reshape[:,0]) == False)[0]
            dxyzn_reshape = dxyzn_reshape[idx_nonan,:]
            dxyzn[range(dxyzn_counter,dxyzn_counter+dxyzn_reshape.shape[0]),:] = dxyzn_reshape
            dxyzn_counter = dxyzn_counter + dxyzn_reshape.shape[0]
            count += 1
            del seed_res, dxyzn_res, dxyzn_reshape
        #remove pickle files
        for ii in pkls:
            os.remove(ii)
        pkls=None
    
        #write as compressed HDF file
        print('\tWriting to CSV, GMT, and shapefiles... ', end='', flush=True)
        #write csv
        header_str='1SeedX, 2SeedY, 3SeedZ, 4MeanX, 5MeanY, 6MeanZ, 7Z_min, 8Z_max, 9Dz_max, 10Dz_min, 11Dz_std, 12Dz_range, 13Dz_9010p, 14Dz_7525p, 15Dz_var, 16Slp_p1, 17Slp_p1r, 18Slp_p2, 19Slp_p2r, 20Nr_lidar, 21CurvLSQ, 22StdZ'
        seed_pts_stats_csv = '_seed_pts_stats_raster_%0.2fm.csv'%(current_rstep_size)
        pc_seed_pts_stats_csv_fn = os.path.join(vrt_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + seed_pts_stats_csv)
        seed_pts_stats_vrt = '_seed_pts_stats_raster_%0.2fm.vrt'%(current_rstep_size)
        pc_seed_pts_stats_vrt_fn = os.path.join(vrt_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + seed_pts_stats_vrt)
        seed_pts_stats_shp = '_seed_pts_stats_raster_%0.2fm.shp'%(current_rstep_size)
        pc_seed_pts_stats_shp_fn = os.path.join(vrt_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + seed_pts_stats_shp)
        seed_pts_stats_dd_shp = '_seed_pts_stats_raster_%0.2fm_dd.shp'%(current_rstep_size)
        pc_seed_pts_stats_dd_shp_fn = os.path.join(vrt_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + seed_pts_stats_dd_shp)
        seed_pts_stats_gmt = '_seed_pts_stats_raster_%0.2fm.gmt'%(current_rstep_size)
        pc_seed_pts_stats_gmt_fn = os.path.join(vrt_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + seed_pts_stats_gmt)
        seed_pts_stats_dd_gmt = '_seed_pts_stats_raster_%0.2fm_dd.gmt'%(current_rstep_size)
        pc_seed_pts_stats_dd_gmt_fn = os.path.join(vrt_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + seed_pts_stats_dd_gmt)
        idxnan = np.where(np.isnan(pts_seed_stats))
        if os.path.exists(pc_seed_pts_stats_csv_fn) == False:
            #before writing to CSV file, remove all lines with np.nan in pts_seed_stats
            pts_seed_stats_nonan = np.copy(pts_seed_stats)
            idx_nodata = np.where(np.isnan(pts_seed_stats_nonan))
            rows2remove = np.unique(idx_nodata[0])
            pts_seed_stats_nonan = np.delete(pts_seed_stats_nonan, (rows2remove), axis=0)
            np.savetxt(pc_seed_pts_stats_csv_fn, pts_seed_stats_nonan, fmt='%.4f', delimiter=',', header=header_str, comments='')
        pts_seed_stats_nonan = None
        idxnan = None
    
        # write VRT for shapefile generation
        vrt_f = open(pc_seed_pts_stats_vrt_fn,'w')
        vrt_f.write('<OGRVRTDataSource>\n')
        vrt_f.write('\t<OGRVRTLayer name="%s">\n'%os.path.basename(pc_seed_pts_stats_vrt_fn))
        vrt_f.write('\t\t<SrcDataSource>%s</SrcDataSource>\n'%os.path.join(vrt_dir, os.path.basename(pc_seed_pts_stats_csv_fn)))
        vrt_f.write('\t\t<SrcLayer>%s</SrcLayer>\n'%'.'.join(os.path.basename(pc_seed_pts_stats_csv_fn).split('.')[0:-1]))
        vrt_f.write('\t\t<LayerSRS>EPSG:%d</LayerSRS>\n'%inps.epsg_code)
        vrt_f.write('\t\t<GeometryType>wkbPoint</GeometryType>\n')
        vrt_f.write('\t\t<GeometryField encoding="PointFromColumns" x="4MeanX" y="5MeanY"/>\n')
        vrt_f.write('\t\t\t<Field name="1SeedX" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="2SeedY" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="3SeedZ" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="4MeanX" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="5MeanY" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="6MeanZ" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="7Z_min" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="8Z_max" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="9Dz_max" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="10Dz_min" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="11Dz_std" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="12Dz_range" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="13Dz_9010p" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="14Dz_7525p" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="15Dz_var" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="16Slp_p1" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="17Slp_p1r" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="18Slp_p2" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="19Slp_p2r" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="20Nr_lidar" type="Real" width="8"/>\n')
        vrt_f.write('\t\t\t<Field name="21CurvLSQ" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t\t\t<Field name="22StdZ" type="Real" width="8" precision="7"/>\n')
        vrt_f.write('\t</OGRVRTLayer>\n')
        vrt_f.write('</OGRVRTDataSource>\n')
        vrt_f.close()
    
        # Generate shapefile from VRT
        if inps.create_shapefiles == 1 and os.path.exists(pc_seed_pts_stats_shp_fn) == False:
            cwd=os.getcwd()
            os.chdir(vrt_dir)
            cmd = ['ogr2ogr', pc_seed_pts_stats_shp_fn, pc_seed_pts_stats_vrt_fn]
            logfile_fname = os.path.join(log_dir,  'ogr2ogr_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '.txt')
            logfile_error_fname = os.path.join(log_dir, 'ogr2ogr_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '_err.txt')
            with open(logfile_fname, 'w') as out, open(logfile_error_fname, 'w') as err:
                subprocess_p = subprocess.Popen(cmd, stdout=out, stderr=err)
                subprocess_p.wait()
            os.chdir(cwd)
        
        if inps.create_shapefiles == 1 and os.path.exists(pc_seed_pts_stats_dd_shp_fn) == False:
            cwd=os.getcwd()
            os.chdir(vrt_dir)
            cmd = ['ogr2ogr', '-t_srs', 'epsg:4326', pc_seed_pts_stats_dd_shp_fn, pc_seed_pts_stats_shp_fn]
            logfile_fname = os.path.join(log_dir, 'ogr2ogr_reproject_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '.txt')
            logfile_error_fname = os.path.join(log_dir, 'ogr2ogr_reproject_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '_err.txt')
            with open(logfile_fname, 'w') as out, open(logfile_error_fname, 'w') as err:
                subprocess_p = subprocess.Popen(cmd, stdout=out, stderr=err)
                subprocess_p.wait()
            os.chdir(cwd)
    
        # Convert VRT to GMT dataset for plotting
        if inps.create_gmt == 1 and os.path.exists(pc_seed_pts_stats_gmt_fn) == False:
            cwd=os.getcwd()
            os.chdir(vrt_dir)
            cmd = ['ogr2ogr', '-f', 'GMT', pc_seed_pts_stats_gmt_fn, pc_seed_pts_stats_vrt_fn]
            logfile_fname = os.path.join(log_dir, 'ogr2ogr_gmt_utm_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '.txt')
            logfile_error_fname = os.path.join(log_dir,  'ogr2ogr_gmt_utm_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '_err.txt')
            with open(logfile_fname, 'w') as out, open(logfile_error_fname, 'w') as err:
                subprocess_p = subprocess.Popen(cmd, stdout=out, stderr=err)
                subprocess_p.wait()
            os.chdir(cwd)

        if inps.create_gmt == 1 and os.path.exists(pc_seed_pts_stats_dd_gmt_fn) == False:
            cwd=os.getcwd()
            os.chdir(vrt_dir)
            cmd = ['ogr2ogr', '-f', 'GMT', '-t_srs', 'epsg:4326', pc_seed_pts_stats_dd_gmt_fn, pc_seed_pts_stats_gmt_fn]
            logfile_fname = os.path.join(log_dir,  'ogr2ogr_gmt_dd_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '.txt')
            logfile_error_fname = os.path.join(log_dir, 'ogr2ogr_gmt_dd_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '_err.txt')
            with open(logfile_fname, 'w') as out, open(logfile_error_fname, 'w') as err:
                subprocess_p = subprocess.Popen(cmd, stdout=out, stderr=err)
                subprocess_p.wait()
            os.chdir(cwd)
    
        print('done.')
    
        ### Interpolate to equally-spaced grid and generate GeoTIFF output
        if inps.create_geotiff == 1:
            print('\tInterpolating seed points (mean-X, mean-Y) and writing geotiff raster... ',end='', flush=True)
            idx_nonan = np.where(np.isnan(pts_seed_stats[:,3])==False)
        
            xres=current_rstep_size
            yres=current_rstep_size
            ncols=len(x_coords)
            nrows=len(y_coords)
            geotransform = (x_coords.min() - (current_rstep_size / 2), xres, 0 , y_coords.min() - (current_rstep_size / 2),0, yres) 

            srs = osr.SpatialReference()
            srs.ImportFromEPSG(inps.epsg_code)
            
            nrlidar_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_nrlidarpts.tif'%(current_rstep_size))
            z_std_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_z_stddev.tif'%(current_rstep_size))
            dz_std_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_dz_stddev.tif'%(current_rstep_size))
            dz_range9010_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_dz_range9010.tif'%(current_rstep_size))
            dz_iqr_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_dz_iqr.tif'%(current_rstep_size))
            p1_slope_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_slope_p1_lstsq.tif'%(current_rstep_size))
            p1_slope_res_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_slope_p1_lstsq_res.tif'%(current_rstep_size))
            p2_slope_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_slope_p2_lstsq.tif'%(current_rstep_size))
            p2_slope_res_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_slope_p2_lstsq_res.tif'%(current_rstep_size))
            mean_curv_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_mean_curv.tif'%(current_rstep_size))
            z_mean_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_z_mean.tif'%(current_rstep_size))

            #interpolate nr_lidar_measurements
            #We have experimented with the following options:
            ## First, use scipy.interpolate (see function griddata_clip_geotif). While this works, it's incredible slow
            #nr_lidar = griddata_clip_geotif(nrlidar_tif_fn, points, pts_seed_stats[idx_nonan,17][0], xxyy=(xx,yy), ncols=ncols, nrows=nrows, geotransform=geotransform)
            #
            ## Second, the API of gdal.Grid, but for some reason this wasn't as efficient (but should work well)
            #gridoptions = gdal.GridOptions(options=[], format="GTiff", outputType=gdal.GDT_Float32, width=nrows, height=ncols, creationOptions=None, \
            #                               outputBounds=[ np.min(x_coords), np.min(y_coords), np.max(x_coords) , np.max(y_coords)], \
            #                               outputSRS=srs, noData=-9999, algorithm='nearest:radius1=3.0:radius2=3.0:angle=0.0:nodata=-9999', \
            #                               layers='Pozo_USGS_UTM11_NAD83_all_color_cl_cat1_seed_pts_stats_raster_1.00m', SQLStatement=None, where=None, \
            #                               spatFilter=None, zfield='18Nr_lidar', z_increase=None, z_multiply=None, callback=None, callback_data=None)
            #ds = gdal.Grid(nrlidar_tif_fn, pc_seed_pts_stats_vrt_fn, options=gridoptions)
            #ds.SetGeoTransform(geotransform)
            #ds.SetProjection( srs.ExportToWkt() )
            #ds.GetRasterBand(1).WriteArray(datai) 
            #ds.FlushCache()
            #ds=None
            #
            ## Third, the command line version of gdal_grid. This appears to be the fastest option
            #                    nr_lidar=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
            #                                                   zfield='18Nr_Lidar', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=nrlidar_tif_fn,\
            #                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
            #
            #For determining the interpolation radius, we use twice the grid spacing
            radius1=current_rstep_size * 2
            radius2=radius1 * 2
            interpolation_method='lin'
            #interpolation_method='invdist'
            
            if os.path.exists(nrlidar_tif_fn) == False:
                print('nr_lidar_measurements, ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    nr_lidar=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='20Nr_Lidar', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=nrlidar_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Int16', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    nr_lidar=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='20Nr_Lidar', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=nrlidar_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Int16', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(nrlidar_tif_fn)
                nr_lidar = np.array(ds.GetRasterBand(1).ReadAsArray()).astype(float)
                nr_lidar[np.where(nr_lidar == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None
            
            if os.path.exists(z_std_tif_fn) == False:
                print('z std. dev., ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    z_std=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='22StdZ', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=z_std_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    z_std=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='22StdZ', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=z_std_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(z_std_tif_fn)
                z_std = np.array(ds.GetRasterBand(1).ReadAsArray())
                z_std[np.where(z_std == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None

            if os.path.exists(dz_std_tif_fn) == False:
                print('Dz std. dev., ', end='', flush=True)                
                if os.path.exists(inps.shapefile_clip):
                    dz_std=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='11Dz_std', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=dz_std_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    dz_std=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='11Dz_std', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=dz_std_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(dz_std_tif_fn)
                dz_std = np.array(ds.GetRasterBand(1).ReadAsArray())
                dz_std[np.where(dz_std == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None

            if os.path.exists(z_mean_tif_fn) == False:
                print('z mean, ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    z_mean=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='6MeanZ', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=z_mean_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    z_mean=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='6MeanZ', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=z_mean_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(z_mean_tif_fn)
                z_mean = np.array(ds.GetRasterBand(1).ReadAsArray())
                z_mean[np.where(z_mean == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None
            
            #interpolate Dz_range 90-10 percentile
            if os.path.exists(dz_range9010_tif_fn) == False:
                print('Dz range (90-10th perc.), ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    dz_range9010=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows, resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='13Dz_9010p', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=dz_range9010_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    dz_range9010=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='13Dz_9010p', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=dz_range9010_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(dz_range9010_tif_fn)
                dz_range9010 = np.array(ds.GetRasterBand(1).ReadAsArray())
                dz_range9010[np.where(dz_range9010 == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None
            
            #interpolate Dz_range 75-25 percentile
            if os.path.exists(dz_iqr_tif_fn) == False:
                print('Dz range (75-25th perc.), ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    dz_iqr=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='14Dz_7525p', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=dz_iqr_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    dz_iqr=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='14Dz_7525p', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=dz_iqr_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(dz_iqr_tif_fn)
                dz_iqr = np.array(ds.GetRasterBand(1).ReadAsArray())
                dz_iqr[np.where(dz_iqr == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None
                        
            #interpolate p1 slope
            if os.path.exists(p1_slope_tif_fn) == False:
                print('p1 slope, ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    p1_slope=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='16Slp_p1', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=p1_slope_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    p1_slope=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='16Slp_p1', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=p1_slope_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(p1_slope_tif_fn)
                p1_slope = np.array(ds.GetRasterBand(1).ReadAsArray())
                p1_slope[np.where(p1_slope == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None

            #interpolate Plane_slope - residuals
            if os.path.exists(p1_slope_res_tif_fn) == False:
                print('slope p1-least squared residuals, ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    p1_slope_res=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='17Slp_p1r', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=p1_slope_res_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    p1_slope_res=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='17Slp_p1r', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=p1_slope_res_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(p1_slope_res_tif_fn)
                p1_slope_res = np.array(ds.GetRasterBand(1).ReadAsArray())
                p1_slope_res[np.where(p1_slope_res == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None

            #interpolate Lst-sq Plane_slope - p2
            if os.path.exists(p2_slope_tif_fn) == False:
                print('LST-sq slope, ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    p2_slope=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='18Slp_p2', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=p2_slope_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    p2_slope=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='18Slp_p2', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=p2_slope_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(p2_slope_tif_fn)
                p2_slope = np.array(ds.GetRasterBand(1).ReadAsArray())
                p2_slope[np.where(p2_slope == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None

            #interpolate Lst-sq Plane_slope - p2 - residuals
            if os.path.exists(p2_slope_res_tif_fn) == False:
                print('slope p2-least squared residuals, ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    p2_slope_res=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='17Slp_p1r', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=p2_slope_res_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    p2_slope_res=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='17Slp_p1r', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=p2_slope_res_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(p2_slope_res_tif_fn)
                p2_slope_res = np.array(ds.GetRasterBand(1).ReadAsArray())
                p2_slope_res[np.where(p2_slope_res == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None
            
            #interpolate LST-Square curvate
            if os.path.exists(mean_curv_tif_fn) == False:
                print('LSTSQ curvature, ', end='', flush=True)
                if os.path.exists(inps.shapefile_clip):
                    plane_curv=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='21CurvLSQ', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=mean_curv_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer=inps.shapefile_clip, interpolation_method=interpolation_method)
                else:
                    plane_curv=gdal_grid_interpolate(x_coords=x_coords, y_coords=y_coords, ncols=ncols, nrows=nrows,resolution_m=current_rstep_size, layer_in=os.path.basename(inps.shapefile_clip)[:-4], \
                                                   zfield='21CurvLSQ', input_vrt=pc_seed_pts_stats_vrt_fn, output_grid=mean_curv_tif_fn,\
                                                   radius1=radius1, radius2=radius2, grid_datatype='Float32', cliplayer='', interpolation_method=interpolation_method)
            else:
                ds = gdal.Open(mean_curv_tif_fn)
                plane_curv = np.array(ds.GetRasterBand(1).ReadAsArray())
                plane_curv[np.where(plane_curv == ds.GetRasterBand(1).GetNoDataValue())] = np.nan
                ds = None
            print(' done.')

        print('\tWriting to HDF file ... ', end='', flush=True)
        hdf5_fn = str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_seed_pts_stats_raster_%0.2fm.h5'%current_rstep_size
        pc_results_fn = os.path.join(hdf_dir, hdf5_fn)
        hdf_out = h5py.File(pc_results_fn,'w')
        hdf_out.attrs['help'] = 'Results from slope normalization: radius %0.2fm'%( current_rstep_size) 
        pts_seed_stats_fc = hdf_out.create_dataset('pts_seed_stats',data=pts_seed_stats, chunks=True, compression="gzip", compression_opts=7)
        pts_seed_stats_fc.attrs['help'] = '''Nr. of seed pts: %d,  pts_seed_stats shape: %d x %d, with col: name 
        0: Seed-X, 1: Seed-Y, 2: Seed-Z, 3: Mean-X, 4: Mean-Y, 5: Mean-Z, 6: Z-min, 7: Z-max, 8: Dz-max, 9: Dz-min,  
        10: Dz-std.dev, 11: Dz-range, 12: Dz-90-10th percentile range, 13: Dz-75-25th percentile range, 14: variance dz, 
        15: slope-p1, 16: Slope-p1 residuals, 17: slope-p2, 18: Slope-p2 residuals, 19: nr. of lidar points, 
        20: curvature_lstsq, 21:Std. Z'''\
        %(nr_of_seed_points, pts_seed_stats.shape[0], pts_seed_stats.shape[1])

        if inps.create_geotiff == 1:
            nr_lidar_fc = hdf_out.create_dataset('nr_lidar',data=nr_lidar, chunks=True, compression="gzip", compression_opts=7)
            nr_lidar_fc.attrs['help'] = 'Nr. of lidar measurements per grid cell for seed point'
            dz_std_fc = hdf_out.create_dataset('dz_std',data=dz_std, chunks=True, compression="gzip", compression_opts=7)
            dz_std_fc.attrs['help'] = 'Dz standard deviation'
            dz_range9010_fc = hdf_out.create_dataset('dz_range9010',data=dz_range9010, chunks=True, compression="gzip", compression_opts=7)
            dz_range9010_fc.attrs['help'] = 'dz_range9010'
            dz_iqr_fc = hdf_out.create_dataset('dz_iqr',data=dz_iqr, chunks=True, compression="gzip", compression_opts=7)
            dz_iqr_fc.attrs['help'] = 'dz_iqr'
            p2_slope_fc = hdf_out.create_dataset('p2_slope',data=p2_slope, chunks=True, compression="gzip", compression_opts=7)
            p2_slope_fc.attrs['help'] = 'p2 Slope-Least squared'
            p1_slope_fc = hdf_out.create_dataset('p1_slope',data=p1_slope, chunks=True, compression="gzip", compression_opts=7)
            p1_slope_fc.attrs['help'] = 'p1 Slope-Least squared'
            p2_slope_res_fc = hdf_out.create_dataset('p2_slope_res',data=p2_slope_res, chunks=True, compression="gzip", compression_opts=7)
            p2_slope_res_fc.attrs['help'] = 'p2 Slope-Least squared residuals'
            p1_slope_res_fc = hdf_out.create_dataset('p1_slope_res',data=p1_slope_res, chunks=True, compression="gzip", compression_opts=7)
            p1_slope_res_fc.attrs['help'] = 'p1 Slope-Least squared residuals'
            plane_curv_fc = hdf_out.create_dataset('plane_curv',data=plane_curv, chunks=True, compression="gzip", compression_opts=7)
            plane_curv_fc.attrs['help'] = 'plane_curv'
            z_std_fc = hdf_out.create_dataset('z_std',data=z_std, chunks=True, compression="gzip", compression_opts=7)
            z_std_fc.attrs['help'] = 'z_std'
            z_mean_fc = hdf_out.create_dataset('z_mean',data=z_mean, chunks=True, compression="gzip", compression_opts=7)
            z_mean_fc.attrs['help'] = 'z_mean'
            geotransform_fc = hdf_out.create_dataset('geotransform',data=geotransform)
            geotransform_fc.attrs['help'] = 'geotransform'
            epsg_code_fc = hdf_out.create_dataset('epsg_code',data=inps.epsg_code)
            epsg_code_fc.attrs['help'] = 'epsg_code'
            nr_lidar = None
            dz_std = None
            dz_range9010 = None
            dz_iqr = None
            p2_slope = None
            p1_slope = None
            p2_slope_res = None
            p1_slope_res = None
            plane_curv = None
            z_std = None
            z_mean = None
        hdf_out.close()
        print('done.')
    
        ### Write to LAS/LAZ file (writing to LAZ file not yet supported by laspy)
        output_las_fn = '_seed_pts_%0.2fm_radius_xyzmean.las'%(current_rstep_size)
        output_las_fn = os.path.join(las_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + output_las_fn)
        if inps.create_las == 1 and os.path.exists(output_las_fn) == False:
            print('\tWriting mean X-Y-Z of seed points to LAS file: %s... '%os.path.basename(output_las_fn), end='', flush=True)    
            pts2write = pts_seed_stats[:,3:6]
            mask = np.all(np.isnan(pts2write) | np.equal(pts2write, 0), axis=1)                
            pts2write = pts2write[~mask]
            mask = None
            #normalize input and generate colors using colormap
            v = pts2write[:,2]
            #stretch to 10-90th percentile
            v_1090p = np.nanpercentile(v, [10, 90])
            v_rescale = exposure.rescale_intensity(v, in_range=(v_1090p[0], v_1090p[1]))
            colormap_PuOr = mpl.cm.PuOr
            rgb = colormap_PuOr(v_rescale)
            #remove last column - alpha value
            rgb = (rgb[:, :3] * (np.power(2,16)-1)).astype('uint16')
        
            outFile = File(output_las_fn, mode='w', header=inFile.header)
            new_header = copy.copy(outFile.header)
            #setting some variables
            new_header.created_year = datetime.datetime.now().year
            new_header.created_day = datetime.datetime.now().timetuple().tm_yday
            new_header.x_max = pts2write[:,0].max()
            new_header.x_min = pts2write[:,0].min()
            new_header.y_max = pts2write[:,1].max()
            new_header.y_min = pts2write[:,1].min()
            new_header.z_max = pts2write[:,2].max()
            new_header.z_min = pts2write[:,2].min()
            new_header.point_records_count = pts2write.shape[0]
            new_header.point_return_count = 0
            outFile.header.count = v.shape[0]
            outFile.X = pts2write[:,0]
            outFile.Y = pts2write[:,1]
            outFile.Z = pts2write[:,2]
            outFile.Red = rgb[:,0]
            outFile.Green = rgb[:,1]
            outFile.Blue = rgb[:,2]    
            outFile.close()    
            print('done.')

        output_las_fn = '_seed_pts_%0.2fm_radius_iqrdZ.las'%(current_rstep_size)
        output_las_fn = os.path.join(las_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + output_las_fn)
        if inps.create_las == 1 and os.path.exists(output_las_fn) == False:
            print('\tWriting IQR dZ of seed points to LAS file: %s... '%os.path.basename(output_las_fn), end='', flush=True)    
            pts2write = pts_seed_stats[:, 3:6] #[:, [3, 4, 13]]
            mask = np.all(np.isnan(pts2write) | np.equal(pts2write, 0), axis=1)                
            pts2write = pts2write[~mask]
            #normalize input and generate colors using colormap
            v = pts_seed_stats[:,13]
            v = v[~mask]
            mask = None
            #stretch to 10-90th percentile
            v_1090p = np.nanpercentile(v, [10, 90])
            v_rescale = exposure.rescale_intensity(v, in_range=(v_1090p[0], v_1090p[1]))
            colormap_PuOr = mpl.cm.YlGn
            rgb = colormap_PuOr(v_rescale)
            #remove last column - alpha value
            rgb = (rgb[:, :3] * (np.power(2,16)-1)).astype('uint16')
        
            outFile = File(output_las_fn, mode='w', header=inFile.header)
            new_header = copy.copy(outFile.header)
            #setting some variables
            new_header.created_year = datetime.datetime.now().year
            new_header.created_day = datetime.datetime.now().timetuple().tm_yday
            new_header.x_max = pts2write[:,0].max()
            new_header.x_min = pts2write[:,0].min()
            new_header.y_max = pts2write[:,1].max()
            new_header.y_min = pts2write[:,1].min()
            new_header.z_max = pts2write[:,2].max()
            new_header.z_min = pts2write[:,2].min()
            new_header.point_records_count = pts2write.shape[0]
            new_header.point_return_count = 0
            outFile.header.count = v.shape[0]
            outFile.X = pts2write[:,0]
            outFile.Y = pts2write[:,1]
            outFile.Z = pts2write[:,2]
            outFile.Red = rgb[:,0]
            outFile.Green = rgb[:,1]
            outFile.Blue = rgb[:,2]    
            outFile.close()    
            print('done.')

        output_las_fn = '_seed_pts_%0.2fm_radius_slope_p2_lstsq.las'%(current_rstep_size)
        output_las_fn = os.path.join(las_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + output_las_fn)
        if inps.create_las == 1 and os.path.exists(output_las_fn) == False:
            print('\tWriting least-squared p2 slopes of seed points to LAS file: %s... '%os.path.basename(output_las_fn), end='', flush=True)    
            pts2write = pts_seed_stats[:, 3:6] #[:, [3, 4, 13]]
            mask = np.all(np.isnan(pts2write) | np.equal(pts2write, 0), axis=1)                
            pts2write = pts2write[~mask]
            #normalize input and generate colors using colormap
            v = pts_seed_stats[:,17]
            v = v[~mask]
            mask = None
            #stretch to 10-90th percentile
            v_1090p = np.nanpercentile(v, [10, 90])
            v_rescale = exposure.rescale_intensity(v, in_range=(v_1090p[0], v_1090p[1]))
            colormap_PuOr = mpl.cm.plasma
            rgb = colormap_PuOr(v_rescale)
            #remove last column - alpha value
            rgb = (rgb[:, :3] * (np.power(2,16)-1)).astype('uint16')
        
            outFile = File(output_las_fn, mode='w', header=inFile.header)
            new_header = copy.copy(outFile.header)
            #setting some variables
            new_header.created_year = datetime.datetime.now().year
            new_header.created_day = datetime.datetime.now().timetuple().tm_yday
            new_header.x_max = pts2write[:,0].max()
            new_header.x_min = pts2write[:,0].min()
            new_header.y_max = pts2write[:,1].max()
            new_header.y_min = pts2write[:,1].min()
            new_header.z_max = pts2write[:,2].max()
            new_header.z_min = pts2write[:,2].min()
            new_header.point_records_count = pts2write.shape[0]
            new_header.point_return_count = 0
            outFile.header.count = v.shape[0]
            outFile.X = pts2write[:,0]
            outFile.Y = pts2write[:,1]
            outFile.Z = pts2write[:,2]
            outFile.Red = rgb[:,0]
            outFile.Green = rgb[:,1]
            outFile.Blue = rgb[:,2]    
            outFile.close()    
            print('done.')

        output_las_fn = '_seed_pts_%0.2fm_radius_slope_p1_lstsq.las'%(current_rstep_size)
        output_las_fn = os.path.join(las_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + output_las_fn)
        if inps.create_las == 1 and os.path.exists(output_las_fn) == False:
            print('\tWriting least-squared p1 slopes of seed points to LAS file: %s... '%os.path.basename(output_las_fn), end='', flush=True)    
            pts2write = pts_seed_stats[:, 3:6] #[:, [3, 4, 13]]
            mask = np.all(np.isnan(pts2write) | np.equal(pts2write, 0), axis=1)                
            pts2write = pts2write[~mask]
            #normalize input and generate colors using colormap
            v = pts_seed_stats[:,15]
            v = v[~mask]
            mask = None
            #stretch to 10-90th percentile
            v_1090p = np.nanpercentile(v, [10, 90])
            v_rescale = exposure.rescale_intensity(v, in_range=(v_1090p[0], v_1090p[1]))
            colormap_PuOr = mpl.cm.plasma
            rgb = colormap_PuOr(v_rescale)
            #remove last column - alpha value
            rgb = (rgb[:, :3] * (np.power(2,16)-1)).astype('uint16')
        
            outFile = File(output_las_fn, mode='w', header=inFile.header)
            new_header = copy.copy(outFile.header)
            #setting some variables
            new_header.created_year = datetime.datetime.now().year
            new_header.created_day = datetime.datetime.now().timetuple().tm_yday
            new_header.x_max = pts2write[:,0].max()
            new_header.x_min = pts2write[:,0].min()
            new_header.y_max = pts2write[:,1].max()
            new_header.y_min = pts2write[:,1].min()
            new_header.z_max = pts2write[:,2].max()
            new_header.z_min = pts2write[:,2].min()
            new_header.point_records_count = pts2write.shape[0]
            new_header.point_return_count = 0
            outFile.header.count = v.shape[0]
            outFile.X = pts2write[:,0]
            outFile.Y = pts2write[:,1]
            outFile.Z = pts2write[:,2]
            outFile.Red = rgb[:,0]
            outFile.Green = rgb[:,1]
            outFile.Blue = rgb[:,2]    
            outFile.close()    
            print('done.')
    
        output_las_fn = '_seed_pts_%0.2fm_radius_mean_curv.las'%(current_rstep_size)
        output_las_fn = os.path.join(las_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + output_las_fn)
        if inps.create_las == 1 and os.path.exists(output_las_fn) == False:
            print('\tWriting mean curvature of seed points to LAS file: %s... '%os.path.basename(output_las_fn), end='', flush=True)    
            pts2write = pts_seed_stats[:, 3:6] #[:, [3, 4, 13]]
            mask = np.all(np.isnan(pts2write) | np.equal(pts2write, 0), axis=1)                
            pts2write = pts2write[~mask]
            #normalize input and generate colors using colormap
            v = pts_seed_stats[:,20]
            v = v[~mask]
            mask = None
            #stretch to 10-90th percentile
            v_1090p = np.nanpercentile(v, [10, 90])
            bounds = np.round(np.mean(np.abs(v_1090p)), decimals=2)
            v_rescale = exposure.rescale_intensity(v, in_range=(-bounds, bounds))
            colormap_PuOr = mpl.cm.PiYG
            rgb = colormap_PuOr(v_rescale)
            #remove last column - alpha value
            rgb = (rgb[:, :3] * (np.power(2,16)-1)).astype('uint16')
        
            outFile = File(output_las_fn, mode='w', header=inFile.header)
            new_header = copy.copy(outFile.header)
            #setting some variables
            new_header.created_year = datetime.datetime.now().year
            new_header.created_day = datetime.datetime.now().timetuple().tm_yday
            new_header.x_max = pts2write[:,0].max()
            new_header.x_min = pts2write[:,0].min()
            new_header.y_max = pts2write[:,1].max()
            new_header.y_min = pts2write[:,1].min()
            new_header.z_max = pts2write[:,2].max()
            new_header.z_min = pts2write[:,2].min()
            new_header.point_records_count = pts2write.shape[0]
            new_header.point_return_count = 0
            outFile.header.count = v.shape[0]
            outFile.X = pts2write[:,0]
            outFile.Y = pts2write[:,1]
            outFile.Z = pts2write[:,2]
            outFile.Red = rgb[:,0]
            outFile.Green = rgb[:,1]
            outFile.Blue = rgb[:,2]    
            outFile.close()    
            print('done.')

        if inps.create_gmt_maps != '':
            if os.path.exists(inps.create_gmt_maps) == True:
                os.chdir(os.path.join(inps.basedir, 'maps'))
                TITLE = inps.gmt_title + ' ' + str(current_rstep_size) + 'm'
                TITLE = '"' + TITLE + '"'
                POSTSCRIPT_BASENAME = inps.gmt_basename + '_' + str(current_rstep_size) + 'm'
                SHAPEFILE = inps.shapefile_clip
                DEM_MN_GRD = z_mean_tif_fn
                DEM_INTERP_GRD = inps.dem_fname
                SLP_LSTSQ_P1_GRD = p2_slope_tif_fn
                SLP_LSTSQ_P2_GRD = p1_slope_tif_fn
                SLP_LSTSQ_P1_RES_GRD = p1_slope_res_tif_fn
                CURV_GRD = mean_curv_tif_fn
                NRLIDARPTS_GRD = nrlidar_tif_fn
                DZ_STDDEV_GRD = dz_std_tif_fn
                DZ_IQR_GRD = dz_iqr_tif_fn
                DZ_R9010_GRD = dz_range9010_tif_fn
                z_delta_tif_fn = os.path.join(geotif_dir, str(os.path.basename(inps.inlas).split('.')[:-1][0]) + '_%0.2fm_z_mean_delta.tif'%(current_rstep_size))
                DELTA_DEM = z_delta_tif_fn
                DEM_RESOLUTION=str(current_rstep_size)
                #print('\tCreating GMT maps with %s... '%inps.create_gmt_maps, end='', flush=True)
                cmd = ['bash', inps.create_gmt_maps, TITLE, POSTSCRIPT_BASENAME, SHAPEFILE, DEM_MN_GRD, \
                       DEM_INTERP_GRD, SLP_LSTSQ_P1_GRD, SLP_LSTSQ_P2_GRD, SLP_LSTSQ_P1_RES_GRD, CURV_GRD, \
                       NRLIDARPTS_GRD, DZ_STDDEV_GRD, DZ_IQR_GRD, DZ_R9010_GRD, DELTA_DEM, DEM_RESOLUTION]
                print(' '.join(cmd))
                logfile_fname = os.path.join(inps.basedir, 'log') + '/gmt_maps_' + str(current_rstep_size) + '_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '.txt'
                logfile_error_fname = os.path.join(inps.basedir, 'log') + '/gmt_maps_' + datetime.datetime.now().strftime('%Y%b%d_%H%M%S') + '_err.txt'
                with open(logfile_fname, 'w') as out, open(logfile_error_fname, 'w') as err:
                    subprocess_p = subprocess.Popen(cmd, stdout=out, stderr=err)
                    subprocess_p.wait()
                print('done.')
                os.chdir(os.path.join(inps.basedir))

        #remove data from memory
        pc_xyz_distance = None
        pc_xyz_distance_id = None
        pc_xyz_rstep_seed = None        
        points = None
        pts2write = None
        v = None
        rgb = None
        idx_nonan = None
        xx = None
        yy = None
        pts_seed_stats = None
        
        print('\ttime: %0.2fs or %0.2fm\n'%(time.time() - ts, (time.time() - ts)/60))        
        
    print('total time: %0.2fs or %0.2fm'%(time.time() - ts0, (time.time() - ts0)/60))        

