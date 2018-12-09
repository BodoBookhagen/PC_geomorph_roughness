#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 07:02:31 2018

@author: bodo
"""

from laspy.file import File
import numpy as np, os, argparse, h5py, time, gdal, sys
from pykdtree.kdtree import KDTree
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

### Command Line Parsing
parser = argparse.ArgumentParser(description='PointCloud (PC) processing for geomorphologic purposes. Estimates slopes, curvature, local roughness, and other parameters. B. Bookhagen (bodo.bookhagen@uni-potsdam.de), May 2018.')
# Important and required:
parser.add_argument('--inlas', type=str, default='',  help='LAS/LAZ file with point-cloud data. Ideally, this file contains only ground points (class == 2)')
# Optional / additional parameters
parser.add_argument('--nr_random_sampling', type=int, default=10,  help='Number of random-point cloud sampling iteration (how many iterations of bootstraping to estimate slope/curvature calculations).')
parser.add_argument('--pc_out_points_factor', type=float, default=2.0,  help='Factor to devide total point number by to obtain randomly sampled point cloud. If --pc_out_points_factor 2 the randomly sampled point cloud will have half the number of input points, with --pc_out_points_factor 4 only a fourth of input points are selected based on neighborhood probability. Default values is --pc_out_points_factor 2 (half of the points).')
parser.add_argument('--subsampling_type', type=int, default=0,  help='Defines how to generate the subsampled point cloud. Either by density (--subsampling_type 0) or by number of neighbors (--subsampling_type 1).')
parser.add_argument('--k_nr_of_neighbors', type=int, default=40,  help='Number of neighbors for dynamic density estimation. args.k_nr_of_neighbors = 40 by default, change to lower number for lower-density point clouds.')
parser.add_argument('--k_sampling', type=int, default=10,  help='Number of neighbors for subsampling. args.k_sampling = 10 will generate a point cloud with 10 neighbors for each seed point.')
parser.add_argument('--ground_only', type=bool, default=1,  help='Filter by point-cloud class to use ground class only. Set --ground_only 0 if you want to use all points (or if point cloud has been filtered before (default = 1, i.e. will use only ground-classified points).')
parser.add_argument('--outputdir', type=str, default='',  help='Output directory to store plots and pickle files. Default is directory containing LAS/LAZ file.')
parser.add_argument('--dem_fname', type=str, default='',  help='Filename of DEM to extract seed-point spacing. Used to identify seed-point coordinates')
args = parser.parse_args()

#args.inlas = '/raid-cachi/bodo/Dropbox/California/SCI/Pozo/cat1/Pozo_USGS_UTM11_NAD83_all_color_cl_cat1.laz'
#args.pc_out_points_factor = 2
#args.k_nr_of_neighbors = 40
#args.k_sampling = 10
#args.nr_random_sampling = 10
#args.ground_only = 1
#args.subsampling_type = 1
#args.dem_fname = '/raid-cachi/bodo/Dropbox/California/SCI/Pozo/cat1/Pozo_cat1_UTM11_NAD83_1m.tif'


def pc_seed_random_k_subsampling(pc_xyzg_k, pts_xyz, k):
    '''
    Sub-samples pointcloud (pts_xyz) for seed point list (pc_xyzg_radius) with k random neighbors
    Subsampling is random, not determined by probability.
    Expects a pykdtree list in pc_xyzg_k
    
    pc_xyzg_k_random = pc_k_subsampling(pc_xyzg_radius, pts_xyz, k=10)
    '''
    
    #iterate through n number of points (length of seed  points)
    n = pc_xyzg_k[0].shape[0]
    pc_xyzg_k_random = np.empty((n*k,4))
    pc_xyzg_k_random.fill(np.nan)
    counter = 0
    for i in range(n):
        current_indices = np.array(pc_xyzg_k[1][i])
        random_numbers = np.random.randint(0,len(current_indices),size=(k))
        if len(current_indices) >= k:
            pc_xyzg_k_random[counter:counter+k,0] = pts_xyz[ current_indices[ random_numbers.astype('int8') ], 0]
            pc_xyzg_k_random[counter:counter+k,1] = pts_xyz[ current_indices[ random_numbers.astype('int8') ], 1]
            pc_xyzg_k_random[counter:counter+k,2] = pts_xyz[ current_indices[ random_numbers.astype('int8') ], 2]
            pc_xyzg_k_random[counter:counter+k,3] = pts_xyz[ current_indices[ random_numbers.astype('int8') ], 3]
            counter = counter + k
        else:
            pc_xyzg_k_random[counter:counter+len(current_indices),0] = pts_xyz[ current_indices, 0]
            pc_xyzg_k_random[counter:counter+len(current_indices),1] = pts_xyz[ current_indices, 1]
            pc_xyzg_k_random[counter:counter+len(current_indices),2] = pts_xyz[ current_indices, 2]
            pc_xyzg_k_random[counter:counter+len(current_indices),3] = pts_xyz[ current_indices, 3]
            counter = counter + len(current_indices)
        
    pc_xyzg_k_random = pc_xyzg_k_random[~np.isnan(pc_xyzg_k_random).any(axis=1)]
    return pc_xyzg_k_random

def pc_DynamicDensity_pykdtree(pc_xyzg_pykdtree, pts_seed_xyz, k):
    radii, _ = pc_xyzg_pykdtree.query(pts_seed_xyz, k = k)
    dens =  k / np.pi / radii[:, -1]**2
    dens.shape = pts_seed_xyz.shape[0]
    disk = np.pi * radii[:, -1]**2
    probability = dens / disk
    radii_stats = np.percentile(radii,[10, 25, 50, 75, 90], axis=0)
    return dens, probability, radii_stats

def pc_random_p_subsampling(pts_xyz, p, nr_of_out_points):
    '''
    Sub-samples indices of PC pc_xyzg_radius with probability weight p based 
    on point density of each point. Will result in greatly reduced point cloud. 
    Give nr_of_out_points for subsampled point cloud, usually len(p)/2
    
    call with a probability
    #pc_xyzg_radius_equal_nr_random = pc_random_p_subsampling(pc_xyzg_radius, pts_xyz, nr_of_out_points)
    '''
    
    #iterate through n number of points (length of seed  points)
    n = len(p)
    if pts_xyz.shape[1] > 3:
        pc_xyzg_p_random = np.empty((int(nr_of_out_points),4))
    elif pts_xyz.shape[1] == 3:
        pc_xyzg_p_random = np.empty((int(nr_of_out_points),3))
    pc_xyzg_p_random.fill(np.nan)
    i = np.random.choice(n, size = int(nr_of_out_points), replace = False, p = p)
    # i = np.random.choice(n, size = int(nr_of_out_points), replace = False)
    pc_xyzg_p_random[:,0] = pts_xyz[i,0]
    pc_xyzg_p_random[:,1] = pts_xyz[i,1]
    pc_xyzg_p_random[:,2] = pts_xyz[i,2]
    if pts_xyz.shape[1] > 3:
        pc_xyzg_p_random[:,3] = pts_xyz[i,3]
    return pc_xyzg_p_random

if args.inlas == '':
    print('No input LAS/LAZ file given. Use --help or rerun with --in_las for input LAS/LAZ file. Exiting.')
    exit()

if args.outputdir == '':
    if len(os.path.dirname(args.inlas)) > 0:
        args.outputdir = os.path.dirname(args.inlas)
    else:
        args.outputdir = os.getcwd()
        
if os.path.exists(os.path.join(args.outputdir, 'log')) == False:
    os.mkdir(os.path.join(args.outputdir, 'log'))

bootstrap_dir = os.path.join(args.outputdir, 'bootstrap')
if os.path.exists(bootstrap_dir) == False:
    os.mkdir(bootstrap_dir)

figure_dir = os.path.join(args.outputdir, 'figure')
if os.path.exists(figure_dir) == False:
    os.mkdir(figure_dir)

pickle_dir = os.path.join(args.outputdir, 'pickle')
if os.path.exists(pickle_dir) == False:
    os.mkdir(pickle_dir)


### Loading data and filtering
print('\nLoading input file: %s'%args.inlas)
inFile = File(args.inlas, mode='r')
pc_xyzic = np.vstack((inFile.get_x()*inFile.header.scale[0]+inFile.header.offset[0], inFile.get_y()*inFile.header.scale[1]+inFile.header.offset[1], inFile.get_z()*inFile.header.scale[2]+inFile.header.offset[2], inFile.get_intensity(), inFile.get_classification())).transpose()
#pc_xyzic is now a point cloud with x, y, z, intensity, and classification
#if args.store_color == True:
#    pc_i = inFile.get_intensity().copy()
#    pc_blue = inFile.get_blue().copy()
#    pc_green = inFile.get_green().copy()
#    pc_red = inFile.get_red().copy()
print('Loaded %s points'%"{:,}".format(pc_xyzic.shape[0]))

if args.ground_only == True:
    print('\nFiltering points to only work with ground points (class == 2)... ',end='\n')
    #get only ground points:
    idx_ground = np.where(pc_xyzic[:,4] == 2)[0]
    pc_xyzig = pc_xyzic[idx_ground,0:4]
    pc_xyzg = pc_xyzic[idx_ground,0:3]
    #pc_xyzg is a point cloud with x, y, z, and for class == 2 only
    pc_xyg = pc_xyzic[idx_ground,0:2]
    #if args.store_color == True:
    ##    pc_i = pc_i[idx_ground,:]
    #    pc_blue = pc_blue[idx_ground,:]
    #    pc_green = pc_green[idx_ground,:]
    #    pc_red = pc_red[idx_ground,:]
else:
    print('\nUsing all points from file... ',end='\n')
    #get only ground points:
    idx_ground = np.where(pc_xyzic[:,4] >= 0)[0]
    pc_xyzig = pc_xyzic[idx_ground,0:4]
    pc_xyzg = pc_xyzic[idx_ground,0:3]
    pc_xyg = pc_xyzic[idx_ground,0:2]
    
total_points = pc_xyzig.shape[0]
nr_of_out_points = total_points / args.pc_out_points_factor

if np.max(pc_xyzic[:,4]) > 2:
    idx_vegetation = np.where((pc_xyzic[:,4] >= 3) & (pc_xyzic[:,4] <= 5))[0] # and pc_xyzic[:,4] == 4 and pc_xyzic[:,4] == 5)[0]
    #getting vegetation indices
    vegetation_intensity= pc_xyzic[idx_vegetation,3]
    vegetation_intensity_mean = np.mean(vegetation_intensity)
    vegetation_intensity_std = np.std(vegetation_intensity)
    print('\nNumber of vegetation points (class==[3,4,5]): %s'%"{:,}".format(idx_vegetation.shape[0]))
    print('Vegetation intensity mean: %2.1f +-std.dev.: %2.1f, 10th perc.: %2.1f, 90th perc.: %2.1f'%(vegetation_intensity_mean, vegetation_intensity_std, np.percentile(vegetation_intensity, 10), np.percentile(vegetation_intensity, 90)) )

#getting ground values
ground_intensity = pc_xyzic[idx_ground,3]
ground_intensity_mean = np.mean(ground_intensity)
ground_intensity_std = np.std(ground_intensity)
print('\nNumber of ground points (class==2): %s'%"{:,}".format(idx_ground.shape[0]))
print('Ground intensity mean: %2.1f +-std.dev.: %2.1f, 10th perc.: %2.1f, 90th perc.: %2.1f'%(ground_intensity_mean, ground_intensity_std, np.percentile(ground_intensity,10), np.percentile(ground_intensity,90)) )

# Determine density of EACH (ALL) point for a given number of neighbors
ts_total=time.time()
print('\nCalculating pyKDTree... ',end='',flush=True)
ts = time.time()
pc_xyzg_pykdtree = KDTree(pc_xyzg, leafsize=32)
if args.subsampling_type == 1:
    pc_xyg_pykdtree = KDTree(pc_xyg, leafsize=32)
total_time_KDTree = (time.time() - ts, (time.time() - ts)/60)       
print('done (%02.1fs or %02.1fm)'%(total_time_KDTree[0],total_time_KDTree[1]) ,end='\n',flush=True)

if args.subsampling_type == 1:
    # Loading coordinates from GeoTIFF file
    if args.dem_fname != '':
        try:
            src_ds = gdal.Open( args.dem_fname )
        except RuntimeError:
            print('Unable to open {}'.format(args.dem_fname))
            sys.exit(1)
    
        bandnr = 1
        try:
            src_band = src_ds.GetRasterBand(bandnr)
        except RuntimeError:
            print('Band ( {} ) not found in {}'.format(bandnr, args.dem_fname))
            sys.exit(1)
    
        #print("{}: [ MIN ] = {}, [ MAX ] = {}".format(os.path.basename(args.dem_fname), src_band.GetMinimum(), src_band.GetMaximum()))
        cols = src_band.XSize
        rows = src_band.YSize
        geo_transform = src_ds.GetGeoTransform()
        x_min = geo_transform[0] 
        x_max = geo_transform[0] + geo_transform[1] * cols
        y_min = geo_transform[3] 
        y_max = geo_transform[3] + geo_transform[5] * rows
        raster_m = geo_transform[1]
        
        x_elements = len(np.arange(x_min, x_max, raster_m))
        if y_min > y_max:
            y_elements = len(np.arange(y_max, y_min, raster_m))
            y_coords = np.arange(y_max, y_min, raster_m) + raster_m / 2
        else:
            y_elements = len(np.arange(y_min, y_max, raster_m))
            y_coords = np.arange(y_min, y_max, raster_m) + raster_m / 2
        
        #get coordinate range and shift coordinates by half of the step size to make sure rater overlay is centered. 
        #This is not really necessary and only matters for very small point clouds with edge effects or for very large steps sizes:
        x_coords = np.arange(x_min, x_max, raster_m) + raster_m / 2
        
        #create combination of all coordinates (this is using lists and could be optimized)
        xy_coordinates = np.array([(x,y) for x in x_coords for y in y_coords])
    else:
        #no GeoTiff file given, using min/max coordinates to generate equally-spaced grid
            
        ### Search KDTree with points on a regularly-spaced raster
        #generating equally-spaced raster overlay from input coordinates with stepsize rstep_size
        #This will be used to query the point cloud. Step_size should be small enough and likely 1/2 of the output file resolution. 
        #Note that this uses a 2D raster overlay to slice a 3D point cloud.
        raster_m = args.raster_m
        [x_min, x_max] = np.min(pc_xyzg[:,0]), np.max(pc_xyzg[:,0])
        [y_min, y_max] = np.min(pc_xyzg[:,1]), np.max(pc_xyzg[:,1])
        [z_min, z_max] = np.min(pc_xyzg[:,2]), np.max(pc_xyzg[:,2])
        x_elements = len(np.arange(x_min.round(), x_max.round(), raster_m))
        y_elements = len(np.arange(y_min.round(), y_max.round(), raster_m))
        
        #get coordinate range and shift coordinates by half of the step size to make sure rater overlay is centered. 
        #This is not really necessary and only matters for very small point clouds with edge effects or for very large steps sizes:
        x_coords = np.arange(x_min.round(), x_max.round(), raster_m) + raster_m / 2
        y_coords = np.arange(y_min.round(), y_max.round(), raster_m) + raster_m / 2
        
        #create combination of all coordinates (this is using lists and could be optimized)
        xy_coordinates = np.array([(x,y) for x in x_coords for y in y_coords])
        
    #using the 2D KDTree to find the points that are closest to the defined 2D raster overlay
    [pc_xyg_pykdtree_distance, pc_xyg_pykdtree_id] = pc_xyg_pykdtree.query(xy_coordinates, k=1)
    
    #take only points that are within the search radius (may omit some points at the border regions
    pc_distances_lt_rstep_size = np.where(pc_xyg_pykdtree_distance <= raster_m)[0]
    
    #the following list contains all IDs to the actual lidar points that are closest to the raster overlay. 
    #We will use these points as seed points for determining slope, planes, and point-cloud ranges
    pc_xyg_pykdtree_id = pc_xyg_pykdtree_id[pc_distances_lt_rstep_size]
    
    #remove from memory
    pc_xyg_pykdtree = None
    
    #now select these points from the 3D pointcloud with X, Y, Z coordinates.
    #We refer to these as seed points from the rstep part of this script
    pc_xyzg_rstep_seed = pc_xyzg[pc_xyg_pykdtree_id]
    
    ### Query points - use KDTree (pyKDTree)
    pc_xyzg_seed_k = pc_xyzg_pykdtree.query(pc_xyzg_rstep_seed, k=args.k_nr_of_neighbors, sqr_dists=False)

print('Calculating density and probability... ',end='',flush=True)
ts = time.time()
pc_xyzg_density_k, pc_xyzg_density_k_p, pc_xyz_radii_stats = pc_DynamicDensity_pykdtree(pc_xyzg_pykdtree, pc_xyzg, k = args.k_nr_of_neighbors)
total_time_density = (time.time() - ts, (time.time() - ts)/60)       
print('done (%02.1fs or %02.1fm)'%(total_time_density[0],total_time_density[1]) ,end='\n',flush=True)
print('Distance for k=nr. of points (median):\nFor k=5: %02.1f, k=10: %02.1f, k=15: %02.1f, k=20: %02.1f'%(pc_xyz_radii_stats[3,4], pc_xyz_radii_stats[3,9], pc_xyz_radii_stats[3,14], pc_xyz_radii_stats[3,19]) )

print('creating characteristic density and length figure ... ',end='',flush=True)
if args.subsampling_type == 1:
    pc_characteristic_length_png_fn = os.path.join(os.path.join(args.outputdir, figure_dir), 'PC_characteristic_length_raster%02dm_k%02d.png'%(raster_m, args.k_nr_of_neighbors))
elif args.subsampling_type == 0:
    pc_characteristic_length_png_fn = os.path.join(os.path.join(args.outputdir, figure_dir), 'PC_characteristic_length_k%02d.png'%(args.k_nr_of_neighbors))

fig = plt.figure(figsize=(11.69,8.27), dpi=150)
fig.clf()
ax1 = fig.add_subplot(111)
ax1.plot(pc_xyz_radii_stats[0,:], color='0.5', linestyle='-', linewidth=0.5, label='10th perc.')
ax1.plot(pc_xyz_radii_stats[1,:], color='b', linestyle='-', linewidth=0.5, label='25th perc.')
ax1.plot(pc_xyz_radii_stats[2,:], color='k', linestyle='-', linewidth=1, label='50th perc.')
ax1.plot(pc_xyz_radii_stats[3,:], color='b', linestyle='-', linewidth=0.5, label='75th perc.')
ax1.plot(pc_xyz_radii_stats[4,:], color='0.5', linestyle='-', linewidth=0.5, label='90th perc.')
ax1.grid()
ax1.set_xticks(np.arange(0, args.k_nr_of_neighbors+1, 2.0))
ax1.set_xlim(0, args.k_nr_of_neighbors)
ax1.set_ylim(0, 2)
if args.subsampling_type == 1:
    ax1.plot([0, args.k_nr_of_neighbors], [raster_m, raster_m], color='r', linestyle='-', linewidth=2, label='raster resolution (%02.1f m)'%raster_m)
    ax1.set_yticks(np.arange(0, raster_m*2+0.2, 0.2))
    ax1.set_ylim(0, raster_m*2)
ax1.set_title('Characteristic number of points for step sizes spacing for %s points with max. k=%02d'%("{:,}".format(int(total_points)), args.k_nr_of_neighbors), y=1.05)
ax1.set_xlabel('k number of neighbors')
ax1.set_ylabel('Distance (m)')
ax1.legend()
fig.savefig(pc_characteristic_length_png_fn, bbox_inches='tight')
plt.close()

print('done.\n')
for i in range(args.nr_random_sampling):
    if args.subsampling_type == 0:
        pc_bootstrapping_fn = os.path.join(os.path.join(args.outputdir, bootstrap_dir), 'PC_subsample_p_iter%02d_k%02d_%dpts.h5'%(i, args.k_nr_of_neighbors, nr_of_out_points))
        pc_subsampled_png_fn = os.path.join(os.path.join(args.outputdir, figure_dir), 'PC_subsampled_density_iter%02d_k%02d.png'%(i, args.k_nr_of_neighbors))
        pc_characteristic_length_png_fn = os.path.join(os.path.join(args.outputdir, figure_dir), 'PC_characteristic_length_iter%02d_k%02d_%dpts.png'%(i, args.k_nr_of_neighbors, nr_of_out_points))
        print('[%02d/%02d]: Subsampling of PC based on probability... '%(i+1, args.nr_random_sampling),end='',flush=True)
        if os.path.exists(pc_bootstrapping_fn) and os.path.exists(pc_subsampled_png_fn):
            print('file exists, skipping to next iteration.')
            continue
        elif os.path.exists(pc_bootstrapping_fn) and not os.path.exists(pc_subsampled_png_fn):
            print('file exists, ',end='',flush=True)
            hdf_in = h5py.File(pc_bootstrapping_fn,'r')
            pc_xyzg_p_random = np.array(hdf_in['pc_xyzg_p_random'])
            pc_xyzg_random_density_k_p = np.array(hdf_in['pc_xyzg_random_density_k_p'])
            pc_xyzg_random_density_k = np.array(hdf_in['pc_xyzg_random_density_k'])
            p_norm = 1./pc_xyzg_density_k_p
            p_norm /= p_norm.sum()
        else:
            ts = time.time()
            p_norm = 1./pc_xyzg_density_k_p
            p_norm /= p_norm.sum()
            pc_xyzg_p_random = pc_random_p_subsampling(pc_xyzig, p_norm, nr_of_out_points)
            pc_xyzg_p_random_pykdtree = KDTree(pc_xyzg_p_random[:,0:3], leafsize=32)
            pc_xyzg_p_random_density_k, pc_xyzg_random_density_k_p, pc_xyzg_random_radii_stats = pc_DynamicDensity_pykdtree(pc_xyzg_p_random_pykdtree, pc_xyzg_p_random[:,0:3], k = args.k_nr_of_neighbors)

    elif args.subsampling_type == 1:
        pc_bootstrapping_fn = os.path.join(os.path.join(args.outputdir, bootstrap_dir), 'PC_subsample_k_iter%02d_k%02d_%dpts.h5'%(i, args.k_nr_of_neighbors, nr_of_out_points))
        pc_subsampled_png_fn = os.path.join(os.path.join(args.outputdir, figure_dir), 'PC_subsampled_density_iter%02d_k%02d.png'%(i, args.k_nr_of_neighbors))
        pc_characteristic_length_png_fn = os.path.join(os.path.join(args.outputdir, figure_dir), 'PC_characteristic_length_iter%02d_k%02d.png'%(i, args.k_nr_of_neighbors))
        print('[%02d/%02d]: Subsampling of PC based on k = %02d number of neighbors... '%(i+1, args.nr_random_sampling, args.k_nr_of_neighbors),end='',flush=True)
        if os.path.exists(pc_bootstrapping_fn) and os.path.exists(pc_subsampled_png_fn):
            print('file exists, skipping to next iteration.')
            continue
        elif os.path.exists(pc_bootstrapping_fn) and not os.path.exists(pc_subsampled_png_fn):
            print('file exists, ',end='',flush=True)
            hdf_in = h5py.File(pc_bootstrapping_fn,'r')
            pc_xyzg_rstep_seed = np.array(hdf_in['pc_xyzg_rstep_seed'])
            pc_xyzg_p_random = np.array(hdf_in['pc_xyzg_p_random'])
            pc_xyzg_random_density_k_p = np.array(hdf_in['pc_xyzg_random_density_k_p'])
            pc_xyzg_random_density_k = np.array(hdf_in['pc_xyzg_random_density_k'])
            pc_xyzg_random_radii_stats = np.array(hdf_in['pc_xyzg_random_radii_stats'])
        else:
            ts = time.time()
            pc_xyzg_p_random = pc_seed_random_k_subsampling(pc_xyzg_seed_k, pc_xyzig, args.k_sampling)
            #find unique rows in pointcloud pc_xyzg_k_random and only keep unique rows
            pc_xyzg_p_random = np.vstack({tuple(row) for row in pc_xyzg_p_random})        
            pc_xyzg_p_random_pykdtree = KDTree(pc_xyzg_p_random[:,0:3], leafsize=32)
            pc_xyzg_p_random_density_k, pc_xyzg_random_density_k_p, pc_xyzg_random_radii_stats = pc_DynamicDensity_pykdtree(pc_xyzg_p_random_pykdtree, pc_xyzg_p_random[:,0:3], k = args.k_nr_of_neighbors)

#        p_norm_squareroot = 1/pc_xyzg_density_k_p**(1/2)
#        p_norm_squareroot /= p_norm_squareroot.sum()
#        pc_xyzg_p_squareroot_random = pc_random_p_subsampling(pc_xyzig, p_norm_squareroot, nr_of_out_points)
#        pc_xyzg_p_squareroot_random_pykdtree = KDTree(pc_xyzg_p_squareroot_random[:,0:3], leafsize=32)
#        pc_xyzg_p_squareroot_random_density_k, pc_xyzg_random_density_k_p_squareroot = pc_DynamicDensity_pykdtree(pc_xyzg_p_squareroot_random_pykdtree, pc_xyzg_p_squareroot_random[:,0:3], k = args.k_nr_of_neighbors)
#
#        p_norm_squareroot2 = 1/pc_xyzg_density_k_p**(1/1.5)
#        p_norm_squareroot2 /= p_norm_squareroot2.sum()
#        pc_xyzg_p_squareroot2_random = pc_random_p_subsampling(pc_xyzig, p_norm_squareroot2, nr_of_out_points)
#        pc_xyzg_p_squareroot2_random_pykdtree = KDTree(pc_xyzg_p_squareroot2_random[:,0:3], leafsize=32)
#        pc_xyzg_p_squareroot2_random_density_k, pc_xyzg_random_density_k_p_squareroot2 = pc_DynamicDensity_pykdtree(pc_xyzg_p_squareroot2_random_pykdtree, pc_xyzg_p_squareroot2_random[:,0:3], k = args.k_nr_of_neighbors)

    total_time_subsampling = (time.time() - ts, (time.time() - ts)/60)       
    print('done (%02.1fs or %02.1fm) '%(total_time_subsampling[0],total_time_subsampling[1]) ,end='',flush=True)

    if args.subsampling_type == 0:
        ts = time.time()
        print('creating characteristic density and length figure ... ',end='',flush=True)
        fig = plt.figure(figsize=(11.69,8.27), dpi=150)
        fig.clf()
        ax1 = fig.add_subplot(111)
        ax1.plot(pc_xyzg_random_radii_stats[0,:], color='0.5', linestyle='-', linewidth=1, label='10th perc.')
        ax1.plot(pc_xyzg_random_radii_stats[1,:], color='b', linestyle='-', linewidth=1, label='25th perc.')
        ax1.plot(pc_xyzg_random_radii_stats[2,:], color='k', linestyle='-', linewidth=2, label='50th perc.')
        ax1.plot(pc_xyzg_random_radii_stats[3,:], color='b', linestyle='-', linewidth=1, label='75th perc.')
        ax1.plot(pc_xyzg_random_radii_stats[4,:], color='0.5', linestyle='-', linewidth=1, label='90th perc.')
        ax1.grid()
        ax1.set_xticks(np.arange(0, args.k_nr_of_neighbors+1, 2.0))
        ax1.set_xlim(0, args.k_nr_of_neighbors)
        ax1.set_ylim(0, 2)
        if args.subsampling_type == 1:
            ax1.plot([0, args.k_nr_of_neighbors], [raster_m, raster_m], color='r', linestyle='-', linewidth=2, label='raster resolution (%02.1f m)'%raster_m)
            ax1.set_yticks(np.arange(0, raster_m*2+0.2, 0.2))
            ax1.set_ylim(0, raster_m*2)
        ax1.set_title('Characteristic number of points for step sizes (subsampled) for %s points with max. k=%02d, iteration=%02d'%("{:,}".format(int(nr_of_out_points)), args.k_nr_of_neighbors, i), y=1.05)
        ax1.set_xlabel('k number of neighbors')
        ax1.set_ylabel('Distance (m)')
        ax1.legend()
        fig.savefig(pc_characteristic_length_png_fn, bbox_inches='tight')
        plt.close()
        print('done. ',end='',flush=True)
        
        print('creating map-view figure... ',end='',flush=True)
        fig = plt.figure(figsize=(16.53*1.5,11.69*1.5), dpi=150)
        #fig = plt.figure(figsize=(11.69,8.27), dpi=150)
        fig.clf()
        
        ax1 = fig.add_subplot(221)
        ax1.grid()
        cax1 = ax1.scatter(pc_xyzg_p_random[:,0], pc_xyzg_p_random[:,1], c=pc_xyzg_p_random[:,2], s=0.5, cmap=plt.get_cmap('terrain'), vmin=np.percentile(pc_xyzg_p_random[:,2], 2), vmax=np.percentile(pc_xyzg_p_random[:,2], 98), linewidth=0)
        ax1.set_title('Bootstrap: Lidar point elevation for %s points'%"{:,}".format(int(nr_of_out_points)),y=1.05)
        cbar = fig.colorbar(cax1)
        cbar.set_label('Elevation (m)')
        ax1.set_xlabel('UTM-X (m)')
        ax1.set_ylabel('UTM-Y (m)')
        ax1.axis('equal')
        
    #    ax1 = fig.add_subplot(222)
    #    ax1.grid()
    #    cax1 = ax1.scatter(pc_xyzg_p_squareroot_random[:,0], pc_xyzg_p_squareroot_random[:,1], c=pc_xyzg_p_squareroot_random_density_k, s=0.5, cmap=plt.get_cmap('gnuplot'), vmin=np.percentile(pc_xyzg_p_squareroot_random_density_k, 2), vmax=np.percentile(pc_xyzg_p_squareroot_random_density_k, 98), linewidth=0)
    #    ax1.set_title('PC Density of subsampled point cloud probability^(1/1.5) with %s points for k=%02d neighbors'%("{:,}".format(int(nr_of_out_points)), args.k_nr_of_neighbors), y=1.05)
    #    cbar1 = fig.colorbar(cax1)
    #    cbar1.set_label('Point density (pts/m^2)')
    #    ax1.set_xlabel('UTM-X (m)')
    #    ax1.set_ylabel('UTM-Y (m)')
    #    ax1.axis('equal')
    
        ax2 = fig.add_subplot(222)
        ax2.grid()
        cax2 = ax2.scatter(pc_xyzg[:,0], pc_xyzg[:,1], c=pc_xyzg_density_k, s=0.5, cmap=plt.get_cmap('gnuplot'), vmin=np.percentile(pc_xyzg_density_k, 2), vmax=np.percentile(pc_xyzg_density_k, 98), linewidth=0)
        ax2.set_title('PC Density of original point cloud with %s points for k=%02d neighbors'%("{:,}".format(int(total_points)), args.k_nr_of_neighbors), y=1.05)
        cbar2 = fig.colorbar(cax2)
        cbar2.set_label('Point density (pts/m^2)')
        ax2.set_xlabel('UTM-X (m)')
        ax2.set_ylabel('UTM-Y (m)')
        ax2.axis('equal')
        
    #    ax1 = fig.add_subplot(223)
    #    ax1.grid()
    #    cax1 = ax1.scatter(pc_xyzg_p_squareroot_random[:,0], pc_xyzg_p_squareroot_random[:,1], c=pc_xyzg_p_squareroot_random_density_k, s=0.5, cmap=plt.get_cmap('gnuplot'), vmin=np.percentile(pc_xyzg_p_squareroot_random_density_k, 2), vmax=np.percentile(pc_xyzg_p_squareroot_random_density_k, 98), linewidth=0)
    #    ax1.set_title('PC Density of subsampled point cloud probability^(1/2) with %s points for k=%02d neighbors'%("{:,}".format(int(nr_of_out_points)), args.k_nr_of_neighbors), y=1.05)
    #    cbar1 = fig.colorbar(cax1)
    #    cbar1.set_label('Point density (pts/m^2)')
    #    ax1.set_xlabel('UTM-X (m)')
    #    ax1.set_ylabel('UTM-Y (m)')
    #    ax1.axis('equal')
    
        ax3 = fig.add_subplot(223)
        ax3.grid()
        cax3 = ax3.scatter(pc_xyzg_p_random[:,0], pc_xyzg_p_random[:,1], c=pc_xyzg_p_random[:,3], s=0.5, cmap=plt.get_cmap('gray'), vmin=np.percentile(pc_xyzg_p_random[:,3], 2), vmax=np.percentile(pc_xyzg_p_random[:,3], 98), linewidth=0)
        ax3.set_title('Lidar Intensity of subsampled point cloud with %s points for k=%02d neighbors'%("{:,}".format(int(nr_of_out_points)), args.k_nr_of_neighbors), y=1.05)
        cbar3 = fig.colorbar(cax3)
        cbar3.set_label('Lidar Intensity')
        ax3.set_xlabel('UTM-X (m)')
        ax3.set_ylabel('UTM-Y (m)')
        ax3.axis('equal')
    #    cax3 = ax3.scatter(pc_xyzg[:,0], pc_xyzg[:,1], c=p_norm, s=0.5, cmap=plt.get_cmap('seismic'), vmin=np.percentile(p_norm, 2), vmax=np.percentile(p_norm, 98), linewidth=0)
    #    ax3.set_title('Lidar Intensity of subsampled point cloud with %s points for k=%02d neighbors'%("{:,}".format(int(nr_of_out_points)), args.k_nr_of_neighbors), y=1.05)
    #    cbar3 = fig.colorbar(cax3)
    #    cbar3.set_label('PC probability')
    #    ax3.set_xlabel('UTM-X (m)')
    #    ax3.set_ylabel('UTM-Y (m)')
    #    ax3.axis('equal')
    
        ax4 = fig.add_subplot(224)
        ax4.grid()
        cax4 = ax4.scatter(pc_xyzg_p_random[:,0], pc_xyzg_p_random[:,1], c=pc_xyzg_p_random_density_k, s=0.5, cmap=plt.get_cmap('gnuplot'), vmin=np.percentile(pc_xyzg_p_random_density_k, 2), vmax=np.percentile(pc_xyzg_p_random_density_k, 98), linewidth=0)
        ax4.set_title('PC Density of subsampled point cloud with %s points for k=%02d neighbors'%("{:,}".format(int(nr_of_out_points)), args.k_nr_of_neighbors), y=1.05)
        cbar4 = fig.colorbar(cax4)
        cbar4.set_label('Point density (pts/m^2)')
        ax4.set_xlabel('UTM-X (m)')
        ax4.set_ylabel('UTM-Y (m)')
        ax4.axis('equal')
    
        fig.savefig(pc_subsampled_png_fn, bbox_inches='tight')
        plt.close()
        print('done. ',end='',flush=True)
    
        #write as HDF file
        print('writing as HDF file...',end='',flush=True)
        hdf_out = h5py.File(pc_bootstrapping_fn,'w')
        hdf_out.attrs['help'] = 'Bootstrapping results (i=%02d of %02d) based on probability for k=%02d neighboring points (from %s to %s points)'%(i, args.nr_random_sampling, args.k_nr_of_neighbors, "{:,}".format(total_points), "{:,}".format(nr_of_out_points)) 
        pc_xyzg_p_random_fc = hdf_out.create_dataset('pc_xyzg_p_random',data=pc_xyzg_p_random, chunks=True, compression="gzip", compression_opts=7)
        pc_xyzg_p_random_fc.attrs['help'] = 'Randomly subsampled point cloud, based on probability'
        pc_xyzg_p_random_density_k_fc = hdf_out.create_dataset('pc_xyzg_p_random_density_k',data=pc_xyzg_p_random_density_k, chunks=True, compression="gzip", compression_opts=7)
        pc_xyzg_p_random_density_k_fc.attrs['help'] = 'PC density of subsampled point cloud'
        pc_xyzg_random_density_k_p_fc = hdf_out.create_dataset('pc_xyzg_random_density_k_p',data=pc_xyzg_random_density_k_p, chunks=True, compression="gzip", compression_opts=7)
        pc_xyzg_random_density_k_p_fc.attrs['help'] = 'PC probability of subsampled pointcloud'
        pc_xyzg_random_radii_stats_fc = hdf_out.create_dataset('pc_xyzg_random_radii_stats',data=pc_xyzg_random_radii_stats, chunks=True, compression="gzip", compression_opts=7)
        pc_xyzg_random_radii_stats_fc.attrs['help'] = 'Length/radii for each k-step neighbor'
        hdf_out.close()
        print('done. ',end='',flush=True)
        total_time_hdfwriting = (time.time() - ts, (time.time() - ts)/60)       

    if args.subsampling_type == 1:
        ts = time.time()
        print('creating characteristic density and length figure ... ',end='',flush=True)
        fig = plt.figure(figsize=(11.69,8.27), dpi=150)
        fig.clf()
        ax1 = fig.add_subplot(111)
        ax1.plot(pc_xyzg_random_radii_stats[0,:], color='0.5', linestyle='-', linewidth=1, label='10th perc.')
        ax1.plot(pc_xyzg_random_radii_stats[1,:], color='b', linestyle='-', linewidth=1, label='25th perc.')
        ax1.plot(pc_xyzg_random_radii_stats[2,:], color='k', linestyle='-', linewidth=2, label='50th perc.')
        ax1.plot(pc_xyzg_random_radii_stats[3,:], color='b', linestyle='-', linewidth=1, label='75th perc.')
        ax1.plot(pc_xyzg_random_radii_stats[4,:], color='0.5', linestyle='-', linewidth=1, label='90th perc.')
        ax1.grid()
        ax1.set_xticks(np.arange(0, args.k_nr_of_neighbors+1, 2.0))
        ax1.set_xlim(0, args.k_nr_of_neighbors)
        ax1.set_ylim(0, 2)
        if args.subsampling_type == 1:
            ax1.plot([0, args.k_nr_of_neighbors], [raster_m, raster_m], color='r', linestyle='-', linewidth=2, label='raster resolution (%02.1f m)'%raster_m)
            ax1.set_yticks(np.arange(0, raster_m*2+0.2, 0.2))
            ax1.set_ylim(0, raster_m*2)
        ax1.set_title('Characteristic number of points for step sizes (subsampled) with k=%02d, iteration=%02d'%(args.k_sampling, i), y=1.05)
        ax1.set_xlabel('k number of neighbors')
        ax1.set_ylabel('Distance (m)')
        ax1.legend()
        fig.savefig(pc_characteristic_length_png_fn, bbox_inches='tight')
        plt.close()
        print('done. ',end='',flush=True)

        print('creating map-view figure... ',end='',flush=True)
        fig = plt.figure(figsize=(16.53*1.5,11.69*1.5), dpi=150)
        #fig = plt.figure(figsize=(11.69,8.27), dpi=150)
        fig.clf()
        
        ax1 = fig.add_subplot(221)
        ax1.grid()
#        cax1 = ax1.scatter(pc_xyzg_p_random[:,0], pc_xyzg_p_random[:,1], c=pc_xyzg_p_random[:,2], s=0.5, cmap=plt.get_cmap('terrain'), vmin=np.percentile(pc_xyzg_p_random[:,2], 2), vmax=np.percentile(pc_xyzg_p_random[:,2], 98), linewidth=0)
#        cax1 = ax1.scatter(pc_xyzg_p_random[:,0], pc_xyzg_p_random[:,1], c=pc_xyzg_p_random[:,2], s=0.5, cmap=plt.get_cmap('terrain'), vmin=np.percentile(pc_xyzg_p_random[:,2], 2), vmax=np.percentile(pc_xyzg_p_random[:,2], 98), linewidth=0)
#        cax1b = ax1.plot(pc_xyzg_rstep_seed[:,0], pc_xyzg_rstep_seed[:,1], c='k', marker='.', markersize=1, linewidth=0)
        cax1 = ax1.scatter(pc_xyzg_rstep_seed[:,0], pc_xyzg_rstep_seed[:,1], c=pc_xyzg_rstep_seed[:,2], s=0.5, cmap=plt.get_cmap('terrain'), vmin=np.percentile(pc_xyzg_rstep_seed[:,2], 2), vmax=np.percentile(pc_xyzg_rstep_seed[:,2], 98), linewidth=0)
        ax1.set_title('Bootstrap: Lidar point elevation for seed-point locations from %s points'%"{:,}".format(int(nr_of_out_points)),y=1.05)
        cbar = fig.colorbar(cax1)
        cbar.set_label('Elevation (m)')
        ax1.set_xlabel('UTM-X (m)')
        ax1.set_ylabel('UTM-Y (m)')
        ax1.axis('equal')
    
        ax2 = fig.add_subplot(222)
        ax2.grid()
        cax2 = ax2.scatter(pc_xyzg[:,0], pc_xyzg[:,1], c=pc_xyzg_density_k, s=0.5, cmap=plt.get_cmap('gnuplot'), vmin=np.percentile(pc_xyzg_density_k, 2), vmax=np.percentile(pc_xyzg_density_k, 98), linewidth=0)
        ax2.set_title('PC Density of original point cloud with %s points for k=%02d neighbors'%("{:,}".format(int(total_points)), args.k_nr_of_neighbors), y=1.05)
        cbar2 = fig.colorbar(cax2)
        cbar2.set_label('Point density (pts/m^2)')
        ax2.set_xlabel('UTM-X (m)')
        ax2.set_ylabel('UTM-Y (m)')
        ax2.axis('equal')
        
        ax3 = fig.add_subplot(223)
        ax3.grid()
        cax3 = ax3.scatter(pc_xyzg_p_random[:,0], pc_xyzg_p_random[:,1], c=pc_xyzg_p_random[:,3], s=0.5, cmap=plt.get_cmap('gray'), vmin=np.percentile(pc_xyzg_p_random[:,3], 2), vmax=np.percentile(pc_xyzg_p_random[:,3], 98), linewidth=0)
        ax3.set_title('Lidar Intensity of subsampled point cloud with %s points for k=%02d neighbors'%("{:,}".format(int(nr_of_out_points)), args.args.k_sampling), y=1.05)
        cbar3 = fig.colorbar(cax3)
        cbar3.set_label('Lidar Intensity')
        ax3.set_xlabel('UTM-X (m)')
        ax3.set_ylabel('UTM-Y (m)')
        ax3.axis('equal')
    
        ax4 = fig.add_subplot(224)
        ax4.grid()
        cax4 = ax4.scatter(pc_xyzg_p_random[:,0], pc_xyzg_p_random[:,1], c=pc_xyzg_p_random_density_k, s=0.5, cmap=plt.get_cmap('gnuplot'), vmin=np.percentile(pc_xyzg_p_random_density_k, 2), vmax=np.percentile(pc_xyzg_p_random_density_k, 98), linewidth=0)
        ax4.set_title('PC Density of subsampled point cloud from seed points (k=%02d) from %s points and density estimations with k=%02d neighbors'%(args.k_sampling, "{:,}".format(int(total_points)), args.k_nr_of_neighbors), y=1.05)
        cbar4 = fig.colorbar(cax4)
        cbar4.set_label('Point density (pts/m^2)')
        ax4.set_xlabel('UTM-X (m)')
        ax4.set_ylabel('UTM-Y (m)')
        ax4.axis('equal')
    
        fig.savefig(pc_subsampled_png_fn, bbox_inches='tight')
        plt.close()
        print('done. ',end='',flush=True)
    
        #write as HDF file
        print('writing as HDF file...',end='',flush=True)
        hdf_out = h5py.File(pc_bootstrapping_fn,'w')
        hdf_out.attrs['help'] = 'Bootstrapping results (i=%02d of %02d) based on k=%02d neighboring points (from %s to %s points)'%(i, args.nr_random_sampling, args.k_nr_of_neighbors, "{:,}".format(total_points), "{:,}".format(nr_of_out_points)) 
        pc_xyzg_p_random_fc = hdf_out.create_dataset('pc_xyzg_p_random',data=pc_xyzg_p_random, chunks=True, compression="gzip", compression_opts=7)
        pc_xyzg_p_random_fc.attrs['help'] = 'Randomly subsampled point cloud, based on probability'
        pc_xyzg_rstep_seed_fc = hdf_out.create_dataset('pc_xyzg_rstep_seed',data=pc_xyzg_rstep_seed, chunks=True, compression="gzip", compression_opts=7)
        pc_xyzg_rstep_seed_fc.attrs['help'] = 'Seed point locations'
        pc_xyzg_p_random_density_k_fc = hdf_out.create_dataset('pc_xyzg_p_random_density_k',data=pc_xyzg_p_random_density_k, chunks=True, compression="gzip", compression_opts=7)
        pc_xyzg_p_random_density_k_fc.attrs['help'] = 'PC density of subsampled point cloud'
        pc_xyzg_random_density_k_p_fc = hdf_out.create_dataset('pc_xyzg_random_density_k_p',data=pc_xyzg_random_density_k_p, chunks=True, compression="gzip", compression_opts=7)
        pc_xyzg_random_density_k_p_fc.attrs['help'] = 'PC probability of subsampled pointcloud'
        pc_xyzg_random_radii_stats_fc = hdf_out.create_dataset('pc_xyzg_random_radii_stats',data=pc_xyzg_random_radii_stats, chunks=True, compression="gzip", compression_opts=7)
        pc_xyzg_random_radii_stats_fc.attrs['help'] = 'Length/radii for each k-step neighbor'
        hdf_out.close()
        print('done. ',end='',flush=True)
        total_time_hdfwriting = (time.time() - ts, (time.time() - ts)/60)       

    print('done (%02.1fs or %02.1fm) '%(total_time_hdfwriting[0],total_time_hdfwriting[1]) ,end='\n',flush=True)
    #print('')


print('')
total_time = (time.time() - ts_total, (time.time() - ts_total)/60)
print('done. total time: %02.1fs or %02.1fm'%(total_time[0],total_time[1]) ,end='\n',flush=True)
