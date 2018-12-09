# Point Cloud (PC) Geomorphologic roughness and topographic detrending
*This algorithm for this code is describe in:
Bingham, N., Bookhagen, B., Johnson, K., and Chadwick, O. (in review): Use of lidar point cloud data to assess human-induced erosion and loss of vegetation cover on contrasting lithologies*

When using the code, please cite the above paper.

# Point Cloud (PC) Geomorphologic roughness and topographic detrending
Detrending Point Cloud (PC) data with slope and calculating topographic roughness and curvature from PCs.

The code reads in a ground-classified PC from a LAS/LAZ file and calculates several geomorphology-relevant metrics on the PC. Input files can be from a lidar or a Structure-from-Motion (SfM) PC, but should be ground classified (for descriptions on how to ground-classify your data, see the [PDF manual](docs/PC_geomorph_roughness_manual.pdf)). 
The algorithm allows defining a radius which is used to fit a linear plane through the point cloud to detrend the data (i.e., normalize the point cloud with mean elevation of 0). These data are used to calculate deviations from the mean (roughness) and identify rills, arroyos, incised canyons, and other forms of erosion processes. By varying the radius over which the plane is fitted, several scales of the landscape can be analyzed (similar to varying radii of topographic relief).  The algorithm choses seed points from the PC with a user-defined spacing (for example 1m) and calculated statistics for each seed point with that given radius. 

Output includes a set of shapefile and geotiffs that show statistics of the PC within the given radius. Also, CSV and H5 files are created that contain lists of seed point location and statistical results for further analysis in Python or Matlab.

The code is parallized using `multiprocessing` and uses by default all available cores. This significantly speeds up statistic calculation of the point coud. For large point clouds, a significant amount of RAM is required or you will need to split your point cloud into smaller tiles.


The code performs several additional steps that are described in detail in a [PDF manual](docs/PC_geomorph_manual.pdf). In summary, these are:
1. Finding seed points with a given spacing, usually 1m to 5m.

2. For each seed point and its neighborhood (for 1m spacing of seed points points within a radius of 0.5m  are used), statistics are calculated from the PC (and all points). These include, for example slope, curvature, variability of height (Z) values (for a full list and detailed description see the manual). The parameters also allow detrending the points within the seed-point radius by its slope and derive surface-roughness parameters.

3. The code allows to subsample a point cloud either by a max. number of neighborhood points (e.g., k=5) or by defining a fraction of points to use to create a point cloud with approximately similar point-cloud density based on probabilities. This step of point-cloud homogenization can also be performed by other approaches (see for example [PDAL filters](https://pdal.io/stages/filters.html). The subsampled point cloud is written as a new LAS file.

4. The code interpolates the seed points to a grid and writes the output as a geotiff. In addition, a point cloud generates a LAS file of all seed points with the relevant metric.

5. If GMT is installed (and that opion is chosen), a set of output maps is generated for initial visualization. 


# Installation
This is a Python 3.x code that will run on any OS, which supports the packages. It runs and has been tested on Linux (Ubuntu/Debian), Windows 10, and Mac OS X. We are using [conda/miniconda](https://conda.io/docs/) to install the required packages, which can be [downloaded here](https://conda.io/miniconda.html). Follow [these instruction](https://conda.io/docs/user-guide/install/index.html) to get miniconda installed.

You will need several packages for python to run this code. These are standard packages and are included in all distributions. We create an environment called `PC_py3` (PointCloud-python3) in the following way (currently we are using Python 3.6, but  3.7 should work equally well):

```
conda config --prepend channels conda-forge/label/dev
conda config --prepend channels conda-forge
conda create -y -n PC_py3 python=3.6 pip scipy pandas numpy matplotlib \
    scikit-image gdal pdal xarray packaging ipython multiprocess \
    h5py lastools pykdtree spyder gmt=5*
```

You can active this environment on the command line with `source activate PC_py3`.

You don't need ipython or spyder to run this code and you can remove these repositories in the command line above, but they usually come in handy. Also, we are installing GMT5 for visualization purposes. If you don't plan to generate maps and/or use GMT, you can safely remove `gmt=5*` from the line above.

Next, Install a fast and simple LAS/LAZ reader/writer. You can do similar steps through `lastools`, but this interface is fairly simple to use. *Please note thas laspy currently does not support writing LAZ files*:
```
source activate PC_py3
pip install laspy
```
If you have issues with pip, see: [here](https://stackoverflow.com/questions/47955397/pip3-error-namespacepath-object-has-no-attribute-sor).

This code uses [pykdtree](https://github.com/storpipfugl/pykdtree). There are other KDTree implementations, for example [scipy.spatial.cKDTree](https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.spatial.cKDTree.html). But pykdtree is faster (but doesn't allow you to save the results such as cKDTree). Because we aim at very large point clouds, the pyKDTree algorithm is significantly faster for generating and querying the KDtree and will increase processing speed (we have run tests with 1e6 to 1e9 points).

Last, install the repository into your favorite github directory, for example ~/github:
```
cd ~/github
git clone https://github.com/UP-RS-ESP/PC_geomorph_roughness

```
You are now ready to run the code from the command line (see below).


# Command line parameters
The code can be run from the command line and you can invoke a list of parameters and help with:
```
python pc_geomorph_roughness.py -h

usage: pc_geomorph_roughness.py [-h] --inlas INLAS [--raster_m RASTER_M]
                                [--raster_m_range RASTER_M_RANGE]
                                [--subsample_1m_pc_k SUBSAMPLE_1M_PC_K]
                                [--subsample_1m_pc_p SUBSAMPLE_1M_PC_P]
                                [--redo_subsample_1m_pc_p REDO_SUBSAMPLE_1M_PC_P]
                                [--k_nr_of_neighbors K_NR_OF_NEIGHBORS]
                                [--dem_fname DEM_FNAME]
                                [--shapefile_clip SHAPEFILE_CLIP]
                                [--epsg_code EPSG_CODE]
                                [--create_geotiff CREATE_GEOTIFF]
                                [--create_shapefiles CREATE_SHAPEFILES]
                                [--create_gmt CREATE_GMT]
                                [--create_las CREATE_LAS]
                                [--mean_z_only MEAN_Z_ONLY]
                                [--nr_of_cores NR_OF_CORES]
                                [--max_nr_of_neighbors_kdtree MAX_NR_OF_NEIGHBORS_KDTREE]
                                [--pt_lower_threshold PT_LOWER_THRESHOLD]
                                [--create_gmt_maps CREATE_GMT_MAPS]
                                [--gmt_title GMT_TITLE]
                                [--gmt_basename GMT_BASENAME]

PointCloud (PC) processing for DEM statistics. Deriving gridded ground data
(elevation and slope) using centroid coordinates. B. Bookhagen
(bodo.bookhagen@uni-potsdam.de), V0.2 Dec 2018.

optional arguments:
  -h, --help            show this help message and exit
  --inlas INLAS         LAS/LAZ file with point-cloud data. Ideally, this file
                        contains only ground points (class == 2)
  --raster_m RASTER_M   Raster spacing for subsampling seed points on LAS/LAZ
                        PC. Usually 0.5 to 10 m, default = 1. Seed points are
                        selected from half that distances
  --raster_m_range RASTER_M_RANGE
                        Raster spacing for subsampling seed points on LAS/LAZ
                        PC. Uses a list of ranges with spacing, e.g.,
                        --raster_m_range "1 10 1" will create raster files
                        with spatial resolutions of 1 to 10 m in 1 m steps.
  --subsample_1m_pc_k SUBSAMPLE_1M_PC_K
                        Number of points in 0.5m radius that are randomly
                        subsampled from the full point cloud. This is useful
                        if point-cloud density greatly varies, because
                        statistics calculated for seed points with different
                        point numbers may be biased. If subsample_pc_k > 0
                        then the point cloud will be homogenized by selecting
                        k=n neighbors for each 1-m seed point. For example, if
                        subsample_pc_k 10, then each 1-mseed point will have
                        only 10 neighbors.
  --subsample_1m_pc_p SUBSAMPLE_1M_PC_P
                        Factor to subsample point cloud by based on
                        probability. If subsample_1m_pc_p 0.8, a pointcloud
                        with 80% of the input points is generated and sampling
                        of point cloud is based on probability. That is, seed
                        points with a high number of neighbors is sampled less
                        often, than a seed point with fewer neighbors. Will
                        use original points, but creates a reduced point
                        cloud. Calculates probability for 1-m seed-point
                        spacing.
  --redo_subsample_1m_pc_p REDO_SUBSAMPLE_1M_PC_P
                        Flag to redo random subsampling based on probability.
                        By default, an existing file with a probability is
                        loaded, if you set "redo_subsample_1m_pc_p true", the
                        random subsampling based on probability will be rerun
                        and stored in a separate file.
  --k_nr_of_neighbors K_NR_OF_NEIGHBORS
                        Number of neighbors for dynamic density estimation
                        (k_nr_of_neighbors = 10 by default, change to lower
                        number for lower-density point clouds).
  --dem_fname DEM_FNAME
                        Filename of DEM to extract point spacing. Used to
                        identify seed-point coordinates. Useful if a DEM
                        exists and one wants to create point-cloud statistics
                        aligned to the DEM grid.
  --shapefile_clip SHAPEFILE_CLIP
                        Name of shapefile to be used to clip interpolated
                        surfaces. Give full pathname. This is likely the
                        shapefile you have previously generated to subset/clip
                        the point-cloud data.
  --epsg_code EPSG_CODE
                        EPSG code (integer) to define projection information.
                        This should be the same EPSG code as the input data
                        (no re-projection included yet) and can be taken from
                        LAS/LAZ input file. Add this to ensure that output
                        shapefile and GeoTIFFs are properly geocoded.
  --create_geotiff CREATE_GEOTIFF
                        Create interpolated geotif files from PC data (default
                        no: --create_geotiff 0, set to --create_geotiff 1 to
                        generate geotiff files). Note that creating geotiff
                        files may increase processing time.
  --create_shapefiles CREATE_SHAPEFILES
                        Create point shapefiles in UTM (see --epsg_code) and
                        Geographic-DD projection. These contain all attributes
                        calculated during the processing (default no:
                        --create_shapefiles 0, set to --create_shapefiles 1 to
                        generate shapefiles).
  --create_gmt CREATE_GMT
                        Create gmt point or vector files for plotting with GMT
                        shapefiles in UTM (see --epsg_code) and Geographic-DD
                        projection. These contain all attributes calculated
                        during the processing (default no: --create_gmt 0, set
                        to --create_gmt 1 to generate GMT files).
  --create_las CREATE_LAS
                        Create LAS point file from seed points (currently no
                        writing of LAZ files supported). The color shows mean
                        elevation of the seed points. These contain all
                        attributes calculated during the processing (default
                        no: --create_las 0, set to --create_las 1 to generate
                        LAS files).
  --mean_z_only MEAN_Z_ONLY
                        Calculate mean elevation for grid cell size and no
                        other parameters.
  --nr_of_cores NR_OF_CORES
                        Max. number of cores to use for multi-core processing.
                        Default is to use all cores (0), set to --nr_of_cores
                        6 to use 6 cores. For some memory-intensive
                        applications, it may be useful to reduce the number of
                        cores.
  --max_nr_of_neighbors_kdtree MAX_NR_OF_NEIGHBORS_KDTREE
                        Setting the max. number of neighbors for KDTree
                        search. This can remain at 100 points for airborne
                        lidar data. You may want to consider increasing this
                        when using terrestrial lidar data or SfM data.
  --pt_lower_threshold PT_LOWER_THRESHOLD
                        Lower point threshold for performing plane fitting and
                        slope normalization. If there are less than
                        pt_lower_threshold in the seed-point neighborhood, a
                        point fitting is not performed and values are set to
                        NaN.
  --create_gmt_maps CREATE_GMT_MAPS
                        BASH File with GMT commands for plotting maps. Full
                        path and filename is required. Will need to be fine
                        tuned (see example).
  --gmt_title GMT_TITLE
                        GMT title to appear in output map.
  --gmt_basename GMT_BASENAME
                        GMT basename for filename. "gmt_basename" with the
                        following endings are generated: "_DEM", "_DEM", .
```

# Examples

The first command-line example is using the simplest form, without pointcloud subsampling or homogenization. 
0. We change directory into the example01 directory: `cd example01`;
1. We generate raster files from a LAZ file `--inlas Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC12.laz` ;
2. starting at 1m grid spacing and going up to 10m in 1m steps `--raster_m_range "1 10 1"`;
3. We also define a shapefile `--shapefile_clip ~/github/PC_geomorph_roughness/example01/SC12.shp` (it's considered good practice to use the full path to the file - it makes creating maps with GMT easier) to clip the interpolated point cloud with;
4. define the UTM zone and EPSG code `--epsg_code 26911`;
5. use all available cores `--nr_of_cores 0`;
6. create geotiff output files `--create_geotiff 1`;
7. create GMT point files `--create_gmt 1` ; 
8. no shapefiles `--create_shapefiles 0`; 
9. and write an output LAS file `--create_las 1`. 

These commands put together:

```
python -W ignore ~/github/PC_geomorph_roughness/pc_geomorph_roughness.py \
    --inlas Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC12.laz \
    --raster_m_range "1 10 1" \
    --shapefile_clip ~/github/PC_geomorph_roughness/example01/SC12.shp \
    --epsg_code 26911 --nr_of_cores 0 --create_geotiff 1 --create_gmt 1  \
    --create_shapefiles 0 --create_las 1 \
    2>&1 | tee Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC12_pcgr_1_10_1.log
```

Next, we subsample the point cloud to have a maximum of 5 neighborhood points by adding the option `--subsample_1m_pc_k 5`. Note that we also change the name of the log file.
```
python -W ignore ~/github/PC_geomorph_roughness/pc_geomorph_roughness.py \
    --inlas Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC12.laz \
    --raster_m_range "1 10 1" \
    --shapefile_clip ~/github/PC_geomorph_roughness/example01/SC12.shp \
    --epsg_code 26911 --nr_of_cores 0 --create_geotiff 1 --create_gmt 1  \
    --create_shapefiles 0 --create_las 1 \
    --subsample_1m_pc_k 5 \
    2>&1 | tee Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC12_pcgr_subsample_k5_1_10_1.log
```

Alternatively, we can sample the original point cloud down to 80% (p=0.8) to generate a more homogenous point cloud with similar point-cloud density by adding the option `--subsample_1m_pc_p 0.5`:
```
python -W ignore ~/github/PC_geomorph_roughness/pc_geomorph_roughness.py \
    --inlas Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC12.laz \
    --raster_m_range "1 10 1" \
    --shapefile_clip ~/github/PC_geomorph_roughness/example01/SC12.shp \
    --epsg_code 26911 --nr_of_cores 0 --create_geotiff 1 --create_gmt 1 \
    --create_shapefiles 0 --create_las 1 \
    --subsample_1m_pc_p 0.8 \
    2>&1 | tee Pozo_USGS_UTM11_NAD83_all_color_cl2_cat1_pc_geomorph_roughness_subsample_p05_1_10_1.log
``` 

If you would like to generate GMT figures with some outputs, edit the bash GMT shell file `example01_create_map_view_of_PC_geomorph_output_gmt.sh` and define the proper variables on the command line. Here, we set:
1. `--create_gmt_maps ~/github/PC_geomorph_roughness/example01/example01_create_map_view_of_PC_geomorph_output_gmt.sh` to the name of the bash file (edit this file first!); 
2. Set the title to contain catchment name and subsampling probability `--gmt_title "Blanca in Pozo (p=0.8)"`;
3. Set the prefix of all filenames that are generated in the `maps` subdirectory `--gmt_basename "Pozo_Blanca_cl2_p08"`

```
python -W ignore ~/Dropbox/soft/github/PC_geomorph_roughness/pc_geomorph_roughness.py \
    --inlas Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC12.laz \
    --raster_m_range "1 10 1" \
    --shapefile_clip ~/github/PC_geomorph_roughness/example01/SC12.shp \
    --epsg_code 26911 --nr_of_cores 0 --create_geotiff 1 --create_gmt 1  \
    --create_shapefiles 0 --create_las 1 \
    --subsample_1m_pc_p 0.8 \
    --create_gmt_maps ~/github/PC_geomorph_roughness/example01/example01_create_map_view_of_PC_geomorph_output_gmt.sh \
    --gmt_title "Blanca in Pozo (p=0.8)" \
    --gmt_basename "Pozo_Blanca_cl2_p08" \
    2>&1 | tee Blanca_in_Pozo_USGS_UTM11_NAD83_all_color_cl2_SC1_pcgr_subsample_p08_1_10_1.log
'''
