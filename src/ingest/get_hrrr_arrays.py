# setup environment
from __future__ import absolute_import
from __future__ import print_function
from utils import Dict, hash2
import subprocess
import json
from datetime import datetime, timedelta
import pickle
import os
import os.path as osp
from pathlib import Path
import pandas as pd
import numpy as np
import sys
import ast
from osgeo import osr, ogr, gdal

# Script used to get array of HRRR data for all of conus associated with an FMDA dictionary 

# Variables used throughout
sys_cfg = Dict(json.load(open('etc/conf.json')))
hrrrpath = "/data001/projects/hirschij/data/hrrr/geotiff_files/" # path for atmospheric data stash

# NOTE: choosing to exclude solar bands 'DLWRF', 'USWRF', 'ULWRF'
# Downward shortwave is expected theoretically to be the most useful solar field 
# RAWS have downward shortwave sensors, so these could be compared to model fields
band_df_hrrr = pd.DataFrame({
    'Band': [616, 620, 624, 628, 629, 661, 561, 612, 643],
    'hrrr_name': ['TMP', 'RH', "WIND", 'PRATE', 'APCP',
                  'DSWRF', 'SOILW', 'CNWAT', 'GFLUX'],
    'dict_name': ["temp", "rh", "wind", "rain", "precip_accum",
                 "solar", "soilm", "canopyw", "groundflux"],
    'descr': ['2m Temperature [K]', 
              '2m Relative Humidity [%]', 
              '10m Wind Speed [m/s]'
              'surface Precip. Rate [kg/m^2/s]',
              'surface Total Precipitation [kg/m^2]',
              'surface Downward Short-Wave Radiation Flux [W/m^2]',
              'surface Total Precipitation [kg/m^2]',
              '0.0m below ground Volumetric Soil Moisture Content [Fraction]',
              'Plant Canopy Surface Water [kg/m^2]',
              'surface Ground Heat Flux [W/m^2]']
})

def check_times(d, t):
    for k in d.keys():
        assert np.array_equal(d[k]["HRRR"]["time"],t), f"Time array mismatch for key {k}"
    print("~"*50)
    print(f"Time array for dictionary: ")
    print(f"Start Time: {t.min()}")
    print(f"End Time: {t.max()}")
    print(f"Number of Time Steps: {len(t)}")

def check_f(d, f):
    for k in d.keys():
        assert all(key in d[k]["HRRR"] for key in f)
    print(f"F times for dictionary: {f}")


## Will save large 3d arrays
def save_hrrr_arrays(d):
    print(f"Building HRRR data arrays associated with dictionary, from stash {hrrrpath}")
    
    # Get Times and check same for all keys
    k = [*d.keys()]
    times = d[k[0]]["HRRR"]["time"]
    fts = [*d[k[0]]["HRRR"].keys()]
    fts = fts[1:len(fts)]
    check_times(d, times)
    check_f(d, fts)
    
    # Set up return
    r = dict({"time": times})

    for ft in fts:
        r[ft]={}
        for index, row in band_df_hrrr.iterrows():
            print("~"*50)
            band = row["Band"]
            dict_name = row["dict_name"]
            arrays = []
            for t in times:
                t = datetime.strptime(t, '%Y-%m-%dT%H:%M:%SZ')
                filename = f"hrrr.t{t.strftime('%H')}z.wrfprs{ft.lower()}.{band}.tif"
                tdir = osp.join(hrrrpath, t.strftime('%Y%m%d'))
                tpath = osp.join(tdir, filename)
                print(f"Getting array for band: {band}, {row['descr']}")
                print(f"Reading file {tpath}")
                assert osp.exists(tpath)
                ds = gdal.Open(tpath)
                data = ds.GetRasterBand(1).ReadAsArray()
                arrays.append(data)
            r[ft][dict_name]=np.stack(arrays, axis=2)
    return r


def get_hrrr_files(d):
    # Get Times and check same for all keys
    k = [*d.keys()]
    times = d[k[0]]["HRRR"]["time"]
    fts = [*d[k[0]]["HRRR"].keys()]
    fts = fts[1:len(fts)]
    check_times(d, times)
    check_f(d, fts)

     # Set up return
    r = dict({"time": times})

    for ft in fts:
        r[ft]={}
        for index, row in band_df_hrrr.iterrows():
            band = row["Band"]
            dict_name = row["dict_name"]
            files = []
            for t in times:
                t = datetime.strptime(t, '%Y-%m-%dT%H:%M:%SZ')
                filename = f"hrrr.t{t.strftime('%H')}z.wrfprs{ft.lower()}.{band}.tif"
                tdir = osp.join(hrrrpath, t.strftime('%Y%m%d'))
                tpath = osp.join(tdir, filename)
                assert osp.exists(tpath)
                files.append(tpath)
            r[ft][dict_name]=files
    return r

# Given formatted dictionary of file names, read 3d array into memory
def read_hour(dfiles, fhour, t_index):
    # Inputs:
    # dfiles: formatted dictionary of file names, output of get_hrrr_files
    # fhour: (str) forecast hour, e.g. 'f01'
    # t_index: (int) index of hour to read, index is realtive to start and stop time of time array
    # Set up return dictionary, key for each HRRR band, then 2d array within each
    r = {}
    
    for k in dfiles[fhour]:
        tpath = dfiles[fhour][k][0]
        print(f"Opening {tpath}")
        ds = gdal.Open(tpath)
        band = ds.GetRasterBand(1)
        data = band.ReadAsArray()
        r[k] = data

    return r

# Set of Functions to get Lat/Lon from HRRR bands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``
## Test for these functions lives at: https://github.com/openwfm/notebooks/blob/main_jm/fmda/data/handle_geotiff.ipynb
## TODO: implement test in wrfxpy


# Generate Arrays of pixel indices
def pixel_to_world(geo_matrix, x, y):
    # Given geotransform info of a geotiff file and an (x,y) pixel coord pair, return the coord pair that matches the geotiff in meters
    # Inputs:
    # geomatrix: output of ds.GetGeoTransform() for given geotiff file
    # tuple of length 6 contains:
    # A geotransform consists in a set of 6 coefficients
    # GT(0) x-coordinate of the upper-left corner of the upper-left pixel.
    # GT(1) w-e pixel resolution / pixel width.
    # GT(2) row rotation (typically zero).
    # GT(3) y-coordinate of the upper-left corner of the upper-left pixel.
    # GT(4) column rotation (typically zero).
    # GT(5) n-s pixel resolution / pixel height (negative value for a north-up image).
    # x: pixel index x coord (1)
    # y: pixel index y coord (1)
    # Return: coordinates of same point as given x,y as offset from UL (m)
    # Example: pixel_to_world(mat, 0, 0) returns UL x,y from geotiff

    ul_x = geo_matrix[0]
    ul_y = geo_matrix[3]
    x_dist = geo_matrix[1]
    y_dist = geo_matrix[5]
    _x = x * x_dist + ul_x
    _y = y * y_dist + ul_y
    return _x, _y


def build_transform_inverse(dataset, EPSG):
    # Given gdal dataset and target EPSG, return transformation function that transforms meter coord pairs to pixel coord pairs
    # Inputs:
    # dataset: geotiff file
    # EPSG: integer
    source = osr.SpatialReference(wkt=dataset.GetProjection())
    target = osr.SpatialReference()
    target.ImportFromEPSG(EPSG)
    return osr.CoordinateTransformation(source, target)

def world_to_epsg(wx, wy, trans):
    # Inputs:
    # wx, wy: output of build_transform_inverse
    # wx: x coordinate (m) related to geotiff reference point
    # wy: y coordinate (m) related to geotiff reference point
    # transform: function to transform to given epsg, function type is osgeo.osr.CoordinateTransformation
    # Return:
    # point from osgeo Geometry object
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(wx, wy)
    point.Transform(trans)
    return point

def find_spatial_coordinate_from_pixel(dataset, x, y, transform=None, epsg=4326):
    # Given gdal dataset, target x y pixel pair, and EPSG, return the EPSG defined coordinate pair
    # dataset: gdal dataset, from geotiff file
    # x (int): pixel x index
    # y (int): pixel y index
    ## Upper left corner is often (0,0)
    # transform: transform inverse. output of build_transform_inverse, default none and it calculates from epsg
    # supply transform to save computational time
    # epsg: default 4326 (WGS84)
    # Return: coord pair in given epsg, eg lat/lon (floats)
    if transform is None:
        transform = build_transform_inverse(ds, epsg)
    world_x, world_y = pixel_to_world(dataset.GetGeoTransform(), x, y)
    point = world_to_epsg(world_x, world_y, transform)
    return point.GetX(), point.GetY()

# Return arrays of lat/lons for HRRR grid
# Assuming this is same for all HRRR grib files and all geotiff files extracted from them
tpath = '/home/hirschij/data/hrrr/geotiff_files/20240107/hrrr.t23z.wrfprsf01.643.tif'
def get_lonlats():
    ds = gdal.Open(tpath)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    gt = ds.GetGeoTransform()
    gp = ds.GetProjection()    
    # Initialize empty arrays
    lons=np.zeros(np.shape(data))
    lats=np.zeros(np.shape(data))

    # get transformation once and reuse
    transform = build_transform_inverse(ds, EPSG=4326)
    # Loop over indices and fill
    for i in range(0, np.shape(lons)[0]): # iterate i over x coord (longitude)
        for j in range(0, np.shape(lons)[1]): # iterate j over y coord (latitude)
            coord = find_spatial_coordinate_from_pixel(ds, j, i, transform=transform) # note order flip is intentional
            lats[i,j]=coord[0]
            lons[i,j]=coord[1]
    return lons, lats
    


# Executed Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print(f"Invalid arguments. {len(sys.argv)} was given but 1 expected")
        print(('Usage: %s <path_to_dict> ' % sys.argv[0]))
        print("Example: python src/ingest/get_hrrr_arrays.py ~/data/rnn_data_dicts/test_CA_202401.pkl")
        sys.exit(-1)
    
    dpath = sys.argv[1]
    
    # Read FMDA Dictionary
    dict1 = pd.read_pickle(dpath)
    
    # Get files associated with given fmda dictionary
    print(f"Searching for files at: {hrrrpath}")
    files = get_hrrr_files(dict1)

    print("HRRR Files:")

    # Get 2d arrays of lats and lons
    print("Getting Lon/Lats for HRRR Grid")
    xypath = "/data001/projects/hirschij/data"
    print(f"Checking {xypath}")
    if not osp.exists(osp.join(xypath, "hrrr_lons.pkl")):
        print(f"Calculating lon/lat from {tpath}")
        lons, lats = get_lonlats()
        print(f"Saving to {xypath}")
        with open(osp.join(xypath, "hrrr_lons.pkl"), 'wb') as handle:
            pickle.dump(lons, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(osp.join(xypath, "hrrr_lats.pkl"), 'wb') as handle:
            pickle.dump(lats, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        print(f"Reading from {xypath}")
        lons = pd.read_pickle(osp.join(xypath, "hrrr_lons.pkl"))
        lats = pd.read_pickle(osp.join(xypath, "hrrr_lats.pkl"))
    print(f"Longitudes Shape: {lons.shape}")
    print(f"Latitudes Shape: {lats.shape}")

    print("~"*50)
    print("Test of Reading into memory")
    
    # Get time array
    k = [*dict1.keys()] # list of dictionary keys
    times = dict1[k[0]]["HRRR"]["time"] # get time array for first key, NOTE: a test that the times are all equal runs above in get_hrrr_files, so execution would have stopped if false
    print(f"Reading forecast hour F01 for time {times[0]}")
    
    data = read_hour(files, "f01", 0)
    print("Summary of hourly data")
    print(f"Top Level Keys: {data.keys()}")
    print(f"Hourly arrays inside")
    for k in data:
        print(f"{k} subdict data type and dimensions")
        print(type(data[k]))
        print(data[k].shape)

