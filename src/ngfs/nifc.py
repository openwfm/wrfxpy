from __future__ import absolute_import
from __future__ import print_function
import os, sys, glob
import subprocess
import time
from datetime import datetime
import json
import re
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
sys.path.insert(1, 'src/ngfs')
from datetime import datetime, timedelta
from urllib.parse import urlencode
from pathlib import Path
from ngfs import ngfs_api as api
from ngfs import ngfs_ftp as ftp
from ngfs import config_manager 
import requests
import folium
import folium.plugins
import fire_init
import utils
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import shapely
import pandas as pd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from pyproj import Proj, transform
import map_utils
#from map_utils import fmt_time_string, compute_end_time, parse_time
try:
   from pyproj import Transformer
except:
   print('pyroj Transformer is not available')


WFIGS_LOCATIONS_QUERY = ("https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/"
                         "WFIGS_Incident_Locations_Current/FeatureServer/0/query")
NIFC_LOCS_CSV = '/data/jhaley/wrfxpy/ngfs/perims/WFIGS_Incident_Locations_Current_-8500334455305135272.csv'


def current_nifc(output_filename=NIFC_LOCS_CSV):
    #Refresh the local copy of WFIGS Incident Locations (Current) that nifc_locs reads.
    #This used to replay a browser download -- a signed ArcGIS export URL plus copied
    #cookies -- which is valid for about an hour, so it expired (2026-08-24) and every
    #call since returned 403 while nifc_locs kept reading the August file.  The public
    #feature-service query below needs no token and does not expire.  It writes the same
    #CSV: the layer's attributes plus x, y (lon/lat, EPSG:4326) from the point geometry.
    #Dates come as epoch milliseconds rather than the export's text; nothing here reads
    #them, and the per-incident REST fallback in get_nifc_incident returns that form too.
    #Written to a temporary file and renamed, so a failed refresh keeps the previous copy.
    rows, offset = [], 0
    try:
        while True:
            r = requests.get(WFIGS_LOCATIONS_QUERY, timeout=120,
                             params={'where': '1=1', 'outFields': '*', 'outSR': 4326, 'f': 'json',
                                     'resultOffset': offset, 'resultRecordCount': 2000})
            r.raise_for_status()
            page = r.json()
            if 'error' in page:
                raise RuntimeError(page['error'])
            for f in page.get('features', []):
                row = dict(f['attributes'])
                g = f.get('geometry') or {}
                row['x'], row['y'] = g.get('x'), g.get('y')
                rows.append(row)
            if not page.get('exceededTransferLimit'):
                break
            offset += len(page.get('features', []))
        if not rows:
            raise RuntimeError('the query returned no incidents')
        tmp = output_filename + '.tmp'
        pd.DataFrame(rows).to_csv(tmp, index=False)
        os.replace(tmp, output_filename)
        print(f"WFIGS incident locations: {len(rows)} written to '{output_filename}'")
    except Exception as err:
        print(f"WFIGS incident locations not refreshed ({WFIGS_LOCATIONS_QUERY}): {err}; "
              f"keeping the previous copy")



def nifc_locs():    ### remove this, put in nifc.py module

###url and geo_url seem to have changed
    url = (
        "https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/"
        "WFIGS_Incident_Locations_Current/FeatureServer/replicafilescache/"
        "WFIGS_Incident_Locations_Current_-8500334455305135272.csv"
    )
    geo_url = (
        "https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/"
        "WFIGS_Incident_Locations_Current/FeatureServer/replicafilescache/"
        "WFIGS_Incident_Locations_Current_-698092179933823454.geojson"
    )


    f = NIFC_LOCS_CSV
    if not os.path.exists(f):
        current_nifc()
    file_time = os.path.getmtime(f)
    #download once each 30 minutes
    now = time.time()
    if now-file_time > 1800:  ###alway dlownload?
        current_nifc()
    #read file
    try:
        df = pd.read_csv(f)
    except:
        print('Error reading the NIFC data')
        df = pd.DataFrame()
    return df

def get_nifc_incident(irwin_id,feature_type='loc'):
    #current and year-to-date data is vailable
    #open locally stored file
    nifc_collections = {
        'locations' : {
            'current' : 'WFIGS_Incident_Locations_Current',
            'ytd' : 'WFIGS_Incident_Locations_YearToDate'
        },
        'perimeters' : {
            'current' : 'WFIGS_Interagency_Perimeters_Current',
            'ytd' : 'WFIGS_Interagency_Perimeters_YearToDate'
        }
    }
    
    if feature_type=='loc':
        #load latest nifc locatio data                                              ### <<<------------- this block is giving errors, why? 
        nifc_csv = nifc_locs()    
        if len(nifc_csv) > 0:   #can be an emptry dataframe                                                  ###       maybe load the csv once and access globally?
            #search for id string in csv file
            inc = nifc_csv[nifc_csv['IrwinID'] == irwin_id]
            if len(inc) > 0:
                print('Found loc via nifc_locs()')
                return nifc_geojson(inc)  #pandas Dataframe --->>> geojson
    
    #if that's not found or a perimeter is needed, try API call
    try:
        if feature_type =='loc':
            collection = 'WFIGS_Incident_Locations_Current'
            irwin_key = 'IrwinID'
        else:
            collection = 'WFIGS_Interagency_Perimeters_YearToDate'  #WFIGS_Current_Interagency_Fire_Perimeters
            irwin_key = 'attr_IrwinID'
        url = (
            f"https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/{collection}/FeatureServer/0/query"
        )
        params = {
            "where": f"{irwin_key}='{irwin_id}'",
            "outFields": "*",
            "f": "geojson"
        }
        r = requests.get(url, params=params)
        r.raise_for_status()
        return r.json()
    except:
        print('error in get_nifc_incident')
        fake = {
            "type": "FeatureCollection",
            "features": []
        }

        return fake

#https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/WFIGS_Interagency_Perimeters_Current/FeatureServer/0/query?outFields=*&where=1%3D1
def query_nearby_incidents(
    lon,
    lat,
    radius_km=20,
    feature_type="loc"
):

    if feature_type == "loc":
        collection = (
            "WFIGS_Incident_Locations_Current"
        )

    else:
        collection = (
            "WFIGS_Interagency_Perimeters_Current" #"WFIGS_Interagency_Perimeters_YearToDate" 
        )

    url = (
        "https://services3.arcgis.com/"
        "T4QMspbfLg3qTGWY/"
        f"arcgis/rest/services/{collection}/"
        "FeatureServer/0/query"
    )

    params = {
        "geometry": f"{lon},{lat}",
        "geometryType": "esriGeometryPoint",
        "spatialRel": "esriSpatialRelIntersects",
        "distance": radius_km,
        "units": "esriSRUnit_Kilometer",
        "inSR": 4326,
        "outFields": "*",
        "f": "geojson"
    }

    r = requests.get(url, params=params)
    r.raise_for_status()
    return r.json()

def nifc_geojson(nifc_csv):
    #turns csv file like nifc_csv = nifc_locs into geojson
    features = []
    for r,row in nifc_csv.iterrows():
        prop = row.to_dict()
        geo = {
            "type" : "Point",
            "coordinates" : [row['x'],row['y']]
        }
        f = {
            "type" : "Feature",
            "geometry" : geo,
            "properties" : prop
        }
        features.append(f)
    g = {
        "type" : "FeatureCollection",
        "features" : features
    }
    return g


if __name__ == "__main__":
   pass