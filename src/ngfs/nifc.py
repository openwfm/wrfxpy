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


def current_nifc():

    # The complete URL containing all SAS token parameters
    url = ("https://stg-arcgisazurecdataprod3.az.arcgis.com/exportfiles-2532-182237/WFIGS_Incident_Locations_Current_-8500334455305135272.csv?"
        "sv=2025-05-05&st=2026-08-24T11%3A15%3A10Z&se=2026-08-24T12%3A20%3A10Z&sr=b&sp=r&sig=JLl%2B1cAFiR7yIHIhRavv7VBtxuCur5j%2F9SIZki9MXc0%3D")

    # HTTP headers from your curl command
    headers = {
        'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:136.0) Gecko/20100101 Firefox/136.0',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Accept-Encoding': 'gzip, deflate, br, zstd',
        'Referer': 'https://data-nifc.opendata.arcgis.com/',
        'DNT': '1',
        'Connection': 'keep-alive',
        'Upgrade-Insecure-Requests': '1',
        'Sec-Fetch-Dest': 'document',
        'Sec-Fetch-Mode': 'navigate',
        'Sec-Fetch-Site': 'same-site',
        'Sec-Fetch-User': '?1',
        'Priority': 'u=0, i'
    }

    # Cookie data included in your curl request
    cookies = {
        'mbox': 'PC#b13e0eba0dff409096354136d77a326d.35_0#1850401438|session#8ebf919a9806485d9f3b12c8c760e3ef#1787158498',
        'OptanonConsent': 'isGpcEnabled=0&datestamp=Wed+Aug+19+2026+10%3A23%3A58+GMT-0600+(Mountain+Daylight+Time)&version=202605.1.0&isIABGlobal=false&hosts=&consentId=ed41683d-2f99-47be-ac0d-2fd5bb403c31&interactionCount=6&landingPath=NotLandingPage&groups=1%3A1%2C4%3A0%2C2%3A1%2C3%3A1&AwaitingReconsent=false&geolocation=US%3BCO&browserGpcFlag=0&isAnonUser=1&isDntEnabled=1&prevHadToken=0&crTime=1783950092973',
        'AMCV_ED8D65E655FAC7797F000101@AdobeOrg': '179643557|MCIDTS|20593|MCMID|38351715925472910283663038632625124206|MCAID|NONE|MCOPTOUT-1779212530s|NONE|MCAAMLH-1750947629|7|MCAAMB-1779205329|j8Odv6LonN4r3an7LhD3WZrU1bUpAkFkkiY1ncBR96t2PTI|MCSYNCSOP|411-19636|vVersion|5.5.0',
        'at_check': 'true',
        'AMCVS_ED8D65E655FAC7797F000101@AdobeOrg': '1',
        's_cc': 'true',
        's_sq': '[[B]]',
        's_ppv': 'developers.arcgis.com%253A%2520python%253A%2520api-reference,32,2,13356',
        's_tp': '41450',
        'esri_locale': 'en'
    }

    output_filename = '/data/jhaley/wrfxpy/ngfs/perims/WFIGS_Incident_Locations_Current_-8500334455305135272.csv'

    try:
        # Send GET request with headers and cookies, streaming the response
        response = requests.get(url, headers=headers, cookies=cookies, stream=True)
        
        # Raise error if response is not 200 OK
        response.raise_for_status()
        
        # Write chunks to file to safely manage memory
        with open(output_filename, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    file.write(chunk)
                    
        print(f"File successfully downloaded and saved as '{output_filename}'")

    except requests.exceptions.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except Exception as err:
        print(f"An error occurred: {err}")



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


    f = '/data/jhaley/wrfxpy/ngfs/perims/WFIGS_Incident_Locations_Current_-8500334455305135272.csv'
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