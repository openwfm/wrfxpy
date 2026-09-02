#api handling for ngfs detections
import requests
import pandas as pd
import os, sys
import statistics as stats
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime, timedelta, timezone
import geopandas as gpd
from shapely.geometry import shape
from io import StringIO
import utils

def today_datetime(date):
    #accepts a string either 'today' or 'YYY-MM-DD' and returns a datetime object
    if date=='today':
        today = datetime.now(timezone.utc).date()
    else:
        today = datetime.strptime(date, "%Y-%m-%d").date()
    return today

def get_csv_date_str(df):
    #returns a date string like '2026_02_04' for a dataframe
    t_max = max(df["acq_date_time"])
    return df.pixel_date_time.max().date().isoformat()

def time_paramaters(df=None,date='today',limit=1000):
    #returns a dictionary to query for the time from the current calendar day or day passed as 
    # string like '2026-01-22', UTC time
    # optionally take the latest time from an existing dataframe
    if df is None:
        today = today_datetime(date)
        hour = '00'
        minute = '00'
    else:
        t_max = max(df["acq_date_time"])
        if type(t_max) is str:   #sometimes this can have strings instead of datetime objects
            df = parse_times(df)
            t_max = max(df["acq_date_time"])
        today = t_max.date()
        hour = str(t_max.hour).zfill(2)
        minute = str(t_max.minute).zfill(2)
    tomorrow = today + timedelta(days=1)
    params = {
        "datetime": (
        f"{today.isoformat()}T{hour}:{minute}:00Z/"
        f"{tomorrow.isoformat()}T{hour}:{minute}:00Z"
        ),
        "limit": limit
    }
    return params

def datetime_to_api_time(dt):
    #converts a datetime object to something OK for api call
    return dt.isoformat(timespec='seconds')[:-6]+'Z'

def duration_parameters(start_time,duration_hours,limit = 1000):
    #returns basic time paramerts given a start time and duration
    #start_time can be a string like 2023-06-28_22:00:00 or datetime object
    if isinstance(start_time,str):
        start_time = utils.esmf_to_utc(start_time)
    end_time = start_time + timedelta(hours = duration_hours)

    start_string = datetime_to_api_time(start_time)
    end_string = datetime_to_api_time(end_time)

    params = {
        "datetime": f"{start_string}/{end_string}",
        "limit": limit
    }
    return params

def add_ogc_bbox(params, bbox, bbox_crs=None):
    """
    Add an OGC API - Features compliant bbox parameter.

    Parameters
    ----------
    params : dict
        Existing API parameters (copied, not modified in-place)
    bbox : tuple
        2D: (min_x, min_y, max_x, max_y)
        3D: (min_x, min_y, min_z, max_x, max_y, max_z)
    bbox_crs : str, optional
        CRS URI or EPSG code. If omitted:
          - 2D bbox defaults to CRS84
          - 3D bbox defaults to CRS84h
    Returns
    -------
    dict
        Updated parameter dictionary with bbox (and bbox-crs if provided)
    """
    new_params = params.copy()
    if len(bbox) == 4:
        min_x, min_y, max_x, max_y = bbox
        new_params["bbox"] = f"{min_x},{min_y},{max_x},{max_y}"
    elif len(bbox) == 6:
        min_x, min_y, min_z, max_x, max_y, max_z = bbox
        new_params["bbox"] = f"{min_x},{min_y},{min_z},{max_x},{max_y},{max_z}"
    else:
        raise ValueError(
            "bbox must be a 4-tuple (2D) or 6-tuple (3D)"
        )
    if bbox_crs is not None:
        new_params["bbox-crs"] = bbox_crs
    return new_params




def add_optional_params(params,key,value):
    #adds paramter key with to the params dictionary, needs to be a queryable field from the API
    params[key] = value
    return params

def get_api_url(sat,domain):
    #available collections
    base_url = "https://fire.data.nesdis.noaa.gov/api/ogc/detections/collections/"
    r = requests.get(base_url)
    r.raise_for_status()
    collections = r.json()["collections"]
    sat_dict = {
        'GOES-18' : 'west',
        'GOES-19' : 'east'
    }
    collections_names = [ ### <<<<<<<-----------------------unused
    'ngfs_schema.ngfs_detections_scene_east_conus',
    'ngfs_schema.ngfs_detections_scene_east_mesoscale1',
    'ngfs_schema.ngfs_detections_scene_east_mesoscale2',
    'ngfs_schema.ngfs_detections_scene_west_conus',
    'ngfs_schema.ngfs_detections_scene_west_mesoscale1',
    'ngfs_schema.ngfs_detections_scene_west_mesoscale2'
    ]
    for c in collections:
        if (sat_dict[sat] in c['id']) and (domain.lower() in c['id']):
            api_url = f"https://fire.data.nesdis.noaa.gov/api/ogc/detections/collections/{c['id']}/items"
            print(api_url)
            return api_url
        
def url_from_dataframe(df):
    #finds appropriate url from an existing datframe
    sat = df.satellite.unique()[0]
    domain = df.scan_domain.unique()[0]
    return get_api_url(sat,domain)

def api_call(url,params):
    #uses the requests package to assemble data from the server
    #empty list for each data feature from the server
    features = []
    returned = 0
    #loop through the pages, each page will contain a link to the next page to be obtained
    print('API params:')
    for k in params.keys():
        print('\t',k,params[k])
    while url:
        r = requests.get(url, params=params)
        r.raise_for_status()
        data = r.json()
        #add data to features
        features.extend(data["features"])
        # Find next link
        print(f'Number of matching features: {data["numberMatched"]}')
        returned = returned + data["numberReturned"]
        print(f'Total returned features: {returned}')
        url = next(
            (link["href"] for link in data.get("links", [])
            if link.get("rel") == "next"),
            None
        )
        params = None #the next link will include the params already
    return features

def features_to_geojson(features):
    return {
            "type" : "FeatureCollection",
            "features" : features
    }

def geodataframe_to_geojson(df):
    features = []
    for i in df.index:
        pt = [df.loc[i]["longitude"],df.loc[i]["latitude"]]
        feat = {
            "type" : "feature",
            "geometry" : pt,
            "properties" : df.loc[i].drop(columns=["geometry"]).astype(str).to_dict()
        }
        features.append(feat)
        print(feat)
    return {
        "type" : "FeatureCollection",
        "features" : features
    }

def dataframe_to_geojson(df):
    features = []
    for i in df.index:
        feat = {
            "type" : "feature",
            "geometry" : {"type":"Point","coordinates": [df.loc[i]['longitude'],df.loc[i]['latitude']]},
            "properties" : df.loc[i].astype(str).to_dict()
        }
        features.append(feat)
    return {
        "type" : "FeatureCollection",
        "features" : features
    }


def features_to_geodataframe(features):
    #convert to geodataframe
    '''
    gdf = gpd.GeoDataFrame(
    [
        {**f["properties"], "geometry": shape(f["geometry"])}
        for f in features
    ],
    crs="EPSG:4326"
    )
    '''
    gdf = gpd.GeoDataFrame.from_features(features=features,crs="EPSG:4326")  #probably what the block above does
    gdf = gdf.drop_duplicates()
    return gdf

def features_to_dataframe(features):
    """
    Convert a list of GeoJSON features into a pandas DataFrame
    using only the properties.
    """
    if not len(features):
        print('Empty features list. Returning empty DataFrame')
        return pd.DataFrame()
    df =  pd.DataFrame(
        f["properties"] for f in features
    )
    df = df.drop_duplicates()
    return df

def parse_times(df):
    time_keys = ["acq_date_time","pixel_date_time"]
    for k in df.keys():
        if 'time' in k and k not in time_keys:
            time_keys.append(k)
    time_keys = list(set(time_keys))
    for tk in time_keys:
        if tk in df.keys():
            df[tk] = pd.to_datetime(df[tk], utc=True)
        else:
            print(f'{tk} is an unknown key in DatFrame, skipping...')
    return df

def make_ngfs_dataframe(sat='GOES-19',domain='CONUS',date='today',params=None,parse_t = True):
    url = get_api_url(sat,domain)
    if not params:
        params = time_paramaters(date=date)
    features = api_call(url=url,params=params)
    df = features_to_dataframe(features)
    if parse_t:
        df = parse_times(df)
    return df

def make_ngfs_geodataframe(sat='GOES-19',domain='CONUS',date='today',params=None,parse_t = True):
    url = get_api_url(sat,domain)
    if not params:
        params = time_paramaters(date=date)
    features = api_call(url=url,params=params)
    df = features_to_geodataframe(features)
    if parse_t:
        df = parse_times(df)
    return df

def get_known_incident_id(known_incident_id,date_str):
    params = time_paramaters(date=date_str)
    
def get_ngfs_data(date='today',ngfs_cfg=None,data = pd.DataFrame()):
    #downloads NGFS detections via api
    if len(data) > 0: #passing an existing dataframe
        df = update_dataframe(data)
    else:
        df = data
        if ngfs_cfg is None:
            days_to_get = 2
            sat = ['GOES-18','GOES-19']
        else:
            days_to_get = ngfs_cfg['goes_cfg'].get('days_to_get',2)
            sat = []
            for ss in list(ngfs_cfg['goes_cfg']['sat_sectors'].keys()):
                sat.append(f'GOES-{ss}')
        #empty datframe 
        today = today_datetime(date)
        for i in range(days_to_get+1):
            api_date = today-timedelta(days=i)
            for s in sat:
                new_df = make_ngfs_dataframe(sat = s, date = api_date.isoformat())
                df = pd.concat([df,new_df],ignore_index = True)
                df = df.drop_duplicates()
    df.reset_index(drop=True, inplace=True)
    print(df)
    print(f'Dataframe size: {len(df)}')
    return df
                                         

def update_dataframe(df):
    #updates current data with api call
    for sat in df.satellite.unique():
        tmp = df[df.satellite == sat]
        url = url_from_dataframe(tmp)
        params = time_paramaters(tmp)
        print(params)
        new_features = api_call(url=url,params=params)
        new_df = features_to_dataframe(new_features)
        df = pd.concat([df,new_df],ignore_index = True)
        df = df.drop_duplicates()
    return df

def save_dataframe_as_csv(df,ingest_directory='ingest/NGFS'):
    #saves datframe with format like: NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2026_01_20_020.csv
    t_max = max(df.acq_date_time)
    year = str(t_max.year)
    month = str(t_max.month).zfill(2)
    day = str(t_max.day).zfill(2)
    day_of_year = str(t_max.day_of_year).zfill(3)
    sat = df.satellite.unique()[0]
    domain = df.scan_domain.unique()[0]
    save_name = f'NGFS_FIRE_DETECTIONS_{sat}_ABI_{domain}_{year}_{month}_{day}_{day_of_year}.csv'
    save_path = f'{ingest_directory}/{save_name}'
    print(f'Saving as {save_path}')
    df.to_csv(save_path,index=False)



if __name__ == "__main__":
        pass
