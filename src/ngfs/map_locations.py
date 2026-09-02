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
from ngfs import nifc
from ngfs import watch_duty
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


def get_new_forecasts(cutoff_hours = 48):
    wksp = '/data/jhaley/wrfxpy/wksp'
    wksp_dir = Path(wksp)
    forecast_dirs = [d for d in wksp_dir.iterdir() if d.is_dir()]
    #filter the newsest
    now = time.time()
    cutoff = now - cutoff_hours*3600
    now0 = pd.Timestamp.now('UTC')
    today = f'{now0.year}-{str(now0.month).zfill(2)}-{str(now0.day).zfill(2)}'
    now1 = now0-timedelta(hours=24)
    yesterday = today = f'{now1.year}-{str(now1.month).zfill(2)}-{str(now1.day).zfill(2)}'
    rf = []
    for f in forecast_dirs:
        nme = f.name
        if today in nme or yesterday in nme:
            rf.append(f)
    
    recent_forecasts = [
        d for d in forecast_dirs
        if d.stat().st_ctime > cutoff
    ]
    
    #recent_forecasts = rf
    print(f'Found {len(recent_forecasts)} forecasts made in the previous {cutoff_hours} hours')
    return recent_forecasts

def get_ngfs_data(map_time=None,ngfs_cfg=None):  #  <<<<------------------------------- use ftp scene downloads
    if not map_time:
        map_time = datetime.now()
    ingest_dir = 'ingest/NGFS/'
    base_str = 'NGFS_FIRE_DETECTIONS_GOES-{}_ABI_{}_{}_{}_{}_{}.csv'
    days_to_get = 3   #### <<<<----------------------------------------------------- should be tied to the cutoff time in the algorithm
    goes_list = []
    for i in range(days_to_get):
        csv_day = map_time - timedelta(days=i)
        day_of_year = csv_day.timetuple().tm_yday
        yyyy, mm, dd = csv_day.year, csv_day.month, csv_day.day
        sat_sectors = [(18, 'CONUS'), (19,'CONUS'), (18, 'Full-Disk')]
        for sat,sector in sat_sectors:
            csv_str = base_str.format(str(sat),sector,yyyy, str(mm).zfill(2), str(dd).zfill(2), str(day_of_year).zfill(3))
            print(csv_str)
            goes_list.append(f'{ingest_dir}{csv_str}')
    print(goes_list)
    df = pd.DataFrame()
    for g in goes_list:
        try:
            data_read = pd.read_csv(g)
            df = pd.concat([df,data_read],ignore_index=True)
            print(g,len(df))
        except:
            print(f'{g} not available yet')
    #df = api.parse_times(df)     ##<<<-------------------------------------------------------make sure this doesn't have any timestamps in data 
    #df = api.update_dataframe(df=df)
    if len(df) == 0:
        return pd.DataFrame()
    df = df[df['country'] == 'United States']
    df = df.drop_duplicates(subset=['latitude','longitude'],keep='first')
    return df

def get_ngfs_viirs(ngfs_cfg,days_to_get = None):
    #dowlonads the NGFS VIIRS detections
    #returns a pandas dataframe
    #download the data
    if not days_to_get:
        ngfs_cfg['viirs_cfg']['days_to_get'] = days_to_get
    #add data to the workspace
    end_time = pd.Timestamp.now('UTC')
    start_time = pd.Timestamp.now('UTC')-timedelta(hours=days_to_get*24)
    df = ftp.add_ngfs_scene(ngfs_cfg,parse_times = False,sat='viirs',start_time=start_time,end_time=end_time)
    print(f'Found {len(df)} VIIRS detections. Will Filter by country and type')
    print(df)
    #filter this
    types = ['Possible Wildland Fire', 'Known Wildland Fire Incident']
    #types = ['Known Wildland Fire Incident','other_type']
    if len(df) == 0:
        return pd.DataFrame()
    df = df[df['country'] == 'United States']
    df = df[df['type_description'].isin(types)]
    
    return df

def make_geojson(df,keep_keys = None):
    #return a minimal geojson for an empty datframe
    if len(df) == 0:
        return { "type": "FeatureCollection", "features": [] }
    if keep_keys:
        df = df[keep_keys]
    #covert the time columns to something OK for geojson
    # Identify all datetime columns (including timezone-aware ones)
    date_cols = df.select_dtypes(include=['datetime64', 'datetimetz']).columns
    # Convert those columns to string
    df[date_cols] = df[date_cols].astype(str)
    features = []
    if 'x' in df.keys():  #for handling nifc location csv files
        df['longitude'] = df['x']
        df['latitude'] = df['y']
    for i in df.index:
        geo = {
            'type' : 'Point',
            'coordinates': [df.loc[i]['longitude'],df.loc[i]['latitude']]
        }
        prop = df.loc[i].to_dict()   ###### <<<--------------------------------------------maybe keep only a subset?
        f = {
            "type" : "Feature",
            "geometry" : geo,
            "properties" : prop
        }
        features.append(f)
    return {"type":"FeatureCollection","features":features}


def resetter(top=42,left = 60,home_lat =39.5,home_lon = -98.35,home_zoom = -98.35):
    
    reset_js = f"""
            <script>
            function resetMapView() {{
                window.ngfsMap.setView(
                [{home_lat}, {home_lon}],
                {home_zoom}
                );
            }}
            </script>
        """
    #print(reset_js)
    reset_button = """
        <div style="
            position: fixed;
            top: 42px;
            left: 60px;
            z-index: 9999;
        ">
        <button 
            onclick="resetMapView()"
            class="ngfs-button"
        ">
            Reset Map View
        </button>
        </div>
        """
    return reset_js, reset_button

def map_name():

    name_js = """
    window.addEventListener("load", function () {

    // -------------------------------------------------
    // Find Leaflet map
    // -------------------------------------------------

    const map = Object.values(window).find(
        v => v instanceof L.Map
    );

    if (!map) {
        console.warn("Map not found");
        return;
    }

    // Stable global reference
    window.ngfsMap = map;
    // -------------------------------------------------
    // Save map state
    // -------------------------------------------------
    map.on("moveend", () => saveMapState(map));

    // wait briefly for LayerControl creation
    setTimeout(function() {
        restoreMapState(map);
        restoreLayerState();
        attachLayerListeners();
        console.log("Map initialization complete");
    }, 500);

    });
    """
    return name_js


def save_map_states():
    map_states_js = """

        function saveMapState(map) {
            console.log("saveMapState CALLED");
            const center = map.getCenter();
            localStorage.setItem(
                "mapState",
                JSON.stringify({
                    lat: center.lat,
                    lon: center.lng,
                    zoom: map.getZoom()
                })
            );
        }

        function restoreMapState(map) {
            const saved = localStorage.getItem("mapState");
            if (saved) {
                const state = JSON.parse(saved);
                map.setView(
                    [state.lat, state.lon],
                    state.zoom
                );
            }
        }

        function saveLayerState() {
            const state = {};
            const checkboxes = document.querySelectorAll(
                '.leaflet-control-layers-selector'
            );
            checkboxes.forEach(cb => {

                // Find nearby label text
                const label =
                    cb.parentElement.textContent.trim();

                state[label] = cb.checked;
            });
            localStorage.setItem(
                "layerState",
                JSON.stringify(state)
            );
            console.log("Saved layer state:", state); 
        }

        function restoreLayerState(map) {
            const saved = localStorage.getItem("layerState");
            if (!saved) return;
            const state = JSON.parse(saved);
            const checkboxes = document.querySelectorAll(
                '.leaflet-control-layers-selector'
            );
            checkboxes.forEach(cb => {
                const label =
                    cb.parentElement.textContent.trim();
                const shouldBeChecked = state[label];
                // Only toggle if different
                if (
                    shouldBeChecked !== undefined &&
                    cb.checked !== shouldBeChecked
                ) {

                    // Simulate user click
                    cb.click();
                }
            });
            console.log("Restored layer state");
        }

        function attachLayerListeners() {
            const checkboxes = document.querySelectorAll(
                '.leaflet-control-layers-selector'
            );
            checkboxes.forEach(cb => {
                cb.addEventListener(
                    "change",
                    saveLayerState
                );
            });
        }

    """
    return map_states_js

def parse_primary_ignition(data):
    try:
        entry = data["ignitions"]["1"][0]  #for single domain only
        lat, lon = entry["latlon"]
        time_str = entry["time_utc"]
        return {
            "lat": lat,
            "lon": lon,
            "ign_utc": time_str,
        }
    except (KeyError, IndexError, TypeError):
        return None

'''
def make_perim(wrfout_path):
    #makes a perimeter around the forecast fire area, based of the FIRE_AREA variable
    #return a geometry dictionary and its area
    def close_ring(coords):
        if not np.allclose(coords[0], coords[-1]):
            coords = np.vstack([coords, coords[0]])
        return coords
    with nc4.Dataset(wrfout_path, 'r') as wrfout:
        fire_area = np.array(wrfout['FIRE_AREA'][:][0])
        fxlong = np.array(wrfout['FXLONG'][:][0])
        fxlat = np.array(wrfout['FXLAT'][:][0])
    msk = (fire_area > 0).astype(float)     ### <<<< ------------------------------------------- remove this block
    if np.sum(msk) < 5:
        return None, 0.0
    cs = plt.contour(fxlong, fxlat, msk, levels=[0.5])
    polygons = []
    for collection in cs.collections:
        for path in collection.get_paths():
            coords = path.vertices
            polygons.append(coords)
    plt.close()
    clean_polys = [close_ring(p).tolist() for p in polygons if len(p) >= 4]
    if not clean_polys:
        return None, 0.0
    from shapely.geometry import Polygon, MultiPolygon
    polys = [Polygon(p) for p in clean_polys]
    if len(polys) == 1:
        poly = polys[0]
    else:
        poly = MultiPolygon(polys)
    if not poly.is_valid:
        poly = poly.buffer(0)
    if poly.is_empty:
        return None, 0.0
    poly_area = map_utils.perim_area(poly)
    smooth = poly.simplify(0.0001)
    if smooth.is_empty:
        return None, 0.0
    from shapely.geometry import mapping
    geom = mapping(smooth)
    return geom, round(poly_area/(1000**2)*247.1)  ###acres
'''

def make_fuels(wrfout_path):
    fuels_dict = {
        1: 'Short grass (1 ft)',
        2: 'Timber (grass and understory)',
        3: 'Tall grass (2.5 ft)',
        4: 'Chaparral (6 ft)',
        5: 'Brush (2 ft)',
        6: 'Dormant brush, hardwood slash',
        7: 'Southern rough',
        8: 'Closed timber litter',
        9: 'Hardwood litter',
        10: 'Timber (litter + understory)',
        11: 'Light logging slash',
        12: 'Medium logging slash',
        13: 'Heavy logging slash',
        14: 'No fuel'
    }

    with nc4.Dataset(wrfout_path, 'r') as wrfout:
        fire_area = np.array(
            wrfout['FIRE_AREA'][0]
        )
        nfuel_cat = np.array(
            wrfout['NFUEL_CAT'][0]
        )
    fire_msk = fire_area > 0
    fire_count = np.sum(fire_msk)
    if fire_count == 0:
        return []

    fuel_stats = []
    for fuel_num, fuel_name in fuels_dict.items():
        fuel_msk = nfuel_cat == fuel_num
        burnt_count = np.sum(
            fuel_msk & fire_msk
        )
        burn_fraction = burnt_count / fire_count
        fuel_stats.append({
            "fuel_num": fuel_num,
            "fuel_name": fuel_name,
            "fraction": burn_fraction
        })
    # sort descending
    fuel_stats.sort(
        key=lambda x: x["fraction"],
        reverse=True
    )
    # remove zero entries
    fuel_stats = [
        f for f in fuel_stats
        if f["fraction"] > 0
    ]
    return fuel_stats[:3]




def transform_xy(x,y,inP,outP):
    if type(inP) == int:
        inP = f'epsg:{inP}'
        inProj = Proj(init=inP)
    else:
        inProj = Proj(inP)
    if type(outP) == int:
        outP = f'epsg:{outP}'
        outProj = Proj(init=outP)
    else:
        outProj = Proj(outP)
    #print(inP,outP)
    if 'Transformer' not in dir():
        #inProj = Proj(init=inP)
        #outProj = Proj(init=outP)
        xp,yp = transform(inProj,outProj,np.array(x),np.array(y))
    else:
    #if correct pyproj version is available and has Transformer
        tf = Transformer.from_crs(inP,outP)
        xp,yp = tf.transform(y,x)
    return xp,yp


def transform_poly(poly,inP,outP):
    if type(poly) == pd.Series: #for sending in a rown of GeoDataFrame
        y = np.array(poly.geometry.exterior.coords.xy[1])
        x = np.array(poly.geometry.exterior.coords.xy[0])
    else: #for sending in a polygon
        y = np.array(poly.exterior.coords.xy[1])
        x = np.array(poly.exterior.coords.xy[0])
    xp,yp = transform_xy(x,y,inP,outP)
    return Polygon(zip(xp,yp))


def perim_area(poly):
    #compute area from polygon of lat-lon coordinates using lambert conformal projection
    y = np.array(poly.exterior.coords.xy[1])
    x = np.array(poly.exterior.coords.xy[0])
    try:
        lon0, lat0 = poly.centroid.x, poly.centroid.y
    except:
        lon0 = np.mean(x)
        lat0 = np.mean(y)
    poly = Polygon(zip(x,y))
    proj_string = f"+proj=lcc +lat_0={lat0} +lon_0={lon0} +lat_1={lat0-0.5} +lat_2={lat0+0.5} +x_0=0 +y_0=0.0 +datum=WGS84 +units=m +no_defs"
    #proj_aea = Proj(proj='aea', lat_0=lat0, lon_0=lon0, datum='WGS84', units='m')
    transPoly = transform_poly(poly,4326,proj_string)
    return transPoly.area

def validate_geom(geom):
    if geom is None:
        return False
    if not isinstance(geom, dict):
        return False
    if "type" not in geom or "coordinates" not in geom:
        return False
    return True

def make_geo_json_set(info,make_new = False):
    #makes a geojson for all the wrfouts in a directory
    #make_new - True will fore recreation of perimeter set
    collection = []
    for w in info["wrfouts"]:
        wksp = w.parent.parent
        wrf_time_str = fmt_time_string(extract_wrfout_time(w)) #string from datetime object
        name = f'{info["incident_name"]} {wrf_time_str}'
        perim_file = f'{info["grid_code"]}_{wrf_time_str.replace(" ","_")}.geojson'
        perim_path = wksp / perim_file
        if perim_path.exists() and not make_new: #load perim instead of making new
            print(f'Found perim file, loading {perim_file}')
            try:
                with open(perim_path, "r") as f:
                    feature = json.load(f)
                collection.append(feature)
            except:
                print('Error reading feature file',perim_path)
            continue
        print(f'Processing {w}')
        print(f'Making {perim_file}')
        geom,area = map_utils.make_perim(w)
        feature = {
            "type": "Feature",
            "geometry": geom,
            "properties": {
                "name" : name,
                "fire_name": info["incident_name"],
                "irwin_id": info["irwin_id"],
                "timestamp": wrf_time_str,
                "perimeter_type": "forecast_timestep",
                "area_acres" : round(area)
            }
        }
        #save the individual feature
        if geom: #and len(feature['geometry']['coordinates'][0]) > 4:
            with open(perim_path, "w") as f:
                json.dump(feature,f,indent=2)
            collection.append(feature) # dont append None type geometries into the collection
    geojson = {
        "type": "FeatureCollection",
        "features": collection
    }
    return geojson


def parse_forecast_dir(d):
    import json

    input_file = d / "input.json"
    if not input_file.exists():
        return None

    with open(input_file) as f:
        jobfile = json.load(f)

    #get ignition information
    ign = parse_primary_ignition(jobfile)
    if not ign:
        return None

    #get name and irwin id from direcorty name
    name_info = parse_forecast_name(d.name)
    if not name_info:
        name_info = {"incident_name": d.name}

    #processing outputs
    wrf_dir = d / 'wrf'
    wrfouts = sorted(
        wrf_dir.glob('wrfout*'),
        key=extract_wrfout_time
    )
    wrfout_count = len(wrfouts)  
    time_now = datetime.utcnow()
    if wrfout_count == 0:
        last_time = time_now
    else:
        last_time = extract_wrfout_time(wrfouts[-1])
    print('Last wrfout time: ',last_time)
    print(f'Found {wrfout_count} wrfout files for forecast')
    if (wrfout_count == 1 and (time_now > last_time)) or (wrfout_count >= 2*int(name_info.get("duration_hr"))):
        forecast_complete = True
    else:
        forecast_complete = False
    try:
        last_geom, last_area = map_utils.make_perim(wrfouts[-1])
        try:
            fuels = map_utils.make_fuels(wrfouts[-1])
            fuel_string = ", ".join(
                f'{f["fuel_name"]} ({100*f["fraction"]:.0f}%)'
                for f in fuels[:3]
            )
            print("fuel_string",fuel_string)
        except:
            print('Error in the fuels module')
            fuel_string = 'Not computed'
    except:
        last_geom, last_area = None, 0.0
        fuel_string = 'Not computed'

    try:                                                             #### <<<<------------------ probably remove this, 
        elevation = map_utils.get_elevation(wrfout_path=wrfouts[-1])
    except:
        elevation = None
    
    wksp_path = str(d)
    info_name = f'{jobfile["grid_code"]}.json'
    info = {
        "name": str(d.name),
        "lat": ign["lat"],
        "lon": ign["lon"],
        "ign_time_string": fmt_time_string(parse_time(ign["ign_utc"])),
        "ign_utc" : parse_time(ign["ign_utc"]), #datetime object this cant be serialized in mapping functions
        "grid_code" : jobfile["grid_code"],
        "info_name" : info_name,
        "info_file" : f'{wksp_path}/{info_name}',
        #forecast information add size, etc
        "forecast_complete" : forecast_complete,
        "wrfouts": wrfouts, #posix path, this can't be serialized in mapping functions
        # parsed metadata
        "wksp_path" : wksp_path,
        "incident_name": name_info.get("incident_name"),
        "irwin_id": name_info.get("irwin_id"),
        "forecast_start": name_info.get("forecast_start"),
        "duration_hr": name_info.get("duration_hr"),
        "elevation" : elevation,
        "geometry" : last_geom,
        "final_area" :last_area,
        "fuels" : fuel_string
    }
    '''
    print('Checking info construction')
    for k in info.keys():
        print(k,info[k])
    
    info['info_file'] = str(info_file)
    info_name = f'{info["grid_code"]}.json'
    info_file = f'{info["wksp_path"]}/{info["info_name"]}'
    info['info_file'] = str(info_file)

    '''

    #perimeters information 
    perimeters_name = f'{info["grid_code"]}.geojson' #'forecast_perimeters.geojson'
    geojson_file = d / perimeters_name
    print(geojson_file)
    if geojson_file.exists():
        with open(geojson_file, 'r', encoding='utf-8') as file:
            geojson = json.load(file)
        if len(geojson['features']) < 2:
            geojson = make_geo_json_set(info)
    else:
        if not forecast_complete:
            geojson_file = d / 'incomplete_perimeters.geojson'
        geojson = make_geo_json_set(info)
        #geom = info["geometry"]
        with open(geojson_file, "w") as f:
            json.dump(geojson, f, indent=2)
    print()
    info['geojson'] = geojson
    return info

def parse_forecast_name(name):
    if name[-3] == '-':
        duration = name[-2:]
    elif name[-2] == '-': # single digit forecast duration
        duration = name[-1:]
        name = name + 'x' #add a character at the end
    elif name[-3] == '-': # triple digit forecast duration
        duration = name[-3:]
        name = name[:-1] #ignore final character
    forecast_start = name[-22:-3]
    #forecast_start = name[-20:-3].replace('_',' ')                 <<<------------------------------------------- remove
    irwin_id = '{'+name[-59:-23]+'}'
    fire_name_raw = name[4:-80]
    # Clean fire name (replace underscores)
    incident_name = fire_name_raw.replace("_", " ")
    return {
        "incident_name": incident_name,
        "forecast_start": forecast_start,    #f'{forecast_start[:-3]} UTC'    <<<------------------------------------------- remove
        "irwin_id": irwin_id,
        "duration_hr": int(duration)
    }

def build_forecast_url(info):
    base_url = "https://www.engr.colostate.edu/~jhaley03/wrfxweb/"
    lat = round(info["lat"], 2)
    lon = round(info["lon"], 2)
    rasters = '&rasters=WINDVEC,FIRE_AREA,PM25_SFC_D'
    params = {
        "zoom": 11,
        "pan": f"{lat},{lon}",
        "job_id": info["name"]
    }
    return f"{base_url}?{urlencode(params)}{rasters}"

def extract_wrfout_time(p):
    # assumes format: '_''
    tstr = p.name.split('_', 2)[-1]
    return datetime.strptime(tstr, "%Y-%m-%d_%H:%M:%S")

def compute_end_time(info):
    start = parse_time(info["forecast_start"])
    return start + timedelta(hours=info["duration_hr"])

def parse_time(time_str):
    return datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")  #from job json files in WRFXPY

def fmt_time_string(dt): #convert datetime object
    return dt.strftime("%Y-%m-%d %H:%M UTC")

def build_popup(info):          #<<<-------------------------------------- move into helper module
    wrfx_url = build_forecast_url(info)
    inc_url = f'https://www.engr.colostate.edu/~jhaley03/NGFS/incidents/{info["grid_code"]}.html'
    incident_name = (info.get("incident_name") or "").replace("_", " ")
    html = f"""
    <div class="fire-label" style="font-family: Arial; font-size: 13px; width: 240px;">
        <h4 style="margin-bottom:6px;">
            {incident_name}
        </h4>

        <b>IRWIN ID:</b><br>
        <div style="margin-bottom:6px;">
            {info.get("irwin_id", "N/A")}
        </div>

        <b>Forecast Status:</b><br>
        <div style="margin-bottom:6px;">
            {info.get('status','unknown')}
        </div>

        <b>Estimated Ignition Time:</b><br>
        <div style="margin-bottom:6px;">
            {info.get("ign_time_string", "N/A")}
        </div>

        <b>Forecast End Time:</b><br>
        <div style="margin-bottom:6px;">
            {fmt_time_string(compute_end_time(info))}
        </div>

        <b>Fire Size (Acres):</b><br>
        <div style="margin-bottom:6px;">
            {info["final_area"]}
        </div>

        <a href="{inc_url}" target="_blank">
            Open incident details page
        </a><br>

        <a href="{wrfx_url}" target="_blank">
            Open forecast visualization
        </a>
    </div>
    """

    return folium.Popup(html, max_width=300), html



def add_legend(map_time,cutoff_hours):
    legend_html = f"""
    <div class="map-overlay" style="
        position: fixed;
        top: 80px;
        left: 10px;
        z-index: 500;
        background-color: white;
        padding: 10px 12px;
        border-radius: 6px;
        font-size: 14px;
        max-width: 275px;
        box-shadow: 0 0 6px rgba(0,0,0,0.2);
        transition: transform 0.2s ease;
    ">
    <h3><b>Wildfire Forecast Monitor</b></h3><br>
    <b>Information</b><br>
    Displaying locations and perimeters of wildfire forecasts made during the previous {cutoff_hours} hours <br>
    Map Created: {fmt_time_string(map_time)}<br>  
    Webpage auto-refreshes every 10 minutes<br>
    <br>

    <b>Legend</b><br>
    <i class="fa fa-fire" style="color:red;"></i>
    <b>Active Forecast</b>: Forecast ends after map creation time<br>
    <i class="fa fa-fire" style="color:orange;"></i>
    <b>Expired Forecast</b>: Forecast ended before map creation time<br>
    <br>

    <b>Links</b><br>
    <a href="https://www.engr.colostate.edu/~hilburn/ngfs/index.html"
        target="_blank"
        rel="noopener noreferrer"
    >Forecast Dashboard</a><br>
    <a href="https://www.engr.colostate.edu/~jhaley03/wrfxweb" 
        target="_blank"
        rel="noopener noreferrer"
        >WRFXWEB Visualization Server</a><br>
    <a href="https://www.engr.colostate.edu/~jhaley03/NGFS/viirs_tracks.html" 
        target="_blank"
        rel="noopener noreferrer"
        >VIIRS Satellite locations</a><br>
    </div>
    
    """
    return legend_html



def return_geojsons(info):
    geojson = info['geojson']
    if len(geojson['features']) < 1:
        return None, None
    last_perim = geojson['features'][-1]
    now = pd.Timestamp.now('UTC')
    for f in geojson['features']:
        #print(f['properties']['timestamp'])
        p_time = pd.Timestamp(f['properties']['timestamp'])
        if p_time > now:
            break
        first_perim = f
    return first_perim,last_perim

def save_incident_info(info):
    #save the whole thing
    #remove non-serializable elemnents from dictionary
    info['name'] = str(info['name'])
    info['wksp_path'] = str(info['wksp_path'])
    for i,w in enumerate(info['wrfouts']):
        info['wrfouts'][i]=str(w)
    #info['wrfouts'] = str(info['wrfouts'])
    info.pop("ign_utc",None)
    #info.pop("geojson",None)
    #print(info)
    info_name = f'{info["grid_code"]}.json'
    info_file = f'{info["wksp_path"]}/{info["info_name"]}'
    info['info_file'] = str(info_file)
    with open(info_file,'w') as file2:
            json.dump(info, file2, indent=2)

def make_incident_webpage(info): 
    #makes an individual webpage for the incident    
    #print('Making incident webpage')
    py_cmd = f'python src/ngfs/incident_webpage.py {info["info_file"]} &'
    py_cmd = f'./incident_webpage.sh {info["info_file"]}'
    sleep_time = 2
    cmd = f'sleep {sleep_time}; {py_cmd}'
    print(cmd)
    subprocess.Popen(cmd,shell=True)
    #print(f'Starting job after {sleep_time} delay: {py_cmd}')





'''
    #loop through our incidents and matCh with current locations file
    features = []
    for info in info_list:
        id = info['irwin_id']
        ss = df[df['IrwinID']==id].copy()
        print(len(ss))
        if len(ss):
            prop = ss.to_dict(orient='records')[0]
            geometry = {
                "type" : "Point",
                "coordinates" : [prop['InitialLongitude'],prop['InitialLatitude']]
            }
            
            p = {
                'IncidentName':ss['IncidentName'].item(),
                'IrwinID': ss['IrwinID'].item(),
                'FireDiscoveryDateTime': ss['FireDiscoveryDateTime'].item()
            }
            
            feat = {
                "type" : "Feature",
                "geometry" : geometry,
                "properties" : prop
            }
            features.append(feat)

    geojson = {
        "type": "FeatureCollection",
        "features": features
    }
    return geojson
'''

#'FireDiscoveryDateTime'



 

if __name__ == "__main__":
    print('Finding recent forecasts')
    #print(dir(map_utils))
    cutoff_hours = 2*24
    forecasts = get_new_forecasts(cutoff_hours = cutoff_hours)

    map_time = datetime.utcnow()
    try:
        ngfs_cfg, wrfxpy_cfg = config_manager.load_cfgs()
        ngfs_cfg['viirs_cfg']['ingest_directory'] = "/data/jhaley/new_wrfxpy/wrfxpy/ingest/NGFS/VIIRS"
    except:
        ngfs_cfg = {
            "viirs_cfg" : {
                "data_source" : "ftp",
                "data_format" : "csv",
                "host" : "bin.ssec.wisc.edu",
                "remote_directory": "pub/volcat/fire_csv/NGFS_scene/VIIRS/SSEC-DB",
                "api_directory" : " ",
                "ingest_directory" : "/data/jhaley/new_wrfxpy/wrfxpy/ingest/NGFS/VIIRS",
                "days_to_get" : 2,
                "sats" : ["NOAA-20","NOAA-21","SNPP"]
                }
        }

    
    #basic map object
    home_lat = 39.5
    home_lon = -98.35
    home_zoom = 3
    m = folium.Map(
        location=[home_lat,home_lon],
        zoom_start=home_zoom,
        tiles = None,
        prefer_canvas=True,
        control_scale=True
    )

    
    #mapping tile layers from https://leaflet-extras.github.io/leaflet-providers/preview/
    #folium.Map(tiles='https://{s}.tiles.example.com/{z}/{x}/{y}.png', attr='My Data Attribution')
    layer_image = folium.TileLayer(tiles= "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",               
	                attr = "Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community",
                    name = "Esri.WorldImagery"
    )
    layer_topo = folium.TileLayer(tiles =  "https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}.png",
                     attr = """Tiles &copy; Esri &mdash; Esri, DeLorme, NAVTEQ, TomTom, Intermap, iPC, USGS, FAO, NPS, NRCAN, 
                               GeoBase, Kadaster NL, Ordnance Survey, Esri Japan, METI, Esri China (Hong Kong), and the GIS User Communit""",
                     name = "Esri.WorldTopoMap"
    )
    layer_osm = folium.TileLayer('OpenStreetMap') #default layer  

    folium.plugins.MousePosition().add_to(m)
    mini_map = folium.plugins.MiniMap(toggle_display=True)
    mini_map.add_to(m)
    folium.plugins.Draw(export=True,position="bottomright").add_to(m)  
    #from folium.plugins import measure_control
    #m.add_child(measure_control)  <<<<<----------------not working?

    #map layers
    layer_osm.add_to(m)    
    layer_image.add_to(m)
    layer_topo.add_to(m)   

    #add detections and locations data
    geojson = watch_duty.load_watchduty()
    nifc_csv = nifc.nifc_locs()
    print(f"Found {len(nifc_csv)} current NIFC locations")
    df = get_ngfs_data()
    #layer control group
    #page_info = folium.FeatureGroup("Information and Title").add_to(m)
    expired = folium.FeatureGroup("Expired Forecasts").add_to(m)
    active = folium.FeatureGroup("Active Forecasts").add_to(m)
    detections = folium.FeatureGroup("GOES FIRE DETECTIONS",show=False).add_to(m)
    viirs = folium.FeatureGroup("VIIRS FIRE DETECTIONS",show=False).add_to(m)
    if len(nifc_csv) > 0:
        nifcGroup = folium.FeatureGroup("NIFC Data",show=False).add_to(m)
    watchDuty = folium.FeatureGroup("WatchDuty Data",show=False).add_to(m)
    folium.LayerControl().add_to(m)         


    #try df[tk], utc=True)
    #dt = pd.to_datetime(df['acq_date_time'],utc=True)
    #df = df[dt > (pd.to_datetime(map_time)-timedelta(hours=cutoff_hours))]   ##<<<<<----------------------- FIX LATER
    types = ['Possible Wildland Fire', 'Known Wildland Fire Incident']
    for t in types:
        popup_keys = [
            'latitude',
            'longitude',
            'bright_t7',
            'pixel_area',
            'acq_date_time',
            'pixel_date_time',
            'satellite',
            'scan_domain',
            'confidence',
            'version',
            'bright_t13',
            'frp',
            'quality_flag',
            'type_description',
            'daynight',
            'country',
            'state',
            'county',
            'gacc_id',
            'nws_region',
            'nws_wfo_code',
            'nws_wfo_name',
            'nws_fire_wx_code',
            'known_incident_name',
            'known_incident_type',
            'known_incident_id',
            'land_cover',
            'fuel',
            'feature_tracking_id',
            'feature_frp'
        ]
        geo = make_geojson(df[df['type_description']==t],keep_keys=popup_keys)
        

        if 'Possible' in t:
            color = 'orange'
        else:
            color = 'red'
        if geo.get('features'):
            folium.GeoJson(
                geo,
                marker=folium.Circle(radius=1000, fill_color=color, color="black", weight=0.2),
                tooltip=folium.GeoJsonTooltip(fields=["satellite", "acq_date_time","type_description","known_incident_name"]),
                popup=folium.GeoJsonPopup(fields=popup_keys,max_width=600)
            ).add_to(detections)
        else:
            print('Skipping empty GOES GeoJson')
        

    #add viirs detections
    #filter now
    #compute the times
    days = int(round(cutoff_hours/24.0) + 1)
    print(f'Getting {days} days of VIIRS data')
    viirs_df = get_ngfs_viirs(ngfs_cfg,days_to_get=days)
    if len(viirs_df) > 0:
        viirs_df.to_csv('map_viirs_detections.csv',index=False)
        for t in types:
            popup_keys = [
                    'latitude',
                    'longitude',
                    'pixel_area',
                    'acq_date_time',
                    'satellite',
                    'confidence',
                    'version',
                    'frp',
                    'state',
                    'county',
                    'gacc_id',
                    'nws_fire_wx_code',
                    'known_incident_name',
                    'known_incident_id',
                    'feature_tracking_id',
                    'feature_frp',
                    'type_description'
                ]
            geo = make_geojson(viirs_df[viirs_df['type_description']==t],keep_keys=popup_keys)
            #popup_keys = list(viirs_df.keys())
            if 'Possible' in t:
                color = 'orange'
                continue   #don't add these right now, too many........... need to filter for area around fires
            else:
                color = 'red' 
            ##### this was breaking with NGFS version 3.9.14   why??  
            if geo.get('features'):
                folium.GeoJson(
                    geo,
                    marker=folium.Circle(radius=250, fill_color=color, color="black", weight=0.4),
                    tooltip=folium.GeoJsonTooltip(fields=["satellite","acq_date_time","type_description","known_incident_name"]),
                    popup=folium.GeoJsonPopup(fields=popup_keys,max_width=600)
                ).add_to(viirs)
            else:
                print('Skipping empty VIIRS GeoJson')
        
        
    ### html elements
    disclaimer_html = map_utils.add_disclaimer(left=75)
    m.get_root().html.add_child(folium.Element(disclaimer_html))
    m.get_root().header.add_child(folium.Element("""
    <title>Wildfire Forecast Monitor</title>
    """))

    
    
    

    #auto refresh

    ##page resher options
    #add button style
    button_string = map_utils.ngfs_button_class()
    m.get_root().html.add_child(folium.Element(button_string))
    '''
    m.get_root().html.add_child(  ###<<<<---------------------------------------- remove this block if what's above works
        folium.Element("""
            <style>
            .ngfs-button {
                background-color: #2c7be5;
                color: white;
                border: none;
                padding: 6px 10px;
                border-radius: 6px;
                font-size: 10px;
                cursor: pointer;
                box-shadow: 0 2px 6px rgba(0,0,0,0.2);
                width: 160px;
                text-align: center;
            }

            .ngfs-button:hover {
                background-color: #1b5fcc;
            }
            </style>
            """
        )
    )
    '''
    #refresh botton
    refresh_js, refresh_button = map_utils.refresher()
    m.get_root().html.add_child(
        folium.Element(refresh_js)
    )
    m.get_root().html.add_child(
        folium.Element(refresh_button)
    )

    #reset view to original map view button
    reset_js, reset_button = map_utils.resetter(home_lat=home_lat,home_lon = home_lon,home_zoom = home_zoom)
    m.get_root().html.add_child(
        folium.Element(reset_js)
    )
    m.get_root().html.add_child(
        folium.Element(reset_button)
    )


    #get map name and attach to window
    map_name_js = map_name()
    m.get_root().script.add_child(
        folium.Element(map_name_js)
    )
    #put key for may amd layer states into the page
    map_keys_js = map_utils.state_keys()
    m.get_root().script.add_child(
        folium.Element(map_keys_js)
    )

    ##map state storage
    map_states_js = map_utils.save_map_states()
    m.get_root().script.add_child(
        folium.Element(map_states_js)
    )

    info_list = []
    id_list = []
    ### add the forecast markers

    for f in forecasts:
        print(f)
        info = parse_forecast_dir(f)
        if info is None:
            print('Forecast failed or has just started')
            continue
        ###changing files in the wksp directory will let older forecasts be put in the list
        if info["ign_utc"] + timedelta(hours = cutoff_hours) < map_time:
            continue
        #print(info)
        info_list.append(info)
        id_list.append(info['irwin_id'])
        url = build_forecast_url(info)
        end_time = compute_end_time(info)
        
        #forecast status    <<<<----------------------------------------------put in its own function
        is_active = end_time > map_time
        status = "ACTIVE" if is_active else "EXPIRED"
        if is_active and not info['forecast_complete']:
            status = "ACTIVE, INCOMPLETE FORESCAST"
        print(f'Forecast end time: {end_time}, forecast is {status}')
        color = "red" if is_active else "orange"
        group = active if is_active else expired
        icon = folium.Icon(color=color, icon="fire")  ### <<<<----------- needed?
        info['status'] = status

        #add satellites data to info
        inc_df = df[df.known_incident_id==info['irwin_id']].copy()
        feature_id = inc_df.feature_tracking_id.unique()
        #print(feature_id)
        for fi in feature_id:
            add_df = df[df.feature_tracking_id==fi].copy()
            inc_df = pd.concat([inc_df,add_df],ignore_index=True)
        inc_df = inc_df.drop_duplicates(keep='first')                    ### <<<<----- make last, to have most recent data on dashboard?
        if len(inc_df)>0:
            print(f'Adding {len(inc_df)} detections to the info["goes_geojson"] file')
            info['goes_geojson'] = make_geojson(inc_df)
        else:
            print('No data for',info['irwin_id'])
            '''
            info['goes_geojson'] = {
                "type":"FeatureCollection",
                "features" : [
                        {
                        "type" : "Feature",
                        "geometry" : {
                            "type" : "Point",
                            "coordinates" : []

                        },
                        "properties" : []
                    }
                ]
            }
            '''
        info['nifc_data'] = []  ###<<<< switch to dictionary   or remove entirely

        #adding NIFC information, if the locations file is found
        if len(nifc_csv) > 0:
            #add nifc perim and location data
            perim_element, loc_element = map_utils.add_nifc_data(info,nifc_csv=nifc_csv)
            #print(f'add debug {type(perim_element)} {type(loc_element)}')
            if perim_element:
                perim_element.add_to(nifcGroup)
            if loc_element:
                loc_element.add_to(nifcGroup)

        
        #scan the geojson file for the incident
        watch_feature = watch_duty.find_watchduty_incident(info,geojson)
        if watch_feature:
            geoElement = map_utils.watchduty_geoElement(watch_feature)
            geoElement.add_to(watchDuty)


    
        #add incident ignition marker
        popup,tooltip_html = map_utils.build_popup(info)
        big_mark = folium.Marker(
            location=[info["lat"], info["lon"]],
            popup=popup,
            zoom_on_click=True,
            tooltip=folium.Tooltip(
                info["incident_name"],
                permanent=False,
                direction="right",   # label appears to the right of marker
                offset=(8, 0),      # small spacing from marker
                style=f"""
                    background-color: rgba(255,255,255,0.4);
                    font-size: 10px;
                """
            ),
            icon=folium.Icon(color=color, icon="fire")
        )
        big_mark.add_to(group)
        
        geojson_poly = {'type': 'Feature',    ### <<<< ----------------------------- does this get used anywhere?
                    'geometry': info["geometry"],               
                  'properties': {"Perimeter": "Forecast Perimeter"}
        }

        #put two forecast perims on the map, last perim and the perim closest to cuurent time
        first_perim,last_perim = return_geojsons(info)
        if first_perim:
            print('Adding returned perims to map...')
            last_perimElement = map_utils.add_perim(perim=last_perim,perim_type='incident') 
            if last_perimElement:
                last_perimElement.add_to(group) 
            else:
                print(last_perim,info)######<<<<--------------------------------------------last_perim,frst_perim can be None if bad dictionary is passed

            first_perimElement = map_utils.add_perim(perim=first_perim,perim_type='incident') 
            if first_perimElement:
                first_perimElement.add_to(group)
            else:
                print(first_perim,info)
        else:
            print('No perims returned :(  ....')
        

        #save the info file and make webpage for the incident
        save_incident_info(info)
        #print('info keys:',info.keys())
        make_incident_webpage(info)
             
    ######## end forecast loop block
    #add nifc locations


    ### legend
    legend_html = map_utils.add_legend(map_time,cutoff_hours)
    m.get_root().html.add_child(folium.Element(legend_html))

    
    m.save("wildfire_map.html")
    map_utils.engr_scp("wildfire_map.html","~/NGFS/.",user='jhaley03')

    #cmd = 'scp  wildfire_map.html jhaley03@linux2.engr.colostate.edu:~/NGFS/.'
   # os.system(cmd)

    



