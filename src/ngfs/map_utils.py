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
import requests
import json
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
try:
   from pyproj import Transformer
except:
   print('pyroj Transformer is not available')


#### General helper function ##########
def fmt_time_string(dt): #convert datetime object
    return dt.strftime("%Y-%m-%d %H:%M UTC")

def compute_end_time(info):
    start = parse_time(info["forecast_start"])
    return start + timedelta(hours=info["duration_hr"])

def parse_time(time_str):
    return datetime.strptime(time_str, "%Y-%m-%d_%H:%M:%S")  #from job json files in WRFXPY

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
    try:
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
    except:
        print('Error getting perimeter area')
        return np.NaN

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
    msk = (fire_area > 0).astype(float)
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
    poly_area = perim_area(poly)
    smooth = poly.simplify(0.0001)
    if smooth.is_empty:
        return None, 0.0
    from shapely.geometry import mapping
    geom = mapping(smooth)    #shapely object -->> geojson dictionary
    if np.isnan(poly_area):
        poly_area = 0.0
    return geom, round(poly_area/(1000**2)*247.1)  ###acres

def get_elevation(wrfout_path):
    with nc4.Dataset(wrfout_path, 'r') as wrfout:
        fire_area = np.array(
            wrfout['FIRE_AREA'][0]
        )
        elevation = np.array(
            wrfout['ZSF'][0]
        )
    fire_msk = fire_area > 0
    fire_count = np.sum(fire_msk)
    if fire_count == 0:
        return None
    
    elev = elevation[fire_msk].mean()/1000.0
    print(f'Mean eleavtion in burned area is {elev}')
    return round(elev,3)
    


    

        

def make_fuels(wrfout_path):
    #returns a dictionary of fuels in the burnt area of the fire domain from a wrfout_file
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
        elevation = np.array(
            wrfout['ZSF'][0]
        )

    #fine where the firespread and its area
    fire_msk = fire_area > 0
    fire_count = np.sum(fire_msk)
    if fire_count == 0:
        return []

    #find the areas of fuels in burnt area
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

    #find the elevation
    
    return fuel_stats[:3]


############################### html, js, css utils  ######################

## js to use with save_map_states to save map properties like boxes ticked and vieing pan and zoom
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
        console.log("Map initialization complete");
    }, 500);

    });
    """
    return name_js

def state_keys(IrwinID = None):
    #assign unique state variables 
    if IrwinID:
        IrwinID = IrwinID.replace('{','').replace('}','')
        map_state = f'map_state_{IrwinID}'
        layer_state = f'layer_state_{IrwinID}'
    else:
        map_state = f'map_state'
        layer_state = f'layer_state'
    #js string to be injected via folium
    js = f"""
        const MAP_STATE_KEY = "{map_state}";
        const LAYER_STATE_KEY = "{layer_state}";
    """
    return js
    

   
def save_map_states():

    map_states_js = """
        function saveMapState(map) {
            console.log("saveMapState CALLED");
            const center = map.getCenter();
            localStorage.setItem(
                MAP_STATE_KEY,
                JSON.stringify({
                    lat: center.lat,
                    lon: center.lng,
                    zoom: map.getZoom()
                })
            );
        }
        function restoreMapState(map) {
            const saved = localStorage.getItem(MAP_STATE_KEY);
            if (saved) {
                const state = JSON.parse(saved);
                map.setView(
                    [state.lat, state.lon],
                    state.zoom
                );
            }
            console.log("Restored map state");
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
                LAYER_STATE_KEY,
                JSON.stringify(state)
            );
            console.log("Saved layer state:", state); 
        }

        function restoreLayerState(map) {
            const saved = localStorage.getItem(LAYER_STATE_KEY);
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

def ngfs_button_class(width = 160):
    button_string = f"""
            <style>
            .ngfs-button {{
                background-color: #2c7be5;
                color: white;
                border: black;
                padding: 6px 10px;
                border-radius: 6px;
                font-size: 10px;
                cursor: pointer;
                box-shadow: 0 2px 6px rgba(0,0,0,0.2);
                width: {width}px;
                text-align: center;
            }}

            .ngfs-button:hover {{
                background-color: #1b5fcc;
            }}
            </style>
            """
    return button_string
   
def refresher(refresh_minutes = 10,top=14,left = 60,z = 9999):
    refresh_js = """
        <script>
        let refreshEnabled =
            localStorage.getItem("autoRefresh") !== "false";
        let refreshTimer = null;
        function startRefresh() {
            stopRefresh();
            refreshTimer = setTimeout(() => {
                location.reload();
            }, 10 * 60 * 1000); // 10 minutes
        }
        function stopRefresh() {
            if (refreshTimer) {
                clearTimeout(refreshTimer);
                refreshTimer = null;
            }
        }
        function toggleRefresh() {
            refreshEnabled = !refreshEnabled;
            localStorage.setItem(
                "autoRefresh",
                refreshEnabled
            );
            updateRefreshButton();
            if (refreshEnabled) {
                startRefresh();
            } else {
                stopRefresh();
            }
        }
        function updateRefreshButton() {
            const btn = document.getElementById(
                "refresh-toggle"
            );
            btn.innerHTML = refreshEnabled
                ? "Auto-refresh: ON"
                : "Auto-refresh: OFF";
        }
        window.addEventListener("load", () => {
            updateRefreshButton();
            if (refreshEnabled) {
                startRefresh();
            }
        });
        </script>
    """
    refresh_button = f"""
        <div style="
            position: fixed;
            top: {top}px;
            left: {left}px;
            z-index: {z};
        ">
        <button id="refresh-toggle"
                onclick="toggleRefresh()"
                class="ngfs-button">
        </button>
        </div>
        """
    return refresh_js, refresh_button

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


def add_disclaimer(bottom = 10,left=60,max_width= 1100):
     ### disclaimer
    disclaimer_html = f"""
        <div style="
            position: fixed;
            bottom: {bottom}px;
            left: {left}px;
            z-index: 450;
            background-color: yellow;
            padding: 6px 10px;
            font-size: 8px;
            border-radius: 4px;
            box-shadow: 0 0 4px rgba(0,0,0,0.2);
            max-width: {max_width}px;
        ">
        <b>Disclaimer:</b> THIS DATA IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
          BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
            IN NO EVENT SHALL THE DATA PROVIDER(S) OR CONTRIBUTOR(S) BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
            EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
            LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
          WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
          ARISING IN ANY WAY OUT OF THE USE OF THIS DATA, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        </div>
        """
    return disclaimer_html

def add_legend(map_time,cutoff_hours,top = 80, left = 10):
    legend_html = f"""
    <div class="map-overlay" style="
        position: fixed;
        top: {top}px;
        left: {left}px;
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


##### Perim and location Elements #########################
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

def nifc_popup(geojson):
    incident_name = geojson.get('IncidentName',None)
    nifc_url = 'https://data-nifc.opendata.arcgis.com/'
    html = f"""
    <div class="fire-label" style="font-family: Arial; font-size: 13px; width: 240px;">
        <h4 style="margin-bottom:6px;">
            {incident_name}
        </h4>

        <b>IRWIN ID:</b><br>
        <div style="margin-bottom:6px;">
            {geojson.get("irwin_id", "N/A")}
        </div>

        <b>Incident Size:</b><br>
        <div style="margin-bottom:6px;">
            {geojson.get('IncidentSize','N/A')}
        </div>

        <b>FireDiscoveryDateTime:</b><br>
        <div style="margin-bottom:6px;">
            {geojson.get("FireDiscoveryDateTime", "N/A")}
        </div>

        <a href="{nifc_url}" target="_blank">
            NIFC Open Data Site
        </a><br>
    </div>
    """
    return html

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

def add_watchduty_loc(info,nifc_csv = None):
    #add location marker
    loc = nifc.get_nifc_incident(info["irwin_id"],feature_type='loc') #geojson object
    #print('found nifc location', loc)
    if loc and len(loc['features']) >0:
        if hasattr(loc,'FeatureCollection'): #a geojson was passed
            lon,lat = loc['features'][0]['geometry']['coordinates']
        else: # a feature was passed
            lon,lat = loc['geometry']['coordinates']
        loc_element = folium.Marker(
            location=[lat,lon],
            tooltip=f'WatchDuty ignition point {info["incident_name"]}',
            icon=folium.Icon(icon="w",prefix = "fa",color = "blue")  #<i class="fa-solid fa-n"></i>
        )
        return loc_element
    else:
        print(f'No Watch Duty data found for {info["irwin_id"]}')

def watchduty_geoElement(watch_feature):
    #takes a feature from GeoJson feature collection and returns a folium.GeoJson element
    wdg = {
                'type': 'FeatureCollection',
                'features': [watch_feature]
            }
    wd_popup = folium.GeoJsonPopup(
        fields=[
            "name",
            "is_active",
            "date_created",
            "date_modified",   
            "containment",
            "acreage",
            "url"
        ],
        aliases=[
            "Name",
            "Active",
            "Created",
            "Modified",
            "Containment",
            "Acreage",
            "Watch Duty"
        ],
        localize=True,
        labels=True,
        style="""
            background-color: white;
            border: 1px solid black;
            border-radius: 4px;
        """
        )

    geoElement = folium.GeoJson(
        wdg,
        popup=wd_popup,
        tooltip=folium.GeoJsonTooltip(
            fields=["name"]
        )
    )
    return geoElement

def engr_scp(file,remote_dir,user = 'jhaley03'):
    import subprocess
    import time
    #example file = "wildfire_map.html", remote_dir = "~/NGFS/." "jhaley03@linux2.engr.colostate.edu:"

    servers = [f"linuxe{i}.engr.colostate.edu" for i in range(1,5)]
    servers.extend([f"linux{i}.engr.colostate.edu" for i in range(1,15)])

    for s in servers:
        cmd = [
            "scp",
            "-v",                      # Verbose output for debugging log files
            "-i", "~/.ssh/id_ed25519", # Explicitly forces your new Ed25519 key
            "-o", "BatchMode=yes",     # DONT prompt for passwords; fail immediately instead
            "-o", "ConnectTimeout=5",  # Drop connection if a server is o
            file,
            f"{user}@{s}:{remote_dir}"
        ]        
        print(cmd)
        result = subprocess.run(cmd)

        if result.returncode == 0:
            print("Transfer successful.")
            break

        print(f"Attempt {s} failed.")
        time.sleep(2)



def watchduty_popup(feature):   
    #find the feature with watchduty.find_watchduty_incident
    #feature = geojson['features'][0]     
    wd_url = feature['properties']['url']
    incident_name = feature['properties']['name']
    active = feature['properties']['is_active']
    containment = f'{feature["properties"]["containment"]}%'
    acreage = feature['properties']['acreage']
    html = f"""
    <div class="fire-label" style="font-family: Arial; font-size: 13px; width: 240px;">
        <h4 style="margin-bottom:6px;">
            {incident_name}
        </h4>

        <b>Active:</b><br>
        <div style="margin-bottom:6px;">
            {active}
        </div>

        <b>Containment:</b><br>
        <div style="margin-bottom:6px;">
            {containment}
        </div>

        <b>Acreage:</b><br>
        <div style="margin-bottom:6px;">
            {acreage}
        </div>

        <b>Created:</b><br>
        <div style="margin-bottom:6px;">
            {feature['properties']['date_created']}
        </div>

        <b>Modified:</b><br>
        <div style="margin-bottom:6px;">
            {feature['properties']['date_modified']}
        </div>

        <a href="{wd_url}" target="_blank">
            WatchDuty
        </a>
    </div>
    """
    return folium.Popup(html, max_width=300), html



def add_nifc_perim(info):
        ## add nifc perim
    perim_element = None
    perims = glob.glob(f'/data/jhaley/wrfxpy/ngfs/perims/*{info["irwin_id"]}*.geojson')
    if len(perims) > 0:
        info['nifc_data'].extend(perims[0])
        print('Found nifc perims for the incident',perims)
        perim_element = add_perim(perims[0],perim_type='nifc')
    else: #try to get a perim from the API
        perim = nifc.get_nifc_incident(info["irwin_id"],feature_type='perim')
        if perim.get('features',None) and len(perim['features']) > 0: #save the perim locally and and add to map
            '''
            save_name = f'/data/jhaley/wrfxpy/ngfs/perims/{info["incident_name"]}_{info["irwin_id"]}.geojson'.replace(' ','_')
            try:
                with open(save_name,'w') as geo_file:
                    json.dump(perim,save_name,indent=2)
            except:                                                                                                         ###<<<------------------- improve this block
                print('Error saving geojson', perim)
            print('save_name',save_name)
            info['nifc_data'] = nifc_data #add path to the info file
            '''
            perim_element = add_perim(perim,perim_type='nifc')
    return perim_element

def add_nifc_loc(info,nifc_csv = None):
    #add location marker
    loc_element = None
    loc = nifc.get_nifc_incident(info["irwin_id"],feature_type='loc') #geojson object
    if not loc:
        return loc_element
    print('found nifc location', loc)
    if loc.get('features') and len(loc['features']) >0:
        popup_html = nifc_popup(loc['features'][0]['properties'])
        #print('debug 783',popup_html)
        lon,lat = loc['features'][0]['geometry']['coordinates']
        loc_element = folium.Marker(
            popup = popup_html,
            location=[lat,lon],
            tooltip=f'NIFC ignition point {info["incident_name"]}',
            icon=folium.Icon(icon="n",prefix = "fa",color = "blue")  #<i class="fa-solid fa-n"></i>
        )
    else:
        print(f'No NIFC data found for {info["irwin_id"]}')
    return loc_element
    
def add_nifc_data(info,nifc_csv = None):
    #returns folium elelemts to add to maps
    feature_types = ["perim","loc"]
    for feature_type in feature_types:
        ## add nifc perim
        if feature_type == 'perim':
            perim_element = add_nifc_perim(info)
        else:
            loc_element = add_nifc_loc(info)
    return perim_element, loc_element

'''
def nifc_locs():    ### remove this, put in nifc.py module
    url = (
        "https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/"
        "WFIGS_Incident_Locations_Current/FeatureServer/replicafilescache/"
        "WFIGS_Incident_Locations_Current_-5929201850136171392.csv"
    )
    geo_url = (
        "https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/"
        "WFIGS_Incident_Locations_Current/FeatureServer/replicafilescache/"
        "WFIGS_Incident_Locations_Current_-698092179933823454.geojson"
    )
    f = '/data/jhaley/wrfxpy/ngfs/perims/WFIGS_Incident_Locations_Current_-5929201850136171392.csv'
    file_time = os.path.getmtime(f)
    #download once each 30 minutes
    now = time.time()
    if now-file_time > 1800:  ###alway dlownload?
        cmd = f'wget -O {f} {url}'
        os.system(cmd)

    #read file
    df = pd.read_csv(f)
    return df
#
'''

'''
def get_nifc_incident(irwin_id,feature_type='loc'):   ### remove this, put in nifc.py module
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
        #load latest nifc locatio data
        nifc_csv = nifc_locs()
        #search for id string in csv file
        inc = nifc_csv[nifc_csv['IrwinID'] == irwin_id]
        if len(inc) > 0:
            return nifc.nifc_geojson(inc)  #geojson type object
    
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
        return None
'''
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

def add_perim(perim,perim_type = 'incident'):
    #takes a geojson dictionary or path to one and returns folim object to add to a map
    if not isinstance(perim, (dict, str)):
        print(f'weird {perim_type} geojson was sent. Type is {type(perim)}')
        print(perim)
        return None
    if perim_type == 'incident':
        fields=["fire_name","timestamp","area_acres"]
        aliases=["Fire Name","Forecast Perimeter Time","Area [acres]"]
        fillcolor = "red"
    elif perim_type == 'nifc':
        fields=["poly_IncidentName","poly_CreateDate","poly_Acres_AutoCalc"]
        aliases=["Fire Name","NIFC Perimeter Time","Area [acres]"]
        fillcolor = "yellow"
    else:
        print('unkown_perim_type')
    perim_element = folium.GeoJson(
                perim,
                style_function=lambda x: {
                    "color": "darkred",
                    "weight": 2,
                    "fillColor": fillcolor,
                    "fillOpacity": 0.4
                 },
                zoom_on_click=True,
                tooltip=folium.GeoJsonTooltip(
                    fields=fields,
                    aliases=aliases
                )
            )
    return perim_element




if __name__ == "__main__":
   pass
