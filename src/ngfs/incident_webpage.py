#builds a html page for the incident, with forecast perimeters and more
from __future__ import absolute_import
from __future__ import print_function
import os, sys, glob
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
import folium
import folium.plugins
import fire_init
import utils
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import shapely
from ngfs import ngfs_api as api
from ngfs import ngfs_ftp as ftp
from ngfs import config_manager 
from shapely.geometry import Polygon
import pandas as pd
from shapely.geometry import Point, LineString, Polygon, MultiPolygon
from pyproj import Proj, transform
from map_locations import return_geojsons, compute_end_time, fmt_time_string
import map_locations as ml
import map_utils
import sat_webpage
import watch_duty
import nifc
try:
   from pyproj import Transformer
except:
   print('pyroj Transformer is not available')
from jinja2 import Template




def add_geojson(geojson_path):
    pass

def padded_geojson_bounds(geojson,margin = 0.3):
    #expand boundary around the geometry of a geojson
    a = np.array(geojson['features'][-1]['geometry']['coordinates'][0])
    min_lat = min(a[:,1])
    max_lat = max(a[:,1])
    min_lon = min(a[:,0])
    max_lon = max(a[:,0])
    center_lat = (min_lat+max_lat)/2.
    center_lon = (min_lon+max_lon)/2.
    pad_lat = (center_lat-min_lat)*margin
    pad_lon = (center_lon-min_lon)*margin
    return [[min_lat-pad_lat,min_lon-pad_lon],[max_lat+pad_lat,max_lon+pad_lon]]

def api_detections(info,unique = False):   ##<---------------------------------------------------------probably get rid of this and used cached data from ftp module

    #assemble parameters for api call
    start_time = utils.esmf_to_utc(info['forecast_start'])-timedelta(hours=2)
    #end_time = start_time + timedelta(hours=info['duration_hr'])
    try:
        bounds = padded_geojson_bounds(info['geojson'],margin=1.2)
        bbox = (bounds[0][1],bounds[0][0],bounds[1][1],bounds[1][0])
    except:
        lat = info['lat']
        lon = info['lon']
        pad = 0.15
        bbox = [lon-0.15,lat-0.15,lon+0.15,lat+0.15]
    params = api.duration_parameters(start_time=start_time,duration_hours=int(info['duration_hr']))
    params = api.add_ogc_bbox(params,bbox=bbox)
    print(params)

    #
    sat = ['GOES-18','GOES-19']
    df = pd.DataFrame()
    for s in sat:
        df = pd.concat([df,api.make_ngfs_dataframe(sat=s,params = params,parse_t=False)])
        df.reset_index(inplace=True,drop=True)
    if unique: #first will keep the last time,
        df = df.drop_duplicates(subset=['latitude','longitude'],keep='last')

    print(f'API data type: {type(df)}, length: {len(df)}')
    df.reset_index(drop=True,inplace=True)
    return df

def make_geojson(df):
    #covert the time columns to something OK for geojson
    # Identify all datetime columns (including timezone-aware ones)
    date_cols = df.select_dtypes(include=['datetime64', 'datetimetz']).columns
    # Convert those columns to string
    df[date_cols] = df[date_cols].astype(str)
    features = []
    for i in df.index:
        geometry = {
            'type' : 'Point',
            'coordinates': [df['longitude'].loc[i],df['latitude'].loc[i]]
        }
        #print(geometry)

        prop = df.loc[i].to_dict()
        #prop = {"property":"test"}
        f = {
            "type" : "Feature",
            "geometry" : geometry,
            "properties" : prop
        }
        #print(i,type(f['geometry']['coordinates']))
        if isinstance(f['geometry']['coordinates'],list):
            features.append(f)
        else:
            print(df.loc[i])
    return {"type":"FeatureCollection","features":features}


def make_info_panel(info):

    #popup,tooltip_html = build_popup(info) 
    #tooltip_html = tooltip_html.replace('Open incident details page','Refresh page').replace('target="_blank"','')
    incident_name = (info.get("incident_name") or "").replace("_", " ")
    if 'fuels' in info.keys():
        fuels = info['fuels']
    else:
        fuels = 'Not computed yet'
    html = f"""
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

        <b>Fire Size [Acres]:</b><br>
        <div style="margin-bottom:6px;">
            {info["final_area"]}
        </div>

        <b>Primary fuels in burnt area:</b><br>
        <div style="margin-bottom:6px;">
            {fuels}
        </div>


    </div>
    """
    return html
'''
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
'''
def dynamic_url(info):
    job_id = info["name"]
    rasters = '&rasters=WINDVEC,FIRE_AREA,PM25_SFC_D'
    base_url = f"https://www.engr.colostate.edu/~jhaley03/wrfxweb/?job_id={job_id}{rasters}"
    dynamic_js = f"""
        function makeVisURL() {{
            const map = window.ngfsMap;
            const center = map.getCenter();
            const lat = center.lat.toFixed(4);
            const lon = center.lng.toFixed(4);
            const zoom = Math.max(map.getZoom(),3);
            console.log("FOUND LAT:", lat);
            console.log("FOUND LON:", lon);
            console.log("FOUND ZOOM:", zoom);
            return (
                "{base_url}" +
                `&zoom=${{zoom}}` +
                `&pan=${{lat}},${{lon}}`
            );
        }}

        function openVisualization() {{
            const url = makeVisURL();
            console.log("FOUND URL:", url);
            window.open(url, "_blank");
        }}
    """
    button_html = """
        <div style="
            position: fixed;
            top: 42px;
            left: 60px;
            z-index: 9999;
        ">
        <button class="ngfs-button"
                onclick="openVisualization()">
            Open Forecast Visualization
        </button>
        </div>
    """
    return dynamic_js,button_html

    
def make_map(info,split_screen=False):

    output_path = f'{info["wksp_path"]}/{info["grid_code"]}.html'
    bounds = None
    #check if geojson file has perimeters
    try:
        bounds = padded_geojson_bounds(info['geojson']) #geojson with all perims for wrfouts
    except: #make a geojson on the fly
        prop = {
            "fire_name" : info["incident_name"],
            "irwin_id" : info["irwin_id"],
            "name" : f'{info["incident_name"]} last perim',
            "timestamp": "forecast end time",
            "perimeter_type": "forecast_timestep",
            "area_acres" : info['final_area']
        }
        try:
            #individual perim features look like
                # CANYON_2026-05-07_18_00_00_64EABB0F-3EA8-47F3-BC7F-19750DF7E579_2026-05-08_03:00_UTC.geojson
            g = glob.glob(f'{info["wksp_path"]}/*UTC.geojson')
            features = []
            if len(g) > 0:
                for perim in g:
                    try:
                        with open(perim,'r') as file:
                            feat = json.load(file)
                        features.append(feat)
                        print(len(features))
                        print(feat['type'])
                    except:
                        print('Error reading ',perim)
                geo = {
                    "type": "FeatureCollection",
                    "features" :features
                }
                print(len(geo['features']))
                info['geojson'] = geo
                bounds = padded_geojson_bounds(info['geojson'])   
        except:
            try: #make a perimeter set of length 1 from last wrout
                feat = {
                    "type":"feature",
                    "geometry":info["geometry"],
                    "properties":prop
                }
                geo = {
                    "type": "FeatureCollection",
                    "features" :feat
                }
                info['geojson'] = geo
                bounds = padded_geojson_bounds(info['geojson'])    
            except:
                print('geojson file has no perimeters')
                html = f"""<!DOCTYPE html>
                            <html>
                            <head>
                                <title>{info['incident_name']}</title>
                                <meta http-equiv="refresh" content="600">
                            </head>
                            <body>
                                <h1>Waiting for forecast output...</h1><br>
                            </body>
                """
                with open(output_path,"w", encoding="utf-8") as file:
                    file.write(html)
                #push to the server
                map_utils.engr_scp(output_path,'~/NGFS/incidents/.')
                #cmd = f'scp {output_path} jhaley03@linux1.engr.colostate.edu:~/NGFS/incidents/.'
                #os.system(cmd)
                return

    #create map
    m = folium.Map(
    location=[info['lat'],info['lon']],
    tiles = None,
    min_zoom = 3,
    zoom_control=True
    )
    if bounds:
        m.fit_bounds(bounds)
  

    
    #mapping tile layers from https://leaflet-extras.github.io/leaflet-providers/preview/
    #folium.Map(tiles='https://{s}.tiles.example.com/{z}/{x}/{y}.png', attr='My Data Attribution')
    layer_image = folium.TileLayer(tiles= "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",               
	                attr = "Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community",
                    name = "Esri.WorldImagery",
                    control=not split_screen
    )
    
    layer_topo = folium.TileLayer(tiles =  "https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}.png",
                     attr = """Tiles &copy; Esri &mdash; Esri, DeLorme, NAVTEQ, TomTom, Intermap, iPC, USGS, FAO, NPS, NRCAN, 
                               GeoBase, Kadaster NL, Ordnance Survey, Esri Japan, METI, Esri China (Hong Kong), and the GIS User Communit""",
                     name = "Esri.WorldTopoMap"
    )
    
    layer_osm = folium.TileLayer('OpenStreetMap',control=not split_screen) #default layer     
    
    if split_screen:
        sbs = folium.plugins.SideBySideLayers(layer_left = layer_osm,layer_right=layer_image)
        layer_osm.add_to(m)
        layer_image.add_to(m)
        sbs.add_to(m)
    else:
        layer_osm.add_to(m)    
        layer_image.add_to(m)
        layer_topo.add_to(m)   

    #button
    m.get_root().html.add_child(
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
    
    ##map state storage
    #put key for map amd layer states into the page
    map_keys_js = map_utils.state_keys(IrwinID=info['irwin_id'])     ##<<----------------------------------put IrwinID in the function
    m.get_root().script.add_child(
        folium.Element(map_keys_js)
    )

    ##map state storage
    map_states_js = map_utils.save_map_states()
    m.get_root().script.add_child(
        folium.Element(map_states_js)
    )
    
    #get map name and attach to window
    map_name_js = ml.map_name()       
    m.get_root().script.add_child(
        folium.Element(map_name_js)
    )
    

    #make url and put button for forecast visualization
    button_js,button_html = dynamic_url(info)
    m.get_root().script.add_child(
        folium.Element(button_js)
    )
    m.get_root().html.add_child(
        folium.Element(button_html)
    )
    
        
    #layer control group
    perims = folium.FeatureGroup("Forecast Perimeters").add_to(m)  
    detections = folium.FeatureGroup("GOES FIRE DETECTIONS").add_to(m)
    viirs = folium.FeatureGroup("VIIRS FIRE DETECTIONS").add_to(m)
    if 'nifc_data' in info.keys():
        nifc_group = folium.FeatureGroup("NIFC Data",show=False).add_to(m)
    watchDuty = folium.FeatureGroup("WatchDuty Data",show=False).add_to(m)   ### maybe add this to info as well

    folium.LayerControl().add_to(m) 

    ### add header
    m.get_root().header.add_child(folium.Element(f"""
    <title>{info['incident_name']} Forecast</title>
    """))
    #auto refre

    ### auto-refresh             ##<<<-------------------------possibly use refresher
    refresh_html = """
        <meta http-equiv="refresh" content="600">
    """
    m.get_root().header.add_child(folium.Element(refresh_html))

    #disclaimer txt
    disclaimer_html = map_utils.add_disclaimer()
    m.get_root().html.add_child(folium.Element(disclaimer_html))

    #add ignition point
    color = "red" if 'ACTIVE' in info['status'] else "orange"
    big_mark = folium.Marker(
            location=[info["lat"], info["lon"]],
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
    big_mark.add_to(m)

    #add information panel on left side  
    info_html = make_info_panel(info)
    m.get_root().html.add_child(folium.Element(info_html))


    #add perimeters
    first_perim,last_perim = return_geojsons(info)
    if first_perim:
        print('Adding returned perims to map...')
        folium.GeoJson(
            last_perim,
            name=f'{info["incident_name"]}_last_perim',
            style_function=lambda x: {
                "color": "darkred",
                "weight": 2,
                "fillColor": "red",
                "fillOpacity": 0.4
                },
            zoom_on_click=True,
            tooltip=folium.GeoJsonTooltip(
                fields=["fire_name","timestamp","area_acres"],
                aliases=["Fire Name","Forecast Perimeter Time","Area [acres]"]                            
            )
            #popup=folium.GeoJsonPopup('Testing popup function')
        ).add_to(perims)
        folium.GeoJson(
            first_perim,
            name=f'{info["incident_name"]}_first_perim',
                style_function=lambda x: {
            "color": "darkred",
            "weight": 2,
            "fillColor": "darkred",
            "fillOpacity": 0.4
            },
            zoom_on_click=True,
            #popup=folium.GeoJsonPopup('Testing popup function'),
            tooltip=folium.GeoJsonTooltip(
               fields=["fire_name","timestamp","area_acres"],
              aliases=["Fire Name","Forecast Perimeter Time","Area [acres]"]   
            )
        ).add_to(perims)
    else:
        print('No perims returned :(  ....')

    #add nifc data  #<<<<-------------------------------------------------------- get nifc location too
    #add nifc perim and location data
    
    perim_element, loc_element = map_utils.add_nifc_data(info)
    #print(f'add debug {type(perim_element)} {type(loc_element)}')
    if perim_element:
        perim_element.add_to(nifc_group)
    if loc_element:
        loc_element.add_to(nifc_group)
    


    #scan the geojson file for the incident, this will only load a saved incident
    watch_feature = watch_duty.find_watchduty_incident(info)
    if watch_feature:
        geoElement = map_utils.watchduty_geoElement(watch_feature)
        geoElement.add_to(watchDuty)
    
    
    ### area plotting
    areas,times = make_area_series(info)
    if True:#areas:
        print('Adding area panel')
        m.get_root().header.add_child(
            folium.Element(
                """
                <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
                """
            )
        )

        chart_data_js = f"""
        const fireTimes = {json.dumps(times)};
        const fireAreas = {json.dumps(areas)};
        """

        m.get_root().script.add_child(
            folium.Element(chart_data_js)
        )

        chart_panel = """
            <div id="fireChartPanel"
                style="
                    display:none;
                    position:fixed;
                    top:80px;
                    right:20px;
                    width:600px;
                    height:350px;
                    background:white;
                    border:1px solid #666;
                    z-index:9999;
                    padding:10px;
                ">

                <button onclick="
                    document.getElementById(
                        'fireChartPanel'
                    ).style.display='none';
                ">
                    Close
                </button>

                <canvas id="fireChart"></canvas>

            </div>
        """

        button_html = """
        <div style="
            position:fixed;
            top:80px;
            left:62px;
            z-index:9999;
        ">
            <button
                class="ngfs-button"
                onclick="showFireGrowth()">

                Fire Growth

            </button>
        </div>
        """

        chart_js = """
            let fireChart = null;

            function showFireGrowth() {

                document.getElementById(
                    "fireChartPanel"
                ).style.display = "block";

                if (fireChart) return;

                const ctx =
                    document.getElementById(
                        "fireChart"
                    );

                fireChart = new Chart(ctx, {

                    type: "line",

                    data: {

                        labels: fireTimes,

                        datasets: [{
                            label: "Burned Area (acres)",
                            data: fireAreas
                        }]
                    },

                    options: {

                        responsive: true,

                        maintainAspectRatio: false,

                        scales: {

                            y: {
                                title: {
                                    display: true,
                                    text: "Area (acres)"
                                }
                            },

                            x: {
                                title: {
                                    display: true,
                                    text: "Forecast Time"
                                }
                            }
                        }
                    }
                });
            }
            """
        
        m.get_root().html.add_child(
            folium.Element(chart_panel)
        )
        
        m.get_root().script.add_child(
            folium.Element(chart_js)
        )





    '''
    if 'nifc_data' in info.keys():
        if len(info['nifc_data']):
            folium.GeoJson(
                    info['nifc_data'][0],
                    style_function=lambda x: {
                        "color": "black",
                        "weight": 2,
                        "fillColor": "yellow",
                        "fillOpacity": 0.4
                    },
                    zoom_on_click=True,
                    tooltip=folium.GeoJsonTooltip(
                        fields=["poly_IncidentName","poly_PolygonDateTime","poly_Acres_AutoCalc"],
                        aliases=["Fire Name","NIFC Perimeter Time","Area [acres]"]
                    )
                ).add_to(nifc)
    '''

    
    #add satellite detections
    #api call to get all detections during the forecast period
    print('Adding detections')
    #compute time parameters
    end_time_padding = 24
    end_time = compute_end_time(info) + timedelta(hours=end_time_padding)
    duration_hr = info['duration_hr'] + end_time_padding
    print(f'Will add detection data from {str(end_time)} back to the previous the previous {duration_hr} hours')
    geo = info.get('goes_geojson',None)
    #print(geo)
    popup_keys = [
                'latitude',
                'longitude',
                'pixel_area',
                'acq_date_time',
                'satellite',
                'scan_domain',
                'confidence',
                'version',
                'bright_t13',
                'frp',
                'quality_flag',
                'type_description',
                'daynight',
                'state',
                'county',
                'gacc_id',
                'nws_region',
                'nws_wfo_code',
                'nws_wfo_name',
                'nws_fire_wx_code',
                'known_incident_name',
                'known_incident_id',
                'land_cover',
                'fuel',
                'feature_tracking_id',
                'feature_frp'
            ]
    if not geo:                              #  <------------  ##### use ftp .add_viirs if older
        print('Getting detections from API call.')
        df = api_detections(info,unique=False)
        print(info['name'],'has ',len(df),' detections')
        #geo = api.dataframe_to_geojson(df)   
        
        #print(geo.keys())
        if len(df) == 0:
            geo = None      
        else:
            geo = make_geojson(df)
            pk = set(popup_keys)
            pf = set(df.keys())
            popup_keys=list(pf.intersection(pk))                                           ## if this is empty
    types = ['Possible Wildland Fire', 'Known Wildland Fire Incident'] 
    if geo.get("features"): #False: #True:#len(geo)>0:
        print(geo)
        for t in types:
                    # # <<<<<------------------------------------------------------- should be subsetting dataaframe here??
            if 'Possible' in t:
                color = 'orange'
            else:
                color = 'red'
            folium.GeoJson(
                geo,
                marker=folium.Circle(radius=1000, fill_color=color, color="black", weight=0.2),
                tooltip=folium.GeoJsonTooltip(
                    fields=["satellite", "acq_date_time","known_incident_name"],
                    aliases=["Satellite","Time","Known Incident Name"]                            
                ),
                popup=folium.GeoJsonPopup(fields=popup_keys,max_width=600)
            ).add_to(detections)
    else:
        print('No detection json found')

    
    #add viirs detections
    #
    days_to_get=duration_hr/24.0
    end_time = pd.Timestamp.now('UTC')
    start_time = pd.Timestamp.now('UTC')-timedelta(hours=days_to_get*24)
    viirs_df = ftp.add_ngfs_scene(ngfs_cfg,parse_times = False,sat='viirs',start_time=start_time,end_time=end_time)
    '''
    viirs_df = ftp.add_ngfs_viirs(
        ngfs_cfg,
        parse_times=False,
        end_time=pd.Timestamp(end_time),
        days_to_get=duration_hr/24.0,
        sat = None
    )
    '''
    #filter
    print(bounds) #use later
    print(f'Loaded {len(viirs_df)} viirs detections')
    if len(viirs_df)>0:
        viirs_df = viirs_df[viirs_df['known_incident_id'] == info['irwin_id']]
        print(f'After filtering, there are {len(viirs_df)} viirs detections for {info["irwin_id"]}')
        for t in types:
            mtype = 'circle'
            if mtype == 'circle':
                geo = make_geojson(viirs_df[viirs_df['type_description']==t])
            else:
                geo = ftp.make_geojson(
                    viirs_df[viirs_df['type_description']==t],
                    geometry_type="polygon"
                )
            
            popup_keys = list(viirs_df.keys())
            '''
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
                    'feature_frp'
                ]
            '''
            if 'Possible' in t:
                color = 'orange'
                continue  ##don't add these right now, too many........... need to filter for area around fires
            else:
                color = 'red'

            #check for empyty gejon file
            if len(geo["features"]) == 0:
                continue
            if mtype == 'circle':
                folium.GeoJson(
                    geo,
                    marker=folium.Circle(radius=250, fill_color=color, color="black", weight=0.4),
                    tooltip=folium.GeoJsonTooltip(
                        fields=["satellite", "acq_date_time","known_incident_name"],
                        aliases=["Satellite","Time","Known Incident Name"]                            
                        ),
                    popup=folium.GeoJsonPopup(fields=popup_keys,max_width=600)
                ).add_to(viirs)
            else:
                folium.GeoJson(geo).add_to(viirs) ###for adding polgons
                    
                '''
                name='VIIRS detection',
                style_function=lambda x: {
                    "color": "black",
                    "weight": 2,
                    "fillColor": "green",
                    "fillOpacity": 0.2
                },
                tooltip=folium.GeoJsonTooltip(
                    fields=["satellite", "acq_date_time","known_incident_name"],
                    aliases=["Satellite","Time","Known Incident Name"]                            
                    ),
                popup=folium.GeoJsonPopup(fields=popup_keys,max_width=600)
                '''


    #testing some plugns
    folium.plugins.MousePosition().add_to(m)
    folium.plugins.MiniMap(toggle_display=True).add_to(m)
    folium.plugins.Draw(export=True,position="bottomright").add_to(m)  
    #from folium.plugins import measure_control
    #m.add_child(measure_control)  <<<<<----------------not working?

    m.save(output_path)
    '''
    try:
        m.save(output_path)
    except:
        print(f'Error saving {output_path}')
    '''

    #push to the server
    #cmd = f'scp {output_path} jhaley03@linux7.engr.colostate.edu:~/NGFS/incidents/.'
    #os.system(cmd)

    map_utils.engr_scp(output_path,'~/NGFS/incidents/.')
               
def make_area_series(info):
    geo_json = f'{info["wksp_path"]}/{info["grid_code"]}.geojson'
    if not os.path.exists(geo_json):
        geo_json = f'{info["wksp_path"]}/incomplete_perimeters.geojson'
    try:
        with open(geo_json,'r') as file:
                g = json.load(file)
    except:
        print('Error opening perimeters geojson file(s)')
        return
    times = []
    areas = []
    for feat in g['features']:
        times.append(feat['properties']['timestamp'])
        areas.append(feat['properties']['area_acres'])
    #print(areas,times)

    return areas,times


if __name__ == "__main__":

    #tempory, while thinking of better way to handle it
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


    #job_args = json.load(open(sys.argv[1]))
    #for n,a in enumerate(sys.argv):
    #  print(n,a)
    with open(sys.argv[1],'r') as info_file:
        info = json.load(info_file)
    
    info['map_time'] = datetime.utcnow()
    #fill_template(info)
    print(f'Building webpage for {info["incident_name"]}')
    if not info['status'] == 'EXPIRED_junk':     ###<<<------------------------------------------remove
        make_map(info,split_screen=False)
    else:
        print('Skipping expired forecast')

    #print satellite informatiom
    lat = info['lat']
    lon = info['lon']
    alt = info['elevation']
    try:
        sat_webpage.sat_overpass(lat,lon,alt,horizon = 15,length = 8)
    except:
        print('Error making satellite calculations')

    print("*******************************************************")
    print()
