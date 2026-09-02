#maps the locations of the JPSS satellites
import folium
import folium.plugins
from urllib.request import urlopen
from pathlib import Path
from datetime import datetime, timedelta
import numpy as np
import os
import pyorbital
from pyorbital.orbital import Orbital
from map_utils import fmt_time_string

def read_tle_files():
    #reads the tles files in  and
    for s,t in sat_dict.items():
        orb = Orbital(sat_dict[s]['name'],tle_file=sat_dict[s]['tle_file'])
        sat_dict[s]['orb'] = orb

def satellite_tle():
    # update tle files twice a day
    # https://celestrak.org/NORAD/elements/gp.php?CATNR=25544&FORMAT=TLE
    now = datetime.utcnow()
    now_str = now.isoformat()
    if ((now.hour == 0) or (now.hour == 12)) and now.mimute < 5:
        for sat_name,sat in sat_dict.items():
            url = f'https://celestrak.org/NORAD/elements/gp.php?NAME={sat["url_name"]}&FORMAT=TLE'
            output_file = Path(sat["tle_file"])
            with urlopen(url) as response:
                text = response.read().decode("utf-8")
            output_file.write_text(text)  
            print(f"Saved {sat_name} -> {output_file}")
            #COPY THE FILE TO SOMETHING WITH TIMESTAMP FOR ARCHIVING
            copy_path = Path(f'{sat["tle_file"]}_{now_str}')
            cp_cmd = f'cp {output_file} {copy_path}'
            os.system(cp_cmd)

def sat_geo_json(sat,hours=0.5,rate =12.0):
    #returns a geojson file wil locations of the satelite for the next hours
    #rate is the positions per hour to record
    now = datetime.utcnow()
    orb = sat['orb']

    features = []
    for i in range(int(hours*rate)+1):
        t = now + timedelta(hours = i / rate)
        lon,lat,alt = orb.get_lonlatalt(t)
        geometry = {
            'type' : 'Point',
            'coordinates': [lon,lat]
        }
        prop = {
                "name" : sat['name'],
                "time" : t.isoformat()
        }
        feat = {
            "type" : "Feature",
            "geometry" : geometry,
            "properties" : prop
        }
        features.append(feat)
    return {
        "type":"FeatureCollection",
        "features":features, 
        "properties" : {"name" : sat['name']}
    }

def line_from_points(points_geo):
    print('Making line geo')
    c = []
    #turn points into line
    for f in points_geo['features']:
        c.append(f['geometry']['coordinates'])
    geometry = {
        "type" : "LineString",
        "coordinates" : c
    }
    prop = { "name" : points_geo['properties']['name']}
    features= {
        "type" : "Feature",
        "geometry" : geometry
    }

    return {
        "type":"FeatureCollection",
        "features":[features], 
        "properties" : prop
    }


def add_legend(map_time):
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
    <h3><b>VIIRS Satellite Tracker</b></h3><br>
    <b>Information</b><br>
    Displaying current locations and paths over the next thrity minutes for the VIIRS satellites.
    The circles have radii of 1500 km, giving a sense of what the instruments onboard will be able to observe. 
       <br>
    Map Created: {fmt_time_string(map_time)}<br>  
    Webpage auto-refreshes every 5 minutes<br>
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
    <a href="https://www.engr.colostate.edu/~jhaley03/NGFS/wildfire_map.html" 
        target="_blank"
        rel="noopener noreferrer"
        >Wildfire Forecast Monitor</a><br>
    </div>
    
    """
    return legend_html

import requests

def get_altitude(lat, lon):
    # Free API endpoint (Open-Elevation)
    url = f"https://open-elevation.com{lat},{lon}"
    
    try:
        response = requests.get(url)
        # Check if the request was successful
        if response.status_code == 200:
            data = response.json()
            # Extract elevation from the results list
            elevation = data['results'][0]['elevation']
            return elevation
        else:
            return f"Error: {response.status_code}"
    except Exception as e:
        return f"An error occurred: {e}"

    # Example: Mount Everest coordinates
    '''
    latitude = 27.9881
    longitude = 86.9250
    altitude = get_altitude(latitude, longitude)

    print(f"The altitude at ({latitude}, {longitude}) is {altitude} meters.")
    '''

def sat_overpass(lat,lon,alt = 500,horizon = 15,length = 8):
    #dictionary with VIIRS satellite information and TLEs
    now_utc = datetime.utcnow()
    starts = []
    for s,t in sat_dict.items():
        satellite_name = sat_dict[s]['name']
        orb = Orbital(satellite_name,tle_file=sat_dict[s]['tle_file'])
        sat_dict[s]['orb'] = orb

        #passes for satellite
        passes = orb.get_next_passes(now_utc, lon=lon, lat=lat, alt=alt, horizon=horizon, length = length)

        print(f"Next passes for {satellite_name}:\n")
        for i, overpass in enumerate(passes):
            rise_time, fall_time, max_elev_time  = overpass
            print(f"Pass {i+1}:")
            print(f"\tStart (Rise): {rise_time.strftime('%Y-%m-%d %H:%M:%S')} UTC")
            print(f"\tPeak (Max Elevation): {max_elev_time.strftime('%Y-%m-%d %H:%M:%S')} UTC")
            print(f"\tEnd (Fall): {fall_time.strftime('%Y-%m-%d %H:%M:%S')} UTC")
            view_duration = (fall_time - rise_time).total_seconds()/60.0
            wait_time = (rise_time - now_utc).total_seconds()/60.0
            print(f'\tFire will be visible to {satellite_name} for {round(view_duration,2)} minutes in about {round(wait_time,2)} minutes')
            if i == 0:
                starts.append((satellite_name,wait_time, rise_time,view_duration))
    
    print(starts)

sat_dict = {
        "noaa_20" : {
            "color" : "red",
            "name" : "NOAA 20 (JPSS-1)",
            "url_name" : "NOAA%2020",
            "tle_file" : "/data/jhaley/wrfxpy/ingest/sat_data/noaa_20.tle"
        },
        "noaa_21" : {
            "color" : "green",
            "name" : "NOAA 21 (JPSS-2)",
            "url_name" : "NOAA%2021",
            "tle_file" : "/data/jhaley/wrfxpy/ingest/sat_data/noaa_21.tle"
        },
        "suomi" : {
            "color" : "blue",
            "name" : "SUOMI NPP",
            "url_name" : "SUOMI%20NPP",
            "tle_file" : "/data/jhaley/wrfxpy/ingest/sat_data/suomi.tle"
        }
    }



if __name__ == "__main__":

    #basic dictionary for datellites
    sat_dict = {
        "noaa_20" : {
            "color" : "red",
            "name" : "NOAA 20 (JPSS-1)",
            "url_name" : "NOAA%2020",
            "tle_file" : "/data/jhaley/wrfxpy/ingest/sat_data/noaa_20.tle"
        },
        "noaa_21" : {
            "color" : "green",
            "name" : "NOAA 21 (JPSS-2)",
            "url_name" : "NOAA%2021",
            "tle_file" : "/data/jhaley/wrfxpy/ingest/sat_data/noaa_21.tle"
        },
        "suomi" : {
            "color" : "blue",
            "name" : "SUOMI NPP",
            "url_name" : "SUOMI%20NPP",
            "tle_file" : "/data/jhaley/wrfxpy/ingest/sat_data/suomi.tle"
        }
    }

    #update the TLE files
    satellite_tle()

    #read the tle file and make Orbital objects
    read_tle_files()

    #make geojson for each each satellite
    geo_jsons = []
    for s,sat in sat_dict.items():
        geo_jsons.append(sat_geo_json(sat))

    map_time = datetime.utcnow()

    #make map
    m = folium.Map(location=[0,-40],zoom_start=2)
    #folium.plugins.Terminator().add_to(m)

    #autorefresh
    refresh_time = 5*60
    refresh_tag = f'<meta http-equiv="refresh" content="{refresh_time}">'

    # Inject meta tag into the map's root HTML
    m.get_root().html.add_child(folium.Element(refresh_tag))

    #add geojsons to map
    colors = ['red','green','blue']
    for n,sat_geo in enumerate(geo_jsons):
        c = colors[n]
        name = f"{sat_geo['features'][0]['properties']['name']} Track"
        folium.GeoJson(
            sat_geo,
            name = name,
            marker=folium.Circle(radius=1500*1000,color=c,fill=True),
            tooltip=folium.GeoJsonTooltip(fields=["name", "time"])
        ).add_to(m)

        #make line
        #path_line = line_from_points(sat_geo)
        #folium.GeoJson(path_line,name=name).add_to(m)
        
    '''
    folium.GeoJson(
    data,
    marker=folium.CircleMarker(radius=10, color="red", fill=True)
).add_to(m)

    radius = 10000
    folium.Circle(
        location=[-27.551667, -48.478889],
        radius=radius,
        color="black",
        weight=1,
        fill_opacity=0.6,
        opacity=1,
        fill_color="green",
        fill=False,  # gets overridden by fill_color
        popup="{} meters".format(radius),
        tooltip="I am in meters",
    ).add_to(m)
    '''
    legend_html = add_legend(map_time)
    m.get_root().html.add_child(folium.Element(legend_html))

    #save the html
    output_path = 'viirs_tracks.html'
    m.save(output_path)

    cmd = f'scp {output_path} jhaley03@linuxe2.engr.colostate.edu:~/NGFS/.'
    os.system(cmd)
