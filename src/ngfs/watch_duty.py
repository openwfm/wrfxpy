#tools for getting information from watchduty.org
import folium
import html
import json, requests, os, sys
import time
from datetime import datetime, timezone
#import map_utils
import cloudscraper #pip install cloudscraper
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
sys.path.insert(1, 'src/ngfs')
from ngfs import nifc



url = "https://api.watchduty.org/api/v1/geo_events/?geo_event_types=*"
link = "https://api.watchduty.org/i/"


#### functions for downloading data ######

def watchduty_scraper():
    #gets recents data from watchduty
    #returns a list of watchduty ditionaries
    scraper = cloudscraper.create_scraper()
    url = (
        "https://api.watchduty.org/"
        "api/v1/geo_events/"
        "?geo_event_types=*"
    )
    r = scraper.get(url)
    r.raise_for_status()
    data = r.json()  
    return data

def watchduty_geojson(watchduty_scrape):
    #goes through the scrape of the watchduty site and makes a geojson from the dictionaries in the list
    features = []
    for f in watchduty_scrape:
        features.append(watchduty_feature(f))
    g = {
        "type" : "FeatureCollection",
        "features" : features,
        "properties" : {
            "created" : time.time()
        }
    }
    return g

def download_watchduty():


    #downloads latest watchduty information, converts to geojon and saves it
    watch_dict = watchduty_scraper()
    geojson = watchduty_geojson(watch_dict)
    save_watchduty(geojson)

def load_watchduty():
    #loads the geojson or downloads fresher data and 
    f = 'ngfs/perims/watch_duty.geojson'
    file_time = os.path.getmtime(f)
    #download once each 30 minutes
    now = time.time()
    if now-file_time > 1800:  ###always dlownload?
        download_watchduty()
    with open(f,'r') as file:
        geojson = json.load(file)
    return geojson

def save_watchduty(geojson):
    f = 'ngfs/perims/watch_duty.geojson'
    with open(f,'w') as file:
        json.dump(geojson,file,indent=2)

#### functions for working with incidents ######
def find_watchduty_incident(info,geojson=None):
    #scans the features in watchduty geojson to find irwin_id of forecast info
    #info can be string or dictionary
    #return a geojson feature
    if isinstance(info,str):
        irwin_id = info
        watchduty_file = None
    else:
        irwin_id = info['irwin_id']
        #try to see if something has been saved before
        watchduty_file = info['info_file'].replace('.json','_watchdutyFeature.json')
        try:
            if os.path.exists(watchduty_file):
                with open(watchduty_file,'r') as file:
                    f = json.load(file)
                    return f
        except:
            print('Error reading file')
    
    #search for id in watchduty geojson file
    if not geojson: #don't search, just return None
        return None
    found_wd = None
    for f in geojson['features']:
        if f['properties']['external_id'] == irwin_id:
             found_wd = f
             break
             #loc.get('feautures') and len(loc['features']) >0
    #if search fails, get the nifc location file and try to match lat,lon with watchduty location
    if not found_wd:
        nifc_loc = nifc.get_nifc_incident(irwin_id,feature_type='loc')   #### <<<----- read thois from local storage, or store nifc coods in info
        if nifc_loc.get('feautures') and len(nifc_loc['features']) >0:
            lon,lat = nifc_loc['features'][0]['geometry']['coordinates']
            for f in geojson['features']:
                lon2,lat2 = f['geometry']['coordinates']
                if (abs(lon-lon2) + abs(lat-lat2)) < 0.05: 
                    found_wd = f
                    break
    #add other strategy like checking name
    if found_wd:
        #save this locally
        with open(watchduty_file,'w') as file:
                json.dump(found_wd,file,indent=2)
        return found_wd

    return None

def watchduty_feature(watch_dict):
    #constructs a geojson feature from  the watchduty scrape
    geo = {
            'type' : 'Point',
            'coordinates': [watch_dict['lng'],watch_dict['lat']]
        }
    wd_url = f"https://app.watchduty.org/i/{watch_dict['id']}"
    prop = {
        "name" : watch_dict.get('name',None),
        "id": watch_dict.get("id",None),                       #### maybe use the whole dictionary here
        "is_active" : watch_dict.get("is_active",None),
        'date_created' : watch_dict.get('date_created',None),
        'date_modified': watch_dict.get('date_modified',None),
        'containment' : watch_dict['data'].get('containment',None),
        'acreage' : watch_dict['data'].get('acreage',None),
        'external_id' : watch_dict.get('external_id',None),
        'url' : wd_url
         }
         
    f = {
            "type" : "Feature",
            "geometry" : geo,
            "properties" : prop
        }
    
    return f




def get_watchduty_reports(geo_event_id):
    url = (
        "https://api.watchduty.org/"
        "api/v1/reports/"
    )
    params = {
        "geo_event_id": geo_event_id
    }
    scraper = cloudscraper.create_scraper()
    r = scraper.get(
        url,
        params=params,
        timeout=30
    )
    r.raise_for_status()
    return r.json()


def make_wd_marker(geojson):
    properties = geojson['features'][0]['properties']
    rows = []
    for key, value in properties.items():
        # clickable link for URLs
        if key == "url":
            value_html = (
                f'<a href="{value}" '
                f'target="_blank">Watch Duty</a>'
            )
        else:
            value_html = html.escape(str(value))
        rows.append(
            f"<tr>"
            f"<td><b>{html.escape(key)}</b></td>"
            f"<td>{value_html}</td>"
            f"</tr>"
        )
    table = (
        "<table style='width:100%'>"
        + "".join(rows) +
        "</table>"
    )
    popup = folium.Popup(
        table,
        max_width=350
    )
    for feature in geojson["features"]:
        coords = feature["geometry"]["coordinates"]
        lon, lat = coords
        props = feature["properties"]
        popup = popup
        
    marker = folium.Marker(
            location=[lat, lon],
            popup=popup,
            tooltip=props.get("name", "incident")
    )

    return marker


'''
def make_wd_map(geojson):

    m = folium.Map(
        location=[39, -105],
        zoom_start=5
    )
    for feature in geojson["features"]:

        coords = feature["geometry"]["coordinates"]

        lon, lat = coords

        props = feature["properties"]

        popup = make_popup(props)

        folium.Marker(
            location=[lat, lon],
            popup=popup,
            tooltip=props.get("name", "incident")
        ).add_to(m)
'''
     

if __name__ == "__main__":
    pass
