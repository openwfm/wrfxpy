from __future__ import absolute_import
from __future__ import print_function
import sys, os
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
import json
import pandas as pd
from datetime import datetime, timedelta
from shapely.geometry import Polygon, Point

def set_red_flags(date_str):
    #returns list of state zones with active red_flag warnings for current the ngfs_day.datestr
    #   each element of the list is a dictionary with key 'geometry' that is shapely object  
    # all active warnings are downloaded from NWA api : https://api.weather.gov/alerts/active?
    #rf_shape is geojson file with warning regions derived from shape file from :
    #  https://www.weather.gov/source/gis/Shapefiles/WSOM/fz05mr24.zip
    # this shape file needs to be updated , see : https://www.weather.gov/gis/
    #convert the time strings in the NWS warnings to UCTC with --> pd.Timestamp.astimezone(pd.Timestamp(ef.properties[2]['effective']),tz='UTC')
    rf_shape = pd.read_json('ingest/red_flag/west_fire_warning_zones.geojson')
    #make list of the zones covered, some WFO may issue red flag warnings for zones not in this western subset
    state_zones = list()
    for i in rf_shape['features']:
        state_zones.append(i['properties']['STATE_ZONE'])
    #rf_xml = 'ingest/red_flag/current_red_flag_'+self.date_str+'.xml'
    #os.system('wget -O ' + rf_xml +' https://api.weather.gov/alerts/urn:oid:2.49.0.1.840.0.7761e4c3d072b6e56fc4dca8513cf8b02487b2ae.001.1.cap ')
    active_json = 'ingest/red_flag/active_'+date_str+'.geojson'
    os.system('wget -O ' + active_json + ' https://api.weather.gov/alerts/active?')
    red_flag_warnings = list()
    redflags = list()
    effective = list()
    ends = list()
    #try:
    active_frame = json.load(open(active_json))
    for f in active_frame['features']:
        if 'Red Flag' in f['properties']['event']:
            red_flag_warnings.append(f)
            for zones in f['properties']['geocode']['UGC']:
                #remove the 'Z' to match codes in shape file
                if zones.replace('Z','') in state_zones:
                    redflags.append(zones.replace('Z',''))
                    effective.append(pd.Timestamp.astimezone(pd.Timestamp(f['properties']['effective']),tz='UTC'))
                    ends.append(pd.Timestamp.astimezone(pd.Timestamp(f['properties']['ends']),tz='UTC'))
    print(len(redflags), ' red flag zones found: ',redflags)
    #except:
           # print('Error reading NWS warning file')
    rf_zones = list()
    for rf in redflags:
        for f in rf_shape['features']:
            if rf == f['properties']['STATE_ZONE']:
                f['geom'] = f['geometry']
                try:
                    f['geometry'] = Polygon(f['geom']['coordinates'][0])
                    rf_zones.append(f)
                except:
                    #print(f['geom'])
                    print('Skipping red flag zone, problem reading object geometry')
                break

    try:
        red_flag_zones = pd.DataFrame(rf_zones)
        red_flag_zones['effective'] = effective
        red_flag_zones['ends'] = ends
    except:
        red_flag_zones = pd.DataFrame()
        print('Error combining red flags zone, returning empty DataFrame')
    '''
    for i in range(len(df.data)):
        viirs_point = Point(df.data.lon_tc.iloc[i],df.data.lat_tc.iloc[i])
        for j in range(len(active_red_flags)):
            if viirs_point.within(active_red_flags.iloc[j].geometry):
                print(df.data.incident_name.iloc[i],active_red_flags.STATE_ZONE.iloc[j])

    '''
    
    return red_flag_warnings, red_flag_zones

def red_flag_incident(lon,lat,ign_utc,red_flag_zones):
    #returns a boolen for whether the point lon,lat is inside a redflag warning area within 
    #red_flag zones is pd.DataFrame with Polygon geometry and timestamp effective and ends columns for each zone
    end_utc = ign_utc + timedelta(hours = 24)
    p = Point(lon,lat)
    for i in range(len(red_flag_zones)):
        rfz = red_flag_zones.iloc[i]
        if not (end_utc < rfz['effective'] or ign_utc > rfz['ends'] ) and p.within(rfz['geometry']):
            #print('The point (',lat,lon,') is within red flag warning zone ',rfz['properties']['NAME'],rfz['properties']['STATE_ZONE'])
            return True
    return False

if __name__ == '__main__':
    pass