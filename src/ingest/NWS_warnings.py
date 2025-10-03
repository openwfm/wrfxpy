from __future__ import absolute_import, print_function
import sys
import os
import json
import time
from datetime import datetime, timedelta

import pandas as pd
from shapely.geometry import Polygon, Point, MultiPolygon

sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')


def subset_wfo_data(full_data, data, wfo_list):
   """
   Finds possible wildfire detections in a NWS WFO region withouth NIFC designation and assigns 
   incident_id_string and incident_name derived fron the feature_tracking_id
   """
   temp_data = full_data[full_data.type_description == 'Possible Wildland Fire']
   
   # Filter data by wfo_list and concatenate results
   wfo_data = pd.concat([temp_data[temp_data.wfo_code == w] for w in wfo_list], ignore_index=True)
   print(f"Found {len(wfo_data)} possible wildfire detections in {', '.join(wfo_list)}")

   # Remove detections already associated with known incidents
   known_features = data.feature_tracking_id.unique()
   wfo_data = wfo_data[~wfo_data.feature_tracking_id.isin(known_features)]

   # Assign incident_id_string and incident_name
   def generate_incident_info(row):
      fid = row.feature_tracking_id.replace(':', '').replace('-', '').replace('_', '')
      wfo_code = row.wfo_code
      row.incident_id_string = f'{{{wfo_code}{fid[-4:]}-{fid[2:6]}-{fid[6:10]}-{fid[11:15]}-{fid[-12:]}}}'
      row.incident_name = f"{wfo_code}_{fid[2:6]}_{fid[6:10]}_{fid[-12:]}"
      return row

   wfo_data = wfo_data.apply(generate_incident_info, axis=1)

   return wfo_data


def set_red_flags(date_str):
    """
    Returns a list of state zones with active red_flag warnings for the given date.
    Each element of the list is a dictionary with a key 'geometry' that is a shapely object.
    """

    #shape file is subset of data in https://www.weather.gov/source/gis/Shapefiles/WSOM/fz05mr24.zip
    #  this needs to be updated periodically, see:  https://www.weather.gov/gis/
    rf_shape = pd.read_json('ingest/red_flag/west_fire_warning_zones.geojson')

    # List of state zones to consider
    state_zones = [i['properties']['STATE_ZONE'] for i in rf_shape['features']
                   if 'AK' not in i['properties']['STATE_ZONE'] and 'HI' not in i['properties']['STATE_ZONE']]

    active_json = f'ingest/red_flag/active_{date_str}.geojson'

    # Check if the active_json file is archived or outdated, and download if needed
    is_archived = os.path.exists(active_json) and (time.time() - os.path.getmtime(active_json)) / 3600 > 72
    is_outdated = not os.path.exists(active_json) or (time.time() - os.path.getmtime(active_json)) / 3600 > 1

    if is_outdated and not is_archived:
        os.system(f'wget -O {active_json} https://api.weather.gov/alerts/active?')
    
    print(f'\tLoading {active_json} active warning file')
    active_frame = json.load(open(active_json))

    # Find zones with active red flag warnings
    red_flag_warnings, redflags, effective, ends = [], [], [], []
    for feature in active_frame['features']:
        if 'Red Flag' in feature['properties']['event']:
            red_flag_warnings.append(feature)
            for zone in feature['properties']['geocode']['UGC']:
                zone = zone.replace('Z', '')
                if zone in state_zones:
                    redflags.append(zone)
                    effective.append(pd.Timestamp(feature['properties']['effective']).tz_convert('UTC'))
                    ends.append(pd.Timestamp(feature['properties']['ends']).tz_convert('UTC'))
    
    print(f'{len(redflags)} red flag zones found: {redflags}')

    # Find the shapes associated with the active red flag zones
    rf_zones = []
    to_remove = []
    for i in range(len(redflags)):
        for j in range(len(rf_shape['features'])):
            feature = rf_shape['features'][j]
            try:
                if redflags[i] == feature['properties']['STATE_ZONE']:
                    #print(rf)
                    feature['geom'] = feature['geometry']
                    if feature['geometry']['type'] == 'Polygon': # isinstance(feature['geometry'],Polygon):
                        try:
                            feature['geometry'] = Polygon(feature['geometry']['coordinates'][0])
                        #some zones have the geometry wrong
                        except:
                            feature['geometry']['type'] = 'MultiPolygon'
                            polys = [Polygon(c[0]) for c in feature['geometry']['coordinates']]
                            feature['geometry'] = MultiPolygon(polys)
                    elif feature['geometry']['type'] == 'MultiPolygon': #isinstance(feature['geometry'],MultiPolygon):
                        polys = [Polygon(c[0]) for c in feature['geometry']['coordinates']]
                        feature['geometry'] = MultiPolygon(polys)
                    else:
                        print('Unknown geometry type')
                        to_remove.append(i)
                    rf_zones.append(feature)
                    #dir()print('Have added ',len(rf_zones), 'red flag zones to datframe')
            except:
                 print(type(feature))
                 print('Error with the red flag feature',redflags[i])
                 #redflags.drop(i)
                 to_remove.append(i)

    if len(to_remove):
        r_effective = []
        r_ends = []
        for k in range(len(rf_zones)):
            if k not in to_remove:
                r_effective.append(effective[k])
                r_ends.append(ends[k])
        effective = r_effective
        ends = r_ends

    #try:
    red_flag_zones = pd.DataFrame(rf_zones)
    red_flag_zones['effective'] = effective
    red_flag_zones['ends'] = ends
    #except:
    #    red_flag_zones = pd.DataFrame()
     #   print('Error combining red flags zone, returning empty DataFrame')

    return red_flag_warnings, red_flag_zones.dropna()


def red_flag_incident(lon, lat, ign_utc, red_flag_zones):
    """
    Returns a boolean indicating whether the point (lon, lat) is inside a red flag warning area 
    during the time for which the warning is in effect.
    """
    end_utc = ign_utc + timedelta(hours=24)
    point = Point(lon, lat)
    if len(red_flag_zones):
        for _, rfz in red_flag_zones.iterrows():
            if not (end_utc < rfz['effective'] or ign_utc > rfz['ends']) and point.within(rfz['geometry']):
                return True
        return False
    else:
        return None

def subset_red_flag_data(full_data, data, rf_zones):
    if len(rf_zones) == 0:
        print('No red flag warnings issued')
        return pd.DataFrame()

    temp_data = full_data[full_data.type_description == 'Possible Wildland Fire']
    
    # Get states with active red flag zones and subset data by these states
    rf_states = rf_zones['properties'].apply(lambda x: x['STATE']).unique()
    state_data = temp_data[temp_data.state.isin(rf_states)]
    
    print('State data length:', len(state_data))

    # Identify hotspots within red flag zones
    rf_data = state_data[state_data.apply(lambda row: red_flag_incident(row.lon_tc, row.lat_tc, row.actual_image_time, rf_zones), axis=1)]
    
    print('Red flag detections length:', len(rf_data))

    # Filter out data for known feature tracking IDs
    known_features = data.feature_tracking_id.unique()
    rf_data = rf_data[~rf_data.feature_tracking_id.isin(known_features)]

    # Assign incident_id_string and incident_name
    def generate_incident_info(row):
        fid = row.feature_tracking_id.replace(':', '').replace('-', '').replace('_', '')
        wfo_code = row.wfo_code
        row.incident_id_string = f'{{REDFLAGS-{wfo_code}-{fid[2:6]}-{fid[6:10]}-{fid[-12:]}}}'
        row.incident_name = f'RED_FLAG_{wfo_code}_{fid[2:6]}_{fid[6:10]}_{fid[-12:]}'
        return row

    rf_data = rf_data.apply(generate_incident_info, axis=1)

    return rf_data


if __name__ == '__main__':
    pass

