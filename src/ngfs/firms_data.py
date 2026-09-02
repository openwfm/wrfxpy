#tools for downloading viirs URT and NRT data
from __future__ import absolute_import
from __future__ import print_function
import sys, logging
sys.path.insert(1, '.')
sys.path.insert(1, '../ingest')
sys.path.insert(1, '/data/jhaley/work/wrfxpy/src/ingest/')
import glob
import time
import numpy as np
import pandas as pd
import csv
import json
import os, zipfile
import shapely as shp
import xml.etree.ElementTree as ET
import re
#from ingest.downloader import download_url as du
from shapely.geometry import Polygon, Point, box
from datetime import timedelta, datetime
from ftplib import FTP
from pathlib import Path

###### for get URT data ####################################################
def extract_placemark(data):
    # for processing FIRMS URT kmz/kml files
    # Define patterns for each piece of VIIRS data
    # Landsat with require a different dictionary, based on the following type
    '''
                    <b>Latitude: </b> 45.731552<br/>
                    <b>Longitude: </b> -121.5075<br/>
                    <b>Detection Time: </b> 2024-09-03 18:55 UTC<br/>
                    <b>Sensor: </b> Landsat 8 OLI<br/>
                    <b>Confidence: </b> High<br/>
                    <b>Day/Night: </b> Day<br/>
                    <b>WRS-2 Path: </b> 046<br/>
                    <b>WRS-2 Row: </b> 028
    '''
    patterns = {
        'latitude': r'<b>Latitude: </b> ([\d.-]+)<br/>',
        'longitude': r'<b>Longitude: </b> ([\d.-]+)<br/>',
        'detection_time': r'<b>Detection Time: </b> ([\d-]+) ([\d:]+) UTC<br/>',
        'sensor': r'<b>Sensor: </b> (.+?)<br/>',
        'confidence': r'<b>Confidence: </b> (.+?)<br/>',
        'day_night': r'<b>Day/Night: </b> (.+?)<br/>',
        'scan': r'<b>Scan: </b> ([\d.]+) km<br/>',
        'track': r'<b>Track: </b> ([\d.]+) km<br/>',
        'frp': r'<b>FRP: </b> ([\d.]+) MW<br/>',
        'brightness': r'<b>Brightness: </b> ([\d.]+) K'
    }

    extracted_data = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, data)
        if match:
            value = match.group(1)
            if key in ['latitude', 'longitude', 'scan', 'track', 'frp', 'brightness']:
                # Convert to float
                extracted_data[key] = float(value)
            elif key == 'detection_time':
                # Convert to date and time
                date_str, time_str = match.groups()
                date_obj = datetime.strptime(f"{date_str} {time_str}", '%Y-%m-%d %H:%M')
                ts = pd.Timestamp(date_obj,tz='UTC')
                extracted_data['acq_date'] = ts #.strftime('%Y-%m-%d')
                extracted_data['acq_time'] = date_obj.strftime('%H%M')
            elif key == 'sensor':
                # Extract sensor as satellite
                extracted_data['satellite'] = 'N'  # or another mapping based on your data
            elif key == 'confidence':
                # Convert confidence to lowercase
                extracted_data['confidence'] = value.lower()
            elif key == 'day_night':
                # Map day/night to desired format
                extracted_data['daynight'] = value[0].upper()
            else:
                extracted_data[key] = value
                
    #  Define the CSV headers
    #  headers = ['latitude', 'longitude', 'bright_ti4', 'scan', 'track', 'acq_date', 'acq_time', 'satellite', 'confidence', 'version', 'bright_ti5', 'frp', 'daynight']
                
    csv_row = {
    'latitude': extracted_data.get('latitude', ''),
    'longitude': extracted_data.get('longitude', ''),
    'bright_ti4': extracted_data.get('brightness', ''),  # Assuming brightness maps to bright_ti4
    'scan': extracted_data.get('scan', ''),
    'track': extracted_data.get('track', ''),
    'acq_date': extracted_data.get('acq_date', ''),
    'acq_date_time': extracted_data.get('acq_date', ''),
    'pixel_date_time': extracted_data.get('acq_date', ''),
    'acq_time': extracted_data.get('acq_time', ''),
    'satellite': extracted_data.get('satellite', ''),
    'confidence': extracted_data.get('confidence', ''),
    'version': 'URT',  # Assuming a constant version value
    'bright_ti5': extracted_data.get('brightness', ''),  # Assuming brightness maps to bright_ti5
    'frp': extracted_data.get('frp', ''),
    'daynight': extracted_data.get('daynight', '')
    }

    return pd.DataFrame([csv_row])

def filter_low_confidence(data_read,sat):
        if sat == 'modis':
            data_read = data_read.loc[data_read['confidence'] > 70]
        else:
            data_read = data_read.loc[data_read['confidence']!='low']
        return data_read

def download_url(self,url,local_path,token=False):
        if token:
            print('Trying with tokens')
            str1 = 'wget -e robots=off -m -np -R .html,.tmp -nd '
            tokens = json.load(open('etc/tokens.json'))

            ts = ' --header "Authorization: Bearer {}" '
            token_string = ts.format(tokens['nrt'])
            local_str = ' -O ' + local_path

            cmd = str1 + url + token_string + local_str
        else:
            cmd = 'wget  -O ' +local_path + ' ' + url
        os.system(cmd)

def fix_times(data_read):
        #converts the acq_time to a column with hours and minutes and adds to the acq_date object
        #time_type depends on the csv time_acq format. For the 24_hour nrt, format like 912 = 09:12
        #archived data has format 09:12
        #this function workswith 0912 instead of 912 or 09:12
        acq_int = data_read['acq_time']
        #remove the colon
        acq_str = [str(i).replace(':','') for i in acq_int]
        #now pad zeros to the left of strings like '912'
        acq_str = [str(i).rjust(4,'0') for i in acq_str]

        t_fmt = '%H%M'
        t1 = datetime.strptime('0000',t_fmt)
        for i in range(data_read.shape[0]):
            t2 = datetime.strptime(acq_str[i],t_fmt)

            data_read.loc[i,'acq_date'] = data_read.loc[i,'acq_date'] + (t2-t1)
        
        data_read['acq_date'] = pd.to_datetime(data_read['acq_date'])
        data_read['pixel_date_time'] = data_read['acq_date']
        data_read['acq_date_time'] = data_read['acq_date']
        return data_read

def add_firms_nrt(ngfs_cfg,sat,csv_timestamp=None,df = pd.DataFrame()):
    
        # Dictionary mapping satellite names to directory and file names
        satellite_info = {
            'noaa_20': ('noaa-20-viirs-c2/csv/', 'J1_VIIRS_C2_USA_contiguous_and_Hawaii_48h.csv'),
            'noaa_21': ('noaa-21-viirs-c2/csv/', 'J2_VIIRS_C2_USA_contiguous_and_Hawaii_48h.csv'),
            'suomi': ('suomi-npp-viirs-c2/csv/', 'SUOMI_VIIRS_C2_USA_contiguous_and_Hawaii_48h.csv'),
            'landsat': ('landsat/csv/', 'LANDSAT_USA_contiguous_and_Hawaii_48h.csv'),
            'noaa_20_Alaska': ('noaa-20-viirs-c2/csv/', 'J1_VIIRS_C2_Alaska_48h.csv'),
            'noaa_21_Alaska': ('noaa-21-viirs-c2/csv/', 'J2_VIIRS_C2_Alaska_48h.csv'),
            'suomi_Alaska': ('suomi-npp-viirs-c2/csv/', 'SUOMI_VIIRS_C2_Alaska_48h.csv')
        }

        sat_list = ngfs_cfg['firms_cfg']['sats']
        if sat not in sat_list:
            print('Unknown satellite named', sat)
            return

        # Constructing csv_url and csv_local using satellite_info dictionary
        sat_dir, sat_file = satellite_info[sat]
        csv_url = f'https://firms.modaps.eosdis.nasa.gov/data/active_fire/{sat_dir}{sat_file}'
        ingest_path = ngfs_cfg['firms_cfg']['ingest_directory']
        csv_local = ingest_path + sat_file

        #download if the csv file doesn't exist or is older than 600 seconds = 10 minutes
        if not os.path.exists(csv_local):
                download_url(csv_url,csv_local)
        else:
            if (time.time()- os.path.getmtime(csv_local)) > 600:
                #copy older file and download again
                cp_str = 'cp ' + csv_local + ' ' + csv_local + '.bak'
                os.system(cp_str)
                del_str = 'rm ' + csv_local
                os.system(del_str)
                download_url(csv_url,csv_local)    

        #read the csv file
        data_read = pd.read_csv(csv_local)

        #make the acq_date into a utc timestamp
        data_read['acq_date'] = pd.to_datetime(data_read['acq_date'], utc = True)
        #make sure the "acq_time" column is a string"
        data_read['acq_time'] = data_read['acq_time'].astype(str)

        #make the acq_date have hours , min, and seconds
        data_read = fix_times(data_read)

        #remove this troublesome column
        data_read.drop(columns='satellite',inplace=True)

        #filter out low confidence detections
        data_read = filter_low_confidence(data_read,sat)

        #join new download to the data 
        if len(df) > 0:
            df = pd.merge(df,data_read,how = 'outer')
        else:
            df = data_read
        
        return df

def add_firms_urt(ngfs_cfg,sat,csv_timestamp=None,df = pd.DataFrame()):
        #(region, satellite_name) in API

        #Alaska files seem to be empty, return if requested
        #Landsat files requires a different dictionary to process
        if 'Alaska' in sat or 'landsat' in sat:
            return
        
        satellite_info = {
            'noaa_20': ('usa_contiguous_and_hawaii', 'noaa-20-viirs-c2'),
            'noaa_21': ('usa_contiguous_and_hawaii', 'noaa-21-viirs-c2'),
            'suomi': ('usa_contiguous_and_hawaii', 'suomi-npp-viirs-c2'),
            'landsat': ('usa_contiguous_and_hawaii', 'landsat'),
            'noaa_20_Alaska': ('alaska', 'noaa-20-viirs-c2'),
            'noaa_21_Alaska': ('alaska', 'noaa-21-viirs-c2'),
            'suomi_Alaska': ('alaska', 'suomi-npp-viirs-c2')
        }
        sat_list = ngfs_cfg['firms_cfg']['sats']
        if sat not in sat_list:
            print('Unknown satellite named', sat)
            return

        #download the URT kmz file
        #get time infomration for the timestamp, for saving the data
        if not csv_timestamp:
            csv_timestamp = pd.Timestamp.now('UTC')
        year = csv_timestamp.year
        month = csv_timestamp.month
        day = csv_timestamp.day
        hour = csv_timestamp.hour

        region, satellite = satellite_info[sat]
        url_format = "https://firms.modaps.eosdis.nasa.gov/api/kml_fire_footprints/{0}/24h/{1}/FirespotArea_{0}_{1}_24h.kmz"
        kmz_url = url_format.format(region,satellite)
        ingest_path = ngfs_cfg['firms_cfg']['ingest_directory']
        local_format = ingest_path + 'URT_' + satellite+'_'+region + '_{0}{1:02d}{2:02d}_{3:02d}h.kmz'
        kmz_local = local_format.format(year,month,day,hour)
        print('Downloading ',kmz_url,' to ' ,kmz_local)

        #download if the kmz file doesn't exist or is older than 600 seconds = 10 minutes
        if not os.path.exists(kmz_local):
                os.system('wget -O ' + kmz_local + ' ' + kmz_url)
                #self.download_url(kmz_url,kmz_local)
        else:
            if (time.time()- os.path.getmtime(kmz_local)) > 600:
                #copy older file and download again
                cp_str = 'cp ' + kmz_local + ' ' + kmz_local + '.bak'
                os.system(cp_str)
                del_str = 'rm ' + kmz_local
                os.system(del_str)
                os.system('wget -O ' + kmz_local + ' ' + kmz_url)
                #self.download_url(kmz_url,kmz_local)


        print('Extracting the kmz to a kml file')
        #extract the kmz file to a kml which can be read
        kmzip = zipfile.ZipFile(kmz_local,'r')
        kml_filename = kmzip.infolist()[0].filename
        #extract file and gets path in one step
        kml_path = kmzip.extract(kml_filename,path=ingest_path)

        #parse the kml file
        tree = ET.parse(kml_path)
        root = tree.getroot()
        folders = root[0].findall('{http://earth.google.com/kml/2.1}Folder')

        centroids = [] # probably not needed
        placemarks = []
        for fol in folders:
            #for Landast, it looks like '30m Fire Detection Centroids (Last 0 to 6hrs)'
            if 'Fire Detection Centroids' in fol[0].text:
                centroids.append(fol)
                placemarks.append(fol.findall('{http://earth.google.com/kml/2.1}Placemark'))
        
        # placemarks is now a list of len 4, with all the detections for different time periods
        # placemarks[0][0].findall('{http://earth.google.com/kml/2.1}description') will get the description of the first elemt in the 
        # set of detections 0-6 hours old
        # placemarks[0][1][1].text will get the text description of the second element in the 
        # set of detections 0-6 hours old
        
        if not len(placemarks):
            print('Empty KML file, returning....')
            return

        data_read = pd.DataFrame()
        for pl in placemarks:
            for detection in pl:
                data_read = data_read.append(extract_placemark(detection[1].text),ignore_index=True)

        #now data_read looks like it came from one of the FIRMS csv files
        #make the acq_date into a utc timestamp
        #data_read['acq_date'] = pd.to_datetime(data_read['acq_date'], utc = True)
        #make sure the "acq_time" column is a string"
        #['acq_time'] = data_read['acq_time'].astype(str)

        #make the acq_date have hours , min, and seconds
        #data_read = self.fix_times(data_read)

        #remove this troublesome column
        data_read.drop(columns='satellite',inplace=True)

        #filter out low confidence detections
        data_read = filter_low_confidence(data_read,sat)

        #join new download to the data 
        if len(df) > 0:
            df = pd.merge(df,data_read,how = 'outer')
        else:
            df = data_read

        return df


def fix_times(data_read):
        #converts the acq_time to a column with hours and minutes and adds to the acq_date object
        #time_type depends on the csv time_acq format. For the 24_hour nrt, format like 912 = 09:12
        #archived data has format 09:12
        #this function workswith 0912 instead of 912 or 09:12
        acq_int = data_read['acq_time']
        #remove the colon
        acq_str = [str(i).replace(':','') for i in acq_int]
        #now pad zeros to the left of strings like '912'
        acq_str = [str(i).rjust(4,'0') for i in acq_str]

        t_fmt = '%H%M'
        t1 = datetime.strptime('0000',t_fmt)
        for i in range(data_read.shape[0]):
            t2 = datetime.strptime(acq_str[i],t_fmt)

            data_read.loc[i,'acq_date'] = data_read.loc[i,'acq_date'] + (t2-t1)
        
        data_read['acq_date'] = pd.to_datetime(data_read['acq_date'])
        data_read['pixel_date_time'] = data_read['acq_date']
        data_read['acq_date_time'] = data_read['acq_date']  ## this is the NGFS time variable name
        return data_read


### get rid of sel. pass in the datframe with all viirs detections   #probaly not needed since polar_ignition_estimate is in ngfs_incident
def best_ign_estimate(viirs_data,ignition_pixel,inc_bb):
        #use shapely to help find the first viirs pixel within the boundaries of the goes pixel
        #inc_bb is the incident bounding box (minx, miny, maxx, maxy)
        #GOignition_pixe
        print('\tFinding best viirs pixel for the incident GOES ignition point')
        try:
            print('\tTime of GOES imaging: ',ignition_pixel.loc['observation_time'])
        except:
            print('\tTime of GOES imaging: ',ignition_pixel.loc['pixel_date_time'])

        #setup ignition pixel polygon
        try: #NGFS Version 1.x
            vertices = [(ignition_pixel[f'lon_tc_c{i}'], ignition_pixel[f'lat_tc_c{i}']) for i in range(1, 5)]
        except: #NGFS Version 2+
            vertices = [(ignition_pixel[f'longitude_c{i}'], ignition_pixel[f'latitude_c{i}']) for i in range(1, 5)]
        ignition_polygon = Polygon(vertices)
        
        #setup incident polygon, easier construction for square box
        incident_polygon = box(inc_bb[0],inc_bb[2],inc_bb[1],inc_bb[3])
        
        #brute force way to get at it, loop through list of all viirs detections
        viirs_inc_data = pd.DataFrame()
        viirs_ign_data = pd.DataFrame()
        new_ign_pts = np.empty([0,2])
        new_inc_pts = np.empty([0,2])
        new_ign_times = np.array([])
        new_inc_times = np.array([])
        #print(viirs_data.dtypes)
        #print('Looking for new detection data for ',viirs_data['acq_date'])
        for i in range(viirs_data.shape[0]):
            try:
                #point object
                viirs_point = Point(viirs_data.loc[i,'longitude'],viirs_data.loc[i,'latitude'])
                #array of numbers
                viirs_ign = (viirs_data.loc[i,'latitude'],viirs_data.loc[i,'longitude'])
                #print('i= ',i,viirs_data.data.loc[i])
            except:
                print('Error finding viirs point')
                viirs_point = Point(600,3000) # <<------- not in any polygon
            
            if viirs_point.within(incident_polygon):
                #print('i = ', i,viirs_data.loc[i]['acq_date'])
                #time.sleep(1)
                #print('\tFound additional incident pixel')
                viirs_inc_data = viirs_inc_data.append(viirs_data.loc[i])
                new_inc_pts = np.append(new_inc_pts,[viirs_ign],axis = 0)
                new_inc_times = np.append(new_inc_times,viirs_data.loc[i,'acq_date'])
                if viirs_point.within(ignition_polygon):
                    viirs_ign_data = viirs_ign_data.append(viirs_data.loc[i])
                    new_ign_pts = np.append(new_ign_pts,[viirs_ign],axis = 0)
                    new_ign_times = np.append(new_ign_times,viirs_data.loc[i,'acq_date'])
                    #print('\tFound viirs pixel within earliest GOES pixel:',viirs_ign,viirs_data.loc[i,'acq_date'])
        print('\t',len(new_inc_pts),'viirs incident detections, ',len(new_ign_pts),' viirs ignition polygon detections')
        if (len(new_ign_pts) == 0 and len(new_inc_pts) > 0):
            print('\tUsing earliest VIIRS pixel, but it is not within the earliest GOES Pixel')
            new_ign_pts = new_inc_pts 
            new_ign_times = new_inc_times 
            viirs_ign_data = viirs_ign_data.append(viirs_inc_data)
            print('\t',new_ign_pts.shape[0],' new incident points found')

        if new_ign_pts.shape[0]>0:
            #sort data by time
            srt_idx = np.array(new_ign_times.argsort())
            new_ign_pts = new_ign_pts[srt_idx]
            new_ign_times = new_ign_times[srt_idx]
            min_ign_time = new_ign_times[0]
            #boolean index
            first_viirs_pixels = new_ign_times <= min_ign_time + timedelta(hours=3)
            print('\tCount of earliest viirs pixels: ',sum(first_viirs_pixels))

            #replace with better estimate of lat/lon pf ignition
            mean_viirs_lat = np.mean(new_ign_pts[first_viirs_pixels,0])
            mean_viirs_lon = np.mean(new_ign_pts[first_viirs_pixels,1])

            # return these to be added to the ngfs_incidents
            new_ign_latlon = [mean_viirs_lat,mean_viirs_lon]
            new_ign_utc = new_ign_times[0]
            print('\tBest viirs ignition estimate location: \n\t\t',new_ign_latlon,new_ign_utc)
        else:
            new_ign_latlon = [float('NaN'),float('NaN')]
            new_ign_utc = pd.NaT  # << -------- return something of correct data type

        return new_ign_latlon, new_ign_utc, viirs_ign_data 


if __name__ == "__main__":
    pass