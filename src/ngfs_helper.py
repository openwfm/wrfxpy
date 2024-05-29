#helps ngfs_start by downloading viir and modis NRT data

from __future__ import absolute_import
from __future__ import print_function
import sys
sys.path.insert(1, '.')
sys.path.insert(1, '../ingest')
sys.path.insert(1, '/data/jhaley/work/wrfxpy/src/ingest/')

import datetime, time
import numpy as np
import pandas as pd
import csv
import json
import os
import shapely as shp


#from ingest.downloader import download_url as du
from shapely.geometry import Polygon, Point, box

from datetime import timedelta, datetime


#from downloader import download_url
#import utils
#dictionary to convert between "CA" and "California", et
# get these two files
#https://firms.modaps.eosdis.nasa.gov/data/active_fire/noaa-20-viirs-c2/csv/J1_VIIRS_C2_USA_contiguous_and_Hawaii_24h.csv --> noaa_20
#https://firms.modaps.eosdis.nasa.gov/usfs/api/kml_fire_footprints/usa_contiguous_and_hawaii/24h/noaa-20-viirs-c2/FirespotArea_usa_contiguous_and_hawaii_noaa-20-viirs-c2_24h.kmz
#https://firms.modaps.eosdis.nasa.gov/data/active_fire/suomi-npp-viirs-c2/csv/SUOMI_VIIRS_C2_USA_contiguous_and_Hawaii_24h.csv --> suomi
#https://firms.modaps.eosdis.nasa.gov/usfs/api/kml_fire_footprints/usa_contiguous_and_hawaii/24h/suomi-npp-viirs-c2/FirespotArea_usa_contiguous_and_hawaii_suomi-npp-viirs-c2_24h.kmz
#https://firms.modaps.eosdis.nasa.gov/data/active_fire/modis-c6.1/csv/MODIS_C6_1_USA_contiguous_and_Hawaii_24h.csv --> modis
#https://firms.modaps.eosdis.nasa.gov/usfs/api/kml_fire_footprints/usa_contiguous_and_hawaii/24h/c6.1/FirespotArea_usa_contiguous_and_hawaii_c6.1_24h.kmz
# os.path.exists(path)

class polar_data():
    # get the csv files above 
    # combine into large pandas data frame
    def __init__(self,csv_timestamp):
        self.data = pd.DataFrame()
        self.timestamp = csv_timestamp
        self.csv_date_str = str(csv_timestamp.year)+'_'+str(csv_timestamp.month)+'_'+str(csv_timestamp)
        self.firms_path = 'ingest/FIRMS/'
        self.satlist = ['noaa_20','noaa_21','suomi','modis','landsat','noaa_20_Alaska','noaa_21_Alaska','suomi_Alaska']
    
    #add nrt data, filter out low confidence detections
    def sort_data_bytime(self):
        self.data = self.data.sort_values(by='acq_time')

    def download_url(self,url,local_path,token=False):
        
        if token:
            print('Trying with tokens')
            str1 = 'wget -e robots=off -m -np -R .html,.tmp -nd '
            tokens = json.load(open('etc/tokens.json'))

            ts = ' --header "Authorization: Bearer {}" '
            token_string = ts.format(tokens['nrt'])
            local_str = ' -P ' + local_path

            cmd = str1 + url + token_string + local_str
        else:
            cmd = 'wget  -P ' +local_path + ' ' + url

        os.system(cmd)
    def filter_low_confidence(self,data_read,sat):
        if sat == 'modis':
            data_read = data_read.loc[data_read['confidence'] > 70]
        else:
            data_read = data_read.loc[data_read['confidence']!='low']
        return data_read
    
    def add_firms_24(self,sat,csv_timestamp):
    
        # Dictionary mapping satellite names to directory and file names
        satellite_info = {
            'noaa_20': ('noaa-20-viirs-c2/csv/', 'J1_VIIRS_C2_USA_contiguous_and_Hawaii_48h.csv'),
            'noaa_21': ('noaa-21-viirs-c2/csv/', 'J2_VIIRS_C2_USA_contiguous_and_Hawaii_48h.csv'),
            'suomi': ('suomi-npp-viirs-c2/csv/', 'SUOMI_VIIRS_C2_USA_contiguous_and_Hawaii_48h.csv'),
            'modis': ('modis-c6.1/csv/', 'MODIS_C6_1_USA_contiguous_and_Hawaii_48h.csv'),
            'landsat': ('landsat/csv/', 'LANDSAT_USA_contiguous_and_Hawaii_48h.csv'),
            'noaa_20_Alaska': ('noaa-20-viirs-c2/csv/', 'J1_VIIRS_C2_Alaska_48h.csv'),
            'noaa_21_Alaska': ('noaa-21-viirs-c2/csv/', 'J2_VIIRS_C2_Alaska_48h.csv'),
            'suomi_Alaska': ('suomi-npp-viirs-c2/csv/', 'SUOMI_VIIRS_C2_Alaska_48h.csv')
        }

        if sat not in self.satlist:
            print('Unknown satellite named', sat)
            return

        # Constructing csv_url and csv_local using satellite_info dictionary
        sat_dir, sat_file = satellite_info[sat]
        csv_url = f'https://firms.modaps.eosdis.nasa.gov/data/active_fire/{sat_dir}{sat_file}'
        csv_local = self.firms_path + sat_file

        #download if the csv file doesn't exist or is older than 7200 seconds = 2 hours
        if not os.path.exists(csv_local):
                self.download_url(csv_url,self.firms_path)
        else:
            if (time.time()- os.path.getmtime(csv_local)) > 1800:
                #copy older file and download again
                cp_str = 'cp ' + csv_local + ' ' + csv_local + '.bak'
                os.system(cp_str)
                del_str = 'rm ' + csv_local
                os.system(del_str)
                self.download_url(csv_url,self.firms_path)

        #read the csv file
        data_read = pd.read_csv(csv_local)
        #make the acq_date into a utc timestamp
        data_read['acq_date'] = pd.to_datetime(data_read['acq_date'], utc = True)
        #make sure the "acq_time" column is a string"
        data_read['acq_time'] = data_read['acq_time'].astype(str)


        #make the acq_date have hours , min, and seconds
        data_read = self.fix_times(data_read)

        #remove this troublesome column
        data_read.drop(columns='satellite',inplace=True)

        #filter out low confidence detections
        data_read = self.filter_low_confidence(data_read,sat)

        #join new download to the data 
        if self.data.shape[0] > 0:
            self.data = pd.merge(self.data,data_read,how = 'outer')
        else:
            self.data = data_read

    def add_firms_dates(self,sat,csv_timestamp,days_to_get):
        #sat is string with sateliite name, modis is not included 
        #csv_datestamp is a pandas.timestamp object
        #gets firms nrt for specific dates, maybe requires tokens
        tokens = json.load(open('etc/tokens.json'))
        #probably should go back at least 24 hours
        #days_to_get = 2

        #get day of the year from the date_str
        day_of_year = csv_timestamp.day_of_year
        year = csv_timestamp.year

        if sat not in self.satlist:
            print('Unknownn satellite named', sat)
            return
        
        #same for all satellites
        url_base = 'https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/'
        #example file path
        #https://nrt3.modaps.eosdis.nasa.gov/api/v2/content/archives/FIRMS/landsat/USA_contiguous_and_Hawaii/LANDSAT_USA_contiguous_and_Hawaii__2022312.txt
        satellite_info = {
            'noaa_20': {
                'csv_file_pattern': 'J1_VIIRS_C2_USA_contiguous_and_Hawaii_VJ114IMGTDL_NRT_{0}{1:03d}.txt',
                'csv_dir': 'noaa-20-viirs-c2/USA_contiguous_and_Hawaii/'
            },
            'noaa_21': {
                'csv_file_pattern': 'J2_VIIRS_C2_USA_contiguous_and_Hawaii_VJ214IMGTDL_NRT_{0}{1:03d}.txt',
                'csv_dir': 'noaa-21-viirs-c2/USA_contiguous_and_Hawaii/'
            },
            'suomi': {
                'csv_file_pattern': 'SUOMI_VIIRS_C2_USA_contiguous_and_Hawaii_VNP14IMGTDL_NRT_{0}{1:03d}.txt',
                'csv_dir': 'suomi-npp-viirs-c2/USA_contiguous_and_Hawaii/'
            },
            'noaa_20_Alaska': {
                'csv_file_pattern': 'J1_VIIRS_C2_Alaska_VJ114IMGTDL_NRT_{0}{1:03d}.txt',
                'csv_dir': 'noaa-20-viirs-c2/Alaska/'
            },
            'noaa_21_Alaska': {
                'csv_file_pattern': 'J2_VIIRS_C2_Alaska_VJ214IMGTDL_NRT_{0}{1:03d}.txt',
                'csv_dir': 'noaa-21-viirs-c2/Alaska/'
            },
            'suomi_Alaska': {
                'csv_file_pattern': 'SUOMI_VIIRS_C2_Alaska_VNP14IMGTDL_NRT_{0}{1:03d}.txt',
                'csv_dir': 'suomi-npp-viirs-c2/Alaska/'
            },
            'landsat': {
                'csv_file_pattern': 'LANDSAT_USA_contiguous_and_Hawaii__{0}{1:03d}.txt',
                'csv_dir' : 'landsat/USA_contiguous_and_Hawaii/'
            }
        }

        #particulars for individual satellites
        for i in range(days_to_get):
            print('Getting polar data from ',day_of_year - i,'th day of ',year)
            
            if sat in satellite_info:
                satellite_data = satellite_info[sat]
                csv_file = satellite_data['csv_file_pattern'].format(year, day_of_year - i)
                csv_dir = satellite_data['csv_dir']
            else:
                print('Unknown or unprogrammed satellite named', sat)
                return


            #remote and local locations of the csv file
            csv_url = url_base + csv_dir + csv_file
            csv_local = self.firms_path + csv_file
            #viirs data fields
            #latitude,longitude,bright_ti4,scan,track,acq_date,acq_time,satellite,confidence,version,bright_ti5,frp,daynight
            #landsat data fields
            #latitude,longitude,path,row,scan,track,acq_date,acq_time,satellite,confidence,daynight

            #downloading not working, files are all in place
            if not os.path.exists(csv_local):
                print('updating VIIRS data cache',csv_local)
                cmd_str = 'wget -e robots=off -m -np -R .html,.tmp -nd "{}" --header "Authorization: Bearer {}" -P ingest/FIRMS'
                #   print(cmd_str.format(url_base,csv_dir[:-1],tokens['nrt']))
                os.system(cmd_str.format(url_base+csv_dir+csv_file,tokens['nrt']))
                
            #if downloads fails
            if not os.path.exists(csv_local):
                print('Missing the file: ',csv_local)
            #else, read the csv file
            else:
                print('\tReading',csv_local)
                data_read = pd.read_csv(csv_local)
                #make the acq_date into a utc timestampjobs 
                data_read['acq_date'] = pd.to_datetime(data_read['acq_date'], utc = True)
                #print(data_read['acq_date'].iloc[0])
                #make sure the "acq_time" column is a string"
                data_read['acq_time'] = data_read['acq_time'].astype(str)

                #remove this troublesome column
                data_read.drop(columns='satellite',inplace=True)

                #turn the weird column of 'acq_time' into houurs, min,ecs and add to the acq_date
                data_read = self.fix_times(data_read)

                #filter out low confidence detections
                data_read = self.filter_low_confidence(data_read,sat)

                #join new download to the data 
                if self.data.shape[0] > 0:
                    self.data = pd.merge(self.data,data_read,how = 'outer')
                else:
                    self.data = data_read


    def fix_times(self,data_read):
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

        return data_read
  
    #convert to work on self
    def best_ign_estimate(self,ignition_pixel,inc_bb):
        #use shapely to help find the first viirs pixel within the boundaries of the goes pixel
        #inc_bb is the incident bounding box (minx, miny, maxx, maxy)
        #GOignition_pixe
        print('\tFinding best viirs pixel for the incident GOES ignition point')
        print('\tTime of GOES imaging: ',ignition_pixel.loc['observation_time'])

        #setup ignition pixel polygon

        vertices = [(ignition_pixel[f'lon_tc_c{i}'], ignition_pixel[f'lat_tc_c{i}']) for i in range(1, 5)]
        ignition_polygon = Polygon(vertices)
        
        #setup incident polygon, easier construction for square box
        incident_polygon = box(inc_bb[0],inc_bb[2],inc_bb[1],inc_bb[3])
        
        #brute force way to get at it, loop through list of all viirs detections
        new_ign_pts = np.empty([0,2])
        new_inc_pts = np.empty([0,2])
        new_ign_times = np.array([])
        new_inc_times = np.array([])
        #print(self.data.dtypes)
        #print('Looking for new detection data for ',self.data['acq_date'])
        for i in range(self.data.shape[0]):
            try:
                #point object
                viirs_point = Point(self.data.loc[i,'longitude'],self.data.loc[i,'latitude'])
                #array of numbers
                viirs_ign = (self.data.loc[i,'latitude'],self.data.loc[i,'longitude'])
                #print('i= ',i,self.data.data.loc[i])
            except:
                print('Error finding viirs point')
                viirs_point = Point(600,3000) # <<------- not in any polygon
            
            if viirs_point.within(incident_polygon):
                #print('i = ', i,self.data.loc[i]['acq_date'])
                #time.sleep(1)
                #print('\tFound additional incident pixel')
                new_inc_pts = np.append(new_inc_pts,[viirs_ign],axis = 0)
                new_inc_times = np.append(new_inc_times,self.data.loc[i,'acq_date'])
                if viirs_point.within(ignition_polygon):
                    new_ign_pts = np.append(new_ign_pts,[viirs_ign],axis = 0)
                    new_ign_times = np.append(new_ign_times,self.data.loc[i,'acq_date'])
                    #print('\tFound viirs pixel within earliest GOES pixel:',viirs_ign,self.data.loc[i,'acq_date'])
                
        if (new_ign_pts.shape[0] == 0 & new_inc_pts.shape[0] != 0):
            print('\tUsing earliest VIIRS pixel, but it is not within the earliest GOES Pixel')
            new_ign_pts = new_inc_pts
            new_ign_times = new_inc_times
            print('\t',new_ign_pts.shape[0],' new incident points found')

        if new_ign_pts.shape[0]>0:
            #sort data by time
            srt_idx = np.array(new_ign_times.argsort())
            new_ign_pts = new_ign_pts[srt_idx]
            new_ign_times = new_ign_times[srt_idx]
            min_ign_time = new_ign_times[0]
            #boolean index
            first_viirs_pixels = new_ign_times == min_ign_time
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

        return new_ign_latlon, new_ign_utc            
