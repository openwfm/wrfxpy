#handles downloading and preparing NGFS data via csv file from ftp site
import pandas as pd
import os, sys
sys.path.insert(1, 'src/')
sys.path.insert(1, 'src/ingest')
from datetime import datetime, timedelta, timezone
from ingest.downloader import download_url
import ngfs_dictionary as nd
from ngfs import config_manager as con_man
from ftplib import FTP
from pathlib import Path
import re
import glob
import geopandas as gpd

def today_datetime(date):
    #accepts a string either 'today' or 'YYY-MM-DD' and returns a datetime object
    if date=='today':
        today = datetime.now(timezone.utc).date()
    else:
        today = datetime.strptime(date, "%Y-%m-%d").date()
    return today

def parse_times(df):
    time_keys = ["acq_date_time","pixel_date_time"]
    for k in df.keys():
        if 'time' in k and k not in time_keys:
            time_keys.append(k)
    for tk in time_keys:
        if tk in df.keys():
            df[tk] = pd.to_datetime(df[tk], utc=True)
        else:
            print(f'{tk} is an unknown key in DatFrame, skipping...')
    return df

def get_ngfs_csv(days_previous,sat,domain):
   #downloads daily NGFS csv file, returns file name and local paths
   #example csv name:
   #   NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2023_05_25_145.csv

   #name the csv file
   # Get the date components
   csv_day = datetime.now() - timedelta(days=days_previous)
   day_of_year = csv_day.timetuple().tm_yday
   yyyy, mm, dd = csv_day.year, csv_day.month, csv_day.day

 
   #join strings to have name and url of the csv file
   # Define the base of the csv_str
   #base_str = 'NGFS_FIRE_DETECTIONS_GOES-{}_ABI_CONUS_{}_{}_{}_{}.csv'
   base_str = 'NGFS_FIRE_DETECTIONS_GOES-{}_ABI_{}_{}_{}_{}_{}.csv'

   # Determine the satellite
   satellite = str(sat)

   # Format the csv_str using satellite and other variables
   csv_str = base_str.format(satellite,domain,yyyy, str(mm).zfill(2), str(dd).zfill(2), str(day_of_year).zfill(3))

   # Define the csv_url using csv_str
   csv_url = f'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/GOES-{"WEST" if sat == 18 else "EAST"}/{domain}/{csv_str}'

  
   print('\tDownloading: ',csv_str)
   csv_path = 'ingest/NGFS/'+csv_str
   # today's and yesterday's
   if days_previous < 2:
      download_url(csv_url,csv_path)
   #older file should be the same as what is cached
   else:
      if os.path.exists(csv_path):
         print('\tUsing older data in ingest path')
      else:
         download_url(csv_url,csv_path)
         
   return csv_str, csv_path

def get_csv_lists(days_to_get,sat_sectors = None):
   #returns list of csv files to be obtained
   #empy lists for storing path names
   csv_str = list()
   csv_path = list()
   if sat_sectors is None:
       sat_sectors = [(18, 'CONUS'), (19,'CONUS'), (18, 'Full-Disk')]
   #number of days of data
   #days_to_get = 
   for i in range(days_to_get):
      print('NGFS data from today' if i == 0 else 'NGFS data from yesterday or before') 
      for satellite, sector in sat_sectors:
         # fix it so download error doesn't crash, also this list should come from the ngfs.json config file
         try:
            c_str, c_path = get_ngfs_csv(i, satellite, sector)
            csv_str.append(c_str)
            csv_path.append(c_path)
         except:
            print('Error getting the csv for GOES-',satellite,sector)
      
   return csv_str, csv_path


def setup_csv_download(ngfs_cfg):
    # Determine if working with today's data
    today = 'now' in sys.argv or ngfs_cfg['run_cfg'].get('today_forecasts', False)
    #check if a csv file was used as system argument
    for sa in sys.argv:
        if '.csv' in sa:
            csv_file = [sa]
            today = False
            print('Using existing CSV file:')
            print(f'\t{csv_file[0]}')
            return today, csv_file
    #determine the satelite sectors to get
    sat_sectors = []
    try:
        for i in ngfs_cfg['goes_cfg']['sat_sectors']:
            for j in ngfs_cfg['goes_cfg']['sat_sectors'][i]:
                sat_sectors.append((int(i),j))
    except:
        print('Error finding satelite sectors from configuation')
        sat_sectors = [(18, 'CONUS'), (19,'CONUS'), (18, 'Full-Disk')]
    if today:
        print('Downloading latest CSV files:')
        # Configure download path and create directory if needed
        ngfs_ingest = ngfs_cfg['goes_cfg'].get('ingest_dir_daily','ingest/NGFS') #for daily files
        days_to_get = ngfs_cfg['goes_cfg'].get('days_to_get', 2)
        print(f'\tWill download to: {ngfs_ingest}')
        
        # Get list of files to download
        csv_str, csv_file = get_csv_lists(days_to_get,sat_sectors=sat_sectors)
    
    return today, csv_file

def read_NGFS_csv_data(csv_file):
    """
    Reads and merges CSV file(s), assigning a date to them too.
    
    Parameters:
    - csv_file (str or list): A string path or a list of strings representing the CSV files.
    
    Returns:
    - data (pd.DataFrame): Merged data from CSV files.
    - csv_date_str (str): Date string extracted from the first file name.
    """
    
    # Convert csv_file to list if it's not already
    csv_file = [csv_file] if not isinstance(csv_file, list) else csv_file


    def get_alaska_wfo(data):
      borough_to_wfo = {
      'Anchorage Municipality': 'KAFC',   # Anchorage Forecast Office
      'Fairbanks North Star Borough': 'KAFG',  # Fairbanks Forecast Office
      'Juneau City and Borough': 'KAJK',  # Juneau Forecast Office
      'Matanuska-Susitna Borough': 'KAFC',
      'Kenai Peninsula Borough': 'KAFC',
      'North Slope Borough': 'KAFG',
      'Nome Census Area': 'KAFG',
      'Ketchikan Gateway Borough': 'KAJK',
      'Yakutat City and Borough': 'KAJK',
      'Bethel Census Area': 'KAFG',
      'Dillingham Census Area': 'KAFC',
      'Kusilvak Census Area': 'KAFG',
      'Valdez-Cordova Census Area': 'KAFC',
      'Wrangell City and Borough': 'KAJK',
      'Haines Borough': 'KAJK',
      'Petersburg Borough': 'KAJK',
      'Sitka City and Borough': 'KAJK',
      'Prince of Wales-Hyder Census Area': 'KAJK',
      'Kodiak Island Borough': 'KAFC',
      'Bristol Bay Borough': 'KAFC',
      'Aleutians East Borough': 'KAFC',
      'Aleutians West Census Area': 'KAFC',
      'Denali Borough': 'KAFG',
      'Southeast Fairbanks Census Area': 'KAFG',
      'Yukon-Koyukuk Census Area': 'KAFG',
      'Lake and Peninsula Borough': 'KAFC',
      'Hoonah-Angoon Census Area': 'KAJK',
      'Skagway Municipality': 'KAJK',
      'Copper River Census Area': 'KAFC',
      'Northwest Arctic Borough': 'KAFG'
      # Add more as needed
      }
      data = data[data.state == 'AK']
      wfo_list = []
      for i in range(len(data)):
         try:
            wfo_list.append(borough_to_wfo[data['county'].iloc[i]])
         except:
            wfo_list.append('Unknown')
      data['nws_wfo_code'] = wfo_list
      return data
          

    # Initialize an empty DataFrame
    data = pd.DataFrame()

    # Loop through each CSV file
    for file_path in csv_file:
        if not os.path.exists(file_path):
            print(f'File {file_path} not found')
            continue
        # Modify feature_tracking_id string for GOES-16 or GOES-19 files
        if 'GOES-16' in file_path or 'GOES-19' in file_path:
            os.system(f"sed -i 's/,ID-20/,I6-20/g' {file_path}")
        try:
            # Attempt to read as v2 CSV
            data_read = pd.read_csv(file_path)
            data_read = parse_times(data_read)
            if 'Full' in file_path:
                print(f'\nReading full-disk data: {file_path}')
                print(f'\tNumber of all Full-Disk detections: {len(data_read)}')
                data_read = data_read[data_read.country == 'United States']
                print(f'\tNumber of USA Full-Disk detections: {len(data_read)}')
                data_read = get_alaska_wfo(data_read)

            # Update fields to match legacy NGFS fields
            '''
            data_read['initial_observation_time'] = data_read[time_cols_v2[1]]
            data_read['incident_start_time'] = data_read[time_cols_v2[1]]
            data_read.rename(columns=nd.v2_to_v1, inplace=True)
            data_read = data_read.astype(nd.v2_dict)
            '''
        except Exception:
            # If reading as v2 fails, try reading as v1
            data_read = pd.read_csv(file_path)
            data_read = parse_times(data_read)
            data_read['actual_image_time'] = data_read['observation_time']
            try:
                data_read = data_read.astype(nd.ngfs_dictionary)
            except Exception as e:
                print('Trouble parsing data types, probably early ngfs version')
                print(e)
                data_read = pd.DataFrame()

        # Merge the current data with existing data
        if data.empty:
            print(f'Reading first CSV: {file_path}')
            data = data_read
        else:
            print(f'Merging with CSV file: {file_path}')
            try:
                data = pd.concat([data, data_read],ignore_index=True)
            except Exception as e:
                print(f"Failed to merge with {file_path}: {str(e)}")

        print(f'\tNumber of detections in CSV: {len(data_read)}')
        print(f'\tTotal number of detections: {len(data)}')

    #give a new index to data
    data.reset_index(drop=True, inplace=True)
    # Extract date string from the first file name
    csv_date_str = csv_file[0][-18:-8]

    # Ensure times are in UTC
    '''
    try:
        data[time_cols_v1[1]] = pd.to_datetime(data[time_cols_v1[1]], utc=True)
    except Exception:
        print('\tData column has correct timezone already')
    '''

    return data, csv_date_str

def get_ngfs_data(date='today',ngfs_cfg=None,data=pd.DataFrame()):
    #for downlaoding using the ngfs daily files
    if ngfs_cfg is None:
        ngfs_cfg, wrfxpy_cfg = con_man.load_cfgs()
    
    today, csv_file = setup_csv_download(ngfs_cfg)
    data, csv_date_str = read_NGFS_csv_data(csv_file)
    if len(data) > 0:
        print('Merging with existing data')
        data = pd.concat([data,data],ignore_index=True)
        data = data.drop_duplicates()
    print(data)
    print(f'Dataframe size: {len(data)}')
    return data, csv_date_str

def parse_ts(name):
            #for file names like NGFS_FIRE_DETECTIONS_NOAA-20_VIIRS_2026_03_11_070_06_17_03.csv
            pattern = re.compile(r'(\d{4})_(\d{2})_(\d{2})_\d{3}_(\d{2})_(\d{2})_(\d{2})')
            m = pattern.search(name)
            #return datetime(*map(int, m.groups())) if m else None
            return pd.Timestamp(*map(int, m.groups())).tz_localize('UTC') if m else None


def download_ngfs_scene(ngfs_cfg):
    #for file names like NGFS_FIRE_DETECTIONS_GOES-19_ABI_CONUS_2026_05_14_134_12_56_17.csv
    #a little difference between GOES and VIIRS
    #https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_scene/GOES-WEST/CONUS/NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2026_05_13_133_02_01_17.csv
    #https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_scene/VIIRS/SSEC-DB/NGFS_FIRE_DETECTIONS_NOAA-20_VIIRS_2026_05_13_133_06_35_29.csv
    #all viirs files are in the same directory, rerdless of satellite. GOES files are in different directories
    #scene files look like this : NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2026_05_13_133_02_01_17.csv
    #daily files look like this: NGFS_FIRE_DETECTIONS_GOES-18_ABI_CONUS_2026_05_10_130.csv


    #get subdirectory based on satellite, will be one of the following three:
    ##https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_scene/VIIRS/SSEC-DB/
    ##https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_scene/GOES-WEST/CONUS/
    ##https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_scene/GOES-EAST/CONUS/
    ##https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_scene/GOES-WEST/Subdomains/Alaska-Full-Disk/

    sat_dirs = ['VIIRS/SSEC-DB','GOES-WEST/CONUS','GOES-EAST/CONUS','GOES-WEST/Subdomains/Alaska-Full-Disk/']

    
    for s in sat_dirs:
        if 'VIIRS' in s:
            sat_key = 'viirs_cfg'
        else:
            sat_key = 'goes_cfg'

        #remote ftp storage
        host = ngfs_cfg[sat_key].get('host','bin.ssec.wisc.edu')
        remote_directory = f"pub/volcat/fire_csv/NGFS_scene/{s}"    ####<<<<--------- fix this so it's in the ngfs_cfg file
        print(f'Downloading latest file from {remote_directory}')
        data_format = ngfs_cfg[sat_key]['data_format']

        #local storage
        ingest_dir = ngfs_cfg[sat_key]['ingest_directory']
        local_dir = Path(ingest_dir)
        local_dir.mkdir(parents=True, exist_ok=True)
        
        # find the files already cached and get the time of the last
        local_times = [parse_ts(p.name) for p in local_dir.glob(data_format)]
        last_ts = max([t for t in local_times if t],default = (pd.Timestamp.now("UTC")-timedelta(days=5)))
        #print(last_ts,type(last_ts))

        #ftp setup
        ftp = FTP(host); ftp.login(); ftp.cwd(remote_directory)

        #find available file on the ftp site
        files = sorted((parse_ts(f), f) for f in ftp.nlst() if f.endswith(data_format))
        print(f'There are {len(files)} {data_format} files on the FTP site')

        down_count = []
        for ts, fname in files:
            if ts and ts > last_ts:
                if not os.path.exists(f'{ingest_dir}/{fname}'):
                    down_count.append(fname)
                    print("\tDownloading", fname)
                    with open(local_dir / fname, "wb") as f:
                       ftp.retrbinary(f"RETR {fname}", f.write)
                else:
                    pass
                    #print('Skipping ', fname, ' file already downloaded')
        print(f'Downloaded {len(down_count)} {data_format} files')
        ftp.quit()


def add_ngfs_scene(ngfs_cfg,
                   data = pd.DataFrame(),
                   parse_times = True,
                   start_time = None,
                   end_time = None,
                   sat = 'goes'
    ):
    #update the cached data
    download_ngfs_scene(ngfs_cfg)

    if sat == 'goes':
        sat_key = 'goes_cfg'
    else:
        sat_key = 'viirs_cfg'
    #determine which data to add
    if len(data) > 0: #updating data
        #check which satelite the dataframe is from
        data_sat = data.satellite.unique()[0]
        if data_sat in ['SNPP', 'NOAA-20', 'NOAA-21']:
            sat_key = 'viirs_cfg'
        else:
            sat_key = 'goes_cfg'

        start_time = data.acq_date_time.max()
        end_time = pd.Timestamp.now('UTC')
        print(f'Udating data of length {len(data)} ')
    if end_time and not start_time:
        start_time = end_time - timedelta(hours = ngfs_cfg[sat_key]['days_to_get']*24)

    ingest_dir = ngfs_cfg[sat_key]['ingest_directory']
    data_format = ngfs_cfg[sat_key]['data_format'] 
    local_dir = Path(ingest_dir)
    print(f'Adding data in {ingest_dir}')


    g = glob.glob(f'{ingest_dir}/*{data_format}')
    print('Length of NGFS GOES file list :',len(g))
    print(f'Will load file between {start_time} and  {end_time}')

    #list of dataframes #make a list of pd.DataFrames
    dfs = []
    for f in local_dir.glob(f'*.{data_format}'):
        ts = parse_ts(f.name)    ##### <<<---------------------------------------  make these all pd.Timestamp, everywhere
        if ts is None or ts < start_time or ts > end_time:
            continue
        print(f'File to read: {f}')
        dfs.append(pd.read_csv(f))
    #join the data_fames together
    if dfs:
        #print(dfs[0])
        print(f'Read {len(dfs)} csv files')
        #join all the dataframes
        data_read = pd.concat(dfs, ignore_index=True)
        data_read = data_read.drop_duplicates()
        print(f'Length of concat NGFS  dataframe {len(data_read)}')
        if parse_times:
            data_read["pixel_date_time"] = pd.to_datetime(data_read["pixel_date_time"], utc=True, errors="coerce")
            data_read["acq_date_time"] = pd.to_datetime(data_read["acq_date_time"], utc=True, errors="coerce")
        data_read['acq_date'] = data_read['acq_date_time']   # <<<------------------------------------------------------ probably remove this
    else: 
        print('No data read, returning empty datframe')
        data_read = pd.DataFrame()

    print(data_read)
    return data_read

def make_geometry(feature,geometry_type = "point"):
    #constructs the geometry for a feature (row in csv) of NFS data
    #features is a row from a dataframe
    #"geometry_type = "point" or (geometry_type=="polygon")
    if ('latitude_c1' not in feature.keys()) and (geometry_type=="polygon"):
        print('Pixel corners not in data, will return point geometry')
        geometry_type = "point"
    if geometry_type == "point":
        gt = "Point"
        c = [feature['longitude'],feature['latitude']]
    else:
        gt = "Polygon"
        c = []
        # Fill corners array with latitude and longitude values
        #xcorners seem to be counterclockwise from NW corner
        for i in range(1,5):
            lat = feature[f'latitude_c{i}']
            lon = feature[f'longitude_c{i}']
            c.append([lat,lon])
        c.append(c[0]) #close on itself
    g = {
        "type" : gt,
        "coordinates" : c
    }
    return g


def make_geojson(df,geometry_type = "point"):
    #makes a geojson file from an NGFS dataframe
    if ('latitude_c1' not in df.keys()) and (geometry_type=="polygon"):
        print('Pixel corners not in data, will return point geometry')
        geometry_type = "point"

    #reset df.index to avoid problems
    df.reset_index(drop=True, inplace=True)

    #covert the time columns to something OK for geojson
    # Identify all datetime columns (including timezone-aware ones)
    date_cols = df.select_dtypes(include=['datetime64', 'datetimetz']).columns
    # Convert those columns to string
    df[date_cols] = df[date_cols].astype(str)


    #loop through the rows in the dataframe and make a feature from each
    features = []
    for i in df.index:
        f = df.loc[i].copy()
        prop = f.to_dict()
        g = make_geometry(f,geometry_type=geometry_type)
        feat = {
            "type" : "Feature",
            "properties" : prop,
            "geometry" : g
        }
        features.append(feat)
    #assemble features into a geojson file
    gj = {
        "type" : "FeatureCollection",
        "features": features
    }
    return gj

def load_ngfs_data(start_time = None,end_time = None,sat_type = 'goes',ngfs_cfg = None,df=pd.DataFrame()):
    #loads the data between two times
    #get file list
    if len(df) > 0:
        print('Updating DataFrame')
        satrt_time = df.acq_date_time.max()
        end_time = pd.Timestamp.utcnow()




if __name__ == "__main__":
    print('Testing the ngfs_ftp module')
    
    ngfs_cfg, wrfxpy_cfg = con_man.load_cfgs()  #in __main__ for testing

    today, csv_file = setup_csv_download(ngfs_cfg = ngfs_cfg)
    #print(today)
    #print(csv_file)
    data, csv_date_str = read_NGFS_csv_data(csv_file)
    print(data)
    print(f'Dataframe size: {len(data)}')
