# reads ngfs csv file and lists incidents and starts a simulation from user selection
# borrows heavily from simple_forecast.py
# usage python ngfs_start.py /path/to/ngsf_csv_file.csv

import pandas as pd
import csv
import sys
import json
from datetime import timedelta, datetime
sys.path.insert(1, 'src/')
import utils
import simple_forecast as sf

class ngfs_incident():
   def __init__(self,name):
      self.name = name
   #add a pandas dataframe subset from csv
   #def add_process_date(self,df):
   #   self.process_date = 
   def add_incident_id_string(self,id_string):
      self.incident_id_string = id_string
   def add_detections(self,df):
      self.data = df
   def add_incident_start_time(self):
      if self.data.shape[0] > 0:
         self.incident_start_time = self.data.iloc[0,3]
      else:
         self.incident_start_time = 'Unset start time'
   #def add_ignition(self):
      



print()
print()
print('Starting ngfs script')

base_cfg = sf.questionnaire()

#read the data and filter
csv_file = sys.argv[1]
print('Loading the file : ',csv_file)
#columns with critical time information
time_cols = ['incident_start_time','observation_time','initial_observation_time']
data = pd.read_csv(csv_file, parse_dates=time_cols)
#make sure the times are UTC
data[time_cols[0]] = pd.DatetimeIndex(pd.to_datetime(data[time_cols[0]])).tz_localize('UTC')
#assumes standard  filename for all detections in day
csv_date_str = csv_file[-18:-8]
print('Detection data from ',csv_date_str)
csv_year = int(csv_date_str[0:4])
csv_month = int(csv_date_str[5:7])
csv_day = int(csv_date_str[8:10])
csv_timestamp = pd.Timestamp(year=csv_year,month=csv_month,day=csv_day)

#filter out "null incidents"
data = data.dropna(subset=['incident_name'])
print(data)

#print('Column headers from csv file')
for i in data.columns:
   print(i)

print('Incidents in csv file')
incident_names = data['incident_name'].unique()
num_incidents = incident_names.shape[0]
print('Number of incidents in csv file: ',num_incidents)
print(incident_names)


      
#array of incident objects
incidents  = [ngfs_incident(incident_names[i]) for i in range(num_incidents)]

for i in range(num_incidents):
#    print(incident_names[i])
#    incidents[i].name = incident_names[i]
    print(incidents[i].name)
    incident_subset = data[data.incident_name == incident_names[i]]
    #sort by observation time	
    incident_subset = incident_subset.sort_values(by='observation_time')
    #incident_subset_idx = data[data.incident_name == incident_names[i]].index
    det_count = incident_subset.shape[0]
    print('      ',det_count,'  detections') 
    incidents[i].add_detections(incident_subset)
    #print(incidents[i].data)
    incidents[i].add_incident_start_time()
    incidents[i].add_incident_id_string(incident_subset['incident_id_string'][incident_subset.index[0]])
    print('       Incident start time : ', incidents[i].incident_start_time)
    print('       Incident start day : ', incidents[i].incident_start_time.day)
    if (incidents[i].incident_start_time.day == csv_day and incidents[i].incident_start_time.month == csv_month):
       print('       New Incident From Today')
       print('       Finding start time and location')
       idx = incident_subset.index[0]
       ign_latlon = incident_subset.lat[idx], incident_subset.lon[idx]
       ign_utc = incident_subset.observation_time[idx]
       start_utc = utils.round_time_to_hour(ign_utc - timedelta(minutes=30))
       end_utc = start_utc + timedelta(hours=5)
       time_utc = utils.utc_to_esmf(ign_utc)
       print('       location: ',ign_latlon)
       print('       utc ignition time', time_utc)

       cfg = base_cfg
       gc_string  = (incident_subset['incident_name'][idx]+'_'+utils.utc_to_esmf(start_utc)+
                                         '_'+incident_subset['incident_id_string'][idx][1:-1])
       cfg['grid_code'] = gc_string.replace(' ','_')
       cfg['domains']['1']['center_latlon'] = ign_latlon
       cfg['domains']['1']['truelats'] = (ign_latlon[0], ign_latlon[0])
       cfg['domains']['1']['stand_lon'] = ign_latlon[1]
       ignitions = cfg['ignitions']
       ign_dur = ignitions['1'][0]['duration_s']
       cfg['ignitions'] = { '1' : [ { 'time_utc' : time_utc,
                                   'duration_s' : ign_dur,
                                   'latlon' : ign_latlon } ] }
       cfg['start_utc'] = utils.utc_to_esmf(start_utc)
       cfg['end_utc'] = utils.utc_to_esmf(end_utc)
       #print(cfg)

       filename = 'jobs/' + cfg['grid_code'] + '.json'
#       json.dump(cfg, open(filename, 'w'), indent=4, separators=(',', ': '))
       print(('INT to start the simulation, execute ./forecast.sh %s' % filename))
       
    
