import json
import sys
import csv
import math


def read_csv(csv_path):
    """
    Read flat csv file, first row are headers
    Return file contents as json {1:{header1:value1,...},2:}
    """

    print('Reading file %s\n' % csv_path)
    with open(csv_path,'rU') as csv_file:
        reader = csv.reader(csv_file, delimiter=',', quotechar='"')
        row_num=0
        out={}
        for row in reader:
            if row_num == 0: 
                head=row;
            else:
                out[row_num]={}
                for header, value in zip (head, row):
                     out[row_num][header]=value
            row_num = row_num + 1

        return(out)

def copyto(partial_path,kml):
    with open(partial_path,'r') as part:
        for line in part:
            kml.write(line)

def json2kml(s,kml_path):

    frp_style={-1:'modis_frp_no_data',40:'modis_frp_gt_400'}
    for i in range(0,40):
        frp_style[i]='modis_frp_%s_to_%s' % (i*10, i*10+10)
 
    with open(kml_path,'w') as kml:

        copyto('src/vis/partial1.kml',kml)

        r = 6378   # Earth radius
        km_lat = 180/(math.pi*r)  # 1km in degrees latitude
        sq_size_km = 1.0   # square size

        type = {'AF':'Active Fires','FRP':'Fire Radiative Power'}
        #type = {'AF':'Active Fires'}

        for t in type:

            kml.write('<Folder>\n')
            kml.write('<name>%s</name>\n' % type[t])

            if t=='FRP':
                copyto('src/vis/partial2.kml',kml)
     

            for p in s:
                latitude=s[p]['latitude']
                longitude=s[p]['longitude']
                acq_date=s[p]['acq_date']
                acq_time=s[p]['acq_time']
                satellite=s[p].get('satellite','Not available')
                instrument=s[p].get('instrument','Not available')
                confidence=s[p].get('confidence','Not available')
                frp=s[p].get('frp',0)
    
                latitude=float(latitude)
                longitude=float(longitude)
                confidence=float(confidence)
                frp=float(frp)
                timestamp=acq_date + 'T' + acq_time[0:2] + ':' + acq_time[2:4] + '-06:00'
                timedescr=acq_date + ' ' + acq_time[0:2] + ':' + acq_time[2:4] + ' UTC'
    
                print([longitude,latitude,acq_date,acq_time,satellite,instrument,confidence,frp])
    
                kml.write('<Placemark>\n<name>Fire detection square</name>\n')
                kml.write('<description>\nlongitude:  %s\n' % longitude 
                                      +  'latitude:   %s\n' % latitude
                                      +  'time:       %s\n' % timedescr
                                      +  'satellite:  %s\n' % satellite
                                      +  'instrument: %s\n' % instrument
                                      +  'confidence: %s\n' % confidence
                                      +  'FRP:        %s\n' % frp 
                        + '</description>\n')
                kml.write('<TimeStamp><when>%s</when></TimeStamp>\n' % timestamp)
                if t == 'AF':
                    if confidence < 30:
                        kml.write('<styleUrl> modis_conf_low </styleUrl>\n')
                    elif confidence < 80: 
                        kml.write('<styleUrl> modis_conf_med </styleUrl>\n')
                    else:
                        kml.write('<styleUrl> modis_conf_high </styleUrl>\n')
                elif t=='FRP':
                    frpx = min(40,math.ceil(frp/10.)-1)
                    kml.write('<styleUrl> %s </styleUrl>\n' % frp_style[frpx] )

                kml.write('<Polygon>\n<outerBoundaryIs>\n<LinearRing>\n<coordinates>\n')
    
    	        km_lon=km_lat/math.cos(latitude*math.pi/180)  # 1 km in longitude
                print km_lat

                if 'track' in s[p].keys():
                    sq_track_size_km=float(s[p]['track'])
                    sq2_lat=km_lat * sq_track_size_km/2
                else:
                    sq2_lat=km_lat * sq_size_km/2

                if 'scan' in s[p].keys():
                    sq_scan_size_km=float(s[p]['scan'])
                    sq2_lon=km_lon * sq_scan_size_km/2
                else:
                    sq2_lon=km_lon * sq_size_km/2

                kml.write('%s,%s,0\n' % (longitude - sq2_lon, latitude-sq2_lat))
                kml.write('%s,%s,0\n' % (longitude - sq2_lon, latitude+sq2_lat))
                kml.write('%s,%s,0\n' % (longitude + sq2_lon, latitude+sq2_lat))
                kml.write('%s,%s,0\n' % (longitude + sq2_lon, latitude-sq2_lat))
                kml.write('%s,%s,0\n' % (longitude - sq2_lon, latitude-sq2_lat))
    
                kml.write('</coordinates>\n</LinearRing>\n</outerBoundaryIs>\n</Polygon>\n</Placemark>\n')
    
            kml.write('</Folder>\n')

        kml.write('</Document>\n</kml>\n')
    
    print('Created file %s\n' % kml_path)
            
if __name__ == '__main__':


    if len(sys.argv) < 3 or not sys.argv[2]:
        print('usage: csv2kml.sh csv_file_path kml_file_path')
        #print('usage: python csv2kml.py csv_file_path kml_file_path')
        sys.exit(1)

    csv_file_path=sys.argv[1]
    kml_file_path=sys.argv[2]

    s = read_csv(csv_file_path) 
    json2kml(s,kml_file_path)

    exit(0) 
