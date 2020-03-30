from scipy import spatial

def deg2str(deg,islat):
    # convert decimal degrees to degress minutes seconds.xx E/W/S/N string
    deg = round(deg*3600.0,2)/3600.0
    d = int(deg)
    m = int((deg - d) * 60)
    s = (deg - d - m/60.0) * 3600.00
    z= round(s, 2)
    NSEW = [['E', 'W'], ['N', 'S']]
    return '%3.0fd%2.0f\'%5.2f\"%s' % (abs(d), abs(m), abs(z), NSEW[islat][d<0])

def coord2str(lon_deg,lat_deg):
    # convert lon-lat coordinates to degress minutes seconds.xx E/W/S/N strings
    # compute longitude
    lon_str = deg2str(lon_deg,False)
    # compute latitude
    lat_str = deg2str(lat_deg,True)
    return '%s, %s' % (lon_str,lat_str)

#def fill_novalues(array,,fill_value):
