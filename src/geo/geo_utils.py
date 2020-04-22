from scipy import spatial
import numpy as np

def deg2str(deg,islat):
    """
    Convert decimal degree coordinate to degress minutes seconds.xx E/W/S/N string
    
    :param deg: degree coordinate
    :param islat: boolean if coordinate is from latitude or not
    :return: string format coordinate 
    """
    deg = round(deg*3600.0,2)/3600.0
    d = int(deg)
    m = int((deg - d) * 60)
    s = (deg - d - m/60.0) * 3600.00
    z= round(s, 2)
    NSEW = [['E', 'W'], ['N', 'S']]
    return '%3.0fd%2.0f\'%5.2f\"%s' % (abs(d), abs(m), abs(z), NSEW[islat][d<0])

def coord2str(lon_deg,lat_deg):
    """
    Convert longitude and latitude degrees to degress minutes seconds.xx E/W/S/N string
    
    :param lon_deg: longitude degree coordinate
    :param lat_deg: latotide degree coordinate
    :return: string format coordinates
    """
    # compute longitude
    lon_str = deg2str(lon_deg,False)
    # compute latitude
    lat_str = deg2str(lat_deg,True)
    return '%s, %s' % (lon_str,lat_str)

#def fill_novalues(array,,fill_value):

def wrf_gc_distance(lon1,lon2,lat1,lat2,plot=False):
    """
    Compute great circle distance in meters from lon-lat coordinates
    
    :param lon1,lat1: lon-lat coordinates to compute distance from
    :param lon2,lat2: lon-lat coordinates to compute distance to
    :param plot: optional, visualize great circle distance in meters
    :return: great circle distance in meters
    """
    # Earth radius in WRF
    R = 6370000.
    # Compute coordinates in radians
    to_rads = np.pi/180.
    lon1 = np.array(lon1*to_rads)
    lon2 = np.array(lon2*to_rads)
    lat1 = np.array(lat1*to_rads)
    lat2 = np.array(lat2*to_rads)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)
    sdlon = np.sin(abs(lon2-lon1))
    cdlon = np.cos(abs(lon2-lon1))
    
    # Central angle in radians times radius in meters
    diff = 2*R*np.arctan(np.sqrt((clat2*sdlon)**2 + (clat1*slat2-slat1*clat2*cdlon)**2)/(slat1*slat2+clat1*clat2*cdlon))
    print('Great-circle distance --> min_distance=%fm -- mean_distance=%fm -- max_distance=%fm' % (diff.min(),diff.mean(),diff.max()))
    
    if plot:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import axes3d
        from matplotlib import cm
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        fig.suptitle("Great-circle distance in meters")
        ax.plot_surface(lon1,lat1,diff,cmap=cm.coolwarm)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_zlabel("Difference")
        plt.show()

    return diff 
