import numpy as np
import logging
from utils import Dict
from geo.geo_utils import fill_categories

def print_stats(varname,v,unit):
    """
    print some statistics on a variable
    :param varname: string with variable name
    :param v: array with values
    """
    logging.info('%s %s shape %s nans %s  min %s max %s mean %s %s' 
        %(varname,repr(type(v)),repr(v.shape),np.sum(np.isnan(np.ma.compressed(v))),
        np.nanmin(v),np.nanmax(v),np.nanmean(v),unit))
    h = np.histogram(v)
    logging.info("%s n: %s %s" % (varname,h[0],unit))
    logging.info("%s bins: %s" % (varname,h[1],))

def interpolate2height_old(var,height,level):
    """
    Interpolate 3d variable to a given height
    :param var: the 3d array to be interpolated, 1st axis is vertical
    :param height: the 3d array of heights of the nodes on which v is defined
    :param level: the target height to interpolate to
    :return: interpolated value, or
    """
    maxlayer = var.shape[0]-1
    z = np.zeros([var.shape[1],var.shape[2]])
    r = z
    for i in range(0, var.shape[1]):
        for j in range(0, var.shape[2]):
            k = np.searchsorted(height[:,i,j],level)
            if k==0 or k>maxlayer:
                r[i,j] = np.nan
                #logging.error("Need height[0,%s,%s]=%s < level=%s <= height[%s,%s,%s]=%s" \
                #   % (i,j,height[0,i,j],level,maxlayer,i,j,height[maxlayer,i,j]))
                #return z
            else:
                # interpolate in the interval height[k-1,i,j] to height[k,i,j]
                r[i,j]=var[k-1,i,j]+(var[k,i,j]-var[k-1,i,j]) \
                    * (level - height[k-1,i,j])/(height[k,i,j] - height[k-1,i,j])
    return r

def interpolate2height(var,height,level):
    """
    Interpolate 3d variable to a given height
    :param var: the 3d array to be interpolated, 1st axis is vertical
    :param height: the 3d array of heights of the nodes on which v is defined
    :param level: the target height to interpolate to
    :return: interpolated value, or
    """
    ix, tx = index8height(height,level)
    r = np.zeros([var.shape[1],var.shape[2]])
    maxlayer = var.shape[0]-1
    for i in range(0, var.shape[1]):
        for j in range(0, var.shape[2]):
            k = ix[i,j]
            t = tx[i,j]
            if k==0 or k>maxlayer:
                #r[i,j] = np.nan
                r[i,j] = 0
            else:
                #r[i,j]=var[k,i,j]+(var[::k+1,i,j]-var[k,i,j])*tx[i,j] 
                r[i,j] = var[k,i,j]*(1.0-t) + var[k+1,i,j]*t 
    logging.info('interpolated to %s: min %s max %s' % (level,np.min(r),np.max(r)))
    return r

def integrate_ratio_to_level(d,t,ratio,height,level):
    """
    Interpolate 3d mixing ratio to a given level
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :param ratio: the 3d array to be integrated, 1st asix is is vertical (X/kg)
    :param height: the 3d array of heights of the nodes on which v is defined
    :param level: the target height to interpolate to (m)
    :return: vertically integrated mass, (x/m^2)
    """
    # logging.info('mixing ratio min %s max %s /m^2' % (np.min(ratio),np.max(ratio)))
    rho = density(d,t)      # air density (kg/m^3)
    dz = dz8w(d,t)      # vertical mesh steps (m)
    var = ratio * rho * dz
    ratio_int = sum_to_level(var,height,level)
    logging.info('integrated to %s: min %s max %s /m^2' % (level,np.min(ratio_int),np.max(ratio_int)))
    return ratio_int

def sum_to_level(var,height,level):
    """
    Interpolate 3d variable to a given height
    :param var: the 3d array to be integrated, 1st axis is vertical
    :param height: the 3d array of heights of the nodes on which v is defined
    :param level: the target height to interpolate to
    :return: interpolated value, or
    """
    ix, tx = index8height(height,level)
    r = np.zeros([var.shape[1],var.shape[2]])
    maxlayer = var.shape[0]-1
    for i in range(0, var.shape[1]):
        for j in range(0, var.shape[2]):
            k = ix[i,j]
            t = tx[i,j]
            if k<=0 :
                r[i,j] = var[k,i,j]*t 
            elif k>maxlayer:
                r[i,j] = np.sum(var[:,i,j],axis=0)
            else:
                #r[i,j]=var[k,i,j]+(var[k+1,i,j]-var[k,i,j])*tx[i,j] 
                r[i,j] = np.sum(var[0:k-1,i,j],axis=0)+var[k,i,j]*t 
    return r

def index8height(height,level):
    """
    Find index and fraction at given height
    :param height: the 3d array of heights of the nodes on which v is defined
    :param level: the target height to interpolate to
    :return: index of level as integer part and fractional part
    If level is <= smallest height or > largest g=height, return all zeros
    """
    maxlayer = height.shape[0]-1
    ix = np.zeros([height.shape[1],height.shape[2]],dtype=np.int_)
    tx = np.zeros([height.shape[1],height.shape[2]])
    for i in range(0, height.shape[1]):
        for j in range(0, height.shape[2]):
            k = np.searchsorted(height[:,i,j],level)
            if k == 0:
                ix[i,j] = -1 if level < height[0,i,j] else 0 # leaving tx[i,j]=0
            elif k < height.shape[0]:
                ix[i,j]=k-1    # integer partheight[:,i,j]
                # interpolation in the interval height[k-1,i,j] to height[k,i,j]
                tx[i,j]= (level - height[k-1,i,j])/(height[k,i,j] - height[k-1,i,j])
            else:
                ix[i,j]=height.shape[0]-1 #  max, leaving tx[i,j]=0
                #logging.warning("Need height[0,%s,%s]=%s < level=%s <= height[%s,%s,%s]=%s, got index k=%s" \
                #   % (i,j,height[0,i,j],level,maxlayer,i,j,height[maxlayer,i,j],k))
    logging.info('grid layer min %s max %s'  % ( np.min(ix), np.max(ix)))
    return ix,tx 

def pressure(d,t):
    """
    Compute pressure at mesh centers
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    """
    P = d.variables['P'][t,:,:,:] + d.variables['PB'][t,:,:,:]
    print_stats('pressure',P,'Pa')
    return P

def smoke_concentration(d,t):
    """
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :returns smoke ug/m^3
    """
    s=d.variables['tr17_1'][t,:,:,:]*density(d,t)
    print_stats('smoke concentration',s,'ug/m^3')
    return s

def pressure8w(d,t):
    """
    Compute pressure height at mesh cell bottoms a.k.a. w-points 
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    """
    ph = pressure(d,t)
    ph8w = ph 
    # average from 2nd layer up 
    ph8w[1:,:,:] = 0.5*(ph[1:,:,:] + ph[0:ph.shape[0]-1,:,:])
    return ph8w 

def u8p(d,t):
    """
    Compute horizontal wind u at mesh cell centers a.k.a. p-points
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    """
    u = d.variables['U'][t,:,:,:]
    return 0.5*(u[:,:,0:u.shape[2]-1]+u[:,:,1:])

def v8p(d,t):
    """
    Compute horizontal wind v at mesh cell centers a.k.a. p-points
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    """
    v = d.variables['V'][t,:,:,:]
    return 0.5*(v[:,0:v.shape[1]-1,:]+v[:,1:,:])

def w8p(d,t):
    """
    Compute vertical wind w at mesh cell centers a.k.a. p-points
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    """
    w = d.variables['W'][t,:,:,:]
    return 0.5*(w[0:w.shape[0]-1,:,:]+w[1:,:,:])

def height8w(d,t):
    """
    Compute height at mesh bottom a.k.a. w-points 
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    """
    ph = d.variables['PH'][t,:,:,:]  
    phb = d.variables['PHB'][t,:,:,:]
    return (phb + ph)/9.81 # geopotential height at W points

def height8w_terrain(d,t):
    """
    Compute height of mesh centers (p-points) above terrain
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    """
    h = height8w(d,t)
    terrain_height = h[0,:,:] 
    for i in range(0, h.shape[1]):
        for j in range(0, h.shape[2]):
            h[:,i,j] -= terrain_height[i,j]
    return h

def height8p(d,t):
    """
    Compute height of mesh centers (p-points)
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    """
    z8w = height8w(d,t)
    return 0.5*(z8w[0:z8w.shape[0]-1,:,:]+z8w[1:,:,:])

def height8p_terrain(d,t):
    """
    Compute height of mesh centers (p-points) above terrain
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    """
    z8w = height8w(d,t)
    h =  0.5*(z8w[0:z8w.shape[0]-1,:,:]+z8w[1:,:,:])
    for i in range(0, h.shape[1]):
        for j in range(0, h.shape[2]):
            h[:,i,j] -= z8w[0,i,j]
    return h

def dz8w(d,t):
    z8w = height8w(d,t)             # height of mesh bottom (m)
    return z8w[1:,:,:]-z8w[0:z8w.shape[0]-1,:,:] # mesh cell heights (m)

def hPa_to_m(p):
    """
    Compute pressure altitude
    :param p: pressure (hPa)
    :return: altitude (ft) 
    """
    # https://www.weather.gov/media/epz/wxcalc/pressureAltitude.pdf
    return (0.3048*145366.45)*(1 - (p/1013.25)**0.190284)

def density(d,t):
    """
    Get air density 
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :return: air density (kg/m^3)
    """
    P = pressure(d,t)   # dry air pressure (Pa )
    T = d.variables['T'][t,:,:,:] + d.variables['T00'][t]  # temperature (K)
    r_d = 287                       # specific gas constant (J/kg/K)
    rho = P/(r_d * T)               # dry air density  (kg/m^3)
    print_stats('temperature',T,'K')
    print_stats('density',rho,'kg/m^3')
    return rho

def cloud_to_level_hPa(d,t,level_hPa):
    """
    Integrate cloud water density from the ground to given pressure height
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :param level_hPa: pressure level
    :return: cloud water intensity to given pressure level (kg/m^2)
    """
    qcloud = d.variables['QCLOUD'][t,:,:,:] # cloud water mixing ratio (kg water/kg dry air)
    p8w = pressure8w(d,t) # pressure at cell bottoms (Pa)
    h8w_m = hPa_to_m(p8w*0.01)  # pressure height at cell bottoms (m)
    level_m = hPa_to_m(level_hPa) # desired pressure level (m)
    cloud_int = integrate_ratio_to_level(d,t,qcloud,h8w_m,level_m)
    # logging.info('integrated cloud water to pressure height %sm min %s max %s kg/m^2' % (level_m,np.min(cloud_int),np.max(cloud_int)))
    return cloud_int
    
def smoke_to_height_terrain(d,t,level):
    """
    Integrate smoke from the ground to given height above terrain
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :param level: height above terrain (m)
    :return: smoke integrated to given level (mg/m^2)
    """
    smoke= d.variables['tr17_1'][t,:,:,:] # smoke mixing ratio (ug/kg dry air)
    print_stats('tr17_1',smoke,'ug/kg')
    # logging.info('smoke tr17_1 min %s max %s ug/kg' % (np.min(smoke),np.max(smoke)))
    h_terrain = height8w_terrain(d,t)  # height above the terrain
    htw = h_terrain[0:h_terrain.shape[0]-1,:,:] # get rid of extra end stagger points at the top
    smoke_int = integrate_ratio_to_level(d,t,smoke,htw,level) 
    logging.info('integrated smoke to %sm min %s max %s ug/m^2' % (level,np.min(smoke_int),np.max(smoke_int)))
    return smoke_int
    
def transform_goes(d):
    """
    Transform GOES mask data into similar to MODIS-VIIRS AF mask
    :param d: open NetCDF4 dataset
    """
    fill = Dict({(150., 151., 152., 153.): 3, (100.): 5, (15., 35.): 7, (14., 34.): 8, (10., 11., 12., 13., 30., 31., 32., 33.): 9})
    array = d.variables['Mask'][:]
    return fill_categories(array, fill)

def fire_front(d,t,var):
    """
    Plot variable only at the fire front
    """
    dt = 15.*60.
    v = d.variables[var][t,:,:]
    lfn = d.variables['LFN'][t,:,:]
    m = np.logical_and(lfn > 60.,lfn <= dt)
    v[~m] = np.nan
    return v

def calculate_vp(q, P):
    """
    Calculate Vapor Pressure from Specific Humidity [kg/kg] and Pressure [Pa]
    :param q: specific humidity [kg/kg]
    :param P: pressure [Pa]
    :return pw: vapor pressure [Pa]
    """
    # Molecular weight of water (18.02 g/mol) to molecular weight of dry air (28.97 g/mol)
    epsilon = 0.622 
    # Vapor pressure [Pa]
    pw = q * P / (epsilon + (1-epsilon) * q)
    return pw

def calculate_svp(T):
    """
    Calculate Saturation Vapor Pressure [Pa] from Temperature [K]
    Following Murphy and Koop, Q.J.R. Meteorol. Soc (2005) 131 1539-1565 eq. (10)
    :param T: temperature [K]
    :return pws: saturation vapor pressure [Pa]
    """ 
    # Saturation vapor pressure [Pa] 
    pws = np.exp( 54.842763 - 6763.22/T - 4.210 * np.log(T) + 0.000367 * T +
        np.tanh(0.0415*(T - 218.8)) * (53.878 - 1331.22 / T - 9.44523 * np.log(T) + 0.014025*T))
    return pws

def calculate_rh(T, q, P):
    """
    Calculate the relative humidity
    ES=610.78*exp(17.269*(T-273.161)/(T-35.861))
    QRS=0.622*ES/(P-0.378*ES)
    RH = Q/QRS
    From function rh_from_q from Adam Kochanski 
    :param T: temperature [K]
    :param q: specific humidity [kg/kg]
    :param P: pressure [Pa]
    :return rh: relative humidity [1]
    """
    # Vapor pressure [Pa]
    pw = calculate_vp(q, P) 
    # Saturation vapor pressure [Pa]
    pws = calculate_svp(T)
    # Realtive humidity [1]
    rh = pw / pws
    return np.clip(rh, 0, 1)

def calculate_ffwi(T, q, P, u, v):
    """
    Calculate the Fosberg Index
    :param T: temperature [K]
    :param q: specific humidity [kg/kg]
    :param P: pressure [Pa]
    :param u: u-component of wind [m/s]
    :param v: v-component of wind [m/s] 
    :return: Fosberg Index (unitless)
    """
    # Wind speed m/s -> mph
    u *= 2.23694
    v *= 2.23694
    ws = np.sqrt(u**2 + v**2) # Wind speed (mph)
    # Temperature K -> F
    T_f = ((T - 273.15) * 1.8) + 32
    # Find relative humidity [%]
    rh = calculate_rh(T, q, P) * 100
    # Initialize equilibrium moisture content
    emc = np.zeros(rh.shape)
    # Case 1: RH >= 50
    rh_mask = rh >= 50
    emc[rh_mask] = 21.0606 + 0.005565 * rh[rh_mask]**2 - 0.00035 * rh[rh_mask] * T_f[rh_mask] - 0.483199 * rh[rh_mask] # Equilibrium Moisture Content (%)
    # Case 2: 10 <= RH < 50
    rh_mask = np.logical_and(rh >= 10, rh < 50)
    emc[rh_mask] = 2.22749 + 0.160107 * rh[rh_mask] - 0.014784 * T_f[rh_mask] # Equilibrium Moisture Content (%)
    # Case 3: RH < 10
    rh_mask = rh < 10
    emc[rh_mask] = 0.03229 + 0.281073 * rh[rh_mask] - 0.000578 * rh[rh_mask] * T_f[rh_mask] # Equilibrium Moisture Content (%)
    # Find eta
    eta = 1 - 2 * (emc/30) + 1.5 * (emc / 30)**2 - 0.5 * (emc / 30)**3
    # Find Fosberg Index
    ffwi = (eta * np.sqrt(1 + ws**2)) / 0.3002 # Fosberg Index Value (unitless)
    return ffwi

def ffwi(d, t):
    """
    Interface to calculate the FFWI from WRF data
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :return: Fosberg Index (unitless)
    """
    # Get the needed variables
    q2 = d.variables['Q2'][t,:,:]        # Specific Humidity (kg/kg)
    psfc = d.variables['PSFC'][t,:,:]    # Surface Pressure (Pa)
    t2 = d.variables['T2'][t,:,:]        # Temperature (K)
    u10 = d.variables['U10'][t,:,:]      # U-Component of the Wind at 10m (m/s)
    v10 = d.variables['V10'][t,:,:]      # V-Component of the Wind at 10m (m/s)
    # Calculate FFWI
    ffwi = calculate_ffwi(t2, q2, psfc, u10, v10)
    return ffwi

def calculate_hdw(T, q, P, u, v):
    """
    Calculate the Hot Dry and Windy (HDW) Index
    :param T: temperature [K]
    :param q: specific humidity [kg/kg]
    :param P: pressure [Pa]
    :param u: u-component of wind [m/s]
    :param v: v-component of wind [m/s] 
    :return: HDW Index (hPa m s-1, but units should be ignored)
    """
    # Find vapor pressure and saturated vapor pressure in Pa
    pw = calculate_vp(q, P)
    pws = calculate_svp(T)
    # Find vapor pressure deficit (hPa)
    vpd = (pws - pw) / 100
    # Wind speed (m/s)
    ws = np.sqrt(u**2 + v**2) 
    # Find Hot, Dry, & Windy Index (hPa m s-1) 
    hdw = ws * vpd
    return hdw

def hdw(d,t):
    """
    Interface to calculate the HDW from WRF data
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :return: HDW Index (hPa m s-1, but units should be ignored)
    """
    # Get the needed variables
    q2 = d.variables['Q2'][t,:,:]        # Specific Humidity (kg/kg)
    psfc = d.variables['PSFC'][t,:,:]    # Surface Pressure (Pa)
    t2 = d.variables['T2'][t,:,:]        # Temperature (K)
    u10 = d.variables['U10'][t,:,:]      # U-Component of the Wind at 10m (m/s)
    v10 = d.variables['V10'][t,:,:]      # V-Component of the Wind at 10m (m/s)
    hdw = calculate_hdw(t2, q2, psfc, u10, v10)
    return hdw