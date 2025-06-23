import numpy as np
import logging

########### FWI TOOLS ###########

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
    :return Pws: saturation vapor pressure [Pa]
    """ 
    # Saturation vapor pressure [Pa] 
    Pws = np.exp( 54.842763 - 6763.22/T - 4.210 * np.log(T) + 0.000367 * T +
        np.tanh(0.0415*(T - 218.8)) * (53.878 - 1331.22 / T - 9.44523 * np.log(T) + 0.014025*T))
    return Pws

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
    Pw = calculate_vp(q, P) 
    # Saturation vapor pressure [Pa]
    Pws = calculate_svp(T)
    # Realtive humidity [1]
    rh = Pw / Pws
    return np.clip(rh, 0, 1)

def calculate_emc(T, rh):
    """
    Calculate equilibrium moisture content (emc) from temperature [K] and relative humidity [%]
    :param T: temperature [K]
    :param rh: relative humidity [1]
    """
    # Temperature K -> F
    T_f = ((T - 273.15) * 1.8) + 32
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
    return emc

def calculate_eta(fmc, fmce=30):
    """
    Calculate moisture damping coefficient (eta)
    :param fmc: fuel moisture [%]
    :param fmce: fuel moisture of extinction [%]
    :return eta: eta value
    """ 
    # Find moisture ratio
    ratio = fmc/fmce
    # Find moisture damping coefficient
    eta = 1 - 2.59 * ratio + 5.11 * ratio**2 - 3.52 * ratio**3
    return np.clip(eta, 0, 1)

def calculate_air_temp(theta, P):
    """
    Calculate air temperature from potential temperature based on Poisson's equation
    :param theta: potential temperature [K]
    :param P: pressure [Pa]
    :return T: air temperature [K]
    """
    P0 = 100000.  # reference pressure [Pa]
    rd = 287.04   # gas constant for dry air in [(J/kg)K]
    cp = 1004.    # specific heta at constant pressure [(J/kg)K]
    T = theta * ((P / P0) ** (rd/cp))
    return T

def calculate_dew_point(q, P):
    """
    Calculate dew point temperature using Magnus-Tetens equation 
    :param q: potential temperature [K]
    :param P: pressure [Pa]
    :return Td: dew points temperature [C]
    """
    # Find vapor pressure in hPa
    Pw = calculate_vp(q, P) / 100.
    # Compute dew point temperature in Celsius
    Td = (243.5 * np.log(Pw / 6.112)) / (17.67 - np.log(Pw / 6.112))
    return Td

########### FWI CALCULATIONS ###########
    
def calculate_ffwi(ws, eta):
    """
    Calculate the Fosberg Index
    :param ws: wind speed [mph]
    :param eta: moisture damping coefficient 
    :return: Fosberg Index (unitless)
    """
    # Fosberg Index Value (unitless)
    ffwi = (eta * np.sqrt(1 + ws**2)) / 0.3002 
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
    Pw = calculate_vp(q, P)
    Pws = calculate_svp(T)
    # Find vapor pressure deficit (kPa)
    vpd = (Pws - Pw) / 1000
    # Wind speed (m/s)
    ws = np.sqrt(u**2 + v**2) 
    # Find Hot, Dry, & Windy Index (kPa m s-1) 
    hdw = ws * vpd
    return hdw

def calculate_lfp(T, Td, ws):
    """
    Calculate the Large Fire Potential (LFP) Index
    :param T: temperature [K]
    :param Td: dew point temperature [K]
    :param ws: sustained wind speed [m/s] 
    :return: LFM Index (unitless)
    """
    # Temperature conversions
    T = (T - 273.15) * 9/5 + 32   # K to F
    Td = (Td - 273.15) * 9/5 + 32 # K to F
    # Compute dew point depression in F
    dp = (T - Td) 
    # Find Large Fire Potential Index (unitless) 
    lfp = 0.001 * np.power(ws, 2) * dp
    return lfp

def calculate_haines(t_950, t_850, t_700, t_500, t_dew850, t_dew700, height):
    """
    Calculate the Haines Index
    :param t_950: Temprature at 950 atm [K]
    :param t_850: Temprature at 850 atm [K]
    :param t_700: Temprature at 700 atm [K]
    :param t_500: Temprature at 500 atm [K]
    :param t_dew850: Dew Point Temprature at 850hPa [K]
    :param t_dew700: Dew Point Temprature at 700hPa [K]
    :param height: Height above sea level [m]
    """
    
    def haines_low(t_950, t_850, t_dew850):
        '''Low Elevation Haines Index'''
        A = np.where(
            t_950 - t_850 <= 3, 
            1, np.where(
                t_950 - t_850 <= 7, 
                2, 3
            )
        ) # stability score
        B = np.where(
            t_850 - t_dew850 <= 5, 
            1, np.where(
                t_850 - t_dew850 <= 9, 
                2, 3
            )
        ) # moisture score
        return A + B
    
    def haines_mid(t_850, t_700, t_dew850):
        '''Mid Elevation Haines Index'''
        A = np.where(
            t_850 - t_700 <= 5, 
            1, np.where(
                t_850 - t_700 <= 10, 
                2, 3
            )
        ) # stability score 
        B = np.where(
            t_850 - t_dew850 <= 5, 
            1, np.where(
                t_850 - t_dew850 <= 12, 
                2, 3
            )
        ) # moisture score
        return A + B
    
    def haines_high(t_700, t_500, t_dew700):
        '''High Elevation Haines Index'''
        A = np.where(
            t_700 - t_500 <= 17, 
            1, np.where(
                t_700 - t_500 <= 21, 
                2, 3
            )
        ) # stability score 
        B = np.where(
            t_700 - t_dew700 <= 14, 
            1, np.where(
                t_700 - t_dew700 <= 20, 
                2, 3
            )
        ) # moisture score
        return A + B
    
    # Determine Haines category based on elevation
    haines_index = np.where(
        height < 304.8,  # <1000 ft
        haines_low(t_950, t_850, t_dew850),
        np.where(
            height < 914.4,  # <3000 ft
            haines_mid(t_850, t_700, t_dew850),
            haines_high(t_700, t_500, t_dew700),
        )
    )
    return haines_index
 
########### FWI INTERFACES ###########

def FFWI(d, t):
    """
    Interface to calculate the Forsberg Fire Weather Index (FFWI) from WRF data
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :return: FFWI (unitless)
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

def HDW(d, t):
    """
    Interface to calculate the HDW from WRF data
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :return: HDW (hPa m s-1, but units should be ignored)
    """
    # Get the needed variables
    q2 = d.variables['Q2'][t,:,:]        # Specific Humidity (kg/kg)
    psfc = d.variables['PSFC'][t,:,:]    # Surface Pressure (Pa)
    t2 = d.variables['T2'][t,:,:]        # Temperature (K)
    u10 = d.variables['U10'][t,:,:]      # U-Component of the Wind at 10m (m/s)
    v10 = d.variables['V10'][t,:,:]      # V-Component of the Wind at 10m (m/s)
    # Wind speed mph
    ws = np.sqrt(u10**2 + v10**2) * 2.23694
    # Find relative humidity [%]
    rh = calculate_rh(t2, q2, psfc) * 100
    # Calculate equilibrium moisture content (EMC)
    emc = calculate_emc(t2, rh)
    # Calculate moisture damping coefficient
    eta = calculate_eta(emc)
    # Calculate HDW
    hdw = calculate_hdw(ws, eta)
    return hdw

def LFP(d, t):
    """
    Interface to calculate the Large Fire Potential (LFP) from WRF data
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :return: LFP (unitless)
    """ 
    # Get the needed variables
    q2 = d.variables['Q2'][t,:,:]        # Specific Humidity (kg/kg)
    psfc = d.variables['PSFC'][t,:,:]    # Surface Pressure (Pa)
    t2 = d.variables['T2'][t,:,:]        # Temperature (K)
    u10 = d.variables['U10'][t,:,:]      # U-Component of the Wind at 10m (m/s)
    v10 = d.variables['V10'][t,:,:]      # V-Component of the Wind at 10m (m/s)
    # Calculate dew point temperature [Celsius]
    Td = calculate_dew_point(q, P)
    # Conversion Celsius to Kelvin
    Td = Td + 273.15
    # Wind speed (m/s)
    ws = np.sqrt(u**2 + v**2) 
    # Calculate LFP
    lfp = calculate_lfp(T, Td, ws)
    return lfp

def Haines(d, t):
    """
    Interface to calculate the Haines Index from WRF data
    :param d: open NetCDF4 dataset
    :param t: number of timestep
    :return: Haines (unitless)
    """
    def nearest_index(array, value):
        '''Convert hPa to Pa and find index'''
        return (np.abs(array - value * 100)).argmin()
    # Get the needed variables
    w = d.variables['QVAPOR'][t,:,:]                            # Water vapor mixing ratio (kg/kg)
    P = d.variables['P'][t,:,:] + d.variables['PB'][t,:,:]      # Full pressure (Pa)
    theta = d.variables['T'][t,:,:] + 300                       # Potential Temperature (K)
    gph = d.variables['PH'][t,:,:] + d.variables['PHB'][t,:,:]  # Geopotential height (m)
    # Calculate air temperature in K
    T = calculate_air_temp(theta, P)
    # Calculate specific humidity from water vapor mixing ratio
    q = w/(1+w)
    # Calculate dew point temperature in K
    Td = calculate_dew_point(q, P) + 273.15
    # Extract values at required pressure levels
    T_950 = T[nearest_index(P, 950)]
    T_850 = T[nearest_index(P, 850)]
    T_700 = T[nearest_index(P, 700)]
    T_500 = T[nearest_index(P, 500)]
    Td_850 = Td[nearest_index(P, 850)]
    Td_700 = Td[nearest_index(P, 700)]
    # Extract surface height (first level of geopotential height)
    surface_height = gph[0] / 9.81  # Convert geopotential height to meters
    # Calculate Haines index
    Haines = calculate_haines(T_950, T_850, T_700, T_500, Td_850, Td_700, surface_height)
    return Haines
