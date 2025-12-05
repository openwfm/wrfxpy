import numpy as np
import f90nml
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

def get_namelist(namelist_fire_path):
    """
    Read namelist
    :param namelist_fire_path: path to namelist file
    :return nml: namelist object with parsed information
    """
    nml = f90nml.read(namelist_fire_path)
    return nml

def calculate_weighting_factors(nml):
    """
    Calculate fuel moisture weighting factors
    :param nml: namelist object with parsed information
    :return weights: weights in a dictionary 
    """
    ## 1) Reading all the parameters
    # Fuel density
    fueldens = np.array(nml["fuel_categories"]["fueldens"], dtype=float)
    # Surface area to volume ratio
    savr_1 = np.array(nml["fuel_categories"]["savr_gc01"], dtype=float)
    savr_10 = np.array(nml["fuel_categories"]["savr_gc02"], dtype=float)
    savr_100 = np.array(nml["fuel_categories"]["savr_gc03"], dtype=float)
    savr_herba = np.array(nml["fuel_categories"]["savr_gc05"], dtype=float)
    savr_woody = np.array(nml["fuel_categories"]["savr_gc06"], dtype=float)
    # Fuel load
    fuelload_1 = np.array(nml["fuel_categories"]["fuelload_gc01"], dtype=float)
    fuelload_10 = np.array(nml["fuel_categories"]["fuelload_gc02"], dtype=float)
    fuelload_100 = np.array(nml["fuel_categories"]["fuelload_gc03"], dtype=float)
    fuelload_herba = np.array(nml["fuel_categories"]["fuelload_gc05"], dtype=float) 
    fuelload_woody = np.array(nml["fuel_categories"]["fuelload_gc06"], dtype=float)
    ## 2) Compute surface area per unit of each size class within each category
    A_1      = savr_1 * fuelload_1 / fueldens
    A_10     = savr_10 * fuelload_10 / fueldens
    A_100    = savr_100 * fuelload_100 / fueldens
    A_herba  = savr_herba * fuelload_herba / fueldens
    A_woody  = savr_woody * fuelload_woody / fueldens
    # mean total surface area per unit fuel cell of the dead and live categories
    A_dead       = A_1 + A_10 + A_100
    A_live       = A_herba + A_woody
    # mean total surface area of the fuel
    AT = A_dead + A_live 
    ## 3) Calculate weighting factors of each size class within each category
    f_1 = np.zeros_like(A_1, dtype=float)
    f_10 = np.zeros_like(A_1, dtype=float)
    f_100 = np.zeros_like(A_1, dtype=float)
    f_herba = np.zeros_like(A_1, dtype=float)
    f_woody = np.zeros_like(A_1, dtype=float)
    f_dead = np.zeros_like(A_1, dtype=float)
    f_live = np.zeros_like(A_1, dtype=float)
    mask_dead = A_dead != 0
    mask_live = A_live != 0
    f_1[mask_dead] = A_1[mask_dead] / A_dead[mask_dead]
    f_10[mask_dead] = A_10[mask_dead] / A_dead[mask_dead]
    f_100[mask_dead] = A_100[mask_dead] / A_dead[mask_dead]
    f_herba[mask_live] = A_herba[mask_live] / A_live[mask_live]
    f_woody[mask_live] = A_woody[mask_live] / A_live[mask_live]
    # weighting factors of the dead and live categories
    mask_tot = AT != 0
    f_dead[mask_tot] = A_dead[mask_tot] / AT[mask_tot]
    f_live[mask_tot]  = A_live[mask_tot] / AT[mask_tot]
    ## 4) Dictionary with results
    weights = {
        "f_1": f_1, "f_10": f_10, "f_100": f_100,
        "f_herba": f_herba, "f_woody": f_woody,
        "f_dead": f_dead, "f_live": f_live
    }
    return weights

def calculate_eta_effective(
    fmc_1, fmc_10, fmc_100, fmc_herba, fmc_woody,
    n_fuel_cat, weights, ratios
):
    """
    Compute a single effective moisture damping coefficient (eta_eff)
    for a given fuel category, combining dead and live fuels.

    Returns
    -------
    eta_eff : float or array
        Effective moisture damping coefficient for the fuel bed.
    eta_dead : float or array
        Dead-fuel moisture damping coefficient.
    eta_live : float or array
        Live-fuel moisture damping coefficient.
    """
    n_fuel_cat = np.asarray(n_fuel_cat, dtype=int)
    idx = n_fuel_cat[None, ...] - 1  # same shape as grid
    # 1) Dead and live FMC (surface-area weighted within group)
    f_1     = weights["f_1"][idx]
    f_10    = weights["f_10"][idx]
    f_100   = weights["f_100"][idx]
    f_herba = weights["f_herba"][idx]
    f_woody = weights["f_woody"][idx]
    fmc_dead = f_1 * fmc_1 + f_10 * fmc_10 + f_100 * fmc_100
    fmc_live = f_herba * fmc_herba + f_woody * fmc_woody
    # 2) Dead & live moisture extinction
    fmce_dead = ratios["fmce_dead"][idx]
    fmce_live = calculate_fmce_live(fmc_1, fmc_10, fmc_100, ratios, n_fuel_cat)
    # 3) Group etas
    eta_dead = calculate_eta(fmc_dead, fmce_dead)
    eta_live = calculate_eta(fmc_live, fmce_live)
    # 4) Group weights (dead vs live) for this fuel category
    f_dead = weights["f_dead"][idx]
    f_live = weights["f_live"][idx]      
    numerator = (
        eta_dead * f_dead + 
        eta_live * f_live
    )
    denom = f_dead + f_live
    mask_denom = denom > 0
    eta_eff = np.where(
        mask_denom,
        numerator / denom,
        0.0
    )
    return eta_eff, eta_dead, eta_live, fmce_live

def calculate_effective_loading_ratios(nml):
    fc = nml["fuel_categories"]
    # Dead fuel moisture of extinction
    fmce_dead = np.array(fc["fuelmce"], dtype=float)
    # Dead loads and SAVR
    w1   = np.array(fc["fuelload_gc01"], dtype=float)
    w10  = np.array(fc["fuelload_gc02"], dtype=float)
    w100 = np.array(fc["fuelload_gc03"], dtype=float)
    savr1   = np.array(fc["savr_gc01"], dtype=float)
    savr10  = np.array(fc["savr_gc02"], dtype=float)
    savr100 = np.array(fc["savr_gc03"], dtype=float)
    # Live loads and SAVR
    w_herb  = np.array(fc["fuelload_gc05"], dtype=float)
    w_woody = np.array(fc["fuelload_gc06"], dtype=float)
    savr_herb  = np.array(fc["savr_gc05"], dtype=float)
    savr_woody = np.array(fc["savr_gc06"], dtype=float)
    # Effective heating numbers
    hd_1 = np.zeros_like(w1, dtype=float)
    hd_10 = np.zeros_like(w1, dtype=float)
    hd_100 = np.zeros_like(w1, dtype=float)
    hl_herb = np.zeros_like(w1, dtype=float)
    hl_woody = np.zeros_like(w1, dtype=float)
    mask_s1   = savr1   != 0
    mask_s10  = savr10  != 0
    mask_s100 = savr100 != 0
    mask_sh   = savr_herb  != 0
    mask_sw   = savr_woody != 0
    hd_1[mask_s1]    = np.exp(-138.0 / savr1[mask_s1])
    hd_10[mask_s10]  = np.exp(-138.0 / savr10[mask_s10])
    hd_100[mask_s100]= np.exp(-138.0 / savr100[mask_s100])
    hl_herb[mask_sh] = np.exp(-500.0 / savr_herb[mask_sh])
    hl_woody[mask_sw]= np.exp(-500.0 / savr_woody[mask_sw])
    # Dead and live effective loads
    sum_dead_eff = w1*hd_1 + w10*hd_10 + w100*hd_100
    sum_live_eff = w_herb*hl_herb + w_woody*hl_woody
    # Dead-to-live effective loading ratio W
    W = np.zeros_like(w1, dtype=float)
    mask_live = sum_live_eff > 0
    W[mask_live] = sum_dead_eff[mask_live] / sum_live_eff[mask_live]
    # Dictionary output
    ratios = {
        "w1": w1, "w10": w10, "w100": w100,
        "hd_1": hd_1, "hd_10": hd_10, "hd_100": hd_100,
        "sum_dead_eff": sum_dead_eff, "sum_live_eff": sum_live_eff,
        "fmce_dead": fmce_dead, "W": W
    }
    return ratios
    
def calculate_fmce_live(
    fmc_1, fmc_10, fmc_100,
    ratios, n_fuel_cat
):
    """
    Compute live fuel moisture of extinction (fmce_live, in %) using
    the Albini/NFDRS formulation.

    Parameters
    ----------
    nml : dict
        Namelist with fuel_categories.
    n_fuel_cat : int
        1-based index of the fuel category.
    fmc_1, fmc_10, fmc_100 : float
        1-h, 10-h, 100-h dead FMC [%].
    fmce_dead : float
        Dead fuel moisture of extinction [%] for this fuel model.

    Returns
    -------
    fmce_live : float
        Live fuel moisture of extinction [%].
    """
    fmc_1   = np.asarray(fmc_1,   dtype=float)
    fmc_10  = np.asarray(fmc_10,  dtype=float)
    fmc_100 = np.asarray(fmc_100, dtype=float)
    # Index into fuel-model arrays (time, y, x) via broadcasting
    n_fuel_cat = np.asarray(n_fuel_cat, dtype=int)
    idx = n_fuel_cat[None, ...] - 1
    # Get information from ratios
    w1 = ratios["w1"][idx]
    w10 = ratios["w10"][idx]
    w100 = ratios["w100"][idx]
    hd_1 = ratios["hd_1"][idx]
    hd_10 = ratios["hd_10"][idx]
    hd_100 = ratios["hd_100"][idx]
    sum_dead_eff = ratios["sum_dead_eff"][idx]
    fmce_dead = ratios["fmce_dead"][idx]
    W = ratios["W"][idx]
    # Weighted fine dead FMC (fraction)
    numerator = (
        fmc_1   * w1   * hd_1 +
        fmc_10  * w10  * hd_10 +
        fmc_100 * w100 * hd_100
    )
    mask_dead = sum_dead_eff > 0
    mf_av = np.where(
        mask_dead,
        numerator / sum_dead_eff,
        0.0
    )
    # Albini expression (percent)
    mxd = fmce_dead
    fmce_live = 2.9 * W * (1.0 - mf_av / mxd) - 0.226
    # where there is no effective dead load, just use dead Mx
    fmce_live = np.where(mask_dead, fmce_live, fmce_dead)
    # Enforce constraints: cannot be < dead Mx, and clip to [0, 3] say.
    fmce_live = np.maximum(fmce_live, fmce_dead)
    fmce_live = np.clip(fmce_live, 0.0, 3.0)
    return fmce_live
    

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

def calculate_fmc(
    fmc_1, fmc_10, fmc_100, fmc_herba, fmc_woody,
    n_fuel_cat, weights, gsi
):
    n_fuel_cat = np.asarray(n_fuel_cat, dtype=int)
    idx = n_fuel_cat[None, ...] - 1  # same shape as grid
    # 1) Dead and live FMC (surface-area weighted within group)
    f_1     = weights["f_1"][idx]
    f_10    = weights["f_10"][idx]
    f_100   = weights["f_100"][idx]
    f_herba = weights["f_herba"][idx]
    f_woody = weights["f_woody"][idx]
    fmc_dead = f_1 * fmc_1 + f_10 * fmc_10 + f_100 * fmc_100
    fmc_live = f_herba * fmc_herba + f_woody * fmc_woody
    fmc_live[fmc_live == 0] = 0.001
    fmc = (0.1 * ((fmc_dead / fmc_live) - 1))
    fmc = fmc + gsi
    fmc = np.power(np.abs(fmc), 1.7)
    return fmc

def calculate_sawti(LFP, FMC):
    '''Calculates Santa Anna Wildfire Threat Index
    VARIBLES
    --------
    LFP: LargeFire Potential
    FMC: Fuel Moisture Component'''
    return LFP*FMC
 
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
    # Wind speed mph
    ws = np.sqrt(u10**2 + v10**2) * 2.23694
    # Find relative humidity [%]
    rh = calculate_rh(t2, q2, psfc) * 100
    # Calculate equilibrium moisture content (EMC)
    emc = calculate_emc(t2, rh)
    # Calculate moisture damping coefficient
    eta = calculate_eta(emc)
    # Calculate FFWI
    ffwi = calculate_ffwi(ws, eta)
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
    # Calculate HDW
    hdw = calculate_hdw(t2, q2, psfc, u10, v10)
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
    Td = calculate_dew_point(q2, psfc)
    # Conversion Celsius to Kelvin
    Td = Td + 273.15
    # Wind speed (m/s)
    ws = np.sqrt(u10**2 + v10**2) 
    # Calculate LFP
    lfp = calculate_lfp(t2, Td, ws)
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
