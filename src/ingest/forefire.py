import netCDF4 as nc
import numpy as np
from pyproj import Proj
import glob
import os, sys
import json
import pandas as pd
import subprocess
import geopandas as gpd
import fiona
from datetime import timedelta, datetime
import signal
import re
import shutil

#---------------------------------------------------------------------------
# configuration
#
# Every host path this module used to hardcode lives in etc/forefire.json now,
# with etc/forefire.json.initial tracked as the example -- the same convention
# etc/conf.json and etc/ngfs.json use. The path is relative, so the wrfxpy
# install the process runs from picks its own configuration.
#---------------------------------------------------------------------------

FF_CFG_PATH = 'etc/forefire.json'

def load_ff_cfg(cfg_path=FF_CFG_PATH):
    #reads the ForeFire configuration, falling back to the tracked example
    path = cfg_path
    if not os.path.exists(path):
        if os.path.exists(f"{path}.initial"):
            print(f"No {path}, using {path}.initial -- copy and edit it for this install")
            path = f"{path}.initial"
        else:
            raise IOError(f"No ForeFire configuration at {path} or {path}.initial")
    with open(path,'r') as openfile:
        return json.load(openfile)

def container_path(host_path,cfg):
    #where the apptainer bind makes a host path visible inside the container
    #forefire is handed bare filenames and resolves them against the working
    #directory, so a run directory outside the bind fails with files "missing"
    #inside a container that simply cannot see them -- catch that here instead
    host_bind = cfg['apptainer']['bind_host'].rstrip('/')
    cont_bind = cfg['apptainer']['bind_container'].rstrip('/')
    abs_path = os.path.abspath(host_path).rstrip('/')
    if abs_path != host_bind and not abs_path.startswith(f"{host_bind}/"):
        raise ValueError(
            f"{abs_path} is outside the apptainer bind {host_bind}, so ForeFire "
            f"would not see it. Put it under bind_host, or widen bind_host in {FF_CFG_PATH}.")
    return cont_bind + abs_path[len(host_bind):]

def bind_flag(cfg):
    #the --bind argument, host:container
    return f"{cfg['apptainer']['bind_host'].rstrip('/')}/:{cfg['apptainer']['bind_container'].rstrip('/')}/"

def ff_nc_name(step,grid_code,utc_str):
    #the FF netcdf filename for one time step, used when building and when linking
    return f"FF_{str(step).zfill(2)}_{grid_code}_{utc_str}.nc"

def read_fmda_fmc(geo_dir,lat,lon,radius_km=10,classes=(0,),weights=None):
    """Weighted mean fuel moisture from an FMDA geogrid tile, around one point.

    classes lists the dead size classes to combine (0=1h, 1=10h, 2=100h) and
    weights their relative contribution. Rothermel characterises dead fuel by a
    surface-area weighted moisture, so the default weights are the standard
    size-class SAVs -- which make it overwhelmingly a 1h quantity. Passing a
    single class reproduces the old behaviour exactly.

    The geogrid is written before wrf.exe runs, so this is the earliest moisture
    estimate available -- unlike FMC_GC_F in the wrfout, which is all zeros in
    wrfinput_d01 and spins up over the first couple of hours.

    Rows are stored south to north. The negative dy in index.json is a
    projection sign, not the array order: reading it as north-to-south makes the
    Sonoran desert wetter than the Olympic rainforest.

    Returns the mean as a fraction, or None if the tile cannot be read.
    """
    try:
        idx = json.load(open(f"{geo_dir}/index.json"))['FMC_GC']
        tiles = [t for t in glob.glob(f"{geo_dir}/FMC_GC/*") if os.path.basename(t) != 'index']
        if len(tiles) != 1:
            print(f"FMDA geogrid: expected one tile in {geo_dir}/FMC_GC, found {len(tiles)}")
            return None
        nx,ny,nz = idx['tile_x'],idx['tile_y'],idx['tile_z']
        classes = list(classes)
        weights = [1.0]*len(classes) if weights is None else list(weights)
        if len(weights) != len(classes):
            print(f"FMDA geogrid: {len(weights)} weights for {len(classes)} classes")
            return None
        bad = [c for c in classes if c >= nz or c < 0]
        if bad:
            print(f"FMDA geogrid: classes {bad} outside the {nz} available")
            return None
        if sum(weights) <= 0:
            print('FMDA geogrid: weights sum to zero')
            return None
        dt = ('>i%d' if idx.get('endian') == 'big' else '<i%d') % idx['wordsize']
        raw = np.fromfile(tiles[0],dtype=dt)
        if raw.size != nx*ny*nz:
            print(f"FMDA geogrid: {raw.size} values, expected {nx*ny*nz}")
            return None
        cube = raw.reshape(nz,ny,nx)

        proj = Proj(proj='lcc',lat_1=idx['truelat1'],lat_2=idx['truelat2'],lat_0=idx['truelat1'],
                    lon_0=idx['stdlon'],a=idx['radius'],b=idx['radius'])
        x0,y0 = proj(idx['known_lon'],idx['known_lat'])
        x,y = proj(lon,lat)
        i = int(round(idx['known_x'] - 1 + (x-x0)/idx['dx']))
        j = int(round(idx['known_y'] - 1 + (y-y0)/abs(idx['dy'])))
        if not (0 <= i < nx and 0 <= j < ny):
            print(f"FMDA geogrid: {lat},{lon} falls outside the tile")
            return None
        r = max(1,int(radius_km*1000/idx['dx']))
        total = 0.0
        for c,w in zip(classes,weights):
            band = cube[c].astype(np.float64) * idx['scale_factor']
            win = band[max(0,j-r):j+r+1, max(0,i-r):i+r+1]
            if not win.size:
                return None
            total += w * float(win.mean())
        return total / sum(weights)
    except Exception as e:
        print(f"FMDA geogrid: could not read {geo_dir}: {e}")
        return None


def resolve_md(wksp_dir,ign_latlon,cfg):
    """Dead fuel moisture for a run, with where it came from.

    Returns (md, provenance). md is None when the fuel table should be left
    alone, which is the behaviour whenever the feature is off or the source
    cannot be trusted -- so a broken FMDA instance degrades to today's run
    rather than to a silently wrong fire.
    """
    mc = cfg.get('moisture',{})
    if not mc.get('enabled'):
        return None,'moisture disabled, fuel table unchanged'
    lo,hi = mc.get('valid_range',[0.02,0.50])
    fallback = mc.get('fallback_md')
    #which dead size classes to combine, and how. 'single' keeps the old
    #class_index behaviour; 'sav' weights by size-class surface-area-to-volume,
    #which is what Rothermel's single dead moisture is meant to represent
    if mc.get('weighting','single') == 'sav':
        dead_classes = tuple(mc.get('dead_classes',[0,1,2]))
        sav_weights = tuple(mc.get('sav_weights',[2000.0,109.0,30.0]))
    else:
        dead_classes = (mc.get('class_index',0),)
        sav_weights = (1.0,)

    md = None
    src = mc.get('source','fmda_geogrid')
    if src == 'constant':
        md = fallback
        if md is not None:
            return md,f"constant Md={md:.4f} from config"
    elif src == 'fmda_geogrid':
        try:
            geo = json.load(open(f"{wksp_dir}/input.json")).get('fmda_geogrid_path')
        except Exception as e:
            geo = None
            print(f"could not read input.json for fmda_geogrid_path: {e}")
        if not geo:
            print('no fmda_geogrid_path in input.json')
        elif not os.path.isdir(geo):
            print(f"fmda_geogrid_path does not exist: {geo}")
        else:
            md = read_fmda_fmc(geo,ign_latlon[0],ign_latlon[1],mc.get('radius_km',10),
                               dead_classes,sav_weights)
    else:
        print(f"unknown moisture source {src!r}")

    if md is None or not (lo <= md <= hi):
        if md is not None:
            print(f"FMDA moisture {md:.4f} outside valid_range [{lo}, {hi}], rejecting")
        if fallback is None:
            return None,'FMDA moisture unusable, fuel table unchanged'
        return fallback,f"FMDA moisture unusable, falling back to Md={fallback:.4f}"
    how = (f"SAV-weighted classes {list(dead_classes)}" if len(dead_classes) > 1
           else f"class {dead_classes[0]}")
    return md,f"Md={md:.4f} from {src} ({how}, {mc.get('radius_km',10)} km)"


def write_fuel_table(md,cfg):
    """Copies the base fuel table with Md replaced in every fuel class.

    Plain text on purpose: the table is ';' separated with '\n' endings, and
    rewriting it through csv.writer appends '\r' to the last column, which makes
    ForeFire fail to find 'me' and silently stop the fire from spreading.
    """
    base = cfg['fuels_table']
    out_dir = f"{cfg['run_dir']}/fueltables"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out = f"{out_dir}/fuels_Md_{md:.4f}.csv"
    lines = open(base,newline='').read().split('\n')
    col = lines[0].split(';').index('Md')
    keep = [lines[0]]
    for ln in lines[1:]:
        if not ln.strip():
            keep.append(ln); continue
        f = ln.split(';'); f[col] = f"{md:.4f}"; keep.append(';'.join(f))
    open(out,'w',newline='').write('\n'.join(keep))
    return out


def moisture_params(wksp_dir,ign_latlon,cfg):
    """{'fuelsTableFile': <container path>} for a run, or {} to leave it alone."""
    md,why = resolve_md(wksp_dir,ign_latlon,cfg)
    print(f"ForeFire fuel moisture: {why}")
    if md is None:
        return {}
    return {'fuelsTableFile': container_path(write_fuel_table(md,cfg),cfg)}


def apply_ff_params(text,params):
    #overrides setParameter[key=value] lines in a .ff script; a key the template
    #does not already set is added next to the other tuning parameters
    for key,value in params.items():
        new_line = f"setParameter[{key}={value}]"
        pattern = re.compile(r'^([ \t]*)setParameter\[' + re.escape(str(key)) + r'=[^\]]*\]',re.M)
        text,n = pattern.subn(lambda m: m.group(1) + new_line, text)
        if n == 0:
            anchor = '# --- Load Data & Domain ---'
            if anchor in text:
                text = text.replace(anchor,f"{new_line}\n\n{anchor}",1)
            else:
                text = f"{new_line}\n{text}"
    return text

def param_tag(params):
    #directory-safe name for a parameter set, e.g. windReductionFactor_0.6
    return '_'.join(f"{k}_{v}" for k,v in sorted(params.items()))

def read_input(wksp_dir):
    #gets the ignition point and time from the input.json file
    input = glob.glob(f'{wksp_dir}/input.json')
    if not input:
        print('No input file found...')
        return None
    with open(input[0],'r') as openfile:
            cfg = json.load(openfile) 
            start_utc = cfg['start_utc']
            ign_utc = cfg['ignitions']['1'][0]['time_utc']          #ignitions is a list   
            ign_latlon = cfg['ignitions']['1'][0]['latlon'] 
            grid_code = cfg['grid_code']
            print(ign_utc,ign_latlon)
    return start_utc,ign_utc,ign_latlon,grid_code

def make_timing_table(wksp_dir,ign_utc):
    #ign utc is the string from the json input file for wrfxpy job
    g = glob.glob(f'{wksp_dir}/wrf/wrfout*')
    if len(g) == 0:
        g = glob.glob(f'{wksp_dir}/wrfout*')
    if len(g) < 5:
        g = glob.glob(f'{wksp_dir}/wrf/saveout*')
    if len(g) < 2: #return empty DF
        return pd.DataFrame()
    #t_step and every seconds column below are measured off g[0], so order matters
    g = sort_by_time(g)
    #compute the time step between wrfouts
    utc0 = make_UTC_string(g[0][-19:])
    utc1 = make_UTC_string(g[1][-19:])
    t_step = int((pd.Timestamp(utc1) - pd.Timestamp(utc0)).total_seconds())

    
    #compute the seconds aft first wrfout for the ignition
    ig_time = pd.Timestamp(make_UTC_string(ign_utc))
    ignition_seconds = int((ig_time-pd.Timestamp(utc0)).total_seconds())
    
    wrf_time = []
    wrf_seconds = []
    ign_seconds = []
    UTC_str = []
    
    for i,gg in enumerate(g):
        UTC_str.append(make_UTC_string(gg[-19:]))
        wt = pd.Timestamp(UTC_str[-1])
        wrf_time.append(wt)
        ws = (wt-wrf_time[0]).total_seconds() #seconds of the start of the wrfout file
        wrf_seconds.append(int(ws))

        ig_s = (wt-ig_time).total_seconds()
        if ws + t_step < ignition_seconds:
            ig_s = -9999 # negative numbers indicate spinup period
        elif (ws < ignition_seconds) and (ws + t_step > ignition_seconds):
            ig_s = ignition_seconds - ws
        else:
            ig_s = ws

        ign_seconds.append(int(ig_s))
    df = pd.DataFrame({
         "wrfout" : g,
         'UTC_str': UTC_str,
         "wrf_seconds" : wrf_seconds,
         "ign_seconds" : ign_seconds
    }
    )
    return df
          
     
     

def find_wrfouts(wksp_dir,ign_utc):
    #strips out the spinup wrfouts
    g = glob.glob(f'{wksp_dir}/wrf/wrfout*')
    if len(g) < 5:
        g = glob.glob(f'{wksp_dir}/wrf/saveout*')
    g = sort_by_time(g)

    strip = []
    keep = []
    ## gg = '/data/jhaley/wrfxpy/wksp/wfc-Wichita_Hills_2026-08-14_17_00_00_9CE196E0-4EE8-45B4-B5BC-279677D1FC29-2026-08-14_15:00:00-30/wrf/wrfout_d01_2026-08-15_21:00:00'
    ## gg[-19:] = '2026-08-15_21:00:00'
    for i,gg in enumerate(g):
         if gg[-19:] < ign_utc:
              strip.append(gg)
         else:
              keep.append(gg)
    
    #the last file stripped is the first file before ingition ime
    keep.insert(0,strip[-1])
         
    return keep

def sort_by_time(g):
     #glob returns directory order, not time order, so sort on the
     #YYYY-MM-DD_HH:MM:SS stamp that ends every wrfout/saveout name
     return sorted(g, key=lambda gg: gg[-19:])

def make_UTC_string(ign_utc):
     #makes a forefire time string like YYYY-MM-DDTHH:MM:SSZ from something like 2026-08-12_03:26:48
                                       #2026-08-12_03:26:48
     if ign_utc[:4] == '0001':
         ign_utc = ign_utc.replace('0001','2001')

     ign_utc = ign_utc.replace('_','T')
     return f'{ign_utc}Z'

def make_script_set(forefire_dir,timing_table,ign_latlon,grid_code,cfg=None,params=None):
    #params overrides ForeFire setParameter values, e.g. {'windReductionFactor':0.6}
    #make individual ignition or restart scripts and save them, addtheir names to the timming_table to restart o
    #forefire_dir = '/home/jhaley/forefire/tests/rosenbaum'

    #write geojson files like
    #sim_GRID_CODE_final_STEP_NUMBER.geojson sim_Wichita_Hills_2026-08-14_17_00_00_9CE196E0-4EE8-45B4-B5BC-279677D1FC29_final_12.geojson

    #sim_GRID_CODE_STEP_NUMBER_END_STEP.geojson sim_Wichita_Hills_2026-08-14_17_00_00_9CE196E0-4EE8-45B4-B5BC-279677D1FC29_12_10059.geojson
    #in the workspace forefire directory

    if cfg is None:
        cfg = load_ff_cfg()
    restart_template = cfg['templates']['restart']
    ignition_template = cfg['templates']['ignition']
    #time of the last wrfout, ign_seconds
    end_time = timing_table['ign_seconds'].max()
    #number of time steps
    steps = len(timing_table)
    ignition_step = True

    n = len(timing_table)
    nc_file = [None]*n  #WRF data reformatted for FOREFIRE
    step_result = [None]*n #FireFire native perimeter
    script_file = [None]*n
    step_geo = [None]*n
    final_geo = [None]*n

    for i,row in timing_table.iterrows():
        if i < steps - 1: #don't go past last steps
            if row['ign_seconds'] > 0:
                row_string = str(i).zfill(2)   #this is the step number
                out_file = f"{forefire_dir}/{ff_nc_name(i,grid_code,row['UTC_str'])}"
                ignition_file = f"{forefire_dir}/ignition_{row_string}_{grid_code}.ff" 
                restart_file = f"{forefire_dir}/restart_{row_string}_{grid_code}.ff"                      
                
                if not os.path.exists(out_file):
                    print('Making FF netcdf file',out_file)
                    print(row)
                    print()
                    try:
                        make_FF_nc(row['wrfout'],out_file)
                    except:
                        ff_ideal_nc(row['wrfout'],out_file)
                else:
                    print('FF netcf exists already',out_file)
                
                #make script files, IGN_LAT and IGN_LONG not in restart scripts
                start_step_string =  str(row['ign_seconds']) #time in seconds after ignition
                end_step_string = str(timing_table.loc[i+1]['ign_seconds']) 
                
                change_dict = {
                    "FF_NC" : os.path.basename(out_file),
                    "UTC_STRING" : str(row['UTC_str']),
                    "START_STEP" : start_step_string,
                    "END_STEP" : end_step_string,
                    "END_TIME" : str(end_time),
                    "STEP_NUMBER" : row_string,
                    "IGN_LONG" : str(ign_latlon[1]),
                    "IGN_LAT" : str(ign_latlon[0]),
                    "GRID_CODE" : grid_code
                }
                if ignition_step:
                    #os.system('cp {ignition_template} {ignition_file}')
                    #copy initial script file and change values
                    ignition_step = False # will do restart steps after this
                    template = ignition_template
                    FF_script = ignition_file
                else:
                    template = restart_template
                    FF_script = restart_file

                #open template file
                with open(template, "r") as f:
                    file_content = f.read()

                #replace the placehoilder values
                for placeholder, value in change_dict.items():
                    print(placeholder,value)
                    file_content = file_content.replace(placeholder, value)

                #override tuning parameters for a parameter-variation run
                if params:
                    file_content = apply_ff_params(file_content,params)
                
                #save new file
                if True: #not os.path.exists(restart_file):
                    with open(FF_script,"w") as f:
                        f.write(file_content)

                #update and return timing_table with inpuit and output files for forfire
                ## inputs
                nc_file[i]= os.path.basename(out_file)
                
                script_file[i] = os.path.basename(FF_script)
                ## output(s)
                #res_file = 'simulation_result_END_STEP.ff'.replace('END_STEP',end_step_string)
                res_file = f"sim_{grid_code}_{end_step_string}.ff"
                step_result[i] = res_file

                #add expected geojson files
                final_geo[i] = f"sim_{grid_code}_final_{row_string}.geojson"
                step_geo[i] = f"sim_{grid_code}_{row_string}_{end_step_string}.geojson"




    #update timining_table
    timing_table['nc_file'] = nc_file
    timing_table['script_file'] = script_file
    timing_table['result_file'] = step_result
    timing_table['step_geojson'] = step_geo
    timing_table['final_geojson'] = final_geo

    return timing_table
    

def ff_ideal_nc(nc_path,out_path):
    wrf_in = nc.Dataset(nc_path, "r") #get winds from these
    wrfinput = f"{os.path.dirname(nc_path)}/wrfinput_d01" #het fuel and top from wrfinput_d01 file
    wrfinput_d01 = nc.Dataset(wrfinput,"r")

    dx = wrf_in.DX
    dy = wrf_in.DY
    nx_coarse = wrf_in.dimensions["west_east"].size
    ny_coarse = wrf_in.dimensions["south_north"].size
    sub_grid_ratio = wrf_in.dimensions["west_east_subgrid"].size / wrf_in.dimensions["west_east_stag"].size
    srx = int(wrf_in.dimensions['west_east_subgrid'].size/(wrf_in.dimensions['west_east'].size+1))      #  <<<<------ removes the strip around the edg
    sry = int(wrf_in.dimensions['south_north_subgrid'].size/(wrf_in.dimensions['south_north'].size+1))  #  <<<<------ removes the strip around the edg

    
    x_subgrid_size = wrf_in.dimensions["west_east_subgrid"].size  - srx
    y_subgrid_size = wrf_in.dimensions["south_north_subgrid"].size  - sry
    #subgrid_size = wrf_in.dimensions["west_east"].size*sub_grid_ratio

    # Compute total domain width in meters based on the atmospheric scale
    total_x_meters = float(nx_coarse * dx)  # 30 cells * 1000m = 30000.0
    total_y_meters = float(ny_coarse * dy)


    bbox_wsen = f"0,0,{total_x_meters},{total_y_meters}"

    fuel_matrix = np.array(wrfinput_d01.variables['NFUEL_CAT'][0,:-sry,:-srx])    #  <<<<------ removes the strip around the edg
    #topo_matrix = wrf_in.variables["ZSF"][0 , :, :]
    topo_matrix = np.array(wrfinput_d01.variables['ZSF'][0,:-sry,:-srx])    #  <<<<------ removes the strip around the edg
    #write 82 into the fuel array tempoarily
    #fuel_matrix[:,:] = 82
    u_matrix = np.array(wrf_in.variables['UF'][0,:-sry,:-srx]) 
    v_matrix = np.array(wrf_in.variables['VF'][0,:-sry,:-srx]) 

    # --- WRITE THE FOREFIRE COMPATIBLE DATASET ---
    #put in block to extract the p-ath of the original nc file
    with nc.Dataset(out_path, "w", format="NETCDF4") as ff_out:
        # Match grid resolutions perfectly by using subgrid sizes
        ff_out.createDimension("nx", x_subgrid_size)
        ff_out.createDimension("ny", y_subgrid_size)
        ff_out.createDimension("nz", 1)
        ff_out.createDimension("nt", 1)

        ff_out.createDimension("fx", x_subgrid_size)
        ff_out.createDimension("fy", y_subgrid_size)
        ff_out.createDimension("fz", 1)
        ff_out.createDimension("ft", 1)
        
        wind_x = x_subgrid_size
        wind_y = y_subgrid_size
        ff_out.createDimension('wind_rows',wind_y)
        ff_out.createDimension('wind_columns',wind_x)
        ff_out.createDimension('wind_dimensions',1)
        ff_out.createDimension('wind_directions',1)
        
        # Inject domain positioning strings
        domain = ff_out.createVariable("domain", str)
        domain.type = "domain"
        domain.BBoxWSEN = bbox_wsen
        domain.SWx = np.float32(0.0)
        domain.SWy = np.float32(0.0)
        domain.Lx =  np.float32(total_x_meters)   # Sets Lx = 30000.0 np.float32(12000.)
        domain.Ly =  np.float32(total_y_meters)   # Sets Ly = 30000.0 np.float32(12000.)
        domain.Lz = np.float32(0.0)
        domain.t0 = np.float32(0.0)
        domain.Lt = np.float32(np.inf)    # Matches the 'Inf' tracking scale
        domain.SWz = np.float32(0.0)


        # Write 4D Fuel Layer (Cast to 16-bit Integer)
        fuel = ff_out.createVariable(
            "fuel", "i2", ("ft", "fz", "fy", "fx"), fill_value=False
        )
        fuel.type = "fuel"
        #fuel[:, :, 0, 0] = fuel_matrix.astype(np.int16)
        four_d_fuel = fuel_matrix[np.newaxis, np.newaxis,:,:].astype(np.int16)
        # Ensure memory is ordered properly for the C++ library
        four_d_fuel = np.ascontiguousarray(four_d_fuel)
        fuel[:] = four_d_fuel


        # Write 4D Elevation Layer (Cast to 16-bit Integer)
        altitude = ff_out.createVariable(
            "altitude", "i2", ("nt", "nz", "ny", "nx"), fill_value=False
        )
        altitude.type = "data"
        #altitude[:, :, 0, 0] = topo_matrix.astype(np.int16)
        four_d_alt = topo_matrix[np.newaxis, np.newaxis,:,:].astype(np.int16)
        four_d_alt = np.ascontiguousarray(four_d_alt)
        altitude[:] = four_d_alt

        
        
        wind_u = u_matrix[np.newaxis,np.newaxis,:,:]
        wind_v = v_matrix[np.newaxis,np.newaxis,:,:]

        w = np.concatenate((wind_u,wind_v),axis=1)
        w = np.ascontiguousarray(w)
        windU = ff_out.createVariable(
            "windU","f4",("wind_dimensions","wind_directions","wind_rows","wind_columns"),fill_value=np.nan
        )
        windU.type = "wind"
        windU[:] = w[:,0,:,:]

        windV = ff_out.createVariable(
            "windV","f4",("wind_dimensions","wind_directions","wind_rows","wind_columns"),fill_value=np.nan
        )
        windV.type = "wind"
        windV[:] = w[:,1,:,:]


def make_FF_nc(nc_path,out_path):
    #takes in a wrf file (wrfout) and writes something compatible with ForeFire
        # Open the WRF Input file containing pre-sliced Landfire data
    #f = '/data/jhaley/wrfxpy/wksp/wfc-TITAN_2023-06-28_22:00:00_4A7E0BD9-1D57-4263-B49D-547B19BB2A38-2023-06-28_21:00:00-27/wrf/wrfout_d01_2023-06-30_00:00:00'
    #should get the following variables NFUEL_CAT, ZSF, UF, VF
    wrf_in = nc.Dataset(nc_path, "r") #get winds from these
    wrfinput = f"{os.path.dirname(nc_path)}/wrfinput_d01" #het fuel and top from wrfinput_d01 file
    wrfinput_d01 = nc.Dataset(wrfinput,"r")

    # 1. Dynamically read global projection metadata from WRF
    lat_0 = wrf_in.CEN_LAT
    lon_0 = wrf_in.CEN_LON
    lat_1 = wrf_in.TRUELAT1
    lat_2 = wrf_in.TRUELAT2

    # 2. Extract physical metrics and dimensions
    dx = wrf_in.DX
    dy = wrf_in.DY
    nx_coarse = wrf_in.dimensions["west_east"].size
    ny_coarse = wrf_in.dimensions["south_north"].size
    sub_grid_ratio = wrf_in.dimensions["west_east_subgrid"].size / wrf_in.dimensions["west_east_stag"].size
    srx = int(wrf_in.dimensions['west_east_subgrid'].size/(wrf_in.dimensions['west_east'].size+1))      #  <<<<------ removes the strip around the edg
    sry = int(wrf_in.dimensions['south_north_subgrid'].size/(wrf_in.dimensions['south_north'].size+1))  #  <<<<------ removes the strip around the edg

    
    x_subgrid_size = wrf_in.dimensions["west_east_subgrid"].size  - srx
    y_subgrid_size = wrf_in.dimensions["south_north_subgrid"].size  - sry
    #subgrid_size = wrf_in.dimensions["west_east"].size*sub_grid_ratio

    # Compute total domain width in meters based on the atmospheric scale
    total_x_meters = float(nx_coarse * dx)  # 30 cells * 1000m = 30000.0
    total_y_meters = float(ny_coarse * dy)

    # 3. Use PyProj to map the geographic bounding box (WSEN)
    proj_lcc = Proj(
        proj="lcc",
        lat_1=lat_1,
        lat_2=lat_2,
        lat_0=lat_0,
        lon_0=lon_0,
        a=6370000,
        b=6370000,
    )
    half_x = total_x_meters / 2
    half_y = total_y_meters / 2

    west, south = proj_lcc(-half_x, -half_y, inverse=True)
    east, north = proj_lcc(half_x, half_y, inverse=True)

    # 4. Generate the exact strings ForeFire uses to anchor spatial logic
    bbox_wsen = f"{west},{south},{east},{north}"
    wsenlbrt_str = f"{west},{south},{east},{north},0.0,0.0,{total_x_meters},{total_y_meters}"

    # 5. Extract raw data blocks (dropping the trailing Time dimension)
    # Shape goes from (1240, 1240, 1) -> (1240, 1240)
    #fuel_matrix = wrf_in.variables["NFUEL_CAT"][0, :, :]

    fuel_matrix = np.array(wrfinput_d01.variables['NFUEL_CAT'][0,:-sry,:-srx])    #  <<<<------ removes the strip around the edg
    #topo_matrix = wrf_in.variables["ZSF"][0 , :, :]
    topo_matrix = np.array(wrfinput_d01.variables['ZSF'][0,:-sry,:-srx])    #  <<<<------ removes the strip around the edg
    #write 82 into the fuel array tempoarily
    #fuel_matrix[:,:] = 82
    u_matrix = np.array(wrf_in.variables['UF'][0,:-sry,:-srx]) 
    v_matrix = np.array(wrf_in.variables['VF'][0,:-sry,:-srx]) 

    # --- WRITE THE FOREFIRE COMPATIBLE DATASET ---
    #put in block to extract the p-ath of the original nc file
    with nc.Dataset(out_path, "w", format="NETCDF4") as ff_out:
        # Match grid resolutions perfectly by using subgrid sizes
        ff_out.createDimension("nx", x_subgrid_size)
        ff_out.createDimension("ny", y_subgrid_size)
        ff_out.createDimension("nz", 1)
        ff_out.createDimension("nt", 1)

        ff_out.createDimension("fx", x_subgrid_size)
        ff_out.createDimension("fy", y_subgrid_size)
        ff_out.createDimension("fz", 1)
        ff_out.createDimension("ft", 1)
        
        wind_x = x_subgrid_size
        wind_y = y_subgrid_size
        ff_out.createDimension('wind_rows',wind_y)
        ff_out.createDimension('wind_columns',wind_x)
        ff_out.createDimension('wind_dimensions',1)
        ff_out.createDimension('wind_directions',1)
        
        # Inject domain positioning strings
        domain = ff_out.createVariable("domain", str)
        domain.type = "domain"
        domain.BBoxWSEN = bbox_wsen
        domain.WSENLBRT = wsenlbrt_str
        domain.SWx = np.float32(0.0)
        domain.SWy = np.float32(0.0)
        domain.Lx =  np.float32(total_x_meters)   # Sets Lx = 30000.0 np.float32(12000.)
        domain.Ly =  np.float32(total_y_meters)   # Sets Ly = 30000.0 np.float32(12000.)
        domain.Lz = np.float32(0.0)
        domain.t0 = np.float32(0.0)
        domain.Lt = np.float32(np.inf)    # Matches the 'Inf' tracking scale
        domain.SWz = np.float32(0.0)


        # Write 4D Fuel Layer (Cast to 16-bit Integer)
        fuel = ff_out.createVariable(
            "fuel", "i2", ("ft", "fz", "fy", "fx"), fill_value=False
        )
        fuel.type = "fuel"
        #fuel[:, :, 0, 0] = fuel_matrix.astype(np.int16)
        four_d_fuel = fuel_matrix[np.newaxis, np.newaxis,:,:].astype(np.int16)
        # Ensure memory is ordered properly for the C++ library
        four_d_fuel = np.ascontiguousarray(four_d_fuel)
        fuel[:] = four_d_fuel


        # Write 4D Elevation Layer (Cast to 16-bit Integer)
        altitude = ff_out.createVariable(
            "altitude", "i2", ("nt", "nz", "ny", "nx"), fill_value=False
        )
        altitude.type = "data"
        #altitude[:, :, 0, 0] = topo_matrix.astype(np.int16)
        four_d_alt = topo_matrix[np.newaxis, np.newaxis,:,:].astype(np.int16)
        four_d_alt = np.ascontiguousarray(four_d_alt)
        altitude[:] = four_d_alt

        
        
        wind_u = u_matrix[np.newaxis,np.newaxis,:,:]
        wind_v = v_matrix[np.newaxis,np.newaxis,:,:]

        w = np.concatenate((wind_u,wind_v),axis=1)
        w = np.ascontiguousarray(w)
        windU = ff_out.createVariable(
            "windU","f4",("wind_dimensions","wind_directions","wind_rows","wind_columns"),fill_value=np.nan
        )
        windU.type = "wind"
        windU[:] = w[:,0,:,:]

        windV = ff_out.createVariable(
            "windV","f4",("wind_dimensions","wind_directions","wind_rows","wind_columns"),fill_value=np.nan
        )
        windV.type = "wind"
        windV[:] = w[:,1,:,:]

def stage_ff_nc(wksp_dir,timing_table,grid_code,cfg):
    """Puts one copy of a fire's FF netcdf files under the apptainer bind.

    The netcdf files carry the WRF weather, fuels and terrain, none of which
    depend on the ForeFire tuning parameters, so a sweep builds them once here
    and every parameter directory links to them. Files an earlier run already
    moved into the workspace are copied in rather than rebuilt -- copied and not
    linked, because the workspace filesystem is not bound into the container and
    a link there would dangle. Returns the staging directory.
    """
    stage_dir = f"{cfg['run_dir']}/nc_{grid_code}"
    container_path(stage_dir,cfg)
    if not os.path.exists(stage_dir):
        os.makedirs(stage_dir)

    steps = len(timing_table)
    for i,row in timing_table.iterrows():
        if i >= steps - 1 or row['ign_seconds'] <= 0:
            continue
        name = ff_nc_name(i,grid_code,row['UTC_str'])
        target = f"{stage_dir}/{name}"
        if os.path.exists(target):
            continue
        from_wksp = f"{wksp_dir}/forefire/{name}"
        if os.path.exists(from_wksp):
            print(f"Staging {name} from the workspace")
            shutil.copy2(from_wksp,target)
        else:
            print(f"Building {name}")
            try:
                make_FF_nc(row['wrfout'],target)
            except Exception:
                ff_ideal_nc(row['wrfout'],target)
    return stage_dir


def sweep_ff_params(wksp_dir,param_sets,cfg=None,overwrite=False,reuse_nc=True,tags=None):
    """Runs one fire once per ForeFire parameter set, each in its own directory.

    param_sets is a list of dicts of setParameter names to values, e.g.
    [{'windReductionFactor':0.4},{'windReductionFactor':0.6}]. Each set gets a
    run directory named for its parameters under run_dir, and its results land in
    <wksp>/forefire/<tag>. Keeping them apart matters: the restart template
    chains time steps through include[sim_<grid>_<start>.ff], so two parameter
    sets sharing a directory would read each other's perimeters.

    tags optionally names the run directories instead of deriving them from the
    parameters, for a sweep whose varying quantity is not itself a ForeFire
    parameter -- fuel moisture, say, which is varied by pointing fuelsTableFile
    at a different table.

    Returns {tag: timing_table}.
    """
    if cfg is None:
        cfg = load_ff_cfg()
    if tags is not None and len(tags) != len(param_sets):
        raise ValueError(f"{len(tags)} tags for {len(param_sets)} parameter sets")
    start_utc,ign_utc,ign_latlon,grid_code = read_input(wksp_dir)
    timing_table = make_timing_table(wksp_dir,ign_utc)
    if len(timing_table) <= cfg.get('min_timesteps',20):
        print('Not running ForeFire forecasts, too few wrfout files')
        return {}

    stage_dir = stage_ff_nc(wksp_dir,timing_table,grid_code,cfg) if reuse_nc else None
    #resolved once: every parameter set in the sweep shares the same moisture
    base_params = moisture_params(wksp_dir,ign_latlon,cfg)

    results = {}
    for i,params in enumerate(param_sets):
        tag = tags[i] if tags is not None else param_tag(params)
        run_dir = f"{cfg['run_dir']}/{tag}"
        container_path(run_dir,cfg)
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)

        #link the shared netcdf files in, relatively, so the paths resolve both
        #on the host and under the container's view of the bind
        if stage_dir:
            for name in sorted(os.listdir(stage_dir)):
                link = f"{run_dir}/{name}"
                if not os.path.exists(link):
                    os.symlink(os.path.relpath(f"{stage_dir}/{name}",run_dir),link)

        print(f"\n=== ForeFire parameter set {tag} in {run_dir}")
        set_params = dict(base_params); set_params.update(params or {})
        tt = make_script_set(run_dir,timing_table.copy(),ign_latlon,grid_code,cfg=cfg,params=set_params)
        tt.to_csv(f"{run_dir}/timing_table.csv",index=False)
        run_timing_table(tt,overwrite=overwrite,cfg=cfg,run_dir=run_dir,
                         dest_dir=f"{wksp_dir}/forefire/{tag}")

        #drop the netcdf links before the results are moved; they point into the
        #staging directory and would dangle in the workspace
        for name in os.listdir(run_dir):
            link = f"{run_dir}/{name}"
            if os.path.islink(link):
                os.unlink(link)

        cleanup_ff_run(run_dir,wksp_dir,grid_code,cfg=cfg,
                       dest_dir=f"{wksp_dir}/forefire/{tag}",suffix=tag)
        results[tag] = tt
    return results


def vary_ff_params(wksp_dir,param,param_value,cfg=None,overwrite=False):
    #one parameter, one value: vary_ff_params(wksp_dir,'windReductionFactor',0.6)
    return sweep_ff_params(wksp_dir,[{param:param_value}],cfg=cfg,overwrite=overwrite)

        
def run_ff_set(ff_set):
    #runs all the forefire scripts in the list of dataframe that contains the indivdual ff files and nc files
    pass

def run_timing_table(df,overwrite=False,cfg=None,run_dir=None,dest_dir=None):
    #run_dir is where the scripts and netcdf files were written; it defaults to
    #the configured run_dir and must match what make_script_set was given
    #dest_dir is where finished results live, and decides which steps count as
    #already done; a parameter sweep must pass its own or it sees the baseline
    #run's perimeters and skips every step

    #look for the files like
    #sim_GRID_CODE_final_STEP_NUMBER.geojson sim_Wichita_Hills_2026-08-14_17_00_00_9CE196E0-4EE8-45B4-B5BC-279677D1FC29_final_12.geojson
    #sim_GRID_CODE_STEP_NUMBER_END_STEP.geojson sim_Wichita_Hills_2026-08-14_17_00_00_9CE196E0-4EE8-45B4-B5BC-279677D1FC29_12_10059.geojson
    #in the workspace forefire directory
    if dest_dir is None:
        dest_dir = f"{os.path.dirname(os.path.dirname(df.wrfout.iloc[0]))}/forefire"
    print(dest_dir)

    # 1. Configuration, from etc/forefire.json
    if cfg is None:
        cfg = load_ff_cfg()
    HOST_BIND_DIR = run_dir if run_dir else cfg['run_dir']
    CONTAINER_PWD = container_path(HOST_BIND_DIR,cfg)
    SIF_IMAGE = cfg['apptainer']['sif_image']
    BIND_FLAG = bind_flag(cfg)
    timeout_s = cfg['apptainer'].get('timeout_s',240)

    # Load your dataframe (replace with your actual data source)
    # df = pd.read_csv("your_data.csv")

    # 2. Filter dataframe for rows where required fields are not empty/NaN
    valid_runs = df.dropna(subset=["nc_file", "script_file", "result_file"])

    print(f"Found {len(valid_runs)} valid simulation runs to process.")

    # 3. Iterate through valid entries
    for row in valid_runs.itertuples():
        final_geo = row.final_geojson
        if not overwrite and os.path.exists(f"{dest_dir}/{final_geo}"):
            print('Time step completed, skipping to next row')
            continue
        script_name = row.script_file
        result_name = row.result_file
        
        # Define exact host paths to check file existence
        host_script_path = os.path.join(HOST_BIND_DIR, script_name)
        host_nc_path = os.path.join(HOST_BIND_DIR, row.nc_file)

        # Create a unique log file path in the same host folder
        log_file_name = f"sim_{os.path.splitext(result_name)[0]}.log"
        host_log_path = os.path.join(HOST_BIND_DIR, log_file_name)
        
        # Verify files exist on host before launching the container
        if not (os.path.exists(host_script_path) and os.path.exists(host_nc_path)):
            print(f"Skipping row {row.Index}: Missing host files for {script_name}")
            continue
            
        # Construct the exact Apptainer command array
        cmd = [
            "apptainer", "exec",
            "--pwd", CONTAINER_PWD,
            "--writable-tmpfs",
            "--bind", BIND_FLAG,
            SIF_IMAGE,
            "forefire", "-i", script_name
        ]
        
        print(f"\nExecuting: forefire -i {script_name} for timestamp {row.UTC_str}")

        # Open the log file and stream the container output into it
        with open(host_log_path, "w") as log_file:
            try:
                # Start Apptainer in a new process group session
                proc = subprocess.Popen(
                    cmd, 
                    stdout=log_file, 
                    stderr=subprocess.STDOUT, 
                    text=True,
                    start_new_session=True
                )
                
                proc.wait(timeout=timeout_s)
                
                # Check if the process exited with an error
                if proc.returncode != 0:
                    raise subprocess.CalledProcessError(proc.returncode, cmd)
                    
                print(f"Success: Process completed for {script_name}")
                
            except subprocess.TimeoutExpired:
                print(f"Error: Process timed out after {timeout_s} seconds for {script_name}")
                
                # Kill the entire process group (Apptainer + inner forefire process)
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except ProcessLookupError:
                    pass # Process already exited right at the deadline
                    
                print(f"Check logs for partial output: {host_log_path}")
                print("Halting automation loop due to timeout.")
                break
                
            except subprocess.CalledProcessError as e:
                print(f"Error: Command failed on {script_name} with exit code {e.returncode}")
                print(f"Check logs for details: {host_log_path}")
                print("Halting automation loop for troubleshooting.")
                break

        '''
        with open(host_log_path, "w") as log_file:
            try:
                # stdout=log_file and stderr=subprocess.STDOUT merges all output into the log file
                subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, check=True, text=True)
                print(f"Success: Process completed for {script_name}")
                
            except subprocess.CalledProcessError as e:
                print(f"Error: Command failed on {script_name} with exit code {e.returncode}")
                print(f"Check logs for details: {host_log_path}")
                
                # Halt loop on failure to protect data integrity
                print("Halting automation loop for troubleshooting.")
                break
        '''


def cleanup_ff_run(forefire_dir,wksp_dir,grid_code,cfg=None,dest_dir=None,suffix=None):
    #destination for files; a parameter-variation run passes its own so the
    #results of one parameter set do not overwrite another's
    #suffix tags every geojson and kml name, e.g. windReductionFactor_0.6, so
    #several parameter sets can be opened together in Google Earth
    if cfg is None:
        cfg = load_ff_cfg()
    if dest_dir is None:
        dest_dir = f"{wksp_dir}/forefire"
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    #move all the files in the forefire_dir to the destination_dir
    
    os.system(f"mv {forefire_dir}/* {dest_dir}/.")

    finalize_ff_outputs(dest_dir,grid_code,cfg=cfg,suffix=suffix)


def finalize_ff_outputs(dest_dir,grid_code,cfg=None,suffix=None):
    """Cleans the geojson perimeters in a finished run directory and builds the kml.

    With a suffix the parameter set is appended to every output name before the
    extension -- sim_<grid>_final_07_windReductionFactor_0.6.geojson -- so files
    from different parameter sets stay apart when they are opened side by side.
    The kml placemarks take their names from the geojson filenames, so tagging
    the files labels the layers in Google Earth as well.

    Split out of cleanup_ff_run so it can be re-run over results already in a
    workspace without repeating the simulation.
    """
    if cfg is None:
        cfg = load_ff_cfg()

    #rewrite the geojson files in the dest_dir to remove the z-coordiate
    g = glob.glob(f"{dest_dir}/*.geojson")
    for gg in g:
        geo_in = json.load(open(gg))
        geo_out = remove_z_coordinate(geo_in)
        geo_out['Name'] = f"{grid_code}_{suffix}" if suffix else grid_code
        if "final" in os.path.basename(gg):
            geo_out['Description'] = 'Final Perimeters'
        else:
            geo_out['Description'] = 'Time Step Perimeters'
        out_path = gg
        if suffix and not gg.endswith(f"_{suffix}.geojson"):
            out_path = f"{gg[:-len('.geojson')]}_{suffix}.geojson"
        with open(out_path,"w") as w_file:
            json.dump(geo_out,w_file,indent = 2)
        if out_path != gg:
            os.remove(gg)

    #make kml files
    tag = f"_{suffix}" if suffix else ""
    steps_file = f"{dest_dir}/{grid_code}_steps{tag}.kml"
    merge_geojson_to_kml(dest_dir,steps_file,final = False,cfg = cfg,grid_code = grid_code)
    final_file = f"{dest_dir}/{grid_code}_final{tag}.kml"
    merge_geojson_to_kml(dest_dir,final_file,final = True,cfg = cfg,grid_code = grid_code)

def prob_forecast(json_directory):
    import glob
    import geopandas as gpd
    import pandas as pd
    from shapely.ops import polygonize



    # 1. Load and combine all 50 GeoJSON files
    geojson_files = glob.glob(f"{json_directory}/*final*.geojson")
    n_ensembles = len(geojson_files)

    gdf_list = []
    for file in geojson_files:
        temp_gdf = gpd.read_file(file)
        gdf_list.append(temp_gdf)

    gdf = pd.concat(gdf_list, ignore_index=True)

    # 2. Explicitly define the input CRS for v0.6.1
    gdf.crs = {'init': 'epsg:4326'}

    # 3. Project to a planar coordinate system for accurate simplification
    # Using EPSG:3857 (Web Mercator) as a generic backup, or replace with your local UTM projection
    gdf = gdf.to_crs({'init': 'epsg:3857'})

    # 4. Explode MultiPolygons into single Polygons (v0.6.1 syntax)
    # In 0.6.1, explode returns a multi-index Series/DataFrame that we must reset
    gdf = gdf.explode()
    gdf = gdf.reset_index(drop=True)

    # 5. Simplify the geometries to speed up line intersections
    # Tolerance is now in meters (e.g., 10 meters)
    gdf["geometry"] = gdf["geometry"].simplify(tolerance=10, preserve_topology=True)



    # --- NEW: VALIDATION & REPAIR STEP ---
    # Check for invalid geometries
    invalid_mask = ~gdf.is_valid
    if invalid_mask.any():
        print(f"Found {invalid_mask.sum()} invalid geometries. Repairing with buffer(0)...")
        
        # Fix invalid polygons by applying a zero-buffer
        gdf.loc[invalid_mask, "geometry"] = gdf.loc[invalid_mask, "geometry"].buffer(0)
        
        # Drop any geometries that failed to repair or became empty
        gdf = gdf[gdf.is_valid & ~gdf.is_empty].reset_index(drop=True)
    # -------------------------------------




    # 6. Extract boundaries and combine them
    boundaries = gdf.boundary.unary_union

    # 7. Rebuild the atomic puzzle pieces (fragments)
    atomic_polygons = list(polygonize(boundaries))
    fragments_gdf = gpd.GeoDataFrame(geometry=atomic_polygons, crs=gdf.crs)

    # 8. Keep only fragments that sit inside the actual fire zones
    total_coverage = gdf.unary_union
    fragments_gdf = fragments_gdf[fragments_gdf.intersects(total_coverage)].copy()

    # 9. Count how many original perimeters cover each fragment
    # Generate an internal representative point for each fragment to perform a robust point-in-polygon join
    fragments_gdf["sample_point"] = fragments_gdf.geometry.representative_point()

    # Workaround for v0.6.1 spatial join using the point geometry
    points_gdf = fragments_gdf.copy()
    points_gdf = points_gdf.set_geometry("sample_point")

    # Spatial join: match sample points inside the original perimeters
    joined = gpd.sjoin(points_gdf, gdf, how="inner", op="within")

    # 10. Map counts back to fragments and calculate final probability
    counts = joined.index.value_counts()
    fragments_gdf["count"] = fragments_gdf.index.map(counts).fillna(0)
    fragments_gdf["probability"] = fragments_gdf["count"] / n_ensembles

    # Clean up temporary point geometry column
    fragments_gdf = fragments_gdf.drop(columns=["sample_point"])

    # 11. Classify the probabilities into discrete operational bands
    # Define your bin edges and matching clear labels
    bins = [-1, 0.5, 0.75, 0.9, 1.01]  # -1 handles exactly 0, 1.01 handles exactly 1.0
    labels = ["P < 0.50", "0.50 <= P < 0.75", "0.75 <= P < 0.90", "P >= 0.90"]

    # Assign each fragment to a probability band
    fragments_gdf["risk_zone"] = pd.cut(fragments_gdf["probability"], bins=bins, labels=labels)

    # Group and merge (union) all adjacent or separated geometries sharing the same risk zone
    prob_map = fragments_gdf.dissolve(by="risk_zone").reset_index()

    # Clean up and reorder columns for clarity
    prob_map = prob_map[["risk_zone", "geometry"]]

    # --- FIX: Convert Pandas Category to String for Fiona compatibility ---
    prob_map["risk_zone"] = prob_map["risk_zone"].astype(str)
    # ----------------------------------------------------------------------


    # 12. Convert back to EPSG:4326 before saving to GeoJSON
    prob_map = prob_map.to_crs({'init': 'epsg:4326'})
    # Save final probabilistic contour map
    prob_file = f"{json_directory}/probabilistic_wildfire_map.geojson"
    prob_map.to_file(prob_file, driver="GeoJSON")

    # --- STEP 13: STYLING FOR KML EXPORT ---
    # Define a color scheme (HEX values) from low risk (yellow) to high risk (dark red)
    color_mapping = {
        "P < 0.50": "#FFFF00",          # Yellow
        "0.50 <= P < 0.75": "#FFA500",  # Orange
        "0.75 <= P < 0.90": "#FF4500",  # Orange-Red
        "P >= 0.90": "#8B0000"          # Dark Red
    }

    # Add a name column (Google Earth uses the 'Name' or 'description' field for labels)
    prob_map["Name"] = prob_map["risk_zone"]
    prob_map["description"] = prob_map["risk_zone"]

    # Inject standard SimpleStyle properties that many modern KML/GeoJSON viewers read
    prob_map["fill"] = prob_map["risk_zone"].map(color_mapping)
    prob_map["stroke"] = "#FFFFFF"  # White borders for contrast
    prob_map["fill-opacity"] = 0.5   # Semi-transparent so map layers underneath are visible

    # Force file into EPSG:4326 (KML requires latitude/longitude)
    prob_map_kml = prob_map.to_crs({'init': 'epsg:4326'})

    # Enable the KML driver in Fiona (disabled by default in older setups)
    import fiona
    fiona.drvsupport.supported_drivers['KML'] = 'rw'
    fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'

    # Determine the best available driver
    kml_driver = 'KML'

    #if 'LIBKML' in fiona.available_drivers:
    #    kml_driver = 'LIBKML'

    kml_file = prob_file.replace('.geojson','.kml')
    try:
        prob_map_kml.to_file(kml_file, driver=kml_driver)
        print(f"Successfully saved styled KML using {kml_driver} driver.")
    except Exception as e:
        print(f"KML Driver failed: {e}. Falling back to styled GeoJSON.")





def merge_geojson_to_kml(input_folder, output_kml_path,final=False,cfg=None,grid_code=None):


    # Enable KML driver in fiona (it is disabled by default for writing)
    fiona.drvsupport.supported_drivers['KML'] = 'rw'

    geojson_files = [f for f in os.listdir(input_folder) if f.endswith('.geojson')]
    
    if not geojson_files:
        print("No GeoJSON files found in the folder.")
        return
        
    gdf_list = []
    
    # Read each GeoJSON file
    for file in geojson_files:
        if final and not "final" in file:
            continue
        if not final and "final" in file:
            continue
        file_path = os.path.join(input_folder, file)
        try:
            gdf = gpd.read_file(file_path)
            # Ensure everything uses the standard KML coordinate system (WGS84)
            if gdf.crs is None:
                gdf.set_crs(epsg=4326, inplace=True)
            else:
                gdf = gdf.to_crs(epsg=4326)
            
            # Optional: Add a column to track which file the feature came from
            gdf['source_file'] = file
            gdf_list.append(gdf)
            print(f"Loaded: {file}")
        except Exception as e:
            print(f"Error loading {file}: {e}")

    #the guard above only catches an empty folder; this catches a folder whose
    #geojsons are all the other kind, which happens when a run produced step
    #perimeters but no finals
    if not gdf_list:
        print(f"No {'final' if final else 'time step'} GeoJSON files to merge, skipping {output_kml_path}")
        return

    # Combine all layers into one Geodataframe
    combined_gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True), crs=gdf_list[0].crs)
    
    # Name each placemark for the file it came from, so the perimeters stay
    # identifiable in the Google Earth sidebar. Everything else is dropped: the
    # kml driver names a placemark from the first attribute field, which used to
    # be ForeFire's numberOfPolygons and made every placemark read as "1".
    names = combined_gdf['source_file'].str.replace('.geojson','',regex=False)
    if grid_code:
        names = names.str.replace(f"sim_{grid_code}_",'',regex=False)
    combined_gdf['Name'] = names
    combined_gdf = combined_gdf[['Name','geometry']]
        
    print(f"Exporting {len(combined_gdf)} features to {output_kml_path}...")
    
    # Write to KML
    #delete any file that is there already
    if os.path.exists(output_kml_path):
        os.system(f"rm {output_kml_path}")
    combined_gdf.to_file(output_kml_path, driver='KML')
    #copy to the perimeter collection directory, if one is configured
    if cfg is None:
        cfg = load_ff_cfg()
    perim_dir = cfg.get('perim_copy_dir')
    if perim_dir:
        if os.path.isdir(perim_dir):
            os.system(f"cp {output_kml_path} {perim_dir}/.")
        else:
            print(f"perim_copy_dir {perim_dir} does not exist, not copying the kml there")
    print("Success!")







def remove_z_coordinate(geojson_data):
  """Removes the z-coordinate (altitude) from a GeoJSON object in-place or returns a clean copy."""

  def clean_coords(coords):
    if not coords:
      return coords
    # Check if the first item is a number (a single coordinate pair/triple)
    if isinstance(coords[0], (int, float)):
      return coords[:2]
    # Otherwise, recurse deeper into nested coordinate arrays
    return [clean_coords(c) for c in coords]

  if "features" in geojson_data:
    for feature in geojson_data["features"]:
      geom = feature.get("geometry")
      if geom and "coordinates" in geom:
        geom["coordinates"] = clean_coords(geom["coordinates"])
  elif "coordinates" in geojson_data:
    geojson_data["coordinates"] = clean_coords(geojson_data["coordinates"])

  return geojson_data


def run_forecasts(wksp_dir,forefire_dir = None, overwrite = False, cfg = None, params = None):
    #forefire_dir defaults to run_dir in etc/forefire.json
    #params overrides ForeFire setParameter values for this run
    if cfg is None:
        cfg = load_ff_cfg()
    if forefire_dir is None:
        forefire_dir = cfg['run_dir']
    #fail now, not after an hour of netcdf building, if the container cannot see it
    container_path(forefire_dir,cfg)
    #read information about the fire
    t0 = pd.Timestamp.now()
    start_utc,ign_utc,ign_latlon,grid_code = read_input(wksp_dir)

    #find list of wrfouts and pair with times
    timing_table = make_timing_table(wksp_dir,ign_utc)     #<<<<-------------------------- add ability to make forecasts from the weather files

    #clean the forefire drectory
    if not os.path.exists(forefire_dir):
        os.makedirs(forefire_dir)
    #remove the files, leaving any parameter-sweep subdirectories in place
    for stale in glob.glob(f"{forefire_dir}/*"):
        if os.path.isfile(stale) or os.path.islink(stale):
            os.remove(stale)
    ##### this should run in an update mode, always making forecasts that don't already exist
    if len(timing_table) > cfg.get('min_timesteps',20): # and not os.path.exists(f"{wksp_dir}/forefire"):
        #dead fuel moisture from FMDA, if configured; an explicit params entry wins
        run_params = dict(moisture_params(wksp_dir,ign_latlon,cfg))
        run_params.update(params or {})
        #generate set of NC files and scripts for ForeFire
        timing_table = make_script_set(forefire_dir,timing_table,ign_latlon,grid_code,cfg=cfg,params=run_params)
        timing_table.to_csv(f"{forefire_dir}/timing_table.csv",index=False)
        #run the forecasts
        run_timing_table(timing_table,overwrite=overwrite,cfg=cfg,run_dir=forefire_dir)
        #make forefire directory in the wksp_dir and move all input files there
        cleanup_ff_run(forefire_dir,wksp_dir,grid_code,cfg=cfg)
    else:
        print('Not runnning ForeFire forecasts, too few wrfout files or forecasts already made')

    
    t1 = pd.Timestamp.now()
    print(t1)
    dt = t1-t0
    print(f"Finished running in {dt.total_seconds()/60} minutes")


def run_days(days2run=1,overwrite=False,cfg=None):
    #overwrite=True re-runs time steps whose final geojson is already in the wksp
    if cfg is None:
        cfg = load_ff_cfg()
    wksp_root = cfg.get('wksp_root','/data/jhaley/wrfxpy/wksp')

    #empty list of wksp directories
    wksp_dirs = []
    now = pd.Timestamp.now()
    for i in range(days2run):
        dt = now - timedelta(days=i)
        date_string = dt.strftime('%Y-%m-%d')
        g = glob.glob(f'{wksp_root}/*{date_string}*')
        wksp_dirs.extend(g)
        print(f"{len(g)} wrksp directoriess for {date_string}")

    for gg in g:
        print(f"Running forecasts for {gg}")
        run_forecasts(gg,overwrite=overwrite,cfg=cfg)

    
    '''
    # Use "w" mode to overwrite, and always specify the encoding
    with open("make_forefire_perims.sh", "w", encoding="utf-8") as file:
        file.write('#!/usr/bin/env bash \n')
        file.write('source /home/jhaley/.bashrc \n')
        file.write('conda activate wrf_test \n')
        file.write('PYTHONPATH=src \n')
        for wd in wksp_dirs:
            file.write(f"python src/ingest/forefire.py {wd} \n")  # Adds a newline character at the end of each item
    '''




if __name__ == "__main__":
    run_days(days2run=40)
    '''
    l = [
        '/data/jhaley/wrfxpy/wksp/wfc-LITTLE_GIANT_2026-07-16_09_00_00_091081ED-BD23-4610-AE4A-270F95D1711E-2026-07-16_09:00:00-27',
        '/data/jhaley/wrfxpy/wksp/wfc-SISI_2026-08-19_09_00_00_DC4342D9-B479-44F1-906C-8ABD42E1F59C-2026-08-19_09:00:00-27'
        ]
    for ll in l:
        run_forecasts(ll)
    #hill = "/home/jhaley/work/WRF/test/em_fire/hill/wrfout_d01_0001-01-01_00:00:00"
    #outfile = "/home/jhaley/forefire/tests/hill.nc"
    #make_FF_nc(hill,outfile)

    
    t0 = pd.Timestamp.now()
    print('Starting FF script',t0)
    overwrite = False
    print(sys.argv)

    if len(sys.argv) > 1:
        print(f"Will make forecasts for {sys.argv[1]}")
        wksp_dir = sys.argv[1]
        if "overwrite" in sys.argv:
            overwrite = True
    else:
        #wksp_dir = '/data/jhaley/wrfxpy/wksp/wfc-China_2026-08-12_02_00_00_C0313CBB-4DEF-40AC-81A7-CF9483683F36-2026-08-12_00:00:00-30'
        #wksp_dir = '/data/jhaley/wrfxpy/wksp/wfc-0722_ROSENBAUM_2026-08-14_02_00_00_2059A645-6D4B-4EF9-8D1A-8EA49732A456-2026-08-14_00:00:00-30'
        wksp_dir = '/data/jhaley/wrfxpy/wksp/wfc-CHURCH_2_2026-08-17_20_00_00_112B9B87-43FE-49C0-98C3-CBC6277625C4-2026-08-17_18:00:00-30'
        #wksp_dir = '/data/jhaley/wrfxpy/wksp/wfc-Keyserville_2026-08-17_20_00_00_F5F1F398-0B0F-42C7-8771-FD8D6324F73A-2026-08-17_18:00:00-30'
        #forefire_dir = '/home/jhaley/forefire/tests/china'

    #maybe change this, this directory gets bound to the apptainer files system
    forefire_dir = '/home/jhaley/forefire/tests/ffwksp'

    run_forecasts(wksp_dir,forefire_dir,overwrite=overwrite)
    '''

    '''
    #read information about the fire
    start_utc,ign_utc,ign_latlon,grid_code = read_input(wksp_dir)

    #find list of wrfouts and pair with times
    timing_table = make_timing_table(wksp_dir,ign_utc)     #<<<<-------------------------- add ability to make forecasts from the weather files

    ##### this should run in an update mode, always making forecasts that don't already exist
    if len(timing_table) > 20: # and not os.path.exists(f"{wksp_dir}/forefire"):
        #generate set of NC files and scripts for ForeFire
        timing_table = make_script_set(forefire_dir,timing_table,ign_latlon,grid_code)
        timing_table.to_csv(f"{forefire_dir}/timing_table.csv",index=False)
        #run the forecasts
        run_timing_table(timing_table)
        #make forefire directory in the wksp_dir and move all input files there
        cleanup_ff_run(forefire_dir,wksp_dir,grid_code)
    else:
        print('Not runnning ForeFire forecasts, too few wrfout files or forecasts already made')

    
    t1 = pd.Timestamp.now()
    print(t1)
    dt = t1-t0
    print(f"Finished running in {dt.total_seconds()/60} minutes")
    '''
    
