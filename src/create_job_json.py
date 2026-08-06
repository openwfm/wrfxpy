# Libraries #
from datetime import datetime, timedelta, UTC
from utils import load_sys_cfg
import os.path as osp
import numpy as np
import json
import os

sys_cfg = load_sys_cfg()
clusters_path = osp.join(sys_cfg["sys_install_path"], "etc/clusters.json")
clusters = json.load(open(clusters_path))

grib_sources = {
    "GFSF": 28, "RAP": 13, "RRFSNA": 3,
    "NAM218": 12, "NAM227": 5, "RRFS": 3, 
    "HRRR": 3
} # resolutions in km (roughly)

profile_sizes = {
    "S": 72, "M": 96, "L": 144, "XL": 192
}

# Functions #
# Replaces the template json values with the user defined values
#
# @Param tmp_dict - The original json dictionary values
# @Param usr_dict - The user"s json dictionary values
#
def replace_values(tmp_dict, usr_dict):
    for key, value in tmp_dict.items():
        if isinstance(value, dict):
            replace_values(value, usr_dict)
        elif key in usr_dict:
            if isinstance(usr_dict[key], (dict, list)):  # Check for nested values
                tmp_dict[key] = usr_dict[key]
            else:
                tmp_dict[key] = usr_dict[key]
        elif isinstance(value, list):  # Check for lists
            for item in value:
                if isinstance(item, dict):
                    replace_values(item, usr_dict)
                    
def make_profile_table(profile_sizes, inner_res):
    headers = []
    grid_sizes = []
    physical_sizes = []

    for label, n in profile_sizes.items():
        size_km = round(n * inner_res / 1000)

        headers.append(label)
        grid_sizes.append(f"{n}x{n}")
        physical_sizes.append(f"{size_km}x{size_km}km")

    rows = [headers, grid_sizes, physical_sizes]

    # Compute column widths dynamically
    col_widths = [
        max(len(str(row[i])) for row in rows) + 2
        for i in range(len(headers))
    ]

    def format_row(row):
        return "|" + "|".join(
            f"{str(cell):^{col_widths[i]}}"
            for i, cell in enumerate(row)
        ) + "|"

    return "\n".join(format_row(row) for row in rows)

OUTER_POSTPROC = [
    "T2_F", "PSFC", "WINDSPD_mph", "WINDVEC", "WINDVEC_mph_D",
    "PM25_SFC_D", "SMOKE1000FT_AGL_D", "SMOKE4000FT_AGL_D",
    "SMOKE6000FT_AGL_D", "VR_SFC", "VR1000FT_AGL", "VR4000FT_AGL",
    "VR6000FT_AGL", "PLUME_HEIGHT_kft", "VR_SFC",
]

INNER_POSTPROC = [
    "T2_F", "PSFC", "WINDSPD_mph", "WINDVEC", "WINDVEC_mph_D",
    "RH_FIRE", "ZSF", "NFUEL_CAT", "1HR_FM", "10HR_FM", "100HR_FM",
    "PM25_SFC_D", "PLUME_HEIGHT_kft", "SMOKE1000FT_AGL_D",
    "SMOKE4000FT_AGL_D", "SMOKE6000FT_AGL_D", "VR_SFC",
    "VR1000FT_AGL", "VR4000FT_AGL", "VR6000FT_AGL", "HDW", "FOSBERG",
    "FIRE_AREA", "FGRNHFX", "FCANHFX"
]

def build_domain_conf(grib_res, n_domains, profile_size, center_latlon, ref_ratio=3, subgrid_res=30):    
    domain_size = profile_size + 1
    parent_start = profile_size // ref_ratio + 1
    parent_end = parent_start + profile_size // ref_ratio - 1
    cell_size = int(round(grib_res * 1000))
    subgrid_ratio = int(np.ceil(cell_size / ref_ratio**n_domains / subgrid_res / 2) * 2)

    cell_sub_size = int(round(cell_size / ref_ratio))
    domains = {
        "1": {
            "cell_size": [cell_sub_size, cell_sub_size],
            "domain_size": [domain_size, domain_size],
            "center_latlon": center_latlon,
            "truelats": [center_latlon[0], center_latlon[0]],
            "stand_lon": center_latlon[1],
            "time_step": max(1, int(6 * cell_sub_size / 1000) - 1),
            "history_interval": 60,
            "geog_res": "30s",
            "subgrid_ratio": [0, 0],
        }
    }
    for dom in range(2, n_domains + 1):
        cell_sub_size /= ref_ratio
        domains[str(dom)] = {
            "parent_id": dom - 1,
            "parent_cell_size_ratio": ref_ratio,
            "parent_time_step_ratio": ref_ratio,
            "geog_res": "30s", #".3s" if cell_sub_size < 1000. else "30s",
            "subgrid_ratio": [
                subgrid_ratio, subgrid_ratio
            ] if dom == n_domains else [0, 0],
            "parent_start": [parent_start, parent_start],
            "parent_end": [parent_end, parent_end],
            "history_interval": 15 if dom == n_domains else 60,
        }

    return domains

def build_job_json(cfg, grib_res, n_domains, profile_size):
    outer_dom = 1
    inner_dom = n_domains
    fmda_utc = datetime.strptime(cfg["start_utc"], "%Y-%m-%d_%H:%M:%S").replace(tzinfo=UTC)
    patch_size = 8
    n_cores = (profile_size / patch_size)**2
    ppn = clusters[sys_cfg["qsys"]].get("ppn", 64)
    n_nodes = int(np.floor(n_cores / ppn))
    if n_nodes == 0:
        n_nodes = 1
        ppn = n_cores

    return {
        "run_wrf": cfg["run_wrf"],
        "use_realtime": True,
        "iofields": True,
        "use_wgrib2": True,
        "grid_code": cfg["grid_code"],
        "grib_source": cfg["grib_source"],
        "wps_namelist_path": "etc/nlists/default.wps",
        "wrf_namelist_path": "etc/nlists/default.input",
        "fire_namelist_path": "etc/nlists/default.fire",
        "emissions_namelist_path": "etc/nlists/default.fire_emissions",
        "geo_vars_path": cfg["geo_vars_path"],
        "num_nodes": n_nodes,
        "ppn": ppn,
        "wall_time_hrs": 12,
        "start_utc": cfg["start_utc"],
        "end_utc": cfg["end_utc"],
        "cycle_start_utc": cfg["cycle_start_utc"],
        "domains": build_domain_conf(
            grib_res,
            n_domains,
            profile_size,
            cfg["center_latlon"],
        ),
        "ignitions": {},
        "postproc": {
            str(outer_dom): OUTER_POSTPROC,
            str(inner_dom): INNER_POSTPROC,
            "shuttle": "incremental",
            "description": cfg["description"],
        },
        "fmda_geogrid_path": fmda_utc.strftime(
            f"wksp_fmda/CONUS/%Y%m/fmda-CONUS-%Y%m%d-%H/fmda-CONUS-%Y%m%d-%H.geo"
        ),
    }

################################################################################
# Create dictionary to hold user data and define current datetime in Local and UTC #
cfg = {} # dictionary to hold user inputs

################################################################################
# Get the grib source from the user #
print(
    "Long-range forecast lead-time limitations by forcing:\n"
    " - GFSF (00/06/12/18 UTC): 384 h: 5h latency\n"
    " - RAP (03/09/15/21 UTC): 51 h: 1h latency\n"
    " - NAM218 (00/06/12/18 UTC): 84 h: 3h latency\n"
    " - NAM227 (00/06/12/18 UTC): 60 h: 3h latency\n"
    " - RRFSNA (00/06/12/18 UTC): 84 h: 4h latency\n"
    " - RRFS (00/06/12/18 UTC): 84 h: 4h latency\n"
    " - HRRR (00/06/12/18 UTC): 48 h: 2h latency\n"
)
while True:
    grib_sources_input = input("Enter a grib source {} (HRRR): ".format(grib_sources.keys()))
    if grib_sources_input == "":
        grib_sources_input = "HRRR"
    grib_sources_input = grib_sources_input.upper()
    if not grib_sources_input in grib_sources.keys():
        print("Please enter a valid grib source.")
        continue
    grib_res = grib_sources[grib_sources_input]
    print(f"You have selected {grib_sources_input} ({grib_res} km)")
    print()
    cfg["grib_source"] = grib_sources_input
    break

################################################################################
# Get domain configuration and resolution based on target resolution
while True:
    inner_res = input("Select innermost atmospheric resolution in meters (300 m): ")
    if inner_res == "":
        inner_res = 300
    try:
        inner_res = float(inner_res)
    except:
        print("Please select a valid resolution")
        continue
    break
    
ref_ratio = 3
n_domains = round(np.log(grib_res * 1000. / inner_res) / np.log(ref_ratio))
inner_res = grib_res * 1000. / ref_ratio**n_domains
print(f"Domain configuration with {n_domains} domains and inner most resolution of {inner_res} m\n")

################################################################################
# Get the profile size input from the user #
while True:
    print("Select size:")
    print(make_profile_table(profile_sizes, inner_res))
    profile_size_input = input("Enter a valid profile size [S, M, L, XL] (M): ")
    if profile_size_input == "":
        profile_size_input = "M"
    profile_size_input = profile_size_input.upper()
    if not profile_size_input in profile_sizes:
        print("Please select a valid profile size")
        continue
    print(f"You have selected the {profile_size_input} profile size")
    print()
    break

################################################################################
# Get the simulation start time from the user 
while True:
    utc_now = datetime.now(UTC)
    utc_now_str = utc_now.strftime("%Y-%m-%d_%H:%M:%S")
    start_utc = input('Current UTC time is {}. Enter start time in UTC with format "%Y-%m-%d_%H:%M:%S": '.format(utc_now_str))
    try:
        start_utc = datetime.strptime(start_utc, "%Y-%m-%d_%H:%M:%S").replace(tzinfo=UTC)
    except:
        print(f"Enter a valid UTC time, not {start_utc}")
        continue
    if abs(utc_now - start_utc) > timedelta(days=7):
        print("WARNING: this is for forecasting operations, you might need to consider modifying some flags.")
    break

print(start_utc.strftime("Simulation Start Time: %Y-%m-%d_%H:%M:%S"))
print()
cfg["start_utc"] = start_utc.strftime("%Y-%m-%d_%H:%M:%S")

################################################################################
# Get the end time and ignition time #
while True:
    time_input = input("How many hours would you like to run the simulation for (48): ")
    if time_input == "":
        time_input = 48
    try:
        hours = int(time_input)
    except:
        print(f"Please enter a valid hours value, not {time_input}")
        continue
    break

end_utc = start_utc + timedelta(hours=hours)
cfg["end_utc"] = end_utc.strftime("%Y-%m-%d_%H:%M:%S")
print(f"The end time of the simulation: {end_utc}")
print()

################################################################################
# Get the cycle start time from the user 
while True:
    utc_now = datetime.now(UTC)
    utc_now_str = utc_now.strftime("%Y-%m-%d_%H:%M:%S")
    start_utc_str = start_utc.strftime("%Y-%m-%d_%H:%M:%S")
    cycle_start_utc = input('Enter the cycle start time in UTC with format "%Y-%m-%d_%H:%M:%S" ({}): '.format(start_utc_str))
    if cycle_start_utc == '':
        cycle_start_utc = start_utc_str
    try:
        cycle_start_utc = datetime.strptime(cycle_start_utc, "%Y-%m-%d_%H:%M:%S").replace(tzinfo=UTC)
    except:
        print(f"Enter a valid UTC time, not {cycle_start_utc}")
        continue
    if cycle_start_utc > start_utc:
        print(f"Cycle start time {cycle_start_utc} needs to be <= start time {start_utc}")
        continue
    break

print(cycle_start_utc.strftime("Simulation Cycle Start Time: %Y-%m-%d_%H:%M:%S"))
print()
cfg["cycle_start_utc"] = cycle_start_utc.strftime("%Y-%m-%d_%H:%M:%S")

################################################################################
# Get the central latitude and longitude values
while True:
    coord_input = input(
        "Enter longitude, latitude of the center of the fire as decimal coordinates: "
    )
    # Remove all spaces, then split on comma
    coord_input_clean = coord_input.replace(" ", "")
    try:
        lon_input, lat_input = map(float, coord_input_clean.split(","))
    except ValueError:
        print(f"Enter valid coordinate values, not {coord_input}")
        continue
    # Validate coordinate ranges
    if not (-90.0 <= lat_input <= 90.0):
        print(f"Latitude must be between -90 and 90 degrees, not {lat_input}")
        continue
    if not (-180.0 <= lon_input <= 180.0):
        print(f"Longitude must be between -180 and 180 degrees, not {lon_input}")
        continue
    break
    
cfg["center_latlon"] = [lat_input, lon_input]
cfg["truelats"] = [lat_input, lat_input]
cfg["stand_lon"] = lon_input
cfg["latlon"] = [lat_input, lon_input]
print()

################################################################################
# Get the fire name
fire_name_input = input("Enter the name of the fire (Test): ")
if fire_name_input == "":
    fire_name_input = "Test"
print(f"You have entered {fire_name_input} as the fire name")
print()

################################################################################
# Get if run or not the whole workflow
while True:
    run_conus = input("Is the fire from CONUS? [Y/N] (Y): ")
    if run_conus == "":
        run_conus = "Y"
    if run_conus.upper() in ["Y", "N"]:
        print(f"You have entered {run_conus} to run a CONUS wildfire")
        break
    print(f"Enter a valid answer (Y/N), not {run_conus}")
print()

################################################################################
# Get if run or not the whole workflow
while True:
    run_wrf = input("Do you want to run the whole workflow? [Y/N] (Y): ")
    if run_wrf == "":
        run_wrf = "Y"
    if run_wrf.upper() in ["Y", "N"]:
        print(f"You have entered {run_wrf} to run the whole workflow")
        break
    print(f"Enter a valid answer (Y/N), not {run_wrf}")
print()

################################################################################
# Get if generate overnight json
while True:
    generate_overnight = input("Do you want to generate an overnight JSON? [Y/N] (Y): ")
    if generate_overnight == "":
        generate_overnight = "Y"
    if generate_overnight.upper() in ["Y", "N"]:
        print(f"You have entered {generate_overnight} to generate an overnight JSON")
        break
    print(f"Enter a valid answer (Y/N), not {generate_overnight}")
print()

################################################################################
# Get the code for the fire
cfg["grid_code"] = datetime.now(UTC).strftime(f"%Y-%m-%d_%H-{fire_name_input.replace(" ", "_").upper()}")

################################################################################
# Construct the description and fmda path for the json #
cfg["description"] = start_utc.strftime(f"{fire_name_input.title()} {profile_size_input} %Y-%m-%d %Hz")
cfg["fmda_geogrid_path"] = start_utc.strftime("wksp_fmda/CONUS/%Y%m/fmda-CONUS-%Y%m%d-%H.geo")

################################################################################
# Modify other flags #
if run_conus.upper() == "Y":
    cfg["geo_vars_path"] = "etc/vtables/geo_vars.json"
else:
    cfg["geo_vars_path"] = "etc/vtables/geo_vars_canada.json"
    
if run_wrf.upper() == "Y":
    cfg["run_wrf"] = True
else:
    cfg["run_wrf"] = False

################################################################################
# Generate the json #
fire_name_modified = fire_name_input.replace(" ", "_")
job_id = start_utc.strftime("{}_{}_%y%m%d_%Hz".format(fire_name_modified, profile_size_input))
destination_path = osp.join("jobs", job_id + ".json")

json_data = build_job_json(
    cfg,
    grib_res,
    n_domains,
    profile_sizes[profile_size_input],
)

# Finalize the changes to the json
os.makedirs(os.path.dirname(destination_path), exist_ok=True)
with open(destination_path, "w") as json_file:
    json.dump(json_data, json_file, indent=4)
    
# Create overnight json
if generate_overnight.upper() == "Y":
    overnight_hour = 8
    json_data_overnight = json_data.copy()
    start_utc_overnight = (start_utc + timedelta(days=1)).replace(hour=overnight_hour)
    json_data_overnight["grid_code"] = start_utc_overnight.strftime(f"%Y-%m-%d_%H-{fire_name_input.replace(' ', '_').upper()}")
    json_data_overnight["start_utc"] = start_utc_overnight.strftime("%Y-%m-%d_%H:%M:%S")
    json_data_overnight["end_utc"] = end_utc.replace(hour=6).strftime("%Y-%m-%d_%H:%M:%S")
    json_data_overnight["cycle_start_utc"] = start_utc_overnight.replace(hour=6).strftime("%Y-%m-%d_%H:%M:%S")
    json_data_overnight["fmda_geogrid_path"] = start_utc_overnight.strftime("wksp_fmda/CONUS/%Y%m/fmda-CONUS-%Y%m%d-%H/fmda-CONUS-%Y%m%d-%H.geo")
    json_data_overnight["postproc"]["description"] = start_utc_overnight.strftime(f"{fire_name_input.title()} {profile_size_input} %Y-%m-%d %Hz")
    overnight_job_id = start_utc_overnight.strftime("{}_{}_%y%m%d_%Hz".format(fire_name_modified, profile_size_input))
    overnight_destination_path = osp.join("jobs", overnight_job_id + ".json")
    with open(overnight_destination_path, "w") as json_file:
        json.dump(json_data_overnight, json_file, indent=4)

print(f"JSON is ready here: {destination_path}")
print(f"./forecast.sh {destination_path} >& logs/{job_id}.log")
