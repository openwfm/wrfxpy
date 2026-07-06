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

# Functions #
# Replaces the template json values with the user defined values
#
# @Param tmpDict - The original json dictionary values
# @Param usrDict - The user"s json dictionary values
#
def replace_values(tmpDict, usrDict):
    for key, value in tmpDict.items():
        if isinstance(value, dict):
            replace_values(value, usrDict)
        elif key in usrDict:
            if isinstance(usrDict[key], (dict, list)):  # Check for nested values
                tmpDict[key] = usrDict[key]
            else:
                tmpDict[key] = usrDict[key]
        elif isinstance(value, list):  # Check for lists
            for item in value:
                if isinstance(item, dict):
                    replace_values(item, usrDict)
                    
def make_profile_table(profileSizes, innerRes):
    headers = []
    grid_sizes = []
    physical_sizes = []

    for label, n in profileSizes.items():
        size_km = round(n * innerRes / 1000)

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
    "VR6000FT_AGL", "PLUME_HEIGHT_ft", "VR_SFC",
]

INNER_POSTPROC = [
    "T2_F", "PSFC", "WINDSPD_mph", "WINDVEC", "WINDVEC_mph_D",
    "RH_FIRE", "ZSF", "NFUEL_CAT", "1HR_FM", "10HR_FM", "100HR_FM",
    "PM25_SFC_D", "PLUME_HEIGHT_ft", "SMOKE1000FT_AGL_D",
    "SMOKE4000FT_AGL_D", "SMOKE6000FT_AGL_D", "VR_SFC",
    "VR1000FT_AGL", "VR4000FT_AGL", "VR6000FT_AGL", "HDW", "FOSBERG",
    "FIRE_AREA", "FGRNHFX", "FLINEINT_btupftps", "ROS_chsph",
]

def build_domain_conf(gribRes, n_domains, profileSize, center_latlon, ref_ratio=3, subgrid_res=30):    
    domain_size = profileSize + 1
    parent_start = profileSize // ref_ratio + 1
    parent_end = parent_start + profileSize // ref_ratio - 1
    cell_size = int(round(gribRes * 1000))
    subgrid_ratio = int(np.ceil(cell_size / 3**(n_domains - 1) / subgrid_res / 2) * 2)

    domains = {
        "1": {
            "cell_size": [cell_size, cell_size],
            "domain_size": [domain_size, domain_size],
            "center_latlon": center_latlon,
            "truelats": [center_latlon[0], center_latlon[0]],
            "stand_lon": center_latlon[1],
            "time_step": max(1, int(6 * gribRes)),
            "history_interval": 60,
            "geog_res": "30s",
            "subgrid_ratio": [0, 0],
        }
    }
    dom_size = cell_size / ref_ratio
    for dom in range(2, n_domains + 1):
        domains[str(dom)] = {
            "parent_id": dom - 1,
            "parent_cell_size_ratio": ref_ratio,
            "parent_time_step_ratio": ref_ratio,
            "geog_res": ".3s" if dom_size < 1000. else "30s",
            "subgrid_ratio": [
                subgrid_ratio, subgrid_ratio
            ] if dom == n_domains else [0, 0],
            "parent_start": [parent_start, parent_start],
            "parent_end": [parent_end, parent_end],
            "history_interval": 15 if dom == n_domains else 60,
        }
        dom_size /= ref_ratio

    return domains

def build_job_json(cfg, gribRes, n_domains, profileSize, profileSizeInput, fireNameInput):
    outer_dom = 1
    inner_dom = n_domains
    fmdaUTC = datetime.strptime(cfg["start_utc"], "%Y-%m-%d_%H:%M:%S").replace(tzinfo=UTC)
    patch_size = 8
    n_cores = (profileSize / patch_size)**2
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
        "num_nodes": n_nodes,
        "ppn": ppn,
        "wall_time_hrs": 12,
        "start_utc": cfg["start_utc"],
        "end_utc": cfg["end_utc"],
        "domains": build_domain_conf(
            gribRes,
            n_domains,
            profileSize,
            cfg["center_latlon"],
        ),
        "ignitions": {},
        "postproc": {
            str(outer_dom): OUTER_POSTPROC,
            str(inner_dom): INNER_POSTPROC,
            "shuttle": "incremental",
            "description": cfg["description"],
        },
        "fmda_geogrid_path": fmdaUTC.strftime(
            f"wksp_fmda/CONUS/%Y%m/fmda-CONUS-%Y%m%d-%H/fmda-CONUS-%Y%m%d-%H.geo"
        ),
    }

################################################################################
# Create dictionary to hold user data and define current datetime in Local and UTC #
cfg = {} # dictionary to hold user inputs
now = datetime.now() # current datetime value
utcNow = datetime.now(UTC) # current utc datetime value
if utcNow.hour > 0 and utcNow.hour < 12:
    currHour = 0
else:
    currHour = 12

################################################################################
# Get the grib source from the user #
while True:
    gribSources = {
        "GFSF": 28, "HRRR": 3, 
        "NAM218": 12, "NAM227": 5
    } # resolutions in km (roughly)
    gribSourceInput = input("Enter a grib source {} (HRRR): ".format(gribSources.keys()))
    if gribSourceInput == "":
        gribSourceInput = "HRRR"
    gribSourceInput = gribSourceInput.upper()
    if not gribSourceInput in gribSources.keys():
        print("Please enter a valid grib source.")
        continue
    gribRes = gribSources[gribSourceInput]
    print(f"You have selected {gribSourceInput} ({gribRes} km)")
    print()
    cfg["grib_source"] = gribSourceInput
    break

################################################################################
# Get domain configuration and resolution based on target resolution
while True:
    innerRes = input("Select innermost atmospheric resolution in meters (300 m): ")
    if innerRes == "":
        innerRes = 300
    try:
        innerRes = float(innerRes)
    except:
        print("Please select a valid resolution")
        continue
    break
    
ref_ratio = 3
n_domains = 1 + round(np.log(gribRes * 1000. / innerRes) / np.log(ref_ratio))
innerRes = gribRes * 1000. / ref_ratio**(n_domains - 1)
print(f"Domain configuration with {n_domains} domains and inner most resolution of {innerRes} m\n")

################################################################################
# Get the profile size input from the user #
profileSizes = {
    "S": 72, "M": 96, "L": 144, "XL": 192
}
while True:
    print("Select size:")
    print(make_profile_table(profileSizes, innerRes))
    profileSizeInput = input("Enter a valid profile size [S, M, L, XL] (M): ")
    if profileSizeInput == "":
        profileSizeInput = "M"
    profileSizeInput = profileSizeInput.upper()
    if not profileSizeInput in profileSizes:
        print("Please select a valid profile size")
        continue
    print(f"You have selected the {profileSizeInput} profile size")
    print()
    break

################################################################################
# Get the simulation start time from the user 
while True:
    utcnow = datetime.now(UTC).strftime("%Y-%m-%d_%H:%M:%S")
    startUTC = input('Current UTC time is {}. Enter start time in UTC with format "%Y-%m-%d_%H:%M:%S": '.format(utcnow))
    try:
        startUTC = datetime.strptime(startUTC, "%Y-%m-%d_%H:%M:%S").replace(tzinfo=UTC)
    except:
        print(f"Enter a valid UTC time, not {startUTC}")
        continue
    if abs(utcNow - startUTC) > timedelta(days=7):
        print("WARNING: this is for forecasting operations, you might need to consider modifying some flags.")
    break

print(startUTC.strftime("Simulation Start Time: %Y-%m-%d_%H:%M:%S"))
print()
cfg["start_utc"] = startUTC.strftime("%Y-%m-%d_%H:%M:%S")

################################################################################
# Get the end time and ignition time #
while True:
    timeInput = input("How many hours would you like to run the simulation for: ")
    try:
        hours = int(timeInput)
    except:
        print(f"Please enter a valid hours value, not {timeInput}")
        continue
    break

endUTC = startUTC + timedelta(hours=hours)
cfg["end_utc"] = endUTC.strftime("%Y-%m-%d_%H:%M:%S")
# Make the ingition time the halfway point between the start and end times
timeDiff = (endUTC-startUTC)/2
ignUTC = startUTC + timeDiff
cfg["time_utc"] = ignUTC.strftime("%Y-%m-%d_%H:%M:%S")
print(f"The end time of the simulation: {endUTC}")
print()

################################################################################
# Get the central latitude and longitude values
while True:
    latInput = input("Enter latitude of the center of the fire as a decimal: ")
    try:
        latInput = float(latInput)
    except:
        print(f"Enter a valid latitude value as a decimal, not {latInput}")
        continue
    break

while True:
    lonInput = input("Enter longitude of the center of the fire as a decimal: ")
    try:
        lonInput = float(lonInput)
    except:
        print(f"Enter a valid longitude value as a decimal, not {lonInput}")
        continue
    break
    
cfg["center_latlon"] = [latInput,lonInput]
cfg["truelats"] = [latInput,latInput]
cfg["stand_lon"] = lonInput
cfg["latlon"] = [latInput,lonInput]
print()

################################################################################
# Get the fire name
fireNameInput = input("Enter the name of the fire (Test): ")
if fireNameInput == "":
    fireNameInput = "Test"
print(f"You have entered {fireNameInput} as the fire name")
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
# Get the code for the fire
hexCode = fireNameInput.encode("utf-8").hex()
cfg["grid_code"] = datetime.now(UTC).strftime(f"%Y-%m-%d_%H-{fireNameInput.replace(" ", "_").upper()}")

################################################################################
# Construct the description and fmda path for the json #
cfg["description"] = startUTC.strftime(f"{fireNameInput.title()} {profileSizeInput} %Y-%m-%d %Hz")
cfg["fmda_geogrid_path"] = startUTC.strftime("wksp_fmda/CONUS/%Y%m/fmda-CONUS-%Y%m%d-%H.geo")

################################################################################
# Modify other flags #
if run_wrf.upper() == "Y":
    cfg["run_wrf"] = True
else:
    cfg["run_wrf"] = False

################################################################################
# Generate the json #
fireNameModified = fireNameInput.replace(" ", "_")
jobID = startUTC.strftime("{}_{}_%y%m%d_%Hz".format(fireNameModified, profileSizeInput))
destinationPath = osp.join("jobs", jobID + ".json")

jsonData = build_job_json(
    cfg,
    gribRes,
    n_domains,
    profileSizes[profileSizeInput],
    profileSizeInput,
    fireNameInput,
)

# Finalize the changes to the json
os.makedirs(os.path.dirname(destinationPath), exist_ok=True)
with open(destinationPath, "w") as jsonFile:
    json.dump(jsonData, jsonFile, indent=4)

print(f"JSON is ready here: {destinationPath}")
print(f"./forecast.sh {destinationPath} >& logs/{jobID}.log")