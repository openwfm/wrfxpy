"""
Constants and configuration values for NGFS incident processing.

This module centralizes all magic numbers, default values, and static
configurations used throughout the NGFS processing pipeline.
"""

# =============================================================================
# TIME CONSTANTS
# =============================================================================

# Default forecast length and timing
DEFAULT_FORECAST_HOURS = 24
IGNITION_LOOKBACK_HOURS = 2
IGNITION_LOOKAHEAD_HOURS = 4
DEFAULT_LOOKBACK_TIME_HOURS = 24

# Time adjustments for ignition detection
IGNITION_TIME_BUFFER_MINUTES = 60
EARLY_DETECTION_WINDOW_HOURS = 3.0
VIIRS_WAIT_HOURS = 4

# Job management
JOB_SLEEP_SECONDS = 150  # Pause between job starts to avoid metgrid crashes

# Data retrieval windows
DEFAULT_DAYS_TO_GET = 2
FIRMS_DEFAULT_DAYS = 3
PICKLE_LOOKBACK_DAYS = 7

# =============================================================================
# SPATIAL CONSTANTS
# =============================================================================

# Detection buffer and domain sizing
DETECTION_BUFFER_RADIUS_KM = 15
DETECTION_BUFFER_RADIUS_M = DETECTION_BUFFER_RADIUS_KM * 1000

# Domain size limits
MAX_DOMAIN_SIZE_DOUBLINGS = 5
MIN_DOMAIN_SIZE = 31

# Geographic thresholds
ALASKA_LATITUDE_THRESHOLD = 54.0
SWIR_MAX_DISTANCE_DEGREES = 0.04  # Max lat/lon diff between SWIR and nominal

# =============================================================================
# CLUSTERING PARAMETERS (DBSCAN)
# =============================================================================

# DBSCAN clustering for incident detection
DBSCAN_MIN_SAMPLES = 10  # Min points to form a cluster
DBSCAN_EPS_MULTIPLIER = 2  # Multiplier for average pixel resolution

# =============================================================================
# FIRE RADIATIVE POWER (FRP) THRESHOLDS
# =============================================================================

# FRP filtering for incident prioritization
FRP_CUTOFF_THRESHOLD = 5e3  # Filter out low-FRP incidents when many exist
FRP_CUTOFF_WHEN_FEW = 0  # No filtering when few incidents

# =============================================================================
# WRF/FIRE MODEL CONFIGURATION
# =============================================================================

# Default ignition duration
DEFAULT_IGNITION_DURATION_S = 240  # 4 minutes

# Time step adjustments for different models
BEHAVE_13_TIME_STEP = 6

# =============================================================================
# NWS WEATHER FORECAST OFFICE (WFO) CODES
# =============================================================================

# Default WFO list for unknown incident detection
DEFAULT_WFO_LIST = [
    "KBOU",  # Boulder, CO
    "KSLC",  # Salt Lake City, UT
    "KPDT",  # Pendleton, OR
    "KMFR",  # Medford, OR
    "KBOI",  # Boise, ID
    "KUNR",  # Rapid City, SD
    "KBYZ",  # Billings, MT
    "KGJT",  # Grand Junction, CO
    "KPUB",  # Pueblo, CO
]

# Alaska borough to WFO mapping
ALASKA_BOROUGH_TO_WFO = {
    'Anchorage Municipality': 'KAFC',
    'Fairbanks North Star Borough': 'KAFG',
    'Juneau City and Borough': 'KAJK',
    'Matanuska-Susitna Borough': 'KAFC',
    'Kenai Peninsula Borough': 'KAFC',
    'North Slope Borough': 'KAFG',
    'Nome Census Area': 'KAFG',
    'Ketchikan Gateway Borough': 'KAJK',
    'Yakutat City and Borough': 'KAJK',
    'Bethel Census Area': 'KAFG',
    'Dillingham Census Area': 'KAFC',
    'Kusilvak Census Area': 'KAFG',
    'Valdez-Cordova Census Area': 'KAFC',
    'Wrangell City and Borough': 'KAJK',
    'Haines Borough': 'KAJK',
    'Petersburg Borough': 'KAJK',
    'Sitka City and Borough': 'KAJK',
    'Prince of Wales-Hyder Census Area': 'KAJK',
    'Kodiak Island Borough': 'KAFC',
    'Bristol Bay Borough': 'KAFC',
    'Aleutians East Borough': 'KAFC',
    'Aleutians West Census Area': 'KAFC',
    'Denali Borough': 'KAFG',
    'Southeast Fairbanks Census Area': 'KAFG',
    'Yukon-Koyukuk Census Area': 'KAFG',
    'Lake and Peninsula Borough': 'KAFC',
    'Hoonah-Angoon Census Area': 'KAJK',
    'Skagway Municipality': 'KAJK',
    'Copper River Census Area': 'KAFC',
    'Northwest Arctic Borough': 'KAFG'
}

# =============================================================================
# GRIB DATA SOURCES BY REGION
# =============================================================================

# Special GRIB sources for specific regions
GRIB_SOURCE_ALASKA = 'NAM198'
GRIB_SOURCE_HAWAII = 'NAM196'

# GRIB sources requiring special cycle handling
CYCLE_BASED_GRIB_SOURCES = ['HRRR', 'HRRR_AK', 'NAM198', 'NAM196']
CYCLE_INTERVAL_HOURS = 6

# =============================================================================
# LANDFIRE / GEOSPATIAL DATA PATHS
# =============================================================================

# States using updated Landfire maps
UPDATED_LANDFIRE_STATES = ['CA', 'AZ', 'NV', 'UT', 'NM']

# Landfire configuration file paths
GEO_VARS_DEFAULT = 'etc/vtables/geo_vars.json'
GEO_VARS_2024 = 'etc/vtables/geo_vars.json_2024'
GEO_VARS_ALASKA = 'etc/vtables/geo_vars.json_alaska'
GEO_VARS_HAWAII = 'etc/vtables/geo_vars.json_hawaii'

# =============================================================================
# FMDA (FUEL MOISTURE DATA ASSIMILATION) PATHS
# =============================================================================

# FMDA base directories for different fire models
FMDA_BASE_CAWFE = '/data/WRFXPY/wksp_fmda/CONUS/'
FMDA_BASE_BEHAVE = '/data/jhaley/wrfxpy/wksp_fmda/CONUS/'

# States excluded from FMDA (non-CONUS)
NON_CONUS_STATES = ['Alaska', 'Hawaii']

# =============================================================================
# SATELLITE DATA CONFIGURATION
# =============================================================================

# GOES satellite sectors
GOES_18_SECTORS = ['CONUS', 'Full-Disk']
GOES_19_SECTORS = ['CONUS']

# Satellite pairs to check
DEFAULT_SAT_SECTOR_PAIRS = [
    (18, 'CONUS'),
    (19, 'CONUS'),
    (18, 'Full-Disk')
]

# FIRMS satellite names
FIRMS_SATELLITES = [
    'noaa_20',
    'noaa_21',
    'suomi',
    'landsat',
    'noaa_20_Alaska',
    'noaa_21_Alaska',
    'suomi_Alaska'
]

# =============================================================================
# CSV DATA PARSING
# =============================================================================

# Time column names for different CSV versions
TIME_COLS_V1 = ['incident_start_time', 'observation_time', 'initial_observation_time']
TIME_COLS_V2 = ['acq_date_time', 'pixel_date_time']

# Null column defaults for CSV reading
NULL_COLUMNS = {
    'incident_name': 'string',
    'incident_conf': 'string',
    'incident_type': 'string'
}

# =============================================================================
# FILE PATHS AND NAMING
# =============================================================================

# Default file paths
DEFAULT_POPULATION_FILE = 'ingest/NGFS/Population_by_US_County_July_2023.txt'
DEFAULT_NGFS_INGEST_DIR = 'ingest/NGFS'
DEFAULT_URBAN_GEOJSON_4326 = 'landfire/urban_contours_4326.geojson'
DEFAULT_URBAN_GEOJSON_5070 = 'landfire/urban_contours_5070.geojson'

# Configuration file paths
BASE_NGFS_CONFIG = 'jobs/base_ngfs_cfg.json'
NGFS_CONFIG_PATH = 'etc/ngfs.json'

# Output directories
JOBS_DIRECTORY = 'jobs/'
LOGS_DIRECTORY = 'logs/'
NGFS_OUTPUT_DIRECTORY = 'ngfs/'

# State pickle serialization.
# Measured on a 70.9 MB state pickle: gzip level 6 gives 7.1x for 1.65 s, while
# level 9 (the pandas default) costs 7.23 s for only 3% more. Decompression is
# 0.23 s either way and is on the critical path, since every run reads the most
# recent pickle back. xz reaches 12.1x but takes 14.6 s to write.
PICKLE_COMPRESSION = {'method': 'gzip', 'compresslevel': 6}
PICKLE_SUFFIX = '.pkl.gz'

# Read side accepts the historical uncompressed files alongside compressed ones,
# so no migration of existing state is required.
PICKLE_PATTERNS = ('*.pkl', '*.pkl.gz', '*.pkl.xz')

# =============================================================================
# INCIDENT NAMING AND FILTERING
# =============================================================================

# Characters to replace in incident names for file system compatibility
FILE_SAFE_REPLACE_CHARS = ['#', '(', ')', ':']
FILE_SAFE_REPLACEMENT = '_'

# Incident type filters
RX_BURN_INDICATOR = 'RX'  # Prescribed burn identifier

# =============================================================================
# PRIORITIZATION PARAMETERS
# =============================================================================

# Default number of incidents to auto-start
DEFAULT_NUM_STARTS = 25
UNLIMITED_STARTS = -1  # Value indicating no limit

# Population-based prioritization
PRIORITY_BY_POPULATION = False  # Default: prioritize by FRP instead

# =============================================================================
# MAP GENERATION
# =============================================================================

# Basemap projections and coordinates
# CONUS map boundaries
CONUS_LLCRNRLON = -119
CONUS_LLCRNRLAT = 22
CONUS_URCRNRLON = -64
CONUS_URCRNRLAT = 49
CONUS_LAT_1 = 33
CONUS_LAT_2 = 45
CONUS_LON_0 = -95

# Alaska map boundaries
ALASKA_LLCRNRLON = -164
ALASKA_LLCRNRLAT = 54
ALASKA_URCRNRLON = -130
ALASKA_URCRNRLAT = 73
ALASKA_LAT_1 = 63
ALASKA_LAT_2 = 68
ALASKA_LON_0 = -151

# Map marker sizes
INCIDENT_MARKER_SIZE_ONGOING = 15
INCIDENT_MARKER_SIZE_NEW = 30
INCIDENT_MARKER_SIZE_STARTED = 30

# Map image processing
MAP_IMAGE_SLEEP_SECONDS = 20  # Wait time for image file operations

# =============================================================================
# PROJECTION SYSTEMS
# =============================================================================

# EPSG codes for coordinate transformations
EPSG_WGS84 = 4326  # Standard lat/lon
EPSG_WEB_MERCATOR = 3857  # Web Mercator (meters)
EPSG_NAD83_ALBERS = 5070  # NAD83 Albers Equal Area

# =============================================================================
# DOMAIN SIZING PARAMETERS
# =============================================================================

# Processor and node scaling
PROCESSOR_DOUBLING_EXPONENT = 2  # How to scale processors with domain size
MAX_PROCESSORS = 400

# Subgrid ratio adjustments
SUBGRID_RATIO_REDUCTION_FACTOR = 0.5

# =============================================================================
# DATA TYPE HANDLING
# =============================================================================

# Fill values
LON_SWIR_FILL_VALUE = -999  # Fill value for missing SWIR lon data
MAX_VALID_LONGITUDE = 180.0  # Used to filter invalid coordinates

# =============================================================================
# PIXEL DETECTION
# =============================================================================

# Number of pixel corners
NUM_PIXEL_CORNERS = 4

# Lat/lon column naming patterns
LAT_CORNER_PATTERN = 'lat_c{}'  # Format: lat_c1, lat_c2, etc.
LON_CORNER_PATTERN = 'lon_c{}'
LAT_TC_CORNER_PATTERN = 'lat_tc_c{}'  # Terrain-corrected
LON_TC_CORNER_PATTERN = 'lon_tc_c{}'

# =============================================================================
# URL TEMPLATES
# =============================================================================

# NGFS CSV download URL template
NGFS_CSV_URL_TEMPLATE = (
    'https://bin.ssec.wisc.edu/pub/volcat/fire_csv/NGFS_daily/'
    'GOES-{goes_direction}/{sector}/{csv_filename}'
)

# GOES direction mapping
GOES_DIRECTION_MAP = {
    18: 'WEST',
    19: 'EAST'
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_lat_corner_columns():
    """Return list of latitude corner column names."""
    return [LAT_CORNER_PATTERN.format(i) for i in range(1, NUM_PIXEL_CORNERS + 1)]

def get_lon_corner_columns():
    """Return list of longitude corner column names."""
    return [LON_CORNER_PATTERN.format(i) for i in range(1, NUM_PIXEL_CORNERS + 1)]

def get_lat_tc_corner_columns():
    """Return list of terrain-corrected latitude corner column names."""
    return [LAT_TC_CORNER_PATTERN.format(i) for i in range(1, NUM_PIXEL_CORNERS + 1)]

def get_lon_tc_corner_columns():
    """Return list of terrain-corrected longitude corner column names."""
    return [LON_TC_CORNER_PATTERN.format(i) for i in range(1, NUM_PIXEL_CORNERS + 1)]

def is_cycle_based_grib_source(grib_source):
    """Check if a GRIB source requires cycle-based handling."""
    return grib_source in CYCLE_BASED_GRIB_SOURCES

def get_alaska_wfo_for_borough(borough_name):
    """Get WFO code for an Alaska borough."""
    return ALASKA_BOROUGH_TO_WFO.get(borough_name, 'Unknown')

def should_use_updated_landfire(state_abbrev):
    """Check if state uses updated Landfire maps."""
    return state_abbrev in UPDATED_LANDFIRE_STATES
