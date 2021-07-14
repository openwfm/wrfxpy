import numpy as np
import struct
import re
import collections
import six

def iterable(arg):
    return (isinstance(arg, collections.Iterable) 
        and not isinstance(arg, six.string_types))

class WPSRecordError(Exception):
    pass

class WPSFormatError(Exception):
    pass

class WPSRecord(dict):
    def __init__(self, args):
        if not isinstance(args, dict):
            raise WPSRecordError('WPSRecord: need a dictionary object instead of {}.\nTry using: WPSRecord.from_file or GeoDriver.from_elements' % type(gdal_ds))
        super(WPSRecord, self).__init__(args)  
        self.ensure_args()
        field = self['field'].decode() if isinstance(self['field'],(bytes, bytearray)) else self['field']
        desc = self['desc'].decode() if isinstance(self['desc'],(bytes, bytearray)) else self['desc']
        units = self['units'].decode() if isinstance(self['units'],(bytes, bytearray)) else self['units']
        hdate = self['hdate'].decode() if isinstance(self['hdate'],(bytes, bytearray)) else self['hdate']
        self.repr = "{}:{}:{}:{}:level {:d} Pa:fcst time {:d}:from {}".format(
                        field,desc,units,_projection_map[self['iproj']][0],
                        int(self['xlvl']),int(self['xfcst']),
                        hdate.replace('_','').replace('-','').replace(':',''))

    def __str__(self):
        return self.repr

    def __repr__(self):
        return self.repr

    def ensure_args(self):
        for k in self.keys():
            if isinstance(self[k],str):
                self[k] = self[k].encode('ascii')

class WPSFormat(list):
    # Reference in https://www2.mmm.ucar.edu/wrf/OnLineTutorial/Basics/IM_files/IM_wps.php
    def __init__(self, record_list):
        if not isinstance(record_list, list) and sum([isinstance(elem,WPSRecord) for elem in record_list]) != len(record_list):
            raise WPSFormatError('WPSFormat: need to be a list of WPSRecord objects.\nTry using: WPSFormat.from_file or WPSFormat.from_params')
        super(WPSFormat, self).__init__(record_list)

    def __repr__(self):
        return '{} {}'.format(self.__class__, self.repr)

    @staticmethod
    def params_info():
        print(_params_info)

    @classmethod
    def from_params(cls, **params):
        _test_params(params)
        n_records = len(params["field"])
        record_list = []
        for r in range(n_records):
            record_params = {}
            for k,v in params.items():
                if iterable(v):
                    record_params.update({k: v[r]})
                else:
                    record_params.update({k: v})
            record_list.append(WPSRecord(record_params))
        wps_format = cls(record_list)
        wps_format.repr = "from_params {}".format(wps_format.__len__())
        return wps_format

    @classmethod
    def from_file(cls, wps_path):
        with open(wps_path, "rb") as f:
            wps_params = {}
            record_list = []
            raw_data = f.read()
            idx_current = 0
            idx_end_file = len(raw_data) - 1
            while True:
                iproj = None
                nx = None
                ny = None
                for idx_iter in range(5):
                    # Each record has the size (in SIZE_BYTES bytes) at the start and end of each record. This might be compiler 
                    # dependent though, so this might need to be modified. Also, the WPS files are stored big endian.
                    record_start = idx_current + SIZE_BYTES
                    record_size = struct.unpack(">i", raw_data[idx_current:record_start])
                    record_end = record_start + record_size[0]
                    record_data = raw_data[record_start:record_end]

                    fmt,fields = _parse_map[idx_iter](iproj)
                    parsed_record = _parse_data(record_data, fmt, fields, nx, ny)

                    iproj = parsed_record.get("iproj", iproj)
                    nx = parsed_record.get("nx", nx)
                    ny = parsed_record.get("ny", ny)

                    wps_params.update(parsed_record)
                    idx_current = record_end + SIZE_BYTES

                record_list.append(WPSRecord(wps_params)) 
                if record_end + SIZE_BYTES > idx_end_file:
                    break

        wps_format = cls(record_list)
        wps_format.wps_path = wps_path
        wps_format.repr = "{} {}".format(wps_format.wps_path,wps_format.__len__())
        return wps_format

    def to_file(self, path):
        with open(path, "wb") as f:
            for record in self:
                for idx_iter in range(5):
                    fmt,fields = _parse_map[idx_iter](record['iproj'])
                    vv = [record[field] for field in fields]
                    if idx_iter < 4:
                        size = None
                    else:
                        size = record['nx']*record['ny']
                    packed = _pack_data(fmt,vv,size=size)
                    np = struct.pack(">i",len(packed))
                    f.write(np + packed + np)

def _test_params(params):
    if 'iproj' not in params:
        raise WPSFormatError('WPSFormat.test_params - invalid or not existend parameter iproj')
    else:
        iproj = params['iproj']
        for v in _parse_map.values():
            param = v(iproj)[1]
            for p in param:
                if p not in params:
                    raise WPSFormatError('WPSFormat.test_params - invalid or not existend parameter {}'.format(p))

def _pack_data(fmt, vv, size=None):
    if size is not None:
        fmt = fmt.format(size)
        vv = list(vv[0].reshape(size,order="F"))
    else:
        fs = re.findall(_fmt_pattern,fmt)
        if len(fs) == len(vv):
            for i,f in enumerate(fs):
                if 's' in f:
                    num = f.split('s')[0]
                    if num == '':
                        vv[i] = vv[i].ljust(1)
                    else:
                        vv[i] = vv[i].ljust(int(num))
    return struct.pack(fmt, *vv)

def _parse_data(data, fmt, fields, nx, ny):
    r = {}
    if "slab" in fields:
        size = nx * ny
        d = struct.unpack(fmt.format(size), data)
        arr = np.array(d, dtype=np.float32)
        parsed = tuple([arr.reshape((nx, ny), order="F")])
    else:
        parsed = struct.unpack(fmt, data)
    if len(parsed) != len(fields):
        raise WPSFormatError('WPSFormat._parse_data - number of fields parsed {} does not match {} field(s)'.format(len(parsed),len(fields)))
    for i,field in enumerate(fields):
        p = parsed[i]
        if isinstance(p,(bytes, bytearray, str)):
            p = p.strip()
        r[field] = p
    return r

_parse_map = {
    0: lambda iproj: (">i", ["ifv"]),
    1: lambda iproj: (">24sf32s9s25s46sfiii", ["hdate","xfcst","map_source",
                                "field","units","desc","xlvl",
                                "nx","ny","iproj"]),
    2: lambda iproj: (_projection_map[iproj][1], _projection_map[iproj][2]),
    3: lambda iproj: (">i", ["is_wind_earth_rel"]),
    4: lambda iproj: (">{}f", ["slab"])
}
_projection_map = {
    0: ("Cylindrical equidistant projection",">8sfffff",
        ["startloc","startlat","startlon","deltalat","deltalon","earth_radius"]),
    1: ("Mercator projection",">8sffffff",
        ["startloc","startlat","startlon","dx","dy","truelat1","earth_radius"]),
    3: ("Lambert conformal projection",">8sffffffff",
        ["startloc","startlat","startlon","dx","dy","xlonc","truelat1","truelat2","earth_radius"]),
    4: ("Gaussian projection",">8sfffff",
        ["startloc","startlat","startlon","nlats","deltalon","earth_radius"]),
    5: ("Polar-stereographic projection",">8sfffffff",
        ["startloc","startlat","startlon","dx","dy","xlonc","truelat1","earth_radius"]),
}
_fmt_pattern = re.compile(r'[0-9]*[a-zA-Z]{1}')   
# Number of bytes used at the start and end of a fortran record to indicate the record size
SIZE_BYTES = 4
# Explanation of all the parameters
_params_info = ('VARNAME\t\t\tFORMAT\t\t\tDESCRIPTION\n'
    'ifv\t\t\tinteger\t\t\tintermediate-format version number for WPS is 5.\n'
    'hdate\t\t\tcharacter (LEN=24)\tThe time, in format "YYYY-MM-DD_HH:mm:ss" (only the first 19 characters are used).\n'
    'xfcst\t\t\treal\t\t\tForecast time (in hours) of the data in the slab.\n'
    'map_source\t\tcharacter (LEN=32)\tSource of data.\n'
    'field\t\t\tcharacter (LEN=9)\tA field name. Names with special meaning are described below.\n'
    'units\t\t\tcharacter (LEN=25)\tUnits describing the field in the slab.\n'
    'dec\t\t\tcharacter (LEN=46)\tText description of the field in the slab.\n'
    'xlvl\t\t\treal\t\t\tPressure-level (Pa) of the data. 200100 Pa indicates surface data; 201300 Pa indicates sea-level pressure.\n'
    'nx\t\t\tinteger\t\t\tSlab dimension in the X direction.\n'
    'ny\t\t\tinteger\t\t\tSlab dimension in the Y direction.\n'
    'iproj\t\t\tinteger\t\t\tFlag denoting the projection: 0-Lat/lon, 1-Mercator, 3-Lambert, 4-Gaussian, 5-Polar-stereographic.\n'
    'startloc\t\tcharacter (LEN=8)\tStart location of data. Could be "CENTER" or "SWCORNER". "SWCORNER" is typical.\n'
    'startlat\t\treal\t\t\tStarting latitude (degrees north).\n'
    'startlon\t\treal\t\t\tStarting longitude (degrees east).\n'
    'deltalat\t\treal\t\t\tLatitude increment (degrees) for lat/lon grid.\n'
    'deltalon\t\treal\t\t\tLongitude increment (degrees) for lat/lon grid.\n'
    'nlats\t\t\treal\t\t\tNumber of latitudes north of equator (for Gaussian grids).\n'
    'dx\t\t\treal\t\t\tGrid-spacing in x (km at TRUELAT1 (and TRUELAT2 as appropriate)).\n'
    'dy\t\t\treal\t\t\tGrid-spacing in y (km at TRUELAT1 (and TRUELAT2 as appropriate)).\n'
    'xlonc\t\t\treal\t\t\tCenter longitude of the projection.\n'
    'truelat1\t\treal\t\t\tExtra latitude (degrees north) used for defining Mercator, Polar Stereographic, and Lambert conformal projections.\n'
    'truelat2\t\treal\t\t\tA second extra latitude (degrees north) used for defining Lambert conformal projection.\n'
    'earth_radius\t\treal\t\t\tRadius of the earth.\n'
    'is_wind_earth_rel\tlogical\t\t\tLogical flag use to indicate if Lambert projected data has Earth or model rotated winds.\n'
    'slab\t\t\treal, dimension(NX,NY)\tTwo-dimensional array of data.')
