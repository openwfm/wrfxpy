# write_geogrid.py
# Jan Mandel, June 2019

from __future__ import absolute_import
from __future__ import print_function
import os.path as osp
import numpy as np
import sys, os, logging, json
from utils import inq, addquotes, Dict
from geo.var_wisdom import get_wisdom
from geo.geo_utils import fill_categories
import six
import pandas as pd

def write_divide(file,divide='=',count=31):
    """
    Write a dividing line
    """
    f=open(file,'a')
    f.write(divide * count + '\n')
    f.close()
    
def write_table(file,lines_dict,mode='w',divider_char='=',divider_count=31,divider_before=False,divider_after=False):
    """
    write table 
    """
    f=open(file,mode)
    if divider_before:
        f.write(divider_char * divider_count + '\n')
    for key,value in six.iteritems(lines_dict):
        f.write("%s = %s\n" % (key,value))
    if divider_after:
        f.write(divider_char * divider_count + '\n')
    f.close()

def wisdom_to_table(varname, wisdom):
    # initialize table
    table = {
        'name': wisdom.get('name', varname),
        'dest_type': wisdom.get('type','continuous'),
        'interp_option': wisdom.get('interp_option','default:average_gcell(4.0)+four_pt+average_4pt'),
        'priority': wisdom.get('priority',1)
    }
    # resolve path to data
    if 'abs_path' in wisdom:
        table.update({'abs_path': wisdom['abs_path']})
    elif 'rel_path' in wisdom:
        table.update({'rel_path': wisdom['rel_path']})
    # some adds to the table
    if 'fill_missing' in wisdom:
        table['fill_missing'] = wisdom['fill_missing']
    if 'smooth_option' in wisdom:
        table['smooth_option'] = wisdom['smooth_option']
    if 'subgrid' in wisdom:
        table['subgrid'] = wisdom['subgrid']
    if 'add_opts' in wisdom:
        for key in wisdom['add_opts'].keys():
            table[key] = wisdom['add_opts'][key]
    return table

def write_geogrid_var(path_dir,var,array,index,bits=32,coord=None):
    """
    write geogrid dataset and index 
    """
    path_dir=osp.abspath(path_dir)
    logging.info('write_geogrid_var path_dir=%s varname=%s array=%s index=%s' % (path_dir, var, inq(array), str(index)))
    if not osp.exists(path_dir):
        os.makedirs(path_dir)

    # get information from src/geo/var_wisdom.py
    wisdom = get_wisdom(var).copy()
    varname = wisdom.get('name', var)

    # write geogrid dataset
    geogrid_ds_path = osp.join(path_dir,varname)
    index['description'] = addquotes(wisdom.get('description',''))
    index['units'] = addquotes(wisdom.get('units',''))
    index['type'] = wisdom.get('type','continuous')
    index['signed'] = wisdom.get('signed','yes')
    bits = wisdom.get('bits',bits)
    scale = wisdom.get('scale',None)
    uscale = wisdom.get('unit_scale',None)

    # some adds to index
    if 'category_range' in wisdom:
        index['category_min'] = wisdom['category_range'][0]
        index['category_max'] = wisdom['category_range'][1]
    if 'missing_value' in wisdom:
        index['missing_value'] = wisdom['missing_value']
    if 'tile_bdr' in wisdom:
        index['tile_bdr'] = wisdom['tile_bdr']

    # categorical substitution and interpolation
    if index['type'] == 'categorical':
        fill = wisdom.get('fill',{})
        if isinstance(fill,str):
            fill_path = fill
            fill = {}
            if osp.exists(fill_path):
                try:
                    fill = Dict({})
                    df = pd.read_csv(
                        fill_path, names=['from','to'], index_col=False
                    )
                    cfrom = np.array(df.loc[:,'from'])
                    cto = np.array(df.loc[:,'to'])
                    # add filling of non-specified categories
                    if np.isnan(cto).sum():
                        near_inds = tuple(
                            cfrom[np.where(np.isnan(cto))].astype(int)
                        )
                        fill.update({near_inds: 'nearest'})
                    # add nearest neighbor interpolation
                    if np.isnan(cfrom).sum():
                        rest_val = int(cto[np.where(np.isnan(cfrom))[0][0]])
                        unique = np.unique(array)
                        rest_ind = np.array([
                            u for u in unique if u not in cfrom
                        ])
                        fill.update({tuple(rest_ind): rest_val})
                    # add all the rest
                    df = df.dropna().reset_index()
                    for k in range(len(df)):
                        fill.update({
                            int(df.loc[k,'from']): int(df.loc[k,'to'])
                        })
                except Exception as e:
                    logging.warning('write_geogrid_var fail reading fill CSV file {}'.format(fill_path))
                    logging.warning('with exception {}'.format(e))
            else:
                logging.warning('write_geogrid_var fail reading fill CSV file {}'.format(fill_path))
        array = fill_categories(array, fill, coord)

    write_geogrid(geogrid_ds_path,array,index,bits=bits,scale=scale,uscale=uscale)
 
    # write also the index as json entry to modify later
    index_json_path = osp.join(path_dir,'index.json')
    try:
        index_json = json.load(open(index_json_path,'r'))
    except:
        index_json = {}
    index_json[varname]=index
    json.dump(index_json,open(index_json_path,'w'), indent=4, separators=(',', ': ')) 

    # add abs_path to wisdom and create table from wisdom
    wisdom['abs_path'] = geogrid_ds_path
    geogrid_tbl_var = wisdom_to_table(varname, wisdom)

    # write a segment of GEOGRID.TBL
    geogrid_tbl_path=osp.join(path_dir,'GEOGRID.TBL')
    write_table(geogrid_tbl_path, geogrid_tbl_var, mode='a', divider_after=True)

    # write also as json 
    geogrid_tbl_json_path = osp.join(path_dir,'geogrid_tbl.json')
    try:
        geogrid_tbl = json.load(open(geogrid_tbl_json_path,'r'))
    except:
        geogrid_tbl = {}
    geogrid_tbl[varname]=geogrid_tbl_var 
    json.dump(index,open(geogrid_tbl_json_path,'w'), indent=4, separators=(',', ': ')) 
    
    json.dump(geogrid_tbl,open(geogrid_tbl_json_path,'w'), indent=4, separators=(',', ': ')) 


def write_geogrid(path,array,index,bits=32,scale=None,uscale=None):
    """
    Write geogrid dataset 
    :param path: the directory where the dataset is to be stored
    :param array: numpy array of real values, 2d or 3d
    :param index: json with geogrid index, with geolocation and description already set 
    :param bits: 16 or 32 (default)
    :param scale: numeric scale or None (default)
    :param uscale: numeric transform scale to change units or None (default)
    :param data_type: 'categorical' or 'continuous' (default)
    """

    logging.info('write_geogrid path=%s array=%s index=%s' % (path, inq(array), str(index)))
    if not osp.exists(path):
        os.makedirs(path)
    # write binary data file
    a = np.array(array)
    if uscale != None:
        a = a*float(uscale)
    dims = a.shape
    if len(dims) < 3:
        dims = dims + (1,)
        a = np.reshape(a,dims)
    xsize, ysize, zsize = dims
    if scale is None:
        scale = 2**(np.ceil(np.log2(np.max(np.abs(a))))-bits+1)
    if scale != 1.:
        a = np.round(a/scale)
    if bits == 32:
        a = np.int32(a) 
    elif bits == 16:
        a = np.int16(a) 
    else:
        print('unsupported word size')
        sys.exit(1) 
    a = a.transpose(2,0,1)
    logging.info('write_geogrid array min=%f max=%f avg=%f' % (a.min(), a.max(), a.mean()))
    zsize, ysize, xsize = a.shape
    data_file = "00001-%05i.00001-%05i" % (xsize, ysize)
    data_path = osp.join(path,data_file)
    a.flatten().tofile(data_path)
    # write index
    index.update({'scale_factor': scale,
                 'wordsize': bits // 8,
                 'tile_x': xsize,
                 'tile_y': ysize,
                 'tile_z': zsize,
                 'endian': sys.byteorder})
    index_path = osp.join(path,'index')
    write_table(index_path,index)
    
def read_table(path):
    """
    Read WPS geogrid-style index file into a dict.
    """
    out = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line:
                continue
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip().strip('"')
    return out

def _auto_cast(v):
    """
    Convert index strings to int/float when possible.
    """
    if isinstance(v, (int, float)):
        return v
    try:
        if any(c in str(v).lower() for c in [".", "e"]):
            return float(v)
        return int(v)
    except ValueError:
        return v

def read_geogrid(path, uscale=None, squeeze=True):
    """
    Reverse write_geogrid().

    Parameters
    ----------
    path : str
        Directory containing geogrid binary tile and index file.
    uscale : float or None
        If write_geogrid used `a = a * uscale` before writing, then this
        function will divide by uscale to recover the original units.
    squeeze : bool
        If True, return 2D array when tile_z == 1.

    Returns
    -------
    array : np.ndarray
        Reconstructed array, typically float64 after rescaling.
    index : dict
        Parsed index metadata.
    """
    index_path = osp.join(path, "index")
    index = read_table(index_path)
    index = {k: _auto_cast(v) for k, v in index.items()}
    xsize = int(index["tile_x"])
    ysize = int(index["tile_y"])
    zsize = int(index.get("tile_z", 1))
    wordsize = int(index["wordsize"])
    scale = float(index.get("scale_factor", 1.0))
    endian = index.get("endian", sys.byteorder)
    if wordsize == 4:
        dtype = np.int32
    elif wordsize == 2:
        dtype = np.int16
    elif wordsize == 1:
        dtype = np.int8
    else:
        raise ValueError(f"Unsupported wordsize: {wordsize}")
    if endian == "big":
        dtype = np.dtype(dtype).newbyteorder(">")
    elif endian == "little":
        dtype = np.dtype(dtype).newbyteorder("<")
    else:
        raise ValueError(f"Unsupported endian value: {endian}")
    data_file = f"00001-{xsize:05d}.00001-{ysize:05d}"
    data_path = osp.join(path, data_file)
    a = np.fromfile(data_path, dtype=dtype)
    expected = xsize * ysize * zsize
    if a.size != expected:
        raise ValueError(
            f"File size mismatch, expected {expected} values, got {a.size}"
        )
    # written shape was (zsize, ysize, xsize)
    a = a.reshape((zsize, ysize, xsize))
    # invert transpose(2,0,1) from the writer
    # saved: original (x,y,z) -> (z,x,y)? let's verify carefully below
    a = a.transpose(1, 2, 0)
    # undo integer scaling
    a = a.astype(np.float64) * scale
    # undo unit scaling, if requested
    if uscale is not None:
        a = a / float(uscale)
    if squeeze and a.shape[2] == 1:
        a = a[:, :, 0]
    return a, index

if __name__ == '__main__':
    test_name = 'test_geo'
    print('testing write_geogrid to',test_name)
    write_geogrid(test_name,[[1.0,2.0],[-1.0,-2.0]],{'name':addquotes(test_name)},bits=32)
