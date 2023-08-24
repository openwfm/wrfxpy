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

def write_geogrid_var(path_dir,varname,array,index,bits=32,coord=None):
    """
    write geogrid dataset and index 
    """
    path_dir=osp.abspath(path_dir)
    logging.info('write_geogrid_var path_dir=%s varname=%s array=%s index=%s' % (path_dir, varname, inq(array), str(index)))
    if not osp.exists(path_dir):
        os.makedirs(path_dir)

    # get information from src/geo/var_wisdom.py
    wisdom = get_wisdom(varname).copy()

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
                    df = pd.read_csv(fill_path, names=['from','to'],index_col=False)
                    cfrom = np.array(df.loc[1:,'from'])
                    cto = np.array(df.loc[1:,'to'])
                    rest_val = df.loc[0,'from']
                    unique = np.unique(array)
                    rest_ind = np.array([u for u in unique if u not in cfrom])
                    fill = Dict({tuple(rest_ind): rest_val})
                    for k,key in enumerate(cfrom):
                        fill.update({key: cto[k]})
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

    geogrid_tbl_var = {'name': varname,
                   'dest_type': wisdom.get('type','continuous'),
                   'interp_option': wisdom.get('interp_option','default:average_gcell(4.0)+four_pt+average_4pt'),
                   'abs_path': geogrid_ds_path,
                   'priority': wisdom.get('priority',1)}

    # some adds to geogrid_tbl_var
    if 'fill_missing' in wisdom:
        geogrid_tbl_var['fill_missing'] = wisdom['fill_missing']
    if 'smooth_option' in wisdom:
        geogrid_tbl_var['smooth_option'] = wisdom['smooth_option']
    if 'subgrid' in wisdom:
        geogrid_tbl_var['subgrid'] = wisdom['subgrid']
    if 'add_opts' in wisdom:
        for key in wisdom['add_opts'].keys():
            geogrid_tbl_var[key] = wisdom['add_opts'][key]

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

if __name__ == '__main__':
    test_name = 'test_geo'
    print('testing write_geogrid to',test_name)
    write_geogrid(test_name,[[1.0,2.0],[-1.0,-2.0]],{'name':addquotes(test_name)},bits=32)
