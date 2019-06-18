# write_geogrid.py
# Jan Mandel, June 2019

import os.path as osp
import numpy as np
import sys, os, logging, json
from utils import inq

def write_divide(file,divide='=',count=25):
    """
    Write a dividing line
    """
    f=open(file,'a')
    f.write(divide * count + '\n')
    f.close()
    
def write_table(file,lines_dict,mode='w',divider_char='=',divider_count=25,divider_before=False,divider_after=False):
    """
    write table 
    """
    f=open(file,mode)
    if divider_before:
        f.write(divider_char * divider_count + '\n')
    for key,value in lines_dict.iteritems():
        f.write("%s = %s\n" % (key,value))
    if divider_after:
        f.write(divider_char * divider_count + '\n')
    f.close()

def write_geogrid_var(path_dir,varname,array,description,index,bits=32):
    """
    write geogrid dataset and index 
    """
    path_dir=osp.abspath(path_dir)
    logging.info('write_geogrid_var path_dir=%s varname=%s array=%s index=%s' % (path_dir, varname, inq(array), str(index)))
    if not osp.exists(path_dir):
        os.makedirs(path_dir)

    # write geogrid dataset
    geogrid_ds_path = osp.join(path_dir,varname)
    index['description']=description
    write_geogrid(geogrid_ds_path,array,index,bits=bits)
 
    # write also the index as json entry to modify later
    index_json_path = osp.join(path_dir,'index.json')
    try:
        index_json = json.load(open(index_json_path,'r'))
    except:
        index_json = {}
    index_json[varname]=index
    json.dump(index_json,open(index_json_path,'w'), indent=4, separators=(',', ': ')) 

    geogrid_tbl_var = {'name':varname,
                   'dest_type':'continuous',
                   'interp_option':'average_4pt',
                   'abs_path':geogrid_ds_path,
                   'priority':1}

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

    
def write_geogrid(path,array,index,bits=32):
    """
    Write geogrid dataset 
    :param path: the directory where the dataset is to be stored
    :param array: numpy array of real values, 2d or 3d
    :param index: json with geogrid index, with geolocation and description already set 
    :param bits: 16 or 32 (default)
    """

    logging.info('write_geogrid_var path=%s array=%s index=%s' % (path, inq(array), str(index)))
    if not osp.exists(path):
        os.makedirs(path)
    # write binary data file
    a=np.array(array)
    dims=a.shape
    if len(dims) < 3:
        dims = dims + (1,)
    xsize, ysize, zsize = dims
    scale = 2**(np.ceil(np.log2(np.max(np.abs(a))))-bits+1)
    aa = np.round(a/scale)
    if bits == 32:
        aa = np.int32(aa) 
    elif bits == 16:
        aa = np.int16(aa) 
    else:
        print 'unsupported word size'
        sys.exit(1) 
    a = aa.transpose().flatten()
    data_file = "00001-%05i.00001-%05i" % (xsize, ysize)
    data_path = osp.join(path,data_file)
    a.tofile(data_path)
    
    # write index
    index.update({'type':'continuous',
                 'signed':'yes',
                 'scale_factor':scale,
                 'wordsize':bits/8,
                 'tile_x':xsize,
                 'tile_y':ysize,
                 'tile_z':zsize,
                 'endian':sys.byteorder})
    index_path = osp.join(path,'index')
    write_table(index_path,index)

def addquotes(s):
    """
    add quotes to string
    """
    return '"'+s+'"'

if __name__ == '__main__':
    test_name = 'test_geo'
    print 'testing write_geogrid to',test_name
    write_geogrid(test_name,[[1.0,2.0],[-1.0,-2.0]],{'name':addquotes(test_name)},bits=32)
