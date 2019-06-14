# write_geogrid.py
# Jan Mandel, June 2019

import os.path as osp
import numpy as np
import sys, os

def write_divide(file,divide='=',count=25):
    """
    Write a dividing line
    """
    f=open(file,'a')
    f.write(divide * count + '\n')
    f.close()
    
def write_table(file,lines,mode='w'):
    """
    write table with lines key = value
    """
    f=open(file,mode)
    for key,value in lines.iteritems():
        f.write("%s = %s\n" % (key,value))
    f.close()
    
def write_geogrid(path,array,index,bits=32):
    """
    Write geogrid dataset 
    :param path: the directory where the dataset is to be stored
    :param array: numpy array of real values, 2d or 3d
    :param index: json with geogrid index, with geolocation and description already set 
    :param bits: 16 or 32 (default)
    """

    if not osp.exists(path):
        os.makedirs(path)

    # write binary data file
    a=np.array(array)
    dims=a.shape
    if len(dims) < 3:
        dims = dims + (1,)
    xsize, ysize, zsize = dims
    scale = 2**(np.ceil(np.log2(np.max(np.abs(a))))-bits+1)
    a = np.round(a.flatten()/scale)
    if bits == 32:
        a = np.int32(a) 
    elif bits == 16:
        a = np.int16(a) 
    else:
        print 'unsupported word size'
        sys.exit(1) 
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
