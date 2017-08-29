# Copyright (C) 2013-2016 Martin Vejmelka, CU Denver
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR
# A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import pygrib
import numpy as np


class GribFile(object):
    """
    Represents the content of one GRIB file.
    """

    def __init__(self, path, def_msg = 1):
        """
        Initializes metadata information.

        :param path: path to the file
        :param def_msg: the default message to process
        """
        self.grbf = pygrib.open(path)
        
    
    def __getitem__(self, msg_id):
        """
        Retrieve the a GRIB message by either index (if msg_id is int) or
        by name if msg_id is a string.
        
        :param msg_id: an int or a string
        :return: the message or raises ValueError if no such message
        """
        if type(msg_id) == str:
            msgs = self.grbf.select(name=msg_id)
            if len(msgs) == 0:
                raise ValueError('No message with name %s found.' % msg_id)
            return GribMessage(msgs[0])
        elif type(msg_id) == int:
            try:
                return GribMessage(self.grbf.message(msg_id))
            except:
                raise ValueError('GRIB file does not have message with index %d' % msg_id)
        
    
    def __iter__(self):
        """
        Implements iterator interface.
        
        :return: self
        """
        self.cur_msg = 0
        return self
        
        
    def next(self):
        """
        Implementation of iterator interface.
        
        :return: next object or raise StopIteration
        """
        self.cur_msg += 1
        try:
            return GribMessage(self.grbf.message(self.cur_msg))
        except:
            raise StopIteration


    def close(self):
        """
        Close the current file.
        """
        self.grbf.close()

        
        
class GribMessage(object):
    """
    Represents one GRIB message.
    """
    
    def __init__(self, grb):
        """
        Initialize a GRIB2 message.
        
        :param grb: The grib message from pygrib
        """
        self.grb = grb
        
        
    def name(self):
        """
        Retrieve the GRIB2 message name.
        
        :return: the name of the grib message
        """
        return self.grb.name
        

    def values(self):
        """
        Returns the values of observations in the message
        in the topology helpfully reconstructed by pygrib.
        
        :return: values in grid
        """
        return self.grb.values
        
        
    def latlons(self):
        """
        Return the latitudes and longitudes of the grid points of the field.
        
        :return: tuple latitudes in grid and longitudes in grid
        """
        return self.grb.latlons()
        
        
    def __str__(self):
        """
        Returns a string representation of the message as provided by pygrib.
        
        :return: descriptive string
        """
        return str(self.grb)



if __name__ == '__main__':
    import sys
    
    if len(sys.argv) < 2:
        print('usage: %s list <grib_file>' % sys.argv[0])
        print('usage: %s to_netcdf <grib_file> <message-name> <netcdf_file>' % sys.argv[0])
        sys.exit(1)
        
    if sys.argv[1] == 'list':
        gf = GribFile(sys.argv[2])
        print('Listing messages in GRIB file')
        for gm in gf:
            print(gm)
        gf.close()
        
    if sys.argv[1] == 'to_netcdf':
        import netCDF4
        gf = GribFile(sys.argv[2])
        msg_name = sys.argv[3]
        gm = gf[msg_name]
        vals = gm.values()
        lats, lons = gm.latlons()
        nc_path = sys.argv[4]
        d = netCDF4.Dataset(nc_path, 'w', format='NETCDF4')
        d.createDimension('south_north', vals.shape[0])
        d.createDimension('west_east', vals.shape[1])
        dlats = d.createVariable('XLAT', 'f4', ('south_north', 'west_east'))
        dlons = d.createVariable('XLONG', 'f4', ('south_north', 'west_east'))
        dvals = d.createVariable('VALUES', 'f4', ('south_north', 'west_east'))
        dlats[:,:] = lats
        dlons[:,:] = lons
        dvals[:,:] = vals
        d.close()
