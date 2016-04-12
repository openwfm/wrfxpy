# Copyright (C) 2013-2016 Martin Vejmelka, UC Denver
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

import gribapi
import numpy as np


class GribFile(object):
    """
    Represents the content of one GRIB file.
    """

    def __init__(self, path):
        """
        Initializes metadata information.

        :param path: path to the file
        """
        self.f = open(path)
        self.gid = gribapi.grib_new_from_file(self.f)

        # extract keys
        self.keyval = {}
        iterid = gribapi.grib_keys_iterator_new(self.gid,None)
        while gribapi.grib_keys_iterator_next(iterid):
            name = gribapi.grib_keys_iterator_get_name(iterid)
            val = gribapi.grib_get_string(iterid,name)
            self.keyval[name] = val
        gribapi.grib_keys_iterator_delete(iterid)

        # extract lat/lon/vals
        iterid = gribapi.grib_iterator_new(self.gid,0)

        self.lats = []
        self.lons = []
        self.values = []
        while 1:
            result = gribapi.grib_iterator_next(iterid)
            if not result: break
            [lat,lon,value] = result
            self.lats.append(lat)
            self.lons.append(lon)
            self.values.append(value)
                          
        gribapi.grib_iterator_delete(iterid)


    def close(self):
        """
        Close the current file.
        """
        gribapi.grib_release(self.gid)
        self.f.close()



if __name__ == '__main__':
    gf = GribFile('ds.terrainh.bin')
    print(np.max(gf.values))
    print(gf.keyval.keys())
    print(gf.keyval['gridType'])
    gf.close()

