#!/usr/bin/env python
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

import logging
import sys
import os
# Postprocessor must be imported first since it selects the Agg backend
from vis.postprocessor import Postprocessor
from vis.rasterizer import make_colorbar, basemap_raster_mercator, basemap_barbs_mercator
from vis.var_wisdom import convert_value, get_wisdom


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print(('usage: %s <wrfout_path> <var_instr> <prefix> [skip]' % sys.argv[0]))
        sys.exit(1)

    wrf_path = sys.argv[1]
    var_instr = None
    if sys.argv[2][0] == '@':
        var_instr = json.load(open(sys.argv[2][1:]))
    else:
        var_instr = {x:{} for x in sys.argv[2].split(',')}
        
    prefix = sys.argv[3]
    skip = 1
    if len(sys.argv) == 5:
        skip = int(sys.argv[4])

    p = Postprocessor(os.path.dirname(prefix), os.path.basename(prefix), var_instr)
    p.process_file(wrf_path, list(var_instr.keys()), skip)


