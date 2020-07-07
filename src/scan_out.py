from __future__ import absolute_import
from __future__ import print_function
import sys
import re


def process_out(file):
    p = re.compile('Timing for')
    with open(file) as f:
        for line in f:
            if p.match(line):
                 print(line)
              
        

file = sys.argv[1]
process_out(file)
