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
