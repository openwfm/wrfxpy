# Helper module for Operational RNN class
# The core RNN functionality is in a distributable package ml_fmda
#     Building model from weights path, cyclical prediction with states,  
#     weight warping for transfer, scaling 3d arrays, ...
# This module is for functionality specific to wrfxpy

import copy
import warnings
import numpy as np

def read_twarp_file(path):
    import ast
    import re
    from pathlib import Path

    lines = Path(path).read_text().splitlines()
    line = next(x for x in lines if x.startswith("Twarp Params:"))

    params = line.removeprefix("Twarp Params: ")
    params = re.sub(r"np\.float64\(([^()]*)\)", r"\1", params)

    return ast.literal_eval(params)


