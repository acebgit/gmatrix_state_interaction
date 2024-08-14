
__author__ = 'Antonio Cebreiro'

import json
import sys 
import numpy as np
from numpy import linalg, sqrt
import pandas as pd
import math 
from scipy import constants
import matplotlib
matplotlib.use('TkAgg') # To be used in visualcode
import matplotlib.pyplot as plt

import timeit
from functools import partial

# Get the absolute path to local folder 'PyQChem'
import sys
import os
parent_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), 'C:\\Users\\HP.LAPTOP-F127N3L3\\Desktop\\my_programs\\PyQChem'))
if parent_directory not in sys.path: # Add 'folder2' to the Python path
    sys.path.append(parent_directory)
from pyqchem.parsers.parser_rasci import parser_rasci