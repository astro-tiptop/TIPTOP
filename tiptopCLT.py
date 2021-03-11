#!/usr/bin/env python

import os
import sys
import argparse
from tiptop import *

my_parser = argparse.ArgumentParser(description='TIPTOP Command Line Tool')

# Add the arguments

my_parser.add_argument('parametersFile',
                       metavar='parametersFile',
                       type=str,
                       help='the parameters file name (full path)')

my_parser.add_argument('windPsdFile',
                       metavar='windPsdFile',
                       type=str,
                       help='the windPsdFile file')

my_parser.add_argument('outputFile',
                       metavar='outputFile',
                       type=str,
                       help='the output file name (full path)')

# Execute the parse_args() method
args = my_parser.parse_args()

InputPath, parametersFile = os.path.split( args.parametersFile )
windPsdFile = args.windPsdFile
baseOutpuPath, outputFile = os.path.split( args.outputFile )

overallSimulation( InputPath, parametersFile, windPsdFile, baseOutpuPath, outputFile, doConvolve=True, pitchScaling = 0.9)

# test command: 
# ./tiptopCLT.py localTest/params1GPU MASTSEL/data/windpsd_mavis.fits localTest/test