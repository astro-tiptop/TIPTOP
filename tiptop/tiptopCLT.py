# Ok combinare sottosezioni dell'ini


# Note: this is for perf, so we do not consider interesting cases which have the same complexity
# i.e. different atmosphere with the same number of layers
# 3 - tipi di atmosfera: JQ1, 35 layer, 1 layer a caso, [ JQ4, 35 layer, 1 layer a caso ]
# 2 - sources_science: 1 position, 9 positions
# 3 - system: SCAO ( 1 HO, 1 LO) / LTAO (6 HO, 1 LO) / MCAO_ (6 HO, 3 LO)
# 2 - magnitudine: 100/1000

# tabella atmosfera: motiplicare per alpha le velocita'



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

my_parser.add_argument('outputFile',
                       metavar='outputFile',
                       type=str,
                       help='the output file name (full path)')

# Execute the parse_args() method
args = my_parser.parse_args()

InputPath, parametersFile = os.path.split( args.parametersFile )
baseOutpuPath, outputFile = os.path.split( args.outputFile )

overallSimulation( InputPath, parametersFile, baseOutpuPath, outputFile, doConvolve=True, pitchScaling = 0.9)

# test command: 
# ./tiptopCLT.py localTest/params1GPU localTest/test
