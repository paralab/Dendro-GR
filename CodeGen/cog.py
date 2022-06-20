#!/usr/bin/python3
""" 
Cog code generation tool.
http://nedbatchelder.com/code/cog
Copyright 2004-2015, Ned Batchelder.

usage: cog.py inputfile     
    
"""

import time
start = time.perf_counter()

import sys
from cogapp import Cog

if(len(sys.argv)>1):
    ret = Cog().main(sys.argv)
    print("//Time to generate %s : %.2f sec" % (sys.argv[1], (time.perf_counter() - start)))
    sys.exit(ret)
