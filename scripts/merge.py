#! /usr/bin/env python

import string
from sys import argv, exit

## Usage is intended as ...
##
##  merge prefix npes saveFile

if len(argv) < 2:
  print "Usage: merge prefix npes output"
  exit(0);

p = int(argv[2]);
## first compute the total number of nodes.
p_temp=2

if len(argv)>3:
 fileOut=argv[3]+"_ws_"+str(p)+".csv"
else:
 fileOut=argv[1]+"_ws_"+str(p)+".csv"

fout=open(fileOut,'w');

while p_temp<=p:
  fname = argv[1]+"_"+str(p_temp)+".stat"
  fin = open (fname, "r")
  #print "reading filename "+fname
  for line in fin:
	fout.write(line);
  fin.close();
  p_temp=p_temp<<1
fout.close();
