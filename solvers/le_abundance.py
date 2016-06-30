import csv
import sys
import numpy
from gurobipy import *
from datetime import datetime
from numpy import array,linspace,zeros,eye,concatenate,sum as SUM,linalg

# directory contains Fb
indir1 = sys.argv[1]
# number of reference genomes
ref = sys.argv[2]
# result file
outfile = sys.argv[3]

def list_files(path):
    # returns a list of names (with extension, without full path) of all files 
    # in folder path
    files = []
    for name in os.listdir(path):
        if os.path.isfile(os.path.join(path, name)):
            files.append(str(os.path.join(path, name)))
    return files         


files = list_files(indir1)

namex = []
x = []
p = 1

total_F = numpy.zeros(shape=(1,ref))
total_c = numpy.zeros(shape=(1,ref))
for f in files:
    fi = open(f, "r")
    filename = str(f).split("/")
    
    print(p, filename[1])
    
    F = 0
    c = 0
    for line in fi:
        line = line.strip()
        part = line.split(",")
        
        F = F + int(part[1])
        c = c + int(part[2])

    namex.append(filename[1])
    total_F[0][p-1] = F
    total_c[0][p-1] = c
    
    p = p + 1
    fi.close()

#print(len(total_F))    
#print(len(total_c))
#print(len(namex))

out = open(outfile,"w")
for i in range(0, len(namex)):
	out.write("%s, %f\n" %(namex[i], (total_c[0][i]*1.0)/total_F[0][i]))
    
out.close()


