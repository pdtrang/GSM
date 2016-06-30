import scipy.optimize
import sys
import csv
import numpy
from numpy import array,linspace,zeros,eye,concatenate,sum as SUM,linalg
from datetime import datetime
from scipy.optimize import nnls

infile_F = sys.argv[1]
infile_b = sys.argv[2]
outfile = sys.argv[3]
#logfile = sys.argv[4]

f1 = open(infile_F, "r")
reader = csv.reader(f1)

# g is number of genomes
g = 0
li = f1.readline()
g = len(li.split(","))
#print g

# n is number of k-mers
n = 0
for row in reader:
    n += 1
n = n + 1
#print n
f1.close()

f = open(infile_F, "r")
s1 = datetime.now()

print "Reading matrix F: %s" %infile_F
F1 = numpy.zeros(shape=(n-1,g-1))
name = []

i = 0
for line in f.readlines():
    line = line.strip()
    #print(line)
    part = line.split(",")
    #print part 
    if (part[0] != "kmer"):
    	for j in range(1,len(part)):
    		F1[i][j-1] = float(part[j])
    	i = i + 1
    else:
        for idx in range(1,len(part)):
            name.append(part[idx])
    

s2 = datetime.now()

f2 = open(infile_b, "r")
print "Reading vector b: %s" %infile_b
b1 = []

for line2 in f2.readlines():
    line2 = line2.strip()
    #print(line)
    part2 = line2.split(",")
    #print part
    if (part2[0] != "K-mer"):
        b1.append(float(part2[len(part2)-1]))

F = F1
b = b1

sol = scipy.optimize.nnls(F, b)

out = open(outfile,'w')
for i in range(0, len(sol[0])):
	out.write("%s, %.30f\n" %(name[i], round(sol[0][i], 30)))
out.close()
