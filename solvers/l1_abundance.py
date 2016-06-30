import csv
import sys
from gurobipy import *
from datetime import datetime

# matrix F
infile1 = sys.argv[1]
# b vector
infile2 = sys.argv[2]

logfile = sys.argv[3]
outfile = sys.argv[4]

print "Read F matrix: %s" %infile1
f1 = open(infile1, "r")

F1 = {}
name = []

# n is number of k-mers
s1 = datetime.now()
n = 0
for line1 in f1.readlines():
    line1 = line1.strip()
    #print(line)
    part1 = line1.split(",")
    #print part
    #if (part1[0].isdigit() != 0):
    if (part1[0] != "kmer"):
        for j in range(1,len(part1)):
        	F1[n,j-1] = float(part1[j])
        n = n + 1
    else:
        for idx in range(1,len(part1)):
            name.append(part1[idx])

# g is number of genomes
g = len(part1)-1

print("Number of K-mer = %d\n" %n)
print("Number of Genomes = %d\n" %g)
s2 = datetime.now()

print "Read b vector: %s" %infile2
f2 = open(infile2, "r")

b1 = []

for line2 in f2.readlines():
    line2 = line2.strip()
    #print(line)
    part2 = line2.split(",")
    #print part
    if (part2[0] != "kmer"):
        b1.append(float(part2[len(part2)-1]))   

F = F1
b = b1

"""
The l1-norm approximation problem is given by

minimize | Pu - q |_1 

with variable u in R^n and problem data P in R^(m x n) and q in R^m. The problem is equivalent to an LP

minimize (1)^T*v 
subject to  |P    -I| * (u, v) <= (q, -q)
            |-P  - I|   

"""

try:
	print "Building model......"
	s3 = datetime.now()
	m = Model("lpfroml1")

	# Variables
	x = [None] * (g)
	y = [None] * (n)
	for i in range(g):
            x[i] = m.addVar(lb = 0, vtype = GRB.CONTINUOUS, name=name[i])

        for j in range(n):
            y[j] = m.addVar(lb = 0, vtype = GRB.CONTINUOUS, name="y"+str(j))
	
	m.update()

	# Objective
	obj = LinExpr()	
	for i in range(n):
		obj += y[i]
	
	m.setObjective(obj, GRB.MINIMIZE)

        print len(b)
	# Constraints
        nc = 0
	for i in range(n):
            constr1 = LinExpr()
            constr2 = LinExpr()
            for j in range(g):
                constr1 += F[i,j] * x[j]
                constr2 += -1.0 * F[i,j] * x[j]
            nc += 2
            m.addConstr(constr1 - y[i] <= b[i])
            m.addConstr(constr2 - y[i] <= -1.0 * b[i])
		
	print "Finish building model"
	s4 = datetime.now()
	print s4 - s3

print "Number of constr = %d\n" %nc

	s5 = datetime.now()
        print "\n"
	m.optimize()
        m.write("file.lp")

        s6 = datetime.now()

        if m.SolCount > 0:
            out = open(outfile,'w')
            info = ""
            for v in m.getVars():    
                out.write('%s, %f \n' %(v.varName, float(v.x)))
            out.close()

        log = open(logfile,'w')
        log.write("Number of genomes = %d\n" %g)
        log.write("Number of k-mers = %d\n" %n)
        log.write("Read input = " + str(s2 - s1) + "\n")
        log.write("Build model = " + str(s4 - s3) + "\n")
        log.write("Solve LP = " + str(s4 - s3))
        log.close()

except GurobiError as e:
	print('Error reported')
        print e.errno
