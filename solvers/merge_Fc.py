import csv
import sys
from datetime import datetime
from numpy import array,linspace,zeros,eye,concatenate,sum as SUM,linalg

# folder for F
indir1 = sys.argv[1]
# folder for b 
indir2 = sys.argv[2]
# output folder
outdir = sys.argv[3]

def list_files(path):
    # returns a list of names (with extension, without full path) of all files 
    # in folder path
    files = []
    for name in os.listdir(path):
        if os.path.isfile(os.path.join(path, name)):
            files.append(str(os.path.join(path, name)))
    return files         

def sort_count(kmer, freq, count):
    for i in range(0, len(count)-1):
        for j in range(i+1, len(count)):
            if (count[i] < count[j]):
                temp = count[i]
                count[i] = count[j]
                count[j] = temp

                temp = freq[i]
                freq[i] = freq[j]
                freq[j] = temp

                temp = kmer[i]
                kmer[i] = kmer[j]
                kmer[j] = temp


files = list_files(indir1)

namex = []
x = []
p = 1
for f in files:
    fi = open(f, "r")
    filename = str(f).split("/")
    b = indir2 + filename[1]
    fc = open(b, "r")

    print(p, filename[1])
    F = []
    nameF = []
    for line in fi:
        line = line.strip()
        part = line.split(",")
        nameF.append(part[0])
        F.append(int(part[1]))

    b = []
    nameb = []
    for lineb in fc:
        lineb = lineb.strip()
        partb = lineb.split(",")
        if (partb[0] != "kmer"):
            nameb.append(partb[0])
            b.append(int(partb[1]))

    p = p + 1
    fi.close()
    fc.close()

    outfile = outdir + filename[1]
    out = open(outfile,"w")
    for i in range(0, len(F)):
        if (nameF[i] == nameb[i]):
            out.write("%s, %d, %d\n" %(nameF[i], F[i], b[i]))
        else:
            out.write("#,%s,%s\n" %(nameF[i], nameb[i]))
    out.close()
