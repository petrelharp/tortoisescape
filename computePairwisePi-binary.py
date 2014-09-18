#!/usr/bin/env python
#PBS -q cmb
#PBS -l walltime=300:00:00
#PBS -l nodes=1:sl230s:ppn=1
#PBS -l mem=2000mb,pmem=2000mb,vmem=2000mb

import gzip
#import sys
#import os
import struct

pathIn = '/home/rcf-40/pralph/panfs/tortoises/first100scaffolds/first100scaffolds_pval1e-6_posteriorDump.geno.gz'
pathOut = '/home/rcf-40/pralph/panfs/tortoises/first100scaffolds/pairwisePi/first100scaffolds_pval1e-6_posteriorDump.pwp.gz'

def computePi(aR, aH, aA, bR, bH, bA):
    pi = aR*bA + 0.5*(aR*bH + aH*bR + aH*bA + aA*bH + aH*bH) + aA*bR
    return pi

nInd = 180

def generateHeader(nInd):
    header = []
    #header.append('{0}   {1}'.format('chrom', 'position'))
    x = 2
    for i in range(1, nInd+1):
        for j in range(x, nInd+1):
            header.append('{0}-{1}'.format(i, j))
        x += 1
    return header

# break binary file into chunks by line
def chunkBin(fin):
    chunk = 'x'
    while chunk:
        chunk=fin.read(3*nInd*8)
        if chunk: yield chunk

#each line contains 3*nInd doubles
genoStruct = '{0}{1}{2}'.format('=', 3*nInd, 'd') # looks like '=36d'

# set up output
fout = gzip.open(pathOut, 'w')
header = generateHeader(nInd)
print len(header)
fout.writelines('\t'.join(header) + '\n')

fin = gzip.open(pathIn,"rb")
            
for chunk in chunkBin(fin):
    line = struct.unpack(genoStruct, chunk)
    #print line
    rowPis = []
    #chrom = line[0]
    #pos = line[1]
    #rowPis.append('{0}   {1}'.format(str(line[0]), str(line[1])))
    x = 1
    for i in range(0, nInd):
        for j in range(x, nInd):
            aR = float(line[i*3])
            aH = float(line[i*3+1])
            aA = float(line[i*3+2])
            bR = float(line[j*3])
            bH = float(line[j*3+1])
            bA = float(line[j*3+2])
            #print aR, aH, aA, bR, bH, bA
            rowPis.append(str(computePi(aR, aH, aA, bR, bH, bA)))
        x += 1
    #print rowPis
    fout.writelines('\t'.join(rowPis) + '\n')

fout.close()            
            
