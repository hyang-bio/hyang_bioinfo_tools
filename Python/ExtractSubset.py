#!/usr/bin/env python
import sys

"""
	extract content which contains seedF in column c[0 based] of fullF
	usage: python ExtractSubset.py seedF fullF c outF [ExtractSubset function]
	source: copy from gaolab
"""

def ExtractSubset(seedF,fullF,c,outF):
	seed = []
	seedF = open(seedF)
	for line in seedF:
		seed.append(line.strip())
	seedF.close()

	fullF = open(fullF)
	outF = open(outF,'w')
	for line in fullF:
		if line.strip().split()[c] in seed:
			outF.write(line)
	outF.close()
	fullF.close()

def ExtractNonSubset(seedF,fullF,c,outF):
	seed = []
	seedF = open(seedF)
	for line in seedF:
		seed.append(line.strip())
	seedF.close()

	fullF = open(fullF)
	outF = open(outF,'w')
	for line in fullF:
		if line.strip().split('\t')[c] not in seed:
			outF.write(line)
	outF.close()
	fullF.close()

def ExtractSubsetWithHeader(seedF,fullF,c,outF):
	seed = []
	seedF = open(seedF)
	for line in seedF:
		seed.append(line.strip())
	seedF.close()

	fullF = open(fullF)
	outF = open(outF,'w')
	n = 0
	for line in fullF:
		if n == 0:
			outF.write(line)
		if line.strip().split()[c] in seed:
			outF.write(line)
		n += 1
	outF.close()
	fullF.close()


def ExtractSubsetBed(seedF,fullF,outF):
	'''
	extract content according to bed location[fully overlap]
	usage: python ExtractSubsetBed.py seedF fullF outF
	'''
	seed = {}
	seedFH = open(seedF)
	for line in seedFH:
		ele = line.strip().split('\t')
		seed['@'.join(ele[0:3])] = ele[3]
		# print(line.strip())
	seedFH.close()

	fullF = open(fullF)
	outF = open(outF,'w')
	for line in fullF:
		if '@'.join(line.strip().split('\t')[0:3]) in seed:
			outF.write(line.strip()+'\t'+str(seed['@'.join(line.strip().split('\t')[0:3])])+'\n')
	outF.close()
	fullF.close()


def ExtractSubsetBedRegion(seedF,fullF,outF):
	seed = []
	seedFH = open(seedF)
	for line in seedFH:
		ele = line.strip().split('\t')
		seed.append('@'.join(ele[0:3]))
	seedFH.close()

	fullF = open(fullF)
	outF = open(outF,'w')
	for line in fullF:
		if '@'.join(line.strip().split('\t')[0:3]) in seed:
			outF.write(line)
	outF.close()
	fullF.close()

def main():
	if sys.argv[1] == "ExtractSubset":
		ExtractSubset(sys.argv[2],sys.argv[3],int(sys.argv[4]),sys.argv[5])
	if sys.argv[1] == "ExtractNonSubset":
		ExtractNonSubset(sys.argv[2],sys.argv[3],int(sys.argv[4]),sys.argv[5])
	if sys.argv[1] == "ExtractSubsetWithHeader":
		ExtractSubsetWithHeader(sys.argv[2],sys.argv[3],int(sys.argv[4]),sys.argv[5])
	if sys.argv[1] == "ExtractSubsetBed":
		ExtractSubsetBed(sys.argv[2],sys.argv[3],sys.argv[4])
	if sys.argv[1] == "ExtractSubsetBedRegion":
		ExtractSubsetBedRegion(sys.argv[2],sys.argv[3],sys.argv[4])

main()
