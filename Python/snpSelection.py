#!/usr/bin/env python
# Modified from Chenfei Wang
import sys
import pysam
import twobitreader
import math
import subprocess
import pandas as pd
from scipy.stats import binom_test
from scipy.stats import fisher_exact

def split_bam_file(bamF,snpF,strain1,strain2):
	snp = {}
	for line in open(snpF,'r'):
		line = line.strip().split('\t')
		if line[0] not in snp:
				snp[line[0]] = {}
		snp[line[0]][int(line[1])] = [line[3],line[4]]

	votes = {}
	match, unmatch = 0,0
	infile = pysam.Samfile(bamF,'rb')
	for line in infile:
		if line.qname not in votes:
				votes[line.qname] = [0,0]
		# print(line.rname)
		chr = line.reference_name
		start = line.pos + 1
		if line.cigar:
			read = ''
			score = ''
			index = 0
			for cg in line.cigar:
				if cg[0] == 0:
					read += line.seq[index:(index+cg[1])]
					score += line.qual[index:(index+cg[1])]
					index += cg[1]
				elif cg[0] == 2:
					read += 'N'*cg[1]
					score += '!'*cg[1]
					index += cg[1]
				else:
					index += cg[1]
			for pos in range(0,len(read)):
				if chr in snp and start+pos in snp[chr] and ord(score[pos])-33 >= 30:
				# if chr in snp and start+pos in snp[chr]: # don't require quality
					if read[pos] in snp[chr][start+pos][0].split(';'): # A;T
						votes[line.qname][0] += 1
						match += 1
					elif read[pos] in snp[chr][start+pos][1].split(';'): # C;G
						votes[line.qname][1] += 1
						match += 1
					else:
						unmatch += 1

	outf1 = pysam.Samfile(bamF[:-4]+'_'+strain1+'.bam','wb',template=infile)
	outf2 = pysam.Samfile(bamF[:-4]+'_'+strain2+'.bam','wb',template=infile)
	outf3 = pysam.Samfile(bamF[:-4]+'_mixed.bam','wb',template=infile)
	outf4 = pysam.Samfile(bamF[:-4]+'_other.bam','wb',template=infile) # 
	maternal, paternal, mixed, other = 0, 0, 0, 0
	for line in pysam.Samfile(bamF,'rb'):
		if line.qname in votes:
			if sum(votes[line.qname]) > 0:
				if votes[line.qname][0]*1.0/sum(votes[line.qname]) >= 2.0/3:
					maternal += 1
					outf1.write(line)
				elif votes[line.qname][0]*1.0/sum(votes[line.qname]) <= 1.0/3:
					paternal += 1
					outf2.write(line)
				else:
					mixed += 1
					outf3.write(line)
		else:
			other += 1
			outf4.write(line)
	outf1.close()
	outf2.close()
	outf3.close()
	outf4.close()
	infile.close()

	print(match,unmatch)
	print(maternal, paternal, mixed, other)


def split_bam_file_methyl(bamF,snpF,strain1,strain2):
	snp = {}
	for line in open(snpF,'r'):
		line = line.strip().split('\t')
		if line[0] not in snp:
				snp[line[0]] = {}
		snp[line[0]][int(line[1])] = [line[3],line[4]]

	votes = {}
	match, unmatch = 0,0
	infile = pysam.Samfile(bamF,'rb')
	for line in infile:
		if line.qname not in votes:
				votes[line.qname] = [0,0]
		# print(line.rname)
		chr = line.reference_name
		start = line.pos + 1
		if line.cigar:
			read = ''
			score = ''
			index = 0
			for cg in line.cigar:
				if cg[0] == 0:
					read += line.seq[index:(index+cg[1])]
					score += line.qual[index:(index+cg[1])]
					index += cg[1]
				elif cg[0] == 2:
					read += 'N'*cg[1]
					score += '!'*cg[1]
					index += cg[1]
				else:
					index += cg[1]
			for pos in range(0,len(read)):
				if chr in snp and start+pos in snp[chr] and ord(score[pos])-33 >= 30:
					if read[pos] in snp[chr][start+pos][0].split(';'):
						votes[line.qname][0] += 1
						match += 1
					elif read[pos] in snp[chr][start+pos][1].split(';'):
						votes[line.qname][1] += 1
						match += 1
					else:
						unmatch += 1

	outf1 = pysam.Samfile(bamF[:-4]+'_'+strain1+'.bam','wb',template=infile)
	outf2 = pysam.Samfile(bamF[:-4]+'_'+strain2+'.bam','wb',template=infile)
	outf3 = pysam.Samfile(bamF[:-4]+'_mixed.bam','wb',template=infile)

	maternal, paternal, mixed = 0, 0, 0
	for line in pysam.Samfile(bamF,'rb'):
		if line.qname in votes:
			if sum(votes[line.qname]) > 0:
				if votes[line.qname][0]*1.0/sum(votes[line.qname]) >= 2.0/3:
					maternal += 1
					outf1.write(line)
				elif votes[line.qname][0]*1.0/sum(votes[line.qname]) <= 1.0/3:
					paternal += 1
					outf2.write(line)
				else:
					mixed += 1
					outf3.write(line)
	outf1.close()
	outf2.close()
	outf3.close()
	infile.close()

	print(match,unmatch)
	print(maternal, paternal, mixed)


def allelic_MPBN(infile, outfile, c1, c2, valid_count = 10):
	df = pd.read_table(infile, sep ='\t', header = None, skiprows = 0)
	df.columns = ['#Chr','Start','End','RefSeq','GeneSymbol','Strand','txStart','txEnd','Score','exonCount','exonStarts','exonEnds', '2cell.C57BL', '2cell.PWK', 'Morula.C57BL', 'Morula.PWK', 'ICM.C57BL', 'ICM.PWK', 'TE.C57BL', 'TE.PWK']
	outf = open(outfile, 'w')
	print('RefSeq\tGeneSymbol\tMat\tPat\tp-value\tAS\tclass', file = outf)

	with open(infile, 'r') as f:
		# next(f)
		# next(f)
		for line in f:
			line = line.strip().split('\t')
			if int(line[c1]) + int(line[c2]) < valid_count:
				# continue
				print('\t'.join([line[3],line[4],line[c1],line[c2],'-','-','N']), file=outf)
			else:
				# pBase = float(df.iloc[:,c1].sum())/(float(df.iloc[:,c1].sum()) + float(df.iloc[:,c2].sum()))
				# pvalue = binom_test(int(line[c1]),int(line[c1]) + int(line[c2]), pBase)+1E-300
				pvalue = binom_test(int(line[c1]),int(line[c1]) + int(line[c2]), 0.5)+1E-300
				if pvalue <= 0.001:
					if int(line[c1]) > int(line[c2]):
						print('\t'.join([line[3],line[4],line[c1],line[c2],str(pvalue),str(-math.log10(pvalue)),'M']), file=outf)
					else:
						print('\t'.join([line[3],line[4],line[c1],line[c2],str(pvalue),str(math.log10(pvalue)),'P']), file=outf)
				else:
					if int(line[c1]) > int(line[c2]):
						print('\t'.join([line[3],line[4],line[c1],line[c2],str(pvalue),str(-math.log10(pvalue)),'B']), file=outf)
					else:
						print('\t'.join([line[3],line[4],line[c1],line[c2],str(pvalue),str(math.log10(pvalue)),'B']), file=outf)
		outf.close()


def chip_allelic(infile, outfile, c1, c2, valid_count = 10):

	outf = open(outfile, 'w')
	print('#chr\tstart\tend\tC57BL\tPWK\tp-value\tAS\tclass', file=outf)
	for line in open(infile,'r'):
		line = line.strip().split('\t')		
		if int(line[c1]) + int(line[c2]) < valid_count:
			continue
			# print('\t'.join(line[:3]+[line[c1],line[c2],str(pvalue),str(-math.log10(pvalue)),'N']), file=outf)
		else:
			pvalue = binom_test(int(line[c1]),int(line[c1]) + int(line[c2]),0.5)+1E-300
			if pvalue <= 0.001:
				if int(line[c1]) > int(line[c2]):
					print('\t'.join(line[:3]+[line[c1],line[c2],str(pvalue),str(-math.log10(pvalue)),'M']), file=outf)
				else:
					print('\t'.join(line[:3]+[line[c1],line[c2],str(pvalue),str(math.log10(pvalue)),'P']), file=outf)
			else:
				if int(line[c1]) > int(line[c2]):
					print('\t'.join(line[:3]+[line[c1],line[c2],str(pvalue),str(-math.log10(pvalue)),'B']), file=outf)
				else:
					print('\t'.join(line[:3]+[line[c1],line[c2],str(pvalue),str(math.log10(pvalue)),'B']), file=outf)
	outf.close()


def chip_allelic_on_gene(infile, outfile, c1, c2):

	chip = {}	
	for line in open(infile, 'r'):
		line = line.strip().split('\t')
		chip['_'.join(line[:3])] = line[3:]

	promoter = {}
	outf = open(outfile, 'w')
	print('#chr\tstart\tend\trefseq\tsymbol\tstrand\tC57BL\tPWK\tp-value\tAS', file=outf)
	for line in open('/mnt/Storage3/Embryo/TEMP/mm9.promoter.bed','r'):
		line = line.strip().split('\t')
		m, p = 0, 0
		for i in range(int(line[1])/1000*1000, int(line[2])/1000*1000, 1000):
			if line[0]+'_'+str(i)+'_'+str(i+1000) in chip:
				m += int(chip[line[0]+'_'+str(i)+'_'+str(i+1000)][c1-3])
				p += int(chip[line[0]+'_'+str(i)+'_'+str(i+1000)][c2-3])
		pvalue = binom_test(m,m+p,0.5)+1E-300
		if m > p:
			print('\t'.join(line[:6]+[str(m),str(p),str(pvalue),str(-math.log10(pvalue))]), file=outf)
		else:
			print('\t'.join(line[:6]+[str(m),str(p),str(pvalue),str(math.log10(pvalue))]), file=outf)
	outf.close()

def fisher_single_site(rl1, cl1, rl2, cl2):

	m1,um1,m2,um2 = 0,0,0,0
	if isinstance(rl1,str):
		m1 = int(float(rl1)*float(cl1))
		um1 = int(float(cl1)-m1)
		m2 = int(float(rl2)*float(cl2))
		um2 = int(float(cl1)-m2)
	else:	
		for i in range(len(rl1)):
			if not rl1[i] == 'NA':
				m1 += int(float(rl1[i])*float(cl1[i]))
				um1 += int((1-float(rl1[i]))*float(cl1[i]))
			if not rl2[i] == 'NA':
				m2 += int(float(rl2[i])*float(cl2[i]))
				um2 += int((1-float(rl2[i]))*float(cl2[i]))

	score,pvalue = fisher_exact([[m1,um1],[m2,um2]])
	if m1*1.0/(m1+um1+1) > m2*1.0/(m2+um2+1):
		return [str(m1),str(um1),str(m2),str(um2),str(pvalue),str(-math.log10(pvalue+1E-300))]
	else:
		return [str(m1),str(um1),str(m2),str(um2),str(pvalue),str(math.log10(pvalue+1E-300))]

def methyl_allelic(value_file, CpG_file, c1, c2, outfile):
	outf = open(outfile,'w')
	methyl_bins = {}
	coverage_bins = {}

	for line in open(value_file,'r'):
		line = line.strip().split("\t")
		if not line[0] == 'chr':
			methyl_bins[line[0]+'_'+line[1]] = line[4:len(line):2]
			coverage_bins[line[0]+'_'+line[1]] = line[3:len(line):2]
		else:
			continue

	region = {}
	for line in open(CpG_file,'r'):
		line = line.strip().split('\t')
		if '@'.join(line[3:]) not in region:
			region['@'.join(line[3:])] = [line[0]+'_'+line[1]]
		else:
			region['@'.join(line[3:])].append(line[0]+'_'+line[1])

	for k in list(region.keys()):
		pos = set(region[k])
		value,coverage = [],[]
		output = []
		for i in pos:
			if i in methyl_bins:
				value.append(methyl_bins[i])
				coverage.append(coverage_bins[i])

		if value:
			value = list(zip(*value))
			coverage = list(zip(*coverage))
			output = fisher_single_site(value[c1], coverage[c1], value[c2], coverage[c2])

		if output:
			print(k.replace('@','\t')+'\t'+'\t'.join(output), file=outf)

	outf.close()


def chip_allelic_simple(readcnt, outname):
	# D: Jan-27-2022 21:17 Thu
	# N: readcnt format: #Chrom	Start	End	ReadCnt_mat	ReadCnt_pat
	# N: readcnt was counted using multiBamCov through multiBamCov snp_mat.bam snp_pat.bam -bed region.bed -D (including duplicate reads): need bai index
	# N: U: undetermined
	outf_detail = open(outname + '.details.txt', 'w')
	outf_basic = open(outname + '.basic.txt', 'w')
	print('#Chrom\tStart\tEnd\tMat\tPat\tP-value\tAS\tClass', file = outf_detail)
	print('#Chrom\tStart\tEnd\tClass', file = outf_basic)
	for line in open(readcnt,'r'):
		line = line.strip().split('\t')		
		pvalue = binom_test(int(line[3]),int(line[3]) + int(line[4]),0.5)+1E-300
		if pvalue <= 0.001:
			if int(line[3]) > int(line[4]):
				print('\t'.join(line[:3]+[line[3],line[4],str(pvalue),str(-math.log10(pvalue)),'M']), file=outf_detail)
				print('\t'.join(line[:3]+['M']), file=outf_basic)
			else:
				print('\t'.join(line[:3]+[line[3],line[4],str(pvalue),str(math.log10(pvalue)),'P']), file=outf_detail)
				print('\t'.join(line[:3]+['P']), file=outf_basic)
		else:
			if int(line[3]) > int(line[4]):
				print('\t'.join(line[:3]+[line[3],line[4],str(pvalue),str(-math.log10(pvalue)),'U']), file=outf_detail)
				print('\t'.join(line[:3]+['U']), file=outf_basic)
			else:
				print('\t'.join(line[:3]+[line[3],line[4],str(pvalue),str(math.log10(pvalue)),'U']), file=outf_detail)
				print('\t'.join(line[:3]+['U']), file=outf_basic)
	outf_detail.close()
	outf_basic.close()


def fisher_single_site_simple(rl1, cl1, rl2, cl2):

	m1,um1,m2,um2 = 0,0,0,0
	if isinstance(rl1,str):
		m1 = int(float(rl1)*float(cl1))
		um1 = int(float(cl1)-m1)
		m2 = int(float(rl2)*float(cl2))
		um2 = int(float(cl1)-m2)
	else:	
		for i in range(len(rl1)):
			if not rl1[i] == 'NA':
				m1 += int(float(rl1[i])*float(cl1[i]))
				um1 += int((1-float(rl1[i]))*float(cl1[i]))
			if not rl2[i] == 'NA':
				m2 += int(float(rl2[i])*float(cl2[i]))
				um2 += int((1-float(rl2[i]))*float(cl2[i]))

	score,pvalue = fisher_exact([[m1,um1],[m2,um2]])
	if pvalue <= 0.001:
		if m1*1.0/(m1+um1+1) > m2*1.0/(m2+um2+1):
			return [str(m1),str(um1),str(m2),str(um2),str(pvalue),str(-math.log10(pvalue+1E-300)), 'M']
		else:
			return [str(m1),str(um1),str(m2),str(um2),str(pvalue),str(math.log10(pvalue+1E-300)), 'P']
	else:
		if m1*1.0/(m1+um1+1) > m2*1.0/(m2+um2+1):
			return [str(m1),str(um1),str(m2),str(um2),str(pvalue),str(-math.log10(pvalue+1E-300)), 'U']
		else:
			return [str(m1),str(um1),str(m2),str(um2),str(pvalue),str(math.log10(pvalue+1E-300)), 'U']


def methyl_allelic_simple(value_file, region_CpG_file, outname):
	# D: Jan-28-2022 11:43 Fri
	# N: value_file format: #Chrom	Start	Mat_coverage	Mat_methratio	Pat_coverage	Pat_methratio
	# N: region_CpG_file format: #Chrom_CpG	Start_CpG	End_CpG	#Chrom_region	Start_region	End_region 
	# N: U: undetermined
	outf_detail = open(outname + '.details.txt','w')
	outf_basic = open(outname + '.basic.txt', 'w')
	print('#Chrom\tStart\tEnd\tMat_M\tMat_U\tPat_M\tPat_U\tP-value\tAS\tClass', file = outf_detail)
	print('#Chrom\tStart\tEnd\tClass', file = outf_basic)
	methyl_bins = {}
	coverage_bins = {}
	
	for line in open(value_file,'r'):
		line = line.strip().split("\t")
		if not line[0] == '#Chrom':
			methyl_bins[line[0]+'_'+line[1]] = line[3:len(line):2]
			coverage_bins[line[0]+'_'+line[1]] = line[2:len(line):2]
		else:
			continue
	
	region = {}
	for line in open(region_CpG_file,'r'):
		line = line.strip().split('\t')
		if '@'.join(line[3:]) not in region:
			region['@'.join(line[3:])] = [line[0]+'_'+line[1]]
		else:
			region['@'.join(line[3:])].append(line[0]+'_'+line[1])
	
	for k in list(region.keys()):
		pos = set(region[k])
		value,coverage = [],[]
		output = []
		for i in pos:
			if i in methyl_bins:
				value.append(methyl_bins[i])
				coverage.append(coverage_bins[i])
	
		if value:
			value = list(zip(*value))
			coverage = list(zip(*coverage))
			output = fisher_single_site_simple(value[0], coverage[0], value[1], coverage[1])
	
		if output:
			print(k.replace('@','\t')+'\t'+'\t'.join(output), file=outf_detail)
			print(k.replace('@','\t')+'\t'+output[6], file=outf_basic)
	
	outf_detail.close()
	outf_basic.close()


def methyl_allelic_count(value_file, region_CpG_file, outname):
	# D: Jan-29-2022 15:50 Sat
	# N: value_file format: #Chrom	Start	Mat_totalC	Mat_methC	Pat_totalC	Pat_methC
	# N: region_CpG_file format: #Chrom_CpG	Start_CpG	End_CpG	#Chrom_region	Start_region	End_region 
	outf_detail = open(outname + '.details.txt','w')
	print('#Chrom\tStart\tEnd\tMat_M\tMat_U\tPat_M\tPat_U', file = outf_detail)
	methyl_bins = {}
	coverage_bins = {}
	
	for line in open(value_file,'r'):
		line = line.strip().split("\t")
		if not line[0] == '#Chrom':
			methyl_bins[line[0]+'_'+line[1]] = line[3:len(line):2]
			coverage_bins[line[0]+'_'+line[1]] = line[2:len(line):2]
		else:
			continue
	
	region = {}
	for line in open(region_CpG_file,'r'):
		line = line.strip().split('\t')
		if '@'.join(line[3:]) not in region:
			region['@'.join(line[3:])] = [line[0]+'_'+line[1]]
		else:
			region['@'.join(line[3:])].append(line[0]+'_'+line[1])
	
	for k in list(region.keys()):
		pos = set(region[k])
		value,coverage = [],[]
		output = []
		for i in pos:
			if i in methyl_bins:
				value.append(methyl_bins[i])
				coverage.append(coverage_bins[i])
	
		if value:
			value = list(zip(*value))
			coverage = list(zip(*coverage))
			output = fisher_single_site_count(value[0], coverage[0], value[1], coverage[1])
	
		if output:
			print(k.replace('@','\t')+'\t'+'\t'.join(output), file=outf_detail)
	
	outf_detail.close()


def main():
	selection = sys.argv[1]
	if selection == "split":
		split_bam_file(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]) # bamF, snpF, strain1, strain2
	if selection == "splitMethyl":
		split_bam_file_methyl(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]) # bamF, snpF, strain1, strain2
	if selection == "alleleMPBN": # allele_negLogPValue
		allelic_MPBN(sys.argv[2],sys.argv[3], int(sys.argv[4]), int(sys.argv[5])) #infile, outfile, c1, c2, valid_count
	if selection == 'chip_allelic_simple':
		chip_allelic_simple(sys.argv[2], sys.argv[3]) # readcnt, outname
	if selection == 'methyl_allelic_simple':
		methyl_allelic_simple(sys.argv[2], sys.argv[3], sys.argv[4]) # value_file, region_CpG_file, outname
	if selection == 'methyl_allelic_count':
		methyl_allelic_count(sys.argv[2], sys.argv[3], sys.argv[4]) # value_file, region_CpG_file, outname

main()