#!/usr/bin/env python
# Jan-24-2018

import sys
from twobitreader import TwoBitFile

USAGE = "USAGE: %prog <bedFile> <outputFile> <2bit>"



def getCpGInfo(sequence):
    """
    sequence Input, CpGInfo output
    """

    try:
        CpGRatio = 1.0 * len(sequence) * sequence.count('CG') / sequence.count('C') / sequence.count('G')
        CGContent = float(sequence.count('C')+sequence.count('G'))/float(len(sequence))
        ATContent = float(sequence.count('A')+sequence.count('T'))/float(len(sequence))
    except:
        CpGRatio = float("nan")
        CGContent = float("nan")
        ATContent = float("nan")

    return CpGRatio,CGContent,ATContent




genome = TwoBitFile(sys.argv[3])
with open(sys.argv[1]) as fhd, open(sys.argv[2], 'w') as rfhd:
    rfhd.write('#chrom\tstart\tend\tlen\tnum\tCpG_Ratio\tCGContent\tATContent\n')
    for line in fhd:
        line = line.strip().split()
        try:
            length = int(line[2])-int(line[1])
            number = genome[line[0]][int(line[1]):int(line[2])].upper().count("CG")
            CpG_Ratio,CGContent,ATContent = getCpGInfo(genome[line[0]][int(line[1]):int(line[2])].upper())
        except:
            length = 0
            number = 0
            CGContent = float('nan')
            ATContent = float('nan')
            CpG_Ratio = float("nan")
        rfhd.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(line[0],line[1],line[2],length,number,CpG_Ratio,CGContent,ATContent))
