#!/usr/bin/env bash
# usage: bash FRiP.sh input.bam peak.bed label outFile

jarPATH=~/scripts/bin/picard.jar
inputBam=${1}
peakBed=${2}
label=${3}
outFile=${4}


readsCount=`samtools view -@ 24 -c ${inputBam}`
if [ "${readsCount}" -ge "4000000" ]
then
	if [ ! -f sd4m_$(basename ${inputBam}) ]
	then
		# step1. sample down to 4M
		p=`bc -l <<< 4000000/${readsCount}`
		java -Xmx5000m -jar ${jarPATH} DownsampleSam I=${inputBam} O=sd4m_$(basename ${inputBam}) P=${p}
		
		# step2. bam to bed
		bamToBed -i sd4m_$(basename ${inputBam}) > sd4m.reads.bed
	fi

	
	# step3. fraction in peaks
	intersectBed -u -a sd4m.reads.bed -b ${peakBed} > sd4m.reads_peak.bed
	peaksCount=`cat sd4m.reads_peak.bed | wc -l`
	sd4mReadsCount=`cat sd4m.reads.bed | wc -l`
	FRiP=$(echo -e "${peaksCount}\t${sd4mReadsCount}" | awk 'BEGIN{FS="\t"}{printf "%.2f", $1/$2*100}')
	echo -e "${label}\t${FRiP}%" >> ${outFile}
fi
