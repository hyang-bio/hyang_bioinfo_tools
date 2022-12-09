#!/usr/bin/env bash
anPATH=~/annotations
pyPATH=~/scripts/python

dirPATH=${1}
datasource=${2}
inputInfo=${3}
type=${4}
sampleN=${5}
genomeVersion=${6}


# prepare
mkdir -p ${dirPATH}/${sampleN};cd ${dirPATH}/${sampleN}

step0_prepare(){
	python ${pyPATH}/download.py -d ${datasource} -i ${inputInfo} -t ${type} -o ${sampleN}
	fastqc -t 2 ${sampleN}*.fq.gz &
}
step0_prepare


# mapping to genome by bowtie
step1_mappingUsingBowtie2(){
	if [ -f ${sampleN}_R2.fq.gz ]
	then
		bowtie2 -p 30 -x ${anPATH}/${genomeVersion}/Bowtie2/${genomeVersion} -1 ${sampleN}_R1.fq.gz -2 ${sampleN}_R2.fq.gz | samtools sort -@ 20 - -o ${sampleN}.bam
	else
		bowtie2 -p 30 -x ${anPATH}/${genomeVersion}/Bowtie2/${genomeVersion} -U ${sampleN}.fq.gz | samtools sort -@ 20 - -o ${sampleN}.bam
	fi
}
step1_mappingUsingBowtie2