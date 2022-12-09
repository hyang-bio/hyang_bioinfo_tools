#!/usr/bin/env bash
anPATH=~/annotations
pyPATH=~/scripts/python

dirPATH=${1}
dataSource=${2} # GEO | EBI | DDBJ
dataId=${3}
dataType=${4}
sampleN=${5}
genomeVersion=${6}


# prepare
mkdir -p ${dirPATH}/${sampleN};cd ${dirPATH}/${sampleN}

step0_prepare(){
	python ${pyPATH}/download.py -d ${dataSource} -i ${dataId} -t ${dataType} -o ${sampleN}
	fastqc -t 2 ${sampleN}*.fq.gz &
	
	if [ -f ${sampleN}_part1_R1.fq.gz ]
	then
		cat ${sampleN}_part*_R1.fq.gz > ${sampleN}_R1.fq.gz
		if [ -f ${sampleN}_part1_R2.fq.gz ]
		then
			cat ${sampleN}_part*_R2.fq.gz > ${sampleN}_R2.fq.gz
		fi
	fi
}
step0_prepare


# mapping to genome by bowtie
step1_mappingUsingBowtie2(){
	if [ -f ${sampleN}_R2.fq.gz ]
	then
		hisat2 -p 30 --dta -x ${anPATH}/${genomeVersion}/HISAT2/${genomeVersion} -1 ${sampleN}_R1.fq.gz -2 ${sampleN}_R2.fq.gz | samtools sort - -o ${sampleN}.bam
	else
		hisat2 -p 30 --dta -x ${anPATH}/${genomeVersion}/HISAT2/${genomeVersion} -U ${sampleN}.fq.gz | samtools sort - -o ${sampleN}.bam
	fi
}
step1_mappingUsingBowtie2

step2_quantify(){
	stringtie ${sampleN}.bam -o ${sampleN}.gtf -p 30 -G ${anPATH}/${genomeVersion}/${genomeVersion}.refGene.withGeneSymbol.gtf -A ${sampleN}.gene_abund.tab -B -e
	mv t_data.ctab ${sampleN}.t_data.ctab
	mv e_data.ctab ${sampleN}.e_data.ctab
	mv i_data.ctab ${sampleN}.i_data.ctab
	mv e2t.ctab ${sampleN}.e2t.ctab
	mv i2t.ctab ${sampleN}.i2t.ctab
}
step2_quantify