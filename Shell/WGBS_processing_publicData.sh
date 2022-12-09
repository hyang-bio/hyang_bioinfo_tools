#!/usr/bin/env bash

anPATH=~/annotations
pyPATH=~/scripts/python

dirPATH=${1}
datasource=${2}
inputInfo=${3}
type=${4}
sampleN=${5}
genomeVersion=${6}


mkdir -p ${dirPATH}/${sampleN};cd ${dirPATH}/${sampleN}


step0_prepare(){
	python ${pyPATH}/download.py -d ${datasource} -i ${inputInfo} -t ${type} -o ${sampleN}
	fastqc -t 20 ${sampleN}*.fq.gz &
}
step0_prepare


step1_mappingUsingBSMap_F(){
	if [ -f ${sampleN}_R2.fq.gz ]
	then
		trim_galore --phred33 --clip_R1 9 --clip_R2 9 --paired ${sampleN}_R1.fq.gz ${sampleN}_R2.fq.gz --cores 8 --fastqc > ${sampleN}.trim_galore.log 2>&1
		cat ${sampleN}_R1_val_1.fq.gz ${sampleN}_R2_val_2.fq.gz > ${sampleN}.fq.gz
	else
		trim_galore --phred33 --clip_R1 9 ${sampleN}.fq.gz --cores 8 --fastqc > ${sampleN}.trim_galore.log 2>&1
		mv ${sampleN}_trimmed.fq.gz ${sampleN}.fq.gz
	fi
	bsmap -a ${sampleN}.fq.gz -d ${anPATH}/${genomeVersion}/${genomeVersion}.fa -s 16 -v 0.1 -n 1 -R -o ${sampleN}.sam -p 20 >> mapping.sample.log 2>&1
	bsmap -a ${sampleN}.fq.gz -d ${anPATH}/Enterobacteria_phage_lambda.fa -s 16 -v 0.1 -n 1 -R -o ${sampleN}.lambda.sam -p 20 >> mapping.lambda.log 2>&1
}
step1_mappingUsingBSMap_F


step2_calculateMethylationLevel_F(){
	# mcall 
	mcall -m ${sampleN}.sam -p 20 -r ${anPATH}/${genomeVersion}/${genomeVersion}.fa --outputDir ./ >> mapping.sample.log 2>&1
	
	# methratio
	python ~/scripts/bin/methratio.py -o ${sampleN}.lambda.txt -d ${anPATH}/Enterobacteria_phage_lambda.fa -z ${sampleN}.lambda.sam >> mapping.lambda.log 2>&1
	
	# generate bw files
	grep -v '#' ${sampleN}.sam.G.bed | awk 'BEGIN{FS=OFS="\t"}{if($5>=3){print $1,$2,$3,$4}}' | sort -k1,1 -k2,2n > ${sampleN}.sam.G.bed.tmp
	bedGraphToBigWig ${sampleN}.sam.G.bed.tmp ${anPATH}/${genomeVersion}/${genomeVersion}.chrom.sizes ${sampleN}.bw
}
step2_calculateMethylationLevel_F


step3_removeTemporaryFiles_F(){
	rm ${sampleN}*.fq.gz
	rm ${sampleN}.sam.G.bed.tmp
	samtools view -@ 20 -bht ${anPATH}/${genomeVersion}/${genomeVersion}.chrom.sizes ${sampleN}.sam > ${sampleN}.bam && rm ${sampleN}.sam
}
step3_removeTemporaryFiles_F