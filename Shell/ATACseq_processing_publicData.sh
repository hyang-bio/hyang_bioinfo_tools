#!/usr/bin/env bash
# usage: bash ATACseq_processing_publicData.sh dataPATH dataSource dataId dataType sampleName genomeVersion
# Note: Refer to Wen Wang
anPATH=~/annotations
pyPATH=~/scripts/python

dirPATH=${1}
dataSource=${2}
dataId=${3}
dataType=${4}
sampleN=${5}
genomeVersion=${6}


# prepare
mkdir -p ${dirPATH}/${sampleN};cd ${dirPATH}/${sampleN}

step0_prepare(){
	python ${pyPATH}/download.py -d ${dataSource} -i ${dataId} -t ${dataType} -o ${sampleN}
	fastqc -t 2 ${sampleN}*.fq.gz &

	if [ -f ${sampleN}_R2.fq.gz ]
	then
		trim_galore --paired --fastqc --stringency 3 --length 30 ${sampleN}_R1.fq.gz ${sampleN}_R2.fq.gz --cores 4 > ${sampleN}.trim_galore.log 2>&1
	else
		trim_galore --fastqc --stringency 3 --length 30 ${sampleN}_R1.fq.gz --cores 4 > ${sampleN}.trim_galore.log 2>&1
	fi
}
step0_prepare


# mapping to genome by bowtie
step1_mappingUsingBowtie2(){
	if [ -f ${sampleN}_R2.fq.gz ]
	then
		bowtie2 -p 20 --trim-to 3:40 -x ${anPATH}/${genomeVersion}/Bowtie2/${genomeVersion} --no-mixed --no-unal -1 ${sampleN}_R1_val_1.fq.gz -2 ${sampleN}_R2_val_2.fq.gz | samtools view -@ 20 -f 0x2 -bSq 30 - -o ${sampleN}.bam 2> mapping.log
	else
		bowtie2 -p 20 --trim-to 3:40 -x ${anPATH}/${genomeVersion}/Bowtie2/${genomeVersion} --no-mixed --no-unal -U ${sampleN}_R1_trimmed.fq.gz | samtools view -@ 20 -bSq 30 - -o ${sampleN}.bam 2> mapping.log
	fi
}
step1_mappingUsingBowtie2


step2_filtering(){
	bamToBed -bedpe -i ${sampleN}.bam | awk 'BEGIN{FS=OFS="\t"}{if($9=="+" && $10=="-"){print $1, $2, $6, $7, $8, "."};if($9=="-" && $10=="+"){print $1, $5, $3, $7, $8, "."}}' | awk '$1 !~ /_/{if($3>$2 && $1!="chrM") print}' | sort -S 1% -k1,1 -k2,2n > ${sampleN}.fragments.bed &
	samtools view -@ 20 -f 0x40 ${sampleN}.bam | cut -f 3 | sort -S 1% | uniq -c | sort -k1,1rg -S 1% | awk 'BEGIN{FS=OFS="\t";print "#Chrom", "Number";}{print $2, $1;}' > ${sampleN}.chromosome_distribution.txt
	wait
}
step2_filtering


step3_pileup(){
	cut -f 1-3 ${sampleN}.fragments.bed | uniq > ${sampleN}.fragments_uniq.bed

	# Fragment length
	awk 'BEGIN{FS=OFS="\t"}{print $3-$2}' ${sampleN}.fragments_uniq.bed | sort -S 1% | uniq -c | sort -k2,2g -S 1% | awk 'BEGIN{FS=OFS="\t";print "FragmentLength", "Number";}{print $2, $1;}' > ${sampleN}.fragments_length.txt &

	# Generate pesudo single-end reads
	awk 'BEGIN{FS=OFS="\t";srand(1007)}{if(rand()<0.5){print $1, $2, $2+50;} else if($3-50<0){print $1, 0, $3} else{print $1, $3-50, $3}}' ${sampleN}.fragments_uniq.bed | sort -k1,1 -k2,2g -S 1% > ${sampleN}.fragments_uniq_SE_reads.bed
	n=$(cat ${sampleN}.fragments_uniq_SE_reads.bed | wc -l)
	c=$(bc -l <<< 1000000/${n})
	genomeCoverageBed -bga -scale ${c} -i ${sampleN}.fragments_uniq_SE_reads.bed -g ${anPATH}/${genomeVersion}/${genomeVersion}_main.chrom.sizes | grep -v '#' | sort -k1,1 -k2,2g -S 1% > ${sampleN}.fragments_uniq_SE_reads.bdg
	bedGraphToBigWig ${sampleN}.fragments_uniq_SE_reads.bdg ${anPATH}/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${sampleN}.fragments_uniq_SE_reads.bw

}
step3_pileup

step4_OCR(){
	awk 'BEGIN{FS=OFS="\t";srand(1007)}{if($3-$2<=100){\
		if(rand()<0.5){print $1, $2, $2+50, "r"NR, 0, "+";} \
		else if($3-50<0){print $1, 0, $3, "r"NR, 0, "-";} \
		else{print $1, $3-50, $3, "r"NR, 0, "-";}}}' ${sampleN}.fragments_uniq.bed | sort -k1,1 -k2,2g > ${sampleN}.fragments_uniq_OCR_SE_reads.bed

	n=$(cat ${sampleN}.fragments_uniq_OCR_SE_reads.bed | wc -l)
	c=$(bc -l <<< 1000000/${n})
	genomeCoverageBed -bga -scale ${c} -i ${sampleN}.fragments_uniq_OCR_SE_reads.bed -g ${anPATH}/${genomeVersion}/${genomeVersion}_main.chrom.sizes | grep -v '#' | sort -k1,1 -k2,2g -S 1% > ${sampleN}.fragments_uniq_OCR_SE_reads.bdg
	bedGraphToBigWig ${sampleN}.fragments_uniq_OCR_SE_reads.bdg ${anPATH}/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${sampleN}.fragments_uniq_OCR_SE_reads.bw
}
step4_OCR

step5_nucleosome(){
	awk 'BEGIN{FS=OFS="\t";srand(1007)}{if($3-$2>=180){\
		if(rand()<0.5){print $1, $2, $2+50, "r"NR, 0, "+";} \
		else if($3-50<0){print $1, 0, $3, "r"NR, 0, "-";} \
		else{print $1, $3-50, $3, "r"NR, 0, "-";}}}' ${sampleN}.fragments_uniq.bed | sort -k1,1 -k2,2g > ${sampleN}.fragments_uniq_nucleosome_SE_reads.bed

	n=$(cat ${sampleN}.fragments_uniq_nucleosome_SE_reads.bed | wc -l)
	c=$(bc -l <<< 1000000/${n})
	genomeCoverageBed -bga -scale ${c} -i ${sampleN}.fragments_uniq_nucleosome_SE_reads.bed -g ${anPATH}/${genomeVersion}/${genomeVersion}_main.chrom.sizes | grep -v '#' | sort -k1,1 -k2,2g -S 1% > ${sampleN}.fragments_uniq_nucleosome_SE_reads.bdg
	bedGraphToBigWig ${sampleN}.fragments_uniq_nucleosome_SE_reads.bdg ${anPATH}/${genomeVersion}/${genomeVersion}_main.chrom.sizes ${sampleN}.fragments_uniq_nucleosome_SE_reads.bw

}
step5_nucleosome