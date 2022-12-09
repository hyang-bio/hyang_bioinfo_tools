#!/usr/bin/env bash
# D: Mar-12-2021 16:08 Fri [Update]

anPATH=~/annotations
rPATH=~/scripts/R

dirPATH=${1}
genomeVersion=${2}
genomeF=${3}
bedF=${4}
regionL=(${5//;/ }) # 
outF=${6}


# overlap with promoter, genebody, intergenic, exon[5UTR, 3UTR, CDS], intron, lncRNA, DHS, LINE, SINE, LTR, CGI, 
cd ${dirPATH}
echo -e "Region\tOverlapedRegionCount\tOverlapedRegionRatio" > ${outF}.overlaped.txt
totalCount=`cat ${bedF} | wc -l`
for region in ${regionL[@]}
do
	# echo ${region}
	overlapedRegionCount=`intersectBed -u -a ${bedF} -b ${region} | wc -l`
	overlapedRegionRatio=`bc -l <<< ${overlapedRegionCount}/${totalCount} | awk '{printf "%.2f",$1}'`
	label=$(echo ${region#*${genomeVersion}.})
	label=$(echo ${label%.merged.bed*})
	echo -e "${label}\t${overlapedRegionCount}\t${overlapedRegionRatio}" >> ${outF}.overlaped.txt
done
Rscript ${rPATH}/barplot_args.r ${outF}.overlaped.txt "Overlaped Count in Different Regions" "${outF}.overlapedRegions.pdf" 6 4

# enrichment in promoter, genebody, intergenic, exon[5UTR, 3UTR, CDS], intron, lncRNA, DHS, LINE, SINE, LTR, CGI, 
genomeLength=`cat ${genomeF} | awk 'BEGIN{OFS="\t";FS="\t";TL=0}{TL+=$2}END{print TL}'`
echo -e "Region\tEnrichScore" > ${outF}.EnrichScore.txt
for region in ${regionL[@]}
do
	overlapLengthWithRegion=`intersectBed -a ${bedF} -b ${region} | awk 'BEGIN{OFS="\t";FS="\t";TL=0}{TL+=$3-$2}END{print TL}'`
	regionLength=`cat ${bedF} | awk 'BEGIN{OFS="\t";FS="\t";TL=0}{TL+=$3-$2}END{print TL}'`
	promoterLength=`cat ${region} | awk 'BEGIN{OFS="\t";FS="\t";TL=0}{TL+=$3-$2}END{print TL}'`
	overlapRatio=`bc -l <<< ${overlapLengthWithRegion}/${regionLength}`
	regionRatio=`bc -l <<< ${promoterLength}/${genomeLength}`
	enrichScore=`bc -l <<< ${overlapRatio}/${regionRatio} | awk '{printf "%.2f",log($1)/log(2)}'`
	label=$(echo ${region#*${genomeVersion}.})
	label=$(echo ${label%.merged.bed*})
	echo -e "${label}\t${enrichScore}" >> ${outF}.EnrichScore.txt
done
Rscript ${rPATH}/barplot_primary.r ${outF}.EnrichScore.txt "Enrichment Score in Different Regions" "${outF}.enrichment.pdf" -8 8 4 6 4