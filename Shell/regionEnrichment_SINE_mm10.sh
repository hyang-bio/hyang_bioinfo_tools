anPATH=~/annotations
rPATH=~/scripts/R

dirPATH=${1}
bedF=${2}
outF=${3}

# overlap with SINE
cd ${dirPATH}

totalCount=`cat ${bedF} | wc -l`
for fml in 5S-Deu-L2 Alu B2 B4 ID MIR tRNA tRNA-Deu tRNA-RTE
do
	echo -e "Region\tOverlapedRegionCount\tOverlapedRegionRatio" > ${outF}.overlaped_SINE_${fml}.txt
	for sub in $(ls ${anPATH}/mm10/Repeats/SINE/${fml}/*.merged.bed)
	do
		overlapedRegionCount=$(intersectBed -u -a ${bedF} -b ${sub} | wc -l)
		overlapedRegionRatio=$(bc -l <<< ${overlapedRegionCount}/${totalCount} | awk '{printf "%.2f",$1}')
		region=$(echo ${sub#*mm10.})
		region=$(echo ${region%.merged.bed*})
		echo -e "${region}\t${overlapedRegionCount}\t${overlapedRegionRatio}" >> ${outF}.overlaped_SINE_${fml}.txt
	done # for sub end
	Rscript ${rPATH}/barplot_args.r ${outF}.overlaped_SINE_${fml}.txt "Overlaped Count in Different Regions" "${outF}.overlapedRegions_SINE_${fml}.pdf" 6 4
done # for fml end



# enrichment in SINE
genomeLength=`cat ${anPATH}/mm10/mm10_main.chrom.sizes | awk 'BEGIN{OFS="\t";FS="\t";TL=0}{TL+=$2}END{print TL}'`
for fml in 5S-Deu-L2 Alu B2 B4 ID MIR tRNA tRNA-Deu tRNA-RTE
do
	echo -e "Region\tEnrichScore" > ${outF}.EnrichScore_SINE_${fml}.txt
	for sub in $(ls ${anPATH}/mm10/Repeats/SINE/${fml}/*.merged.bed)
	do
		overlapLengthWithRegion=`intersectBed -a ${bedF} -b ${sub} | awk 'BEGIN{OFS="\t";FS="\t";TL=0}{TL+=$3-$2}END{print TL}'`
		regionLength=`cat ${bedF} | awk 'BEGIN{OFS="\t";FS="\t";TL=0}{TL+=$3-$2}END{print TL}'`
		promoterLength=`cat ${sub}| awk 'BEGIN{OFS="\t";FS="\t";TL=0}{TL+=$3-$2}END{print TL}'`
		overlapRatio=`bc -l <<< ${overlapLengthWithRegion}/${regionLength}`
		regionRatio=`bc -l <<< ${promoterLength}/${genomeLength}`
		enrichScore=`bc -l <<< ${overlapRatio}/${regionRatio} | awk '{printf "%.2f",log($1)/log(2)}'`
		region=$(echo ${sub#*mm10.})
		region=$(echo ${region%.merged.bed*})
		echo -e "${region}\t${enrichScore}" >> ${outF}.EnrichScore_SINE_${fml}.txt
	done # for sub end
	Rscript ${rPATH}/barplot_primary.r ${outF}.EnrichScore_SINE_${fml}.txt "Enrichment Score in Different Regions" "${outF}.enrichment_SINE_${fml}.pdf" -5 5 2.5 6 4
done # for fml end