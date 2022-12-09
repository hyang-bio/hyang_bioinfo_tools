#!/usr/bin/env bash
pyPATH=~/scripts/python
shPATH=~/scripts/shell

bedF=${1}
samGbedF=${2}
nBin=${3}
outF=${4}
label=${5}


declare -a myarray
myarray[${#myarray[*]}]=${label}


grep '^chr' ${bedF} | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3}' > ${bedF}.tmp # rm

# divide into multiple bins
bedtools makewindows -b ${bedF}.tmp -n ${nBin} > ${bedF}.splited.tmp # rm

# obtain average for each bin
nBin_for=`bc <<< ${nBin}-1`
for n in `seq 1 ${nBin_for}` 0
do
	cat ${bedF}.splited.tmp | awk -v NB=${nBin} -v N=${n} 'BEGIN{OFS="\t";FS="\t"}{if(NR%NB==N){print $1,$2,$3}}' > ${bedF}.splited.tmp_${n}
	bash ${shPATH}/averageMethylInRegionMultipleThreads.sh ${bedF}.splited.tmp_${n} ${samGbedF} ${bedF}.splited.tmp_${n}.methyl 8
	myarray[${#myarray[*]}]=`cat ${bedF}.splited.tmp_${n}.methyl | awk 'BEGIN{OFS="\t";FS="\t";S=N=0;}{if($4!="NA"){S+=$4;N+=1;}}END{if(N==0){print "NA"}else{print S/N}}'`
done # n end

# write into file
echo ${myarray[*]} | tr " " "\t" >> ${outF}

# rm
rm ${bedF}.tmp ${bedF}.splited.tmp*