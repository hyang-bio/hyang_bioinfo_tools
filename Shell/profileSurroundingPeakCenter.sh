#!/usr/bin/env bash
shPATH=~/scripts/shell
rPATH=~/scripts/R


# usage: bash profileSurroudingregion.sh
# Function: generate average signal for profile surrouding regions
# Note: Check average signal to exclude those with specific high values

regionPATH=${1}
len_up=${2}
len_down=${3}
bin=${4}
bwPATH=${5}
label=${6}
outF=${7}

declare -a myarray


# label
myarray[${#myarray[*]}]=${label}

# extend
cat ${regionPATH} | awk -v Len_up=${len_up} -v Len_down=${len_down} 'function max(a,b){return a > b ? a:b}BEGIN{FS=OFS="\t"}{\
	print $1,max(int(($2+$3)/2)-Len_up, 0), int(($2+$3)/2)+Len_down}' > ${outF}.${label}.extend.bed # rm ${outF}.${label}.extend.bed

# region
bash ${shPATH}/averageSignalInRegion.sh ${outF}.${label}.extend.bed ${bwPATH} ${bin} ${outF}.${label}.extend.signal # rm ${outF}.${label}.extend.signal
myarray[${#myarray[*]}]=`cut -f 4- ${outF}.${label}.extend.signal | awk 'BEGIN{FS="\t";OFS="\t"}{for(i=1;i<=NF;++i){if($i!="NA"){a[i]+=$i;b[i]+=1;}}}END{for(x in a){print a[x]/b[x]}}'`

# write into file
echo ${myarray[*]} | tr " " "\t"  >> ${outF}

# rm
rm ${outF}.${label}.extend.signal
rm ${outF}.${label}.extend.bed