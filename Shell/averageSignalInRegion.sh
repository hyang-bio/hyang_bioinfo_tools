#!/usr/bin/env bash
pyPATH=~/scripts/python

bedF=${1}
bwF=${2}
n=${3}
outF=${4}

grep '^chr' ${bedF} | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3}' > ${outF}.bed.tmp2
python ${pyPATH}/getBigWigValue.py ${outF}_tmp ${outF}.bed.tmp2 ${n} ${bwF} # ${outF}_tmp_siteprof${n}[temp] # from Wen Wang()
gunzip -f ${outF}_tmp_siteprof1.gz
sort -S100G --parallel=24 -k1,1 -k2,2n ${outF}_tmp_siteprof1 > ${outF}
rm ${outF}_tmp_siteprof1 ${outF}.bed.tmp2