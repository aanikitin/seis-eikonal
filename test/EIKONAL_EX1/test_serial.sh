#!/bin/sh
# bin Nseq BIseq BJseq BKseq
bin_name="`basename ${2}`"
dir="${1}/${bin_name}"
mkdir -p --verbose ${dir}
for i in `seq ${3} ${4} ${5}`;
do
    for ntest in `seq ${6} ${7} ${8}`
    do
	TNAME="${bin_name}-dim-${i}-ntest-${ntest}.csv"
	echo ${TNAME}
	exec 3>&2
	exec 2>> ${dir}/${TNAME}
	OMP_NUM_THREADS=1 ./${2} ${i} ${i} ${i} >> ${dir}/${TNAME}
	exec 2>&3
    done
done
