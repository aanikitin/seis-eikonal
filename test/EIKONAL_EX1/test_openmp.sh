#!/bin/sh
# bin Nseq BIseq BJseq BKseq
bin_name="`basename ${2}`"
dir="${1}/${bin_name}"
mkdir -p --verbose ${dir}
for i in `seq ${3} ${4} ${5}`;
do
    for t in `seq ${6} ${7} ${8}`
    do
	for ntest in `seq ${9} ${10} ${11}`
	do
	    TNAME="${bin_name}-dim-${i}-thr-${t}-ntest-${ntest}.csv"
	    echo ${TNAME}
	    exec 3>&2
	    exec 2>> ${dir}/${TNAME}
	    OMP_NUM_THREADS=${t} ./${2} ${i} ${i} ${i} >> ${dir}/${TNAME}
	    exec 2>&3
	done
    done
done
