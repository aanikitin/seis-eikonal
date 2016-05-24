#!/bin/sh
# bin Nseq BIseq BJseq BKseq
bin_name="`basename ${2}`"
dir="${1}/${bin_name}"
mkdir -p --verbose ${dir}
for i in `seq ${3} ${4} ${5}`;
do
    for t in `seq ${6} ${7} ${8}`
    do
	for bi in `seq ${9} ${10} ${11}`
	do
	    for bj in `seq ${12} ${13} ${14}`
	    do
		for bk in `seq ${15} ${16} ${17}`
		do
		    for ntest in `seq ${18} ${19} ${20}`
		    do
			TNAME="${bin_name}-dim-${i}-thr-${t}-bsize-${bi}-${bj}-${bk}-ntest-${ntest}.csv"
			echo $TNAME
			exec 3>&2
			exec 2>> ${dir}/${TNAME}
			OMP_NUM_THREADS=$t ./${2} ${i} ${i} ${i} ${bi} ${bj} ${bk} >> ${dir}/${TNAME}
			exec 2>&3
		    done
		done
	    done
	done
    done
done
