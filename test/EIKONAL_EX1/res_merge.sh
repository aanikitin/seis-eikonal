#!/bin/bash
# merge results into csv files
printf "" > ${PWD##*/}.csv
flag=0
for dir in `find . -maxdepth 1 -mindepth 1 -type d`
do
  echo ${dir}
  FILES="${dir}/*"
  arr=($FILES)
  head -n -1 ${arr[0]} > ${dir}.csv
  if [ ${flag} -eq 0 ]; then
    head -n -1 ${arr[0]} >> ${PWD##*/}.csv
    flag=1
  fi
  for f in ${FILES}
  do
    # echo "Processing $f file..."
    # take action on each file. $f store current file name
    tail -n1 ${f} >> ${dir}.csv
    tail -n1 ${f} >> ${PWD##*/}.csv
  done
done
