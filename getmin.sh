#!/bin/bash
IFS=$'\n'
t=$(printf '\t')
#input is smth like: 0:04.97
#this script converts from hh:mm:ss or mm:ss format to only minutes(mm)
#awk -v with variables is used to be able to do math operations with double numbers (bc only works with ints)

timefile=$1
while read line; do
	step=$(echo $line | awk -F $'\t' '{print $1}')
	tme=$(echo $line | awk -F $'\t' '{print $2}')
	if [[ $(grep -o ':' <<< $(echo "${tme}") | grep -c .) -eq 1 ]]; then
		min=$(echo "${tme}" | awk -F ":" '{print $1}')
		sec=$(echo "${tme}" | awk -F ":" '{print $2}')
		write=$(awk -v m=${min} -v s=${sec} 'BEGIN {result=m+(s/60); print result}')
		echo "${step}${t}${write}" >> ${timefile}.min
	else
		hr=$(echo "${tme}" | awk -F ":" '{print $1}')
		min=$(echo "${tme}"| awk -F ":" '{print $2}')
		sec=$(echo "${tme}"| awk -F ":" '{print $3}')
		write=$(awk -v h=${hr} -v m=${min} -v s=${sec} 'BEGIN {result=(h*60)+m+(s/60); print result}')
		echo "${step}${t}${write}" >> ${timefile}.min
	fi
done < $timefile
