#!/bin/bash
IFS=$'\n'

last=0
t=$(printf '\t')
step_t="Physlrstep"
interval_t="Elapsedtime"
begin_t="Begin"
end_t="End"
timefile=$1

echo "${step_t}${t}${interval_t}${t}${begin_t}${t}${end_t}" >> ${timefile}.tsv

while read line
do
	step=$(echo $line | awk -F $'\t' '{print $1}')
	interval=$(echo $line | awk -F $'\t' '{print $2}')
	begin=$last
	last=$(awk -v b=${begin} -v i=${interval} 'BEGIN {result=b+i; print result}')
	end=$last
	echo "${step}${t}${interval}${t}${begin}${t}${end}" >> ${timefile}.tsv
done < $timefile
