#!/bin/bash
IFS=$'\n'

last=0
t=$(printf '\t')
step_t="Physlrstep"
mem_t="Memory"
memfile=$1

echo "${step_t}${t}${mem_t}" >> ${memfile}.tsv

while read line
do
        step=$(echo $line | awk -F $'\t' '{print $1}')
		mem=$(echo $line | awk -F $'\t' '{print $2}')        
	    echo "${step}${t}${mem}" >> ${memfile}.tsv
done < $memfile
