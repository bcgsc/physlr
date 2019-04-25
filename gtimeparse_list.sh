IFS=$'\n'
t=$(printf '\t')
#dir=$1 #this is the data directory that holds the output files
#key=$2 #this is the key identifier of the dataset (ie. f1chr2R or f1chr4)
#this script takes the output of a "gtime" command -- extracts the step name(need to be fixed), the process time in either hh:mm:ss or mm:ss (handles both), and the memory usage.
#it later on uses getMinutes.sh to convert this time format to readable minutes
#it writes both time and memory into separate files
#those files then need to be converted into tsv

#ls -d : to not include subdirectories, but simply the directory itself
#$dir* : to access the relative path and not just the file names --> because we need them in "cat $file" --> if we don't use the relative path, we look for the files in this script's directory & can't find them
#for the directory: you need to specify the last backward slash (ex: ~/PycharmProjects/physlrtiming/physlr/data/)

files=$1
while read f;
do
	key=$(echo $(basename $f) | awk -F "." '{print $1}')
	#step=$(cat $f | grep "Command" | grep -oP '../[a-zA-Z]+/physlr-* *[a-zA-Z]+-*[a-zA-Z]*-*[a-zA-Z]*' | awk '{ if(NF==2) {print $2;} else { split($0,a,"physlr-"); print a[2]}}')
	step=$(cat $f | grep "Command" | grep -Eo '/bin/physlr\s[a-zA-Z]+-*[a-zA-Z]*|physlr-*[a-zA-Z]+-*[a-zA-Z]*')
	clocktime=$(cat $f | grep "Elapsed" | grep -Po '[0-9]*:*[0-9]+:[0-9]+.[0-9]+')
	memory=$(cat $f | grep "Maximum" | awk '{print $6}')
	echo "${step}${t}${clocktime}" >> ${key}.t
	echo "${step}${t}${memory}" >> ${key}.m
done < $files
. ./getmin.sh ${key}.t
