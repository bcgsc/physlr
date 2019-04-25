IFS=$'\n'
t=$(printf '\t')
files=$1
while read f;
do
	key=$(echo $(basename $f) | awk -F "." '{print $1}')
	step=$(cat $f | grep "Command" | grep -Eo '/bin/physlr\s[a-zA-Z]+-*[a-zA-Z]*|physlr-*[a-zA-Z]+-*[a-zA-Z]*')
	clocktime=$(cat $f | grep "Elapsed" | grep -Po '[0-9]*:*[0-9]+:[0-9]+.[0-9]+')
	memory=$(cat $f | grep "Maximum" | awk '{print $6}')
	echo "${step}${t}${clocktime}" >> ${key}.t
	echo "${step}${t}${memory}" >> ${key}.m
done < $files
