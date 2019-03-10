#!/bin/bash
set -eu
IFS=$'\n'
cpp=$1
py=$2

cut -f1 $cpp | sort -h > barcodes.$cpp
cut -f1 $py | sort -h > barcodes.$py

diff=$(diff -q barcodes.$cpp barcodes.$py)
if [[ "$diff" != "" ]]; then
	echo "Barcodes are not identical!"
	exit 1
fi

while read bx; do
	grep "${bx}" $cpp | cut -f2 | tr " " "\n" | sort -n > ${bx}.${cpp}
	grep "${bx}" $py | cut -f2 | tr " " "\n" | sort -n > ${bx}.${py}
	diff=$(diff -q ${bx}.${cpp} ${bx}.${py})
	if [[ "$diff" != "" ]]; then
		echo "${bx} has different minimizers!"
		exit 1
	else
		rm ${bx}.${cpp} ${bx}.${py}
	fi
done < barcodes.${cpp}

rm barcodes.${cpp} barcodes.${py}
