#!/usr/bin/bash

# Description:
# 	Compares "*.mol.tsv" of the physlr's barcode to molecule stage, or any other algorithm for barcode to molecule separation of linked reads,
#			 with a reference file reporting the number of molecules per barcode.	
#	- this script only compares the reported number of molecules per barcode, and consistency just means equal number of molecules per barcode.
#	- undersplitting: less number of molecules per barcode compared to reference (golden truth ) of the molecule per barcodes.
#	- oversplitting: more number of molecules per barcode compared to reference (golden truth ) of the molecule per barcodes.
#
# Inputs: 
# 	$1: the name of the ".mol.tsv" query file in "../" directory (without ".mol.tsv") - or the name of the file for molecule counts in current directory (molCount_"NAME").
# 	$2: the reference (Tigmint results or results of another physlr run with the same barcodes - (molCount_"NAME")).
#	$3: If == "keep", keeps the files for comparison, otherwise removes the created files.
#	$4: If present (any string), the script skips the stats for 2+ mol/barcode as well as the stats for the whole data.
#
# Output:
# 	prints the stats of Consistent/underplsit/oversplit barcodes.
# 
# Written by Amirhossein Afshinfard

if [[ "$#" -eq 0 ]] || [[ "$#" -eq 2 ]] || [[ "$#" -eq 3 ]] || [[ "$#" -eq 4 ]]; then
	echo ""
else
	printf "\n Error: unxpected number of parameters: $#."
	printf "\n Expecting 1, 3 or 4."
	printf "\n Exiting.\n"
	exit
fi
if [ -n "$1" ];then
	toCheck=$1
	ref=$2
	keepfiles=$3	
	plus2=$4
else
	toCheck=f1_n50_st2 #n50.t75.st2e3  #test1.n50 #n50_st2
	#toCheck=GT_m4d100k./
	ref=GT_m8d100k
	plus2=""
fi

file1=1
file2=1
if [ ! -e molCount_${toCheck} ]; then
	if [ -e ../${toCheck}.mol.tsv ]; then
		printf "\n ##################################################\n ## Comparing ${toCheck}.mol.tsv with ref being molCount_${ref}:\n"		
		printf "\n Making molCount_${toCheck} from ${toCheck}.mol.tsv ... "
		file1=2
		header=$(cat ../${toCheck}.mol.tsv | awk '{if(gsub(/\t/,"")<2) print gsub(/\t/,"")}' | wc -l )
		headerLenght=$(($header-1))
		head ../${toCheck}.mol.tsv -n $headerLenght > ${toCheck}.header.mol.tsv
		awk -F\| '$1 { print substr($1,1,18)}' ${toCheck}.header.mol.tsv | uniq -c | sort -k 2 > molCount_${toCheck}
		rm ${toCheck}.header.mol.tsv
	else
		printf "\n Error in accessing the query; no file named ${toCheck}.mol.tsv or molCount_${toCheck} ."
		printf "\n You need to make ${toCheck}.mol.tsv using physlr's barcode to molecule stage, (or any other similar file)."
		printf "\n Exiting..."
		exit
	fi
else
	printf "\n ##################################################\n ## Comparing molCount_${toCheck} with ref being molCount_${ref}:\n"		
fi

if [ ! -e molCount_${ref} ]; then
	if [ -e ../${ref}.mol.tsv ]; then
		printf "\n Making molCount_${ref} from ${ref}.mol.tsv ... "
		file2=2
		header=$(cat ../${ref}.mol.tsv | awk '{if(gsub(/\t/,"")<2) print gsub(/\t/,"")}' | wc -l )
		headerLenght=$(($header-1))
		head ../${ref}.mol.tsv -n $headerLenght > ${ref}.header.mol.tsv
		awk -F\| '$1 { print substr($1,1,18)}' ${ref}.header.mol.tsv | uniq -c | sort -k 2 > molCount_${ref}
		rm ${toCheck}.header.mol.tsv
	else
		#cut -f4 fly/fly.f1.a0.65.d100k.m4.n5.q1.s2000.molecule.bed | uniq -c | sort -k 2  > compare/molCount_GT_m4d100k
		printf "\n Error in accessing the reference; no file named ${ref}.mol.tsv or molCount_${ref} ."
		printf "\n If you're comparing two physlr outputs you need to make ${ref}.mol.tsv for both using physlr,"
		printf "\n And if you're comparing your results with a golden truth like Tigmint, you need to make molCount_${ref} from your .bed file."
		printf "\n\t-Hint for making file from Tigmint's .bed:"
		printf "\n\t cut -f4 tigmintOutput.molecule.bed | uniq -c | sort -k 2  > molCount_nameHere"
		printf "\n Exiting..."
		exit
	fi
fi
awk 'NR==FNR{seen[$2]; next} $2 in seen' molCount_${toCheck} molCount_${ref} > molCount_${ref}_filt${toCheck} 
awk 'NR==FNR{seen[$2]; next} $2 in seen' molCount_${ref} molCount_${toCheck} > molCount_${toCheck}_filt${ref}
diff -y molCount_${toCheck}_filt${ref}  molCount_${ref}_filt${toCheck}   | grep "\-.*\-" > diff_${toCheck}_${ref} 

### All the barcodes 
printf "\n Number of barcodes for query (molCount_${toCheck}):"
num_barcodes=$(wc -l molCount_${toCheck} | awk {'print $1'})
printf " $num_barcodes"
printf "\n Number of barcodes for reference (molCount_${ref}):"
num_barcodes=$(wc -l molCount_${ref} | awk {'print $1'})
printf " $num_barcodes"
printf "\n\n - For all shared barcodes:"
printf "\n\n\tNumber of shared barcodes bewteen the two files:"
num_barcodes=$(wc -l diff_${toCheck}_${ref} | awk {'print $1'})
printf "\n\t $num_barcodes\t|  %% 100 \n"

printf "\n\tConsistent:"
num_cons=$(grep -v "|" diff_${toCheck}_${ref} | wc -l)
printf "\n\t $num_cons\t|  %% "
printf %.1f\\n "$((1000 * $num_cons/$num_barcodes))e-1"

printf "\n\tUndersplitting:"
num_under=$(grep "|" diff_${toCheck}_${ref} | awk {'print $1-$4'} | sort -g | uniq -c | awk {'if($2<0) print $1'} | awk '{s+=$1} END {print s}')
printf "\n\t $num_under\t|  %% "
printf %.1f\\n "$((1000 * $num_under/$num_barcodes))e-1"

printf "\n\tOversplitting:"
num_over=$(grep "|" diff_${toCheck}_${ref} | awk {'print $1-$4'} | sort -g | uniq -c | awk {'if($2>0) print $1'} | awk '{s+=$1} END {print s}')
printf "\n\t $num_over\t|  %% "
printf %.1f\\n "$((1000 * $num_over/$num_barcodes))e-1"

### Barcodes of more than 1 molecule reported by reference
if [ ! -n "$plus2" ]; then
	printf "\n - Only for barcodes of 2+ molecules reported by molCount_${ref}:"
	num_m_cons=$(grep -v "|" diff_${toCheck}_${ref} | awk {'if($3>1) print'} | wc -l)
	num_m_under=$(grep "|" diff_${toCheck}_${ref} | awk {'if($4>1) print'} | awk {'print $1-$4'} | sort -g | uniq -c | awk {'if($2<0) print $1'} | awk '{s+=$1} END {print s}')
	num_m_over=$(grep "|" diff_${toCheck}_${ref} | awk {'if($4>1) print'} | awk {'print $1-$4'} | sort -g | uniq -c | awk {'if($2>0) print $1'} | awk '{s+=$1} END {print s}')
	num_m_barcodes=$(($num_m_cons+$num_m_under+$num_m_over))
	printf "\n\n\tNumber of shared barcodes (with #mol > 1):"
	printf "\n\t $num_m_barcodes\t|  %% 100 \n"
	printf "\n\tConsistent (with #mol > 1):"
	printf "\n\t $num_m_cons\t|  %% "
	printf %.1f\\n "$((1000 * $num_m_cons/$num_m_barcodes))e-1"
	
	printf "\n\tUndersplitting (with #mol > 1):"
	printf "\n\t $num_m_under\t|  %% "
	printf %.1f\\n "$((1000 * $num_m_under/$num_m_barcodes))e-1"

	printf "\n\tOversplitting (with #mol > 1):"
	printf "\n\t $num_m_over\t|  %% "
	printf %.1f\\n "$((1000 * $num_m_over/$num_m_barcodes))e-1"
fi
printf "\n"
rm molCount_${ref}_filt${toCheck}
rm molCount_${toCheck}_filt${ref}

if [ keepfiles == "keep" ]; then
	printf "\n Removing created files"
	if [ file1 == 2 ]; then
		rm molCount_${toCheck}
		printf "\n Removed: molCount_${toCheck}"
	fi
	if [ file2 == 2 ]; then
		rm molCount_${ref}
		printf "\n Removed: molCount_${ref}"
	fi
	rm diff_${toCheck}_${ref}
	printf "\n Removed: diff_${toCheck}_${ref}"
fi
printf "\n Finished.\n"
