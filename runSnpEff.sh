#!/usr/bin/env bash

# This is a naive script to run snpEff on multiple isolate genomes against the reference. It will format the output (so long as utility scripts
# are in the same directory) for analysis with the bioconductor package "GOstats" in R.
# Paths to different inputs need to be customised below (lines appended with '#customise' need to be customised).

EFFPATH=/home/mark/local/src/snpEff #customise
INPUTS=($(ls /home/mark/Rdrive/CCDM_Prog_10_Share-HANE0J-SE00128/Sclerotinia/pan_genome/outputs/Pawsey_02-17/*.var/*.highqual.vcf)) #customise
INTERPROPATH=/home/mark/Rdrive/CCDM_Prog_10_Share-HANE0J-SE00128/Sclerotinia/localBackup/Research/2016/03.16/out/interpro_scan/version_2.interpro.txt #customise

# Function to get commands for running snpEff in parallel on an array of inputs from the inputs directory.
function get_commands()

{

	array=("${INPUTS[@]}")
	effcommands=()

	for ((i=0; i<${#array[@]}; i++)); do
		stub=${array[$i]}
		stub=${stub##*/}
		stub=${stub%.*}
		effcommands+=("java -jar $EFFPATH/snpEff.jar eff Ssclerotiorum_v2 -stats $stub ${array[$i]} > $stub.vcf")
	done

	tabcommands=()

	for ((i=0; i<${#array[@]}; i++)); do
		stub=${array[$i]}
		stub=${stub##*/}
		stub=${stub%.*}
		tabcommands+=("perl onePerLine.pl $stub.vcf > $stub.opl.vcf")
	done

}

# Function to get interpro terms and format for GOstats analysis.
function interpro_formatting()

{

	stub=$1
	stub=${stub##*/}
	stub=${stub%.*}

	array=($( awk '{print $3}' $1 ))
	newarray=()

	for name in ${array[@]}; do
		newarray+=("$name	")
	done

	echo ${newarray[*]} | tr ' ' '\n' |  fgrep -f - $INTERPROPATH > $stub.interpro.txt
	perl InterProScanToGOstats.pl -i $stub.interpro.txt -p $stub.gos

}

# Get the commands.
get_commands

# Check if script has been run in current directory, if not do the appropriate process. First snpEff, then reformatting commands.
if ! ls *vcf > /dev/null 2>&1; then
	parallel --no-notice --jobs 6 ::: "${effcommands[@]}"
fi

if ! ls *opl.vcf > /dev/null 2>&1; then
	parallel --no-notice --jobs 6 ::: "${tabcommands[@]}"
fi

if ! ls *mf.txt > /dev/null 2>&1; then
	for file in ${INPUTS[@]}; do
		stub=$file
		stub=${stub##*/}
		stub=${stub%.*}
		grep "HIGH"	$stub.opl.vcf | awk 'BEGIN{OFS="\t"}{split($0,a,"|");print a[2],a[3],a[4]}' > $stub.HIGH.txt
		grep "MODERATE"	$stub.opl.vcf | awk 'BEGIN{OFS="\t"}{split($0,a,"|");print a[2],a[3],a[4]}' > $stub.MODERATE.txt
		grep "LOW"	$stub.opl.vcf | awk 'BEGIN{OFS="\t"}{split($0,a,"|");print a[2],a[3],a[4]}' > $stub.LOW.txt
	done
fi

# Now get interpro domains using previous function and reformat for GOstats.
if ! ls *interpro.txt >/dev/null 2>&1; then

	lowarray=($(ls *LOW.txt))
	modarray=($(ls *MODERATE.txt))
	higharray=($(ls *HIGH.txt))

	for file in ${lowarray[@]}; do
		interpro_formatting $file
	done
	for file in ${modarray[@]}; do
		interpro_formatting $file
	done
	for file in ${higharray[@]}; do
		interpro_formatting $file
	done
fi

