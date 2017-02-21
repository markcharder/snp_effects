#!/usr/bin/env bash

# This is a naive script to run snpEff on multiple isolate genomes against the reference. It will format the output (so long as utility scripts
# are in the same directory) for analysis with the bioconductor package "GOstats" in R.
# Paths to different inputs need to be customised below (lines appended with '#customise' need to be customised).

EFFPATH=/home/mark/local/src/snpEff #customise
INPUTS=($(ls /home/mark/Rdrive/CCDM_Prog_10_Share-HANE0J-SE00128/Sclerotinia/pan_genome/outputs/Pawsey_02-17/*.var/*.highqual.vcf)) #customise
INTERPROPATH=/home/mark/Analyses/2017/02.17/filter_genes/interpro.filtered.txt #customise
OUTDIR=$1

if [ ! -d $1 ]; then 
	mkdir $1
fi

# Function to get commands for running snpEff in parallel on an array of inputs from the inputs directory.
function get_commands()

{

	array=("${INPUTS[@]}")
	effcommands=()

	for ((i=0; i<${#array[@]}; i++)); do
		stub=${array[$i]}
		stub=${stub##*/}
		stub=${stub%.*}
		effcommands+=("java -jar $EFFPATH/snpEff.jar eff Ssclerotiorum_v2 -stats $OUTDIR/$stub ${array[$i]} > $OUTDIR/$stub.vcf")
	done

	tabcommands=()

	for ((i=0; i<${#array[@]}; i++)); do
		stub=${array[$i]}
		stub=${stub##*/}
		stub=${stub%.*}
		tabcommands+=("perl onePerLine.pl $OUTDIR/$stub.vcf > $OUTDIR/$stub.opl.vcf")
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

	echo ${newarray[*]} | tr ' ' '\n' |  fgrep -f - $INTERPROPATH > $OUTDIR/$stub.interpro.txt
	perl InterProScanToGOstats.pl -i $OUTDIR/$stub.interpro.txt -p $OUTDIR/$stub.gos

}

# Get the commands.
get_commands

# Check if script has been run in current directory, if not do the appropriate process. First snpEff, then reformatting commands.
if ! ls $OUTDIR/*vcf > /dev/null 2>&1; then
	parallel --no-notice --jobs 6 ::: "${effcommands[@]}"
fi

if ! ls $OUTDIR/*opl.vcf > /dev/null 2>&1; then
	parallel --no-notice --jobs 6 ::: "${tabcommands[@]}"
fi

if ! ls $OUTDIR/*mf.txt > /dev/null 2>&1; then
	for file in ${INPUTS[@]}; do
		stub=$file
		stub=${stub##*/}
		stub=${stub%.*}
		grep "HIGH"	$OUTDIR/$stub.opl.vcf | awk 'BEGIN{OFS="\t"}{split($0,a,"|");print a[2],a[3],a[4]}' > $OUTDIR/$stub.HIGH.txt
		grep "MODERATE"	$OUTDIR/$stub.opl.vcf | awk 'BEGIN{OFS="\t"}{split($0,a,"|");print a[2],a[3],a[4]}' > $OUTDIR/$stub.MODERATE.txt
		grep "LOW"	$OUTDIR/$stub.opl.vcf | awk 'BEGIN{OFS="\t"}{split($0,a,"|");print a[2],a[3],a[4]}' > $OUTDIR/$stub.LOW.txt
	done
fi

# Now get interpro domains using previous function and reformat for GOstats.
if ! ls *interpro.txt >/dev/null 2>&1; then

	lowarray=($(ls $OUTDIR/*LOW.txt))
	modarray=($(ls $OUTDIR/*MODERATE.txt))

	for file in ${lowarray[@]}; do
		interpro_formatting $file
	done
	for file in ${modarray[@]}; do
		stub=$file
		stub=${stub##*/}
		stub=${stub%.*}
		stub=${stub%/MODERATE/HIGH}
		array=($(awk '{print $3}' $file))
		for name in ${array[@]}; do
			newarray+=("$name\n")
		done
		echo ${newarray[*]} | tr ' ' '\n' | fgrep -v -f - $OUTDIR/$stub.txt > $OUTDIR/$stub.temp
		mv $OUTDIR/$stub.temp $OUTDIR/$stub.txt
		interpro_formatting $file
	done
	higharray=($(ls $OUTDIR/*HIGH.txt))
	for file in ${higharray[@]}; do
		interpro_formatting $file
	done
fi

