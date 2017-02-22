#!/usr/bin/env python

import re
import optparse

parser	= optparse.OptionParser()
parser.add_option(	"-i",
			"--input",
			dest="input",
			metavar="INPUT FILE",
			help="Output file generated by getDomainInfoInterpro.py."	)

parser.add_option(	"-v",
			"--vcf",
			dest="vcf",
			metavar="VCF FILE",
			help="Vcf file output from snpEff."	)

parser.add_option(	"-o",
			"--outfile",
			dest="outfile",
			metavar="OUTPUT FILE",
			help="Output file to write to (defaults to 'out.txt'"		)

(options, args) = parser.parse_args()

if not options.input:
	parser.error("\nPlease provide input file.\n")

if not options.vcf:
	parser.error("\nPlease provide vcf file.\n") 

if not options.outfile:
	options.outfile	= "out.txt"

with open(options.input, "r") as f:
	lines	= f.readlines()

names = dict()
for line in lines:
	fields	= line.split("\t")
	name	= fields[0]
	start	= fields[1]
	end	= fields[2]
	dstart	= fields[3]
	dend	= fields[4]
	domain	= fields[5]
	desc	= fields[6]
	if name not in names:
		names[name]	= [start, end, dstart, dend, domain, desc]
	else:
		names[name]	= names[name] + [start, end, dstart, dend, domain, desc]

output	= open(options.outfile, "w")

with open(options.vcf, "r") as f:
	lines	= f.readlines()
	for line in lines:
		fields	= line.split("\t")
		if not re.match("^\t", line):
			chrom	= fields[0]
			pos	= fields[1]
		else:
			if fields[7] == "ANN":
				desc	= fields[8].split("|")
				if desc[3] in  names:
					if pos > names[desc[3]][2] and pos < names[desc[3]][3]:
						out	= [desc[3], desc[1], names[desc[3]][4], names[desc[3]][5]]
                                        	out	= "\t".join(out)
						output.write(out)
		

