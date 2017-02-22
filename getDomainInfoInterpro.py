#!/usr/bin/env python

# This script parses the output from InterProScan and retrieves information on positions of domains relative to gene coordinates 
# in a gff3 file.

import optparse
import sys
import re

parser = optparse.OptionParser()

parser.add_option(	"-i",
			"--input",
			dest="interpro_file",
			help="Interpro output containing Pfam domains.",
			metavar="INTERPRO OUTPUT"	)

parser.add_option(	"-g",
			"--gff3",
			dest="gff3",
			help="Gff3 file containing gene predictions.",
			metavar="GFF3"			)

parser.add_option(	"-t",
			"--type",
			dest="predictor",
			help="Type of predictions to keep.",
			metavar="TYPE"			)

parser.add_option(	"-o",
			"--output",
			dest="outfile",
			help="File to write output to (defaults to 'out.txt')",
			metavar="OUTPUT"		)

(options, args) = parser.parse_args()

if not options.interpro_file:
	parser.error("\nPlease specify Interpro output file.\n")
if not options.gff3:
	parser.error("\nPlease specify gff3 annotation file.\n")
if not options.predictor:
	options.predictor	= "Pfam"
if not options.outfile:
	options.outfile		= "out.txt"

with open(options.interpro_file, "r") as f:
	lines	= f.readlines()

names_starts 	= dict()
names_ends 	= dict()
domains		= dict()
identifiers	= dict()

for line in lines:
	fields	= line.split("\t")
	name	= fields[0]
	start	= fields[6]
	end	= fields[7]
	if fields[3] == options.predictor:
		if name not in names_starts:
			names_starts[name] 	= [fields[6]]
			names_ends[name]	= [fields[7]]
			domains[name]		= [fields[5]]
			identifiers[name]	= [fields[4]]
		else:
			names_starts[name]	= names_starts[name] + [fields[6]]
			names_ends[name]	= names_ends[name] + [fields[7]]
			domains[name]		= domains[name] + [fields[5]]
			identifiers[name]	= identifiers[name] + [fields[4]]
if not domains:
	sys.exit("\nDomain type " + options.predictor + " does not exist in file provided.\n")

with open(options.gff3, "r") as f:
	lines = f.readlines()

coords			= dict()

output	= open(options.outfile, "w")

for line in lines:
	if not re.match("#", line):
		fields		= line.split("\t")
		scaffold	= fields[0]
		feature_type	= fields[2]
		start		= int(fields[3])
		end		= int(fields[4])
		attributes	= fields[8].split(";")
		if feature_type == "gene":
			name		= attributes[0]
			name		= name.split("=")
			name		= name[1]
			if name in names_starts:
				for i in range(0, len(names_starts[name])):
					dstart		= int(names_starts[name][i]) * 3 - 2
					dend		= int(names_ends[name][i]) * 3 - 2
					domain		= domains[name][i]
					identifier	= identifiers[name][i]
					tdstart		= start + dstart
					tdend		= start + dend
					out		= [name, str(start), str(end), str(tdstart), str(tdend), identifier, domain]
					out		= "\t".join(out)
					output.write(out + "\n")


