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

def dict_append(a,b):
        x = a.copy()
        b.update(x)
        return b

with open(options.interpro_file, "r") as f:
	lines	= f.readlines()

domains		= dict()
identifiers	= dict()
seqlen		= dict()
starts		= dict()
ends		= dict()

for line in lines:
	fields	= line.split("\t")
	if fields[3] == options.predictor:
		if fields[0] not in domains:
			domains[fields[0]]	= [fields[5]]
			identifiers[fields[0]]	= [fields[4]]
			seqlen[fields[0]]	= [fields[2]]
			starts[fields[0]]	= [fields[6]]
			ends[fields[0]]		= [fields[7]]
		else:
			domains[fields[0]]	= domains[fields[0]] + [fields[5]]
			identifiers[fields[0]]	= identifiers[fields[0]] + [fields[4]]
			seqlen[fields[0]]	= seqlen[fields[0]] + [fields[2]]
			starts[fields[0]]	= starts[fields[0]] + [fields[6]]
			ends[fields[0]]		= ends[fields[0]] + [fields[7]]

with open(options.gff3, "r") as f:
	lines	= f.readlines()

gcoords	= dict()
ccoords	= dict()

for line in lines:
	if not re.match("^#", line):
		fields	= line.split("\t")
		attrs	= fields[8].split(";")
		if fields[2] == "gene":
			(undef, name)	= attrs[0].split("=")
			gcoords[name]	= [fields[3], fields[4], fields[6]]
			
		if fields[2] == "CDS":
			if name in ccoords:
				ccoords[name]	= ccoords[name] + [fields[3], fields[4]]
			else:
				ccoords[name]	= [fields[3], fields[4]]

tcoords		= dict()
translator	= dict()
for key in gcoords.iterkeys():
	for i in xrange(0, len(ccoords[key]), 2):
		if i == 0:
			start	= 1
			end	= int(ccoords[key][i+1]) - int(ccoords[key][i]) + 1
		else:
			start	= end + 1
			end	= int(ccoords[key][i+1]) - int(ccoords[key][i]) + start
		tcoords	= dict_append(tcoords, {start : ccoords[key][i], end : ccoords[key][i+1]})
		translator[key]	= tcoords
	tcoords	= dict()

allpos	= dict()
for key in translator.iterkeys():
	if key in starts:
		count = 1
		for i in range(0, len(starts[key])):
			iikey		= key + "_" + str(i)
			if key not in allpos:
				allpos[key]	= {iikey: [int(starts[key][i]),int(ends[key][i]), translator[key]]}
			else:
				allpos[key]	= dict_append(allpos[key], {iikey : [int(starts[key][i]), int(ends[key][i]), translator[key]]})
dpos	= dict()
sortpos	= list()
for key in allpos.iterkeys():
	if gcoords[key][2] == "+":
		for ikey in allpos[key].iterkeys():
			sortpos	= list()
			for iikey in allpos[key][ikey][2].iterkeys():
				sortpos	= sortpos + [iikey]
			sortpos	= [allpos[key][ikey][0], allpos[key][ikey][1]] + sortpos
			sortpos	= sorted(sortpos)
			for i in range(0, len(sortpos)):
				if sortpos[i] in translator[key]:
					start	= translator[key][sortpos[i]]
					tstart	= sortpos[i]
				else:
					newpos	= int(start) + (int(sortpos[i]) - int(tstart))
					if key not in dpos:
						dpos[key]	= [newpos]
					else:
						dpos[key]	= dpos[key] + [newpos]
	else:
		for ikey in allpos[key].iterkeys():
			sortpos	= list()
			for iikey in allpos[key][ikey][2].iterkeys():
				sortpos	= sortpos +[iikey]
			sortpos	= [allpos[key][ikey][0], allpos[key][ikey][1]] + sortpos
			dstart		= sortpos[0]
			dend		= sortpos[1]
			del		sortpos[0]
			del 		sortpos[0]
			sortpos		= sorted(sortpos)
			dend		= sortpos[-1] - dend + 1
			dstart		= sortpos[-1] - dstart + 1
			sortpos		= sortpos + [dstart, dend]
			sortpos		= sorted(sortpos)
			for i in range(0, len(sortpos)):
				if sortpos[i] in translator[key]:
					start	= translator[key][sortpos[i]]
					tstart	= sortpos[i]
				else:
					newpos	= int(start) + (int(sortpos[i]) - int(tstart))
					if key not in dpos:
						dpos[key]	= [newpos]
					else:
						dpos[key]	= dpos[key] + [newpos]

with open(options.outfile, "w") as f:
	for key in dpos.iterkeys():
		domsort	= zip(starts[key], domains[key])
		idsort	= zip(starts[key], identifiers[key])
		domsort.sort()
		idsort.sort()
		n=0
		for i in xrange(0, len(dpos[key]) - 1, 2):
			output	= [str(key), str(gcoords[key][0]), str(gcoords[key][1]), str(dpos[key][i]), str(dpos[key][i+1]), str(domsort[n][1]), str(idsort[n][1])]
			output	= "\t".join(output)
			f.write(output + "\n")
			n += 1
