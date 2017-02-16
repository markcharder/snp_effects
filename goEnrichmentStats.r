# This R source file can be used to run GOStats hypergeometric tests on GO terms associated with arbitrary gene sets.
# The get_files function is used to read in three different types of files. Each file has the following format:
# GENE_ID	EVIDENCE	GO_TERM
# Where gene ID is the ID of the gene, the evidence is the GO term evidence code for assignment of a particular GO term and
# GO_TERM is the GO term itself. There may be duplicate genes and duplicate GO terms, but not duplicate genes and GO terms that match.
# The 'files.h' object handles this file for genes in which snpEff has identified SNPs with a high impact on coding sequence -
# such as causing a premature stop, frame shift or intron/exon boundary change. The 'files.m' and 'files.l' objects contain genes that,
# according to SNP calling and snpEff analysis contain SNPs that have a moderate and low effect on protein sequence, respectively.
# A moderate effect is a non-synonymous SNP and a low effect is a synonymous SNP, according to snpEff documentation.
# Each of the files can either be specified exactly or with a regex to denote a list of files in a directory.
# However, the lists of each file type must be the same length.

require("ALL")
require("hgu95av2.db")
require("GO.db")
require("annotate")
require("genefilter")
require("GOstats")
require("RColorBrewer")
require("xtable")
require("Rgraphviz")
require("GSEABase")

# Function for reading large number of files.
get_files	<- function(pattern.h, pattern.m, pattern.l)

{

	files.h	<- list.files(pattern	= pattern.h)
	files.m	<- list.files(pattern	= pattern.m)
	files.l	<- list.files(pattern	= pattern.l)

	data.frames.h	<- list()
	data.frames.m	<- list()
	data.frames.l	<- list()

	for (i in 1:length(files.h)){
		data.frames.h[[i]] <- data.frame(read.csv(files.h[i], sep="\t", header=TRUE))
		data.frames.m[[i]] <- data.frame(read.csv(files.m[i], sep="\t", header=TRUE))
		data.frames.l[[i]] <- data.frame(read.csv(files.l[i], sep="\t", header=TRUE))
	}

	data.frames.all	<- list(data.frames.h, data.frames.m, data.frames.l)
	return (data.frames.all)
}

# Function for running GOstats.
run_go		<- function(file, name, gene_file)

{
	
	goFrame		=  GOFrame(file, organism=name)
	goAllFrame	=  GOAllFrame(goFrame)
	gsc 		<- GeneSetCollection(goAllFrame, setType=GOCollection())
	background	=  as.vector(file$GENE_ID)
	genes		<- data.frame(read.csv(gene_file, sep="\t", header=TRUE))
	genes		=  as.vector(genes$GENE_ID)

	
	params		<- suppressWarnings (	GSEAGOHyperGParams(	name = "Sclerotinia sclerotiorum",
									geneSetCollection = gsc,
									geneIds = genes,
									universeGeneIds = background,
									ontology = "MF",
									pvalueCutoff = 0.05,
									conditional = FALSE,
									testDirection = "over"	)	)

	or.test		<- suppressWarnings(	hyperGTest(params)	)
	return(or.test)

}

# Function for funning GOstats on all isolates in the dataframes produced by 'get_files'.
run_go_all	<- function(file_list, name, gene_file)

{

	test.all	<- list()
	test.all.all	<- list()

	for (i in 1:length(file_list)){
		for (n in 1:length(file_list[[i]])){
			input		<- data.frame(file_list[[i]][[n]])
			test		<- run_go(input, name, gene_file)
			test.all[[n]]	<- test
		}

	test.all.all[[i]] <- test.all
	
	}

	return(test.all.all)

}

# Function to create table for ggplot.


# Function for plotting GO terms with GGPlot.
