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
require("reshape2")

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
run_go		<- function(file, name, gene_file, type=NULL, pval=0.05)

{
	
	frame 		=  read.csv(gene_file, header=TRUE, sep="\t")
	goFrame		=  GOFrame(frame, organism=name)
	goAllFrame	=  GOAllFrame(goFrame)
	gsc		<- GeneSetCollection(goAllFrame, setType=GOCollection())
	genes		<- file$GENE_ID
	background	=  frame$GENE_ID
	
	params		<- suppressWarnings (	GSEAGOHyperGParams(	name = "Sclerotinia sclerotiorum",
									geneSetCollection = gsc,
									geneIds = genes,
									universeGeneIds = background,
									ontology = type,
									pvalueCutoff = pval,
									conditional = FALSE,
									testDirection = "over"	)	)

	or.test		<- suppressWarnings(	hyperGTest(params)	)
	return(or.test)

}

# Function for funning GOstats on all isolates in the dataframes produced by 'get_files'.
run_go_all	<- function(file_list, name, gene_file, type=NULL, pval=0.05)

{

	test.all	<- list()
	test.all.all	<- list()

	for (i in 1:length(file_list)){
		for (n in 1:length(file_list[[i]])){
			input		<- data.frame(file_list[[i]][[n]])
			test		<- run_go(input, name, gene_file, type=type, pval=pval)
			test.all[[n]]	<- test
		}

	test.all.all[[i]] <- test.all
	
	}

	return(test.all.all)

}

# Function to create table for ggplot. Will return dataframes for logodds and p-values with rownames as GO term descriptors.
make_plot	<- function(file_list, file_list.all)

{

	domains			<- c()
	domains.all		<- list()
	domains.logodds		<- list()
	domains.pval		<- list()
	domains.count		<- list()
	domains.pval.all	<- list()
	domains.logodds.all	<- list()
	domains.count.all	<- list()
	domains.count.tab	<- list()
	for (i in 1:length(file_list)){
		for (n in 1:length(file_list[[i]])){
			sum <- summary(file_list[[i]][[n]])$Term
			for (v in 1:length(sum)){
				if (!sum[v] %in% domains){
					domains <- c(domains, sum[v])
				}
			}
		}
	domains.all[[i]] 	<- domains
	domains			<- c()
	}
	for (i in 1:length(file_list.all)){
		for (n in 1:length(file_list.all[[i]])){
			logodds		<- summary(file_list.all[[i]][[n]])$OddsRatio
			pval		<- summary(file_list.all[[i]][[n]])$Pvalue
			count		<- summary(file_list.all[[i]][[n]])$Count
			domains.ns	<- summary(file_list.all[[i]][[n]])$Term
			for (v in 1:length(domains.ns)){
				domain	<- domains.ns[v]
				if (domain %in% domains.all[[i]]){
					domains.pval[[domain]]		<- c(domains.pval[[domain]], pval[v])
					domains.logodds[[domain]]	<- c(domains.logodds[[domain]], logodds[v])
					domains.count[[domain]]		<- c(domains.count[[domain]], count[v])
				}	
			}
			for (v in 1:length(domains.all[[i]])){
				domain	<- domains.all[[i]][v]
				if (! domain %in% domains.ns){
					domains.pval[[domain]]		<- c(domains.pval[[domain]], 0)
					domains.logodds[[domain]]	<- c(domains.logodds[[domain]], 0)
					domains.count[[domain]]		<- c(domains.count[[domain]], count[v])
				}
			}
		}
		domains.pval.tab	<- c()
		domains.logodds.tab	<- c()
		domains.count.tab	<- c()
		for (name in names(domains.pval)){
			domains.pval.tab	<- rbind(domains.pval.tab, domains.pval[[name]])
			domains.logodds.tab	<- rbind(domains.logodds.tab, domains.logodds[[name]])
			domains.count.tab	<- rbind(domains.count.tab, domains.count[[name]])
		}
		df.pval				<- data.frame(domains.pval.tab)
		df.logodds			<- data.frame(domains.logodds.tab)
		df.count			<- data.frame(domains.count.tab)
		df.pval				<- cbind(names(domains.pval), df.pval)
		df.logodds			<- cbind(names(domains.logodds), df.logodds)
		df.count			<- cbind(names(domains.count), df.count)
		colnames(df.pval)[1]		<- "NAME"
		colnames(df.logodds)[1]		<- "NAME"
		colnames(df.count)[1]		<- "NAME"
		xymelt.pval 			<- melt(df.pval, id.vars = "NAME")
		xymelt.logodds 			<- melt(df.logodds, id.vars = "NAME")
		xymelt.count 			<- melt(df.count, id.vars = "NAME")
		domains.pval.all[[i]]		<- xymelt.pval
		domains.logodds.all[[i]]	<- xymelt.logodds
		domains.count.all[[i]]		<- xymelt.count
		domains.logodds	<- list()
		domains.pval	<- list()
		domains.count	<- list()
	}
	domains.pval.logodds	<- list(domains.pval.all, domains.logodds.all, domains.count.all)
	return(domains.pval.logodds)

}
