# GRAM
GRAM: A GeneRAlized Model to predict the molecular effect of a non-coding variant in a cell type-specific manner


=Variant Prioritization=

==A. Dependencies==
The following tools are required: <br>
* sed, awk, grep <br>
* [http://tools.genes.toronto.edu/deepbind/ DeepBind] (version deepbind-v0.11) <br>
* [http://code.google.com/p/bedtools/downloads/list bedtools] (version bedtools-2.17.0) <br>
* [http://sourceforge.net/projects/samtools/files/tabix/ tabix] (version tabix-0.2.6 and up) <br>
*  [http://www.r-project.org/ R] (require packages: andomForest, glmnet, reshape2, gplots) <br>

==B. Tool Download==

This is a  Linux/UNIX-based tool. At the command-line prompt, type the following. 
 $ git clone https://github.com/gersteinlab/GRAM.git

==C. Pre-built Data Context==
All of the data can be downloaded under ‘Downloads’ in the Funseq3 web server (http://funseq3.gersteinlab.org/). 

==D. Tool Usage==
 $ cd gram/


To display the usage of tool, type ‘./grammar.sh -h’. <br>
 * Usage : bash ./grammar.sh -i bed -e exp
        Options :
                	-i		[Required] User Input SNVs file (BED format: chr st ed ref mut sample-id rsid)
                	-e	 	[Required] User Input gene expression matrix
                  -h     help message
             
                	
                	NOTE: Please make sure you have sufficient memory, at least 3G.

-i : Required format: chr st ed ref mut sample-id rsid<br>
-e: The rows correspond to genes and columns correspond to samples. Sample ids need to match with those in the variant bed file. <br>

==E. Input files==
* User input SNV file (-i): BED format 
 
In addition to the three required BED fields, please prepare your files as following (5 required fields, tab delimited; 
 the 6th column is reserved for sample names, do not put other information there): 
 chromosome, start position, end position, reference allele, alternative allele, sample id, rsid.
        Chromosome - name of the chromosome (e.g. chr3, chrX)
        Start position - start coordinates of variants. (0-based)
        End position - end coordinates of variants. (end exclusive)
                e.g., chr1   0     100  spanning bases numbered 0-99
        Reference allele - germlime allele of variants
        Alternative allele - mutated allele of variants
        Sample id - the sample id, specifying the input sample or cell line (e.g. "Patient-1", "GM12878")
        RSID -  the id for the variant (e.g. rs9347341)

e.g.
        chr2	242382863	242382864	C	T	Patient-1	rs3771570
        chr2 	242382863	242382864	C	T	Patient-2	rs3771570
        chr6	117210051 	117210052	T	C	Patient-3	rs339331
        …	…		…		…	…	…		…


* User input expression matrix (-e): The gene expression file should be prepared as a matrix with first column stores gene names (use gene symbols) and first row as sample names. Other fields are gene expression data either in RPKM or raw read counts format. Tab delimited. 
e.g.
        Gene	Sample1	Sample2	Sample3	Sample4	…
        A1BG	1	5	40	0	…
        A1CF	20	9	0	23	…
        …	…	…	…	…	…

==F. Output files==
Five output files will be generated: ‘Output.format’, ‘Output.indel.format’, ‘Recur.Summary’, ‘Candidates.Summary’ and ‘Error.log’. Output.format: stores detailed results for all samples; Output.indel.format: contains results for indels; Recur.Summary: the recurrence result when having multiple samples; Candidates.Summary: brief output of potential candidates (coding nonsynonymous/prematurestop variants, non-coding variants with score (>= 5 of un-weighted scoring scheme and >=1.5 for weighted scoring scheme) and variants in or associated with known cancer genes); Error.log: error information. For un-weighted scoring scheme, each feature is given value 1. 

When provided with gene expression files, two additional files will be produced – ‘DE.gene.txt’ stores differentially expressed genes and ‘DE.pdf ’is the differential gene expression plot. 

* Sample GRAM output

        chr6 29897016 29897017 C A Patient-1 6:29897017:C:A_rs7767188 0.58  0.494 -0.500289915458444  0.42  0.439 0.415735691763279
        chr6 29897016 29897017 C A Patient-2 6:29897017:C:A_rs7767188 0.58  0.494 -0.500289915458444  0.377 0.439 0.337216556180445
        chr6 29915765 29915766 C A Patient-3 6:29915766:C:A_rs7767188 0.296 0.286 -0.0699306742423629 0.351 0.666 0.325006854130755
        chr6 29915765 29915766 C A Patient-4 6:29915766:C:A_rs7767188 0.296 0.286 -0.0699306742423629 0.35  0.666 0.323297953713459

Columns:
        1: (chr) name of the chromosome
        2: (st) start coordinates of variants. (0-based)
        3: (ed) end coordinates of variants. (end exclusive)
        4: (ref) reference allele of variants
        5: (mut) mutant allele of variants
        6: (sampleid) the ID of the sample
        7: (snpid) the ID of the SNV
        8: (ref.enhAct) general regulatory activity of the reference allele
        9: (alt.enhAct) general regulatory activity of the mutant allele
        10: (logodds) predicted cell type modifier score
        11: (rank.las)
        12: (db.las)
        13: (score) predicted GRAM score


=Contact=
shaoke DOT lou AT yale DOT edu

This software is freely licensed under the Creative Commons license (Attribution-NonCommerical). The main aspects of this license are that: The work can be made available for non-commercial use Derivatives can be made of the work Derivatives do not have to be made available under the same terms that they were first used, and We should be cited.
