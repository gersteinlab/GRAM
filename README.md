# GRAM
GRAM: A GeneRAlized Model to predict the molecular effect of a non-coding variant in a cell type-specific manner

## A. Dependencies==
The following tools are required: <br>
* sed, awk, grep <br>
* <a href="http://tools.genes.toronto.edu/deepbind/"> DeepBind (version deepbind-v0.11) </a><br>
* <a href="http://code.google.com/p/bedtools/downloads/list"> bedtools (version bedtools-2.17.0)</a> <br>
* <a href="http://sourceforge.net/projects/samtools/files/tabix/"> tabix (version tabix-0.2.6 and up) </a> <br>
* <a href="http://www.r-project.org/"> R </a> (require packages: andomForest, glmnet, reshape2, gplots) <br>

## B. Tool Download

This is a  Linux/UNIX-based tool. At the command-line prompt, type the following. <br>
 $ git clone https://github.com/gersteinlab/GRAM.git

## C. Pre-built Data Context
Before you run grammar, you must provide some information in grammar scripts:<br>
genome=<path of genome fasta file>, chromomosome name should start with 'chr'> <br>
gencode=<gene code cds bed>, we have provide hg19 gencode cds file with this repo, just decompress it <br>
dpath=<deepbind folder >, deepbind fold should include 'db' and 'db/params', which are the parameters for TF binding models <br>
path_funseq=<funseq whole score>. it can be downloaded under ‘Downloads’ in the Funseq3 web server (http://funseq3.gersteinlab.org/). <br>

## D. Tool Usage
 $ cd gram/


To display the usage of tool, type ‘./grammar -h’. <br>
 * Usage : bash ./grammar -i bed -e exp <br>
        Options :<br>
                	-i		[Required] User Input SNVs file (BED format: chr st ed ref mut sample-id rsid)<br>
                	-e	 	[Required] User Input gene expression matrix <br>
                 -h     help message <br>
             
                	
                	NOTE: Please make sure you have sufficient memory, at least 3G.

-i : Required format: chr st ed ref mut sample-id rsid<br> 
-e: The rows correspond to genes and columns correspond to samples. Sample ids need to match with those in the variant bed file. <br>

## E. Input files
* User input SNV file (-i): BED format 
 
In addition to the three required BED fields, please prepare your files as following (5 required fields, tab delimited; <br>
 the 6th column is reserved for sample names, do not put other information there): <br>
 chromosome, start position, end position, reference allele, alternative allele, sample id, rsid.<br>
        Chromosome - name of the chromosome (e.g. chr3, chrX)<br>
        Start position - start coordinates of variants. (0-based)<br>
        End position - end coordinates of variants. (end exclusive)<br>
                e.g., chr1   0     100  spanning bases numbered 0-99<br>
        Reference allele - germlime allele of variants<br>
        Alternative allele - mutated allele of variants<br>
        Sample id - the sample id, specifying the input sample or cell line (e.g. "Patient-1", "GM12878")<br>
        RSID -  the id for the variant (e.g. rs9347341)<br>

e.g.<br>
        chr2	242382863	242382864	C	T	Patient-1	rs3771570<br>
        chr2 	242382863	242382864	C	T	Patient-2	rs3771570<br>
        chr6	117210051 	117210052	T	C	Patient-3	rs339331<br>
        …	…		…		…	…	…		…


* User input expression matrix (-e): <br>The gene expression file should be prepared as a matrix with first column stores gene names (use gene symbols) and first row as sample names. Other fields are gene expression data either in RPKM or raw read counts format. Tab delimited. <br>
e.g.<br>
        Gene	Sample1	Sample2	Sample3	Sample4	…<br>
        A1BG	1	5	40	0	…<br>
        A1CF	20	9	0	23	…<br>
        …	…	…	…	…	…<br>

## F. Output files
Five output files will be generated: ‘Output.format’, ‘Output.indel.format’, ‘Recur.Summary’, ‘Candidates.Summary’ and ‘Error.log’. Output.format: stores detailed results for all samples; Output.indel.format: contains results for indels; Recur.Summary: the recurrence result when having multiple samples; Candidates.Summary: brief output of potential candidates (coding nonsynonymous/prematurestop variants, non-coding variants with score (>= 5 of un-weighted scoring scheme and >=1.5 for weighted scoring scheme) and variants in or associated with known cancer genes); Error.log: error information. For un-weighted scoring scheme, each feature is given value 1. 

When provided with gene expression files, two additional files will be produced – ‘DE.gene.txt’ stores differentially expressed genes and ‘DE.pdf ’is the differential gene expression plot. 

* Sample GRAM output<br>

        chr6 29897016 29897017 C A Patient-1 6:29897017:C:A_rs7767188 0.58  0.494 -0.500289915458444  0.42  0.439 0.415735691763279
        chr6 29897016 29897017 C A Patient-2 6:29897017:C:A_rs7767188 0.58  0.494 -0.500289915458444  0.377 0.439 0.337216556180445
        chr6 29915765 29915766 C A Patient-3 6:29915766:C:A_rs7767188 0.296 0.286 -0.0699306742423629 0.351 0.666 0.325006854130755
        chr6 29915765 29915766 C A Patient-4 6:29915766:C:A_rs7767188 0.296 0.286 -0.0699306742423629 0.35  0.666 0.323297953713459

Columns:
       1: (chrome) name of the chromosome <br>
       2: (start) start coordinates of variants. (0-based)<br>
       3: (end) end coordinates of variants. (end exclusive)<br>
       4: (ref) reference allele of variants<br>
       5: (mut) mutant allele of variants<br>
       6: (sampleid) the ID of the sample<br>
       7: (snpid) the ID of the SNV<br>
       8: (ref.enhAct) general regulatory activity of the reference allele<br>
       9: (alt.enhAct) general regulatory activity of the mutant allele<br>
       10: (logodds) logodds calculated from reference and mutant allele regulatory activity<br>
       11: (expr.modifier) cell type modifier score predicted from TF expression<br>
       12: (binding.modifier) cell type modifier score predicted from TF binding<br>
       13: (gram.prob) predicted GRAM score<br>


## Contact
shaoke DOT lou AT yale DOT edu

This software is freely licensed under the Creative Commons license (Attribution-NonCommerical). The main aspects of this license are that: The work can be made available for non-commercial use Derivatives can be made of the work Derivatives do not have to be made available under the same terms that they were first used, and We should be cited.
