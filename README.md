Program: HCPmiR (Identifying High-Confidence Plant miRNAs Using Multiple sRNA-Seq Data)
Version: 1.0
Contact: Chao Sun <sunchao@caf.ac.cn>
    Research Institute of Forestry, Chinese Academy of Forestry

1. Introduction
Using whole genome sequence and multiple sRNA-Seq data sets, HCPmiR identifies high-confidence
plant miRNAs through expression quantity filtering and results integration.

2. Installation
You must have successfuly installed DNApi, Cutadapt, Bowtie and Mireap.

DNApi (https://github.com/jnktsj/DNApi)
Cutadapt (https://github.com/marcelm/cutadapt)
Bowtie (https://github.com/BenLangmead/bowtie)
Mireap (https://github.com/liqb/mireap)

3. Usage
perl HCPmiR.pl -r genome -s sRNA_list -t TPM -d set_number -o out_prefix [options]*

Prerequiried tools for running: DNApi, cutadapt, Bowtie, mireap

       -r    [required] Fasta format of the reference genome sequence.
       -s    [required] File list of the sRNA-Seq data.
             The suffixes of the sRNA-Seq data should be '.fastq'.
             The file list and the sRNA-Seq data should be placed in the same directory. 
       -t    [required] TPM (transcripts per million) cutoff of miRNA identification. 
       -d    [required] Set number cutoff of miRNA identification.
       -o    [required] The prefix name of the output file.
       -p    [optional] The number of threads used for Bowtie alignment.
             Default: 5.
       -m    [optional] Maximal mapping number of sRNAs on the genome.
             Default: 30.  
       -h    [Optional] Print this usage message.

4. Output format
MIREAP produces one result file (*.list) and multiple directories named by sRNA-Seq id at each run.

*.list
This file contains miRNA genes and miRNAs. The file contains seven columns, which are miRNA 
gene id, chromosomes, miRNA gene starting position, miRNA gene termination position, strand, 
miRNA-5p, and miRNA-3p. Among them, the miRNA gene id is separated by commas, and the record 
represents the miRNA gene obtained by the sRNA-Seq data after running Mireap; the miRNA-5p 
and miRNA-3p information are separated by underline. Before the underline, it is separated 
by semicolons, which represents the location and sequence of the miRNA. After the underline, 
it is also separated by semicolons, indicating the expression values of the miRNA in sRNA-Seq
data. If it is not expressed, it is masked as "-". If miRNA-5p or miRNA-3p is not identified, 
it is represented by "NA".

The directories named by sRNA-Seq id
Their are six files in each directory: 
*.aln, *.gff and *.log are the output files of Mireap;
*TPM* is the gff file of miRNA genes and miRNAs filtered by expression values;
*.3adapt contains the 3' adapter sequences;
*.filtered.fa is the reads file prepared for Mireap.
