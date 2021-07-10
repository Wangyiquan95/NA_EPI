This README describes the scripts used for the sequence analysis in:   
[Antigenic evolution of human influenza H3N2 neuraminidase is constrained by charge balancing](xxx)

## ANALYSIS FOR H3N2 NA ANTIGENIC region of interest by DEEP MUTATIONAL SCANNING
This study aims to understand how epistasis influence NA antigenic evolution and characterize the underlying biophysical constraints. The repository here describes the analysis for the deep mutational scanning experiment that focuses on NA residues 328, 329, 344, 367, 368, 369, 370 in six different genetic backgrounds, namely A/Hong Kong/1/1968 (HK68), A/Bangkok/1/1979 (Bk79), A/Beijing/353/1989 (Bei89), A/Moscow/10/1999 (Mos99),A/Victoria/361/2011 (Vic11), and A/Hong Kong/2671/2019 (HK19). 

### REQUIREMENTS
* [Python](https://www.python.org/) version 3.6
* [R](https://www.r-project.org/) version 4.0
* [Weblogo](https://weblogo.berkeley.edu) version 3.6
* [MAVE-NN](https://github.com/jbkinney/mavenn)

### INPUT FILE
* All raw sequencing reads, which can be downloaded from NIH SRA database [PRJNA742436](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA742436), should be placed in fastq/ folder. The filename for read 1 should match those described in [./data/SampleInfo.tsv](./data/SampleInfo.tsv). The filename for read 2 should be the same as read 1 except "R1" is replaced by "R2"
* [./data/SampleInfo.tsv](./data/SampleInfo.tsv): Describes the sample identity for each fastq file
* [./Fasta/RefSeq.fa](./Fasta/RefSeq.fa): Reference (wild type) nucleotide sequences for the sequencing data
* [./Fasta/N2.fa](./Fasta/N2.fa):Reference (wild type) amino acid sequences for the sequencing data
* [./data/WTseq.tsv](./data/WTseq.tsv): Amino acids for the wild type sequences at residues 328, 329, 344, 367, 368, 369, 370
* 


### ANALYSIS PIPELINE
1. [./script/fastq\_to\_fitness.py](./script/fastq_to_fitness.py): Converts raw reads to variant counts and fitness measures.
    - Input files:
      - Raw sequencing reads in fastq/ folder
      - [./data/SampleInfo.tsv](./data/SampleInfo.tsv)
      - [./Fasta/RefSeq.fa](./Fasta/RefSeq.fa)
    - Output files:
      - result/NA\_Epi\_\*.tsv
2. [./script/complie\_fit\_result.py](./script/compile_fit_result.py): Complie variants info(amino acid;charge;fitness) in six different genetic background
    - Input files:
      - [./data/lib\_variants.tsv](./data/lib_variants.tsv)
      - [./Fasta/N2.fa](./Fasta/N2.fa)
      - [./data/WTseq.tsv](./data/WTseq.tsv)
      - result/NA\_Epi\_\*.tsv
    - Output files:
      - [./result/NA\_compile\_results.tsv](./result/NA_compile_results.tsv)
