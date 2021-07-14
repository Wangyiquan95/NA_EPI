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
## Fitness landscape libraries analysis
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
3. [./script/complie\_fit\_result.py](./script/compile_fit_result.py): Complie variants info(amino acid;charge;fitness) in six different genetic background
    - Input files:
      - [./data/lib\_variants.tsv](./data/lib_variants.tsv)
      - [./Fasta/N2.fa](./Fasta/N2.fa)
      - [./data/WTseq.tsv](./data/WTseq.tsv)
      - result/NA\_Epi\_\*.tsv
    - Output files:
      - [./result/NA\_compile\_results.tsv](./result/NA_compile_results.tsv)
4.[./script/NAEpi\_PrefEvol.py](./script/NAEpi_PrefEvol.py): Amino acid sequences of NA antigenic region of interest in naturally occurring strains were extracted
    - Input files:
      - [./result/NA_compile\_results.tsv](./result/NA_compile_results.tsv)
      - [./Fasta/Human_H3N2\_NA\_2020.aln](./Fasta/Human_H3N2_NA_2020.aln)
    - Output files:
      - [./result/HumanN2_numseq_year.tsv](./result/HumanH3N2_HAecto_year.tsv)
      - [./result/Fits\_ByYear.tsv](./result/Prefs\_ByYear.tsv)
      - [./result/Motif\_ByYear.tsv](./result/Motif\_ByYear.tsv)
## Inference of epistasis and additive fitness effect
1. [./script/GE\_regression.ipynb](./script/GE_regression.ipynb): Model training and robustness validation
    - Input files:
      - [./result/NA\_compile\_results.tsv](./result/NA_compile_results.tsv)
    - Output files:
      - result/*_epi.csv
      - result/*_add.csv
2. [./script/GE\_regression\_v2.ipynb](./script/GE_regression_v2.ipynb): Cross-validation and regularization
    - Input files:
      - [./result/NA\_compile\_results.tsv](./result/NA_compile_results.tsv)
    - Output files:
      - [./result/reg\_r2.csv](./result/reg_r2.csv)
## Natural strains sequence analysis 
1. [./script/ExtractSweep.py](./script/ExtractSweep.py): Calculate all variant frequencies of each residue of H3N2 NA over year
    - Input files:
      - [./Fasta/Human\_H3N2\_NA\_2020.aln](./Fasta/Human_H3N2_NA_2020.aln)
    - Output files:
      - [./result/HumanN2Sweep\_All.tsv](./result/HumanN2Sweep_All.tsv)
2. [./script/analyze_charge_natural_strain.py](./script/analyze_charge_natural_strain.py): Analyze antigenic region local charge in natural circulating strains
    - Input files:
      - [./Fasta/Human_H3N2_NA_2020.aln](./Fasta/Human_H3N2_NA_2020.aln)
    - Output files:
      - [./result/HumanH3N2_NA_charge.tsv](./result/HumanH3N2_NA_charge.tsv)
3. [./script/mut_freq_ByYear.py](./script/mut_freq_ByYear.py): Analyze HA and NA accumulating mutation since 1968
    - Input files:
      - [./Fasta/Human_H3N2_NA_2020.aln](./Fasta/Human_H3N2_NA_2020.aln)
      - [./Fasta/Human_H3N2_HA_2020.aln](./Fasta/Human_H3N2_HA_2020.aln)
    - Output files:
      - [./result/HumanH3N2_mutation_year.tsv](./result/HumanH3N2_mutation_year.tsv)
## Coevolution analysis
1. [./script/Coevolution_analysis_NA.ipynb](./script/Coevolution_analysis_NA.ipynb): Analyze NA antigenic region charge state coevolution
    - Input files:
      - [./result/HumanH3N2_NA_classified.csv](./result/HumanH3N2_NA_classified.csv)
    - Output files:
      - [./result/compare_peak.csv](./result/compare_peak.csv)
      - [./result/Coevols.csv](./result/Coevols.csv)

### Plot
## Natural evolution of an antigenic region in human H3N2 NA
1. [./script/NA\_epi.pml](./script/NA_epi.pml): plot the NA head domain (Fig. 1a)
2. [./script/NA\_epi\_zoom.pml](./script/NA_epi_zoom.pml): plot the antigenic region of interest (Fig. 1b)
3. [./]()
