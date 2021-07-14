This README describes the scripts used for the sequence analysis in:   
[Antigenic evolution of human influenza H3N2 neuraminidase is constrained by charge balancing](https://www.biorxiv.org/content/10.1101/2021.07.10.451918v1?rss=1) Now available in BioRxiv.

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
* [./Fasta/Human_H3N2_HA_2020.aln.gz](./Fasta/Human_H3N2_HA_2020.aln.gz): Full-length HA protein sequences from human H3N2 downloaded from [GISAID](https://www.gisaid.org/)
* [./Fasta/Human_H3N2_NA_2020.aln.gz](./Fasta/Human_H3N2_NA_2020.aln.gz): Full-length NA protein sequences from human H3N2 downloaded from [GISAID](https://www.gisaid.org/)


## ANALYSIS PIPELINE
### Fitness landscape libraries analysis
1. [./script/fastq_to_fitness.py](./script/fastq_to_fitness.py): Converts raw reads to variant counts and fitness measures.
    - Input files:
      - Raw sequencing reads in fastq/ folder
      - [./data/SampleInfo.tsv](./data/SampleInfo.tsv)
      - [./Fasta/RefSeq.fa](./Fasta/RefSeq.fa)
    - Output files:
      - result/NA_Epi_*.tsv
2. [./script/complie_fit_result.py](./script/compile_fit_result.py): Complie variants info(amino acid;charge;fitness) in six different genetic background
    - Input files:
      - [./data/lib_variants.tsv](./data/lib_variants.tsv)
      - [./Fasta/N2.fa](./Fasta/N2.fa)
      - [./data/WTseq.tsv](./data/WTseq.tsv)
      - result/NA_Epi_*.tsv
    - Output files:
      - [./result/NA_compile_results.tsv](./result/NA_compile_results.tsv)
3. [./script/complie_fit_result.py](./script/compile_fit_result.py): Complie variants info(amino acid;charge;fitness) in six different genetic background
    - Input files:
      - [./data/lib_variants.tsv](./data/lib_variants.tsv)
      - [./Fasta/N2.fa](./Fasta/N2.fa)
      - [./data/WTseq.tsv](./data/WTseq.tsv)
      - result/NA_Epi_*.tsv
    - Output files:
      - [./result/NA_compile_results.tsv](./result/NA_compile_results.tsv)
4.[./script/NAEpi_PrefEvol.py](./script/NAEpi_PrefEvol.py): Amino acid sequences of NA antigenic region of interest in naturally occurring strains were extracted
    - Input files:
      - [./result/NA_compile_results.tsv](./result/NA_compile_results.tsv)
      - [./Fasta/Human_H3N2_NA_2020.aln](./Fasta/Human_H3N2_NA_2020.aln)
    - Output files:
      - [./result/HumanN2_numseq_year.tsv](./result/HumanH3N2_HAecto_year.tsv)
      - [./result/Fits_ByYear.tsv](./result/Prefs_ByYear.tsv)
      - [./result/Motif_ByYear.tsv](./result/Motif_ByYear.tsv)
### Inference of additive fitness and pairwise epistasis
1. [./script/GE_regression.ipynb](./script/GE_regression.ipynb): Model training and robustness validation
    - Input files:
      - [./result/NA_compile_results.tsv](./result/NA_compile_results.tsv)
    - Output files:
      - result/*_epi.csv
      - result/*_add.csv
2. [./script/GE_regression_v2.ipynb](./script/GE_regression_v2.ipynb): Cross-validation and regularization
    - Input files:
      - [./result/NA_compile_results.tsv](./result/NA_compile_results.tsv)
    - Output files:
      - [./result/reg_r2.csv](./result/reg_r2.csv)
3. [./script/Distance_CA.py](./script/Distance_CA.py): Calculate C alpha-alpha Distance within NA antigenic region
    - Output files:
      - [./result/CA_distance.tsv](./result/CA_distance.tsv)
### Natural strains evolution analysis 
1. [./script/ExtractSweep.py](./script/ExtractSweep.py): Calculate all variant frequencies of each residue of H3N2 NA over year
    - Input files:
      - [./Fasta/Human_H3N2_NA_2020.aln](./Fasta/Human_H3N2_NA_2020.aln)
    - Output files:
      - [./result/HumanN2Sweep_All.tsv](./result/HumanN2Sweep_All.tsv)
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
4. [./script/Coevolution_analysis_NA.ipynb](./script/Coevolution_analysis_NA.ipynb): Analyze NA antigenic region charge state coevolution
    - Input files:
      - [./result/HumanH3N2_NA_classified.csv](./result/HumanH3N2_NA_classified.csv)
    - Output files:
      - [./result/compare_peak.csv](./result/compare_peak.csv)
      - [./result/Coevols.csv](./result/Coevols.csv)

## PLOT
### Natural evolution of an antigenic region in human H3N2 NA (Figure 1)
1. [./script/Plot_Mutation_year.R](./script/Plot_Mutation_year.R): plot HA and NA accumulating mutation since 1968 ([Supplementary Fig. 1](./graph/H3N2_mutation_year.png))
2. [./script/NA_epi.pml](./script/NA_epi.pml): plot the NA head domain ([Fig. 1a](./graph/NA_epi.png))
3. [./script/NA_epi_zoom.pml](./script/NA_epi_zoom.pml): plot the antigenic region of interest ([Fig. 1b](./graph/NA_epi_zoom.png))
4. [./script/TrackAAFreq.R](./script/TrackAAFreq.R):plot natural occurrence frequencies of the amino acid variants/charge state (Fig. 1c and Supplementary Fig. 13)
    - Input files:
      - [./result/HumanN2Sweep_All.tsv](./result/HumanN2Sweep_All.tsv)
    - Output files:
      - [./graph/NatMutFreq_roi.png](./graph/NatMutFreq_roi.png)
      - [./result/HumanH3N2_NA_classified.csv](./result/HumanH3N2_NA_classified.csv)
### Comparing the local fitness landscapes of the NA antigenic region (Figure 2)
1. [./script/Plot_CompareLib.R](./script/Plot_CompareLib.R): plot fitness distribution and correlations of different background (Fig. 2a and b)
    - Input files:
      - [./result/NA_compile_results.tsv](./result/NA_compile_results.tsv)
    - Output files:
      - [./graph/Lib_fit_sina.png](./graph/Lib_fit_sina.png)
      - [./graph/LibCorPairs.png](./graph/LibCorPairs.png)
2. [./script/Plot_Nat_motif_Freq.R](./script/Plot_Nat_motif_Freq.R):plot naturally occurring variant frequencies over year ([Supplementary Fig. 2](./graph/FreqByYear.png))
3. [./script/Plot_CompareRep.R](./script/Plot_CompareRep.R): plot the biological repeat correlation ([Supplementary Fig. 3](./graph/Compare_Rep.png))
4. [./script/Plot_TrackPref.R](./script/Plot_TrackPref.R): plot naturally occurring variants fitness over year(Fig. 2c)
    - Input files:
      - [./result/Fits_ByYear.tsv](./result/Fits_ByYear.tsv)
    - Output files:
      - [./graph/FitsByYear.png](./graph/FitsByYear.png)
5. [./script/Plot_NA_titer.R](./script/Plot_NA_titer.R): plot virus rescue experiment of WT strains([Supplementary Fig. 4](./graph/NA_titer.png))
### Inference of additive fitness and pairwise epistasis (Figure 3)
1. [./script/Plot_hyperpar_R2.R](./script/Plot_hyperpar_R2.R): plot evaluation of model hyperparameters using repeated k-fold cross-validation ([Supplementary Fig. 5](./graph/reg_r2.png)&(./graph/hyperpar_r.png))
2. [./script/Plot_add_heatmap.R](./script/Plot_add_heatmap.R): plot parameters for additive fitness in different genetic backgrounds (Fig. 3a)
    - Input files:
      - result/*_add.csv
    - Output files:
      - [./graph/add_heatmap.png](./graph/add_heatmap.png)
3. [./script/Plot_epi_heatmap.R](./script/Plot_epi_heatmap.R): plot pairwise epistasis heatmap and epistasis classified by charge states (Fig. 3b and Supplementary Fig. 8; Fig. 4a and Supplementary Fig. 9)
    - Input files:
      - result/*_epi.csv
    - Output files:
      - [./graph/EPI_heatmap.png](./graph/EPI_heatmap.png)
      - [./graph/classified_epi.png](./graph/classified_epi.png)
4. [./script/Plot_CorEPI.R](./script/Plot_CorEPI.R):plot correlation matrices of additive fitness and pairwise epistasis among six genetic backgrounds (Fig. 3c and Supplementary Fig. 6-7)
    - Input files:
      - result/*_add.csv
      - result/*_epi.csv
    - Output files:
      - [./graph/compare_ADD.png](./graph/compare_ADD.png)
      - [./graph/compare_EPI.png](./graph/compare_EPI.png)
      - [./graph/heatmap_epi_add_cor.png](./graph/heatmap_epi_add_cor.png)
### The importance of local net charge in the NA antigenic region (Figure 4)
1. [./script/NA_epi_bind.pml](./script/NA_epi_bind.pml): plot the NA antigenic region interaction ([Fig. 4b](./graph/NA_epi_bind1.png) and [Supplementary Fig. 10](./graph/NA_epi_bind2.png))
2. [./script/Plot_charge_vs_fit.R](./script/Plot_charge_vs_fit.R): Plot variant fitness with local net charge ([Fig. 4c](./graph/Compare_charge_vs_fit_bK79.png) and [Supplementary Fig. 12](./graph/Compare_charge_vs_fit.png))
3. [./script/Plot_Distance_vs_EPI](./script/Plot_Distance_vs_EPI): plot Cα-Cα distances and epistasis
    - Input files:
      - [./result/CA_distance.tsv](./result/CA_distance.tsv)
      - result/*_epi.csv
    - Output files:
      - [./graph/Distance_vs_epi.png](./graph/Distance_vs_epi.png)
### Predicting coevolution of charge states in the NA antigenic region using epistasis (Figure 5)
1. [./script/plot_charge_natural_strain.R](./script/plot_charge_natural_strain.R): plot the evolution of local net charge at the NA antigenic region (Fig. 5a)
    - Input files:
      - [./result/HumanH3N2_NA_charge.tsv](./result/HumanH3N2_NA_charge.tsv)
    - Output files:
      - [./graph/NA_year_vs_charge.png](./graph/NA_year_vs_charge.png)
2. [./script/Plot_epi_heatmap_bycharge.R](./script/Plot_epi_heatmap_bycharge.R): plot pairwise epistasis of charge states (Fig. 5b and Supplementary Figure 14)
    - Input files:
      - result/*_epi.csv
    - Output files:
      - [./graph/epi_charge.png](./graph/epi_charge.png)
3. [./script/Plot_epi_vs_Coevol.R](./script/Plot_epi_vs_Coevol.R):plot relationship between coevolution score and pairwise epistasis (Fig. 5c and Supplementary Fig. 16)
    - Input files:
      - [./result/Coevols.csv](./result/Coevols.csv)
      - result/*_epi.csv
    - Output files:
      - [./graph/compare_natural_epi_bk79.png](./graph/compare_natural_epi_bk79.png)
      - [./graph/compare_natural_epi.png](./graph/compare_natural_epi.png)
