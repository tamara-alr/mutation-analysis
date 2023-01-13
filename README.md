# mutation-analysis
This pipeline can be used to analyse a large number of SARS-CoV-2 sequences and to detect SARS-CoV-2 AA changes. A sample dataset of 10,000 sequences and corresponding metadadta is provided here. 
Steps involve: file chunking, sequence filtering, grouping by geographical regions, epiweek and lineage analysis, mutation detection. Mutational detection here outlined for Nucleocapsid Protein (NP) ORF but coordinates can be changed to suit any ORF.

# Requirements
1. R version 4.2.2
2. R packages: R.utils, seqinr, ape
3. MAFFT v7.49

# Analysis steps for sample dataset
R scripts are passed as arguments from command line

1. Download sample_dataset to desired path (contains sequences, metadata, refseq and all_R_scripts)

2. FILE CHUNKING: Use extract_chunks.R to split 10,000 sequences into 2 chunks of 5000.  
ARGS: Rscript [path]/sample_dataset/all_R_scripts/extract_chunks.R [path]/sample sequences.fasta metadata.tsv 5000  

3. SEQUENCE FILTERING: Use run_filter_analysis.R on each chunk to extract high quality sequences, processing 1000 sequences at a time  
ARGS: Rscript [path]/sample_dataset/all_R_scripts/run_filter_analysis.R [path]/sample_dataset/chunk1 chunk1.fasta chunk1.tsv 1000  
ARGS: Rscript [path]/sample_dataset/all_R_scripts/run_filter_analysis.R [path]/sample_dataset/chunk2 chunk2.fasta chunk2.tsv 1000  

4. SEQUENCE MERGING: Use merge_filtered.R to merge filtered sequences and metadata from all chunks in directory.  
Merged sequences will be stored in new directory: merged_filtered where subsequent analysis takes place.  
ARGS: Rscript [path]/sample_dataset/all_R_scripts/merge_filtered.R [path]/sample_dataset  

4. GROUPING BY GEO: Use group_regions.R to split merged dataset into 6 geographical regions: Africa, Asia, Europe, North America, Oceania, South America  
ARGS: Rscript [path]/sample_dataset/all_R_scripts/group_regions.R [path]/sample_dataset/merged_filtered merged_filtered_data.tsv merged_filtered_seqs.fasta  

5. ASSIGN EPIWEEK: Use assign_epiweek.R to assign epidemiological weeks to each geo dataset. Output csv contains ID, collection date, epiweek and lineage.  
ARGS: Rscript [path]/sample_dataset/all_R_scripts/assign_epiweek.R [path]/sample_dataset/merged_filtered/Africa Africa_data.tsv 10  

6. ANALYSE EPI: Use count_epi_n_lineage.R to count and frequency per epiweek and lineage in each geo dataset.  
ARGS: Rscript [path]/sample_dataset/all_R_scripts/count_epi_n_lineage.R [path]/sample_dataset/merged_filtered/Africa  

7. FIND MUTATIONS: Use find_N_mut.R to detect AA changes in NP (genomic coordinates: 28274 29533)  
ARGS: Rscript [path]/sample_dataset/all_R_scripts/assign_epiweek.R [path]/sample_dataset/merged_filtered/Africa Africa_seqs.fasta [path]/sample/Refseq/MN908947.fasta 28274 29533  

# GISAID EPI IDs
GISAD_EPI_SET file contains EPI_SET for the high-quality dataset (HQR) of 3,051,084 seqeunces used in study.


