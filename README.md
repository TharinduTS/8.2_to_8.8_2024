# 8.2_to_8.8_2024

First I want to extract needed regions INCLUDING INVARIENT regions using vcf and reference fasta

I am working in Beluga
```
/scratch/premacht/8_mb_region_2024
```
I have bam files and reference genome in respective directories

Created a bed file with region to extract

```8.2_to8.8_mb.bed
Chr7    8200000 8800000
```

Then extracted the region with following command
```
module load gatk/4.4.0.0 StdEnv/2023
for i in ../new_project_Apr_2023/new_data_Jun23/*.bam; do gatk HaplotypeCaller -R ../new_project_Apr_2023/new_data_Jun23/reference_genome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -I ${i} -L 8.2_to8.8_mb.bed -ERC BP_RESOLUTION -O ${i##../new_project_Apr_2023/new_data_Jun23/}8.2_to_8.8mb.vcf;done
```
made a directory and moved all the related vcf files there

```
mkdir 8.2_to_8.8mb_vcfs
mv *8.2_to_8.8mb.vcf* 8.2_to_8.8mb_vcfs
```
Copied vcf2phylip from
```
https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py
```
Then converted vcf to fasta
```
for i in *vcf;do python ../vcf2phylip.py --input ${i} --fasta;done
```
Combined all fastas
```
cat *.fasta >8.2_to_8.8mb_combined_fasta.fa
```
Aligned with MAFFT
```
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 mafft-mpi/7.471
mafft --auto 8.2_to_8.8mb_combined_fasta.fa > 8.2_to_8.8mb_combined_fasta.fa_aligned.fasta
```
Then I downloaded the fasta file and analysed it with following R script
```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

#install.packages("seqinr")
library(seqinr)

my_fasta<-read.fasta("Nucleotide alignment.fasta")
my_seq_df<-as.data.frame(my_fasta)



female_fixed <- my_seq_df[my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_Nigeria_EUA0333_combined__sorted.bam 
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_SierraLeone_AMNH17272_combined__sorted.bam
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_SierraLeone_AMNH17274_combined__sorted.bam
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam==my_seq_df$F_IvoryCoast_xen228_combined__sorted.bam
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_Ghana_ZY_BJE4360_combined__sorted.bam
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_Nigeria_EUA0334_combined__sorted.bam
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_Nigeria_EUA0335_combined__sorted.bam
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_SierraLeone_AMNH17271_combined__sorted.bam
                          & my_seq_df$F_Nigeria_EUA0331_combined__sorted.bam!=my_seq_df$M_SierraLeone_AMNH17273_combined__sorted.bam,]

male_fixed <- my_seq_df[my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_Nigeria_EUA0333_combined__sorted.bam 
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_SierraLeone_AMNH17272_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_SierraLeone_AMNH17274_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_IvoryCoast_xen228_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Ghana_ZY_BJE4360_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Nigeria_EUA0334_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Nigeria_EUA0335_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_SierraLeone_AMNH17271_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_SierraLeone_AMNH17273_combined__sorted.bam,]


# Without Nigeria

female_fixed_no_niger <- my_seq_df[my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam==my_seq_df$F_SierraLeone_AMNH17272_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam==my_seq_df$F_SierraLeone_AMNH17274_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam==my_seq_df$F_IvoryCoast_xen228_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$M_Ghana_ZY_BJE4360_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$M_SierraLeone_AMNH17271_combined__sorted.bam
                          & my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam!=my_seq_df$M_SierraLeone_AMNH17273_combined__sorted.bam,]

male_fixed_no_niger <- my_seq_df[my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_Ghana_WZ_BJE4687_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_SierraLeone_AMNH17272_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_SierraLeone_AMNH17274_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam!=my_seq_df$F_IvoryCoast_xen228_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_Ghana_ZY_BJE4360_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_SierraLeone_AMNH17271_combined__sorted.bam
                        & my_seq_df$M_Ghana_WY_BJE4362_combined__sorted.bam==my_seq_df$M_SierraLeone_AMNH17273_combined__sorted.bam,]


# This region is for sequence selecting with previos results

# Change these ******

SNP_location<-44900
upstream<-100
downstream<-900

#********************

extract_begin<-SNP_location-upstream
extract_end<-SNP_location+downstream
```



