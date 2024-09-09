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
Then copied it into a different directory to blast and blasted it
```
module load StdEnv/2023  gcc/12.3 blast+/2.14.1
makeblastdb -in ../../new_project_Apr_2023/new_data_Jun23/reference_genome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -title reference -dbtype nucl -out ref.fa
blastn -db ./reference/ref.fa -query full_sequence_of_selected_sample  -outfmt "6 qstart qend sstart  send length qseq sseq " | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > best_single_hits.blastn

```

#********************************************************

Aligned with MAFFT
```
module load StdEnv/2020  gcc/9.3.0  openmpi/4.0.3 mafft-mpi/7.471
mafft --auto 8.2_to_8.8mb_combined_fasta.fa > 8.2_to_8.8mb_combined_fasta.fa_aligned.fasta
```
Then I downloaded the fasta file and analysed it with following R script
This extract the regions with male or female fixed regions and saves marked and unmarked text files
```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))
library(seqinr)

#************************************

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



#mark locations with male and female fixed SNPs

library(tibble)
male_fixed_no_niger_with_loc <- tibble::rownames_to_column(male_fixed_no_niger, "location")

female_fixed_no_niger_with_loc <- tibble::rownames_to_column(female_fixed_no_niger, "location")

my_seq_df_with_loc <- tibble::rownames_to_column(my_seq_df, "location")

male_location_list<-male_fixed_no_niger_with_loc$location

female_location_list<-female_fixed_no_niger_with_loc$location



#Marking male fixed SNPs

for (i in 1:length(male_location_list)) {
  for (j in 1:length(my_seq_df_with_loc$location)) {
    if (my_seq_df_with_loc$location[j]==male_location_list[i]) {
      my_seq_df_with_loc$F_Ghana_WZ_BJE4687_combined__sorted.bam[j]<-paste("*",my_seq_df_with_loc$F_Ghana_WZ_BJE4687_combined__sorted.bam[j],"*",sep = "")
    }
  }
}


#Marking female fixed SNPs

for (i in 1:length(female_location_list)) {
  for (j in 1:length(my_seq_df_with_loc$location)) {
    if (my_seq_df_with_loc$location[j]==female_location_list[i]) {
      my_seq_df_with_loc$F_Ghana_WZ_BJE4687_combined__sorted.bam[j]<-paste("(",my_seq_df_with_loc$F_Ghana_WZ_BJE4687_combined__sorted.bam[j],")",sep = "")
    }
  }
}


# This region is for extracting sequences as 1000 bp regions

# Creating directory to save text files

dir.create('marked_sequence_files')
dir.create('unmarked_sequence_files')

locations_left<-as.numeric(male_location_list)
while (0<length(locations_left)) {
  print(paste("extracting region "))


SNP_location<-min(as.numeric(locations_left))

# Change these to change sequence length******
upstream<-100
downstream<-900

#********************

extract_begin<-SNP_location-upstream
extract_end<-SNP_location+downstream

extracted_sequence<-my_seq_df_with_loc$F_Ghana_WZ_BJE4687_combined__sorted.bam[extract_begin:extract_end]

selected_sequence<-paste(extracted_sequence,collapse = '')

selected_sequence

print(paste("extracted SNPs from",extract_begin,"to",extract_end))

saving_file_name<-paste('Sex_sp_locations_marked_sequence_from_',extract_begin,'_to_',extract_end,'.txt.',sep = '')

print(paste("saving the sequence as",saving_file_name,"Males fixed SNPs are marked by * and Female fixed SNPs are marked by ()"))

writeLines(selected_sequence,paste('./marked_sequence_files/',saving_file_name,sep = ''))

# remove markings for pure sequence

unmarked_selected_sequence<-stringr::str_replace_all(selected_sequence, '\\*', '')
unmarked_selected_sequence<-stringr::str_replace_all(unmarked_selected_sequence, '\\(', '')
unmarked_selected_sequence<-stringr::str_replace_all(unmarked_selected_sequence, '\\)', '')

saving_file_name_unmarked<-paste('Sex_sp_locations_unmarked_sequence_from_',extract_begin,'_to_',extract_end,'.txt.',sep = '')

writeLines(unmarked_selected_sequence,paste('./unmarked_sequence_files/',saving_file_name_unmarked,sep = ''))

locations_left<-subset(locations_left,locations_left>extract_end)

}

full_sample_sequence<-my_seq_df_with_loc$F_Ghana_WZ_BJE4687_combined__sorted.bam

full_sequence<-paste(extracted_sequence,collapse = '')

# remove markings for pure sequence

unmarked_full_sequence<-stringr::str_replace_all(full_sequence, '\\*', '')
unmarked_full_sequence<-stringr::str_replace_all(unmarked_full_sequence, '\\(', '')
unmarked_full_sequence<-stringr::str_replace_all(unmarked_full_sequence, '\\)', '')

writeLines(unmarked_full_sequence,"./full_sequence_of_selected_sample")


```



