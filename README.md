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
mafft --auto 8.2_to_8.8mb_combined_fasta.fa > Nucleotide_alignment.fasta
```
Then I downloaded the fasta file and analysed it with following R script
This extract the regions with male or female fixed regions and saves marked and unmarked text files
```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))
library(seqinr)


#************************************
#install.packages("RFLPtools")
library(RFLPtools)



my_fasta<-read.fasta("Nucleotide_alignment.fasta")
my_seq_df<-as.data.frame(my_fasta)

#Drop any rows with 'n' in any sample

#my_seq_df<-my_seq_df[rowSums(my_seq_df == "n")==0, , drop = FALSE]

#add location as a column
my_seq_df_with_loc <- tibble::rownames_to_column(my_seq_df, "location")

#export sequence from a selected sample to blast

# ********EDIT SELECTED SAMPLE ACCORDINGLY******

selected_sample<-'F_Ghana_WZ_BJE4687_combined__sorted.bam'

#**************************************************

full_sample_sequence<-my_seq_df_with_loc[[selected_sample]]

full_sequence<-paste(full_sample_sequence,collapse = '')

# remove markings for pure sequence

unmarked_full_sequence<-stringr::str_replace_all(full_sequence, '\\*', '')
unmarked_full_sequence<-stringr::str_replace_all(unmarked_full_sequence, '\\(', '')
unmarked_full_sequence<-stringr::str_replace_all(unmarked_full_sequence, '\\)', '')

writeLines(unmarked_full_sequence,"./full_sequence_of_selected_sample")

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


male_location_list<-male_fixed_no_niger_with_loc$location

female_location_list<-female_fixed_no_niger_with_loc$location

#Marking male fixed SNPs

for (i in 1:length(male_location_list)) {
  for (j in 1:length(my_seq_df_with_loc$location)) {
    if (my_seq_df_with_loc$location[j]==male_location_list[i]) {
      my_seq_df_with_loc[[selected_sample]][j]<-paste("*",my_seq_df_with_loc[[selected_sample]][j],"*",sep = "")
    }
  }
}


#Marking female fixed SNPs

for (i in 1:length(female_location_list)) {
  for (j in 1:length(my_seq_df_with_loc$location)) {
    if (my_seq_df_with_loc$location[j]==female_location_list[i]) {
      my_seq_df_with_loc[[selected_sample]][j]<-paste("(",my_seq_df_with_loc[[selected_sample]][j],")",sep = "")
    }
  }
}


# This region is for extracting sequences as 1000 bp regions

# Creating directory to save text files

dir.create('marked_sequence_files')
dir.create('unmarked_sequence_files')
#dir.create('to_blast')

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

extracted_sequence<-my_seq_df_with_loc[[selected_sample]][extract_begin:extract_end]

selected_sequence<-paste(extracted_sequence,collapse = '')

selected_sequence

print(paste("extracted SNPs from",extract_begin,"to",extract_end))

saving_file_name<-paste('Sex_sp_locations_marked_sequence_from_',extract_begin,'_to_',extract_end,'.txt',sep = '')
saving_file_name_unmarked<-paste('Sex_sp_locations_unmarked_sequence_from_',extract_begin,'_to_',extract_end,'.txt',sep = '')

print(paste("saving the sequence as",saving_file_name,"Males fixed SNPs are marked by * and Female fixed SNPs are marked by () in",selected_sample))

#removing file if already exists to avoid conflicts

file.remove(paste('./marked_sequence_files/',saving_file_name,sep = ''))
file.remove(paste('./unmarked_sequence_files/',saving_file_name_unmarked,sep = ''))
#file.remove(paste('./to_blast/',saving_file_name_unmarked,sep = ''))

#writing marked file only for selected file
first_line<-paste(">",selected_sample,sep = '')
write(first_line,file = paste('./marked_sequence_files/',saving_file_name,sep = ''),append = TRUE)
write(selected_sequence,file = paste('./marked_sequence_files/',saving_file_name,sep = ''),append = TRUE)

#writing unmarked file only for blast

# remove markings for pure sequence

# unmarked_selected_sequence<-stringr::str_replace_all(selected_sequence, '\\*', '')
# unmarked_selected_sequence<-stringr::str_replace_all(unmarked_selected_sequence, '\\(', '')
# unmarked_selected_sequence<-stringr::str_replace_all(unmarked_selected_sequence, '\\)', '')
# 
# write(first_line,file = paste('./to_blast/',saving_file_name_unmarked,sep = ''),append = TRUE)
# write(unmarked_selected_sequence,file = paste('./to_blast/',saving_file_name_unmarked,sep = ''),append = TRUE)

#writing unmarked multifasta for all samples

sample_list<-colnames(my_seq_df_with_loc)
sample_list<-sample_list[2:length(sample_list)]
for (sample in sample_list) {
  current_selected_sample<-sample
  first_line<-paste(">",current_selected_sample,sep = '')
  current_sequence<-my_seq_df_with_loc[[current_selected_sample]][extract_begin:extract_end]
  current_sequence_to_write<-paste(current_sequence,collapse = '')
  
  
  #for unmarked sequence file
  # remove markings for pure sequence
  
  unmarked_current_sequence_to_write<-stringr::str_replace_all(current_sequence_to_write, '\\*', '')
  unmarked_current_sequence_to_write<-stringr::str_replace_all(unmarked_current_sequence_to_write, '\\(', '')
  unmarked_current_sequence_to_write<-stringr::str_replace_all(unmarked_current_sequence_to_write, '\\)', '')
  
  #write unmarked files
  write(first_line,file = paste('./unmarked_sequence_files/',saving_file_name_unmarked,sep = ''),append = TRUE)
  write(unmarked_current_sequence_to_write,file = paste('./unmarked_sequence_files/',saving_file_name_unmarked,sep = ''),append = TRUE)
}

locations_left<-subset(locations_left,locations_left>extract_end)

}

```
From resulting files,
copied all the unmarked files to computecanada and then made a new directory to blast
```
mkdir unmarked_sequence_files_blasted
```
Then
```
module load StdEnv/2023  gcc/12.3 blast+/2.14.1
```

 if you did not create reference db already

```
makeblastdb -in ../../new_project_Apr_2023/new_data_Jun23/reference_genome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -title reference -dbtype nucl -out ref.fa
```

Then Blast files 

```
for i in unmarked_sequence_files/* ;do blastn -query $i -db reference/ref.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq qseq" | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > unmarked_sequence_files_blasted/${i##*/}_blasted;done

```
Then inside the output directory
```
mkdir to_download
mv Sex_sp_locations* to_download/
```

Downloaded
```
rsync -axvH --no-g --no-p premacht@beluga.computecanada.ca:/scratch/premacht/8_mb_region_2024/BLAST/unmarked_sequence_files_blasted/to_download .
```
Then I copied this to_download folder and unmarked_sequence_files folder to a different directory.
Following R script does the final analysis giving you the exact locations of sex fixed SNPs 
```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))
library(seqinr)
library(RFLPtools)
library(insight)
require(plyr)
library(plyr)
install.packages(rowr)

#cbind.fill function
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("msa")

file_list<-list.files("./to_download/")
sequence_data <- data.frame(matrix(ncol = 14, nrow = 0))

#I am getting rid of empty files/not blasted here 

## Get vector of all file names
ff <- dir("./to_download", recursive=TRUE, full.names=TRUE)
## Extract vector of empty files' names
eff <- ff[file.info(ff)[["size"]]==0]
#make a directory for files with too much deletions
dir.create("./unmatched_files_to_check_manually")
file.copy(eff,"./unmatched_files_to_check_manually")
## Remove empty files
unlink(eff, recursive=TRUE, force=FALSE)
print(paste("moved _",eff,"because that file was empty"))

#create file to store sex fixed locations
male_fixed_locs <- data.frame(matrix(ncol = 24, nrow = 0))
female_fixed_locs <- data.frame(matrix(ncol = 24, nrow = 0))


for (file_name in file_list) {
  working_file<-paste("./to_download/",file_name,sep = '')
  print(working_file)

my_data<-read.table(working_file)
colnames(my_data)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","sseq","qseq")

#neme sex fixed dfs properly
other_cols<-c("sloc","qloc","sseq","qseq")
sample_name_columns<-my_data$qseqid
final_col<-c("file")
total_col_names<-append(other_cols,sample_name_columns)
total_col_names<-append(total_col_names,final_col)
colnames(male_fixed_locs)<-total_col_names
colnames(female_fixed_locs)<-total_col_names
# # ******************************************************

# remove rows with different lengths before comparing sequences

# Factor_1<-factor(my_data$length)
# df1<-data.frame(Factor_1)
# common_length<-names(which.max(table(df1$Factor_1)))
# 
# removed_samples<-my_data[my_data$length != common_length, ]
# 
# if (length(removed_samples!=0)) {
# 
#   print_color(paste(removed_samples$qseqid,"was removed from",file_name,"analysis due to length difference"),'red')
#   Sys.sleep(2)
# }
# 
# 
# my_data<-my_data[my_data$length == common_length, ]

#*******This converts blast reads into fasta like read df
#Create empty df

sequence_data <- data.frame(matrix(ncol = 4, nrow = 0))
x <- c("sloc", "qloc", "sseq",'qseq')
colnames(sequence_data) <- x

# This does not need a loop. Just left it coz lazy
# Use a number of a sample you trust better as end number
#select sample no from here
my_data$qseqid
sample_no<-2

for (row_no in 1:sample_no) {
  print(row_no)
  
  selected_sample<-my_data[1,]
  
  #make sequence location list
  s_location_range<-seq(my_data$sstart[row_no],my_data$sstart[row_no]+my_data$length[row_no]-1)
  
  #make query location list
  q_location_range<-seq(my_data$qstart[row_no],my_data$qstart[row_no]+my_data$length[row_no]-1)
  
  #make sequence bp list
  sseq_long<-my_data$sseq[row_no]
  sbp_list<-strsplit(sseq_long,split = '')
  
  #make query bp list
  qseq_long<-my_data$qseq[row_no]
  qbp_list<-strsplit(qseq_long,split = '')
  
  #make df with all the data
  s_loc_col<-data.frame(unlist(s_location_range))
  q_loc_col<-data.frame(unlist(q_location_range))
  sbp_loc_col<-data.frame(unlist(sbp_list))
  qbp_loc_col<-data.frame(unlist(qbp_list))
  #********
  
  
  full_info_ref_df<-data.frame(unlist(s_location_range),unlist(q_location_range),unlist(sbp_list),unlist(qbp_list))
  
  #combine dfs
  sequence_data<-rbind(sequence_data,full_info_ref_df)
}


colnames(full_info_ref_df) <- x


# now adding other samples as columns
for (j in 1:nrow(my_data)) {
  
  selected_sample<-my_data[j,]
  
  current_seq<-selected_sample$qseq
  
  #make sequence location list
  s_location_range<-seq(selected_sample$sstart[1],selected_sample$sstart[1]+selected_sample$length[1]-1)
  
  #make query location list
  q_location_range<-seq(selected_sample$qstart[1],selected_sample$qstart[1]+selected_sample$length[1]-1)
  
  #make sequence bp list
  sseq_long<-selected_sample$sseq[1]
  sbp_list<-strsplit(sseq_long,split = '')
  
  #make query bp list
  qseq_long<-selected_sample$qseq[1]
  qbp_list<-strsplit(qseq_long,split = '')
  qbp_df<-as.data.frame(qbp_list)
  current_sample_name<-selected_sample[1]
  colnames(qbp_df)<-current_sample_name
  
  # make df matching s_loc and qbp
  
  #matching names
  sloc<-s_location_range
  
  qbp_with_locs<-cbind(sloc,qbp_df)
  
  
  
  full_info_ref_df<-merge(full_info_ref_df,qbp_with_locs,by='sloc', all.x=TRUE, all.y = FALSE)
  
} 
# # ******************************************************Now Compare it******************************

#moving full_info_ref_df to  my_seq_df to make it easy

my_seq_df<-full_info_ref_df

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

#drop nonsense columns ue to name issue

#male_fixed_locs <- male_fixed_locs[,colSums(is.na(male_fixed_locs))<nrow(male_fixed_locs)]
#female_fixed_locs <- female_fixed_locs[,colSums(is.na(female_fixed_locs))<nrow(female_fixed_locs)]

#add file name to the end to make it easier to find
#adding new column only on first round

male_fixed_no_niger<-cbind.fill(male_fixed_no_niger,file=file_name)
female_fixed_no_niger<-cbind.fill(female_fixed_no_niger,file=file_name)

male_fixed_no_niger<-as.data.frame(male_fixed_no_niger)
female_fixed_no_niger<-as.data.frame(female_fixed_no_niger)

col_no_for_file<-length(male_fixed_no_niger)
colnames(male_fixed_no_niger)[col_no_for_file]<-"file"
colnames(female_fixed_no_niger)[col_no_for_file]<-"file"

  # col_no_for_file<-length(male_fixed_no_niger)
  # male_fixed_no_niger<-cbind.fill(male_fixed_no_niger,file_name)
  # colnames(male_fixed_no_niger)[col_no_for_file]<-"file"
  # male_fixed_no_niger<-as.data.frame(male_fixed_no_niger)
  # 
  # col_no_for_file<-length(female_fixed_no_niger)
  # female_fixed_no_niger<-cbind.fill(female_fixed_no_niger,file_name)
  # colnames(female_fixed_no_niger)[col_no_for_file]<-"file"
  # female_fixed_no_niger<-as.data.frame(female_fixed_no_niger)


  # male_fixed_no_niger$file<-file_name
  # female_fixed_no_niger$file<-file_name


male_fixed_locs<-rbind.fill(male_fixed_locs,male_fixed_no_niger)
female_fixed_locs<-rbind.fill(female_fixed_locs,female_fixed_no_niger)

}

#drop nonsense columns ue to name issue
# 
# male_fixed_locs <- male_fixed_locs[,colSums(is.na(male_fixed_locs))<nrow(male_fixed_locs)]
# female_fixed_locs <- female_fixed_locs[,colSums(is.na(female_fixed_locs))<nrow(female_fixed_locs)]

write.table(male_fixed_locs, file='./male_fixed_locs.txt', quote=FALSE, sep='\t', col.names = NA)
write.table(female_fixed_locs, file='./female_fixed_locs.txt', quote=FALSE, sep='\t', col.names = NA)

#For more cleaner look, subsetting only needed columns

cleaned_male_fixed_locs<-data.frame(male_fixed_locs$sloc,male_fixed_locs$qloc,male_fixed_locs$sseq,male_fixed_locs$F_Ghana_WZ_BJE4687_combined__sorted.bam,male_fixed_locs$F_SierraLeone_AMNH17272_combined__sorted.bam,male_fixed_locs$F_SierraLeone_AMNH17274_combined__sorted.bam,male_fixed_locs$F_IvoryCoast_xen228_combined__sorted.bam,male_fixed_locs$M_Ghana_WY_BJE4362_combined__sorted.bam,male_fixed_locs$M_Ghana_ZY_BJE4360_combined__sorted.bam,male_fixed_locs$M_SierraLeone_AMNH17271_combined__sorted.bam,male_fixed_locs$M_SierraLeone_AMNH17273_combined__sorted.bam,male_fixed_locs$file)
colnames(cleaned_male_fixed_locs)<-c("sloc","qloc","reference","F_Ghana_WZ_BJE4687_combined__sorted.bam","F_SierraLeone_AMNH17272_combined__sorted.bam","F_SierraLeone_AMNH17274_combined__sorted.bam","F_IvoryCoast_xen228_combined__sorted.bam","M_Ghana_WY_BJE4362_combined__sorted.bam","M_Ghana_ZY_BJE4360_combined__sorted.bam","M_SierraLeone_AMNH17271_combined__sorted.bam","M_SierraLeone_AMNH17273_combined__sorted.bam","file")
cleaned_male_fixed_locs<-cleaned_male_fixed_locs[order(cleaned_male_fixed_locs$sloc),]
write.table(cleaned_male_fixed_locs, file='./cleaned_male_fixed_locs.txt', quote=FALSE, sep='\t', col.names = NA)

cleaned_female_fixed_locs<-data.frame(female_fixed_locs$sloc,female_fixed_locs$qloc,female_fixed_locs$sseq,female_fixed_locs$F_Ghana_WZ_BJE4687_combined__sorted.bam,female_fixed_locs$F_SierraLeone_AMNH17272_combined__sorted.bam,female_fixed_locs$F_SierraLeone_AMNH17274_combined__sorted.bam,female_fixed_locs$F_IvoryCoast_xen228_combined__sorted.bam,female_fixed_locs$M_Ghana_WY_BJE4362_combined__sorted.bam,female_fixed_locs$M_Ghana_ZY_BJE4360_combined__sorted.bam,female_fixed_locs$M_SierraLeone_AMNH17271_combined__sorted.bam,female_fixed_locs$M_SierraLeone_AMNH17273_combined__sorted.bam,female_fixed_locs$file)
colnames(cleaned_female_fixed_locs)<-c("sloc","qloc","reference","F_Ghana_WZ_BJE4687_combined__sorted.bam","F_SierraLeone_AMNH17272_combined__sorted.bam","F_SierraLeone_AMNH17274_combined__sorted.bam","F_IvoryCoast_xen228_combined__sorted.bam","M_Ghana_WY_BJE4362_combined__sorted.bam","M_Ghana_ZY_BJE4360_combined__sorted.bam","M_SierraLeone_AMNH17271_combined__sorted.bam","M_SierraLeone_AMNH17273_combined__sorted.bam","file")
cleaned_female_fixed_locs<-cleaned_female_fixed_locs[order(cleaned_female_fixed_locs$sloc),]
write.table(cleaned_female_fixed_locs, file='./cleaned_female_fixed_locs.txt', quote=FALSE, sep='\t', col.names = NA)
```



