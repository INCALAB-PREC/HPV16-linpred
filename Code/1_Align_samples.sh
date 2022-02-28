#!/bin/bash
#########################################################
###                                                   ###
###       ALIGN SAMPLES ON THE HPV16 REF GENOME       ###
###                    with mafft                     ###
###               lasensio@idibell.cat                ###
###                                                   ###
#########################################################

## HOW TO INSTALL MAFFT?
# https://mafft.cbrc.jp/alignment/software/

## DESCRIPTION: 
# This is a bash code to align the samples of study over the HPV16 reference genome.

## RECOMENDED FOLDER DISTRIBUTION: 
#- PATH
#   - Data
#        - HPV16_REFERENCE.fasta (available the github repository or PAVE webpage)
#        - Samples_to_align.fasta
#   - Aligned (NEW - STEP 1) 
#        - Samples.fasta (NEW - STEP 1)
#        - tmp_samples (Temporary folder - STEP 2) 
#        - tmp_aligned (Temporary folder - STEP 3)
#        - Aligned_samples.fasta (NEW - STEP 3)

## We will ALWAYS use the HPV16_A1 reference genome to keep all positions constant! 
# gi|333031|lcl|HPV16REF.1|HPV16REF or GenBank Accession code: K02718.1

#########################################################

# STEP 1 - OPTIONAL: Concatenate files to be aligned (to align all the samples in a single time if they provide from different fasta files) 

cd PATH
mkdir Aligned
cd Aligned

# Copy the samples to align inside the Aligned folder
cp PATH/HPV16_10_PAVE_REFERENCE.fasta samples1.fasta
cp PATH/HPV16_588_NCBI.fasta samples2.fasta

# Concatenate the samples in a single file
cat samples1.fasta samples2.fasta > Samples.fasta

rm samples1.fasta
rm samples2.fasta

#########################################################

# STEP 2: Create single and temporary files for each HPV16 genome

mkdir tmp_samples/
cd tmp_samples
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        name=${line#>}
        outfile=${name%%.*}.fasta
        echo $outfile
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < ../Samples.fasta

#########################################################

# STEP 3: Align each sample (fasta file) on the reference set with mafft.

## Conditions and parameters used for mafft:
# --add: adding unaligned full-length sequence(s) into an existing alignment
# --keeplength: the alignment length is unchanged.  Insertions at the new sequences are deleted.  

mkdir ../tmp_aligned/
for file in *.fasta
do
  echo "Processing $file file..." 
  mafft --add $file --keeplength ../../Data/HPV16_REFERENCE.fasta > ../tmp_aligned/$file
  name=$(echo $file | sed "s/.fasta//")
  echo $name
  echo $name
  sed -n '/'"$name"'/,$p'  ../tmp_aligned/$file >> ../Aligned_samples.fasta
  rm ../tmp_aligned/$file
done

#All aligned samples are appended in a single fasta file called Aligned_samples.fasta inside the Aligned folder. 

#########################################################

# STEP 4 - OPTIONAL: Remove temporary files that won't be useful anymore and they waste a lot of computer memory. 

cd ../
rm Samples.fasta
rm -r tmp_aligned
rm -r tmp_samples


#####   FINISH    #####
