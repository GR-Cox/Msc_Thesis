## Usearch/Vsearch Pipeline for 18S General Eukaryote Primers

## This is an alternative to Ian Dickie's protocol adapated for 18S micro-eukaryote data by George Cox.
## Designed to work on a UC RCC virtual machine.

## If this is quite new to you I highly reccomend reading Pipeline_set_up.txt as it shows you how to set up the virtual machine/computer so that this runs smoothly.
## That file also has code for how to set up the blast database and how to create a conda environment. 

## Make sure your code is updated with your project names, file paths, database, primers, and sequence lengths before running the following commands. 

##-------------------------------------------------------------------------------------------------------
## 1st step: Transfer files and updated code using Rsync (Change folder names to suit your VM).  

rsync -r /your/files/code rccuser@ip.address:~/code

rsync -r /your/files/data rccuser@ip.address:~/data

##-------------------------------------------------------------------------------------------------------
## 2nd Step: unzipping files

ssh rccuser@ip.address

cd data

gunzip *.gz

##-------------------------------------------------------------------------------------------------------
## 3rd Step: Cleaning up file names

## Remove the _001 from names
for file in *; do [ -f "$file" ] && ( mv "$file" "$(echo $file | sed -e 's/_001//g')" ); done

for file in *; do [ -f "$file" ] && ( mv "$file" "$(echo $file | sed -e 's/forward/R1/g')" ); done

## remove the underscore, which causes sequencing naming problems later, also removed "processed" for simplicity
for file in *; do [ -f "$file" ] && ( mv "$file" "$(echo $file | sed -e 's/processed_//g')" ); done

##-------------------------------------------------------------------------------------------------------
## 4th Step: Turn on conda enviroment

conda activate cutenv

##-------------------------------------------------------------------------------------------------------
## 5th Step: Running the pipeline

#Rscript ~/code/Usearch_pipeline3.R

#Or if you are working on a large dataset that may take a long time to run use this code

nohup Rscript ~/code/Usearch_pipeline3.R > ~/pipeline_stdout.txt

# Inspect the output for errors and make sure it looks alright
# can view the run settings via logStats.log file
nano logStats.log

##-------------------------------------------------------------------------------------------------------
## 6th Step: Blasting the pipeline outputs (make sure you have set up your blast database).

#Rscript ~/code/postUsearch_Blast3.R

nohup Rscript ~/code/postUsearch_Blast3.R > ~/blast_stdout.txt

##-------------------------------------------------------------------------------------------------------
## 7th Step: Data Compliation/Formatting

Rscript ~/code/DataCompilation3.R

##-------------------------------------------------------------------------------------------------------
## 8th Step: Export outputs and begin community analysis

rsync rccuser@ip.address:~/data/ProjectName_OtuOutput_tabulated.Rdat /your/files/outputs
