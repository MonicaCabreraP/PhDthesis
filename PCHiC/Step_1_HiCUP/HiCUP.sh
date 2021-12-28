#!/bin/bash

# Details of the script
############################################################################################################################################################################################################################################################################
#Script Name	: HiCUP.sh                                                                                          
#Description	: Mapping and processing Hi-C data and send the job to Mare Nostrum 4 (were HiCUP is installed)
#Requirements   : Perl, Bowtie2, R (tested with version 3.6.0), SAM tools (version 0.1.18 or later)
#Output         : A results folder will be created containing the bam files of the samples analyzed in the folder where you have this script                                                            
#Author       	: Monica Cabrera Pasadas                                                
#Email         	: monica.cabrera.pasadas@gmail.com     
#Last update    : 20 October 2021                                   
############################################################################################################################################################################################################################################################################

# Part to edit by the user to run HiCUP
#############################################################################################################################################################################################################################################################################
# EDIT START # Put your own paths and variables 
#############################################################################################################################################################################################################################################################################
# Fetching the paths to run the script:
samples="DMSO_WT1 DMSO_WT2 DMSO_KO1 DMSO_KO2 ETO_WT1 ETO_WT2 ETO_KO1 ETO_KO2 NUT_WT1 NUT_WT2 NUT_KO1 NUT_KO2" #sample names to analyze
path_digest_genome="/gpfs/projects/bsc08/bsc08471/p53/PCHiC/Genome/Homo_Sapiens/GRCh37/Digest/GRCh37" #path to the digested genome obtained from 
path_index_genome="/gpfs/projects/bsc08/bsc08471/p53/PCHiC/Genome/Homo_Sapiens/GRCh37/Digest/GRCh37/Digest_Human_GRCh37_92_HindIII_None_10-51-28_27-03-2020.txt.gz"
path_fastq="/gpfs/projects/bsc08/bsc08471/p53/PCHiC/Data" #path were your fastq files are located
path_R="/apps/R/3.5.1/INTEL/bin/R"
path_bowtie=""
path_bowtie2="/apps/BOWTIE2/2.3.2/INTEL/IMPI/bin"
# Article about HiCUP: https://pubmed.ncbi.nlm.nih.gov/26835000/        
# Full documentation: https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html                     

############################################################################################################################################################################################################################################################################
# EDIT END # After here you don't need to edit anything else more
############################################################################################################################################################################################################################################################################
mkdir ${PWD}/results
mkdir ${PWD}/results/configFiles

for sample in $samples
do
	mkdir ${PWD}/results/${sample}

	#### Creating the configuration file for each sample ####
	configs=${PWD}/results/configFiles/HiCUP_config_${sample}.txt
	{
	echo "Outdir: ${PWD}/results/${sample}" #Directory to which output files should be written (optional parameter), set to current working directory by default
	echo "Threads: 10" #Number of threads to use
	echo "Quiet:0" #Suppress progress updates (0: off, 1: on)
	echo "Keep:0" #Retain intermediate pipeline files (0: off, 1: on)
	echo "Zip:1" #Compress outputfiles (0: off, 1: on)"
	#echo "Bowtie: /apps/BOWTIE/1.2.2/GCC/bowtie" #Path to the bowtie alignment program if we use the bowtie indices
	echo "Bowtie2: ${path_bowtie2}" #Path to the alignment program bowtie2 if we use the bowtie2 indices
	echo "R: " # Path to R
	#echo "Index: ${path_index_genome}/Homo_sapiens_GRCh37_92_index" #Path bowtie reference genome indice
	echo "Index: /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/human/GRCh37/current/indices/bt2/Homo_sapiens_GRCh37_92_index" #Path bowtie reference genome indice
	echo "Digest: ${path_digest_genome}/Digest_Human_GRCh37_92_HindIII_None_10-51-28_27-03-2020.txt.gz" #Path to the genome digest file produced by hicup_digester

	#FASTQ format (valid formats: 'Sanger', 'Solexa_Illumina_1.0','Illumina_1.3' or 'Illumina_1.5'). 
	## If not specified, HiCUP will try to determine the format automatically by analysing one of the FASTQ files. 
	### All input FASTQ will assumed to be in that format.
	echo "Format: "

	echo "Longest: 800" # Maximum di-tag length (optional parameter)
	echo "Shortest: 150" # Minimum di-tag length (optional parameter)

	#FASTQ files to be analysed, placing paired files on adjacent lines
	FASTQ1=$(echo ${path_fastq}/${sample}/*R1.fastq.gz)
	FASTQ2=$(echo ${path_fastq}/${sample}/*R2.fastq.gz)
	echo "$FASTQ1"
	echo "$FASTQ2"

	} > $configs

	#### Creating the command containing all the configuration files of each sample to run in a single job ####
	echo "hicup --config ${PWD}/results/configFiles/HiCUP_config_${sample}.txt" >> ${PWD}/results/HiCUP.cmd

done

#### Creating the job command and sending it to the queue ####
mkdir ${PWD}/logs

command=${PWD}/results/HiCUP.cmd
job=${PWD}/results/HiCUP.job
	{
	echo "#!/bin/bash"

	#Job directives
	################
	echo "#SBATCH --job-name=HiCUP" #Name of the job that will apear in the queue
	echo "#SBATCH --workdir=${PWD}" #where the job will run
	echo "#SBATCH --error=${PWD}/logs/HiCUP_%j.err" #File to collect the standard error outputs (stderr) of the job
	echo "#SBATCH --output=${PWD}/logs/HiCUP_%j.out" #File to collect the standard output (stdout) of the job
	#echo "#SBATCH --nodes=NUMBER" #The number of requested nodes
	echo "#SBATCH --ntasks=$(wc -l $command | awk '{print $1}')" #The number of processes to start
	echo "#SBATCH --cpus-per-task=4" #The number of cores assigned to the job will be the total _tasks number*cpus_per_task number
	#echo "#SBATCH --task-per-node=NUMBER" 	#Number of tasks assigned to a node
	#echo "#SBATCH --time=15:00:00" #Limit of wall clok time after the job will be killed after  
	#echo "#SBATCH --qos=debug" #Debug queue for small test
	echo "#SBATCH --constraint=highmem" #HighMem node with 7928 MB of RAM per core. Every core comes with 2GB RAM
	echo "#SBATCH --mail-type=all" #Enable e-mail notifications when a job starts (begin) and (all)/or ends (end) or none - any e-mail
	echo "#SBATCH --mail-user=mcabrera@bsc.es" #Email address

	# Required modules
	###################
	#echo "module load bowtie"
	echo "module load bowtie2"
	echo "module load hicup"
	echo "module load samtools"
	!!!!!! R 3.6.1

	echo "/apps/GREASY/latest/INTEL/IMPI/bin/greasy ${PWD}/results/$command"

	} > $job

chmod +x -R $command $job

#sbatch $job

#########################
# Mónica Cabrera Pasadas
#########################
