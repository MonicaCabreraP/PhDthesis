#!/bin/bash

#__SYNOPSIS__
# HiCUP.sh is used to run the whole HiCUP pipeline and send the job to Mare Nostrum 4
# HiCUP.sh [PATH FOR THE HiCUP RESULTS] [PATH OF THE DIGEST FILE] [PATH OF THE INDEX FILE]

# Activate in case of no permision to give permisions
#chmod +x -R /gpfs/scratch/bsc08/bsc08471/data

# Command to prompt in terminal
if [ $# -lt 3 ]
then
	echo "Usage: HiCUP.sh <path-HiCUP_results-folder> <path-digest> <path-index>"
	exit
fi

results=$1

#############################
# EDIT START #
#############################
# Fetching the sample names
samples = ""

fastq1= ""
fastq2=""

#############################
# EDIT END #
#############################

#rm $1HiCUP.cmd

for sample in $samples
do
	mkdir $1$sample
	
	#### Creating the configuration file for each sample ####
	configs=$1$sample/HiCUP_config_${sample}.txt
	{
	#Directory to which output files should be written (optional parameter)
	#Set to current working directory by default"
	echo "Outdir: $1$sample"

	#Number of threads to use
	echo "Threads: 10"

	#Suppress progress updates (0: off, 1: on)
	echo "Quiet:0"

	#Retain intermediate pipeline files (0: off, 1: on)
	echo "Keep:0"

	#Compress outputfiles (0: off, 1: on)"
	echo "Zip:1"

	#Path to the alignment program (Bowtie or Bowtie2)
	#Bowtie when using Bowtie indices, or Bowtie2 when using Bowtie2 indices.
	echo "Bowtie: /apps/BOWTIE/1.2.2/GCC/bowtie"
	#echo "Bowtie2: /apps/BOWTIE2/2.3.2/INTEL/IMPI/bin"

	#Path to R
	echo "R: /apps/R/3.5.1/INTEL/bin/R"

	#Path to the reference genome indices
	##For bowtie
	#Remember to include the basename of the genome indices
	echo "Index: /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/human/GRCh37/current/indices/ebwt/Homo_sapiens_GRCh37_92_index"
	##For bowtie2
	#echo "Index: /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/human/GRCh37/current/indices/bt2/Homo_sapiens_GRCh37_92_index"

	#Path to the genome digest file produced by hicup_digester
	echo "Digest: /gpfs/projects/bsc08/shared_projects/IJC_3Dchromatin/genomes/human/GRCh37/current/digest_genome/Digest_Human_GRCh37_92_HindIII_None_10-51-28_27-03-2020.txt.gz"

	#FASTQ format (valid formats: 'Sanger', 'Solexa_Illumina_1.0','Illumina_1.3' or 'Illumina_1.5'). 
	# If not specified, HiCUP will #try to determine the format automatically by analysing one of
	# the FASTQ files. All input FASTQ will assumed to be in that format.
	echo "Format: "

	#Maximum di-tag length (optional parameter)
	echo "Longest: 800"

	#Minimum di-tag length (optional parameter)
	echo "Shortest: 150"

	#FASTQ files to be analysed, placing paired files on adjacent lines
	$fastq1
	$fastq2

	} > $configs

	#### Creating the command containing all the configuration files of each sample ####
	echo "hicup --config $1${sample}/HiCUP_config_${sample}.txt" >> $1HiCUP.cmd

done

chmod +x -R $1

#### Creating the job command and sending it to the queue ####
mkdir $1logs

command=$1HiCUP.cmd
job=$1HiCUP.job
	{
	echo "#!/bin/bash"

	#Job directives
	################
	#Name of the job that will apear in the queue
	echo "#SBATCH --job-name=HiCUP"

	#where the job will run
	echo "#SBATCH --workdir=$1"

	#File to collect the standard error outputs (stderr) of the job
	echo "#SBATCH --error=$1logs/HiCUP_%j.err"

	#File to collect the standard output (stdout) of the job
	echo "#SBATCH --output=$1logs/HiCUP_%j.out"

	#The number of requested nodes
	##echo "#SBATCH --nodes=NUMBER"

	#The number of processes to start
	echo "#SBATCH --ntasks=$(wc -l $command | awk '{print $1}')"

	#Optional: specify how many threads each process would open
	#The number of cores assigned to the job will be the total _tasks number*cpus_per_task number
	echo "#SBATCH --cpus-per-task=4"

	#Number of tasks assigned to a node
	##echo "#SBATCH --task-per-node=NUMBER"

	#Limit of wall clok time. Job will be killed after the time has passed- value greater than the real executation time 
	#echo "#SBATCH --time=02:00:00"

	#Debug queue for small test
	#echo "#SBATCH --qos=debug"

	#Limit of wall clok time. Job will be killed after the time has passed- value greater than the real executation time 
	echo "#SBATCH --time=15:00:00"

	#HighMem node with 7928 MB of RAM per core- without this directive the jobs will be sent to standard nodes with 1880 MB of RAM per core.
	##There are limited number of high memory nodes available (216 nodes- 10368 cores) out of 3456 nodes-16588 cores in total. Every core comes with 2GB RAM
	echo "#SBATCH --constraint=highmem"

	#Enable e-mail notifications when a job starts (begin) and (all)/or ends (end) or none - any e-mail
	echo "#SBATCH --mail-type=all"
	echo "#SBATCH --mail-user=monica.cabrera@bsc.es"

	echo "module load bowtie2"
	echo "module load bowtie"
	echo "module load hicup"
	echo "module load samtools"

	echo "/apps/GREASY/latest/INTEL/IMPI/bin/greasy $command"

	} > $job

chmod +x $command $job

sbatch $job


# Mónica Cabrera Pasadas