#!/bin/bash

#BAITS=$1
#BAMS=2

grep -v '#' /gpfs/home/bsc08/bsc08471/data/0.sample_names.txt > /gpfs/scratch/bsc08/bsc08471/tmp/data.txt
SAMPLES=$(</gpfs/scratch/bsc08/bsc08471/tmp/data.txt)

#COUNTER=0
for sample in $SAMPLES
do
	command=/gpfs/projects/bsc08/bsc08471/p53/HiCUP/newRun/HiCUP_capture_efficiency.cmd
	{
	grep -v '#' /gpfs/home/bsc08/bsc08471/data/2.bams.txt | grep -- "${sample}" > /gpfs/scratch/bsc08/bsc08471/tmp/bam_${sample}.txt
	BAM=$(cat /gpfs/scratch/bsc08/bsc08471/tmp/bam_$sample.txt)
	echo "perl /gpfs/home/bsc08/bsc08471/scripts/get_captured_reads --baits /gpfs/home/bsc08/bsc08471/HindIII/baits_coordinates_Human_GRCh37_92.txt $BAM" 
	} >> $command
done

chmod +x -R $command

mkdir /gpfs/projects/bsc08/bsc08471/p53/HiCUP/newRun/logs

command=/gpfs/projects/bsc08/bsc08471/p53/HiCUP/newRun/HiCUP_capture_efficiency.cmd
job=/gpfs/projects/bsc08/bsc08471/p53/HiCUP/newRun/capture_effieciency.job
	{
	echo "#!/bin/bash"

	#Job directives
	################
	#Name of the job that will apear in the queue
	echo "#SBATCH --job-name=capture_effieciency"

	#where the job will run
	echo "#SBATCH --workdir=."

	#File to collect the standard error outputs (stderr) of the job
	echo "#SBATCH --error=/gpfs/projects/bsc08/bsc08471/p53/HiCUP/newRun/logs/capture_effieciency_%j.err"

	#File to collect the standard output (stdout) of the job
	echo "#SBATCH --output=/gpfs/projects/bsc08/bsc08471/p53/HiCUP/newRun/logs/capture_effieciency_%j.out"

	#The number of requested nodes
	##echo "#SBATCH --nodes=NUMBER"

	#The number of processes to start
	echo "#SBATCH --ntasks=12"

	#Optional: specify how many threads each process would open
	#The number of cores assigned to the job will be the total _tasks number*cpus_per_task number
	echo "#SBATCH --cpus-per-task=4"

	#Number of tasks assigned to a node
	##echo "#SBATCH --task-per-node=NUMBER"

	#Limit of wall clok time. Job will be killed after the time has passed- value greater than the real executation time 
	echo "#SBATCH --time=02:00:00"

	#Debug queue for small test
	echo "#SBATCH --qos=debug"

	#Limit of wall clok time. Job will be killed after the time has passed- value greater than the real executation time 
	#echo "#SBATCH --time=15:00:00"

	#HighMem node with 7928 MB of RAM per core- without this directive the jobs will be sent to standard nodes with 1880 MB of RAM per core.
	##There are limited number of high memory nodes available (216 nodes- 10368 cores) out of 3456 nodes-16588 cores in total. Every core comes with 2GB RAM
	#echo "#SBATCH --constraint=highmem"

	#Enable e-mail notifications when a job starts (begin) and (all)/or ends (end) or none - any e-mail
	#echo "#SBATCH --mail-type=all"
	#echo "#SBATCH --mail-user=monica.cabrera@bsc.es"

	echo "module load samtools"
	echo "module load perl"
	echo "module load R/3.6.1"

	echo "/apps/GREASY/latest/INTEL/IMPI/bin/greasy $command"

	} > $job

chmod +x $command $job

sbatch $job


##########################
# Mónica Cabrera Pasadas
##########################


