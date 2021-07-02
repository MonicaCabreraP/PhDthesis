#!/bin/bash

# #path to the scripts where the commands are located
cd /gpfs/projects/bsc08/bsc08471/p53/HiCUP/downsampled

# for sample in $SAMPLES
# do
	command=downsampled.cmd
	{

	echo "samtools view -s 103.642 -b /gpfs/projects/bsc08/bsc08471/p53/HiCUP/MOCK2KO/lane3_NoIndex_L003_R2_1.hicup.captured.bam > MOCK2KO_downsampled.bam"
	echo "samtools view -s 103.837 -b /gpfs/projects/bsc08/bsc08471/p53/HiCUP/MOCK1WT/lane3_NoIndex_L003_R1_2.hicup.captured.bam > MOCK1WT_downsampled.bam"
	echo "samtools view -s 103.610 -b /gpfs/projects/bsc08/bsc08471/p53/HiCUP/MOCK2WT/lane1_NoIndex_L001_R2_1.hicup.captured.bam > MOCK2WT_downsampled.bam"

	echo "samtools view -s 103.503 -b /gpfs/projects/bsc08/bsc08471/p53/HiCUP/NUT1KO/lane3221_HCT116-_-_Nutlin3a_1BR-101.1_NoIndex_L004_R2_1.hicup.captured.bam > NUT1KO_downsampled.bam"
	echo "samtools view -s 103.520 -b /gpfs/projects/bsc08/bsc08471/p53/HiCUP/NUT2KO/lane3223_HCT116-_-_Nutlin3a_2BR-103.1_NoIndex_L006_R1_2.hicup.captured.bam > NUT2KO_downsampled.bam"
	echo "samtools view -s 103.566 -b /gpfs/projects/bsc08/bsc08471/p53/HiCUP/NUT1WT/lane3220_HCT116_Nutlin3a_1BR-100.1_NoIndex_L003_R1_2.hicup.captured.bam > NUT1WT_downsampled.bam"
	echo "samtools view -s 103.653 -b /gpfs/projects/bsc08/bsc08471/p53/HiCUP/NUT2WT/lane3222_HCT116_Nutlin3a_2BR-102.1_NoIndex_L005_R2_1.hicup.captured.bam > NUT2WT_downsampled.bam"

	} >> $command

#done

job=downsampled.job
	{
	echo "#!/bin/bash"

	#Job directives
	################
	#Name of the job that will apear in the queue
	echo "#SBATCH --job-name=downsampled"

	#where the job will run
	echo "#SBATCH --workdir=/gpfs/projects/bsc08/bsc08471/p53/HiCUP/downsampled"

	#File to collect the standard error outputs (stderr) of the job
	echo "#SBATCH --error=/gpfs/projects/bsc08/bsc08471/p53/HiCUP/downsampled/logs/downsampled_%j.err"

	#File to collect the standard output (stdout) of the job
	echo "#SBATCH --output=/gpfs/projects/bsc08/bsc08471/p53/HiCUP/downsampled/logs/downsampled_%j.out"

	#The number of requested nodes
	##echo "#SBATCH --nodes=NUMBER"

	#The number of processes to start
	echo "#SBATCH --ntasks=7"

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
	echo "#SBATCH --constraint=highmem"

	#Enable e-mail notifications when a job starts (begin) and (all)/or ends (end) or none - any e-mail
	echo "#SBATCH --mail-type=all"
	echo "#SBATCH --mail-user=monica.cabrera@bsc.es"

	echo "module load samtools"

	echo "/apps/GREASY/latest/INTEL/IMPI/bin/greasy /gpfs/projects/bsc08/bsc08471/p53/HiCUP/downsampled/$command"

	} > $job

chmod +x $command $job

sbatch $job


# Mónica Cabrera Pasadas