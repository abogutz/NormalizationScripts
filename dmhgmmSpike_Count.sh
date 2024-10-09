#! /bin/bash
#SBATCH --account=def-mlorincz       		   # required
#SBATCH --cpus-per-task=8                        # number of MPI processes
#SBATCH --mem-per-cpu=4G                 # memory; default unit is megabytes
#SBATCH --time=00-20:00                   # time (DD-HH:MM)
#SBATCH --mail-user=aaron.bogutz@ubc.ca
#SBATCH --mail-type=ALL


# Drosophila + Human Spike-in counts for Normalization

# Setup for Graham ( but can be adapted for other servers)

module load samtools

THREADS=$SLURM_CPUS_PER_TASK
MEM=4G
BIN_SIZE=1
SMOOTH_WIN=0
MAPQ=5

TEMP_DIR=$SCRATCH"/"$SLURM_JOB_ID"/"
TEMP="$TEMP_DIR/temp"
HEADER_DM="header_dm.sam"
HEADER_HG="header_hg.sam"
HEADER_MM="header_mm.sam"
VIEW="samtools view -@ $THREADS "
BAM_COVERAGE="/project/def-mlorincz/scripts/utilities/miniconda3/bin/bamCoverage"

echo -e "Name\tDm_Reads\tDm_Reads_MapQ\tDm_Reads_MapQ_NoDup\tHg_Reads\tHg_Reads_MapQ\tHg_Reads_MapQ_NoDup\tMm_Reads\tMm_Reads_MapQ\tMm_Reads_MapQ_NoDup\tUnmapped" >> spike_stats.txt

mkdir $TEMP_DIR

for FILE in *.bam
do
	DM_SAM=$TEMP_DIR/${FILE//.bam/_dm6.sam}
	DM_BAM=${FILE//.bam/_dm6.bam}
	HG_SAM=$TEMP_DIR/${FILE//.bam/_hg38.sam}
	HG_BAM=${FILE//.bam/_hg38.bam}
	MM_SAM=$TEMP_DIR/${FILE//.bam/_mm10.sam}
	MM_BAM=${FILE//.bam/_mm10.bam}
	$VIEW $FILE | awk -v dm=$DM_SAM -v hg=$HG_SAM -v mm=$MM_SAM '{
		if ($3 ~ /dm/) {
			print > dm;
		} else if ($3 ~ /hg/) {
			print > hg;
		} else {
			print > mm;
		}
	}'
	$VIEW -H $FILE | awk -v dm=$HEADER_DM -v hg=$HEADER_HG -v mm=$HEADER_MM '{
		if ($2 ~ /^SN/) {
			if ($2 ~ /dm/) {
				print > dm;
			} else if ($2 ~ /hg/) {
				print > hg;
			} else {
				print > mm;
			}
		} else {
			print > dm;
			print > hg;
			print > mm;
		}
	}'
	cat $HEADER_DM $DM_SAM | sed 's/dm_//g' > $TEMP
	mv $TEMP $DM_SAM
	cat $HEADER_HG $HG_SAM | sed 's/hg_//g' > $TEMP
	mv $TEMP $HG_SAM
	cat $HEADER_MM $MM_SAM > $TEMP
	mv $TEMP $MM_SAM
	rm $HEADER_DM $HEADER_MM $HEADER_HG
	$VIEW -bh -o $DM_BAM $DM_SAM
	samtools index $DM_BAM
	$VIEW -bh -o $HG_BAM $HG_SAM
	samtools index $HG_BAM
	$VIEW -bh -o $MM_BAM $MM_SAM
	samtools index $MM_BAM
	rm $DM_SAM $MM_SAM $HG_SAM
	DM=`$VIEW -c $DM_BAM`
	DM_MQ=`$VIEW -c -q $MAPQ $DM_BAM`
	DM_MQDP=`$VIEW -c -q $MAPQ -F 1024 $DM_BAM`
	
	HG=`$VIEW -c $HG_BAM`
	HG_MQ=`$VIEW -c -q $MAPQ $HG_BAM`
	HG_MQDP=`$VIEW -c -q $MAPQ -F 1024 $HG_BAM`

	MM=`$VIEW -c -F 4 $MM_BAM`
	MM_MQ=`$VIEW -c -q $MAPQ $MM_BAM`
	MM_MQDP=`$VIEW -c -q $MAPQ -F 1024 $MM_BAM`
	UM=`$VIEW -c -f 4 $MM_BAM`
	echo -e "$FILE\t$DM\t$DM_MQ\t$DM_MQDP\t$HG\t$HG_MQ\t$HG_MQDP\t$MM\t$MM_MQ\t$MM_MQDP\t$UM" >> spike_stats.txt
	$BAM_COVERAGE --binSize 1 -p $THREADS --normalizeUsing CPM --smoothLength 0 --outFileFormat bigwig --minMappingQuality 5 --ignoreDuplicates -b $DM_BAM --outFileName ${DM_BAM//.bam/.bw}
	$BAM_COVERAGE --binSize 1 -p $THREADS --normalizeUsing CPM --smoothLength 0 --outFileFormat bigwig --minMappingQuality 5 --ignoreDuplicates -b $HG_BAM --outFileName ${HG_BAM//.bam/.bw}
	$BAM_COVERAGE --binSize 1 -p $THREADS --normalizeUsing CPM --smoothLength 0 --outFileFormat bigwig --minMappingQuality 5 --ignoreDuplicates -b $MM_BAM --outFileName ${MM_BAM//.bam/.bw}
done


rm -r $TEMP_DIR