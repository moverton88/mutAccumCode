# Full pipeline for aligning reads and calling variants for a single sample


module load bowtie2
module load samtools
module load bcftools

# Reference to align reads
export REFPREFIX=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_rnm_UCSD_2020
export REFSEQ=${REFPREFIX}.fna
export REFTYPE=RM

# Set up dir variables
export DATADIR=/oasis/tscc/scratch/mioverto/data/MAseq1
export FQDIR=${DATADIR}/reads/trim
export BAMDIR=${DATADIR}/bam

# Renamed, trimmed reads
export INDEX=N_A05
export R1FILE=${DATADIR}/reads/trim/${INDEX}_1_R1P.trimmed.fastq
export R2FILE=${DATADIR}/reads/trim/${INDEX}_1_R2P.trimmed.fastq
export R1UFILE=${DATADIR}/reads/trim/${INDEX}_1_R1U.trimmed.fastq
export R2UFILE=${DATADIR}/reads/trim/${INDEX}_1_R2U.trimmed.fastq


# VCF file output
export VCFOUT=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/${INDEX}_${REFTYPE}.vcf.gz

# Store de-duplication metrics
export METRICS=/oasis/tscc/scratch/mioverto/data/MAseq1/bam/picard_metrics

# If variants need to be filtered by position before analysis
# export POSFILE=/oasis/tscc/scratch/mioverto/data/refseq/POS_files/RMxBY_variants.vcf

# Start with alignment of renamed, trimmed fastq files and un-indexed reference fasta file

# If reference has not been indexed yet
if [ ! -f ${REFPREFIX}.1.bt2 ]; then
    "Reference sequence needs to be indexed (bowtie2-build)"
    bowtie2-build ${REFSEQ} ${REFPREFIX}
fi

# Align paired and unpaired reads reads with bowtie2. 
# ${INDEX} used for naming alignments. Output is sorted bam alignment
# -X flag gives the maximum valid gap between paired reads
# samtools view -h [include header] -b [bam output] | sort -m [memory allocated] -o [output file] -T [temp files]
bowtie2 --rg-id ${INDEX} --rg SM:${INDEX} -X 1000 -x ${REFPREFIX} -1 ${R1FILE} -2 ${R2FILE} -U "${R1UFILE},${R2UFILE}" \
| samtools view -h -b | samtools sort -m 10000000 -o ${BAMDIR}/${INDEX}.bam -O bam -T ${BAMDIR}/${INDEX}

# Paths to run GATK
GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar
PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH

# Remove duplicate reads
java -jar $GATK MarkDuplicates \
	--INPUT ${BAMDIR}/${INDEX}.bam \
	--OUTPUT ${BAMDIR}/${INDEX}.dm.bam \
	--METRICS_FILE ${METRICS}/${INDEX}.txt \
    --REMOVE_DUPLICATES TRUE \
	--ASSUME_SORTED TRUE \
    --TMP_DIR ${BAMDIR}/temp_dir

# index bam alignment
samtools index ${BAMDIR}/DeDup/${INDEX}.dm.bam

if [ ! -f ${REFPREFIX}.fna.fai ]; then
    "Reference sequence needs to be indexed (samtools faidx)"
    samtools faidx ${REFPREFIX} ${REFPREFIX}.fna
fi

# Call variants. -L [only include alleles in given vcf file]
java -jar $GATK HaplotypeCaller  \
   -R $REFSEQ \
   --max-reads-per-alignment-start 100 \
   # -L $POSVCF \
   -I ${BAMDIR}/${INDEX}.dm.bam \
   -O $VCFOUT


'''
for R1PFILE in ${FQDIR}/R1P.trimmed.fastq; do
    export R1FILE=${R1PFILE}
    export R2FILE=${R1PFILE/R1P/R2P}
    # the "$(cut -d'/' -f9 function truncates off the directory portion 
    # of a file path, leaving only the filename
    # -f9 valid for path DATADIR=/oasis/tscc/scratch/mioverto/data/MAseq1 
    # if path is different, the 9 in -f9 needs to be changed to match the # of "/" + 1
    name_only="$(cut -d'/' -f9 <<<${R1FILE})"
    export INDEX=${name_only/_R1P.trimmed.fastq/}
    echo "Submitting ${INDEX}"
    qsub \
        -V \
        -N align_${INDEX} \
        -o ${LOGDIR}/align_${INDEX}_012920.out \
        -e ${LOGDIR}/align_${INDEX}_012920.err \
        /oasis/tscc/scratch/mioverto/code/bowtie_samtools_012920.sbatch
done
'''