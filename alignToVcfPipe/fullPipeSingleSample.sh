# Alignment

module load bowtie2
module load samtools
module load bcftools

export REFPREFIX=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_rnm_UCSD_2020
export REFSEQ=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/RM_refseq_rnm_UCSD_2020.fna
export DATADIR=/oasis/tscc/scratch/mioverto/data/MAseq1

export FQDIR=${DATADIR}/reads/trim
export BAMDIR=${DATADIR}/bam

export R1FILE=${DATADIR}/reads/trim/N_A05_1_R1P.trimmed.fastq
export R2FILE=${DATADIR}/reads/trim/N_A05_1_R2P.trimmed.fastq
export R1UFILE=${DATADIR}/reads/trim/N_A05_1_R1U.trimmed.fastq
export R2UFILE=${DATADIR}/reads/trim/N_A05_1_R2U.trimmed.fastq

export INDEX=N_A05_RM

export VCFOUT=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/N_A05_RM.vcf.gz
export VCFSLIM=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/N_A05_RM.vcf.gz

export METRICS=/oasis/tscc/scratch/mioverto/data/MAseq1/bam/picard_metrics

export POSFILE=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/POS_files/Bloom_pos_crctd.txt

/oasis/tscc/scratch/mioverto/code/bowtie_samtools_012920.sbatch

# Start with alignment. Assumes renamed, trimmed fastq files and un-indexed reference fasta file

if [ ! -f ${REFPREFIX}.1.bt2 ]; then
    bowtie2-build ${REFSEQ} ${REFPREFIX}
fi

bowtie2 --rg-id ${INDEX} --rg SM:${INDEX} -X 1000 -x ${REFPREFIX} -1 ${R1FILE} -2 ${R2FILE} -U "${R1UFILE},${R2UFILE}" \
| samtools view -h -b | samtools sort -m 10000000 -o ${BAMDIR}/${INDEX}.bam -O bam -T ${BAMDIR}/${INDEX}

GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar
# PATH=/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.252.b09-2.el7_8.x86_64/bin:$PATH

java -jar $GATK MarkDuplicates \
	--INPUT ${BAMDIR}/${INDEX}.bam \
	--OUTPUT ${BAMDIR}/${INDEX}.dm.bam \
	--METRICS_FILE ${METRICS}/picard-${INDEX}.txt \
    --REMOVE_DUPLICATES TRUE \
	--ASSUME_SORTED TRUE \
    --TMP_DIR ${BAMDIR}/temp_dir

samtools index ${BAMDIR}/${INDEX}.dm.bam

# export REFSEQ=/oasis/tscc/scratch/mioverto/data/S288C_R63_refseq.fna
# export BAMDIR=/oasis/tscc/scratch/mioverto/data/MAseq1/bam/DeDup
# export INDEX=full-GM50_1
# export POSFILE=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/POS_files/Bloom_pos_crctd.txt
# export VCFOUT=/oasis/tscc/scratch/mioverto/data/MAseq1/variant_calls/full-GM50_14shft.vcf.gz

# bcftools mpileup -f $REFSEQ ${BAMDIR}/${INDEX}.dm.bam | bcftools call -m -Oz -T $POSFILE -o $VCFOUT

POSFILE=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/bam/vcf/RM_ref.vcf.gz
POSVCF=/oasis/tscc/scratch/mioverto/data/refseq/RM_ref/vcf/RMxBY_variants.vcf

bcftools mpileup -f $REFSEQ ${BAMDIR}/${INDEX}.dm.bam | bcftools call -m -Oz -T $POSVCF -o $VCFOUT



GATK=/home/mioverto/bin/gatk_4/gatk-package-4.1.8.0-local.jar

java -jar $GATK HaplotypeCaller  \
   -R $REFSEQ \
   --max-reads-per-alignment-start 100 \
   -L $POSVCF \
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