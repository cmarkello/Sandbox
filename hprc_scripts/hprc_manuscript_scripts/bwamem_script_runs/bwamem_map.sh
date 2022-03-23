#!/usr/bin/env bash
# bwamem_map.sh: map HG001, NA12891, NA12892, HG002, HG003, HG004, HG005, HG006, HG007 150bp paired reads with BWAMEM to the HS38d1-based graph references
# bwamem_map.sh /data/Udpbinfo/usr/markellocj/hprc_feb3_bwamem_runs HG002 GRCh38-f1g-90-mc-aug11-clip.d9.m1000.D10M.m1000 GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7 .novaseq.pcr-free.30x
module load singularity

set -ex
set -o pipefail

function download() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}"
    fi
}

function wget_download() {
    if [ ! -e "${2}" ] ; then
        wget "${1}" -O "${2}"
    fi
}

function copy() {
    if [ ! -e "${2}" ] ; then
        cp "${1}" "${2}"
    fi
}

WORKDIR=$1
SAMPLE_NAME=$2
INDEX_BASENAME=$3
REF_FASTA_BASENAME=$4
READ_BASENAME=$5
REF_FASTA="${REF_FASTA_BASENAME}.fna"
INPUT_SAM=bwamem_${SAMPLE_NAME}.sam
INPUT_BAM=bwamem_${SAMPLE_NAME}.bam

# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp_${SAMPLE_NAME}"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR

SEQ_DICT="${REF_FASTA_BASENAME}.dict"

# run bwamem mapper
#cd $WORKDIR
#echo ${SAMPLE_NAME}${READ_BASENAME}
#singularity exec \
#--env REF_FASTA=${REF_FASTA},SAMPLE_NAME="${SAMPLE_NAME}",READ_BASENAME="${READ_BASENAME}" \
#-H ${PWD}:${HOME} --pwd ${HOME} docker://biocontainers/bwa:v0.7.17_cv1 \
#bwa mem \
#-t 32 \
#${REF_FASTA} \
#"${SAMPLE_NAME}${READ_BASENAME}.R1.fastq.gz" \
#"${SAMPLE_NAME}${READ_BASENAME}.R2.fastq.gz" \
#> bwamem_${SAMPLE_NAME}.sam


#cd $WORKDIR
#singularity exec \
#--env SAMPLE_NAME=${SAMPLE_NAME} \
#--env INPUT_SAM=${INPUT_SAM} \
#-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/ucsc_cgl/samtools:latest \
#    samtools sort -@ 32 ${INPUT_SAM} -O BAM > positionsorted.bwamem_${SAMPLE_NAME}.bam && rm ${INPUT_SAM}

#singularity exec \
#--env INPUT_BAM=${INPUT_BAM},SEQ_DICT=${SEQ_DICT} \
#-H ${PWD}:${HOME} --pwd ${HOME} docker://broadinstitute/picard:2.21.9 \
#  java -Xmx20g -XX:ParallelGCThreads=16 -jar /usr/picard/picard.jar \
#  ReorderSam \
#  VALIDATION_STRINGENCY=SILENT \
#  INPUT=positionsorted.${INPUT_BAM} \
#  OUTPUT=reordered.positionsorted.${INPUT_BAM} \
#  SEQUENCE_DICTIONARY=${SEQ_DICT} && rm positionsorted.${INPUT_BAM}

# Indel realign input BAM
#cd $WORKDIR
#singularity exec \
#--env INPUT_BAM=${INPUT_BAM},SAMPLE_NAME=${SAMPLE_NAME} \
#-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/ucsc_cgl/samtools:latest \
#samtools addreplacerg \
#-@ 32 -O BAM -r ID:1 -r LB:lib1 -r SM:${SAMPLE_NAME} -r PL:illumina -r PU:unit1 \
#reordered.positionsorted.${INPUT_BAM} > gatk_ready.reordered.positionsorted.${INPUT_BAM}

#singularity exec \
#--env INPUT_BAM=${INPUT_BAM} \
#-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/ucsc_cgl/samtools:latest \
#samtools index -@ 32 gatk_ready.reordered.positionsorted.${INPUT_BAM}

#singularity exec \
#--env SAMPLE_NAME=${SAMPLE_NAME},REF_FASTA=${REF_FASTA},INPUT_BAM=${INPUT_BAM} \
#-H ${PWD}:${HOME} --pwd ${HOME} docker://broadinstitute/gatk3:3.8-1 \
#  java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
#  -R ${REF_FASTA} \
#  -I gatk_ready.reordered.positionsorted.${INPUT_BAM} -o ${SAMPLE_NAME}.intervals

awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' ${SAMPLE_NAME}.intervals > ${SAMPLE_NAME}.intervals.bed
singularity exec \
--env SAMPLE_NAME=${SAMPLE_NAME},INPUT_BAM=${INPUT_BAM},REF_FASTA=${REF_FASTA} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/biocontainers/abra2:2.24--h7d875b9_0 \
  abra2 \
  --targets ${SAMPLE_NAME}.intervals.bed \
  --in gatk_ready.reordered.positionsorted.${INPUT_BAM} \
  --out indel_realigned.${INPUT_BAM} \
  --ref ${REF_FASTA} \
  --threads 32

singularity exec \
--env INPUT_BAM=${INPUT_BAM} \
-H ${PWD}:${HOME} --pwd ${HOME} docker://quay.io/ucsc_cgl/samtools:latest \
    samtools index -@ 32 indel_realigned.${INPUT_BAM}

