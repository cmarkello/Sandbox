#!/usr/bin/env bash
# real_read_deepvariant_call.sh: indel realign input BAMs and run the DeepVariant Genotyper

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

USE_DEFAULT_MODEL=true
WORKDIR=${HOME}/run_deeptrio_genotyping
CHILD_NAME="${1}"
MATERNAL_NAME="${2}"
PATERNAL_NAME="${3}"
BAM_FILE_CHILD="${4}"
BAM_FILE_MATERNAL="${5}"
BAM_FILE_PATERNAL="${6}"
USE_DEFAULT_MODEL=${7}
INPUT_CHILD_BAM="${CHILD_NAME}.bam"
INPUT_MATERNAL_BAM="${MATERNAL_NAME}.bam"
INPUT_PATERNAL_BAM="${PATERNAL_NAME}.bam"
REF_FASTA="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz"
SEQ_DICT="GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict"


# Where should temp files go?
mkdir -p "${WORKDIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
copy ${BAM_FILE_CHILD} "${WORKDIR}/${INPUT_CHILD_BAM}"
copy ${BAM_FILE_MATERNAL} "${WORKDIR}/${INPUT_MATERNAL_BAM}"
copy ${BAM_FILE_PATERNAL} "${WORKDIR}/${INPUT_PATERNAL_BAM}"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz "${WORKDIR}/${REF_FASTA}"
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict -O "${WORKDIR}/${SEQ_DICT}"
if [[ ${CHILD_NAME} == *"HG002"* ]]; then
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
elif [[ ${CHILD_NAME} == *"HG005"* ]]; then
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
    wget_download https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/HG005_GRCh38_1_22_v4.2.1_benchmark.bed "${WORKDIR}/${CHILD_NAME}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
fi
if [[ $USE_DEFAULT_MODEL != true ]]; then 
    wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-child-0806.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0806.tar.gz && tar -xzf ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0806.tar.gz
fi

# Sort and reorder input trio BAMS
for SAMPLE_NAME in "${CHILD_NAME}" "${MATERNAL_NAME}" "${PATERNAL_NAME}"; do
    cd $WORKDIR
    INPUT_BAM=${SAMPLE_NAME}.bam
    docker run \
    -e INPUT_BAM=${INPUT_BAM} \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
    samtools sort -@ 32 ${INPUT_BAM} -O BAM > positionsorted.${INPUT_BAM} && rm ${INPUT_BAM}

    docker run \
    -e INPUT_BAM=${INPUT_BAM} \
    -e SEQ_DICT=${SEQ_DICT} \
    -v ${PWD}:${HOME} -w ${HOME} broadinstitute/picard:2.21.9 \
      java -Xmx20g -XX:ParallelGCThreads=16 -jar /usr/picard/picard.jar \
      ReorderSam \
      VALIDATION_STRINGENCY=SILENT \
      INPUT=positionsorted.${INPUT_BAM} \
      OUTPUT=reordered.positionsorted.${INPUT_BAM} \
      SEQUENCE_DICTIONARY=${SEQ_DICT} && rm positionsorted.${INPUT_BAM}

    # Indel realign input BAM
    cd $WORKDIR
    docker run \
    -e INPUT_BAM=${INPUT_BAM} \
    -e SAMPLE_NAME=${SAMPLE_NAME} \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
    samtools addreplacerg \
    -@ 32 -O BAM -r ID:1 -r LB:lib1 -r SM:${SAMPLE_NAME} -r PL:illumina -r PU:unit1 \
    reordered.positionsorted.${INPUT_BAM} > gatk_ready.reordered.positionsorted.${INPUT_BAM}

    docker run \
    -e INPUT_BAM=${INPUT_BAM} \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/ucsc_cgl/samtools:latest \
    samtools index -@ 32 gatk_ready.reordered.positionsorted.${INPUT_BAM}

    docker run \
    -e SAMPLE_NAME=${SAMPLE_NAME} \
    -e REF_FASTA=${REF_FASTA} \
    -e INPUT_BAM=${INPUT_BAM} \
    -v ${PWD}:${HOME} -w ${HOME} broadinstitute/gatk3:3.8-1 \
      java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
      -R ${REF_FASTA} \
      -I gatk_ready.reordered.positionsorted.${INPUT_BAM} -o ${SAMPLE_NAME}.intervals

    awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' ${SAMPLE_NAME}.intervals > ${SAMPLE_NAME}.intervals.bed
    docker run \
    -e SAMPLE_NAME=${SAMPLE_NAME} \
    -e INPUT_BAM=${INPUT_BAM} \
    -e REF_FASTA=${REF_FASTA} \
    -v ${PWD}:${HOME} -w ${HOME} quay.io/biocontainers/abra2:2.24--h7d875b9_0 \
      --targets ${SAMPLE_NAME}.intervals.bed \
      --in gatk_ready.reordered.positionsorted.${INPUT_BAM} \
      --out indel_realigned.${INPUT_BAM} \
      --ref ${REF_FASTA} \
      --threads 16
done

# Run DeepTrio genotyper on input trio BAMs
cd ${WORKDIR}
mkdir -p ${WORKDIR}/tmp_deeptrio
docker run \
-e REF_FASTA=${REF_FASTA} \
-e CHILD_NAME=${CHILD_NAME} \
-e MATERNAL_NAME=${MATERNAL_NAME} \
-e PATERNAL_NAME=${PATERNAL_NAME} \
-v ${PWD}:${HOME} -w ${HOME} docker://google/deepvariant:deeptrio-1.1.0 \
/bin/bash -c 'seq 0 31 | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/deeptrio/make_examples \
  --mode calling \
  --ref ${REF_FASTA} \
  --reads_parent1 indel_realigned.${PATERNAL_NAME}.bam \
  --reads_parent2 indel_realigned.${MATERNAL_NAME}.bam \
  --reads indel_realigned.${CHILD_NAME}.bam \
  --examples tmp_deeptrio/make_examples.tfrecord@32.gz \
  --sample_name ${CHILD_NAME} \
  --sample_name_parent1 ${PATERNAL_NAME} \
  --sample_name_parent2 ${MATERNAL_NAME} \
  --gvcf gvcf.tfrecord@32.gz \
  --min_mapping_quality 1 \
  --pileup_image_height_child 60 \
  --pileup_image_height_parent 40 \
  --task {}'

# Use the custom deeptrio model if processing VG Pedigre-aligned BAMs
if [[ $USE_DEFAULT_MODEL != true ]]; then
    docker run \
    -v ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0806:${HOME}/dt-giraffe-child-0806 \
    -v ${PWD}:${HOME} -w ${HOME} docker://google/deepvariant:deeptrio-1.1.0 \
    /bin/bash -c '/opt/deepvariant/bin/call_variants \
      --outfile tmp_deeptrio/call_variants_output_child.tfrecord.gz \
      --examples tmp_deeptrio/make_examples_child.tfrecord@32.gz \
      --checkpoint dt-giraffe-child-0806/model.ckpt'
else
    docker run \
    -v ${PWD}:${HOME} -w ${HOME} docker://google/deepvariant:deeptrio-1.1.0 \
    /bin/bash -c '/opt/deepvariant/bin/call_variants \
      --outfile tmp_deeptrio/call_variants_output_child.tfrecord.gz \
      --examples tmp_deeptrio/make_examples_child.tfrecord@32.gz'
fi

docker run \
-e REF_FASTA=${REF_FASTA} \
-e CHILD_NAME=${CHILD_NAME} \
-v ${PWD}:${HOME} -w ${HOME} docker://google/deepvariant:deeptrio-1.1.0 \
/bin/bash -c '/opt/deepvariant/bin/postprocess_variants \
  --ref ${REF_FASTA} \ 
  --infile tmp_deeptrio/call_variants_output_child.tfrecord.gz \
  --outfile ${CHILD_NAME}_DEEPTRIO.abra_gatk_targets.vcf.gz \
  --nonvariant_site_tfrecord_path tmp_deeptrio/gvcf_child.tfrecord@32.gz \
  --gvcf_outfile ${CHILD_NAME}_DEEPTRIO.abra_gatk_targets.g.vcf.gz'


