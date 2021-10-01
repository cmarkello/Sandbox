#!/usr/bin/env bash
# giraffe_map.sh: map HG002, HG003, 150bp and 250bp paired reads with VG GIRAFFE, VG GIRAFFE FAST, and VG GIRAFFE PRIMARY to the HS38d1-based graph references

set -ex
set -o pipefail

function download() {
    if [ ! -e "${2}" ] ; then
        aws s3 cp --no-progress "${1}" "${2}"
    fi
}

function wget_download() {
    if [ ! -e "${2}" ] ; then
        wget_download "${1}" -O "${2}"
    fi
}

function copy() {
    if [ ! -e "${2}" ] ; then
        cp "${1}" "${2}"
    fi
}

WORKDIR=${HOME}/run_giraffe_pedigree_mapping
WORKFLOW_INPUT_DIR=${WORKDIR}/inputs

# Where should temp files go?
mkdir -p "${WORKDIR}" "${WORKFLOW_INPUT_DIR}"
export TMPDIR="${WORKDIR}/tmp"
mkdir -p "${TMPDIR}"

# Download data input data
cd $WORKDIR
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/path_list_whole_genome.txt -O ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.xg -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.xg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.gbwt -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.gbwt
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.gg -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.gg
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.min -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.min
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/liftover_snp1kg_grch38_nosegdup.dist -O ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.dist
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict -O ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/snpEff_v5_0_GRCh38.99.zip -O ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh38.99.zip
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/eagle_data_grch38.tar.gz -O ${WORKFLOW_INPUT_DIR}/eagle_data_grch38.tar.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-child-0806.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-child-0806.tar.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dt-giraffe-parent-0806.tar.gz -O ${WORKFLOW_INPUT_DIR}/dt-giraffe-parent-0806.tar.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/grch38_inputs/dv-giraffe-0507.tar.gz -O ${WORKFLOW_INPUT_DIR}/dv-giraffe-0507.tar.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG002.ped -O ${WORKFLOW_INPUT_DIR}/HG002.ped
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG005.ped -O ${WORKFLOW_INPUT_DIR}/HG005.ped

# Download input reads
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG002.novaseq.pcr-free.35x.R1.fastq.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG002.novaseq.pcr-free.35x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/HG002.novaseq.pcr-free.35x.R2.fastq.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG003.novaseq.pcr-free.35x.R1.fastq.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG003.novaseq.pcr-free.35x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/HG003.novaseq.pcr-free.35x.R2.fastq.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG004.novaseq.pcr-free.35x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG004.novaseq.pcr-free.35x.R1.fastq.gz
wget_download https://storage.googleapis.com/cmarkell-vg-wdl-dev/test_input_reads/HG004.novaseq.pcr-free.35x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/HG004.novaseq.pcr-free.35x.R2.fastq.gz

wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG005.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG005.novaseq.pcr-free.30x.R1.fastq.gz
wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG005.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/HG005.novaseq.pcr-free.30x.R2.fastq.gz
wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG006.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG006.novaseq.pcr-free.30x.R1.fastq.gz
wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG006.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/HG006.novaseq.pcr-free.30x.R2.fastq.gz
wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG007.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG007.novaseq.pcr-free.30x.R1.fastq.gz
wget_download https://storage.googleapis.com/deepvariant/benchmarking/fastq/wgs_pcr_free/30x/HG007.novaseq.pcr-free.30x.R2.fastq.gz ${WORKFLOW_INPUT_DIR}/HG007.novaseq.pcr-free.30x.R2.fastq.gz

git clone --single-branch --branch vg_pedigree_workflow_deepvariant_grch38 https://github.com/vgteam/toil-vg.git 
git clone https://github.com/cmarkello/toil.git
python3 -m venv toilvg_venv
source toilvg_venv/bin/activate
pip install ./toil
pip install ./toil-vg
deactivate


# Run the VG Pedigree workflow for the HG002 cohort
source $WORKDIR/toilvg_venv/bin/activate
OUTSTORE="${WORKDIR}/HG002_vg_pedigree_outstore"
JOBSTORE="${WORKDIR}/HG002_vg_pedigree_jobstore"
LOGFILE="${WORKDIR}/HG002_vg_pedigree.log"
TMPDIR="${WORKDIR}/tmp_HG002_vg_pedigree"
export TOIL_SLURM_ARGS='-t 20:00:00'
export SINGULARITY_CACHEDIR=${WORKDIR}/singularity_cache
rm -fr ${LOGFILE} ${OUTSTORE} ${TMPDIR}
mkdir -p ${OUTSTORE} ${TMPDIR} $SINGULARITY_CACHEDIR
cd ${WORKDIR}

toil clean ${JOBSTORE}

time toil-vg pedigree \
--genome_build "GRCh38" \
--retryCount 0 \
--rotatingLogging \
--setEnv PATH=$PATH \
--disableProgress \
--realTimeLogging \
--batchSystem slurm \
--statePollingWait 120 \
--rescueJobsFrequency 120 \
--container Singularity \
--logInfo \
--logFile ${LOGFILE} \
--workDir ${TMPDIR} \
--cleanWorkDir onSuccess \
--whole_genome_config \
--vg_docker 'quay.io/vgteam/vg:v1.31.0' \
${JOBSTORE} \
${OUTSTORE} \
HG002 \
HG004 \
HG003 \
--mapper giraffe \
--caller deepvariant \
--sibling_genders 0 \
--sibling_affected 0 \
--fastq_proband ${WORKFLOW_INPUT_DIR}/HG002.novaseq.pcr-free.35x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG002.novaseq.pcr-free.35x.R2.fastq.gz \
--fastq_maternal ${WORKFLOW_INPUT_DIR}/HG004.novaseq.pcr-free.35x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG004.novaseq.pcr-free.35x.R2.fastq.gz \
--fastq_paternal ${WORKFLOW_INPUT_DIR}/HG003.novaseq.pcr-free.35x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG003.novaseq.pcr-free.35x.R2.fastq.gz \
--ref_fasta ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz \
--ref_fasta_index ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai \
--ref_fasta_dict ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict \
--deeptrio_child_model /data/Udpbinfo/Scratch/markellocj/toil_vg_workflow_inputs/toil_vg_inputs/grch38_inputs/dt-giraffe-child-0806.tar.gz \
--deeptrio_parent_model /data/Udpbinfo/Scratch/markellocj/toil_vg_workflow_inputs/toil_vg_inputs/grch38_inputs/dt-giraffe-parent-0806.tar.gz \
--deepvariant_model ${WORKFLOW_INPUT_DIR}/dv-giraffe-0507.tar.gz \
--use_haplotypes \
--xg_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.xg \
--gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gbwt \
--graph_gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gg \
--minimizer_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.min \
--distance_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.dist \
--id_ranges ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \
--path_list ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \
--ped_file ${WORKFLOW_INPUT_DIR}/HG002.ped \
--eagle_data ${WORKFLOW_INPUT_DIR}/eagle_data_grch38.tar.gz \
--snpeff_database ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh38.99.zip \
--bam_output \
--use_decoys \
--indel_realign_bams \
--snpeff_annotation

# Run the VG Pedigree workflow for the HG005 cohort
source $WORKDIR/toilvg_venv/bin/activate
OUTSTORE="${WORKDIR}/HG005_vg_pedigree_outstore"
JOBSTORE="${WORKDIR}/HG005_vg_pedigree_jobstore"
LOGFILE="${WORKDIR}/HG005_vg_pedigree.log"
TMPDIR="${WORKDIR}/tmp_HG005_vg_pedigree"
export TOIL_SLURM_ARGS='-t 20:00:00'
export SINGULARITY_CACHEDIR=${WORKDIR}/singularity_cache
rm -fr ${LOGFILE} ${OUTSTORE} ${TMPDIR}
mkdir -p ${OUTSTORE} ${TMPDIR} $SINGULARITY_CACHEDIR
cd ${WORKDIR}

toil clean ${JOBSTORE}

time toil-vg pedigree \
--genome_build "GRCh38" \
--retryCount 0 \
--rotatingLogging \
--setEnv PATH=$PATH \
--disableProgress \
--realTimeLogging \
--batchSystem slurm \
--statePollingWait 120 \
--rescueJobsFrequency 120 \
--container Singularity \
--logInfo \
--logFile ${LOGFILE} \
--workDir ${TMPDIR} \
--cleanWorkDir onSuccess \
--whole_genome_config \
--vg_docker 'quay.io/vgteam/vg:v1.31.0' \
${JOBSTORE} \
${OUTSTORE} \
HG005 \
HG007 \
HG006 \
--mapper giraffe \
--caller deepvariant \
--sibling_genders 0 \
--sibling_affected 0 \
--fastq_proband ${WORKFLOW_INPUT_DIR}/HG005.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG005.novaseq.pcr-free.30x.R2.fastq.gz \
--fastq_maternal ${WORKFLOW_INPUT_DIR}/HG007.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG007.novaseq.pcr-free.30x.R2.fastq.gz \
--fastq_paternal ${WORKFLOW_INPUT_DIR}/HG006.novaseq.pcr-free.30x.R1.fastq.gz ${WORKFLOW_INPUT_DIR}/HG006.novaseq.pcr-free.30x.R2.fastq.gz \
--ref_fasta ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz \
--ref_fasta_index ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai \
--ref_fasta_dict ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict \
--deeptrio_child_model /data/Udpbinfo/Scratch/markellocj/toil_vg_workflow_inputs/toil_vg_inputs/grch38_inputs/dt-giraffe-child-0806.tar.gz \
--deeptrio_parent_model /data/Udpbinfo/Scratch/markellocj/toil_vg_workflow_inputs/toil_vg_inputs/grch38_inputs/dt-giraffe-parent-0806.tar.gz \
--deepvariant_model ${WORKFLOW_INPUT_DIR}/dv-giraffe-0507.tar.gz \
--use_haplotypes \
--xg_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.xg \
--gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gbwt \
--graph_gbwt_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.gg \
--minimizer_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.min \
--distance_index ${WORKFLOW_INPUT_DIR}/liftover_snp1kg_grch38_nosegdup.dist \
--id_ranges ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \
--path_list ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \
--ped_file ${WORKFLOW_INPUT_DIR}/HG005.ped \
--eagle_data ${WORKFLOW_INPUT_DIR}/eagle_data_grch38.tar.gz \
--snpeff_database ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh38.99.zip \
--bam_output \
--use_decoys \
--indel_realign_bams \
--snpeff_annotation
