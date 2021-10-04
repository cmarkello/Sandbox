#!/bin/bash
cd /data/Udpwork/usr/markellocj/sim_map_eval_150bp
module load python/3.7 singularity
source /data/markellocj/test_toil_vg_run/toil_vg_tools/toilvg_venv/bin/activate
WORK_DIR="/data/Udpwork/usr/markellocj/sim_map_eval_150bp/construct_haplotype_graph_workdir"
rm -fr ${WORK_DIR}/vg-construct-actually_v4.2.1_grch38_HG005-outstore
rm -fr ${WORK_DIR}/tmp-vg-construct-actually_v4.2.1_grch38_HG005
mkdir -p ${WORK_DIR}/vg-construct-actually_v4.2.1_grch38_HG005-outstore
mkdir -p ${WORK_DIR}/tmp-vg-construct-actually_v4.2.1_grch38_HG005
toil clean ${WORK_DIR}/vg-construct-actually_v4.2.1_grch38_HG005-jobstore

toil-vg construct \
${WORK_DIR}/vg-construct-actually_v4.2.1_grch38_HG005-jobstore \
${WORK_DIR}/vg-construct-actually_v4.2.1_grch38_HG005-outstore \
--batchSystem singleMachine \
--container Singularity \
--logInfo \
--workDir ${WORK_DIR}/tmp-vg-construct-actually_v4.2.1_grch38_HG005 \
--cleanWorkDir onSuccess \
--fasta /data/markellocj/fasta_references/grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.gz \
--out_name baseline \
--logFile construct_baseline.actually_v4.2.1_grch38_HG005.log \
--xg_index \
--haplo_sample HG005 \
--fasta_regions \
--vcf /data/markellocj/benchmark_data/HG005_cohort/HG005_GRCh38_1_22_draft_2_v4.2.1_benchmark.no_asterisk.vcf.gz \
--vcf_phasing /data/markellocj/benchmark_data/HG005_cohort/HG005_GRCh38_1_22_draft_2_v4.2.1_benchmark.no_asterisk.vcf.gz \
--statePollingWait 120 \
--rescueJobsFrequency 120 

