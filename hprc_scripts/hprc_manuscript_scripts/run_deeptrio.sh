#!/bin/bash
# run_deeptrio.sh HG002 HG002 HG003 HG004 /data/Udpbinfo/usr/markellocj/hprc_feb3_giraffe_runs indel_realigned.bwamem_ GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.hprc_july7.fna
module load bwa samtools/1.13 picard bcftools singularity
COHORT_NAME=$1
SAMPLE_NAME=$2
PATERNAL_NAME=$3
MATERNAL_NAME=$4
WORKDIR=$5
BAM_BASENAME=$6
FASTA_FILE=$7
cd ${WORKDIR}

mkdir ${SAMPLE_NAME}_bam_by_chr
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    samtools view -@ 32 -b ${BAM_BASENAME}${SAMPLE_NAME}.bam GRCh38.chr${CHR} > ${SAMPLE_NAME}_bam_by_chr/raw.GRCh38.chr${CHR}.bam
    samtools index -@ 32 ${SAMPLE_NAME}_bam_by_chr/raw.GRCh38.chr${CHR}.bam
done

mkdir ${PATERNAL_NAME}_bam_by_chr
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    samtools view -@ 32 -b ${BAM_BASENAME}${PATERNAL_NAME}.bam GRCh38.chr${CHR} > ${PATERNAL_NAME}_bam_by_chr/raw.GRCh38.chr${CHR}.bam
    samtools index -@ 32 ${PATERNAL_NAME}_bam_by_chr/raw.GRCh38.chr${CHR}.bam
done

mkdir ${MATERNAL_NAME}_bam_by_chr
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    samtools view -@ 32 -b ${BAM_BASENAME}${MATERNAL_NAME}.bam GRCh38.chr${CHR} > ${MATERNAL_NAME}_bam_by_chr/raw.GRCh38.chr${CHR}.bam
    samtools index -@ 32 ${MATERNAL_NAME}_bam_by_chr/raw.GRCh38.chr${CHR}.bam
done


rm deeptrio_calling.${SAMPLE_NAME}.default.minmapq1.make_examples.swarm
rm deeptrio_calling.${SAMPLE_NAME}.default.minmapq1.call_variants.child.swarm
rm deeptrio_calling.${PATERNAL_NAME}.default.minmapq1.call_variants.parent1.swarm
rm deeptrio_calling.${MATERNAL_NAME}.default.minmapq1.call_variants.parent2.swarm
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    cd ${WORKDIR}
    mkdir tmp_deeptrio_${SAMPLE_NAME}_${CHR}
    echo -e "cd ${WORKDIR}; singularity run -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/deepvar_models:\${HOME}/deepvar_models -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference /data/markellocj/singularity_cache/cache/oci-tmp/eb983d45c44c75d8e558c5c9746128b70b565a4aa6821cb40ab0daf242c578a8 /bin/bash -c 'seq 0 31 | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/deeptrio/make_examples --mode calling --regions GRCh38.chr${CHR} --ref grch38_reference/${FASTA_FILE} --reads_parent1 ${PATERNAL_NAME}_bam_by_chr/raw.GRCh38.chr${CHR}.bam --reads_parent2 ${MATERNAL_NAME}_bam_by_chr/raw.GRCh38.chr${CHR}.bam --reads ${SAMPLE_NAME}_bam_by_chr/raw.GRCh38.chr${CHR}.bam --examples tmp_deeptrio_${SAMPLE_NAME}_${CHR}/make_examples.tfrecord@32.gz --sample_name ${SAMPLE_NAME} --sample_name_parent1 ${PATERNAL_NAME} --sample_name_parent2 ${MATERNAL_NAME} --gvcf tmp_deeptrio_${SAMPLE_NAME}_${CHR}/gvcf.tfrecord@32.gz --min_mapping_quality 1 --pileup_image_height_child 60 --pileup_image_height_parent 40 --keep_legacy_allele_counter_behavior=true --normalize_reads=true --task {}'" >> deeptrio_calling.${SAMPLE_NAME}.default.minmapq1.make_examples.swarm
    echo -e "cd ${WORKDIR}; singularity run --nv -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/deepvar_models:\${HOME}/deepvar_models -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference /data/markellocj/singularity_cache/cache/oci-tmp/ba2a9b6dd3cfd7c01d73a2d0aba90cd4c3f23c12fdb6344b82d468244de6f1bd /bin/bash -c '/opt/deepvariant/bin/call_variants --checkpoint /opt/models/deeptrio/wgs/child/model.ckpt --outfile ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/call_variants_output_child.tfrecord.gz --examples ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/make_examples_child.tfrecord@32.gz'; singularity run --nv -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/deepvar_models:\${HOME}/deepvar_models -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference /data/markellocj/singularity_cache/cache/oci-tmp/ba2a9b6dd3cfd7c01d73a2d0aba90cd4c3f23c12fdb6344b82d468244de6f1bd /bin/bash -c '/opt/deepvariant/bin/postprocess_variants --ref grch38_reference/${FASTA_FILE} --infile ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/call_variants_output_child.tfrecord.gz --outfile ${SAMPLE_NAME}.trianed.GRCh38.chr${CHR}.minmapq1.vcf.gz --nonvariant_site_tfrecord_path ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/gvcf_child.tfrecord@32.gz --gvcf_outfile ${SAMPLE_NAME}.default.GRCh38.chr${CHR}.minmapq1.g.vcf.gz'" >> deeptrio_calling.${SAMPLE_NAME}.default.minmapq1.call_variants.child.swarm
    echo -e "cd ${WORKDIR}; singularity run --nv -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/deepvar_models:\${HOME}/deepvar_models -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference /data/markellocj/singularity_cache/cache/oci-tmp/ba2a9b6dd3cfd7c01d73a2d0aba90cd4c3f23c12fdb6344b82d468244de6f1bd /bin/bash -c '/opt/deepvariant/bin/call_variants --checkpoint /opt/models/deeptrio/wgs/parent/model.ckpt --outfile ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/call_variants_output_parent1.tfrecord.gz --examples ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/make_examples_parent1.tfrecord@32.gz';  singularity run --nv -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/deepvar_models:\${HOME}/deepvar_models -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference /data/markellocj/singularity_cache/cache/oci-tmp/ba2a9b6dd3cfd7c01d73a2d0aba90cd4c3f23c12fdb6344b82d468244de6f1bd /bin/bash -c '/opt/deepvariant/bin/postprocess_variants --ref grch38_reference/${FASTA_FILE} --infile ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/call_variants_output_parent1.tfrecord.gz --outfile ${PATERNAL_NAME}.trianed.GRCh38.chr${CHR}.minmapq1.vcf.gz --nonvariant_site_tfrecord_path ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/gvcf_parent1.tfrecord@32.gz --gvcf_outfile ${PATERNAL_NAME}.default.GRCh38.chr${CHR}.minmapq1.g.vcf.gz'" >> deeptrio_calling.${PATERNAL_NAME}.default.minmapq1.call_variants.parent1.swarm
    echo -e "cd ${WORKDIR}; singularity run --nv -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/deepvar_models:\${HOME}/deepvar_models -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference /data/markellocj/singularity_cache/cache/oci-tmp/ba2a9b6dd3cfd7c01d73a2d0aba90cd4c3f23c12fdb6344b82d468244de6f1bd /bin/bash -c '/opt/deepvariant/bin/call_variants --checkpoint /opt/models/deeptrio/wgs/parent/model.ckpt --outfile ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/call_variants_output_parent2.tfrecord.gz --examples ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/make_examples_parent2.tfrecord@32.gz';  singularity run --nv -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/deepvar_models:\${HOME}/deepvar_models -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference /data/markellocj/singularity_cache/cache/oci-tmp/ba2a9b6dd3cfd7c01d73a2d0aba90cd4c3f23c12fdb6344b82d468244de6f1bd /bin/bash -c '/opt/deepvariant/bin/postprocess_variants --ref grch38_reference/${FASTA_FILE} --infile ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/call_variants_output_parent2.tfrecord.gz --outfile ${MATERNAL_NAME}.trianed.GRCh38.chr${CHR}.minmapq1.vcf.gz --nonvariant_site_tfrecord_path ./tmp_deeptrio_${SAMPLE_NAME}_${CHR}/gvcf_parent2.tfrecord@32.gz --gvcf_outfile ${MATERNAL_NAME}.default.GRCh38.chr${CHR}.minmapq1.g.vcf.gz'" >> deeptrio_calling.${MATERNAL_NAME}.default.minmapq1.call_variants.parent2.swarm
done

swarm_make_examples_jid=$(swarm -f deeptrio_calling.${SAMPLE_NAME}.default.minmapq1.make_examples.swarm -t 32 -g 50 --time=6:00:00 --module singularity)
echo ${swarm_make_examples_jid}
swarm_call_child_jid=$(swarm -f deeptrio_calling.${SAMPLE_NAME}.default.minmapq1.call_variants.child.swarm -t 16 -g 20 --partition=gpu --gres=gpu:k80:1 --time=6:00:00 --dependency=afterany:$swarm_make_examples_jid --module singularity)
echo ${swarm_call_child_jid}
swarm_call_parent1_jid=$(swarm -f deeptrio_calling.${PATERNAL_NAME}.default.minmapq1.call_variants.parent1.swarm -t 16 -g 20 --partition=gpu --gres=gpu:k80:1 --time=6:00:00 --dependency=afterany:$swarm_make_examples_jid --module singularity)
echo ${swarm_call_parent1_jid}
swarm_call_parent2_jid=$(swarm -f deeptrio_calling.${MATERNAL_NAME}.default.minmapq1.call_variants.parent2.swarm -t 16 -g 20 --partition=gpu --gres=gpu:k80:1 --time=6:00:00 --dependency=afterany:$swarm_make_examples_jid --module singularity)
echo ${swarm_call_parent2_jid}

cd ${WORKDIR}
rm merge_${SAMPLE_NAME}_deeptrio_vcfs.sh
echo -e "#!/bin/bash" >> merge_${SAMPLE_NAME}_deeptrio_vcfs.sh
echo -e "module load bcftools; cd ${WORKDIR};
bcftools concat -O v ${SAMPLE_NAME}.trianed.GRCh38.chr{?,??}.minmapq1.vcf.gz > ${SAMPLE_NAME}.merged.deeptrio.indel_realigned.vcf \
&& bcftools sort ${SAMPLE_NAME}.merged.deeptrio.indel_realigned.vcf -O v > ${SAMPLE_NAME}.merged.deeptrio.indel_realigned.sorted.vcf && rm ${SAMPLE_NAME}.merged.deeptrio.indel_realigned.vcf \
&& bgzip -f ${SAMPLE_NAME}.merged.deeptrio.indel_realigned.sorted.vcf \
&& tabix -p vcf ${SAMPLE_NAME}.merged.deeptrio.indel_realigned.sorted.vcf.gz" >> merge_${SAMPLE_NAME}_deeptrio_vcfs.sh

sbatch -c 8 --mem=50GB --time=4:00:00 --dependency=afterany:$swarm_call_child_jid merge_${SAMPLE_NAME}_deeptrio_vcfs.sh

cd ${WORKDIR}
rm merge_${PATERNAL_NAME}_deeptrio_vcfs.sh
echo -e "#!/bin/bash" >> merge_${PATERNAL_NAME}_deeptrio_vcfs.sh
echo -e "module load bcftools; cd ${WORKDIR};
bcftools concat -O v ${PATERNAL_NAME}.trianed.GRCh38.chr{?,??}.minmapq1.vcf.gz > ${PATERNAL_NAME}.merged.deeptrio.indel_realigned.vcf \
&& bcftools sort ${PATERNAL_NAME}.merged.deeptrio.indel_realigned.vcf -O v > ${PATERNAL_NAME}.merged.deeptrio.indel_realigned.sorted.vcf && rm ${PATERNAL_NAME}.merged.deeptrio.indel_realigned.vcf \
&& bgzip -f ${PATERNAL_NAME}.merged.deeptrio.indel_realigned.sorted.vcf \
&& tabix -p vcf ${PATERNAL_NAME}.merged.deeptrio.indel_realigned.sorted.vcf.gz" >> merge_${PATERNAL_NAME}_deeptrio_vcfs.sh

sbatch -c 8 --mem=50GB --time=4:00:00 --dependency=afterany:$swarm_call_parent1_jid merge_${PATERNAL_NAME}_deeptrio_vcfs.sh

cd ${WORKDIR}
rm merge_${MATERNAL_NAME}_deeptrio_vcfs.sh
echo -e "#!/bin/bash" >> merge_${MATERNAL_NAME}_deeptrio_vcfs.sh
echo -e "module load bcftools; cd ${WORKDIR};
bcftools concat -O v ${MATERNAL_NAME}.trianed.GRCh38.chr{?,??}.minmapq1.vcf.gz > ${MATERNAL_NAME}.merged.deeptrio.indel_realigned.vcf \
&& bcftools sort ${MATERNAL_NAME}.merged.deeptrio.indel_realigned.vcf -O v > ${MATERNAL_NAME}.merged.deeptrio.indel_realigned.sorted.vcf && rm ${MATERNAL_NAME}.merged.deeptrio.indel_realigned.vcf \
&& bgzip -f ${MATERNAL_NAME}.merged.deeptrio.indel_realigned.sorted.vcf \
&& tabix -p vcf ${MATERNAL_NAME}.merged.deeptrio.indel_realigned.sorted.vcf.gz" >> merge_${MATERNAL_NAME}_deeptrio_vcfs.sh

sbatch -c 8 --mem=50GB --time=4:00:00 --dependency=afterany:$swarm_call_parent2_jid merge_${MATERNAL_NAME}_deeptrio_vcfs.sh

