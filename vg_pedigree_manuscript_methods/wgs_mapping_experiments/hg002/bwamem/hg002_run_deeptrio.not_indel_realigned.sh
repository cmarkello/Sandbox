#!/bin/bash
module load bwa samtools/1.11 picard bcftools singularity
COHORT_NAME="HG002"
SAMPLE_NAME="HG002"
PATERNAL_NAME="HG003"
MATERNAL_NAME="HG004"
LOW_COHORT_NAME=$(echo $COHORT_NAME | tr '[:upper:]' '[:lower:]')
LOW_SAMPLE_NAME=$(echo $SAMPLE_NAME | tr '[:upper:]' '[:lower:]')
LOW_PATERNAL_NAME=$(echo $PATERNAL_NAME | tr '[:upper:]' '[:lower:]')
LOW_MATERNAL_NAME=$(echo $MATERNAL_NAME | tr '[:upper:]' '[:lower:]')
cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem/
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    samtools index -@ 32 bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr/raw.chr${CHR}.bam
done

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    samtools index -@ 32 bwamem_grch38_${LOW_PATERNAL_NAME}_bam_by_chr/raw.chr${CHR}.bam
done

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    samtools index -@ 32 bwamem_grch38_${LOW_MATERNAL_NAME}_bam_by_chr/raw.chr${CHR}.bam
done


rm deepvtrio_calling.${LOW_SAMPLE_NAME}.default.minmapq1.swarm
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y ; do
    echo -e "cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem/; singularity run --nv -H \${PWD}:\${HOME} --pwd \${HOME} -B /data/markellocj/fasta_references/grch38_reference:\${HOME}/grch38_reference docker://google/deepvariant:deeptrio-1.1.0-gpu /opt/deepvariant/bin/deeptrio/run_deeptrio --make_examples_extra_args 'min_mapping_quality=1' --model_type=WGS --regions chr${CHR} --ref=grch38_reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna --reads_child=bwamem_grch38_${LOW_SAMPLE_NAME}_bam_by_chr/raw.chr${CHR}.bam --reads_parent1=bwamem_grch38_${LOW_PATERNAL_NAME}_bam_by_chr/raw.chr${CHR}.bam --reads_parent2=bwamem_grch38_${LOW_MATERNAL_NAME}_bam_by_chr/raw.chr${CHR}.bam --sample_name_child='${SAMPLE_NAME}' --sample_name_parent1='${PATERNAL_NAME}' --sample_name_parent2='${MATERNAL_NAME}' --output_vcf_child=${SAMPLE_NAME}.bwamem.deeptrio.chr${CHR}.vcf.gz --output_vcf_parent1=${PATERNAL_NAME}.bwamem.deeptrio.chr${CHR}.vcf.gz --output_vcf_parent2=${MATERNAL_NAME}.bwamem.deeptrio.chr${CHR}.vcf.gz --output_gvcf_child=${SAMPLE_NAME}.bwamem.deeptrio.chr${CHR}.g.vcf.gz --output_gvcf_parent1=${PATERNAL_NAME}.bwamem.deeptrio.chr${CHR}.g.vcf.gz --output_gvcf_parent2=${MATERNAL_NAME}.bwamem.deeptrio.chr${CHR}.g.vcf.gz --num_shards=16" >> deepvtrio_calling.${LOW_SAMPLE_NAME}.default.minmapq1.swarm
done

swarm_jid=$(swarm -f deepvtrio_calling.${LOW_SAMPLE_NAME}.default.minmapq1.swarm -t 16 -g 20 --partition=gpu --gres=gpu:k80:1 --time=6:00:00 --module singularity)

cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem/
rm merge_${LOW_SAMPLE_NAME}_deeptrio_vcfs.not_indel_realigned.sh
echo -e "#!/bin/bash" >> merge_${LOW_SAMPLE_NAME}_deeptrio_vcfs.not_indel_realigned.sh
echo -e "cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/${LOW_COHORT_NAME}_cohort/bwamem/;
bcftools concat -O v ${SAMPLE_NAME}.bwamem.deeptrio.chr{?,??}.vcf.gz > merged.deeptrio.vcf \
&& bcftools sort merged.deeptrio.vcf -O v > merged.deeptrio.sorted.vcf && rm merged.deeptrio.vcf \
&& bgzip merged.deeptrio.sorted.vcf \
&& tabix -p vcf merged.deeptrio.sorted.vcf.gz" >> merge_${LOW_SAMPLE_NAME}_deeptrio_vcfs.not_indel_realigned.sh

sbatch -c 8 --mem=50GB --time=4:00:00 --dependency=afterany:$swarm_jid merge_${LOW_SAMPLE_NAME}_deeptrio_vcfs.not_indel_realigned.sh
