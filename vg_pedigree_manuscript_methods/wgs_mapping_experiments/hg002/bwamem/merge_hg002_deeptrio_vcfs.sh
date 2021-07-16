#!/bin/bash
cd /data/Udpwork/usr/markellocj/non_segdup_vg_pedigree_runs/hg002_cohort/bwamem/;
bcftools concat -O v HG002.bwamem.deeptrio.chr{?,??}.abra_indel_realigned.vcf.gz > merged.deeptrio.indel_realigned.vcf && bcftools sort merged.deeptrio.indel_realigned.vcf -O v > merged.deeptrio.indel_realigned.sorted.vcf && rm merged.deeptrio.indel_realigned.vcf && bgzip merged.deeptrio.indel_realigned.sorted.vcf && tabix -p vcf merged.deeptrio.indel_realigned.sorted.vcf.gz
