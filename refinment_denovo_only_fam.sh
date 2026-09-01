#!/bin/bash
#latest GATK GATK_refinment pipeline updated in Oct30-2019
#Ting Liu

#PBS -N GATK_refinment
#PBS -l nodes=1:ppn=8
#PBS -l gres=localhd:80
#PBS -l mem=98g
#PBS -l vmem=98g
#PBS -l walltime=19:00:00

work_path=$PBS_O_WORKDIR
file_path=$PBS_O_WORKDIR
cd $work_path
#outname=combined_CMP
#file=$file_path/gvcf/*.HC.g.vcf.gz.tbi

outdir=$PBS_O_WORKDIR
#echo $file
file_predix="GATK4.0_combined_CMPs_22fam_denovo_gvcfs_chr1-22X.indels_snps.recalibrated.passed_chr1-22"
pedigree_file="fam.ped"
#predix of gzipped vcf file and pedigree file should be specified here

GATK_bundle=/references/homo_sapiens/v37_decoy/gatk
PICARD=/picard-tools/2.18.0
module load java/1.8.0_91
module load picard-tools/2.18.0
module load bwa/0.7.15
module load GATK/4.1.2.0

gatk --java-options "-Xmx68G" CalculateGenotypePosteriors \
-R /references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta \
--supporting /references/homo_sapiens/v37_decoy/gatk_bundle/1000G_phase3_v4_20130502.sites.vcf.gz \
-ped $pedigree_file \
-V $file_path/population/$file_predix.vcf.gz \
-O $file_predix.recaled_v4.1.vcf.gz && echo "CalculateGenotypePosteriors done **"

gatk --java-options "-Xmx68G" VariantFiltration \
-R /references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta \
-V $file_predix.recaled_v4.1.vcf.gz \
-G-filter "GQ < 20.0 || DP < 10" \
-G-filter-name lowGQ \
-O $file_predix.recaled_v4.1.vcf.gz.lowGQfiltered_v4.1.vcf.gz && echo "anoot lowGQ done **"

gatk --java-options "-Xmx68G" VariantAnnotator \
-R /references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta \
-V $file_predix.recaled_v4.1.vcf.gz.lowGQfiltered_v4.1.vcf.gz \
-A PossibleDeNovo \
-A InbreedingCoeff \
-ped $pedigree_file \
-O $file_predix.recaled_v4.1.vcf.gz.lowGQfiltered_v4.1.unfil_deNovos.vcf && echo "annot PossibleDeNovo done **"

awk '{if($1~/^#/||($1~/^chr/)&&($0~/hiConfDeNovo/)) print $0}' $file_predix.recaled_v4.1.vcf.gz.lowGQfiltered_v4.1.unfil_deNovos.vcf > $file_predix.recaled_v4.1.vcf.gz.lowGQfiltered_v4.1.final_awk_deNovo_variants.vcf && echo "to identify dovo variants is done **"
# rm -f $outdir/poplation/${outname}.HC.snps.VQSR.vcf.gz
# $outdir/poplation/${outname}.HC.*.recal
