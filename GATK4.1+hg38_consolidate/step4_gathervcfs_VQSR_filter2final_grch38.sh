#!/bin/bash

#PBS -N GatherVcfs_vqsr2final
#PBS -l nodes=1:ppn=8
#PBS -l gres=localhd:60
#PBS -l vmem=98g
#PBS -l walltime=20:00:00

#assign an prefix for outputfile
prefix=afterVQSR_final_recaled

#set up the work,file and out paths
work_path=$PBS_O_WORKDIR
cd $work_path
file_path=$PBS_O_WORKDIR
file=$file_path/gatk/*.HC.g.vcf.gz.tbi
outdir=$PBS_O_WORKDIR
#echo $file

gatk_bundle="/GRCh38_Resources/GATK_bundle"
PICARD=/picard-tools/2.18.0
module load java/1.8.0_91
module load picard-tools/2.18.0
module load bwa/0.7.15
module load GATK/4.1.2.0
module load samtools/1.9
module load R/3.6.1

if [ ! -d  $work_path/population ]
then mkdir -p $work_path/population
fi

chrs=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" "chrM")

for i in ${chrs[@]}
do
   chr_vcfs=${chr_vcfs}"-I $work_path/population/$i.HC.vcf.gz "
done

#first gather up the vcfs from all chroms, then index the vcf file
echo -e $chr_vcfs
gatk GatherVcfs ${chr_vcfs} -O $work_path/population/$prefix.HC.vcf.gz && echo "** $prefix.HC.vcf.gz done **"
gatk IndexFeatureFile -F $work_path/population/$prefix.HC.vcf.gz  && echo "** $prefix.HC.vcf.gz index done **"

####then perferm the two-step VQSR for SNPs an Indels separately####
##first the snp mode
gatk VariantRecalibrator \
-R "${gatk_bundle}/Homo_sapiens_assembly38.fasta" \
-V "${work_path}/population/${prefix}.HC.vcf.gz" \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 "${gatk_bundle}/hapmap_3.3.hg38.vcf.gz" \
-resource:omni,known=false,training=true,truth=true,prior=12.0 "${gatk_bundle}/1000G_omni2.5.hg38.vcf.gz" \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 "${gatk_bundle}/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
-resource:dbsnp,known=true,training=false,truth=false,prior=2 "${gatk_bundle}/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz" \
-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
--trust-all-polymorphic \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
-mode SNP \
--max-gaussians 6 \
--rscript-file "${outdir}/population/${prefix}.HC.snps.plots.R" \
--tranches-file "${outdir}/population/${prefix}.HC.snps.tranches" \
-O "${outdir}/population/${prefix}.HC.snps.recal" && echo "** SNPs $prefix VariantRecalibrator done **"

gatk ApplyVQSR \
-R "${gatk_bundle}/Homo_sapiens_assembly38.fasta" \
-V "${work_path}/population/$prefix.HC.vcf.gz" \
--recal-file "${outdir}/population/${prefix}.HC.snps.recal" \
--tranches-file "${outdir}/population/${prefix}.HC.snps.tranches" \
--truth-sensitivity-filter-level 90.0 \
--create-output-variant-index true \
-mode SNP \
-O "${outdir}/population/${prefix}.snps.recalibrated.vcf.gz" && echo "** $prefix SNPs ApplyVQSR done **"

## then the Indel mode
gatk VariantRecalibrator \
-R "${gatk_bundle}/Homo_sapiens_assembly38.fasta" \
-V "${work_path}/population/${prefix}.HC.vcf.gz" \
--trust-all-polymorphic \
-resource:mills,known=false,training=true,truth=true,prior=12 "${gatk_bundle}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
-resource:dbsnp,known=true,training=false,truth=false,prior=2 "${gatk_bundle}/beta/Homo_sapiens_assembly38.known_indels.vcf.gz" \
-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
-mode INDEL \
--max-gaussians 4 \
--rscript-file "${outdir}/population/${prefix}.HC.indels.plots.R" \
--tranches-file "${outdir}/population/${prefix}.HC.indels.tranches" \
-O "${outdir}/population/${prefix}.HC.indels.recal" && echo "** Indels $prefix VariantRecalibrator done **"

gatk ApplyVQSR \
-R "${gatk_bundle}/Homo_sapiens_assembly38.fasta" \
-V "${outdir}/population/${prefix}.snps.recalibrated.vcf.gz" \
--recal-file "${outdir}/population/${prefix}.HC.indels.recal" \
--truth-sensitivity-filter-level 90.0 \
--tranches-file "${outdir}/population/${prefix}.HC.indels.tranches" \
--create-output-variant-index true \
-mode INDEL \
-O "${outdir}/population/${prefix}.indels_snps.recalibrated.vcf.gz" && echo "** $prefix SNPs and Indels ApplyVQSR done **"

gatk SelectVariants \
-R "${gatk_bundle}/Homo_sapiens_assembly38.fasta" \
-V "${outdir}/population/${prefix}.indels_snps.recalibrated.vcf.gz" \
--exclude-filtered \
-O "${outdir}/population/${prefix}.indels_snps.recalibrated.passed.vcf.gz" && echo "** $prefix pipeline done **"


