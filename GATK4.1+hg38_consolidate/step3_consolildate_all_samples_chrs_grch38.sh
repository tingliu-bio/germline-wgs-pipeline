#!/bin/bash

#PBS -N consoli_joint
#PBS -l nodes=1:ppn=8
#PBS -l gres=localhd:60
#PBS -l vmem=98g
#PBS -l walltime=200:00:00
#PBS -t 0-24

work_path=$PBS_O_WORKDIR
file_path=/bam2fq_noqc
cd $work_path
outname=geathervcfs
file=$file_path/gatk/*.HC.g.vcf.gz.tbi
outdir=$PBS_O_WORKDIR
#echo $file
chrs=("chrM" "chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

gatk_bundle=/GRCh38_Resources/GATK_bundle
PICARD=/picard-tools/2.18.0
module load java/1.8.0_91
module load picard-tools/2.18.0
module load bwa/0.7.15
module load GATK/4.1.2.0
module load samtools/1.9
module load R/3.5.1

if [ ! -d  $work_path/population ]
then mkdir -p $work_path/population
fi


echo "the gvcf files are"
for i in $file
do
  if [[ $i =~ (.*/.*.HC.g.vcf.gz).tbi ]] ;then
   #echo $i
   prefix=${BASH_REMATCH[1]}
   echo $prefix
   sample_gvcfs=${sample_gvcfs}"-V $prefix "
fi
done

date
echo "Running de novo analysis on ${chrs[$PBS_ARRAYID]} for $file"

#echo -e $sample_gvcfs
gatk GenomicsDBImport ${sample_gvcfs} --genomicsdb-workspace-path ${chrs[$PBS_ARRAYID]} -L $work_path/intervals/${chrs[$PBS_ARRAYID]}.intervals && echo "** ${chrs[$PBS_ARRAYID]}.GenomicsDBImport done **"

gatk GenotypeGVCFs -R $gatk_bundle/Homo_sapiens_assembly38.fasta \
-V gendb://${chrs[$PBS_ARRAYID]} \
-O $outdir/population/${chrs[$PBS_ARRAYID]}.HC.vcf.gz && echo "** ${chrs[$PBS_ARRAYID]}.HC.vcf.gz done **"
