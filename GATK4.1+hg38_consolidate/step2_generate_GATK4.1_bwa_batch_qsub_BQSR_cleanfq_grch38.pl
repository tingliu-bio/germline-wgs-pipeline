#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Basename;
use Cwd;

# set to 1 for test print of script without automated qsub
# set to 0 for automated qsub submission after script generation
my $noqsub = 1;

##setup the name and expected walltime of your job here
my $job_name="GATK4.1_grch38_map2BQSR";
my $walltime="99:59:59";

# set up the paths, or else the current path will be taken as default settings
my $work_path_pl = getcwd();
my $file_path_pl = getcwd();

##check if there is a cleanfq folder from the step1
my $fastq_path_pl = "$file_path_pl/cleanfq";

# Filename of script that will be generated (and optionally submitted)
my $script_path = "./step2_run_bwa_GATK4.1_BQSR_callgvcf4sample_qsub_grch38.sh";

# Read directory specified by $bam_path; get the absolute path to all files ending in .bam
opendir(my $fastq_dh, $fastq_path_pl) or die "Could not open dirhandle to $fastq_path_pl; died";
my @fastq_files = ();

while ( defined (my $fastq_rel = readdir($fastq_dh)) ) {
	chomp $fastq_rel;
	next unless ($fastq_rel =~ /\.r1_paired\.fq\.gz$/);
	my $fastq_abs = "$fastq_path_pl/$fastq_rel";
	push(@fastq_files, $fastq_abs) if (-e $fastq_abs and -f $fastq_abs and -r $fastq_abs);
}

closedir($fastq_dh) or die "Could not close dirhandle to $fastq_path_pl; died";

# Format filename and sanple name array elements to match Bash list syntax
my $fastq_list_pl = join(' ', map {qq("$_")} @fastq_files);

# Create Bash script for running bedtools genomecov
open(my $script_fh, ">", $script_path) or die "Could not write to $script_path: $!";

my $gcov_script = <<"END_SCRIPT";
#!/bin/bash

#PBS -N $job_name
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:40
#PBS -l vmem=64g
#PBS -l walltime=$walltime
#PBS -t 0-$#fastq_files
#PBS -d $work_path_pl

set -o pipefail
set -e

##path to the clean fq data, then read the filename for r1 only
##use the r1 prefix as the output prefix
##please make sure r1 and r2 data are all included in the cleanfq directory
echo "pbspath \$PBS_O_WORKDIR"
#work_path=/tingliu/GATK_testrun3_qsubscript
work_path="\$PBS_O_WORKDIR"
#file=\$file_path/cleanfq/*.r1_paired.fq.gz
gatk_bundle="/GRCh38_Resources/GATK_bundle"

module load java/1.8.0_91
module load picard-tools/2.18.0
module load bwa/0.7.15
module load GATK/4.1.2.0
module load samtools/1.9

##find the cleanreads directory
if [ ! -d "\${work_path}/cleanfq" ]; then
	echo "ERROR: No clean reads folder"
	exit 1
fi

##find the raw_bam directory to get the header line(\@RG) of the raw_bam,rewrite it into the remapped bam files
if [ ! -d "\${work_path}/raw_bam" ]; then
	echo "ERROR: No raw_bam folder"
	exit 1
fi

##output directory
if [ ! -d "\${work_path}/bwa" ]; then
	mkdir -p "\${work_path}/bwa"
fi

if [ ! -d "\${work_path}/gatk" ]; then
	mkdir -p "\${work_path}/gatk"
fi

if [ ! -d "\${work_path}/mito_data" ]; then
	mkdir -p "\${work_path}/mito_data"
fi

files=($fastq_list_pl)

input_fastq="\${files[\$PBS_ARRAYID]##*/}"

###traverse the ./cleanfq file folder, then get the file prefix from the r1 reads
###the script will not check if the r2 reads exists or not, please make sure the paired reads is placed in the same ./cleanfq folder

if [[ "\$input_fastq" =~ (.*).r1_paired.fq.gz ]]
 then
	echo "\$input_fastq"
	prefix="\${BASH_REMATCH[1]}"
	echo "\$prefix"

##get the header first from the raw_bam files
reorg_rg=\$(samtools view -H "\${work_path}/raw_bam/\${prefix}.bam" | grep '\@RG' | sed 's#\\t#\\\\t#g'| head -1)

##use bwa men to re-aln the reads
time bwa mem -t 4 -M -R "\$reorg_rg" "\${gatk_bundle}/Homo_sapiens_assembly38.fasta" "\${work_path}/cleanfq/\${prefix}.r1_paired.fq.gz" "\${work_path}/cleanfq/\${prefix}.r2_paired.fq.gz" | samtools view -Sb - > "\${work_path}/bwa/\${prefix}.clean.bam" && echo "** \$prefix bwa mem done **"

time samtools sort -@ 4 -m 8G -O bam -o "\${work_path}/bwa/\${prefix}.clean.sorted.bam" "\${work_path}/bwa/\${prefix}.clean.bam" && echo "** sorted clean \$prefix bamfile done **"

gatk MarkDuplicates \\
-I "\${work_path}/bwa/\${prefix}.clean.sorted.bam" \\
-M "\${work_path}/bwa/\${prefix}.clean.markdup_matrics.txt" \\
-O "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.bam" && echo "** \$prefix clean markdup done **"

gatk BaseRecalibrator \\
-R "\${gatk_bundle}/Homo_sapiens_assembly38.fasta" \\
-I "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.bam" \\
--known-sites "\${gatk_bundle}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \\
--known-sites "\${gatk_bundle}/beta/Homo_sapiens_assembly38.known_indels.vcf.gz" \\
--known-sites "\${gatk_bundle}/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz" \\
--known-sites "\${gatk_bundle}/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \\
-O "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.recal_data.table" && echo "** \$prefix recal_data.table done **"

gatk ApplyBQSR \\
-bqsr-recal-file "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.recal_data.table" \\
-R "\${gatk_bundle}/Homo_sapiens_assembly38.fasta" \\
-I "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.bam" \\
-O "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.BQSR.bam" && echo "** ApplyBQSR \$prefix done **"

samtools index "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.BQSR.bam" && echo "** \$prefix sorted.markdup.BQSR.bam index done **"

gatk Mutect2 \\
-R "\${gatk_bundle}/Homo_sapiens_assembly38.fasta" \\
-L chrM \\
--mitochondria-mode \\
-I "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.BQSR.bam" \\
-O "\${work_path}/mito_data/\${prefix}.clean.sorted.HC.g.chrM.vcf.gz" && echo "** \$prefix chrM HC.g.vcf.gz done **"

gatk HaplotypeCaller --emit-ref-confidence GVCF \\
-R "\${gatk_bundle}/Homo_sapiens_assembly38.fasta" \\
-I "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.BQSR.bam" \\
-O "\${work_path}/gatk/\${prefix}.clean.sorted.HC.g.vcf.gz" && echo "** \$prefix HC.g.vcf.gz done **"

##option split the chrom by each chrom
#chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM )
#for i in \${chrom[@]}
#do
#echo "\${i}"
#gatk HaplotypeCaller --emit-ref-confidence GVCF \\
#-R "/references/homo_sapiens/v37_decoy/gatk/human_g1k_v37_decoy_mod.fasta" \\
#-I "\${work_path}/gatk/\${prefix}.clean.sorted.markdup.BQSR.bam" \\
#-L "\${i}" \\
#-O "\${work_path}/split_chrs/\${prefix}.HC.\${i}.g.vcf.gz" && echo "** \$prefix.HC.HC.\${i}.g.vcf.gz done **"
#done

fi

END_SCRIPT

print $script_fh "$gcov_script";

close($script_fh) or die "Could not close filehandle to $script_path; $!";

# qsub the script that was just made and capture the job ID assigned to it
# This is done as late as possible in case the "--noqsub" option is set, in order to allow the scripts to be printed out before quitting
exit 0 if ($noqsub);
my $qsub_job_id = `qsub $script_path`;
die "Failed to qsub job; died" unless ($qsub_job_id and $qsub_job_id =~ /^[0-9]+(?:\[\])?/);
chomp $qsub_job_id;

print "Submitted job with job id $qsub_job_id\n";
