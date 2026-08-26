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
my $job_name="GATK4.1_bam2fas_qc";
my $walltime="49:59:59";

# set up the paths, or else the current path will be taken as default settings
my $work_path_pl = getcwd();
my $file_path_pl = getcwd();

#please make sure there is a raw_bam folder under the path, script will take all the .bam files under this folder
my $bam_path_pl = "$file_path_pl/raw_bam";

# filename of script that will be generated (and optionally submitted)
# set up the name of qsub script here
my $script_path = "$work_path_pl/step1_run_GATK4.1_bam2cleanfq_fastp_qc_qsub_grch38.sh";

# Read directory specified by $bam_path; get the absolute path to all files ending in .bam
opendir(my $bam_dh, $bam_path_pl) or die "Could not open dirhandle to $bam_path_pl; died";
my @bam_files = ();

while ( defined (my $bam_rel = readdir($bam_dh)) ) {
	chomp $bam_rel;
	next unless ($bam_rel =~ /\.bam$/);
	my $bam_abs = "$bam_path_pl/$bam_rel";
	push(@bam_files, $bam_abs) if (-e $bam_abs and -f $bam_abs and -r $bam_abs);
}

closedir($bam_dh) or die "Could not close dirhandle to $bam_path_pl; died";

# Format filename and sanple name array elements to match Bash list syntax
my $bam_list_pl = join(' ', map {qq("$_")} @bam_files);

# Create Bash script for running bedtools genomecov
open(my $script_fh, ">", $script_path) or die "Could not write to $script_path: $!";

my $gcov_script = <<"END_SCRIPT";
#!/bin/bash

#PBS -N $job_name
#PBS -l nodes=1:ppn=8
#PBS -l gres=localhd:40
#PBS -l vmem=64g
#PBS -l walltime=$walltime
#PBS -t 0-$#bam_files
#PBS -d $work_path_pl

set -o pipefail
set -e

file_path=$file_path_pl
##the bam files should be put into a separate folder ./raw_bam
work_path=$work_path_pl
echo "\$PBS_O_WORKDIR"
echo "\$work_path"
cd "\$work_path" || exit 1


module load java/1.8.0_91
module load picard-tools/2.18.0
module load bwa/0.7.15
module load fastqc/0.11.5
module load fastp/v0.20.0

if [ ! -d "\$work_path/rawfq" ]; then
	mkdir -p "\$work_path/rawfq"
fi

#if [ ! -d \$work_path/fastqc_dir ]
#then mkdir -p \$work_path/fastqc_dir
#fi

if [ ! -d "\$work_path/cleanfq" ]; then
	mkdir -p "\$work_path/cleanfq"
fi

files=($bam_list_pl)

input_bam="\${files[\$PBS_ARRAYID]##*/}"
echo "\$input_bam"
if [[ "\$input_bam" =~ (.*).bam ]]
then
  echo "\$input_bam"
  num="\${BASH_REMATCH[1]}"
  #echo \$num
  java -jar /picard-tools/2.18.0/picard.jar SamToFastq I="$file_path_pl/raw_bam/\${input_bam}" F="\${work_path}/rawfq/\${num}.r1.fq.gz" F2="\${work_path}/rawfq/\${num}.r2.fq.gz" && echo "\$num bam2fastq done"
# fastqc \${work_path}/rawfq/*.fq -o \${work_path}/fastqc_dir/ && echo "\$num fastqc done"
  fastp --in1 "\${work_path}/rawfq/\${num}.r1.fq.gz" --in2 "\${work_path}/rawfq/\${num}.r2.fq.gz" --out1 "\${work_path}/cleanfq/\${num}.r1_paired.fq.gz" --out2 "\${work_path}/cleanfq/\${num}.r2_paired.fq.gz" --unpaired1 "\${work_path}/cleanfq/\${num}.r1_unpaired.fq.gz" --unpaired2 "\${work_path}/cleanfq/\${num}.r2_unpaired.fq.gz" -j "\${work_path}/cleanfq/\${num}.json" -h "\${work_path}/cleanfq/\${num}.html" > "\${work_path}/cleanfq/\${num}.report" && echo "\$num fastp reads qc done"
#  java -jar trimmomatic/0.36/trimmomatic-0.36.jar PE -threads 3 "\${work_path}/rawfq/\${num}.r1.fq.gz" "\${work_path}/rawfq/\${num}.r2.fq.gz" "\${work_path}/cleanfq/\${num}.r1_paired.fq.gz" "\${work_path}/cleanfq/\${num}.r1_unpaired.fq.gz" "\${work_path}/cleanfq/\${num}.r2_paired.fq.gz" "\${work_path}/cleanfq/\${num}.r2_unpaired.fq.gz" ILLUMINACLIP:/hpf/tools/centos6/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:5:15 MINLEN:50 && echo "\$num trimmomatic done"
fi

#rm -rf \${work_path}/rawfq
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
