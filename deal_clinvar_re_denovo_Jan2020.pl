use strict;
use warnings;
#dealing with the output ouf GATK_refinment modules, filter out the lowGQ samples one de novo mutation for each line
##take in the output of intervar, revert the genomic coordinates into vcf style
#take in sample file which named sample.out, genrate by the script print_samples.pl, put the sample id in the out file

my %all_vcf=();
open (IN,"final_denovo_intervar_out.hg38_multianno.txt.intervar") or die "cannot open vcf file";

my %chd_genes=();
open(CHD,"/genelists/Gene_list_reclassification_20190409.txt");
while(<CHD>){
  chomp;
  my($gene,$catg)=split("\t");
  if(!exists $chd_genes{$gene}){
    $chd_genes{$gene}=$catg;
  }#if
  else{
    print "repeated id $gene error\n";
  }
}#while
my %cmp_genes=();
open(CMP,"/genelists/genelist.txt");
while(<CMP>){
  chomp;
  my($gene,$catg)=split("\t");
  if(!exists $cmp_genes{$gene}){
    $cmp_genes{$gene}=$catg;
  }#if
  else{
    print "repeated id $gene error\n";
  }
}#while
my $cmp_gene_tag=".";
my $chd_gene_tag=".";

open (OUT,">re_mod_clinvar_denovo.out");
print OUT "chr\tstart\tref\talt\ttype\tref_gene\tfunc_g\texoic_Gne\treformed_annot\tintv_annot\tens_gene\tsnpid\taa_ens\taa_refg\tgnomAD_freqall\tfreq_esp\tfreq_1kg_all\tCAD_raw\tCAD_phred\tSIFT\tGERP\n";
my @head_id=();
while(<IN>){
  next if(/^##/); #ignore header
  chomp;
  if (/^#Chr/){ #if chrom, save the sample id
    @head_id=split("\t");
   }#if
  else{
  my ($chr, $start, $end, $ref, $alt, $ref_gene,$func_g,$exoic_Gne,$ens_gene,$snpid,$aa_ens,$aa_refg,$clinvar,$intervar,$gnomAD_freqall,$freq_esp,$freq_1kgall,$CAD_raw,$CAD_phred,$sift,$gerp) = split("\t");
  my @line= split("\t");
#  print OUT "$intervar\n";

   my  @interv= split(":",$intervar);
   my ($intv_annot)= ($interv[1]=~/(.*)PVS1/);
#   print OUT "$intv_annot\t$intervar\n";

   my @INFO=split(" ",$interv[1]);
#   print OUT "$INFO[0]\n";
   my (%INFO_item) =();
   for (my $i= 0; $i<= $#INFO-1; $i++){ #parsing the INFO item of each ID, storing each ID=desc to a pair of hash value
      my $info_id='';
      my $info_desc='';
     if (($info_id,$info_desc) = ($INFO[$i] =~ /(.*)\=(.*)/)){
     if(exists $INFO_item{$info_id}){
      print "error\n";
       }#if info
     else{
         $INFO_item{$info_id}=$info_desc;
              }#else
          }#if
         }#for INFO
# my $line_index={'ref'=>$ref,'alt'=>$alt,'AC'=>$INFO_item{'AC'},'sample'=> $INFO_item{'hiConfDeNovo'}};
# print "$chr\t$pos\t$ref\t$alt\tAC=$INFO_item{'AC'}\tsample=$INFO_item{'hiConfDeNovo'}\n";
# print "$chr\t$start\t$end\t$ref\t$alt\t";
# print "AC_before=$INFO_item{'AC'}\n";
# if  ($INFO_item{'AC'} ne 1) {  #AC>=1, then count genotype
 #print "$interv[1]\t";
     my ($ps_evidences) = ($interv[1] =~ /PS\=\[(.*?)\]/);
     my ($ps_crit)= &check_evd($ps_evidences);
     $ps_crit=$ps_crit+1;  # only works for de novo list
     my ($pm_evidences) = ($interv[1] =~ /PM\=\[(.*?)\]/);
     my ($pm_crit)= &check_evd($pm_evidences);
     my ($pp_evidences) = ($interv[1] =~ /PP\=\[(.*?)\]/);
     my ($pp_crit)= &check_evd($pp_evidences);
     my ($bs_evidences) = ($interv[1] =~ /BS\=\[(.*?)\]/);
     my ($bs_crit)= &check_evd($bs_evidences);
     my ($bp_evidences) = ($interv[1] =~ /BP\=\[(.*?)\]/);
     my ($bp_crit)= &check_evd($bp_evidences);
  my $tag="VUS";
  print "pvs1=$INFO_item{'PVS1'}\tba1=$INFO_item{'BA1'}\tps_evi=$ps_evidences\tps_crit=$ps_crit\tpm_evi=$pm_evidences\tpm_crit=$pm_crit\tpp_evi=$pp_evidences\tpp_crit=$pp_crit\n";
  if (($INFO_item{'PVS1'}==1)&&(($ps_crit>=1)||($pm_crit>=2)||(($pm_crit==1)&&($pp_crit==1))||($pp_crit>=2))){
     $tag="Pathogenic";
    }#if i
  elsif($ps_crit>=2){
    $tag="Pathogenic";
    }#elsif ii
  elsif(($ps_crit==1)&&(($pm_crit>=3)||(($pm_crit==2)&&($pp_crit>=2))||(($pm_crit==1)&&($pp_crit>=4)))) {
   $tag="Pathogenic";
     }#elsif iii
###Likely_pathogenic######################
  elsif (($INFO_item{'PVS1'}==1)&&($pm_crit==1)){
    $tag="Likely_pathogenic";
  }#elsif likely i
  elsif(($ps_crit==1)&&(($pm_crit==1)||($pm_crit==2))) {
   $tag="Likely_pathogenic";
 }#elsif likely ii
  elsif (($ps_crit==1)&&($pp_crit>=2)) {
    $tag="Likely_pathogenic";
  }#elsif likely iii
   elsif ($pm_crit>=3) {
      $tag="Likely_pathogenic";
    }#elsif likely iv
  elsif (($pm_crit==2)&&($pp_crit>=2)) {
      $tag="Likely_pathogenic";
    }#elsif likely v
  elsif (($pm_crit==1)&&($pp_crit>=4)) {
      $tag="Likely_pathogenic";
    }#elsif likely vi
###########benign######################
  elsif ($INFO_item{'BA1'}==1) {
     $tag="Benign";
   }#elsif Benign i
  elsif ($bs_crit>=2) {
    $tag="Benign";
  }#elsif Benign ii
###Likely_benign######################
 elsif (($bs_crit==1)&&($bp_crit==1)) {
     $tag="Likely_benign";
   }#elsif likely benign i
 elsif ($bp_crit>=2) {
     $tag="Likely_benign";
   }#elsif likely benign ii
######################################
#print OUT "$chr\t$start\t$ref\t$alt\t$ref_gene\t$func_g\t$exoic_Gne\t$intervar\n";
#print OUT "$intv_annot\t$tag\n";
my $var_tag="INDEL";
if ((length($ref)==1)&&(length($alt)==1)&&($ref ne "-")&&($alt ne "-")){
  $var_tag="SNP";
  }

my $chd_out=$ref_gene;
if(exists $chd_genes{$ref_gene}){
  $chd_out=$ref_gene."_".$chd_genes{$ref_gene};
  }#if
my $cmp_out=$ref_gene;
 if(exists $cmp_genes{$ref_gene}){
   $cmp_out=$ref_gene."_".$cmp_genes{$ref_gene};
#   print CMP "chr$chr\t$start\t$ref\t$alt\t$var_tag\t$ref_gene\t$func_g\t$exoic_Gne\t$tag\t$intv_annot\t$ens_gene\t$snpid\t$aa_ens\t$aa_refg\n";
  }#if
   print OUT "chr$chr\t$start\t$ref\t$alt\t$var_tag\t$ref_gene\t$func_g\t$exoic_Gne\t$tag\t$intv_annot\t$ens_gene\t$snpid\t$aa_ens\t$aa_refg\t$gnomAD_freqall\t$freq_esp\t$freq_1kgall\t$CAD_raw\t$CAD_phred\t$sift\t$gerp\n";

  }#else not header
}#while

sub check_evd {
 my ($evd_str)=@_;
# print OUT "$evd_str\n";
 my $evd_num=0;
my @item=split(", ",$evd_str);
#print OUT "$evd_num\t$item[1]\n";
 for(my $i=0 ;$i<=$#item; $i++) {
  if($item[$i]==1){
    $evd_num++;
    }#if ne 0
  }#for
return $evd_num;
}#sub
