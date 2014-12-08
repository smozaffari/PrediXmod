#!/usr/bin/perl
#use strict;
use warnings;
use List::Util qw[min max];

####This perl script takes the GTEx imputed vcf file as input, removes ambiguous-strand SNPs (A/T and C/G)
#### and makes several output files for each autosome for future parallel computing:
#### .mlinfo.gz and .mldose.gz MACH files for GCTA
#### .SNPxID matrix for quick scanning into R
#### .bim plink bim file with SNP pos info in .SNPxID

my $dir = "/nas40t2/gtex/GTEx_Analysis_2014-06-13/genotypes/OMNI_arrays/";
my $file = "GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv_allchr_genot_imput_info04_maf01_HEW1E6.vcf.gz";

system("zcat ${dir}${file} > tmp.vcf");

open(VCF, "tmp.vcf");
open(HAP, "/nas40t2/hwheeler/PrediXcan_CV/hapmapSnpsCEU.list"); ##from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/hapmapSnpsCEU.txt.gz

my %hapmapsnps;
while(<HAP>){
    chomp;
    my ($snp) = split(/\n/);
    $hapmapsnps{$snp} = 1;
}

#parse by chr
for(my $i = 1; $i <= 22; $i++){
    my $snpxidhandle = "SNPxID" . $i;
    my $mlinfohandle = "MLINFO" . $i;
    my $bimhandle = "BIM" . $i;
    open($snpxidhandle, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr${i}.SNPxID");
    open($bimhandle, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr${i}.bim");
    open($mlinfohandle, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr${i}.mlinfo");
}
open(INTRO, ">intro");


while(<VCF>){
    chomp;
    my ($chr, $pos, $rs, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/);
    my ($expfreq, $impinfo, $cert) = split(/;/,$info);
    my ($a, $freqalt) = split(/=/,$expfreq);
    my ($freqref) = 1 - $freqalt;
    my ($maf) = min($freqref,$freqalt);
    my ($b, $quality) = split(/=/,$cert);
    my ($c, $rsq) = split(/=/,$impinfo);
    if($chr eq "#CHROM"){
	for(my $i = 1; $i <= 22; $i++){
	    my $mlinfohandle = "MLINFO" . $i;
  	    print $mlinfohandle "SNP\tAl1\tAl2\tFreq1\tMAF\tQuality\tRsq\n";
	}
        foreach my $id (@genos){
	    my ($d,$e) = split(/-/,$id);
	    my $shortid = $d . '-' .$e;	    
	    print INTRO "$shortid->$shortid MLDOSE\n";
        }
    }
    if($ref eq "A" && $alt eq "T"){ ##rm potentially ambiguous strand SNPs
        next;
    }elsif($ref eq "T" && $alt eq "A"){
        next;
    }elsif($ref eq "C" && $alt eq "G"){
        next;
    }elsif($ref eq "G" && $alt eq "C"){
        next;
    }elsif(defined($hapmapsnps{$rs}) && $pos =~ m/\d+/){ ###only pull rsid SNPs at this time & don't print header rows
	my $snpxidhandle = "SNPxID" . $chr;		
	my $bimhandle = "BIM" . $chr;
	my $mlinfohandle = "MLINFO" . $chr;
	print $bimhandle "$chr\t$rs\t0\t$pos\t$ref\t$alt\n";
	print $mlinfohandle "$rs\t$ref\t$alt\t$freqref\t$maf\t$quality\t$rsq\n";
	foreach my $geno (@genos){
	    my ($probs, $gt, $dos) = split(/:/,$geno);
	    print $snpxidhandle "$dos\t";
	}
	print $snpxidhandle "\n";
    }
}

system("rm tmp.vcf");


for(my $i = 1; $i <= 22; $i++){
    open(R, ">runR.R") or die "cant mae runR.R\n";
    print R "dat<-read.table(\"GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr${i}.SNPxID\")\n";
    print R "dat<-t(dat)\n";
    print R "write.table(dat,\"t.dos.chr${i}\",col.names=F,row.names=F,quote=F)\n";
    close(R);
    system("R --vanilla < runR.R");
    system("paste -d\' \' intro t.dos.chr${i} > GTEx_Analysis_2014-06-13.hapmapSnpsCEU.chr${i}.mldose");
}

system("gzip *.mldose");
system("gzip *.mlinfo");
system("rm intro t.dos.chr*");
