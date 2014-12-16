#!/usr/bin/perl
use strict;
use warnings;


my $dir = "/group/im-lab/nas40t2/haky/Data/dbGaP/GTEx/41400/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2014-06-13/genotypes/OMNI_arrays/";
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

#outfiles for lasso:
open(DOS, ">GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv.hapmapSnpsCEU.unamb.dosage");
open(ID, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.ID.list");
open(SNP, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.SNP.list");
open(SNPxID, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.SNPxID");
open(BIM, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.bim");

while(<VCF>){
    chomp;
    my ($chr, $pos, $rs, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/);
    if($chr eq "#CHROM"){
        foreach my $id (@genos){
            print ID "$id\n";
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
	print DOS "$chr\t$rs\t$pos\t$ref\t$alt\t";
	print BIM "$chr\t$rs\t0\t$pos\t$ref\t$alt\n";
	print SNP "$rs\n";
	foreach my $geno (@genos){
	    my ($probs, $gt, $dos) = split(/:/,$geno);
	    print DOS "$dos\t";
	    print SNPxID "$dos\t";
	}
	print DOS "\n";
	print SNPxID "\n";
    }
}

system("rm tmp.vcf");

system("gzip GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv.hapmapSnpsCEU.unamb.dosage");
