#!/usr/bin/perl
use strict;
use warnings;

####This perl script takes the GTEx imputed vcf file as input, removes ambiguous-strand SNPs (A/T and C/G)
#### and makes several output files:
#### .dosage.gz file for PrediXcan (compute.scores.py)
#### .mlinfo.gz and .mldose.gz MACH files for GCTA
#### .SNPxID matrix for quick scanning into R
#### .ID.list colnames of matrix
#### .SNP.list rownames of matrix
#### .bim plink bim file with SNP pos info

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

#outfiles for downstream applications:
open(DOS, ">GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv.hapmapSnpsCEU.unamb.dosage");
open(ID, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.ID.list");
open(SNP, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.SNP.list");
open(SNPxID, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.SNPxID");
open(BIM, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.bim");
open(INTRO, ">intro");
open(MLINFO, ">GTEx_Analysis_2014-06-13.hapmapSnpsCEU.mlinfo");

while(<VCF>){
    chomp;
    my ($chr, $pos, $rs, $ref, $alt, $qual, $filter, $info, $format, @genos) = split(/\t/);
    my ($expfreq, $impinfo, $cert) = split(/;/,$info);
    my ($a, $freqalt) = split(/=/,$expfreq);
    my ($freqref) = 1 - $freqalt;
    my ($b, $quality) = split(/=/,$cert);
    my ($c, $rsq) = split(/=/,$impinfo);
    if($chr eq "#CHROM"){
        foreach my $id (@genos){
	    my ($d,$e) = split(/-/,$id);
	    my $shortid = $d . '-' .$e;	    
            print ID "$shortid\n";
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
	print DOS "$chr\t$rs\t$pos\t$ref\t$alt\t";
	print BIM "$chr\t$rs\t0\t$pos\t$ref\t$alt\n";
	print SNP "$rs\n";
	print MLINFO "$rs\t$ref\t$alt\t$freqref\t$quality\t$rsq\n";
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


open(R, ">runR.R") or die "cant mae runR.R\n";
print R "dat<-read.table(\"GTEx_Analysis_2014-06-13.hapmapSnpsCEU.SNPxID\")\n";
print R "dat<-t(dat)\n";
print R "write.table(dat,\"t.dos\",col.names=F,row.names=F,quote=F)\n";
close(R);
system("R --vanilla < runR.R");
system("paste -d\' \' intro t.dos > GTEx_Analysis_2014-06-13.hapmapSnpsCEU.mldose");


system("gzip GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv.hapmapSnpsCEU.unamb.dosage");
system("gzip GTEx_Analysis_2014-06-13.hapmapSnpsCEU.mldose");
system("gzip GTEx_Analysis_2014-06-13.hapmapSnpsCEU.mlinfo");
system("rm intro t.dos");
