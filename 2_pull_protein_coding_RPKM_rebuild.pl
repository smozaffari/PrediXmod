#!/usr/bin/perl
use strict;
use warnings;
use Compress::Zlib;


####This perl script takes the GTEx RNA-Seq *gene_rpkm.gxt.gz file as input, extracts protein coding genes,
#### and makes several output files:
#### .gene_rpkm.protein-coding.gct
#### .GENExID matrix for quick scanning into R
#### .ID.list colnames of matrix
#### .GENE.list rownames of matrix (Ensembl ID)
#### .GENEname.list matching gene names for Ensembl IDs


my $gtexdir = "/nas40t2/gtex/GTEx_Analysis_2014-06-13/rna-seq/";
my $file = "GTEx_Analysis_2014-06-13_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz";

my $dir = "/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/";

my $pcfile = ${dir} . "gencode.v18.genes.patched_contigs.summary.protein";
open(PCGENES, $pcfile); ##from /nas40t2/gtex/GTEx_Analysis_2014-06-13/reference_files/gencode.v18.genes.patched_contigs.gtf.gz, grep protein_coding

my %pcgenes;
while(<PCGENES>){
    chomp;
    my ($chr, $str, $start, $end, $ensid, $gene) = split(/\t/);
    $pcgenes{$ensid} = 1;
}

#outfiles for lasso:
open(RPKM, ">GTEx_Analysis_2014-06-13_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.protein-coding.gct");
open(ID, ">GTEx_Analysis_2014-06-13.RNA-seq.ID.list");
open(GENE, ">GTEx_Analysis_2014-06-13.RNA-seq.GENE.list");
open(GENExID, ">GTEx_Analysis_2014-06-13.RNA-seq.GENExID");
open(GENENAME, ">GTEx_Analysis_2014-06-13.RNA-seq.GENEname.list");

my $line;
my $rpkmfile = ${gtexdir} . ${file};
my $gzrpkm = gzopen($rpkmfile, "rb") or die "Cannot open $rpkmfile: $gzerrno\n";

while($gzrpkm->gzreadline($line)){
    my ($ensid, $gene, @rest) = split(/\t/, $line);
    if($ensid eq "Name"){
	my $rest = join("\t",@rest);
        print RPKM "$ensid\t$gene\t$rest\n";
	foreach my $id (@rest){
	    print ID "$id\n";
	}
    }
    elsif(defined($pcgenes{$ensid})){
	my $rest = join("\t",@rest);
	print RPKM "$ensid\t$gene\t$rest\n";
	print GENE "$ensid\n";
	print GENENAME "$ensid\t$gene\n";
	print GENExID "$rest\n";
    }
}	

$gzrpkm->gzclose();
system("gzip GTEx_Analysis_2014-06-13_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.protein-coding.gct");
