#!/bin/bash

#PBS -N vcf.4.impute
#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -o $HOME/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e $HOME/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR

$DIR=/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/cis.v.trans.prediction/DGN-WB_genotypes/

####Convert to vcf for UM Imputationserver
##see https://imputationserver.sph.umich.edu/start.html#!pages/help

module load plink/1.09
module load vcftools
module load tabix/0.2.6 

for (( i = 1 ; i <= 22; i++))
do
    plink --bfile ${DIR}DGN.hapmap2.chr1-22.QC --chr $i --recode vcf --out ${DIR}DGN.hapmap2.QC.test.chr$i
    vcf-sort ${DIR}DGN.hapmap2.QC.test.chr$i.vcf | bgzip -c > ${DIR}DGN.hapmap2.QC.test.chr$i.vcf.gz
done

##use sftp to upload *vcf.gz files to https://imputationserver.sph.umich.edu/start.html#!pages/run, see 'vcf.filelist.for.sftp' for paths
##started QC run for all chromosomes, expecting strand flip errors
##download statistics.txt file from Imputationserver and renamed 'statistics.DGN.chr1-22.txt' to get errors
##To fix strand flips:

grep FILTER ${DIR}statistics.DGN.chr1-22.txt> ${DIR}test
grep -o rs[0-9]* ${DIR}test > ${DIR}DGN.chr1-22.strand.switches

for (( i = 1 ; i <= 22; i++))
do
    plink --bfile ${DIR}DGN.hapmap2.chr1-22.QC --chr $i --flip ${DIR}DGN.chr1-22.strand.switches --recode vcf --out ${DIR}DGN.hapmap2.QC.chr$i.strand.switches
    vcf-sort ${DIR}DGN.hapmap2.QC.chr$i.strand.switches.vcf | bgzip -c > ${DIR}DGN.hapmap2.QC.chr$i.strand.switches.vcf.gz
done

##use sftp to upload *strand.switches.vcf.gz files to https://imputationserver.sph.umich.edu/start.html#!pages/run, see 'vcf.filelist.for.sftp' for paths 
##save output in /group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/transfers/PrediXmod/DGN-WB/DGN-imputation/UMich-imputation-results/