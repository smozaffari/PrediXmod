####by Heather E. Wheeler 20140602####
"%&%" = function(a,b) paste(a,b,sep="")
args <- commandArgs(trailingOnly=T)
date <- Sys.Date() 
###############################################
### Directories & Variables
 
cri = F #on cri cluster?
if(cri) precri = "/group/dolan-lab/hwheeler/" else precri = '/'

my.dir <- precri %&% "nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/"
 
#tissue <- "Nerve - Tibial" ###check GTEx_Analysis_2014-06-13.SampleTissue.annot for available tissues###
tissue <- "Adipose - Subcutaneous"
tis <- "GTEx-AdS"
Nk <- 15 ##number of peer factors to calculate, recommend 25% of sample size, but no more than 100, GTEx included 15 in pilot analyses
 
################################################
### Functions & Libraries
 
library(SNPRelate)
library(peer)
library(preprocessCore)
#library(GenABEL)
##if can't install GenABEL
source(my.dir %&% 'GenABEL/R/ztransform.R')
source(my.dir %&% 'GenABEL/R/rntransform.R')
 
################################################
sam <- read.table(my.dir %&% "GTEx_Analysis_2014-06-13.SampleTissue.annot",header=T,sep="\t") 
### above file includes SAMPID and SMTSD from: /nas40t2/gtex/GTEx_Analysis_2014-06-13/sample_annotations/GTEx_Data_2014-06-13_Annotations_SampleAttributesDS.txt
sample <- subset(sam,SMTSD == tissue) ### pull sample list of chosen tissue###
 
expidlist <- scan("GTEx_Analysis_2014-06-13.RNA-seq.ID.list","character")
expgenelist <- scan("GTEx_Analysis_2014-06-13.RNA-seq.GENE.list","character")
exp <- scan("GTEx_Analysis_2014-06-13.RNA-seq.GENExID")
expdata <- matrix(exp, ncol=length(expidlist), byrow=T)
t.expdata <- t(expdata)
rownames(t.expdata) <- expidlist
colnames(t.expdata) <- expgenelist

gencodefile <- my.dir %&% "gencode.v18.genes.patched_contigs.summary.protein"
gencode <- read.table(gencodefile) ##split into 10 files, call each from run_1_CV_GTEx_polyscore_PrediXcan_subset*.sh
rownames(gencode) <- gencode[,5]
t.expdata <- t.expdata[,intersect(colnames(t.expdata),rownames(gencode))] ###pull protein coding gene expression data

tissue.exp <- t.expdata[intersect(rownames(t.expdata),sample$SAMPID),] ###pull expression data for chosen tissue###
expsamplelist <- rownames(tissue.exp) ###samples with exp data###
substr.expsamplelist <- substr(expsamplelist,1,10) ###to match with genotype data###
 
famfile <- my.dir %&% "GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv_Pheno_2for5M_1for2.5M.fam"
fam <- read.table(famfile)
gtsamplelist <- fam$V1
substr.gtsamplelist <- substr(gtsamplelist,1,10) ###to match with exp data###
rownames(fam) <- substr.gtsamplelist
samplelist <- intersect(substr.gtsamplelist,substr.expsamplelist)
nsample <- length(samplelist)
tissue.exp.substr <- tissue.exp
for(id in samplelist){  ####take mean of exp for samples with >1 RNA Seq dataset
        matchexp <- tissue.exp.substr[substr(rownames(tissue.exp.substr),1,10)==id,]
        if(is.array(matchexp)=='TRUE'){
                expmean <- colMeans(matchexp)
                for(i in rownames(matchexp)){
                        tissue.exp.substr[i,] <- expmean
                }
        }
}


rownames(tissue.exp.substr) <- substr.expsamplelist ###change rownames of tissue.exp to match with genotypes###
exp.w.geno <- tissue.exp.substr[samplelist,] ###get expression of samples with genotypes###
explist <- subset(colMeans(exp.w.geno), colMeans(exp.w.geno)>0) ###pull genes with mean expression > 0###
explist <- names(explist)
exp.w.geno <- exp.w.geno[,explist]
 

###get first 3 PCs from genos, as in GTEx
pc.matrix<-read.table(my.dir %&% "GTEx_Analysis_2014-06-13_OMNI_2.5M_5M_451Indiv_PostImput_20genotPCs.txt",header=T)
pcs3 <- pc.matrix[,3:5]
pclist <- substr(pc.matrix[,1],1,10)
rownames(pcs3) <- pclist
pcs <- as.matrix(pcs3[samplelist,])

###pull gender, used as cov in GTEx
famsubset <- fam[samplelist,]
gender <- famsubset$V5
names(gender)<-rownames(famsubset)

###quantile normalize and transform to standard normal exp.w.geno matrix, as in GTEx###
t.exp.w.geno <- t(exp.w.geno)

rowtable<-function(x) length(table(x))>2 ##function to determine if more than 2 exp levels per gene
nonbin<-apply(t.exp.w.geno,1,rowtable) ##apply to matrix
t.exp.w.geno <- t.exp.w.geno[nonbin,] ##remove binary genes from matrix

qn.t.exp.w.geno <- normalize.quantiles(t.exp.w.geno) ##quantile normalize
rn.qn.t.exp.w.geno <- apply(qn.t.exp.w.geno,1,"rntransform") ##rank transform to normality & transposes, not sure why?##

###Now we can create the model object, ### from https://github.com/PMBio/peer/wiki/Tutorial

model = PEER()

###set the observed data,

PEER_setPhenoMean(model,as.matrix(rn.qn.t.exp.w.geno))

dim(PEER_getPhenoMean(model))

###(NULL response means no error here), say we want to infer K=20 hidden confounders,

PEER_setNk(model,Nk)

PEER_getNk(model)

####and perform the inference. ###for Nk=20 and GTEx-NT, it took 323 iterations, for Nk=15 and GTEx-NT, it took 37 iterations

PEER_update(model)

factors = PEER_getX(model)
rownames(factors) <- rownames(exp.w.geno)
write.table(factors,file <- tis %&% "." %&% Nk %&% ".PEER.factors." %&% date %&% ".txt", quote=F)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

png(file= tis %&% "." %&% Nk %&% ".PEER.factors.plotmodel." %&% date %&% ".png")
PEER_plotModel(model)
dev.off()

adj.exp.matrix<-matrix(NA,nrow=dim(rn.qn.t.exp.w.geno)[1],ncol=dim(rn.qn.t.exp.w.geno)[2])

for(i in 1:dim(rn.qn.t.exp.w.geno)[2]){
	res <- summary(lm(rn.qn.t.exp.w.geno[,i] ~ factors + pcs + gender, na.action=na.exclude))
	resid <- residuals(res)
	adj.exp.matrix[,i] <- resid
}

colnames(adj.exp.matrix) <- rownames(t.exp.w.geno)
rownames(adj.exp.matrix) <- colnames(t.exp.w.geno)

write.table(adj.exp.matrix, file= tis %&% ".exp.adj." %&% Nk %&% "PEERfactors.3PCs.gender.IDxGENE", quote=F, row.names=F, col.names=F)
write(colnames(adj.exp.matrix), file = tis %&% ".exp.adj." %&% Nk %&% "PEERfactors.3PCs.gender.GENE.list", ncolumns=1)
write(rownames(adj.exp.matrix), file = tis %&% ".exp.adj." %&% Nk %&% "PEERfactors.3PCs.gender.ID.list", ncolumns=1)
