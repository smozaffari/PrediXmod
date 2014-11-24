####script file found in /nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/
####by Heather E. Wheeler 20140602####
args <- commandArgs(trailingOnly=T)
date <- Sys.Date() 
###############################################
### Directories & Variables
 
my.dir <- "/nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/"
 
tissue <- "Nerve - Tibial" ###check GTEx.SampleTissue.annot for available tissues###
#tissue <- "Whole Blood"
tis <- "GTEx-NT"
k <- 10 ### k-fold CV
n <- 10 #number of k-fold CV replicates
gencodefile <- args[1]
gencodeset <- args[2]
 
################################################
### Functions & Libraries
 
"%&%" = function(a,b) paste(a,b,sep="")
 
library(SNPRelate)
library(glmnet)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(x))
lower <- function(x) quantile(x,0.025,na.rm=TRUE)
upper <- function(x) quantile(x,0.975,na.rm=TRUE)

## convenience function to select best lambda over cv bootstraps for model linear by Keston edited by Heather to get predicted values
glmnet.select <- function(response, covariates, nrep.set = 10, nfold.set = 10, alpha.set = 1, ...) {
  require(glmnet)
  best.lam.sim = vector()
  best.cvm.sim = vector()
  pred.matrix = matrix(0,nrow=dim(covariates)[1],ncol=nrep.set)
  for (i in 1:nrep.set) {
    glmnet.fit = cv.glmnet(covariates, response, nfolds = nfold.set, alpha = alpha.set, keep = TRUE)
    new.df = data.frame(glmnet.fit$cvm, glmnet.fit$lambda, glmnet.fit$glmnet.fit$df, 1:length(glmnet.fit$lambda))
    best.lam = new.df[which.min(new.df[,1]),] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
    cvm.best = best.lam[,1]
    nrow.max = best.lam[,4]
    best.lam.sim[i] = nrow.max
    best.cvm.sim[i] = cvm.best
    pred.matrix[,i] = glmnet.fit$fit.preval[,nrow.max]
    }
  cvm.avg = mean(best.cvm.sim) # average cvm
  nrow.max = as.integer(mean(best.lam.sim)) # best lambda over cv bootstraps
  ret <- as.data.frame(glmnet.fit$glmnet.fit$beta[,nrow.max])
  ret[ret == 0.0] <- NA
  ret.vec = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
  names(ret.vec) = rownames(ret)[which(!is.na(ret))]
  min.lambda <- glmnet.fit$glmnet.fit$lambda[nrow.max]
  pred.avg <- rowMeans(pred.matrix)
  output = list(ret.vec, cvm.avg, nrow.max, min.lambda, pred.avg)
#  cat("avg cvm ->", cvm.avg, "lambda iteration ->", nrow.max, "with", length(ret.vec), "effective degrees of freedom \n")
  gc()
  return(output)
}

 
################################################
###input adjusted expression data###
expidfile <- my.dir %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.ID.list"
expgenefile <- my.dir %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.GENE.list"
expfile <- my.dir %&% tis %&% ".exp.adj.15PEERfactors.3PCs.gender.IDxGENE"

expidlist <- scan(expidfile,"character")
expgenelist <- scan(expgenefile,"character")
exp <- scan(expfile)
expdata <- matrix(exp, ncol=length(expgenelist), byrow=T)
rownames(expdata) <- expidlist
colnames(expdata) <- expgenelist

#gencodefile <- my.dir %&% "gencode.v18.genes.patched_contigs.summary.protein"
gencode <- read.table(gencodefile) ##split into 10 files, call each from run_1_CV_GTEx_polyscore_PrediXcan_subset*.sh
rownames(gencode) <- gencode[,5]
tissue.exp <- expdata[,intersect(colnames(expdata),rownames(gencode))] ###pull protein coding gene expression data

###input genotype data###
gtidfile <- my.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.ID.list"
gtsnpfile <- my.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.SNP.list"
gtfile <- my.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.SNPxID"

gtidlist <- scan(gtidfile,"character")
gtidlist <- substr(gtidlist,1,10) ###to match with exp data###
gtsnplist <- scan(gtsnpfile,"character")
gt <- scan(gtfile)
gtdata <- matrix(gt, ncol=length(gtidlist), byrow=T)
rownames(gtdata) <- gtsnplist
colnames(gtdata) <- gtidlist
X <- t(gtdata)
X <- X[expidlist,]

bimfile <- my.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.bim"
bim <- read.table(bimfile)
rownames(bim) <- bim$V2

exp.w.geno <- tissue.exp
explist <- colnames(exp.w.geno)
 
###create results array
resultsarray <- array(0,c(length(explist),7))
dimnames(resultsarray)[[1]] <- explist
dimnames(resultsarray)[[2]] <- c("gene","mean.cvm","mean.lambda.iteration","lambda.min","n.snps","R2","pval")

###run LASSO CV

set.seed(1001)

for(i in 1:length(explist)){
    cat(i,"/",length(explist),"\n")
    gene <- explist[i]
    geneinfo <- gencode[gene,]
    chr <- geneinfo[1]
    c <- substr(chr$V1,4,5)
    start <- geneinfo$V3 - 1e6 ### 1Mb lower bound for cis-eQTLS
    end <- geneinfo$V4 + 1e6 ### 1Mb upper bound for cis-eQTLs
    chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
    cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
    cisgenos <- X[,intersect(colnames(X),cissnps[,2])] ### pull cis-SNP genotypes
    cisgenos[cisgenos >= 3] <- NA ###.gds files have missing=3
#    if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
    if(is.null(dim(cisgenos))){
        bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else if(dim(cisgenos)[2] == 0){
        bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else{
	minorsnps <- subset(colMeans(cisgenos), colMeans(cisgenos,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
        minorsnps <- names(minorsnps)
        cisgenos <- cisgenos[,minorsnps]
        cisgenos <- scale(cisgenos, center=T, scale=T)
        cisgenos[is.na(cisgenos)] <- 0
	if(is.null(dim(cisgenos))){
          bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
        }else if(dim(cisgenos)[2] == 0){
          bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
        }else{

          exppheno <- exp.w.geno[,gene] ### pull expression data for gene
          exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
          exppheno[is.na(exppheno)] <- 0
          rownames(exppheno) <- rownames(exp.w.geno)

          cv <- glmnet.select(exppheno,cisgenos,nrep.set=n,nfold.set=k,alpha.set=1) ###run lasso k-fold CV n times to determine best lambda & betas

          bestbetas <- cv[[1]] ###how many SNPs in best predictor?
        }
    }
    if(length(bestbetas) > 0){
        pred.lasso <- cv[[5]] ###mean k-fold CV predictions from n reps

        ### calculate correlation between predicted and observed expression
        res <- summary(lm(exppheno~pred.lasso))
        genename <- as.character(gencode[gene,6])
        resultsarray[gene,1] <- genename
        resultsarray[gene,2] <- cv[[2]] ###add mean minimum cvm (cross-validated mean-squared error) to results
        resultsarray[gene,3] <- cv[[3]] ###add mean of best lambda iteration to results
        resultsarray[gene,4] <- cv[[4]] ###add best lambda to results
        resultsarray[gene,5] <- length(bestbetas) ###add #snps in prediction to results
        resultsarray[gene,6] <- res$r.squared ###lm R2
        resultsarray[gene,7] <- res$coefficients[2,4] ###lm p-value

        ### output bestbetas for PrediXcan
        bestbetalist <- names(bestbetas)
        bestbetainfo <- bim[bestbetalist,]
        betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
        betafile<-cbind(betatable[,2],betatable[,6],betatable[,7]) ###middle column: [,6] for GEUVADIS, [,5] for GTEx pilot, [,6] for GTEx 2014-06-13 release
        colnames(betafile) <- c("SNP","eff.allele","beta")
#	write.table(betafile, file="/nas40t2/haky/Signatures/data/betas/LASSO/transcriptome-LCL-GD/" %&% genename %&% "-" %&% tis %&% ".txt",quote=F,row.names=F,sep="\t")
        write.table(betafile, file="/nas40t2/hwheeler/PrediXcan_CV/LASSO/hapmap2/transcriptome-" %&% tis %&% "/" %&% genename %&% "-" %&% tis %&% ".txt",quote=F,row.names=F,sep="\t")
    }else{
	genename <- as.character(gencode[gene,6])
        resultsarray[gene,1] <- genename
        resultsarray[gene,2:7] <- c(NA,NA,NA,0,NA,NA)
    }
}


write.table(resultsarray,file=tis %&% ".adj.exp." %&% k %&% "-fold.CV." %&% n %&% "-reps.lasso.PrediXcan.subset." %&% gencodeset %&% "."   %&% date %&% ".txt",quote=F)



