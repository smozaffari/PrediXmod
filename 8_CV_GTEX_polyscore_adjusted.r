####script file found in /nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/
####by Heather E. Wheeler 20140602####
args <- commandArgs(trailingOnly=T)
date <- Sys.Date() 
"%&%" = function(a,b) paste(a,b,sep="")
###############################################
### Directories & Variables
cri = T #on cri cluster?
if(cri) precri = "/group/dolan-lab/hwheeler/" else precri = '/'

tis <- "GTEx-NT"
k <- 10 ### k-fold CV
n <- 10 #number of k-fold CV replicates
 
my.dir <- precri %&% "nas40t2/hwheeler/PrediXcan_CV/GTEx_2014-06013_release/"
poly.dir <- precri %&% "nas40t2/hwheeler/PrediXcan_CV/Polyscore/hapmap2/transcriptome-" %&% tis %&% "/"
gencodefile <- my.dir %&% "gencode.v18.genes.patched_contigs.summary.protein." %&% args[1]
gencodeset <- args[1]

################################################
### Functions & Libraries
 
library(glmnet)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(x))
lower <- function(x) quantile(x,0.025,na.rm=TRUE)
upper <- function(x) quantile(x,0.975,na.rm=TRUE)

## fit betas by Vasa edited by Heather to get pvals
fit.betas <- function(G, Y) {
require(RcppGSL)
betas <- rep(NA, ncol(G))
stderrs <- rep(NA, ncol(G))
dfs <- rep(NA, ncol(G))
for (i in 1:ncol(G)) {
res <- fastLmPure(y=Y, X=matrix(G[,i]))
betas[i] <- res$coefficients
stderrs[i] <- res$stderr
dfs[i] <- res$df
}
t <- betas/stderrs
pvals <- 2*pt(-abs(t),df=dfs)
output <- cbind(betas,pvals)
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
gtfile <- my.dir %&% "GTEx.IDxSNP.rds"

X <- readRDS(gtfile)
X <- X[expidlist,]

samplelist<-expidlist

bimfile <- my.dir %&% "GTEx_Analysis_2014-06-13.hapmapSnpsCEU.bim"
bim <- read.table(bimfile)
rownames(bim) <- bim$V2

exp.w.geno <- tissue.exp
explist <- colnames(exp.w.geno)
 
###create results array
resultsarray <- array(NA,c(n,6,length(explist)))
dimnames(resultsarray)[[3]] <- explist
finalresults <- matrix(NA,nrow=length(explist),ncol=7)
rownames(finalresults) <- explist
colnames(finalresults) <- c("gene","R2.mean","R2.se","R2.lci","R2.uci","mean.n.snps","P.mean")

###run polyscore CV

set.seed(1001)

for(j in 1:length(explist)){
    cat(j,"/",length(explist),"\n")
    gene <- explist[j]
    genename <- as.character(gencode[gene,6])
    geneinfo <- gencode[gene,]
    chr <- geneinfo[1]
    c <- substr(chr$V1,4,5)
    start <- geneinfo$V3 - 1e6 ### 1Mb lower bound for cis-eQTLS
    end <- geneinfo$V4 + 1e6 ### 1Mb upper bound for cis-eQTLs
    chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
    cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
    cisgenos <- X[,intersect(colnames(X),cissnps[,2])] ### pull cis-SNP genotypes
    cisgenos[cisgenos >= 3] <- NA ###.gds files have missing=3
    if(is.null(dim(cisgenos))){###effectively skips genes with <2 cis-SNPs
                finalresults[gene,1] <- genename
                finalresults[gene,2:6] <- rep(NA,5)
        }else{
              	minorsnps <- subset(colMeans(cisgenos), colMeans(cisgenos,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
                minorsnps <- names(minorsnps)
                cisgenos <- cisgenos[,minorsnps]
                cisgenos <- scale(cisgenos, center=T, scale=T) ##need to scale for fastLmPure to work properly
                cisgenos[is.na(cisgenos)] <- 0
                if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs after removing those with no minor alleles
                        finalresults[gene,1] <- genename
                        finalresults[gene,2:6] <- rep(NA,5)
                }else{
                      	exppheno <- exp.w.geno[,gene] ### pull expression data for gene
                        exppheno <- scale(exppheno, center=T, scale=T)
                        exppheno[is.na(exppheno)] <- 0
                        rownames(exppheno) <- rownames(exp.w.geno)
                        ### run k-fold CV n times
                        for(i in 1:n){
                                ###randomize samples into CV groups
                                samplelength <- length(exppheno)
                                g <- 1:k ##k-fold CV
                                groupid <- sample(g,samplelength,replace=T)
                                newiddata <- data.frame(groupid,samplelist)
                                ### create empty matrix to fill in with polyscore vals
                                cvset <- data.frame(samplelist,matrix(NA,samplelength,1),exppheno)
                                rownames(cvset) <- samplelist
                                names(cvset) <- c("ID","Pred","Exp")
                                for(idsubset in 1:k){
                                        trainset <- data.frame(newiddata$samplelist[newiddata$groupid != idsubset])
                                        testset <- data.frame(newiddata$samplelist[newiddata$groupid == idsubset])
                                        ### pull trainset genos and phenos
                                        traingenos <- cisgenos[intersect(trainset[,1],rownames(cisgenos)),]
                                        trainexp <- exppheno[match(trainset[,1],rownames(exppheno))]  ###if scaled, use rownames
                                        ### calculate betas and pvals from training set
                                        betares <- fit.betas(traingenos,trainexp) ###~30sec/gene
                                        rownames(betares) <- colnames(traingenos)
                                        topres <- subset(betares,betares[,2] < 0.05) ###pull snps with P<0.05
                                        betas <- topres[,1]
                                        if(length(betas) <= 1){
                                                names(betas) <- rownames(topres)
                                        }
                                        ### polyscore
                                        ### multiply betas by test genotypes and take sum
                                        testgenos <- cisgenos[intersect(testset[,1],rownames(cisgenos)),intersect(names(betas),colnames(cisgenos))]
                                        if(is.array(testgenos)=='FALSE'){ ### calcs polyscore if only one individual in test set, happens when expression sample size is low (<100)
                                                if(dim(testset)[1] < 2){ ### calcs polyscore if only one individual in test set, happens when expression sample size is low (<100)
                                                        pred.polyscore <- sum(testgenos * betas)
                                                        names(pred.polyscore)<-testset[,1]
                                                }else{ ### calcs polyscore if only one beta with p<0.05
                                                        pred.polyscore <- testgenos * betas
                                                }
                                                cvset$Pred[match(names(pred.polyscore),rownames(cvset))] <- pred.polyscore
                                        }else{
                                              	testsweep <- sweep(testgenos, 2, betas, "*")
                                                pred.polyscore <- as.vector(rowSums(testsweep))
                                                names(pred.polyscore) <- rownames(testsweep)
                                                ### add predictions to cvset
                                                cvset$Pred[match(names(pred.polyscore),rownames(cvset))] <- pred.polyscore
                                        }
                                }

                                ### calculate correlation between predicted and observed expression
                                res<-summary(lm(cvset$Pred~cvset$Exp))
                                resultsarray[i,1:4,gene] <- coef(res)[2,]
                                resultsarray[i,5,gene] <- res$r.squared
                                resultsarray[i,6,gene] <- length(betas)


                        }
                        finalresults[gene,1] <- genename
                        finalresults[gene,2] <- mean(resultsarray[,5,gene],na.rm=TRUE)
                        finalresults[gene,3] <- stderr(resultsarray[,5,gene])
                      	finalresults[gene,4] <- lower(resultsarray[,5,gene])
                        finalresults[gene,5] <- upper(resultsarray[,5,gene])
                        finalresults[gene,6] <- as.integer(mean(resultsarray[,6,gene],na.rm=TRUE))
                        finalresults[gene,7] <- mean(resultsarray[,4,gene],na.rm=TRUE)

                        ### output bestbetas for PrediXcan
                        allbetares <- fit.betas(cisgenos,exppheno) ###~30sec/gene
                        rownames(allbetares) <- colnames(cisgenos)
                        alltopres <- subset(allbetares,allbetares[,2] < 0.05) ###pull snps with P<0.05
                        allbetas <- alltopres[,1]
                        if(length(allbetas) <= 1){
                                names(allbetas) <- rownames(alltopres)
                        }
                        allbetalist <- names(allbetas)
                        allbetainfo <- bim[allbetalist,]
                        allbetatable <- as.matrix(cbind(allbetainfo,alltopres))
                        betafile <- cbind(allbetatable[,2],allbetatable[,6],allbetatable[,7],allbetatable[,8])
                        colnames(betafile) <- c("SNP","eff.allele","beta","p.value")
                        write.table(betafile, file=poly.dir %&% genename %&% "-" %&% tis %&% ".txt",quote=F,row.names=F,sep="\t")
                }
        }
}

date <- Sys.Date()



write.table(finalresults,file=tis %&% ".adj.exp." %&% k %&% "-fold.CV." %&% n %&% "-reps.polyscore.p0.05.PrediXcan.subset." %&% gencodeset %&% "."   %&% date %&% ".txt",quote=F)



