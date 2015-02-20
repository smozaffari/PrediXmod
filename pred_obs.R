####### Rscript scripts/pred_obs.R /group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx-WB.exp.adj.15PEERfactors.3PCs.gender.txt /group/im-lab/nas40t2/smozaffari/Lasso/GTEx_pilot_predicted_alpha1 /group/im-lab/nas40t2/smozaffari/Observed_GTEx/WB/Observed_WB_pilot

args <- commandArgs(TRUE)

original <- args[1]  #original file of observed gene expression, corrected for PCs
#/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx-WB.exp.adj.15PEERfactors.3PCs.gender.txt

# ID ENSG00000237613.2 ENSG00000186092.4 ENSG00000235249.1
#1 GTEX.PWOO       -0.16662285       0.032003635       -0.14115172
#2 GTEX.PX3G       -0.20667492       0.009758687       -0.08671254
#3 GTEX.PLZ5        0.05458398      -0.052856978        0.08435900
#4 GTEX.OXRK       -0.14251255      -0.116945977        0.15908600
#5 GTEX.QVJO        0.02352867      -0.176810495        0.04375056

predicted <- args[2]  #predicted file output from SNP2GReX.pl
#Lasso: /group/im-lab/nas40t2/smozaffari/Lasso/GTEx_pilot_predicted_alpha1
#Elastic Net: /group/im-lab/nas40t2/smozaffari/Elastic_Net/GTEx_pilot_predicted

#Predicted:
#	 gene   GTEX.P4PP   GTEX.PWO3   GTEX.SNOS   GTEX.PW2O
#1 C9orf152 -0.02820108 -0.54730006 -1.11385679 -0.54492097
#2     AACS  0.34167236  0.10777809  0.07381163  0.43808672
#3    FSTL1 -0.43832337 -0.37899113 -0.01571455 -0.76836501
#4    ELMO2 -0.07581912 -0.16853799 -0.17259455  0.01643681
#5    RPS11 -0.13663735 -0.09192275 -0.20130430 -0.21466426

observed <- args[3]
#/group/im-lab/nas40t2/smozaffari/Observed_GTEx/WB/Observed_WB_pilot

origin <- read.table(original, header = T)
#origin <- read.table("/group/im-lab/nas40t2/hwheeler/PrediXcan_CV/GTEx-WB.exp.adj.15PEERfactors.3PCs.gender.txt", header = T)

rownames(origin) <- origin$ID
origin$ID <- NULL
torigin <- t(origin)
write.table(torigin, observed, quote = F)
ensembl <- paste(observed, "ensembl", sep="_")
system(paste("perl /group/im-lab/nas40t2/smozaffari/scripts/ensemblids.pl", observed,  ensembl, sep=" "))

#ensembl <- "/group/im-lab/nas40t2/smozaffari/Observed_GTEx/WB/Observed_WB_pilot_ensembl"
obs <- read.table(ensembl, header = T)
#predicted <- "/group/im-lab/nas40t2/smozaffari/Elastic_Net/GTEx_pilot_predicted"
pred <- read.table(predicted, header = T)

#Observed:
#		Gene   GTEX.PWOO    GTEX.PX3G   GTEX.PLZ5    GTEX.OXRK
#1    FAM138A -0.16662285 -0.206674917  0.05458398 -0.142512549
#2      OR4F5  0.03200363  0.009758687 -0.05285698 -0.116945977
#3     OR4F29 -0.14115172 -0.086712536  0.08435900  0.159086001
#4 AL669831.1 -0.36049792  0.675464297  0.65251709  0.357673903
#5     SAMD11 -0.08088877 -0.014935383 -0.06005709  0.006045255

rownames(pred)<- pred$gene	       #put genes in rownames
pred$gene <- NULL		       #remove column with gene names - now in rows

pred_1 <- pred[,sort(colnames(pred))]          #sort individuals in columns
pred_2 <- pred_1[sort(rownames(pred_1)),]      #sort genes in rows

obs_1 <- obs[, sort(colnames(obs))]	#sort individuals in columns
obs_2 <- obs_1[order(obs_1$Gene),]	#sort genes in rows

obs_3 <- obs_2[which(obs_2$Gene%in%rownames(pred_2)),]		#keep overlap of predicted genes in new observed table
pred_3 <- pred_2[which(rownames(pred_2)%in%obs_2$Gene),]	#keep overlap of observed genes in new predicted table

duplicates <- c()  #could be duplicate genes in observed file
duplicates <- anyDuplicated(obs_3$Gene)  #count duplicate genes in observed file

vars <- c()

#if duplicate genes exist, make two observed tables with each- will be labeled with genename
if (duplicates > 0) {
   files <- c()
   duplicates <- c( duplicates, anyDuplicated(obs_3$Gene, fromLast=T))					#count the other duplicate
   for (i in 1:length(duplicates)) {
       x <- obs_3[-duplicates[i],]								        #remove duplicate gene
       filename <- paste(ensembl, obs_3$Gene[duplicates[i]], duplicates[i], sep="_");			#new filename to output table
       print (paste("There are 2 Ensembl IDs for 1 gene: ", obs_3$Gene[duplicates[i]], sep = "")) ; 	#notify which gene is duplicated
       var_name <- paste(obs_3$Gene[duplicates[i]], duplicates[i], sep="_");	       	                #new name 
       vars <- c(vars, var_name)
       files <- c(files, filename);		                                                        #put new names of table to call on later
       rownames(x) <- x$Gene						                                #assign genenames to row
       x$Gene <- NULL 								                        #remove column with gene names
       x <- x[,-c(which(!colnames(x)%in%colnames(pred_3)))]						#only include individuals in predicted file
       write.table(x, filename, quote = F, row.names = T)			                        #write table to have for later
    }
} else {
  rownames(obs_3) <- obs_3$Gene
  obs_3$Gene <- NULL
  obs_4 <- obs_3[,-c(which(!colnames(obs_3)%in%colnames(pred_3)))]
    new_name <- paste(ensembl, "unique", sep = "_")
    write.table(obs_4, new_name , row.names = T, quote = F)
}

pred_4 <- pred_3[,-c(which(!colnames(pred_3)%in%colnames(obs_3)))]
new_pname <- paste(predicted, "WB_observed_overlap", sep="_")
write.table(pred_4, new_pname, row.names = T, quote = F)

dim(pred_4)

for (j in 1:length(files)) {
    mypvals <- c()
    mycorvec <- c()
    obs <- read.table(files[j], header = T)
    for (i in 1:dim(pred_4)[1]) {
        pred <- pred_4
	genesum <- summary(lm(as.numeric(obs[i,]) ~ as.numeric(pred[i,])))
        ctest <- cor.test(as.numeric(obs[i,]), as.numeric(pred[i,]))
        cor <- c(ctest$estimate[1])
        mycorvec <- c(mycorvec, cor)
        pval <- genesum$coefficients[8];
    	mypvals <- c(mypvals, pval)
    }
    names(mypvals) <- rownames(pred)
    names(mycorvec) <- rownames(pred)

    write.table(mypvals, paste(predicted, vars[j],  "Pvals", sep = "_"), quote = F)
    write.table(mycorvec, paste(predicted, vars[j],  "Corvec", sep = "_"), quote = F)
}