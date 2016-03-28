# The following code is contributed by Sue McClatchy and Gary Churchill
# and is thus not under the purview of the project license, as described in
# README.md

#------------------------------------------------------------------------------

# Hypothesis: Decreased levels of autophagy gene expression in 
# pancreatic islet will result in severe diabetic traits. 
# Model: 
# decreased Atg genes -> increased insulin -> higher blood glucose

# Load QTL library to do genome scans.
library(qtl)
library(ggplot2)
# Load data.
load(file = "data/in/BTBR.clean.data.Rdata")


### Define functions that we will use in the analysis below

####
# some useful plotting functions
# see documentation for "pairs" function
# note that you can adjust the multiplier on cex.cor
#    and the threshold for gray vs black
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="complete.obs")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor*1.0, col=c("gray60", "black")[(abs(r)>0.5)+1])
}
#
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2],0,1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

####
# fit causal models to a triplet with BIC scoring
# X is a transcript used here as first argument to make "apply" easy
# Y is a clincal
# Q is a genotype (factor or numeric)

triple.fit <- function(X,Y,Q){
  #remove any rows with missing values
  indx <- sort(unique(
    c(which(is.na(X)),which(is.na(Y)),which(is.na(Q)))
  ))
  X <- X[-indx]
  Y <- Y[-indx]
  Q <- Q[-indx]
  
  # fit models and compute scores
  b1 <- BIC(lm(X~Q)) + BIC(lm(Y~Q))	#independent X<-Q->Y
  b2 <- BIC(lm(X~Y)) + BIC(lm(Y~Q))	#reactive	 Q->Y->X
  b3 <- BIC(lm(X~Q)) + BIC(lm(Y~X))	#causal		 Q->X->Y
  b4 <- BIC(lm(X~Q)) + BIC(lm(Y~Q+X)) #complex
  scores <- c(b1,b2,b3,b4)
  names(scores) <- c("independent","reactive","causal","complex")
  scores
}

###################################################################
# Find the parkin gene (gene symbol Park2) in the annotation data.
grep("Park2", annot$gene1)
grep("Park2", annot$gene1, value = TRUE)

# Find the ID number for the Park2 gene in the annotation data.
# Use this ID number to pull out expression data from pancreatic islet.
annot[grep("Park2", annot$gene1),]
annot[grep("Park2", annot$gene1), 1]

# Retrieve Park2 expression data from the islet tissue.
park2_islet <- islet.rz[, annot[grep("Park2", annot$gene1), 1]]

# Repeat this procedure for Iapp and Pink1.
grep("Iapp", annot$gene1, value = TRUE)
grep("Pink1", annot$gene1, value = TRUE)
grep("Atg7", annot$gene1, value = TRUE)

# Retrieve expression data for these genes.
iapp_islet <- islet.rz[, annot[grep("Iapp", annot$gene1), 1]]
pink1_islet <- islet.rz[, annot[grep("Pink1", annot$gene1), 1]]
atg7_islet <- islet.rz[, annot[grep("Atg7", annot$gene1), 1]]
park2_islet <- islet.rz[, annot[grep("Park2", annot$gene1), 1]]
atg7_liver <- liver.rz[, annot[grep("Atg7", annot$gene1), 1]]
atg5_liver <- liver.rz[, annot[grep("Atg5", annot$gene1), 1]]
atg5_islet <- islet.rz[, annot[grep("Atg5", annot$gene1), 1]]

# Move the clinical and gene expression phenotypes in to the cross object.
f2g$pheno <- cbind(f2g$pheno[,c("MouseNum","Sex","pgm")],
                   phenotypes.rz[,c("GLU.10wk","INS.10wk","TRIG.10wk","HOMA.10wk")],
                   park2_islet, iapp_islet, pink1_islet, 
                   atg5_liver, atg7_liver,
                   atg7_islet, atg5_islet)

# look at all pairwise scatterplots of clinical and expression traits
pairs(f2g$pheno[,4:length(f2g$pheno)], upper.panel=panel.cor, diag.panel=panel.hist)
# strongest correlation between insulin and a gene expression trait is Atg5 in islet at 0.29
# Also note correlations between Atg5 and Atg7 in islet and liver, and Park2 and Atg5
# in islet

# Pull out sex as a numeric variable so that it can be used as a covariate
# in genome scans.
sex <- as.numeric(f2g$pheno$Sex)

# Calculate genotype probabilities before running genome scans.
f2g <- calc.genoprob(f2g, step = 1)

# Run genome scans with sex as a covariate. This will allow the average
# phenotype values to differ between the two sexes.
scan1 <- scanone(f2g,
                 pheno.col = c("GLU.10wk","INS.10wk","TRIG.10wk","HOMA.10wk",
                               "park2_islet", "iapp_islet", "pink1_islet",
                               "atg5_liver", "atg7_liver",
                               "atg5_islet", "atg7_islet"), 
                 method = "hk", addcovar = sex)

# Identify LOD significance thresholds for gene expression traits.
perm1 <-scanone(f2g, pheno.col = 8, addcovar = sex, 
                method = "hk", n.perm = 100, perm.Xsp = TRUE)
summary(perm1)

# Plot genome scans for each phenotype.
par(mfrow=c(3,1))
for(i in 1:(length(scan1)-2)){
  plot(scan1, lodcolumn = i)
  add.threshold(scan1,
                perms=perm1, alpha=0.05,
                lty="dashed",lwd=2,col="orange")
  add.threshold(scan1,
                perms=perm1, alpha=0.10,
                lty="dashed",lwd=2,col="purple")
}

# Plot genome scans for Park2 and Atg5 together on the same page. 
# Add a title and LOD thresholds.
par(mfrow=c(2,1)) # Plots in 2 rows, one column.
plot(scan1, lodcolumn = 5, main = "Quantitative Trait Loci for 
     Parkin Expression in Pancreatic Islet")
add.threshold(scan1,
              perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1,
              perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 10, main = "Quantitative Trait Loci for 
     Atg5 Expression in Pancreatic Islet")
add.threshold(scan1,
              perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1,
              perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

# Plot genome scans for Atg5, Atg7, and insulin. Add a title and LOD thresholds.
par(mfrow=c(3,1)) # Plots in 3 rows, 1 column.
plot(scan1, lodcolumn = 2, main = "Quantitative Trait Loci for 
     Insulin at 10 Weeks")
add.threshold(scan1,
              perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1,
              perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 10, main = "Quantitative Trait Loci for 
     Atg5 Expression in Pancreatic Islet")
add.threshold(scan1,
              perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1,
              perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

plot(scan1, lodcolumn = 11, main = "Quantitative Trait Loci for 
     Atg7 Expression in Pancreatic Islet")
add.threshold(scan1,
              perms=perm1, alpha=0.05,lty="dashed",lwd=2,col="orange")
add.threshold(scan1,
              perms=perm1, alpha=0.10,lty="dashed",lwd=2,col="purple")

# View genome scan summaries in different formats.
summary(scan1)
summary(scan1, format = "tabByChr")
summary(scan1, format = "tabByChr", perms=perm1, alpha = 0.05)

# Insulin, Atg5 and Atg7 expression share the same peak on chr2.
# Note that their confidence intervals are nested one inside another.
# What on chr2 drives Atg5, Atg7, or insulin?
# We can download a list of all genes within the largest confidence
# interval (between ci.low and ci.high) from 50.0 cM to 83.6 cM.
# Locate the markers nearest 50.0 cM and 83.6 cM. 
find.marker(cross = f2g, chr = 2, pos = 50.0)
find.marker(cross = f2g, chr = 2, pos = 83.6)

# SNPs rs13476645 and rs8260429 on chromosome 2 define the largest
# confidence interval for the Atg5, Atg7, and insulin peaks.
# Go to ensembl.org to find base pair positions for these SNPs.
# rs13476645 is at 95078782 bp and rs8260429 is at 163751760 bp.
# To download gene symbols, click on BioMart, choose database
# Ensembl Genes, choose dataset Mus musculus genes. 
# Click Filters. Open up Region, select chromosome 2, then enter
# the two base pair positions for gene start and end.
# Click Attributes. Open up Gene, uncheck Ensembl Gene ID and Ensembl
# Transcript ID. Open up External. Check MGI symbol.
# Click Results button at upper left. Export a file as TSV.
# The file should be named mart_export.txt. Change the name
# to something descriptive, like chr2_genes.txt.   
# Read the chromosome 2 gene list into R.

chr2_genes <- scan(file = "chr2_genes.txt", 
                   what = "character", 
                   skip = 1)

##########################################################
### I. conditional genome scans
# looking for a gene expression trait that "blocks" the chr 2 QTL peak
# Atg5 and Atg7 in islet are likely suspects

# scan insulin conditional on Atg gene expression traits
scan_atg5 <- scanone(f2g, pheno.col="INS.10wk", 
                     addcovar=f2g$pheno$atg5_islet, method="hk")
scan_atg7 <- scanone(f2g, pheno.col="INS.10wk", 
                     addcovar=f2g$pheno$atg7_islet, method="hk")


# plot the conditional scans
par(mfrow=c(3,1))
plot(scan1, lodcolumn = 2,
     main = "Quantitative Trait Loci for Insulin at 10 Weeks")
add.threshold(scan1,
              perms=perm1, alpha=0.05,
              lty="dashed", lwd=2, col="orange")
add.threshold(scan1,
              perms=perm1, alpha=0.10,
              lty="dashed", lwd=2, col="purple")

plot(scan_atg5,
     main = "Insulin Scan with Islet Atg5 Expression as Covariate",
     ylab = "INS.10wk")
add.threshold(scan1,
              perms=perm1, alpha=0.05,
              lty="dashed", lwd=2, col="orange")
add.threshold(scan1,
              perms=perm1, alpha=0.10,
              lty="dashed", lwd=2, col="purple")

plot(scan_atg7,
     main = "Insulin Scan with Islet Atg7 Expression as Covariate",
     ylab = "INS.10wk")
add.threshold(scan1,
              perms=perm1, alpha=0.05,
              lty="dashed", lwd=2, col="orange")
add.threshold(scan1,
              perms=perm1, alpha=0.10,
              lty="dashed", lwd=2, col="purple")

# Loop through all genes in the chromosome 2 interval, adding each
# in as covariate, running the scan for insulin, Atg5, and Atg7.

# Add chr2 gene expression traits into cross object.
f2g$pheno <- cbind(f2g$pheno, 
                   islet.rz[,match(chr2_genes,
                                   annot$gene1,
                                   nomatch = 0)])
names(f2g$pheno)

# Column names for gene expression traits are the microarray probe IDs.
# Replace probe IDs with gene symbols.
names(f2g$pheno)[15:ncol(f2g$pheno)] <- annot$gene1[match(names(f2g$pheno)[15:ncol(f2g$pheno)],
                                                          annot$a_gene_id)]
names(f2g$pheno)

# scan insulin, Atg5, and Atg7 conditional on chr2 gene expression traits.
# Scan only for chromosome 2.
scan_cond <- scanone(f2g, 
                     pheno.col=c("INS.10wk", "atg5_islet", "atg7_islet"),
                     addcovar=f2g$pheno$Crnkl1, method="hk")

# note how we used "cbind" to concatenate the results
#  "cbind" calls the specialized function "cbind.scanone"
# Scan only on chromosome 2.
for(i in 16:ncol(f2g$pheno)){
  scan_cond <- cbind(scan_cond,
                  scanone(f2g,
                          pheno.col=c("INS.10wk", "atg5_islet", "atg7_islet"),
                          addcovar=f2g$pheno[,i], method="hk") )
}

# Re-name the conditional scan columns with the gene symbol. Every 3 columns
# will reference a gene scanned conditionally on insulin, Atg5, and Atg7 in
# that order.
head(names(scan_cond))
dim(scan_cond)
for (i in 15:(ncol(f2g$pheno))) {
  names(scan_cond)[(3*(i-14)):(3*(i-14)+2)] <- names(f2g$pheno)[i]
}

par(mfrow=c(6,1), mar=c(3,4,1,4) + 0.1)
for( i in 1:ncol(scan_cond)){
  plot(scan_cond, lodcolumn=i, main=names(scan_cond)[i+2])
  add.threshold(scan1, perms=perm1, alpha=0.05, 
                lty="dashed", lwd=2, col="orange")
  add.threshold(scan1, perms=perm1, alpha=0.10, 
                lty="dashed", lwd=2, col="purple")
}

# Several genes drop the chromosome 2 peaks for all 3 phenotypes
# (insulin, Atg5, Atg7) down below significance. Find these genes in 
# the scan object. Subtract 2 from each index number to account 
# for the first two columns in the scan object, which are 
# chromosome number and cM position.

# Which genes drop the chromosome 2 peak below the 10%
# significance threshold of 3.67?
for (i in 3:length(summary(subset(scan_cond, chr = 2)))) {
  if (summary(subset(scan_cond, chr = 2))[[i]]< 3.67)
    print(summary(subset(scan_cond, chr = 2))[i])
  }

# Look for gene symbols repeated 3 times, indicating that the gene
# influences insulin, Atg5, and Atg7 expression. Note: this is not
# a fail-safe method. Check the plots to verify. Look for flattened
# peaks on chromosome 2 for all three scans - insulin, Atg5, Atg7.
summary(subset(scan_cond, chr = 2))[,which(names(scan_cond) 
                                           %in% c("Pdrg1", "Gatm", 
                                                  "Nphp1", "Cds2", 
                                                  "Ino80", "Dtd1", 
                                                  "Aqr", "Vps39", 
                                                  "Adal", "Cbfa2t2"))]

# Plot each and compare to original scans for insulin, Atg5, and Atg7.
plot(scan1, lodcolumn = 2)
add.threshold(scan1, perms=perm1, alpha=0.05,
              lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,
              lty="dashed", lwd=2, col="purple")

plot(scan1, lodcolumn = 10)
add.threshold(scan1, perms=perm1, alpha=0.05,
              lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,
              lty="dashed", lwd=2, col="purple")

plot(scan1, lodcolumn = 11)
add.threshold(scan1, perms=perm1, alpha=0.05,
              lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,
              lty="dashed", lwd=2, col="purple")

# Pdrg1 conditional scans
plot(scan_cond, lodcolumn = 25) # for insulin
add.threshold(scan1, perms=perm1, alpha=0.05,
              lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,
              lty="dashed", lwd=2, col="purple")

plot(scan_cond, lodcolumn = 26) # for Atg5
add.threshold(scan1, perms=perm1, alpha=0.05,
              lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,
              lty="dashed", lwd=2, col="purple")

plot(scan_cond, lodcolumn = 27) # for Atg7
add.threshold(scan1, perms=perm1, alpha=0.05,
              lty="dashed", lwd=2, col="orange")
add.threshold(scan1, perms=perm1, alpha=0.10,
              lty="dashed", lwd=2, col="purple")


# Plot all other conditional genome scans for Gatm, Nphp1, Cds2, Ino80, 
# Dtd1, Aqr, Vps39, Adal, Cbfa2t2.
# Use the index numbers produced by which(), less
# 2 to account for the first two columns: chr and pos.
which(names(scan_cond) 
      %in% c("Pdrg1", "Gatm", "Nphp1", "Cds2", "Ino80", 
             "Dtd1", "Aqr", "Vps39", "Adal", "Cbfa2t2"))-2

##########################################################
####  II. Mediation analysis
#  use linear models to compute mediation on the peak chr 8 marker
#  this is just another way to get at the same question

###
# get the genotypes at peak marker for insulin
#
# fill in the missing genotypes
f2g <- fill.geno(f2g, method="argmax")
#
#create a three level factor Q2 and add to pheno data
f2g$pheno <- transform(f2g$pheno,
                       Q2 = as.factor(
                         f2g$geno[[2]]$data[,find.marker(
                           f2g, 2, 75.2)]))
levels(f2g$pheno$Q2) <- c("B","H","R")

# plot insulin against Pdrg1 expression
pdrg1_islet <- islet.rz[, annot[grep("Pdrg1", annot$gene1), 1]]
f2g$pheno <- transform(f2g$pheno, pdrg1_islet)
qplot(pdrg1_islet, INS.10wk, color=Q2, shape=Sex, data=f2g$pheno) +
  geom_smooth(aes(group=Q2), method="lm", se=FALSE)
# There's reasonably good correlation between genotype classes.

#####
# check 4 conditions for Pdrg1 gene expression 
# as a mediator of Q2 effect on insulin

####  i) Insulin is linked to Q2
anova(lm(INS.10wk ~ Sex + Q2, data = f2g$pheno))
# significant

####  ii) Pdrg1 gene expression is linked to Q2
anova(lm(pdrg1_islet ~ Sex + Q2, data = f2g$pheno))
# significant

####  iii) Insulin not linked after accounting for Q2
anova(lm(INS.10wk ~ Sex + pdrg1_islet + Q2, data = f2g$pheno))
# not significant

####  iv) Pdrg1 gene expression is still linked after 
# accounting for insulin
anova(lm(pdrg1_islet ~ Sex + INS.10wk + Q2, data = f2g$pheno))
# significant

# all 4 conditions for a mediator are satisfied
##########

##########################################################
#  III. establish mediator using model selection with BIC scoring

# note that missing values will mess up BIC analysis
apply(is.na(f2g$pheno), 2, sum)
#
# easiest to use the triple.fit function that removes missing data
# and fits the three models

# compute BIC scores for the causal triplet models
#   using pdrg1_islet
with(f2g$pheno,
     triple.fit(pdrg1_islet, INS.10wk, Q2))
#  causal model has lowest score
#  suggests that pdrg1_islet is a strong mediator of insulin

# compute BIC scores for the causal triplet models
#   using atg5_islet and atg7_islet
with(f2g$pheno,
     triple.fit(pdrg1_islet, atg5_islet, Q2))
#  causal model has lowest score
#  suggests that pdrg1_islet is a strong mediator of Atg5 expression in islet

with(f2g$pheno,
     triple.fit(pdrg1_islet, atg7_islet, Q2))
#  causal model has lowest score but differs by only 5 from complex model
#  therefore inconclusive


