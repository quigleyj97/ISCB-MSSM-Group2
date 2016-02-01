###########################################################
#  Analysis of Haptoglobin in Attie BTBR eQTL data
#  Dec 5, 2015
#  copyright Gary Churchill
#
# this script carries out a genomescan of haptoglobin in plasma
# as well as haptoglobin expression in six tissues
# we observe a string QTL peak for Haptoglobin on Chr 8
# that is shared by several of the gene expression traits
# we then use three different methods to establish that
# liver gene expression is a causal mediator of the QTL 8 effect
#
# Note: Haptoglobin is a protein that binds free hemoglobin in plasma
#   preventing it from doing damage outside of red blood cells
###########################################################

###########################################################
### load libraries
library(qtl)
# library(dplyr)  #need this?
library(ggplot2)

###########################################################
### Define functions that we will use in the analsysis below

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
# a function to plot genome scans with additive and interactive covariates
# and show the difference in lod scores
# could be rewritten to be general - this version depends on environment
#   variable that will be created in the script below
three.plot <- function(lc)
{
  #plot the genome scans
  par(mfrow=c(3,1))
  #
  plot(scan1a, lodcolumn=lc,
       main=paste(names(scan1a)[lc+2],"w/ sex as addcov") )
  add.threshold(scan1a, perms=f2g.perm1a, alpha=0.05, lty="dashed", lwd=2, col="red")
  add.threshold(scan1a, perms=f2g.perm1a, alpha=0.63, lty="dashed", lwd=2, col="green")
  #
  plot(scan1b, lodcolumn=lc,
       main=paste(names(scan1a)[lc+2],"w/ sex as intcov"))
  add.threshold(scan1b, perms=f2g.perm1b, alpha=0.05, lty="dashed",lwd=2, col="red")
  add.threshold(scan1b, perms=f2g.perm1b, alpha=0.63, lty="dashed",lwd=2, col="green")
  #
  plot(scan1b - scan1a, lodcolumn=lc,
       main=paste(names(scan1a)[lc+2],"sex interaction test") )
  abline(2,0,col="red",lty=2)
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

###########################################################
### import the cross data and get set up

#read cross data into R environment
# may vary depending on your local file configuration
load("BTBR.clean.data.Rdata")

#clean up the f2g cross phenotypes
f2g$pheno <- f2g$pheno[,c("MouseNum","Sex","pgm")]

# pull Haptoglobin phenotype and gene expression data into the 
# f2 cross data
# I did some homework to discover the gene symbol for haptoglobin is Hp
# the gene is on chromosome 8 at 109.6Mb - distal
# note rankZ
f2g$pheno <- transform(f2g$pheno, Haptoglobin=rz.transform(phenotypes$Haptoglobin),
                      Hp.adipose = adipose.rz[,which(annot$gene_symbol=="Hp")],
                      Hp.gastroc = gastroc.rz[,which(annot$gene_symbol=="Hp")],
                      Hp.hypo = hypo.rz[,which(annot$gene_symbol=="Hp")],
                      Hp.islet = islet.rz[,which(annot$gene_symbol=="Hp")],
                      Hp.kidney = kidney.rz[,which(annot$gene_symbol=="Hp")],
                      Hp.liver = liver.rz[,which(annot$gene_symbol=="Hp")] )
dim(f2g$pheno)
names(f2g$pheno)

# look at all pairwise scatterplots of Haptoglobin traits
pairs(f2g$pheno[,4:10], upper.panel=panel.cor, diag.panel=panel.hist)
# strongest correlation with circulating haptoglobin is liver at 0.6
# muscle and adipose also have strong positive correlations
# followed by hypothalamus and kidney

##########################################################
### GENOME SCANS
#########

###set up for genome scans
f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

# need numeric sex for covariate
sex <- as.numeric(f2g$pheno$Sex)-1

f2g.perm1a <- scanone(f2g, pheno.col = 5, addcovar = sex, 
                      method = "hk", n.perm = 100, perm.Xsp = TRUE)
f2g.perm1b <- scanone(f2g, pheno.col = 5, intcovar = sex, 
                      method = "hk", n.perm = 100, perm.Xsp = TRUE)

# Save the permutations so that you can load them later instead
# of re-computing. Un-comment lines 142 and 143 above to load permutations.
save(list="f2g.perm1a", file="f2g_perm1a.Rdata")
save(list="f2g.perm1b", file="f2g_perm1b.Rdata")

#precomputed permutations on .rz variables scan with sex as addcov and intcov
# load(file="f2g_perm1a.Rdata")
# load(file="f2g_perm1b.Rdata")

# run scans
scan1a <- scanone(f2g, pheno.col=4:10, addcovar=sex, method="hk")
scan1b <- scanone(f2g, pheno.col=4:10, intcovar=sex, method="hk")

####
# plot the genome scans

# generate all seven scans into a pdf
pdf("haptoglobin_scans.pdf")
# draw the plot for each trait
for( i in 1:7){
  three.plot(i)
}
dev.off()  #dont forget to close the device

####
# summarize QTL peaks
summary(scan1a, perms=f2g.perm1a, alpha=0.10, format="tabByChr", ci.function="lodint")

##########################################################
### I. conditional genome scans
# which Tissue is driving the serum Hapotoglobin levels?
# looking for a gene expression trait that "blocks" the chr 8 QTL peak

# scan Haptoglobin conditional on gene expression traits
scan1c <- scanone(f2g, pheno.col=4, addcovar=f2g$pheno$Hp.adipose, method="hk")

# note how we used "cbind" to concatenate the results
#  "cbind" calls the specialized function "cbind.scanone"
for(i in 6:10){
  scan1c <- cbind(scan1c,
                 scanone(f2g, pheno.col=4, addcovar=f2g$pheno[,i], method="hk") )
}

names(scan1c)[3:8] <- names(f2g$pheno)[5:10]
head(scan1c)

# plot the conditional scans
pdf("haptoglobin_cond_scans.pdf")
par(mfrow=c(6,1),mar=c(3,4,1,4)+0.1)
for( i in 1:6){
  plot(scan1c, lodcolumn=i, main=names(scan1c)[i+2], xlab="", ylab="")
  add.threshold(scan1a, perms=f2g.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
}
dev.off()
# liver appears to be the primary mediator
# there a residual genetic effect on chr 6 but no obvious tissue driving it

##########################################################
####  II. Mediation analysis
#  use linear models to compute mediation on the peak chr 8 marker
#  this is just another way to get at the same question

###
# get the genotypes at peak marker
#
# fill in the missing genotypes
f2g <- fill.geno(f2g, method="argmax")
#
#create a three level factor Q8 and add to pheno data
f2g$pheno <- transform(f2g$pheno,
      Q8 = as.factor(f2g$geno[[8]]$data[,find.marker(f2g,8,57.9)]) )
levels(f2g$pheno$Q8) <- c("B","H","R")

# plot Haptoglobin against liver expression
# note correlated within genotype group
qplot(Hp.liver, Haptoglobin, color=Q8, shape=Sex, data=f2g$pheno) +
  geom_smooth(aes(group=Q8),method="lm",se=FALSE)
  
#####
# check 4 conditions for liver gene expression as a mediator of Q8 effect on Haptoglobin

####  i) Haptoglobin is linked to Q8
anova(lm(Haptoglobin~Sex+Q8,data=f2g$pheno))
# significant

####  ii) Liver gene expression is linked to Q8
anova(lm(Hp.liver~Sex+Q8,data=f2g$pheno))
# significant

####  iii) Haptoglobin not linked after accounting to Q8
anova(lm(Haptoglobin~Sex+Hp.liver+Q8,data=f2g$pheno))
# not significant

####  iv) Liver gene expression is still linked after accouting for Haptoglobin
anova(lm(Hp.liver~Sex+Haptoglobin+Q8,data=f2g$pheno))
# significant

# all 4 conditions for a mediator are satistified
##########

# repeat mediation analysis for Adipose
qplot(Hp.adipose, Haptoglobin, color=Q8, shape=Sex, data=f2g$pheno) +
  geom_smooth(aes(group=Q8),method="lm",se=FALSE)
#weak correlation within genotype classes
#
anova(lm(Haptoglobin~Sex+Q8,data=f2g$pheno))
anova(lm(Hp.adipose~Sex+Q8,data=f2g$pheno))
anova(lm(Haptoglobin~Sex+Hp.adipose+Q8,data=f2g$pheno))
anova(lm(Hp.adipose~Sex+Haptoglobin+Q8,data=f2g$pheno))
# note that mediation condition iii) fails
# we could go ahead and check other tissues, but...

##########################################################
#  III. establish mediator using model selection with BIC scoring

# note that missing values will mess up BIC analysis
apply(is.na(f2g$pheno),2,sum)
#
# easiest to use the triple.fit function that removes missing data
# and fits the three models

# compute BIC scores for the causal triplet models
#   using Hp.liver
with(f2g$pheno,
     triple.fit(Hp.liver, Haptoglobin, Q8))
# independent    reactive      causal     complex
#    2660.318    2530.023    2466.512    2478.552
#  causal model has lowest score
#  it is better than the next best model (complex) by more than 10 units
#  and so it is significantly best
#  suggests that Hp.liver is a strong mediator


# compute BIC scores for the causal triplet models
#   using Hp.adipose
with(f2g$pheno,
     triple.fit(Hp.adipose, Haptoglobin, Q8))
# independent    reactive      causal     complex
#    2755.024    2734.562    2737.672    2733.997
#  complex model has lowest score but
#  all three models (other than independent) are within 5 units
#  and therefore indistinguishable.








