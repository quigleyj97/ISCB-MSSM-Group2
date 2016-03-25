# The following code is contributed by Gary Churchill and
# is thus not under the purview of the project license
# or style guide, as described in README.md

#----------------------------------------------------------

###########################################################
#  Attie BTBR eQTL data
#  Leptin and Islet GEX first pass
#  Uses the "clean" data. 
#  Feb 5, 2013 - GAC
###########################################################
# 
#set working directory
#setwd("/Users/garyc/Desktop/BTBR_ISCB_Data")

#load qtl library
library(qtl)
library(ggplot2)

# some useful plotting functions
# see documentation for "pairs" function 
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use="complete.obs")
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor*1.0, col=c("gray60", "black")[(abs(r)>0.35)+1])
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

###########################################################
### load the data 
load("BTBR.clean.data.Rdata")
ls()
# smc - Chelse and Larry: I recommend that you reduce the clean
# data set by eliminating some of the objects you won't use in
# this script, including kidney, liver, gastrocnemius, and 
# hypothalamus data. Use the following command to do so.

rm(list=c("kidney.rz", "liver.rz", "gastroc.rz", "hypo.rz"))
ls()

# Your machine should be much happier now.

###########################################################
# select some phenotypes of interest
my.phenos <- phenotypes.rz[,c("GLU.10wk","INS.10wk","TRIG.10wk","Leptin","IL.1beta","NEFA","FGF.basic")]

# smc - use long names for now and return to clean script to figure
# out why short names didn't take
# my.phenos <- phenotypes.rz[,c("GLUCOSE (mg/dl) 10 wk","INSULIN (ng/ml) 10 wk","TRIGLYCERIDE (mg/dl) 10 wk", "Leptin","IL.1beta","NEFA","FGF.basic")]

names(my.phenos)

# scatterplots
# smc - use x11() instead of x11 if you're on a PC instead of a Macintosh

x11()

pairs(my.phenos, upper.panel=panel.cor,diag.panel=panel.hist)

###########################################################
### scan phenotypes for QTL
#
# move phenotypes into the cross structure
f2g$pheno <- cbind(f2g$pheno[,c(2,6,7)],my.phenos)

# setup for scanning
f2g <- calc.genoprob(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)
f2g <- sim.geno(f2g, step=1, stepwidth="fixed", map.function="c-f", err=0.002)

#convenient to keep "sex" handy as a numeric variable
sex <- as.numeric(f2g$pheno$Sex)

#scan with sex as an additive covariate
my.scan1a <- scanone(f2g, pheno.col=c(4:10), addcovar=sex, method="hk")
#
#run permutations
my.perm1a <-scanone(f2g,pheno.col=4,addcovar=sex,method="hk",n.perm=100,perm.Xsp=TRUE)

# plot scans
#x11()
for(i in 1:7){
# x11()
plot(my.scan1a,lodcolumn=i)
add.threshold(my.scan1a, lodcolumn=1,
	perms=my.perm1a, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(my.scan1a, lodcolumn=1,
	perms=my.perm1a, alpha=0.20,lty="dashed",lwd=2,col="green")
}
#report results
summary(my.perm1a)
# 
summary(my.scan1a, perms=my.perm1a, alpha=0.10, format="allpheno")

###########################################################
###look for correlates of Leptin in pancreas gex data
#
#create a correlation "scan" function
# x is a matrix or dataframe with lots varaibles 
# y is a vector that we will assign outisde of the function
# "apply" this function to compute correlation of y with each element of x
mycorr <- function(x){
	cor(x,y,use="complete.obs")
	}
#
x11()
qplot(Sex,Leptin,data=f2g$pheno)

# compute correlations of Srebf1.liver with liver expression
y <- phenotypes$Leptin
islet.cor <- apply(islet.rz,2,"mycorr")
#
x11()
hist(islet.cor, breaks=100)

# what are the most highly correlated genes?
sum(abs(islet.cor)>0.25)
#
indx.islet <- which(abs(islet.cor)>0.25)
annot[indx.islet, c("a_gene_id","gene_symbol")]

###########################################################
# genome scan for islet gene expression
#TBD

f2g$pheno <- cbind(f2g$pheno, islet.rz[,indx.islet])
names(f2g$pheno)
# columns 11 through 40 are GEX
#
names(f2g$pheno)[11:40] <- annot[indx.islet, c("gene_symbol")]
names(f2g$pheno)

# correlation heatmap
x11()
heatmap(cor(f2g$pheno[,-c(1:3)]))
#
# scatterplots
x11()
pairs(f2g$pheno[,11:26],
	upper.panel=panel.cor,diag.panel=panel.hist)
	
#scan with sex as an additive covariate
my.scan1b <- scanone(f2g, pheno.col=c(11:26), addcovar=sex, method="hk")
#
#run permutations
my.perm1b <-scanone(f2g,pheno.col=15,addcovar=sex,method="hk",n.perm=100,perm.Xsp=TRUE)

# plot scans
#x11()
for(i in 1:7){
x11()
plot(my.scan1b,lodcolumn=i)
add.threshold(my.scan1b, lodcolumn=1,
	perms=my.perm1b, alpha=0.05,lty="dashed",lwd=2,col="red")
add.threshold(my.scan1b, lodcolumn=1,
	perms=my.perm1b, alpha=0.20,lty="dashed",lwd=2,col="green")
}
#report results
summary(my.perm1b)
# 
summary(my.scan1b, perms=my.perm1b, alpha=0.10, format="allpheno")

###########################################################
###look for correlates of Leptin in adipose gex data
#
# compute correlations of Srebf1.liver with liver expression
y <- phenotypes$Leptin
adipose.cor <- apply(adipose.rz,2,"mycorr")
#
x11()
hist(adipose.cor, breaks=100)
#
# what are the most highly correlated genes?
sum(abs(adipose.cor)>0.23)
#
adipose.islet <- which(abs(adipose.cor)>0.23)
annot[adipose.islet, c("a_gene_id","gene_symbol")]

###########################################################
#  what is the chance distribution of adipose.cor??
#
#compute corr for 100 permuted versions of Leptin
permuted.cor <- NULL
for(i in 1:100)
{
	y <- sample(phenotypes$Leptin, size = length(phenotypes$Leptin))
	permuted.cor <- cbind(permuted.cor, apply(adipose.rz,2,"mycorr"))
}
#
# find the max corr for each permuted set 
max.cor <- apply(permuted.cor, 2, max)
# find top percentiles
sort(max.cor)[c(80,90,95)]
# it looks like cor > 0.25 is reasonable multiple-test adjust threshold
# nothing that strong in adipose.

###########################################################
###find my Leptin gene expression data 
(gene.names <- grep('^Lep',annot$gene_symbol,value=TRUE))
##Leptin not there, but leptin receptor is
# Hmm...maybe we need to look at other adipokines. 

# Adiponectin and its receptor
(gene.names <- grep('^Adipo',annot$gene_symbol,value=TRUE))

#Visfatin
(gene.names <- grep('^Nampt',annot$gene_symbol,value=TRUE))

#Apelin and its receptor
(gene.names <- grep('^Apln',annot$gene_symbol,value=TRUE))

#Resistin and three related proteins
(gene.names <- grep('^Retn',annot$gene_symbol,value=TRUE))
