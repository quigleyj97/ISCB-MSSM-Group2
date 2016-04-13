##########################################
#  Demo running QTLNET
#  for a small set of traits
#  GAC last updated: March 21, 2015
#  SMC revised August 5, 2015

library(qtl)
library(qtlnet)
library(ggplot2)

#read cross data into R environment
load("BTBR.clean.data.Rdata")

# Pull some traits into the cross object
f2g$pheno <- f2g$pheno[,c("MouseNum","Sex","pgm")]
f2g$pheno <- transform(f2g$pheno,
            Ucp1 = adipose.rz[,which(annot$gene_symbol == "Ucp1")],
            Cidea = adipose.rz[,which(annot$gene_symbol == "Cidea")],
            Plin5 = adipose.rz[,which(annot$gene_symbol == "Plin5")]
            )
#check
names(f2g$pheno)

# set up for genome scans
f2g <- calc.genoprob(f2g, step = 1, stepwidth = "fixed",
                     map.function = "c-f", err = 0.002)

#choose phenotypes for QTLNET analysis
pheno.col <- c(4:length(f2g$pheno))
names(f2g$pheno)[pheno.col]

# limit the number of edges pointing into a trait
max.parents <- min(length(pheno.col) - 1,4)

# adjust the lod threshold to control number of QTL in model
min.lod <- 3.80

# run the monte carlo chain
# this is a single run - repeats with random start
# the number of samples and thinning rate can be varied
# using sex as "addcov" generated an error
# verbose will let you know it is running
mcmc <- mcmc.qtlnet(f2g, pheno.col,
                    threshold = min.lod,
                    max.parents = max.parents,
                    nsamples=10000, thinning=10,
                    rev.method="nbhd",
                    init.edges=NULL,
                    verbose=TRUE)

# edge probabilities
summary(mcmc)

# plot the model, note uses X11
plot(mcmc)

# plot progress of the chain
x11()
qplot(seq(floor(1000)), mcmc$post.bic)

# tabulate the models that were visited in by the chain
table(mcmc$post.model)
names(f2g$pheno)[pheno.col]
