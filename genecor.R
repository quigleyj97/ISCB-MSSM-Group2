# genecor.R
# R gene analysis script
# Joe Quigley 2016 for ISCB 2015-16 at MSSM Group 2 with JAX Labs
# Licensed under GNU General Public License 3.0
#                      Free as in Freedom

# depends on ggplot2 and qtl
library(ggplot2)
library(qtl)

# Ensure that the dataset is in the same wd as this script
load(file="BTBR.data.clean.Rdata")

# ======== UTIL FUNCS ========
# Functions to make things  easier. An attempt has been
# made to ensure they are rigourous to make up for R's
# lackadasical type system. This means they check to
# ensure preconditions are met. The checks may have
# edge cases I didn't account for, but in general if they
# fail then a precondition was not met.
# ============================

# Accepts a regexp string 'gene', returns islet expression if the gene exists
# IMPORTANT: Stops if grep doesn't return one gene
# IMPORTANT: Must be run after the dataset is loaded
# Side Effect: Accesses variables not passed as arguments
# returns: islet gene expression
load_islet <- function(gene) {
  # Initial validity checks. If any of these fail, then
  # the preconditions are not met. See comments above.
  stopifnot(exists(deparse(substitute(annot))))
  gene_id <- grep(gene, annot$gene1)
  stopifnot(length(gene_id) == 1)
  islet.rz[, annot[gene_id, 1]]
}

# Same as load_islet but general case
# IMPORTANT: tissue must be valid tissue
load_expression <- function(gene, tissue) {
  stopifnot(exists(deparse(substitute(annot))))
  stopifnot(exists(deparse(substitute(tissue))))
  gene_id <- grep(gene, annot$gene1)
  stopifnot(length(gene_id) == 1)
  tissue[, annot[gene_id, 1]]
}

# Performs a scanone conditioning for a gene
# IMPORTANT: cross must be prepared for scanning before calling this func
# IMPORTANT: gene must exist in the cross
# Side effect: Accesses variables not passed as arguments
# Returns the scan results
scan_cond_gene <- function(gene, pheno) {
  # Ensure the gene expression exists and has been processed
  for (gene_i in gene)  {
    stopifnot(exists(deparse(substitute(cross$pheno[gene_i]))))
  }
  scanone(cross, pheno.col=pheno, addcovar=cross$pheno[gene])
}

# Again, same, but for 
# Prepares the f2g cross object for QTL scans by
# removing unneeded columns and inserting clinical phenotypes
# and gene phenotypes into the cross object. Returns the resulting
# prepared cross object.
# IMPORTANT: All clinical phenotypes must be in phenotypes.rz
# IMPORTANT: All gene phenotypes must already be loaded before being
#            passed to this function.
# IMPORTANT: f2g must exist and be a cross object
# Side Effect: Accesses variables not passed as arguments
# returns: new f2g$pheno
prepare_for_scan <- function(clin_pheno, gene_pheno) {
  # Precondition checks
  stopifnot(exists(deparse(substitute(f2g))))
  stopifnot(all(clin_pheno %in% names(phenotypes.rz)))
  cbind(
    f2g$pheno[,c("MouseNum","Sex","pgm")],
    phenotypes.rz[,clin_pheno],
    gene_pheno)
}

# =====Plotting functions
# Taken from atg.r written by Dr. Gary Churchill and refactored slightly

# Side Effect: modifies variables not passed as arguments
# returns: void
panel.cor <- function(x, y, sigfigs=2, prefix="",
                      sig_threshold=0.5,  cex.cor, ...) {
  # To ensure non-destructive parameter changing
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1)) # Possible issue could come from an edge case here...
  r <- cor(x, y, use = "complete.obs")
  txt <- format(r, digits = sigfigs)
  txt <- paste(prefix, txt, sep = "")
  if(missing(cex.cor)) {
    cex.cor <- 0.8 / strwidth(txt)
  }
  text(0.5, 0.5, txt, cex = cex.cor,
       col = c("gray60", "black")[(abs(r) > sig_threshold) + 1])
}

# Side Effect: modifies variables not passed as arguments
# returns: void
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nbrks <- length(breaks)
  y <- h$counts
  y <- y / max(y)
  # I *think* what this does is just use the first and last break
  # but I don't know enough about the histogram function to be
  # sure
  rect(breaks[-nbrks], 0, breaks[-1], y, col = "cyan", ...)
}

# =====END Plotting functions

# Fit causal models to a triplet with BIC scoring
# Taken from atg.R written by Dr. Gary Churchill
# (also refactored slightly)
# If R is to be believed, this is a pure function
#    But of course R being R we can't prove that
triple.fit <- function(X, Y, Q) {
  # Remove NA rows
  idx <- sort(unique(c(
    which(is.na(X)),
    which(is.na(Y)),
    which(is.na(Q))
  )))
  X <- X[-idx]
  Y <- Y[-idx]
  Q <- Q[-idx]
  # The resulting X, Y, and Qs should not be empty
  stopifnot(length(X) > 0)
  stopifnot(length(Y) > 0)
  stopifnot(length(Q) > 0)
  scores <- c(
    (BIC(lm(X~Q)) + BIC(lm(Y~Q))),    # Independent X<-Q->Y
    (BIC(lm(X~Y)) + BIC(lm(Y~Q))),    # Reactive    Q->Y->X
    (BIC(lm(X~Q)) + BIC(lm(Y~X))),    # Causal      Q->X->Y
    (BIC(lm(X~Q)) + BIC(lm(Y~Q + X))) # Complex
  )
  names(scores) <- c("Independent", "Reactive", "Causal", "Complex")
  scores
}


# ======== USER SPACE ========
# Enter your code here: genes,
# plotting, etc.
# ============================

pheno <- c("INS.10wk", "GLU.10wk")
genes <- c("Irx3$", "Sirt1$")
cond_genes <- c("Myt1l", "Cmpk2", "Cog5", "Colec11", "Efcab10", "Lpin1")
# and others, there's no upper limit other than practical

islet_expression <- lapply(genes, load_islet)
names(islet_expression) <- genes
cond_ids <- match(cond_genes, annot$gene1, nomatch=0)

# Non-destructive, we don't have to reload f2g$pheno if we screw up
cross <- f2g
cross$pheno <- prepare_for_scan(pheno, islet_expression)

# 4 in 4:length(cross) is a magic number, we don't want the first 3 for pairwise scans
pairs(cross$pheno[, 4:length(cross$pheno)],
      upper.panel = panel.cor, diag.panel = panel.hist)

# Sex as a covariate
sex <- as.numeric(cross$pheno$Sex)
cross_gp <- calc.genoprob(cross, step = 1)

# Perform a QTL scan
scan1 <- scanone(cross_gp, pheno.col = c(pheno, genes),
                 method = "hk", addcovar = sex)

perms_filename <- paste(c("BTBR.perms",
                        sort(c(genes, pheno)),
                        "Rdata"), collapse=".")

if(!file.exists(perms_filename))  {
  # ID LOD significance thresholds for gene expression traits
  # 4+length(pheno) is a magic number, we need to throw out the
  # first 3 columns and ignore clinical phenotypes
  perm1 <- scanone(cross_gp, pheno.col = (4 + length(pheno)),
           addcovar = sex, method = "hk", n.perm = 100, perm.Xsp = TRUE)
  save(perm1, file=perms_filename)
} else {
  perm1 <- load(file=perms_filename)
}

# plot 3 rows, 1 col
par(mfrow = c(3,1))
# Not sure what the magic number 2 is for
# I *think* it's for removing the first 2 columns as they're not meant for plotting
#   First two cols are chromosome and 'pos' ( I *think* peak )
for(i in 1:(length(scan1) - 2))  {
  plot(scan1, lodcolumn = i)
  title(main = paste(names(scan1)[i + 2], "expression"), ylab = "LOD")
  add.threshold(scan1, perms = perm1, alpha = 0.05,
                lty = "dashed", lwd = 2, col = "orange")
  add.threshold(scan1, perms = perm1, alpha = 0.10,
                lty = "dashed", lwd = 2, col = "purple")
}

# Add conditional genes to cross object
cross$pheno <- cbind(cross$pheno,
                     islet.rz[,match(cond_genes,
                                     annot$gene1,
                                     nomatch = 0)])

# Bind our phenotypes of interest first
cond_scans <- lapply(scan_cond_genes, cond_gene, c("INS.10wk", "Sirt1$"))

# Magic number exists for the same reason as before
for(i in 1:length(cond_scans) - 2)  {
  plot(cond_scans, lodcolumn=i)
  add.threshold(scan1, perms = perm1, alpha = 0.05,
                lty = "dashed", lwd = 2, col = "orange")
  add.threshold(scan1, perms = perm1, alpha = 0.10,
                lty = "dashed", lwd = 2, col = "purple")
}
