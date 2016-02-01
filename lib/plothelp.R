# copied from atg.R
# to load, use source()
# or sourceDirectory() the /lib folder
# sourceDirectory() is part of R.utils
# otherwise, you can use:
#	sapply(list.files(pattern="[.]R$", path="lib/", full.names=TRUE), source)

####
# some useful plotting functions
# see documentation for "pairs" function
# note that you can adjust the multiplier on cex.cor
# and the threshold for gray vs black
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r <- cor(x, y, use="complete.obs")
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt):
	text(0.5, 0.5, txt, cex = cex.cor*1.0, col=c("gray60", "black")[(abs(r)>0.5)+1])
}

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
	b4 <- BIC(lm(X~Q)) + BIC(lm(Y~Q+X))	#complex
	scores <- c(b1,b2,b3,b4)
	names(scores) <- c("independent","reactive","causal","complex")
	scores
}
