"coef.fda"<-
function(object, type = c("canonical", "discriminant"), ...)
{
	type <- match.arg(type)
	fit <- object$fit
	Coefs <- fit$coef
	if(is.null(Coefs))
		stop("No explicit coefficients in this formulation")
	dimension <- object$dimension
	Coefs <- Coefs %*% object$theta[, seq(dimension), drop = F]
	lambda <- object$values
	alpha <- sqrt(lambda[seq(dimension)])
	sqima <- sqrt(1 - lambda[seq(dimension)])
	Coefs <- scale(Coefs, F, sqima * alpha)
	if(type == "discriminant")
		Coefs <- Coefs %*% t(object$means)
	Coefs
}
"confusion"<-
function(object, ...)
UseMethod("confusion")
"contr.fda"<-
function(p = rep(1, d[1]), contrast.default = contr.helmert(length(p)))
{
	d <- dim(contrast.default)
	sqp <- sqrt(p/sum(p))
	x <- cbind(1, contrast.default) * outer(sqp, rep(1, d[2] + 1))
	qx <- qr(x)
	J <- qx$rank
	qr.qy(qx, diag(d[1])[, seq(2, J)])/outer(sqp, rep(1, J - 1))
}
"fda"<-
function(formula = formula(data), data = sys.parent(), weights, theta, 
	dimension = J - 1, eps = .Machine$double.eps, method = polyreg, 
	keep.fitted = (n * dimension < 1000), ...)
{
#
# Flexible Discriminant Analysis
# Function for fitting models described in
# Hastie, Tibshirani and Buja, 1994, JASA
# "Flexible Discriminant Analysis by Optimal Scoring"
# Hastie, Buja and Tibshirani, 1995, Annals of Statistics
# "Penalized Discriminant Analysis"
# Modified 2/15/95 by T Hastie
#


  this.call <- match.call()       
  
  m <- match.call(expand = F)
  m[[1]] <- as.name("model.frame")
  m <- m[match(names(m), c("", "formula", "data", "weights"), 0)]
  
  m <- eval(m, sys.parent())
  Terms <- attr(m, "terms")
	g <- model.extract(m, response)
	attr(Terms, "intercept") <- 0
	x <- model.matrix(Terms, m)
	dd <- dim(x)
	n <- dd[1]
	m <- dd[2]	#
	weights <- model.extract(m, weights)
	if(!length(weights))
		weights <- rep(1, n)
	else if(any(weights < 0))
		stop("negative weights not allowed")
	if(length(g) != n)
		stop("g should have length nrow(x)")
	fg <- factor(g)	#
# if some levels are missing, this gets rid of them
# 
	prior <- table(fg)
	prior <- prior/sum(prior)
	cnames <- levels(fg)
	g <- as.numeric(fg)
	J <- length(cnames)	#
#
# construct indicator matrix for response
#
# Initialization uses orthogonal contrasts for theta
# Even if contrasts are supplied, they need to be orthogonalized
	iswt <- F
	if(missing(weights))
		dp <- table(g)/n
	else {
		weights <- (n * weights)/sum(weights)
		dp <- tapply(weights, g, sum)/n
		iswt <- T
	}
	if(missing(theta))
		theta <- contr.helmert(J)
	theta <- contr.fda(dp, theta)
	Theta <- theta[g,  , drop = F]	#
#Theta is now an n x K matrix of responses for the (nonparametric)
# regression, normalized wrt the data (i.e. theta normalized wrt dp)
# where K is min(J-1, ncol of starting theta) 
#
	fit <- method(x, Theta, weights, ...)
	if(iswt)
		Theta <- Theta * weights
	ssm <- t(Theta) %*% fitted(fit)/n	#
# This n could be (n-1) to make it "unbiassed"
#
# now we need the svd of ssm (unweighted)
#
	ed <- svd(ssm, nu = 0)
	thetan <- ed$v
	lambda <- ed$d	
	# Note: the discriminant eigenvalues are a transformation
# of the optimal scaling values 
	lambda[lambda > 1 - eps] <- 1 - eps	#
# If lambda is one we get errors, so we illiminate this possibility
	discr.eigen <- lambda/(1 - lambda)
	pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
	dimension <- min(dimension, sum(lambda > eps))
	if(dimension == 0) {
		warning("degenerate problem; no discrimination")
		return(structure(list(dimension = 0, fit = fit, call = 
			this.call), class = "fda"))
	}
	thetan <- thetan[, seq(dimension), drop = F]	#
	pe <- pe[seq(dimension)]	#Now produce projected centroids
#
	alpha <- sqrt(lambda[seq(dimension)])
	sqima <- sqrt(1 - lambda[seq(dimension)])	#
# package up results
#
	vnames <- paste("v", seq(dimension), sep = "")
	means <- scale(theta %*% thetan, F, sqima/alpha)
	dimnames(means) <- list(cnames, vnames)
	names(lambda) <- c(vnames, rep("", length(lambda) - dimension))	#
	names(pe) <- vnames
	obj <- structure(list(percent.explained = pe, values = lambda, means = 
		means, theta.mod = thetan, dimension = dimension, prior = prior,
		fit = fit, call = this.call, terms = Terms), class = "fda")
	obj$confusion <- confusion(predict(obj), fg)	
	# get rid of the fitted values; these take up too much space
	if(!keep.fitted) obj$fit$fitted.values <- NULL	#
	obj
}
"mars"<-
function(x, y, w = rep(1, nrow(x)), wp, degree = 1, nk = max(21, 2 * ncol(x) + 
	1), penalty = 2, thresh = 0.001, prune = T, trace.mars = F, 
	forward.step = T, prevfit = NULL, ...)
{
	this.call <- match.call()
	if((nk %% 2) != 1)
		nk <- nk - 1
	x <- as.matrix(x)
	np <- dim(x)
	n <- np[1]
	p <- np[2]
	y <- as.matrix(y)
	nclass <- ncol(y)
	if(is.null(np)) {
		np <- c(length(x), 1)
		x <- as.matrix(x)
	}
	if(forward.step) {
		interms <- 1
		lenb <- nk
		bx <- matrix(rep(0, nrow(x) * nk), nrow = n)
		res <- matrix(rep(0, nrow(x) * ncol(y)), nrow = n)
		fullin <- rep(0, nk)
		cuts <- NULL
		factor <- NULL
	}
	else {
		bx <- model.matrix.mars(prevfit, x, full = T)
		interms <- ncol(bx)
		lenb <- prevfit$lenb
		o <- prevfit$all.terms
		fullin <- rep(0, ncol(bx))
		fullin[o] <- 1
		res <- prevfit$res
		factor <- prevfit$factor
		cuts <- prevfit$cuts
		if(missing(penalty))
			penalty <- prevfit$penalty
		degree <- prevfit$degree
		nk <- lenb
		thresh <- prevfit$thresh
	}
	if(missing(penalty) & (degree > 1))
		penalty <- 3
	if(!missing(wp)) {
		if(any(wp <= 0))
			stop("wp should all be positive")
		wp <- sqrt(wp/sum(wp))
		y <- y * outer(rep(1, n), wp)
	}
	else wp <- NULL
	tagx <- x
	storage.mode(tagx) <- "integer"
	for(j in 1:p) {
		tagx[, j] <- order(x[, j])
	}
	bestin <- rep(0, nk)
	flag <- matrix(rep(0, nk * p), nrow = nk, ncol = p)
	if(is.null(cuts))
		cuts <- matrix(rep(0, nk * p), nrow = nk, ncol = p)
	if(is.null(factor)) {
		dir <- matrix(rep(0, nk * p), nrow = nk, ncol = p)
	}
	else {
		dir <- factor
	}
	alpha <- rep(0, nclass)
	beta <- matrix(rep(0, nk * nclass), nrow = nk)
	bestgcv <- 0
	storage.mode(y) <- "double"
	storage.mode(x) <- "double"
	storage.mode(bx) <- "double"
	storage.mode(flag) <- "integer"
	storage.mode(cuts) <- "double"
	storage.mode(dir) <- "double"
	storage.mode(res) <- "double"
	storage.mode(beta) <- "double"
	lenscrat <- 1 + n + 2 * n * nk + 4 * nk * nk + 3 * nk + 3 * nk * nclass +
		3 * nclass + 28 * n + 51	#
#	if(!is.loaded(Fortran.symbol("marss")))
#		dyn.load(.mars.object)
	junk <- .Fortran("marss",
		as.integer(n),
		as.integer(n),
		as.integer(p),
		as.integer(nclass),
		as.matrix(y),
		as.matrix(x),
		as.double(w),
		as.matrix(tagx),
		as.integer(degree),
		as.integer(nk),
		as.double(penalty),
		as.double(thresh),
		as.logical(forward.step),
		as.integer(interms),
		as.logical(prune),
		bx = as.matrix(bx),
		fullin = as.integer(fullin),
		lenb = as.integer(lenb),
		bestgcv = as.double(bestgcv),
		bestin = as.integer(bestin),
		flag = as.matrix(flag),
		cuts = as.matrix(cuts),
		dir = as.matrix(dir),
		res = as.matrix(res),
		alpha = as.double(alpha),
		beta = as.matrix(beta),
		double(lenscrat),
		integer(4 * nk),
		trace.mars)
	lenb <- junk$lenb
	all.terms <- seq(lenb)[junk$fullin[1:lenb] == 1]
	selected.terms <- seq(lenb)[junk$bestin[1:lenb] == 1]
	coefficients <- junk$beta[seq(selected.terms),  , drop = F]
	residuals <- junk$res
	fitted.values <- y - residuals
	if(!is.null(wp)) {
		TT <- outer(rep(1, n), wp)
		residuals <- residuals/TT
		fitted.values <- fitted.values/TT
		coefficients <- coefficients/outer(rep(1, length(selected.terms
			)), wp)
	}
	dir <- junk$dir[seq(lenb),  , drop = F]
	dimnames(dir) <- list(NULL, dimnames(x)[[2]])
	cutss <- junk$cuts[seq(lenb),  , drop = F]
	x <- junk$bx[, selected.terms, drop = F]
	structure(list(call = this.call, all.terms = all.terms, selected.terms
		 = selected.terms, penalty = penalty, degree = degree, nk = nk, 
		thresh = thresh, gcv = junk$bestgcv, factor = dir, cuts = cutss,
		residuals = residuals, fitted.values = fitted.values, lenb = 
		junk$lenb, coefficients = coefficients, x = x), class = "mars")
}
"model.matrix.mars"<-
function(object, x, which = object$selected.terms, full = F, ...)
{
#
# make mars design matrix from output of mars
# x can be any matrix of predictors (without a column of 1s)
#
# if x is missing, the x component of the mars object is returned
# which terms to compute; if not given, the component selected.terms
# from the mars object is used. which is a vector of indexes, ranging
# from 1 thru nrow(object$factor)
# 
# if full =T, a  full matrix (including 0 columns
# or unused terms) is returned
	if(missing(x)) return(object$x)
	x <- as.matrix(x)
	dd <- dim(x)
	n <- dd[1]
	p <- dd[2]
	nterms <- length(which)
	dir <- object$factor
	cut <- object$cuts
	if(full) {
		bx <- matrix(0, nrow = n, ncol = object$lenb)
		bx[, 1] <- 1
	}
	else bx <- matrix(1, nrow = n, ncol = nterms)
	which <- which[-1]
	for(i in seq(along = which)) {
		j <- which[i]
		if(all(dir[j,  ] == 0)) {
			stop("error in factor or which")
		}
		temp1 <- 1
		for(k in 1:p) {
			if(dir[j, k] != 0) {
				temp2 <- dir[j, k] * (x[, k] - cut[j, k])
				temp1 <- temp1 * temp2 * (temp2 > 0)
			}
		}
		if(full)
			bx[, j] <- temp1
		else bx[, i + 1] <- temp1
	}
	bx
}
"polybasis"<-
function(x, degree = 1, monomial = F)
{
	if(degree >= 3)
		warning("This is not a smart polynomial routine. You may get numerical problems with a degree of 3 or more"
			)
	x <- as.matrix(x)
	dn <- dimnames(x)
	dd <- dim(x)
	np <- dd[2]
	if(np == 1)
		monomial <- T
	if(degree > 1) {
		if(monomial) {
			px <- x
			cc <- sapply(split(paste(diag(np)), rep(seq(np), rep(np,
				np))), paste, collapse = "")
			tx <- x
			for(i in 2:degree) {
				px <- px * tx
				x <- cbind(x, px)
				cc <- c(cc, sapply(split(paste(diag(np) * i), 
				  rep(seq(np), rep(np, np))), paste, collapse
				   = ""))
			}
		}
		else {
			matarray <- array(x, c(dd, degree))
			for(i in 2:degree)
				matarray[,  , i] <- x^i
			matarray <- aperm(matarray, c(1, 3, 2))
			x <- matarray[,  , np]
			ad0 <- seq(degree)
			ad <- ad0
			ncol.mat0 <- degree
			ncol.x <- degree
			d0 <- paste(ad0)
			cc <- d0
			for(ii in seq(np - 1, 1)) {
				index0 <- rep(seq(ncol.mat0), ncol.x)
				index <- rep(seq(ncol.x), rep(ncol.mat0, ncol.x
				  ))
				newad <- ad0[index0] + ad[index]
				retain <- newad <= degree
				mat0 <- matarray[,  , ii]
				if(any(retain))
				  newmat <- mat0[, index0[retain], drop = F] * 
				    x[, index[retain], drop = F]
				else newmat <- NULL
				ddn <- paste(d0[index0[retain]], cc[index[
				  retain]], sep = "")
				zeros <- paste(rep(0, nchar(cc[1])), collapse
				   = "")
				cc <- paste(0, cc, sep = "")
				d00 <- paste(d0, zeros, sep = "")
				x <- cbind(mat0, x, newmat)
				cc <- c(d00, cc, ddn)
				ad <- c(ad0, ad, newad[retain])
				ncol.x <- length(ad)
			}
		}
		if(!is.null(dn))
			dn[[2]] <- cc
		else dn <- list(NULL, cc)
		dimnames(x) <- dn
	}
	cbind(Intercept = 1, x)
}
"polyreg"<-
function(x, y, w, degree = 1, monomial = F, ...)
{
	x <- polybasis(x, degree, monomial)
	if(iswt <- !missing(w)) {
		if(any(w <= 0))
			stop("only positive weights")
		w <- sqrt(w)
		y <- y * w
		x <- x * w
	}
	qrx <- qr(x)
	coef <- as.matrix(qr.coef(qrx, y))
	fitted <- qr.fitted(qrx, y)
	if((df <- qrx$rank) < ncol(x))
		coef[qrx$pivot,  ] <- coef
	if(iswt)
		fitted <- fitted/w
	structure(list(fitted.values = fitted, coefficients = coef, degree = 
		degree, monomial = monomial, df = df), class = "polyreg")
}
"pplot"<-
function(x, g, colors, pch, add = F, type = "p", ...)
{
	g <- as.factor(g)
	cc <- codes(g)
	np <- seq(levels(g))
	if(missing(colors))
		colors <- np + 1
	else colors <- rep(colors, length = length(np))
	if(missing(pch))
		pch <- paste(np)
	else pch <- rep(pch, length = length(np))
	if(!add)
		plot(x, type = "n", ...)
	for(i in unique(cc))
		points(x[cc == i,  , drop = F], col = colors[i], pch = pch[i], 
			type = type)
}
"predict.fda"<-
function(object, x, type = c("class", "variates", "posterior", "hierarchical", 
	"distances"), prior, dimension = J - 1)
{
	dist <- function(x, mean, m = ncol(mean))
	(scale(x, mean, F)^2) %*% rep(1, m)
	type <- match.arg(type)
	means <- object$means
	Jk <- dim(means)
	J <- Jk[1]
	k <- Jk[2]	#
# Note for type=="hierarchical" dimension can be a vector
#
	if(type == "hierarchical") {
		if(missing(dimension))
			dimension.set <- seq(k)
		else {
			dimension.set <- dimension[dimension <= k]
			if(!length(dimension.set))
				dimension.set <- k
			dimension <- max(dimension.set)
		}
	}
	else dimension <- min(max(dimension), k)
	if(missing(x))
		y <- predict(object$fit)
	else {
		if(inherits(x, "data.frame") || is.list(x)) {
			Terms <- delete.response(terms(object))
			attr(Terms, "intercept") <- 0
			x <- model.matrix(Terms, x)
		}
		y <- predict(object$fit, x)
	}
	y <- y %*% object$theta[, seq(dimension), drop = F]
	lambda <- object$values
	alpha <- sqrt(lambda[seq(dimension)])
	sqima <- sqrt(1 - lambda[seq(dimension)])
	x <- scale(y, F, sqima * alpha)
	if(missing(prior))
		prior <- object$prior
	else {
		if(any(prior < 0) | round(sum(prior), 5) != 1)
			stop("innappropriate prior")
	}
	means <- means[, seq(dimension), drop = F]
	switch(type,
		variates = return(x),
		class = {
			n <- nrow(x)
			prior <- 2 * log(prior)
			mindist <- dist(x, means[1,  ], dimension) - prior[1]
			pclass <- rep(1, n)
			for(i in seq(2, J)) {
				ndist <- dist(x, means[i,  ], dimension) - 
				  prior[i]
				l <- ndist < mindist
				pclass[l] <- i
				mindist[l] <- ndist[l]
			}
			levels(pclass) <- dimnames(means)[[1]]
			return(factor(pclass))
		}
		,
		posterior = {
			pclass <- matrix(0, nrow(x), J)
			for(i in seq(J))
				pclass[, i] <- exp(-0.5 * dist(x, means[i,  ], 
				  dimension)) * prior[i]
			dimnames(pclass) <- list(dimnames(x)[[1]], dimnames(
				means)[[1]])
			return(pclass/drop(pclass %*% rep(1, J)))
		}
		,
		hierarchical = {
			prior <- 2 * log(prior)
			Pclass <- vector("list", length(dimension.set))
			names(Pclass) <- paste("D", dimension.set, sep = "")
			for(ad in seq(along = dimension.set)) {
				d <- dimension.set[ad]
				dd <- seq(d)
				mindist <- dist(x[, dd, drop = F], means[1, dd, 
				  drop = F], d) - prior[1]
				pclass <- rep(1, nrow(x))
				for(i in seq(2, J)) {
				  ndist <- dist(x[, dd, drop = F], means[i, dd, 
				    drop = F], d) - prior[i]
				  l <- ndist < mindist
				  pclass[l] <- i
				  mindist[l] <- ndist[l]
				}
				levels(pclass) <- dimnames(means)[[1]]
				Pclass[[ad]] <- pclass
			}
			rownames <- dimnames(x)[[1]]
			if(is.null(rownames))
				rownames <- paste(seq(nrow(x)))
			return(structure(Pclass, class = "data.frame", 
				row.names = rownames, dimensions = 
				dimension.set))
		}
		,
		distances = {
			dclass <- matrix(0, nrow(x), J)
			for(i in seq(J))
				dclass[, i] <- dist(x, means[i,  ], dimension)
			dimnames(dclass) <- list(dimnames(x)[[1]], dimnames(
				means)[[1]])
			return(dclass)
		}
		)
}
"predict.mars"<-
function(object, x)
{
#
# computes fitted values for design points x, based on mars fit 
# in object
#
	if(missing(x)) {
		z <- fitted(object)
		if(is.null(z))
			stop("need to supply x")
		else return(z)
	}
	model.matrix.mars(object, x) %*% object$coefficients
}
"predict.polyreg"<-
function(object, x, ...)
{
	if(missing(x)) {
		z <- fitted(object)
		if(is.null(z))
			stop("need to supply x")
		else return(z)
	}
	degree <- object$degree
	monomial <- object$monomial
	polybasis(x, degree, monomial) %*% object$coef
}
"print.fda"<-
function(x, ...)
{
	if(!is.null(cl <- x$call)) {
		cat("Call:\n")
		dput(cl)
	}
	cat("\nDimension:", format(x$dimension), "\n")
	cat("\nPercent Between-Group Variance Explained:\n")
	print(round(x$percent, 2))
	error <- x$confusion
	df <- x$fit
	if(!is.null(df))
		df <- df$df
	if(!is.null(df)) {
		cat("\nDegrees of Freedom (per dimension):", format(sum(df)), 
			"\n")
	}
	if(!is.null(error)) {
		n <- as.integer(sum(error))
		error <- format(round(attr(error, "error"), 5))
		cat("\nTraining Misclassification Error:", error, "( N =", n, 
			")\n")
	}
	invisible(x)
}
"softmax"<-
function(x, gap = F)
{
	d <- dim(x)
	maxdist <- x[, 1]
	pclass <- rep(1, d[1])
	for(i in seq(2, d[2])) {
		l <- x[, i] > maxdist
		pclass[l] <- i
		maxdist[l] <- x[l, i]
	}
	dd <- dimnames(x)[[2]]
	if(gap) {
		x <- abs(maxdist - x)
		x[cbind(seq(d[1]), pclass)] <- drop(x %*% rep(1, d[2]))
		gaps <- do.call("pmin", data.frame(x))
	}
	pclass <- if(is.null(dd) || !length(dd)) pclass else factor(pclass, 
			levels = seq(d[2]), labels = dd)
	if(gap)
		list(class = pclass, gaps = gaps)
	else pclass
}
".bruto.object"<-
"/fs/trevor/docs/discr/BRUTO.o"
".mars.object"<-
"/n/rice/usr/trevor/docs/discr/MARS.o"
".mda.object"<-
"/usr/trevor/docs/mda/mda.o"
"make.dumpdata.mda"<-
function()
dump(c("bruto", "coef.fda", "confusion", "contr.fda", "fda", "mars", 
	"model.matrix.mars", "polybasis", "polyreg", "pplot", "predict.bruto", 
	"predict.fda", "predict.mars", "predict.polyreg", "print.fda", 
	"softmax", ".bruto.object", ".mars.object", ".mda.object", 
	"make.dumpdata.mda", "confusion.default", "confusion.fda", 
	"confusion.list", "fix.theta", "gen.ridge", "predict.gen.ridge", 
	"laplacian", "kmeans.start", "llmult", "lvq.start", "mda", "mda.fit", 
	"mda.means", "mda.start", "mean.penalty", "plot.fda", "pplot4", 
	"predict.mda", "print.mda", "shrink", "shrink.mda", "transform.penalty"
	), "dumpdata.mda")
"confusion.default"<-
function(predict, true, ...)
{
	if(inherits(predict, "data.frame"))
		confusion.list(predict, true)
	else {
		jt <- table(predict, true)
		jd <- dimnames(jt)
		jn <- unlist(jd)
		ju <- jn[duplicated(jn)]
		j1 <- jd[[1]][!match(jd[[1]], ju, 0)]
		j2 <- jd[[2]][!match(jd[[2]], ju, 0)]
		jt <- jt[c(ju, j1), c(ju, j2), drop = F]
		realjt <- jt[ju, ju, drop = F]
		ntot <- sum(jt)
		mismatch <- (ntot - sum(realjt))/ntot
		structure(jt, error = (1 - sum(diag(realjt))/sum(realjt)), 
			mismatch = if(mismatch > 0) mismatch else NULL)
	}
}
"confusion.fda"<-
function(object, data, ...)
{
	if(missing(data))
		return(object$confusion)
	Terms <- terms(object)
	attr(Terms, "intercept") <- 0
	m <- model.frame(Terms, data)
	x <- model.matrix(Terms, m)
	g <- model.extract(m, response)
	confusion.default(predict(object, x, ...), g)
}
"confusion.list"<-
function(pred, truth)
{
	dd <- names(pred)
	y <- seq(dd)
	x <- attr(pred, "dimension")
	if(!length(x))
		x <- seq(dd)
	for(i in y) {
		confi <- confusion(pred[, i], truth)
		y[i] <- attr(confi, "error")
	}
	return(x, y)
}
"fix.theta"<-
function(theta, Q)
{
	M <- t(theta) %*% Q %*% theta
	eM <- eigen(M, sym = T)
	scale(theta %*% eM$vectors, F, sqrt(eM$values))
}
"gen.ridge"<-
function(x, y, weights, lambda = 1, omega, df, ...)
{
#
# ||Y-XB||^2 + lambda*B^T Omega Beta
#    = ||Y-X* B*||^2 + lambda ||B*||^2
# where X*=XU*, Omega=UDU^T,  U*=UD^{-1/2} and B*=D^{1/2}U^T B
# This is a simple ridge problem now
# Let X* = UDV^T
# then H = UDV^T(V^TDV + lambda I)^{-1} VDU^T
#        = UD(D^2 + lambda I)^{-1}DU^T
# with trace sum(d_j^2/(lambda+d_j^2)
# The coefficients B* are given by V(D^2 + lambda I)^{-1}DU^TY
# and then we see that U*B*=B
#
	if(missing(df) && lambda <= .Machine$double.eps) return(polyreg(x, y))
	d <- dim(x)
	mm <- apply(x, 2, mean)
	x <- scale(x, mm, F)
	simple <- if(missing(omega)) T else F
	if(!simple) {
		if(!all(match(c("values", "vectors"), names(omega), F)))
			stop("You must supply an  eigen-decomposed version of omega"
				)
		vals <- pmax(.Machine$double.eps, sqrt(omega$values))
		basis <- scale(omega$vectors, F, vals)
		x <- x %*% basis
	}
	svd.x <- svd(x)
	dd <- svd.x$d
	if(!missing(df)) {
		while(sum(dd^2/(dd^2 + lambda)) > df) lambda <- lambda * 10
		junk <- df.inv(dd^2, df, lambda)
		lambda <- junk$lambda
		df <- junk$df
	}
	else df <- sum(dd^2/(dd^2 + lambda))
	y <- (t(t(y) %*% svd.x$u) * dd)/(dd^2 + lambda)
	coef <- svd.x$v %*% y
	fitted <- x %*% coef
	if(!simple)
		coef <- basis %*% coef
	structure(list(fitted.values = fitted, coefficients = coef, df = df, 
		lambda = lambda, xmeans = mm), class = "gen.ridge")
}
"predict.gen.ridge"<-
function(object, x, ...)
{
	if(missing(x))
		fitted(object)
	else scale(x, object$xmeans, F) %*% object$coef
}
"laplacian"<-
function(size = 16, compose = F)
{
#build gamma
#build gamma
# Here I follow very closely the material on page 635 in JASA 1991
# of O'Sullivan's article on discretized Laplacian Smoothing
	gmat <- matrix(0, size, size)
	xx <- seq(size)
	for(v in xx)
		gmat[, v] <- sqrt(2/size) * cos(((xx - 0.5) * pi * (v - 1))/
			size)
	gmat[, 1] <- gmat[, 1]/sqrt(2)
	lvec <-  - (2 * size^2) * (1 - cos(((xx - 1) * pi)/size))
	gmat <- kronecker(gmat, gmat)
	lvec <- rep(lvec, rep(size, size)) + rep(lvec, size)
	if(compose)
		gmat %*% (lvec^2 * t(gmat))
	else list(vectors = gmat, values = lvec^2)
}
"kmeans.start"<-
function(x, g, subclasses)
{
	cnames <- levels(g <- factor(g))
	J <- length(cnames)
	g <- as.numeric(g)
	weights <- as.list(cnames)
	names(weights) <- cnames
	subclasses <- rep(subclasses, length = length(cnames))	#
	R <- sum(subclasses)
	cl <- rep(seq(J), subclasses)
	cx <- x[seq(R),  , drop = F]
#	if(!is.loaded(symbol.For("kmns")))
#		dyn.load("mda.o")
	for(j in seq(J)) {
		nc <- subclasses[j]
		which <- cl == j
		if(nc <= 1) {
			cx[which,  ] <- apply(x[g == j,  , drop = F], 2, mean)
			wmj <- matrix(1, sum(g == j), 1)
		}
		else {
			xx <- x[g == j,  , drop = F]
			start <- xx[sample(1:nrow(xx), size = nc),  ]
			TT <- kmeans(xx, start)
			cx[which,  ] <- TT$centers
			wmj <- diag(nc)[TT$cluster,  ]
		}
		dimnames(wmj) <- list(NULL, paste("s", seq(nc), sep = ""))
		weights[[j]] <- wmj
	}
	list(x = cx, cl = factor(cl, labels = cnames), weights = weights)
}
"llmult"<-
function(p, g)
{
	index <- cbind(seq(along = g), as.numeric(g))
	p <- p[index]
	-2 * sum(log(p[p > .Machine$double.eps]))
}
"lvq.start"<-
function(x, g, subclasses)
{
	cnames <- levels(fg <- factor(g))
	J <- length(cnames)
	g <- as.numeric(g)
	weights <- as.list(cnames)
	names(weights) <- cnames
	subclasses <- rep(subclasses, length = length(cnames))	#
	size <- sum(subclasses)
#	if(!is.loaded(symbol.For("olvq")))
#		dyn.load("mda.o")
	cb <- lvqinit(x, g, size = size)
	TT <- olvq1(x, g, codebk = cb)
	TT <- lvq3(x, g, codebk = TT)
	cl <- as.numeric(TT$cl)
	R <- length(cl)
	cx <- TT$x
	p <- ncol(cx)	#Compute the weights based on within class assignments
#
	for(j in seq(J)) {
		which <- cl == j
		number <- sum(which)
		if(number == 0) {
			cx <- rbind(cx, apply(x[g == j,  ], 2, mean))
			cl <- c(cl, j)
			wmj <- matrix(1, sum(g == j), 1)
			number <- 1
		}
		else if(number == 1)
			wmj <- matrix(1, sum(g == j), 1)
		else {
			jcx <- cx[which,  ]
			jcl <- seq(number)
			jcluster <- lvqtest(list(x = jcx, cl = jcl), x[g == j,  
				])
			needed <- unique(jcluster)
			rcl <- rep(0, number)
			rcl[needed] <- j
			cl[which] <- rcl
			wmj <- diag(number)[jcluster, needed, drop = F]
			number <- length(needed)
		}
		dimnames(wmj) <- list(NULL, paste("s", seq(number), sep = ""))
		weights[[j]] <- wmj
	}
	TT <- cl > 0
	list(x = cx[TT,  , drop = F], cl = factor(cl[TT], labels = cnames), 
		weights = weights)
}
"mda"<-
function(formula = formula(data), data = sys.parent(), subclasses = 3, sub.df
	 = NULL, tot.df = NULL, dimension = sum(subclasses) - 1, eps = .Machine$
	double.eps, iter = 5, weights = mda.start(x, g, subclasses, trace, ...),
	method = polyreg, keep.fitted = (n * dimension < 1000), trace = F, ...
	)
{
#
# Mixture Discriminant Analysis
# Function for fitting models described in
# Hastie and Tibshirani 1995 JRSSB
# modified by Trevor on 2/15/95
#
	this.call <- match.call()	#
# This extracts the x and g from the formula or data frame
#
# -------< not for human consumption >--------
	m <- match.call(expand = F)
	m[[1]] <- as.name("model.frame")
	m <- m[match(names(m), c("", "formula", "data"), 0)]
	m <- eval(m, sys.parent())
	Terms <- attr(m, "terms")
	g <- model.extract(m, response)
	attr(Terms, "intercept") <- 0
	x <- model.matrix(Terms, m)
	dd <- dim(x)
	n <- dd[1]
	m <- dd[2]	#
#
# Define a function needed later
	rowmin <- function(mat)
	{
		ncc <- ncol(mat)
		if(ncc == 1)
			return(drop(mat))
		rowm <- pmin(mat[, 1], mat[, 2])
		if(ncc == 2)
			return(rowm)
		else {
			for(j in seq(from = 3, to = ncc))
				rowm <- pmin(rowm, mat[, j])
		}
		rowm
	}
# ----------------------------------------------
#
	if(length(g) != n) stop("g should have length nrow(x)")	#
# turn g into a factor
# if some levels are missing, this gets rid of them
	fg <- factor(g)	#	
# weights is a special beast. It is a list of matrices of probabilities
# that are the subclass probabilites. The names correspond to the levels
# of the (implicit) factor g
#
# One can supply an mda object itself as weights
#
	if(inherits(weights, "mda")) {
		if(is.null(weights$weights))
			weights <- predict(weights, x, type = "weights", g = fg
				)
		else weights <- weights$weights
	}
	subclasses <- sapply(weights, ncol)	#
# extract codes and class labels
# I don't use codes since the class labels can get mixed up
#
	prior <- table(fg)
	dim(prior) <- NULL
	prior <- prior/sum(prior)
	cnames <- levels(fg)
	g <- as.numeric(fg)	#
	J <- length(cnames)	#
	Assign <- split(seq(sum(subclasses)), rep(seq(J), subclasses))
	names(Assign) <- cnames	#
#
# see if shrinking was called for
#	
	if(!is.null(tot.df)) {
		if(tot.df >= sum(subclasses))
			tot.df <- NULL
	}
	if(!is.null(sub.df)) {
		sub.df <- rep(sub.df, length = length(prior))
		sub.df <- pmin(sub.df, subclasses)
		if(all(sub.df == subclasses))
			sub.df <- NULL
	}
# generate starting sub-class weight
	for(counter in seq(iter)) {
		fit <- mda.fit(x, g, weights, assign.theta = Assign, Rj = 
			subclasses, sub.df = sub.df, tot.df = tot.df, dimension
			 = dimension, eps = .Machine$double.eps, method = 
			method, trace = trace, ...)	#
# predict.fda works on a mda.fit object
#
		dmat <- predict.fda(fit, type = "distance")
		sub.prior <- fit$sub.prior	#
# Now recompute weights
#
		for(j in seq(J)) {
			TT <- dmat[g == j, Assign[[j]], drop = F]
			TT <- exp(-0.5 * (TT - rowmin(TT)))
			TT <- TT * outer(rep(1, nrow(TT)), sub.prior[[j]])
			weights[[j]][] <- TT/drop(TT %*% rep(1, ncol(TT)))
		}
#
#Compute posterior  probabilities
#
		pclass <- matrix(1, n, J)
		dmat <- exp(-0.5 * (dmat - rowmin(dmat)))
		for(j in seq(J)) {
			priorj <- sub.prior[[j]]
			ass <- Assign[[j]]
			TT <- dmat[, ass, drop = F] * outer(rep(1, n), priorj)
			TTot <- drop(TT %*% rep(1, length(ass)))
			pclass[, j] <- prior[j] * TTot
		}
		pclass <- pclass/drop(pclass %*% rep(1, J))
		if(trace)
			cat("Iteration", counter, "\tDeviance(multinomial)", 
				format(round(ll <- llmult(pclass, g), 5)), "\n"
				)
	}
	if(!trace) ll <- llmult(pclass, g)	
	# get rid of the fitted values; these take up too much space
	if(!keep.fitted) fit$fit$fitted.values <- NULL	#
	dimnames(pclass) <- list(NULL, names(Assign))
	conf <- confusion(softmax(pclass), fg)
	fit <- c(fit, list(weights = weights, prior = prior, assign.theta = 
		Assign, deviance = ll, confusion = conf, terms = Terms))
	fit$call <- this.call
	fit$sub.df <- sub.df
	fit$tot.df <- tot.df
	class(fit) <- c("mda", "fda")
	fit
}
"mda.fit"<-
function(x, g, weights, theta, assign.theta, Rj, sub.df, tot.df, dimension = R - 
	1, eps = .Machine$double.eps, method = polyreg, ...)
{
	this.call <- match.call()
	n <- nrow(x)
	cnames <- names(weights)
	J <- length(cnames)	#
# now extract the lengths and names of the appropriate matrices in
# weights
	R <- sum(Rj)	#
# Get the total weight for each subclass (to be used in weighting the scores)
# Keep it in list form
#
	wtots <- lapply(weights, function(x)
	apply(x, 2, sum))	#
	sub.prior <- lapply(wtots, function(x)
	x/sum(x))	#
#dp is the unconditional probability of being in sublcass r of class j
#
	dp <- unlist(wtots)
	subclass.names <- names(dp)
	dp <- dp/sum(dp)	#
# Initialization uses orthogonal contrasts for theta
# Even if contrasts are supplied, they need to be orthogonalized
	if(missing(theta))
		theta <- contr.helmert(R)
	theta <- contr.fda(dp, theta)	#
# Here is the fix for shrinking
#
# 
	if(!(is.null(sub.df) & is.null(tot.df))) {
		Q <- diag(dp) + transform.penalty(prior = dp, cl = rep(seq(J), 
			Rj), df = sub.df, tot.df = tot.df)
		theta <- fix.theta(theta, Q)
	}
#
# Now compute the weighted score means for each observation
#
	Theta <- matrix(0, n, R - 1)
	obs.weights <- double(n)
	for(j in seq(J)) {
		Theta[g == j,  ] <- weights[[j]] %*% theta[assign.theta[[j]],  ,
			drop = F]
		obs.weights[g == j] <- weights[[j]] %*% rep(1, Rj[j])
	}
#Theta is now an n x K matrix of responses for the regression, normalized wrt
# the weights (i.e. theta normalized wrt dp)
# where K is min(R-1, ncol of starting theta) 
#
#
	fit <- method(x, Theta, obs.weights, ...)
	Theta <- Theta * obs.weights
	ssm <- t(Theta) %*% fitted(fit)/n	#
# now we need the svd of ssm (unweighted)
#
	ed <- svd(ssm, nu = 0)
	thetan <- ed$v
	lambda <- ed$d	
	# Note: the discriminant eigenvalues are a transformation
# of the optimal scaling values 
	lambda[lambda > 1 - eps] <- 1 - eps	#
# If lambda is one we get errors, so we illiminate this possibility
	discr.eigen <- lambda/(1 - lambda)
	pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
	dimension <- min(dimension, sum(lambda > eps))
	if(dimension == 0) {
		warning("degenerate problem; no discrimination")
		return(structure(list(dimension = 0, fit = fit, call = 
			this.call), class = "fda"))
	}
	thetan <- thetan[, seq(dimension), drop = F]	#
	pe <- pe[seq(dimension)]	#Now produce projected centroids
#
	alpha <- sqrt(lambda[seq(dimension)])
	sqima <- sqrt(1 - lambda[seq(dimension)])	#
# package up results
#
	vnames <- paste("v", seq(dimension), sep = "")
	means <- scale(theta %*% thetan, F, sqima/alpha)
	dimnames(means) <- list(subclass.names, vnames)
	names(lambda) <- c(vnames, rep("", length(lambda) - dimension))	#
	names(pe) <- vnames
	structure(list(percent.explained = pe, values = lambda, means = means, 
		theta.mod = thetan, dimension = dimension, sub.prior = 
		sub.prior, fit = fit, call = this.call))
}
"mda.means"<-
function(object, x, y)
{
	weights <- means <- object$weights
	nn <- names(object$weights)
	for(i in nn) {
		xx <- x[y == i,  ]
		ww <- weights[[i]]
		nc <- ncol(ww)
		xm <- matrix(0, ncol(x), nc)
		for(j in seq(nc)) {
			www <- ww[, j]
			www <- www/sum(www)
			xm[, j] <- apply(xx * www, 2, sum)
		}
		means[[i]] <- xm
	}
	means
}
"mda.start"<-
function(x, g, subclasses = 3, trace.mda.start = F, start.method = c("kmeans", 
	"lvq"), tries = 5, criterion = c("misclassification", "deviance"), ...
	)
{
#	if(!length(find("lvqtest")))
#		stop("mda requires functions in the classif collection donated by Brian Ripley to statlib/S"
#			)
	start.method <- match.arg(start.method)
#	if((start.method == "kmeans") && !length(find("kmeans")))
#		stop("mda with start.method=kmeans requires the kmeans() function, currently only available in Splus"
#			)
	criterion <- match.arg(criterion)
	name.criterion <- switch(criterion,
		misclassification = "Misclassification Error",
		deviance = "Deviance(multinomial)")
	starter <- get(paste(start.method, "start", sep = "."), mode = 
		"function")
	fg <- factor(g)
	cnames <- levels(fg)
	prior <- table(fg)
	J <- length(cnames)
	n <- length(g)
	g <- as.numeric(fg)	# Now loop over tries
#
	best.ll <- 1/.Machine$double.eps
	for(try in seq(tries)) {
		start <- starter(x, fg, subclasses)
		weights <- start$weights
		if(criterion == "misclassification") {
			pg <- lvqtest(start, x)
			ll <- attr(confusion(pg, fg), "error")
		}
		else {
			subclasses <- sapply(weights, ncol)	#
			Assign <- split(seq(sum(subclasses)), rep(seq(J), 
				subclasses))
			names(Assign) <- cnames	#
			fit <- mda.fit(x, g, weights, assign.theta = Assign, Rj
				 = subclasses, eps = .Machine$double.eps, 
				method = polyreg, ...)	#
# predict.fda works on a mda.fit object
#
			dmat <- exp(-0.5 * predict.fda(fit, type = "distance"))
			sub.prior <- fit$sub.prior	#
#
#Now compute the weights and the posteriors
#
			pclass <- matrix(1, n, J)
			for(j in seq(J)) {
				priorj <- sub.prior[[j]]
				ass <- Assign[[j]]
				TT <- dmat[, ass, drop = F]
				TT <- TT * outer(rep(1, n), priorj)
				TTot <- drop(TT %*% rep(1, length(ass)))
				wmj <- TT[g == j,  , drop = F]/TTot[g == j]
				pclass[, j] <- prior[j] * TTot
				dimnames(wmj) <- list(NULL, paste("s", seq(
				  along = ass), sep = ""))
				weights[[j]] <- wmj
			}
			pclass <- pclass/drop(pclass %*% rep(1, J))
			ll <- llmult(pclass, g)
		}
		if(trace.mda.start)
			cat(start.method, "start   \t", name.criterion, format(
				round(ll, 5)), "\n")
		if(ll < best.ll) {
			keep.weights <- weights
			best.ll <- ll
		}
	}
	structure(keep.weights, criterion = best.ll, name.criterion = 
		name.criterion)
}
"mean.penalty"<-
function(prior, cl)
{
	n <- length(prior)
	Q <- matrix(0, n, n)
	ll <- unique(cl)
	for(lll in ll) {
		which <- cl == lll
		pp <- prior[which]
		npp <- length(pp)
		pp <- diag(npp) - outer(rep(1, npp), pp/sum(pp))
		Q[which, which] <- t(pp) %*% pp
	}
	attr(Q, "cl") <- cl
	attr(Q, "prior") <- prior
	Q
}
"plot.fda"<-
function(object, data, g, coords = c(1, 2), group = c("true", "predicted"), 
	colors, pch, mcolors = max(colors) + 1, mpch, pcex = 0.5, mcex = 2.5, 
	...)
{
	group <- match.arg(group)
	if(missing(data)) {
		vars <- predict(object, type = "var")
		g <- predict(object)
		group <- "predict"
	}
	else {
		ff <- terms(object)
		attr(ff, "intercept") <- 0
		m <- model.frame(ff, data)
		x <- model.matrix(ff, m)
		vars <- predict(object, x, type = "var")
		if(group == "true")
			g <- model.extract(m, response)
		else g <- predict(object, x)
	}
	means <- object$means
	g <- as.factor(g)
	cc <- as.numeric(g)
	np <- seq(levels(g))
	if(missing(colors))
		colors <- np
	else colors <- rep(colors, length = length(np))
	if(missing(pch))
		pch <- paste(np)
	else pch <- rep(paste(pch), length = length(np))
	mcolors <- rep(mcolors, length = length(np))
	if(missing(mpch))
		mpch <- pch
	else mpch <- rep(paste(mpch), length = length(np))
	assign <- object$assign
	if(is.null(assign))
		assign <- split(seq(np), seq(np))
	if(!is.matrix(coords)) {
		coords <- matrix(coords, length(coords), length(coords))
		tt <- lower.tri(coords)
		coords <- cbind(t(coords)[tt], coords[tt])
	}
	for(ii in seq(nrow(coords))) {
		coord.pair <- coords[ii,  ]
		plot(rbind(vars[, coord.pair], means[, coord.pair]), ..., type
			 = "n", xlab = paste("Discriminant Var", coord.pair[1]),
			ylab = paste("Discriminant Var", coord.pair[2]), main
			 = paste("Discriminant Plot for", group, "classes"))
		for(i in np) {
			which <- cc == i
			if(any(which))
				points(vars[which, coord.pair, drop = F], col
				   = colors[i], pch = pch[i], cex = pcex)
			points(means[assign[[i]], coord.pair, drop = F], col = 
				mcolors[i], pch = 1, cex = mcex)
			points(means[assign[[i]], coord.pair, drop = F], col = 
				mcolors[i], pch = mpch[i], cex = mcex/2)
		}
	}
	invisible()
}
"pplot4"<-
function(x, ...)
{
	par(mfrow = c(3, 2))
	for(i in 1:3) {
		for(j in (i + 1):4)
			pplot(x[, c(i, j)], xlab = paste("var", i), ylab = 
				paste("var", j), ...)
	}
}
"predict.mda"<-
function(object, x, type = c("class", "variates", "posterior", "hierarchical", 
	"weights"), prior = NULL, dimension = R - 1, g, ...)
{
	type <- match.arg(type)	#
# Note for type=="hierarchical" dimension can be a vector
#
	Rk <- dim(object$means)
	R <- Rk[1]
	k <- Rk[2]
	if(type == "hierarchical") {
		if(missing(dimension))
			dimension.set <- seq(k)
		else {
			dimension.set <- dimension[dimension <= k]
			if(!length(dimension.set))
				dimension.set <- k
			dimension <- max(dimension.set)
		}
		Pclass <- vector("list", length(dimension.set))
		names(Pclass) <- paste("D", dimension.set, sep = "")
		for(ad in seq(along = dimension.set)) {
			d <- dimension.set[ad]
			Pclass[[ad]] <- if(missing(x)) Recall(object, prior = 
				  prior, dimension = d, ...) else Recall(object,
				  x, prior = prior, dimension = d, ...)
		}
		rownames <- names(Pclass[[1]])
		if(is.null(rownames))
			rownames <- paste(seq(along = Pclass[[1]]))
		return(structure(Pclass, class = "data.frame", row.names = 
			rownames, dimensions = dimension.set))
	}
	else dimension <- min(max(dimension), k)
	if(is.null(prior))
		prior <- object$prior
	else {
		if(any(prior < 0) | round(sum(prior), 5) != 1)
			stop("innappropriate prior")
	}
	if(type == "variates") return(NextMethod("predict"))	
	# Define a function needed later
	rowmin <- function(mat)
	{
		ncc <- ncol(mat)
		if(ncc == 1)
			return(drop(mat))
		rowm <- pmin(mat[, 1], mat[, 2])
		if(ncc == 2)
			return(rowm)
		else {
			for(j in seq(from = 3, to = ncc))
				rowm <- pmin(rowm, mat[, j])
		}
		rowm
	}
	dmat <- if(missing(x)) predict.fda(object, type = "distances", 
			dimension = dimension, ...) else predict.fda(object, x, 
			type = "distances", dimension = dimension, ...)
	Assign <- object$assign
	sub.prior <- object$sub.prior
	J <- length(Assign)
	if(type == "weights") {
		if(missing(x))
			return(object$weights)
		g <- as.numeric(g)
		weights <- Assign
		for(j in seq(J)) {
			TT <- dmat[g == j, Assign[[j]], drop = F]
			TT <- exp(-0.5 * (TT - rowmin(TT)))
			TT <- TT * outer(rep(1, nrow(TT)), sub.prior[[j]])
			weights[[j]] <- TT/drop(TT %*% rep(1, ncol(TT)))
		}
		return(weights)
	}
	pclass <- matrix(1, nrow(dmat), J)
	dmat <- exp(-0.5 * (dmat - rowmin(dmat)))
	for(j in seq(J)) {
		TT <- dmat[, Assign[[j]], drop = F]
		TT <- TT * outer(rep(1, nrow(TT)), sub.prior[[j]])
		pclass[, j] <- prior[j] * drop(TT %*% rep(1, ncol(TT)))
	}
#	dimnames(pclass) <- list(dimnames(x)[[1]], names(Assign))
	dimnames(pclass) <- list(NULL, names(Assign))
	switch(type,
		class = softmax(pclass),
		posterior = pclass/drop(pclass %*% rep(1, J)))
}
"print.mda"<-
function(x, ...)
{
	NextMethod("print")
	if(!is.null(x$deviance))
		cat("\nDeviance:", format(round(x$deviance, 3)), "\n")
	invisible(x)
}
"shrink"<-
function(object, ...)
UseMethod("shrink")
"shrink.mda"<-
function(object, sub.df = NULL, tot.df = NULL, ...)
{
#
# This function takes an MDA object and shrinks
# it, creating a new mda object
#
# First check for a null condition
#
	if(is.null(sub.df) & is.null(tot.df)) {
		warning("No shrinkage parameters (sub.df or tot.df) given")
		return(object)
	}
#
# First recover theta 
	dimension <- object$dimension
	lambda <- object$values[seq(dimension)]
	theta.mod <- object$theta.mod
	theta <- object$means
	alpha <- sqrt(lambda)
	sqima <- sqrt(1 - lambda)	#
	theta <- scale(theta, F, alpha/sqima)	#
# Construct penalty
#
	sub.prior <- unlist(object$sub.prior)
	prior <- object$prior
	Rj <- sapply(object$assign.theta, length)
	dp <- sub.prior * rep(prior, Rj)
	cl <- rep(seq(Rj), Rj)
	P <- diag(dp) + transform.penalty(prior = dp, cl = cl, df = sub.df, 
		tot.df = tot.df)
	K <- t(theta) %*% P %*% theta
	TT <- chol((K + t(K))/2)
	Tinv <- backsolve(TT, diag(length(lambda)))	#
# lambda*Tinv = diag(lambda)%*% Tinv
#
	M <- t(Tinv) %*% (lambda * Tinv)
	ed <- svd(M)
	thetan <- ed$v
	lambda <- ed$d	
	# Note: the discriminant eigenvalues are a transformation
# of the optimal scaling values 
	discr.eigen <- lambda/(1 - lambda)
	pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
	dimension <- min(dimension, sum(lambda > .Machine$double.eps))
	if(dimension == 0) {
		warning("degenerate problem; no discrimination")
		return(structure(list(dimension = 0, fit = fit, call = 
			this.call), class = "fda"))
	}
	thetan <- thetan[, seq(dimension), drop = F]	#
	pe <- pe[seq(dimension)]	#Now produce projected centroids
#
	alpha <- sqrt(lambda[seq(dimension)])
	sqima <- sqrt(1 - lambda[seq(dimension)])	#
# package up results
#
	dm <- dimnames(theta)
	vnames <- dm[[2]][seq(dimension)]
	means <- scale(theta %*% Tinv %*% thetan, F, sqima/alpha)
	dimnames(means) <- list(dm[[1]], vnames)
	names(lambda) <- c(vnames, rep("", length(lambda) - dimension))	#
	names(pe) <- vnames
	theta.mod <- theta.mod %*% Tinv %*% thetan
	object$confusion <- object$deviance <- NULL
	incl.names <- c("percent.explained", "values", "means", "theta.mod", 
		"dimension")
	rl <- list(pe, lambda, means, theta.mod, dimension)
	names(rl) <- incl.names
	object$sub.df <- sub.df
	object$tot.df <- tot.df
	object[incl.names] <- rl
	object$weights <- NULL
	update(object, sub.df = sub.df, tot.df = tot.df, weights = object, ...)
}
"transform.penalty"<-
function(Q, prior, cl, df = NULL, tot.df = NULL)
{
	if(missing(Q))
		Q <- mean.penalty(prior, cl)
	if(missing(prior))
		prior <- attr(Q, "prior")
	if(missing(cl))
		cl <- attr(Q, "cl")
	transform.pen <- function(Q, prior, df)
	{
		df.inv <- function(d, df, lambda = NULL, iterations = 10)
		{
#
#this function solves for lambda such that sum(1/(1 + d*lambda)) = df
			if(is.null(lambda)) {
				lambda <- 0.10000000000000001
				while(sum(1/(1 + d * lambda)) >= df) lambda <- 
				    lambda * 2
			}
			df.diriv <- function(d, lambda)
 - sum((d * lambda)/(1 + d * lambda)^2)
			current.df <- sum(1/(1 + d * lambda))
			if(abs((df - current.df)/df) < 0.0001 | iterations == 1
				)
				return(lambda, df = current.df)
			else {
#cat(df, lambda, current.df,"\n")
				lambda <- exp(log(lambda) - (current.df - df)/
				  df.diriv(d, lambda))
				Recall(d, df, lambda, iterations - 1)
			}
		}
		pQp <- Q/outer(sqrt(prior), sqrt(prior))
		d <- svd(pQp)$d
		lambda <- df.inv(d, df)$lambda
		lambda * Q
	}
	if(!is.null(tot.df)) {
		if(tot.df >= length(prior))
			return(Q * 0)
		else return(transform.pen(Q, prior, tot.df))
	}
	else {
		ncl <- unique(cl)
		df <- rep(df, length = length(ncl))
		for(i in seq(along = ncl)) {
			which <- cl == ncl[i]
			Q[which, which] <- Recall(Q[which, which, drop = F], 
				prior[which], tot.df = df[i])
		}
		return(Q)
	}
}
