"bruto" <-
function (x, y, w = rep(1, n), wp = rep(1/np, np), dfmax, cost = 2, 
          maxit.select = 20, maxit.backfit = 20, thresh = 1e-04,
          trace.bruto = FALSE, start.linear = TRUE, fit.object, ...)
{
    this.call <- match.call()
    y <- as.matrix(y)
    x <- as.matrix(x)
    np <- ncol(y)
    d <- dim(x)
    nq <- d[2]
    n <- d[1]
    xnames <- dimnames(x)[[2]]
    if (!length(xnames)) 
        xnames <- NULL
    ynames <- dimnames(y)[[2]]
    if (!length(ynames)) 
        ynames <- NULL
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(w) <- "double"
    storage.mode(wp) <- "double"
    storage.mode(cost) <- "double"
    if (missing(fit.object)) {
        nknotl <- function(n) {
            a1 <- log(50)/log(2)
            a2 <- log(100)/log(2)
            a3 <- log(140)/log(2)
            a4 <- log(200)/log(2)
            cx <- as.vector(cut(n, c(0, 50, 200, 800, 3200)))
            if (is.na(cx)) 
                cx <- 5
            floor(switch(cx, n, 2^(a1 + ((a2 - a1) * (n - 50))/150), 
                2^(a2 + ((a3 - a2) * (n - 200))/600), 2^(a3 + 
                  ((a4 - a3) * (n - 800))/2400), 200 + (n - 3200)^0.2) + 
                6)
        }
        check.range <- apply(x, 2, var)
        if (any(check.range < .Machine$double.eps)) 
            stop(paste("A column of x is constant;",
                       "do not include an intercept column"))
        nkmax <- nknotl(n) - 4
        coef <- matrix(double(nkmax * np * nq), ncol = nq)
        ybar <- apply(y * w, 2, sum)/sum(w)
        if (start.linear && (nq * cost > n)) 
            start.linear <- FALSE
        if (start.linear) {
            start.fit <- polyreg(x, y, w)
            eta <- fitted(start.fit)
            coef[seq(from = 2, by = 2, length = np), ] <-
                t(start.fit$coef)[, -1]
            type <- as.integer(rep(2, nq))
            df <- as.double(rep(1, nq))
        }
        else {
            eta <- outer(rep(1, n), ybar)
            type <- integer(nq)
            df <- double(nq)
        }
        nk <- integer(nq)
        knot <- matrix(double((nkmax + 4) * nq), ncol = nq)
        Match <- matrix(integer(n * nq), ncol = nq)
        nef <- integer(nq)
        lambda <- double(nq)
        xrange <- matrix(double(2 * nq), 2, nq)
        df <- double(nq)
        if (missing(dfmax)) 
            dfmax <- (2 * nkmax)/3
        if (length(dfmax) != nq) 
            dfmax <- rep(dfmax, length = nq)
        if (cost > 0) {
            TD <- (n - sum(df))/cost
            TT <- dfmax > TD
            if (any(TT)) 
                dfmax[TT] <- TD
        }
        storage.mode(dfmax) <- "double"
    }
    else {
        this.call <- fit.object$call
        ybar <- fit.object$ybar
        nkmax <- fit.object$nkmax
        dfmax <- fit.object$dfmax
        eta <- fit.object$fitted.values
        if (is.null(eta)) 
            eta <- predict(fit.object, x)
        nk <- fit.object$nk
        knot <- fit.object$knot
        Match <- fit.object$Match
        nef <- fit.object$nef
        lambda <- fit.object$lambda
        coef <- fit.object$coef
        type <- unclass(fit.object$type)
        xrange <- fit.object$xrange
        maxit.select <- 0
        maxit.backfit <- fit.object$nit[2]
        df <- fit.object$df
    }
    maxit <- as.integer(c(maxit.select, maxit.backfit))
    names(df) <- xnames
    names(maxit) <- c("selection", "backfitting")
    gcv.select <- if (maxit.select) 
        matrix(double(maxit.select * nq), nq, maxit.select)
    else double(1)
    gcv.backfit <- if (maxit.backfit) 
        matrix(double(maxit.backfit * nq), nq, maxit.backfit)
    else double(1)
    df.select <- if (maxit.select) 
        matrix(double(maxit.select * nq), nq, maxit.select)
    else double(1)
    names(lambda) <- xnames
    fit <-
        .Fortran("bruto",
                 x,
                 as.integer(n),
                 as.integer(nq), 
                 y,
                 as.integer(np),
                 w,
                 knot = knot,
                 nkmax = as.integer(nkmax), 
                 nk = nk,
                 wp,
                 Match = Match,
                 nef = nef,
                 dfmax = dfmax, 
                 cost = cost,
                 lambda = lambda,
                 df = df,
                 coef = coef,
                 type = type, 
                 xrange = xrange,
                 gcv.select = gcv.select,
                 gcv.backfit = gcv.backfit, 
                 df.select = df.select,
                 maxit = maxit,
                 nit = maxit,
                 fitted.values = eta, 
                 residuals = y - eta,
                 as.double(thresh),
                 double((2 * np + 2) * ((n + 1) + 1) + (2 * np + 16) *
                        (n + 1) + 2 * (n + 1) + np),
                 integer(n),
                 trace.bruto,
                 PACKAGE = "mda")[c("knot", "nkmax", "nk", "Match",
                 "nef", "dfmax", "cost", "lambda", "df", "coef", "type",
                 "xrange", "gcv.select", "gcv.backfit", "df.select",
                 "maxit", "nit", "fitted.values", "residuals")]
    if (TN <- fit$nit[1]) {
        TT <- fit$gcv.select[, seq(TN), drop = FALSE]
        dimnames(TT) <- list(xnames, NULL)
    }
    else TT <- NULL
    fit$gcv.select <- TT
    if (TN) {
        TT <- fit$df.select[, seq(TN), drop = FALSE]
        dimnames(TT) <- list(xnames, NULL)
    }
    else TT <- NULL
    fit$df.select <- TT
    if (TN <- fit$nit[2]) {
        TT <- fit$gcv.backfit[, seq(TN), drop = FALSE]
        dimnames(TT) <- list(xnames, NULL)
    }
    else TT <- NULL
    fit$gcv.backfit <- TT
    TT <- factor(fit$type, levels = 1:3, labels = c("excluded", 
        "linear", "smooth"))
    names(TT) <- xnames
    fit$type <- TT
    fit$ybar <- ybar
    fit$call <- this.call
    structure(fit, class = "bruto")
}
"coef.fda" <-
function (object, type = c("canonical", "discriminant"), ...) 
{
    type <- match.arg(type)
    fit <- object$fit
    Coefs <- fit$coef
    if (is.null(Coefs)) 
        stop("No explicit coefficients in this formulation")
    dimension <- object$dimension
    Coefs <- Coefs %*% object$theta[, seq(dimension), drop = FALSE]
    lambda <- object$values
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    Coefs <- scale(Coefs, FALSE, sqima * alpha)
    if (type == "discriminant") 
        Coefs <- Coefs %*% t(object$means)
    Coefs
}
"confusion" <-
function (object, ...) 
UseMethod("confusion")
"confusion.default" <-
function (object, true, ...) 
{
    if (inherits(object, "data.frame")) 
        confusion.list(object, true)
    else {
        jt <- table(object, true)
        jd <- dimnames(jt)
        jn <- unlist(jd)
        ju <- jn[duplicated(jn)]
        j1 <- jd[[1]][!match(jd[[1]], ju, 0)]
        j2 <- jd[[2]][!match(jd[[2]], ju, 0)]
        jt <- jt[c(ju, j1), c(ju, j2), drop = FALSE]
        realjt <- jt[ju, ju, drop = FALSE]
        ntot <- sum(jt)
        mismatch <- (ntot - sum(realjt))/ntot
        structure(jt, error = (1 - sum(diag(realjt))/sum(realjt)), 
            mismatch = if (mismatch > 0) 
                mismatch
            else NULL)
    }
}
"confusion.fda" <-
function (object, data, ...) 
{
    if (missing(data)) 
        return(object$confusion)
    Terms <- terms(object)
    attr(Terms, "intercept") <- 0
    m <- model.frame(Terms, data)
    x <- model.matrix(Terms, m)
    g <- model.extract(m, response)
    confusion.default(predict(object, x, ...), g)
}
"confusion.list" <-
function (object, truth, ...)
{
    dd <- names(object)
    y <- seq(dd)
    x <- attr(object, "dimension")
    if (!length(x)) 
        x <- seq(dd)
    for (i in y) {
        confi <- confusion(object[, i], truth)
        y[i] <- attr(confi, "error")
    }
    return(x, y)
}
"contr.fda" <-
function (p = rep(1, d[1]), contrast.default = contr.helmert(length(p))) 
{
    d <- dim(contrast.default)
    sqp <- sqrt(p/sum(p))
    x <- cbind(1, contrast.default) * outer(sqp, rep(1, d[2] + 
        1))
    qx <- qr(x)
    J <- qx$rank
    qr.qy(qx, diag(d[1])[, seq(2, J)])/outer(sqp, rep(1, J - 
        1))
}
"fda" <-
function (formula = formula(data), data = sys.frame(sys.parent()), 
    weights, theta, dimension = J - 1, eps = .Machine$double.eps, 
    method = polyreg, keep.fitted = (n * dimension < 1000), ...) 
{
    this.call <- match.call()
    m <- match.call(expand = FALSE)
    m[[1]] <- as.name("model.frame")
    m <- m[match(names(m), c("", "formula", "data", "weights"), 
        0)]
    m <- eval(m, sys.frame(sys.parent()))
    Terms <- attr(m, "terms")
    g <- model.extract(m, response)
    attr(Terms, "intercept") <- 0
    x <- model.matrix(Terms, m)
    dd <- dim(x)
    n <- dd[1]
    weights <- model.extract(m, weights)
    if (!length(weights)) 
        weights <- rep(1, n)
    else if (any(weights < 0)) 
        stop("negative weights not allowed")
    if (length(g) != n) 
        stop("g should have length nrow(x)")
    fg <- factor(g)
    prior <- table(fg)
    prior <- prior/sum(prior)
    cnames <- levels(fg)
    g <- as.numeric(fg)
    J <- length(cnames)
    iswt <- FALSE
    if (missing(weights)) 
        dp <- table(g)/n
    else {
        weights <- (n * weights)/sum(weights)
        dp <- tapply(weights, g, sum)/n
        iswt <- TRUE
    }
    if (missing(theta)) 
        theta <- contr.helmert(J)
    theta <- contr.fda(dp, theta)
    Theta <- theta[g, , drop = FALSE]
    fit <- method(x, Theta, weights, ...)
    if (iswt) 
        Theta <- Theta * weights
    ssm <- t(Theta) %*% fitted(fit)/n
    ed <- svd(ssm, nu = 0)
    thetan <- ed$v
    lambda <- ed$d
    lambda[lambda > 1 - eps] <- 1 - eps
    discr.eigen <- lambda/(1 - lambda)
    pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
    dimension <- min(dimension, sum(lambda > eps))
    if (dimension == 0) {
        warning("degenerate problem; no discrimination")
        return(structure(list(dimension = 0, fit = fit, call = this.call), 
            class = "fda"))
    }
    thetan <- thetan[, seq(dimension), drop = FALSE]
    pe <- pe[seq(dimension)]
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    vnames <- paste("v", seq(dimension), sep = "")
    means <- scale(theta %*% thetan, FALSE, sqima/alpha)
    dimnames(means) <- list(cnames, vnames)
    names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
    names(pe) <- vnames
    obj <- structure(list(percent.explained = pe, values = lambda, 
        means = means, theta.mod = thetan, dimension = dimension, 
        prior = prior, fit = fit, call = this.call, terms = Terms), 
        class = "fda")
    obj$confusion <- confusion(predict(obj), fg)
    if (!keep.fitted) 
        obj$fit$fitted.values <- NULL
    obj
}
"fix.theta" <-
function (theta, Q) 
{
    M <- t(theta) %*% Q %*% theta
    eM <- eigen(M, sym = TRUE)
    scale(theta %*% eM$vectors, FALSE, sqrt(eM$values))
}
"gen.ridge" <-
function (x, y, weights, lambda = 1, omega, df, ...) 
{
    if (missing(df) && lambda <= .Machine$double.eps) 
        return(polyreg(x, y))
    d <- dim(x)
    mm <- apply(x, 2, mean)
    x <- scale(x, mm, FALSE)
    simple <- if (missing(omega)) TRUE else FALSE
    if (!simple) {
        if (!all(match(c("values", "vectors"), names(omega), 
            FALSE))) 
            stop("You must supply an  eigen-decomposed version of omega")
        vals <- pmax(.Machine$double.eps, sqrt(omega$values))
        basis <- scale(omega$vectors, FALSE, vals)
        x <- x %*% basis
    }
    svd.x <- svd(x)
    dd <- svd.x$d
    if (!missing(df)) {
        while (sum(dd^2/(dd^2 + lambda)) > df) lambda <- lambda * 
            10
        junk <- df.inv(dd^2, df, lambda)
        lambda <- junk$lambda
        df <- junk$df
    }
    else df <- sum(dd^2/(dd^2 + lambda))
    y <- (t(t(y) %*% svd.x$u) * dd)/(dd^2 + lambda)
    coef <- svd.x$v %*% y
    fitted <- x %*% coef
    if (!simple) 
        coef <- basis %*% coef
    structure(list(fitted.values = fitted, coefficients = coef, 
        df = df, lambda = lambda, xmeans = mm), class = "gen.ridge")
}
"kmeans.start" <-
function (x, g, subclasses) 
{
    cnames <- levels(g <- factor(g))
    J <- length(cnames)
    g <- as.numeric(g)
    weights <- as.list(cnames)
    names(weights) <- cnames
    subclasses <- rep(subclasses, length = length(cnames))
    R <- sum(subclasses)
    cl <- rep(seq(J), subclasses)
    cx <- x[seq(R), , drop = FALSE]
    for (j in seq(J)) {
        nc <- subclasses[j]
        which <- cl == j
        xx <- x[g == j, , drop = FALSE]
        if ((nc <= 1) || (nrow(xx) <= nc)) {
            cx[which, ] <- apply(xx, 2, mean)
            wmj <- matrix(1, sum(g == j), 1)
        }
        else {
            start <- xx[sample(1:nrow(xx), size = nc), ]
            TT <- kmeans(xx, start)
            cx[which, ] <- TT$centers
            wmj <- diag(nc)[TT$cluster, ]
        }
        dimnames(wmj) <-
            list(NULL, paste("s", seq(dim(wmj)[2]), sep = ""))
        weights[[j]] <- wmj
    }
    list(x = cx, cl = factor(cl, labels = cnames), weights = weights)
}
"laplacian" <-
function (size = 16, compose = FALSE) 
{
    gmat <- matrix(0, size, size)
    xx <- seq(size)
    for (v in xx) gmat[, v] <- sqrt(2/size) * cos(((xx - 0.5) * 
        pi * (v - 1))/size)
    gmat[, 1] <- gmat[, 1]/sqrt(2)
    lvec <- -(2 * size^2) * (1 - cos(((xx - 1) * pi)/size))
    gmat <- kronecker(gmat, gmat)
    lvec <- rep(lvec, rep(size, size)) + rep(lvec, size)
    if (compose) 
        gmat %*% (lvec^2 * t(gmat))
    else list(vectors = gmat, values = lvec^2)
}
"llmult" <-
function (p, g) 
{
    index <- cbind(seq(along = g), as.numeric(g))
    p <- p[index]
    -2 * sum(log(p[p > .Machine$double.eps]))
}
"lvq.start" <-
function (x, g, subclasses) 
{
    cnames <- levels(fg <- factor(g))
    J <- length(cnames)
    g <- as.numeric(g)
    weights <- as.list(cnames)
    names(weights) <- cnames
    subclasses <- rep(subclasses, length = length(cnames))
    size <- sum(subclasses)
    cb <- lvqinit(x, g, size = size)
    TT <- olvq1(x, g, codebk = cb)
    TT <- lvq3(x, g, codebk = TT)
    cl <- as.numeric(TT$cl)
    R <- length(cl)
    cx <- TT$x
    p <- ncol(cx)
    for (j in seq(J)) {
        which <- cl == j
        number <- sum(which)
        if (number == 0) {
            cx <- rbind(cx, apply(x[g == j, ], 2, mean))
            cl <- c(cl, j)
            wmj <- matrix(1, sum(g == j), 1)
            number <- 1
        }
        else if (number == 1) 
            wmj <- matrix(1, sum(g == j), 1)
        else {
            jcx <- cx[which, ]
            jcl <- seq(number)
            jcluster <- lvqtest(list(x = jcx, cl = jcl), x[g == 
                j, ])
            needed <- unique(jcluster)
            rcl <- rep(0, number)
            rcl[needed] <- j
            cl[which] <- rcl
            wmj <- diag(number)[jcluster, needed, drop = FALSE]
            number <- length(needed)
        }
        dimnames(wmj) <- list(NULL, paste("s", seq(number), sep = ""))
        weights[[j]] <- wmj
    }
    TT <- cl > 0
    list(x = cx[TT, , drop = FALSE], cl = factor(cl[TT], labels = cnames), 
        weights = weights)
}
"make.dumpdata.mda" <-
function () 
dump(c("bruto", "coef.fda", "confusion", "contr.fda", "fda", 
    "mars", "model.matrix.mars", "polybasis", "polyreg", "pplot", 
    "predict.bruto", "predict.fda", "predict.mars", "predict.polyreg", 
    "print.fda", "softmax", "make.dumpdata.mda", "confusion.default", 
    "confusion.fda", "confusion.list", "fix.theta", "gen.ridge", 
    "predict.gen.ridge", "laplacian", "kmeans.start", "llmult", 
    "lvq.start", "mda", "mda.fit", "mda.means", "mda.start", 
    "meanPenalty", "plot.fda", "pplot4", "predict.mda", "print.mda", 
    "shrink", "shrink.mda", "transformPenalty"), "dumpdata.mda")
"mars" <-
function (x, y, w = rep(1, nrow(x)), wp, degree = 1, nk = max(21, 
    2 * ncol(x) + 1), penalty = 2, thresh = 0.001, prune = TRUE, 
    trace.mars = FALSE, forward.step = TRUE, prevfit = NULL, ...) 
{
    this.call <- match.call()
    if ((nk%%2) != 1) 
        nk <- nk - 1
    x <- as.matrix(x)
    np <- dim(x)
    n <- np[1]
    p <- np[2]
    y <- as.matrix(y)
    nclass <- ncol(y)
    if (is.null(np)) {
        np <- c(length(x), 1)
        x <- as.matrix(x)
    }
    if (forward.step) {
        interms <- 1
        lenb <- nk
        bx <- matrix(rep(0, nrow(x) * nk), nrow = n)
        res <- matrix(rep(0, nrow(x) * ncol(y)), nrow = n)
        fullin <- rep(0, nk)
        cuts <- NULL
        factor <- NULL
    }
    else {
        bx <- model.matrix.mars(prevfit, x, full = TRUE)
        interms <- ncol(bx)
        lenb <- prevfit$lenb
        o <- prevfit$all.terms
        fullin <- rep(0, ncol(bx))
        fullin[o] <- 1
        res <- prevfit$res
        factor <- prevfit$factor
        cuts <- prevfit$cuts
        if (missing(penalty)) 
            penalty <- prevfit$penalty
        degree <- prevfit$degree
        nk <- lenb
        thresh <- prevfit$thresh
    }
    if (missing(penalty) & (degree > 1)) 
        penalty <- 3
    if (!missing(wp)) {
        if (any(wp <= 0)) 
            stop("wp should all be positive")
        wp <- sqrt(wp/sum(wp))
        y <- y * outer(rep(1, n), wp)
    }
    else wp <- NULL
    tagx <- x
    storage.mode(tagx) <- "integer"
    for (j in 1:p) {
        tagx[, j] <- order(x[, j])
    }
    bestin <- rep(0, nk)
    flag <- matrix(rep(0, nk * p), nrow = nk, ncol = p)
    if (is.null(cuts)) 
        cuts <- matrix(rep(0, nk * p), nrow = nk, ncol = p)
    if (is.null(factor)) {
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
    lenscrat <- 1 + n + 2 * n * nk + 4 * nk * nk + 3 * nk + 3 * 
        nk * nclass + 3 * nclass + 28 * n + 51
    junk <- .Fortran("marss", as.integer(n), as.integer(n), as.integer(p), 
        as.integer(nclass), as.matrix(y), as.matrix(x), as.double(w), 
        as.matrix(tagx), as.integer(degree), as.integer(nk), 
        as.double(penalty), as.double(thresh), as.logical(forward.step), 
        as.integer(interms), as.logical(prune), bx = as.matrix(bx), 
        fullin = as.integer(fullin), lenb = as.integer(lenb), 
        bestgcv = as.double(bestgcv), bestin = as.integer(bestin), 
        flag = as.matrix(flag), cuts = as.matrix(cuts), dir = as.matrix(dir), 
        res = as.matrix(res), alpha = as.double(alpha), beta = as.matrix(beta), 
        double(lenscrat), integer(4 * nk), trace.mars, PACKAGE = "mda")
    lenb <- junk$lenb
    all.terms <- seq(lenb)[junk$fullin[1:lenb] == 1]
    selected.terms <- seq(lenb)[junk$bestin[1:lenb] == 1]
    coefficients <- junk$beta[seq(selected.terms), , drop = FALSE]
    residuals <- junk$res
    fitted.values <- y - residuals
    if (!is.null(wp)) {
        TT <- outer(rep(1, n), wp)
        residuals <- residuals/TT
        fitted.values <- fitted.values/TT
        coefficients <- coefficients/outer(rep(1, length(selected.terms)), 
            wp)
    }
    dir <- junk$dir[seq(lenb), , drop = FALSE]
    dimnames(dir) <- list(NULL, dimnames(x)[[2]])
    cutss <- junk$cuts[seq(lenb), , drop = FALSE]
    x <- junk$bx[, selected.terms, drop = FALSE]
    structure(list(call = this.call, all.terms = all.terms, selected.terms = selected.terms, 
        penalty = penalty, degree = degree, nk = nk, thresh = thresh, 
        gcv = junk$bestgcv, factor = dir, cuts = cutss, residuals = residuals, 
        fitted.values = fitted.values, lenb = junk$lenb, coefficients = coefficients, 
        x = x), class = "mars")
}
"mda" <-
function(formula = formula(data), data = sys.frame(sys.parent()), 
         subclasses = 3, sub.df = NULL, tot.df = NULL, dimension =
         sum(subclasses) - 1, eps = .Machine$double.eps, iter = 5,
         weights = mda.start(x, g, subclasses, trace, ...), method =
         polyreg, keep.fitted = (n * dimension < 1000), trace = FALSE,
         ...) 
{
    this.call <- match.call()
    m <- match.call(expand = FALSE)
    m[[1]] <- as.name("model.frame")
    m <- m[match(names(m), c("", "formula", "data"), 0)]
    m <- eval(m, sys.frame(sys.parent()))
    Terms <- attr(m, "terms")
    g <- model.extract(m, response)
    attr(Terms, "intercept") <- 0
    x <- model.matrix(Terms, m)
    dd <- dim(x)
    n <- dd[1]
    m <- dd[2]
    rowmin <- function(mat) {
        ncc <- ncol(mat)
        if (ncc == 1) 
            return(drop(mat))
        rowm <- pmin(mat[, 1], mat[, 2])
        if (ncc == 2) 
            return(rowm)
        else {
            for (j in seq(from = 3, to = ncc)) rowm <- pmin(rowm, 
                mat[, j])
        }
        rowm
    }
    if (length(g) != n) 
        stop("g should have length nrow(x)")
    fg <- factor(g)
    if (inherits(weights, "mda")) {
        if (is.null(weights$weights)) 
            weights <- predict(weights, x, type = "weights", 
                g = fg)
        else weights <- weights$weights
    }
    subclasses <- sapply(weights, ncol)
    prior <- table(fg)
    dim(prior) <- NULL
    prior <- prior/sum(prior)
    cnames <- levels(fg)
    g <- as.numeric(fg)
    J <- length(cnames)
    Assign <- split(seq(sum(subclasses)), rep(seq(J), subclasses))
    names(Assign) <- cnames
    if (!is.null(tot.df)) {
        if (tot.df >= sum(subclasses)) 
            tot.df <- NULL
    }
    if (!is.null(sub.df)) {
        sub.df <- rep(sub.df, length = length(prior))
        sub.df <- pmin(sub.df, subclasses)
        if (all(sub.df == subclasses)) 
            sub.df <- NULL
    }
    for (counter in seq(iter)) {
        fit <- mda.fit(x, g, weights, assign.theta = Assign, 
            Rj = subclasses, sub.df = sub.df, tot.df = tot.df, 
            dimension = dimension, eps = .Machine$double.eps, 
            method = method, trace = trace, ...)
        dmat <- predict.fda(fit, type = "distance")
        sub.prior <- fit$sub.prior
        for (j in seq(J)) {
            TT <- dmat[g == j, Assign[[j]], drop = FALSE]
            TT <- exp(-0.5 * (TT - rowmin(TT)))
            TT <- TT * outer(rep(1, nrow(TT)), sub.prior[[j]])
            weights[[j]][] <- TT/drop(TT %*% rep(1, ncol(TT)))
        }
        pclass <- matrix(1, n, J)
        dmat <- exp(-0.5 * (dmat - rowmin(dmat)))
        for (j in seq(J)) {
            priorj <- sub.prior[[j]]
            ass <- Assign[[j]]
            TT <- dmat[, ass, drop = FALSE] * outer(rep(1, n), priorj)
            TTot <- drop(TT %*% rep(1, length(ass)))
            pclass[, j] <- prior[j] * TTot
        }
        pclass <- pclass/drop(pclass %*% rep(1, J))
        if (trace) 
            cat("Iteration", counter, "\tDeviance(multinomial)", 
                format(round(ll <- llmult(pclass, g), 5)), "\n")
    }
    if (!trace) 
        ll <- llmult(pclass, g)
    if (!keep.fitted) 
        fit$fit$fitted.values <- NULL
    dimnames(pclass) <- list(NULL, names(Assign))
    conf <- confusion(softmax(pclass), fg)
    fit <- c(fit, list(weights = weights, prior = prior, assign.theta = Assign, 
        deviance = ll, confusion = conf, terms = Terms))
    fit$call <- this.call
    fit$sub.df <- sub.df
    fit$tot.df <- tot.df
    class(fit) <- c("mda", "fda")
    fit
}
"mda.fit" <-
function (x, g, weights, theta, assign.theta, Rj, sub.df, tot.df, 
    dimension = R - 1, eps = .Machine$double.eps, method = polyreg, 
    ...) 
{
    this.call <- match.call()
    n <- nrow(x)
    cnames <- names(weights)
    J <- length(cnames)
    R <- sum(Rj)
    wtots <- lapply(weights, function(x) apply(x, 2, sum))
    sub.prior <- lapply(wtots, function(x) x/sum(x))
    dp <- unlist(wtots)
    subclass.names <- names(dp)
    dp <- dp/sum(dp)
    if (missing(theta)) 
        theta <- contr.helmert(R)
    theta <- contr.fda(dp, theta)
    if (!(is.null(sub.df) & is.null(tot.df))) {
        Q <- diag(dp) + transformPenalty(prior = dp, cl = rep(seq(J), 
            Rj), df = sub.df, tot.df = tot.df)
        theta <- fix.theta(theta, Q)
    }
    Theta <- matrix(0, n, R - 1)
    obs.weights <- double(n)
    for (j in seq(J)) {
        Theta[g == j, ] <- weights[[j]] %*% theta[assign.theta[[j]], 
            , drop = FALSE]
        obs.weights[g == j] <- weights[[j]] %*% rep(1, Rj[j])
    }
    fit <- method(x, Theta, obs.weights, ...)
    Theta <- Theta * obs.weights
    ssm <- t(Theta) %*% fitted(fit)/n
    ed <- svd(ssm, nu = 0)
    thetan <- ed$v
    lambda <- ed$d
    lambda[lambda > 1 - eps] <- 1 - eps
    discr.eigen <- lambda/(1 - lambda)
    pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
    dimension <- min(dimension, sum(lambda > eps))
    if (dimension == 0) {
        warning("degenerate problem; no discrimination")
        return(structure(list(dimension = 0, fit = fit, call = this.call), 
            class = "fda"))
    }
    thetan <- thetan[, seq(dimension), drop = FALSE]
    pe <- pe[seq(dimension)]
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    vnames <- paste("v", seq(dimension), sep = "")
    means <- scale(theta %*% thetan, FALSE, sqima/alpha)
    dimnames(means) <- list(subclass.names, vnames)
    names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
    names(pe) <- vnames
    list(percent.explained = pe, values = lambda, means = means, 
        theta.mod = thetan, dimension = dimension, sub.prior = sub.prior, 
        fit = fit, call = this.call)
}
"mda.means" <-
function (object, x, y) 
{
    weights <- means <- object$weights
    nn <- names(object$weights)
    for (i in nn) {
        xx <- x[y == i, ]
        ww <- weights[[i]]
        nc <- ncol(ww)
        xm <- matrix(0, ncol(x), nc)
        for (j in seq(nc)) {
            www <- ww[, j]
            www <- www/sum(www)
            xm[, j] <- apply(xx * www, 2, sum)
        }
        means[[i]] <- xm
    }
    means
}
"mda.start" <-
function (x, g, subclasses = 3, trace.mda.start = FALSE, start.method = c("kmeans", 
    "lvq"), tries = 5, criterion = c("misclassification", "deviance"), 
    ...) 
{
    if (!require(class, quietly = TRUE)) 
        stop("mda() requires package `class'")
    start.method <- match.arg(start.method)
    criterion <- match.arg(criterion)
    name.criterion <- switch(criterion, misclassification = "Misclassification Error", 
        deviance = "Deviance(multinomial)")
    starter <- get(paste(start.method, "start", sep = "."), mode = "function")
    fg <- factor(g)
    cnames <- levels(fg)
    prior <- table(fg)
    J <- length(cnames)
    n <- length(g)
    g <- as.numeric(fg)
    best.ll <- 1/.Machine$double.eps
    for (try in seq(tries)) {
        start <- starter(x, fg, subclasses)
        weights <- start$weights
        if (criterion == "misclassification") {
            pg <- lvqtest(start, x)
            ll <- attr(confusion(pg, fg), "error")
        }
        else {
            subclasses <- sapply(weights, ncol)
            Assign <- split(seq(sum(subclasses)), rep(seq(J), 
                subclasses))
            names(Assign) <- cnames
            fit <- mda.fit(x, g, weights, assign.theta = Assign, 
                Rj = subclasses, eps = .Machine$double.eps, method = polyreg, 
                ...)
            dmat <- exp(-0.5 * predict.fda(fit, type = "distance"))
            sub.prior <- fit$sub.prior
            pclass <- matrix(1, n, J)
            for (j in seq(J)) {
                priorj <- sub.prior[[j]]
                ass <- Assign[[j]]
                TT <- dmat[, ass, drop = FALSE]
                TT <- TT * outer(rep(1, n), priorj)
                TTot <- drop(TT %*% rep(1, length(ass)))
                wmj <- TT[g == j, , drop = FALSE]/TTot[g == j]
                pclass[, j] <- prior[j] * TTot
                dimnames(wmj) <- list(NULL, paste("s", seq(along = ass), 
                  sep = ""))
                weights[[j]] <- wmj
            }
            pclass <- pclass/drop(pclass %*% rep(1, J))
            ll <- llmult(pclass, g)
        }
        if (trace.mda.start) 
            cat(start.method, "start   \t", name.criterion, format(round(ll, 
                5)), "\n")
        if (ll < best.ll) {
            keep.weights <- weights
            best.ll <- ll
        }
    }
    structure(keep.weights, criterion = best.ll, name.criterion = name.criterion)
}
"meanPenalty" <-
function (prior, cl) 
{
    n <- length(prior)
    Q <- matrix(0, n, n)
    ll <- unique(cl)
    for (lll in ll) {
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
"model.matrix.mars" <-
function (object, x, which = object$selected.terms, full = FALSE, 
    ...) 
{
    if (missing(x)) 
        return(object$x)
    x <- as.matrix(x)
    dd <- dim(x)
    n <- dd[1]
    p <- dd[2]
    nterms <- length(which)
    dir <- object$factor
    cut <- object$cuts
    if (full) {
        bx <- matrix(0, nrow = n, ncol = object$lenb)
        bx[, 1] <- 1
    }
    else bx <- matrix(1, nrow = n, ncol = nterms)
    which <- which[-1]
    for (i in seq(along = which)) {
        j <- which[i]
        if (all(dir[j, ] == 0)) {
            stop("error in factor or which")
        }
        temp1 <- 1
        for (k in 1:p) {
            if (dir[j, k] != 0) {
                temp2 <- dir[j, k] * (x[, k] - cut[j, k])
                temp1 <- temp1 * temp2 * (temp2 > 0)
            }
        }
        if (full) 
            bx[, j] <- temp1
        else bx[, i + 1] <- temp1
    }
    bx
}
"plot.fda" <-
function (x, data, g, coords = c(1, 2), group = c("true", 
    "predicted"), colors, pch, mcolors = max(colors) + 1, mpch, 
    pcex = 0.5, mcex = 2.5, ...) 
{
    object <- x                         # generic/method
    group <- match.arg(group)
    if (missing(data)) {
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
        if (group == "true") 
            g <- model.extract(m, response)
        else g <- predict(object, x)
    }
    means <- object$means
    g <- as.factor(g)
    cc <- as.numeric(g)
    np <- seq(levels(g))
    if (missing(colors)) 
        colors <- np
    else colors <- rep(colors, length = length(np))
    if (missing(pch)) 
        pch <- paste(np)
    else pch <- rep(paste(pch), length = length(np))
    mcolors <- rep(mcolors, length = length(np))
    if (missing(mpch)) 
        mpch <- pch
    else mpch <- rep(paste(mpch), length = length(np))
    assign <- object$assign
    if (is.null(assign)) 
        assign <- split(seq(np), seq(np))
    if (!is.matrix(coords)) {
        coords <- matrix(coords, length(coords), length(coords))
        tt <- lower.tri(coords)
        coords <- cbind(t(coords)[tt], coords[tt])
    }
    for (ii in seq(nrow(coords))) {
        coord.pair <- coords[ii, ]
        plot(rbind(vars[, coord.pair], means[, coord.pair]), 
            ..., type = "n", xlab = paste("Discriminant Var", 
                coord.pair[1]), ylab = paste("Discriminant Var", 
                coord.pair[2]), main = paste("Discriminant Plot for", 
                group, "classes"))
        for (i in np) {
            which <- cc == i
            if (any(which)) 
                points(vars[which, coord.pair, drop = FALSE], col = colors[i], 
                  pch = pch[i], cex = pcex)
            points(means[assign[[i]], coord.pair, drop = FALSE], 
                col = mcolors[i], pch = 1, cex = mcex)
            points(means[assign[[i]], coord.pair, drop = FALSE], 
                col = mcolors[i], pch = mpch[i], cex = mcex/2)
        }
    }
    invisible()
}
"polybasis" <-
function (x, degree = 1, monomial = FALSE) 
{
    if (degree >= 3) 
        warning("This is not a smart polynomial routine. You may get numerical problems with a degree of 3 or more")
    x <- as.matrix(x)
    dn <- dimnames(x)
    dd <- dim(x)
    np <- dd[2]
    if (np == 1) 
        monomial <- TRUE
    if (degree > 1) {
        if (monomial) {
            px <- x
            cc <- sapply(split(paste(diag(np)), rep(seq(np), 
                rep(np, np))), paste, collapse = "")
            tx <- x
            for (i in 2:degree) {
                px <- px * tx
                x <- cbind(x, px)
                cc <- c(cc, sapply(split(paste(diag(np) * i), 
                  rep(seq(np), rep(np, np))), paste, collapse = ""))
            }
        }
        else {
            matarray <- array(x, c(dd, degree))
            for (i in 2:degree) matarray[, , i] <- x^i
            matarray <- aperm(matarray, c(1, 3, 2))
            x <- matarray[, , np]
            ad0 <- seq(degree)
            ad <- ad0
            ncol.mat0 <- degree
            ncol.x <- degree
            d0 <- paste(ad0)
            cc <- d0
            for (ii in seq(np - 1, 1)) {
                index0 <- rep(seq(ncol.mat0), ncol.x)
                index <- rep(seq(ncol.x), rep(ncol.mat0, ncol.x))
                newad <- ad0[index0] + ad[index]
                retain <- newad <= degree
                mat0 <- matarray[, , ii]
                if (any(retain)) 
                  newmat <- mat0[, index0[retain], drop = FALSE] * 
                    x[, index[retain], drop = FALSE]
                else newmat <- NULL
                ddn <- paste(d0[index0[retain]], cc[index[retain]], 
                  sep = "")
                zeros <- paste(rep(0, nchar(cc[1])), collapse = "")
                cc <- paste(0, cc, sep = "")
                d00 <- paste(d0, zeros, sep = "")
                x <- cbind(mat0, x, newmat)
                cc <- c(d00, cc, ddn)
                ad <- c(ad0, ad, newad[retain])
                ncol.x <- length(ad)
            }
        }
        if (!is.null(dn)) 
            dn[[2]] <- cc
        else dn <- list(NULL, cc)
        dimnames(x) <- dn
    }
    cbind(Intercept = 1, x)
}
"polyreg" <-
function (x, y, w, degree = 1, monomial = FALSE, ...) 
{
    x <- polybasis(x, degree, monomial)
    y <- as.matrix(y)                   # just making sure ...
    if (iswt <- !missing(w)) {
        if (any(w <= 0)) 
            stop("only positive weights")
        w <- sqrt(w)
        y <- y * w
        x <- x * w
    }
    qrx <- qr(x)
    coef <- as.matrix(qr.coef(qrx, y))
    fitted <- qr.fitted(qrx, y)
    if ((df <- qrx$rank) < ncol(x)) 
        coef[qrx$pivot, ] <- coef
    if (iswt) 
        fitted <- fitted/w
    structure(list(fitted.values = fitted, coefficients = coef, 
        degree = degree, monomial = monomial, df = df), class = "polyreg")
}
"pplot" <-
function (x, g, colors, pch, add = FALSE, type = "p", ...) 
{
    g <- as.factor(g)
    cc <- as.numeric(g)
    np <- seq(levels(g))
    if (missing(colors)) 
        colors <- np + 1
    else colors <- rep(colors, length = length(np))
    if (missing(pch)) 
        pch <- paste(np)
    else pch <- rep(pch, length = length(np))
    if (!add) 
        plot(x, type = "n", ...)
    for (i in unique(cc)) points(x[cc == i, , drop = FALSE], col = colors[i], 
        pch = pch[i], type = type)
}
"pplot4" <-
function (x, ...) 
{
    par(mfrow = c(3, 2))
    for (i in 1:3) {
        for (j in (i + 1):4) pplot(x[, c(i, j)], xlab = paste("var", 
            i), ylab = paste("var", j), ...)
    }
}
"predict.bruto" <-
function (object, x, type = c("fitted", "terms"), ...)
{
    if (missing(x)) {
        z <- fitted(object)
        if (is.null(z)) 
            stop("need to supply x")
        else return(z)
    }
    d <- as.integer(dim(x))
    type <- match.arg(type)
    nq <- d[2]
    n <- d[1]
    if (nq != length(object$df)) 
        stop(paste("x should have the same number of columns",
                   "as the df component of object"))
    ybar <- object$ybar
    np <- as.integer(length(ybar))
    eta <- matrix(double(n * np), n, np)
    Type <- as.numeric(object$type)
    storage.mode(Type) <- "integer"
    storage.mode(x) <- "double"
    if (type == "fitted") {
        .Fortran("pbruto", x, n, nq, ybar, np, object$knot, object$nkmax, 
            object$nk, object$coef, Type, object$xrange, eta = eta, 
            eta, PACKAGE = "mda")$eta
    }
    else {
        ob <- as.list(seq(nq))
        names(ob) <- dimnames(x)[[2]]
        knot <- object$knot
        nk <- object$nk
        xrange <- object$xrange
        coef <- object$coef
        fitm <- matrix(double(n * np), n, np)
        dimnames(fitm) <- list(dimnames(x)[[1]], names(ybar))
        for (i in seq(nq)) {
            if (Type[i] > 1) 
                fit <- .Fortran("psspl2", x[, i], n, np, knot[, 
                  i], nk[i], xrange[, i], coef[, i], coef[, i], 
                  fit = fitm, as.integer(0), Type[i], PACKAGE = "mda")$fit
            else fit <- fitm
            ob[[i]] <- list(x = x[, i], y = fit)
        }
        ob
    }
}
"predict.fda" <-
function (object, x, type = c("class", "variates", "posterior", 
    "hierarchical", "distances"), prior, dimension = J - 1, ...)
{
    dist <- function(x, mean, m = ncol(mean)) (scale(x, mean, 
        FALSE)^2) %*% rep(1, m)
    type <- match.arg(type)
    means <- object$means
    Jk <- dim(means)
    J <- Jk[1]
    k <- Jk[2]
    if (type == "hierarchical") {
        if (missing(dimension)) 
            dimension.set <- seq(k)
        else {
            dimension.set <- dimension[dimension <= k]
            if (!length(dimension.set)) 
                dimension.set <- k
            dimension <- max(dimension.set)
        }
    }
    else dimension <- min(max(dimension), k)
    if (missing(x)) 
        y <- predict(object$fit)
    else {
        if (inherits(x, "data.frame") || is.list(x)) {
            Terms <- delete.response(terms(object))
            attr(Terms, "intercept") <- 0
            x <- model.matrix(Terms, x)
        }
        y <- predict(object$fit, x)
    }
    y <- y %*% object$theta[, seq(dimension), drop = FALSE]
    lambda <- object$values
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    x <- scale(y, FALSE, sqima * alpha)
    if (missing(prior)) 
        prior <- object$prior
    else {
        if (any(prior < 0) | round(sum(prior), 5) != 1) 
            stop("innappropriate prior")
    }
    means <- means[, seq(dimension), drop = FALSE]
    switch(type, variates = return(x), class = {
        n <- nrow(x)
        prior <- 2 * log(prior)
        mindist <- dist(x, means[1, ], dimension) - prior[1]
        pclass <- rep(1, n)
        for (i in seq(2, J)) {
            ndist <- dist(x, means[i, ], dimension) - prior[i]
            l <- ndist < mindist
            pclass[l] <- i
            mindist[l] <- ndist[l]
        }
        ## 2001-10-27: Need to provide levels or else if we get an error
        ## if the predicted classes do no contain all possible classes.
        ## Reported by Greg Jefferis <jefferis@stanford.edu>, fix by
        ## Bj/orn-Helge Mevik <bjorn-helge.mevik@matforsk.no>.
        return(factor(pclass, levels = seq(J),
                      labels = dimnames(means)[[1]]))
    }, posterior = {
        pclass <- matrix(0, nrow(x), J)
        for (i in seq(J)) pclass[, i] <- exp(-0.5 * dist(x, means[i, 
            ], dimension)) * prior[i]
        dimnames(pclass) <- list(dimnames(x)[[1]], dimnames(means)[[1]])
        return(pclass/drop(pclass %*% rep(1, J)))
    }, hierarchical = {
        prior <- 2 * log(prior)
        Pclass <- vector("list", length(dimension.set))
        names(Pclass) <- paste("D", dimension.set, sep = "")
        for (ad in seq(along = dimension.set)) {
            d <- dimension.set[ad]
            dd <- seq(d)
            mindist <- dist(x[, dd, drop = FALSE], means[1, dd, drop = FALSE], 
                d) - prior[1]
            pclass <- rep(1, nrow(x))
            for (i in seq(2, J)) {
                ndist <- dist(x[, dd, drop = FALSE], means[i, dd, 
                  drop = FALSE], d) - prior[i]
                l <- ndist < mindist
                pclass[l] <- i
                mindist[l] <- ndist[l]
            }
            levels(pclass) <- dimnames(means)[[1]]
            Pclass[[ad]] <- pclass
        }
        rownames <- dimnames(x)[[1]]
        if (is.null(rownames)) 
            rownames <- paste(seq(nrow(x)))
        return(structure(Pclass, class = "data.frame", row.names = rownames, 
            dimensions = dimension.set))
    }, distances = {
        dclass <- matrix(0, nrow(x), J)
        for (i in seq(J)) dclass[, i] <- dist(x, means[i, ], 
            dimension)
        dimnames(dclass) <- list(dimnames(x)[[1]], dimnames(means)[[1]])
        return(dclass)
    })
}
"predict.gen.ridge" <-
function (object, x, ...) 
{
    if (missing(x)) 
        fitted(object)
    else scale(x, object$xmeans, FALSE) %*% object$coef
}
"predict.mars" <-
function (object, x, ...) 
{
    if (missing(x)) {
        z <- fitted(object)
        if (is.null(z)) 
            stop("need to supply x")
        else return(z)
    }
    model.matrix.mars(object, x) %*% object$coefficients
}
"predict.mda" <-
function (object, x, type = c("class", "variates", "posterior", 
    "hierarchical", "weights"), prior = NULL, dimension = R - 
    1, g, ...) 
{
    type <- match.arg(type)
    Rk <- dim(object$means)
    R <- Rk[1]
    k <- Rk[2]
    if (type == "hierarchical") {
        if (missing(dimension)) 
            dimension.set <- seq(k)
        else {
            dimension.set <- dimension[dimension <= k]
            if (!length(dimension.set)) 
                dimension.set <- k
            dimension <- max(dimension.set)
        }
        Pclass <- vector("list", length(dimension.set))
        names(Pclass) <- paste("D", dimension.set, sep = "")
        for (ad in seq(along = dimension.set)) {
            d <- dimension.set[ad]
            Pclass[[ad]] <- if (missing(x)) 
                Recall(object, prior = prior, dimension = d, 
                  ...)
            else Recall(object, x, prior = prior, dimension = d, 
                ...)
        }
        rownames <- names(Pclass[[1]])
        if (is.null(rownames)) 
            rownames <- paste(seq(along = Pclass[[1]]))
        return(structure(Pclass, class = "data.frame", row.names = rownames, 
            dimensions = dimension.set))
    }
    else dimension <- min(max(dimension), k)
    if (is.null(prior)) 
        prior <- object$prior
    else {
        if (any(prior < 0) | round(sum(prior), 5) != 1) 
            stop("innappropriate prior")
    }
    if (type == "variates") 
        return(NextMethod("predict"))
    rowmin <- function(mat) {
        ncc <- ncol(mat)
        if (ncc == 1) 
            return(drop(mat))
        rowm <- pmin(mat[, 1], mat[, 2])
        if (ncc == 2) 
            return(rowm)
        else {
            for (j in seq(from = 3, to = ncc)) rowm <- pmin(rowm, 
                mat[, j])
        }
        rowm
    }
    dmat <- if (missing(x)) 
        predict.fda(object, type = "distances", dimension = dimension, 
            ...)
    else predict.fda(object, x, type = "distances", dimension = dimension, 
        ...)
    Assign <- object$assign
    sub.prior <- object$sub.prior
    J <- length(Assign)
    if (type == "weights") {
        if (missing(x)) 
            return(object$weights)
        g <- as.numeric(g)
        weights <- Assign
        for (j in seq(J)) {
            TT <- dmat[g == j, Assign[[j]], drop = FALSE]
            TT <- exp(-0.5 * (TT - rowmin(TT)))
            TT <- TT * outer(rep(1, nrow(TT)), sub.prior[[j]])
            weights[[j]] <- TT/drop(TT %*% rep(1, ncol(TT)))
        }
        return(weights)
    }
    pclass <- matrix(1, nrow(dmat), J)
    dmat <- exp(-0.5 * (dmat - rowmin(dmat)))
    for (j in seq(J)) {
        TT <- dmat[, Assign[[j]], drop = FALSE]
        TT <- TT * outer(rep(1, nrow(TT)), sub.prior[[j]])
        pclass[, j] <- prior[j] * drop(TT %*% rep(1, ncol(TT)))
    }
    dimnames(pclass) <- list(NULL, names(Assign))
    switch(type, class = softmax(pclass), posterior = pclass/drop(pclass %*% 
        rep(1, J)))
}
"predict.polyreg" <-
function (object, x, ...) 
{
    if (missing(x)) {
        z <- fitted(object)
        if (is.null(z)) 
            stop("need to supply x")
        else return(z)
    }
    degree <- object$degree
    monomial <- object$monomial
    polybasis(x, degree, monomial) %*% object$coef
}
"print.fda" <-
function (x, ...) 
{
    if (!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    cat("\nDimension:", format(x$dimension), "\n")
    cat("\nPercent Between-Group Variance Explained:\n")
    print(round(x$percent, 2))
    error <- x$confusion
    df <- x$fit
    if (!is.null(df)) 
        df <- df$df
    if (!is.null(df)) {
        cat("\nDegrees of Freedom (per dimension):", format(sum(df)), 
            "\n")
    }
    if (!is.null(error)) {
        n <- as.integer(sum(error))
        error <- format(round(attr(error, "error"), 5))
        cat("\nTraining Misclassification Error:", error, "( N =", 
            n, ")\n")
    }
    invisible(x)
}
"print.mda" <-
function (x, ...) 
{
    NextMethod("print")
    if (!is.null(x$deviance)) 
        cat("\nDeviance:", format(round(x$deviance, 3)), "\n")
    invisible(x)
}
"shrink" <-
function (object, ...) 
UseMethod("shrink")
"shrink.mda" <-
function (object, sub.df = NULL, tot.df = NULL, ...) 
{
    if (is.null(sub.df) & is.null(tot.df)) {
        warning("No shrinkage parameters (sub.df or tot.df) given")
        return(object)
    }
    dimension <- object$dimension
    lambda <- object$values[seq(dimension)]
    theta.mod <- object$theta.mod
    theta <- object$means
    alpha <- sqrt(lambda)
    sqima <- sqrt(1 - lambda)
    theta <- scale(theta, FALSE, alpha/sqima)
    sub.prior <- unlist(object$sub.prior)
    prior <- object$prior
    Rj <- sapply(object$assign.theta, length)
    dp <- sub.prior * rep(prior, Rj)
    cl <- rep(seq(Rj), Rj)
    P <- diag(dp) + transformPenalty(prior = dp, cl = cl, df = sub.df, 
        tot.df = tot.df)
    K <- t(theta) %*% P %*% theta
    TT <- chol((K + t(K))/2)
    Tinv <- backsolve(TT, diag(length(lambda)))
    M <- t(Tinv) %*% (lambda * Tinv)
    ed <- svd(M)
    thetan <- ed$v
    lambda <- ed$d
    discr.eigen <- lambda/(1 - lambda)
    pe <- (100 * cumsum(discr.eigen))/sum(discr.eigen)
    dimension <- min(dimension, sum(lambda > .Machine$double.eps))
    if (dimension == 0) {
        warning("degenerate problem; no discrimination")
        return(structure(list(dimension = 0, fit = fit, call = this.call), 
            class = "fda"))
    }
    thetan <- thetan[, seq(dimension), drop = FALSE]
    pe <- pe[seq(dimension)]
    alpha <- sqrt(lambda[seq(dimension)])
    sqima <- sqrt(1 - lambda[seq(dimension)])
    dm <- dimnames(theta)
    vnames <- dm[[2]][seq(dimension)]
    means <- scale(theta %*% Tinv %*% thetan, FALSE, sqima/alpha)
    dimnames(means) <- list(dm[[1]], vnames)
    names(lambda) <- c(vnames, rep("", length(lambda) - dimension))
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
    update(object, sub.df = sub.df, tot.df = tot.df, weights = object, 
        ...)
}
"softmax" <-
function (x, gap = FALSE) 
{
    d <- dim(x)
    maxdist <- x[, 1]
    pclass <- rep(1, d[1])
    for (i in seq(2, d[2])) {
        l <- x[, i] > maxdist
        pclass[l] <- i
        maxdist[l] <- x[l, i]
    }
    dd <- dimnames(x)[[2]]
    if (gap) {
        x <- abs(maxdist - x)
        x[cbind(seq(d[1]), pclass)] <- drop(x %*% rep(1, d[2]))
        gaps <- do.call("pmin", data.frame(x))
    }
    pclass <- if (is.null(dd) || !length(dd)) 
        pclass
    else factor(pclass, levels = seq(d[2]), labels = dd)
    if (gap) 
        list(class = pclass, gaps = gaps)
    else pclass
}
"transformPenalty" <-
function (Q, prior, cl, df = NULL, tot.df = NULL) 
{
    if (missing(Q)) 
        Q <- meanPenalty(prior, cl)
    if (missing(prior)) 
        prior <- attr(Q, "prior")
    if (missing(cl)) 
        cl <- attr(Q, "cl")
    transform.pen <- function(Q, prior, df) {
        df.inv <- function(d, df, lambda = NULL, iterations = 10) {
            if (is.null(lambda)) {
                lambda <- 0.1
                while (sum(1/(1 + d * lambda)) >= df) lambda <- lambda * 
                  2
            }
            df.diriv <- function(d, lambda) -sum((d * lambda)/(1 + 
                d * lambda)^2)
            current.df <- sum(1/(1 + d * lambda))
            if (abs((df - current.df)/df) < 1e-04 | iterations == 
                1) 
                return(list(lambda = lambda, df = current.df))
            else {
                lambda <- exp(log(lambda) - (current.df - df)/df.diriv(d, 
                  lambda))
                Recall(d, df, lambda, iterations - 1)
            }
        }
        pQp <- Q/outer(sqrt(prior), sqrt(prior))
        d <- svd(pQp)$d
        lambda <- df.inv(d, df)$lambda
        lambda * Q
    }
    if (!is.null(tot.df)) {
        if (tot.df >= length(prior)) 
            return(Q * 0)
        else return(transform.pen(Q, prior, tot.df))
    }
    else {
        ncl <- unique(cl)
        df <- rep(df, length = length(ncl))
        for (i in seq(along = ncl)) {
            which <- cl == ncl[i]
            Q[which, which] <- Recall(Q[which, which, drop = FALSE], 
                prior[which], tot.df = df[i])
        }
        return(Q)
    }
}
