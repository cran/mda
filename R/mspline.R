mspline <-
function(x, y, w = rep(1, n), df=5, lambda,
           thresh = 1e-04, ...)
{
    this.call <- match.call()
    y <- as.matrix(y)
    np <- ncol(y)
    n <- length(x)
    wp=rep(1/np,np)
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(w) <- "double"
    storage.mode(wp) <- "double"

     nknotl <- function(n) {
            a1 <- log(50)/log(2)
            a2 <- log(100)/log(2)
            a3 <- log(140)/log(2)
            a4 <- log(200)/log(2)
            cx <- as.numeric(cut(n, c(0, 50, 200, 800, 3200)))
            if (is.na(cx)) 
                cx <- 5
            floor(switch(cx, n, 2^(a1 + ((a2 - a1) * (n - 50))/150), 
                2^(a2 + ((a3 - a2) * (n - 200))/600), 2^(a3 + 
                  ((a4 - a3) * (n - 800))/2400), 200 + (n - 3200)^0.2) + 
                6)
        }
    nkmax <- nknotl(n) - 4
    coef <- double(nkmax * np)
    smat <- matrix(0,n,np)
    storage.mode(smat) <- "double"
    nk <- as.integer(0)
    knot <- double(nkmax + 4)
    Match <- integer(n)
    nef <- as.integer(0)
    if(df==0){
      method=1
      if(missing(lambda))stop("need a lambda")
         }
      else {
        method =2
        lambda=0
      }
    storage.mode(method) <- "integer"
    storage.mode(df) <- "double"
    storage.mode(lambda) <- "double"
    xrange <- double(2)
     fit <-
        .Fortran("sspl0",
                 x,
                 y,
                 w,
                 as.integer(n),
                 as.integer(np),
                 knot = knot,
                 nk=nk,
                 method,
                 tol=as.double(thresh),
                 wp,
                 Match = Match,
                 nef = nef,
                 center=as.integer(0),
                 dfoff=as.double(0),
                 dfmax = as.double(df+2), 
                 cost = as.double(0),
                 lambda = lambda,
                 df = df,
                 cv=double(1),
                 gcv=double(1),
                 coef = coef,
                 s=smat,
                 lev=double(n),
                 xrange = xrange,
                 work=double((2 * np + 2) * ((n + 1) + 1) + (2 * np + 16) *(n + 1) + 2 * (n + 1) + np),
                 iwork=integer(n),
                 integer(1),
                 PACKAGE = "mda")[c("knot",  "nk", "Match",
                 "nef", "lambda", "df", "coef","s",
                 "xrange","lev")]
    fit$call <- this.call
    structure(fit, class = "mspline")
}
