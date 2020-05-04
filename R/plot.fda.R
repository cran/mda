plot.fda <-
  function (x, data, coords = c(1, 2), group = c("true",
                                            "predicted"), colors, pch, mcolors = colors, mpch,
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
    if(group=="predicted"){
       vars <- predict(object, data, type = "var")
       g <- predict(object, data)
     }
    else{
      ff <- terms(object)
      attr(ff, "intercept") <- 0
      m <- model.frame(ff, data)
      x <- model.matrix(ff, m)
      vars <- predict(object, x, type = "var")
      g <- model.extract(m, "response")
    }
  }
  means <- object$means
  if(ncol(means)==1)stop("Only one canonical variate; plot requires at least two")
  g <- as.factor(g)
  cc <- as.numeric(g)
  np <- seq(levels(g))
 mit.colors=c(Orange = "#FF9233", Cyan = "#29D0D0", Lt.Green = "#81C57A",
Dk.Gray = "#575757", Red = "#AD2323", Blue = "#2A4BD7", Green = "#1D6914",
Brown = "#814A19", Purple = "#8126C0", Lt.Gray = "#A0A0A0", Yellow = "#FFEE33",
Pink = "#FFCDF3")
 if (missing(colors)) colors <- mit.colors
 colors <- rep(colors, length = length(np))
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
         ..., type = "n", xlab = paste("Canonical Var",
                            coord.pair[1]), ylab = paste("Canonical Var",
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

