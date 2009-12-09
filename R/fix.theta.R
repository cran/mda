fix.theta <-
function (theta, Q) 
{
    M <- t(theta) %*% Q %*% theta
    eM <- eigen(M, sym = TRUE)
    scale(theta %*% eM$vectors, FALSE, sqrt(eM$values))
}

