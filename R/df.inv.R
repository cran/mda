"df.inv" <-
function(d, df, lambda = 10^4, iterations = 10)
{
    df.diriv <- function(d, lambda)
        - sum((d * lambda)/(d + lambda)^2)
    current.df <- sum(d/(d + lambda))	
    ##	cat("lambda,df", lambda, current.df, "\n")
    if(abs((df - current.df)/df) < 0.0001 | iterations == 1)
        return(list(lambda = lambda, df = current.df))
    else {
        lambda <- exp(log(lambda) - (current.df - df)/df.diriv(d, lambda))
        Recall(d, df, lambda, iterations - 1)
    }
}
