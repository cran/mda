set.seed(20101)

data(ESL.mixture)
mixture.train=data.frame(x=I(ESL.mixture$x),y=ESL.mixture$y)

fda.fit=fda(y~x,data=mixture.train)
print(coef(fda.fit))
mda.fit=mda(y~x,data=mixture.train)
print(coef(mda.fit))
mda.predict = predict(mda.fit,mixture.train$x)
print(mda.predict)
fitb=fda(y~x,data=mixture.train,method=bruto,cost=1)
fitm=fda(y~x,data=mixture.train,method=mars,cost=1)
fitb.confusion = confusion(fitb)
print(fitb.confusion)
fitm.confusion = confusion(fitm)
print(fitm.confusion)
mda.ppred=predict(mda.fit,newdata=ESL.mixture$xnew)
##print(mda.ppred)

objects  <- list(
    fda.fit = fda.fit,
    mda.fit = mda.fit,
    fitb = fitb,
    fitm = fitm,
    fitb.confusion = fitb.confusion,
    fitm.confusion = fitm.confusion,
    mda.ppred = mda.ppred)
##saveRDS(objects, "test_results/mda-0.4-results.RDS")

expected  <- readRDS("test_results/mda-0.4-results.RDS")
for (x in names(objects)) {
    cat(sprintf("Testing %s\n", x))
    if (x == "fitb") {
        expect_equal(objects[[x]], expected[[x]], tol = 1e-7)
    } else {
        expect_equal(objects[[x]], expected[[x]])
    }
}


