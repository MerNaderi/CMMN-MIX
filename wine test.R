library(semTools)
library(gclus)
library(mvShapiroTest)
data("wine")

S = subset(wine, wine[, 1] == 3)
mardiaSkew(S[, -1])
mardiaKurtosis(S[, -1])

mardiaSkew(wine[, -1])
mardiaKurtosis(wine[, -1])

mvShapiro.Test(as.matrix(wine[, -1]))
MVN :: mvn(as.matrix(S[, -1]))
