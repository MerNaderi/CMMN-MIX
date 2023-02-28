rm(list = ls(all.names = TRUE)) 

###--------- Data generation ---------###
CMMNE.MIX.gen = function(n, mu, sigma, lambda, nu, gamma, pi, MMNE = FALSE){
  
  #Mu and lambda are G * p
  #sigma is p * p * G array
  Class = numeric(n)
  Point.Type = numeric(n)
  p = dim(mu)[2]
  x = matrix(NA, nrow = n, ncol = p)
  
  for(i in 1:n){
    Z = rmultinom(1, 1, pi)
    j = Class[i] = which.max(Z)
    b = rbinom(1, 1, gamma[j])
    u = rexp(1, 1)
    
    if(b == 1) {
      Point.Type[i] = "Good"
      z = mvtnorm :: rmvnorm(1, mean = rep(0, p), sigma = sigma[, , j])
      x[i, ] = mu[j, ] + u * lambda[j, ] + z
    }
    
    if(b == 0) {
      Point.Type[i] = "bad"
      z = mvtnorm :: rmvnorm(1, rep(0, p), sigma = nu[j] * sigma[, , j])
      x[i, ] = mu[j, ] + u * sqrt(nu[j]) * lambda[j, ] + z
    }
    
    if(MMNE){
      z = mvtnorm :: rmvnorm(1, mean = rep(0, p), sigma = sigma[, , j])
      x[i, ] = mu[j, ] + u * lambda[j, ] + z
    }
  }
  
  if(MMNE) out = list(X = x, Class = Class)
  else out = list(X = x, Class = Class, Type = Point.Type)
  
  return(out)
}


PATH = paste("...")
source(paste(PATH, 'Initioal function.r', sep=""))
source(paste(PATH, 'Mix-CMMNE.r', sep=""))
source(paste(PATH, 'Mix-MMNE.r', sep=""))
#source(paste(PATH, '/Functions/Mix-MMNE-optim.r', sep=""))

d = 2
g = 3
sigma = array(0, dim = c(d, d, g))
sigma[,,1] = matrix( c(1, -0.18, -0.18, 1), ncol=d, byrow = T)
sigma[,,2] = matrix( c(2.33, 0, 0, 0.88), ncol=d, byrow = T)
sigma[,,3] = matrix( c(1.6 ,  0.58,  0.58, 0.9 ), ncol=d, byrow = T)

mu = matrix(0, nrow = g, ncol = d)
mu[1,] = matrix(c(-3,3) , byrow = T, nrow = 2)
mu[2,] = matrix(c(0,-1) , byrow = T, nrow = 2)
mu[3,] = matrix(c(5,2)  , byrow = T, nrow = 2)

lambda = matrix( c(1, -1, 1, 4, -3, 3), nrow = g)
nu = c(2, 5, 8)
PI = c(0.25, 0.35, 0.4)
gamma = 1 - c(0.15, 0.2, 0.1)


n = 1000
precent = 0.02

data = CMMNE.MIX.gen(n, mu, sigma, lambda, nu, gamma, PI, MMNE = T)
Y = data$X
plot(Y[,1], Y[,2], col = data$Class)

nI = c(rmultinom(1, n * precent, PI))
Y = rbind(Y, cbind(runif(nI[1], max(Y[data$Class == 1, 1]), max(Y[data$Class == 1, 1]) + 3), 
                   runif(nI[1], max(Y[data$Class == 1, 2]), max(Y[data$Class == 1, 2]) + 5)))
Y = rbind(Y, cbind(runif(nI[3], 5, 10), runif(nI[3], 15, 20)))
Y = rbind(Y, cbind(runif(nI[2], min(Y[data$Class == 2, 1]), -5), 
                   runif(nI[2], min(Y[data$Class == 2, 2]) - 5, min(Y[data$Class == 2, 1]))))

index = c(rep(1, n), rep(2, n * precent))
plot(Y[,1], Y[,2], col = c(data$Class, rep(1, nI[1]), rep(3, nI[3]), rep(2, nI[2])))


Type = c(rep(0, n), rep(1, nI[1]), rep(3, nI[3]), rep(2, nI[2]))


PATH = paste("...")
source(paste(PATH, 'Initioal function.r', sep=""))
source(paste(PATH, 'Mix-CMMNE.r', sep=""))
source(paste(PATH, 'Mix-MMNE.r', sep=""))

para  = rgpar(Y, g, method = 'inc-kmed')

para[[1]]$mu = mu[1,]
para[[2]]$mu = mu[2,]
para[[3]]$mu = mu[3,]

para[[1]]$delta = lambda[1, ]
para[[1]]$delta = lambda[2, ]
para[[1]]$delta = lambda[3, ]

para[[1]]$sigma = sigma[,,1]
para[[2]]$sigma = sigma[,,2]
para[[3]]$sigma = sigma[,,3]

para$pi = PI

CMMNE = MIX.CMMNE.EM(Y, para = para, tol = 1e-5,
                     type = "VVV", Stp.rule = "Atiken",
                     max.in.iter = 200, max.iter = 2000, per = 1,
                     print = F)
MMNE = MIX.MMNE.EM(Y, para = para, tol = 1e-5, fix.nu = T,
                   type = "VVV", Stp.rule = "Atiken",
                   max.in.iter = 200, max.iter = 2000, per = 100,
                   print = F)


Class = data$Class
km.clus = CMMNE$post.cate[1:n]
Clustering.CE = c(rand.index(Class, km.clus)[2], error.rate(km.clus, Class))

km.clus = MMNE$post.cate[1:n]
Clustering.E = c(rand.index(Class, km.clus)[2], error.rate(km.clus, Class))


Clustering.CE
Clustering.E



d.MMNE = function(Y,
                  mu,
                  Sigma,
                  lambda,
                  tau = NULL,
                  log = F) {
  if(is.null(tau))
    tau = 1
  inv.sig = matrixcalc::matrix.inverse(Sigma)
  xmu = sweep(Y, 2, mu)
  A = (c(lambda %*% inv.sig %*% t(xmu)) - tau) / sqrt(mahalanobis(lambda, 0, Sigma))
  
  PDF = mvtnorm::dmvnorm(Y,
                         mean = mu,
                         sigma = Sigma,
                         log = T) + pnorm(A, log.p = T) +
    log(tau * sqrt(2 * pi)) + 0.5 * A ^ 2 - log(sqrt(mahalanobis(lambda, 0, Sigma)))
  ifelse(log == T, return(PDF), return(exp(PDF)))
}

d.CMMNE = function(Y,
                   para = NULL,
                   tau = NULL,
                   nu = NULL,
                   gam = NULL,
                   log = F) {
  mu = para$mu
  Sigma = para$sigma
  lambda = para$delta
  
  Part1 = d.MMNE(Y, mu, nu * Sigma, sqrt(nu) * lambda, tau, log = T)
  Part2 = d.MMNE(Y, mu, Sigma, lambda, tau, log = T)
  Out = Part2 + log(gam + (1 - gam) * exp(Part1 - Part2))
  
  ifelse(log == T, return(Out), return(exp(Out)))
}

Mix.CMMNE <- function(Y, para = NULL) {
  fx = matrix(0, nrow = nrow(Y), ncol = length(para$pi))
  
  for (k in 1:length(para$pi)) {
    fx[, k] = d.CMMNE(
      Y,
      para = para[[k]],
      tau = para$tau[k],
      nu = para$nu[k],
      gam = para$gam[k],
      log = F
    )
  }
  
  val = apply(fx, 1, weighted.sum, wt = para$pi)
  val[val == 0] <- .Machine$double.xmin
  return(val)
}


d.MMNE1 = function(Y, para, nu = NULL, log = F) {
  p = ncol(Y)
  Y = as.matrix(Y, ncol = p)
  mu = para$mu
  Sigma = para$sigma
  lambda = para$delta
  if (!is.numeric(nu))
    nu = 1
  inv.sig = matrixcalc::matrix.inverse(Sigma)
  xmu = sweep(Y, 2, mu)
  A = (c(lambda %*% inv.sig %*% t(xmu)) - nu) / sqrt(mahalanobis(lambda, 0, Sigma))
  
  PDF = mvtnorm::dmvnorm(Y,
                         mean = mu,
                         sigma = Sigma,
                         log = T) + pnorm(A, log.p = T) +
    log(nu * sqrt(2 * pi)) + 0.5 * A ^ 2 - log(sqrt(mahalanobis(lambda, 0, Sigma)))
  ifelse(log == T, return(PDF), return(exp(PDF)))
}

Mix.MMNE <- function(Y, para = NULL, nu = NULL) {
  fx = matrix(0, nrow = nrow(Y), ncol = length(para$pi))
  if (!is.numeric(nu))
    nu = rep(1, length(para$pi))
  for (k in 1:length(para$pi))
    fx[, k] = d.MMNE1(Y, para[[k]], nu = nu[k])
  val = apply(fx, 1, weighted.sum, wt = para$pi)
  return(val)
}


para = CMMNE$para
X1 = seq(min(Y[, 1]), max(Y[, 1]), length = 100)
X2 = seq(min(Y[, 2]), max(Y[, 2]), length = 100)

F.out.CE = matrix(0, length(X1), length(X2))
for(i in 1:length(X1))
  for(j in 1:length(X2))
    F.out.CE[i, j] = Mix.CMMNE(cbind(X1[i], X2[j]), para)

para = MMNE$para

F.out = matrix(0, length(X1), length(X2))
for(i in 1:length(X1))
  for(j in 1:length(X2))
    F.out[i, j] = Mix.MMNE(cbind(X1[i], X2[j]), para, para$nu)



postscript(paste(PATH, '/fig1-sim1.eps', sep=''), width=15, height=35)
par(mfrow=c(1, 3), mar=c(4,4,3,0.5))

plot(Y[,1], Y[,2], ylab = expression(bold(X[2])),
     xlab = expression(bold(X[1])), main = "", pch = 21 + Type, 
     bg = c("deeppink3", "lemonchiffon3", "deepskyblue3")[unclass(
       c(data$Class, rep(1, nI[1]), rep(3, nI[3]), rep(2, nI[2])))], 
     cex.axis = 1.7, cex = 1.7, cex.main = 1, cex.lab= 1, font.lab=3)


plot(Y[,1], Y[,2], ylab = expression(bold(X[2])), #col = CMMNE$post.cate,
     xlab = expression(bold(X[1])), main = "CMMNE with ARI=0.99", pch = 21 + Type, 
     bg = c("deeppink3", "lemonchiffon3", "deepskyblue3")[unclass(CMMNE$post.cate)], 
     cex.axis = 1.7, cex = 1.7, cex.main = 1.5, cex.lab=1, font.lab=3)
contour(X1, X2, F.out.CE, add = T, nlevels = 10, lwd = 1.6, lty = 1)


Coll = MMNE$post.cate
Coll[MMNE$post.cate == 2] = 5
Col = ifelse(Coll == 3, 2, ifelse(Coll == 5, 3, 1))

plot(Y[,1], Y[,2], ylab = expression(bold(X[2])), #col = MMNE$post.cate,
     xlab = expression(bold(X[1])), main = "CMMNE with ARI=0.96", pch = 21 + Type, 
     bg = c("deeppink3", "lemonchiffon3", "deepskyblue3")[unclass(Col)], 
     cex.axis = 1.7, cex = 1.7, cex.main = 1.5, cex.lab=1, font.lab=3)
contour(X1, X2, F.out, add = T, nlevels = 10, lwd = 1.6, lty = 1)

dev.off()


