

nSigma.par <- function(Sig.Struc = NULL,
                       p = NULL,
                       G = NULL) {
  if (is.null(p))
    stop("p is null")
  if (is.null(G))
    stop("G is null")
  if (is.null(Sig.Struc))
    stop("Sigma structure is null")
  
  if (Sig.Struc == "EII")
    npar = 1
  else if (Sig.Struc == "VII")
    npar = G
  else if (Sig.Struc == "EEI")
    npar = p
  else if (Sig.Struc == "VEI")
    npar = p + G - 1
  else if (Sig.Struc == "EVI")
    npar = p * G - G + 1
  else if (Sig.Struc == "VVI")
    npar = p * G
  else if (Sig.Struc == "EEE")
    npar = p * (p + 1) / 2
  else if (Sig.Struc == "EEV")
    npar = G * p * (p + 1) / 2 - (G - 1) * p
  else if (Sig.Struc == "VEV")
    npar = G * p * (p + 1) / 2 - (G - 1) * (p - 1)
  else if (Sig.Struc == "VVV")
    npar = G * p * (p + 1) / 2
  else if (Sig.Struc == "EVE")
    npar = p * (p + 1) / 2 + (G - 1) * (p - 1)
  else if (Sig.Struc == "VVE")
    npar = p * (p + 1) / 2 + (G - 1) * p
  else if (Sig.Struc == "VEE")
    npar = p * (p + 1) / 2 + (G - 1)
  else if (Sig.Struc == "EVV")
    npar = G * p * (p + 1) / 2 - (G - 1)
  else
    stop("Sigma structure is not correctly defined")
  
  return(npar)
}


Stop.rule <- function(log.lik) {
  if (length(log.lik) >= 3) {
    n = length(log.lik)
    l.new = log.lik[n]
    l.old  = log.lik[(n - 1)]
    l.old2  = log.lik[(n - 2)]
    ait = (l.new - l.old) / (l.old - l.old2)
    ln.Inf = l.old + (l.new - l.old) / (1 - ait)
    out = ln.Inf - l.new
    if (!is.na(out))
      out = abs(out)
    else
      out = 0
    
  } else
    out = (log.lik[2] - log.lik[1]) / abs(log.lik[1])
  return(out)
}


rgpar <- function(data, g = NULL, clus = NULL, 
                  method = c("k-means", "fast-kmed", "inc-kmed", 
                             "rank-kmed", "Simple-kmed", "trim-kmeans")) {
  w = matrix(1/g, nrow = nrow(data), ncol = g)
  if(is.null(clus)){
    if(method == "k-means") clus = kmeans(data, g)$clus
    
    ###  k-medoid = kmed
    if(method == "fast-kmed"){
      num <- as.matrix(data)
      mrwdist <- kmed :: distNumeric(num, num, method = "mrw")
      result <- kmed :: fastkmed(mrwdist, ncluster = g, iterate = 50)
      clus = result$cluster
    }
    if(method == "inc-kmed"){
      num <- as.matrix(data)
      mrwdist <- kmed :: distNumeric(num, num, method = "mrw")
      result <- kmed :: inckmed(mrwdist, ncluster = g, iterate = 50, alpha = 1.5)
      clus = result$cluster
    }
    if(method == "rank-kmed"){
      num <- as.matrix(data)
      mrwdist <- kmed :: distNumeric(num, num, method = "mrw")
      result <- kmed :: rankkmed(mrwdist, ncluster = g, iterate = 50)
      clus = result$cluster
    }
    if(method == "Simple-kmed"){
      num <- as.matrix(data)
      mrwdist <- kmed :: distNumeric(num, num, method = "mrw")
      result <- kmed :: skm(mrwdist, ncluster = g, seeding = 50)
      clus = result$cluster
    }
    if(method == "trim-kmeans"){
      result <- trimcluster :: trimkmeans(data, g - 1)
      clus = result$classification
    }
  } 
  
  #TT = table(clus)
  #O = order(TT, decreasing = T)
  val = list()
  for (k in 1:g){
    #val[[k]] = rmpar(data = data[clus == O[k],], g, wt = w[clus == O[k], O[k]])
    val[[k]] = rmpar(data = data[clus == k,], g, wt = w[clus == k, k])
  }
  #for (k in 1:g) val$pi[k] = sum(clus == O[k])/dim(data)[1]
  val$pi = table(clus)/dim(data)[1]
  val$nu = rep(1.5, g)
  val$gam = rep(0.9, g)
  val$tau = rep(1, g)
  return(val)
}


rmpar <- function(data , g , wt) {
  library(moments)
  par = list()
  p = ncol(data)
  par$mu = apply(data, 2, mean)
  par$delta = apply(data, 2, skewness)
  if (p == 1) {
    par$sigma = var(data)
  }
  if (p > 1) {
    sigma = ((cov.wt(
      data, wt = wt, method = "ML"
    )$cov))
    for (i in 1:p) {
      if (sigma[i, i] < 0.1) {
        sigma[i, i] = 0.1
      }
    }
    if (any(eigen(sigma)$values <= 0)) {
      par$sigma =  diag(apply(data, 2, var))
    }
    else{
      par$sigma = sigma
    }
  }
  return(par)
}

inv.mat = function(M) {
  eg = eigen(M)
  val = diag(eg$val)
  vec = cbind(eg$vec)
  vec %*% solve(val) %*% solve(vec)
}

weighted.sum <- function(z, wt, ...)
  return(sum(z * wt, ...))


S = function(b, u1, u2, nu) {
  u1 * b + nu * u2 * (1 - b)
}


newD3.MM <-
  function(D = NULL,
           d = NULL,
           G = NULL,
           Wk = NULL,
           Ak = NULL,
           tmax = 100) {
    z = matrix(0, d, d)
    lambda = 0
    for (g in 1:G) {
      lambdak = max(eigen(Wk[, , g])$values)
      z = z + diag(1 / Ak[, g]) %*% t(D) %*% Wk[, , g]  -
        lambdak * (diag(1 /
                          Ak[, g]) %*% t(D))
    }
    z1 = svd(z)
    Xk1 = (z1$v) %*% t(z1$u)
    return(Xk1)
  }

newD4.MM <-
  function(D = NULL,
           d = NULL,
           G = NULL,
           Wk = NULL,
           Ak = NULL,
           tmax = 100) {
    z = matrix(0, d, d)
    lambda = 0
    for (g in 1:G) {
      lambdak = max(1 / Ak[, g])
      #z = z + diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]  -lambdak *(diag(1/Ak[,g])%*% t(D) )
      z = z + Wk[, , g] %*% (D) %*% diag(1 / Ak[, g])  -
        lambdak * (Wk[, , g] %*% (D))
    }
    z1 = svd(z)
    #	Xk1 = (z1$v) %*% t(z1$u)
    #	return( t(Xk1) )
    # OR
    Xk1 = (z1$v) %*% t(z1$u)
    return(t(Xk1))
    
  }

newD <-
  function(D = NULL,
           d = NULL,
           G = NULL,
           Wk = NULL,
           Ak = NULL,
           tmax = 100) {
    D6 = D
    D6 = newD3.MM(
      D = D6,
      d = d,
      G = G,
      Wk = Wk,
      Ak = Ak,
      tmax = 100
    )
    D6 = newD4.MM(
      D = D6,
      d = d,
      G = G,
      Wk = Wk,
      Ak = Ak,
      tmax = 100
    )
    return(D6)
  }


tr = function(A)
  matrixcalc::matrix.trace(A)


Expect.TN = function(mu,
                     sigma2,
                     lower = 0,
                     upper = Inf) {
  alpha1 = (lower - mu) / sqrt(sigma2)
  alpha2 = (upper - mu) / sqrt(sigma2)
  
  if (lower == 0 && upper == Inf) {
    HF = exp(dnorm(alpha1, log = T) - pnorm(alpha1, lower.tail = F, log.p = T))
    E.x = mu + sqrt(sigma2) * HF
    E.x2 = mu ^ 2 + sigma2 + sigma2 * alpha1 * HF + 2 * mu * sqrt(sigma2) * HF
  }
  if (lower == -Inf && upper == 0) {
    HF = exp(dnorm(alpha2, log = T) - pnorm(alpha2, log.p = T))
    E.x = mu - sqrt(sigma2) * HF
    E.x2 = mu ^ 2 + sigma2 - sigma2 * alpha2 * HF - 2 * mu * sqrt(sigma2) * HF
  }
  
  return(E.x = E.x, E.x2 = E.x2)
}


M.step = function(Y,
                  para,
                  zhat,
                  Bhat,
                  what,
                  w2hat,
                  wtild,
                  w2tild,
                  type = NULL,
                  g,
                  max.in.iter) {
  
  ni = colSums(zhat)
  n = dim(Y)[1]
  d = dim(Y)[2]
  
  for (i in 1:g) {
    nu = para$nu[i]
    S1 = S(Bhat[, i], 1, 1, 1 / nu)
    S2 = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu))
    Iner = colSums(zhat[, i] * S1 * Y) - para[[i]]$delta * sum(zhat[, i] * S2)
    para[[i]]$mu = Iner / sum(S1 * zhat[, i])
    
    S1 = S(Bhat[, i], w2hat[, i], w2tild[, i], 1)
    eij = sweep(Y, 2, para[[i]]$mu)
    para[[i]]$delta = colSums(zhat[, i] * S2 * eij) / sum(S1 * zhat[, i])
  }
  
  if (type == "EII") {
    ro = matrix(c(numeric(d)), d, d)
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = ro + t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
    }
    Omeg.hat = tr(ro) / n / d
    for (i in 1:g) {
      para[[i]]$sigma = Omeg.hat * diag(d)
    }
  }
  
  if (type == "VII") {
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      para[[i]]$sigma = diag(d) * tr(ro) / ni[i] / d
    }
  }
  if (type == "EEI") {
    ro = matrix(c(numeric(d)), d, d)
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = ro + t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
    }
    for (i in 1:g) {
      para[[i]]$sigma = diag(diag(ro)) / n
    }
  }
  
  if (type == "VEI") {
    ro = array(0, c(d, d, g))
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro1 = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      ro[, , i] = ro1
    }
    
    In.iter = 0
    B = invB = diag(d)
    BB = 0
    while (In.iter < max.in.iter) {
      In.iter = In.iter + 1
      Lam.K = numeric(g)
      for (i in 1:g) {
        Lam.K[i] = tr(ro[, , i] %*% invB) / d / ni[i]
        BB = BB + ro[, , i] / Lam.K[i]
      }
      B = diag(diag(BB) / (det(diag(diag(
        BB
      )))) ^ (1 / d))
      invB = diag(1 / diag(B))
    }
    for (i in 1:g) {
      para[[i]]$sigma = B * Lam.K[i]
    }
  }
  
  if (type == "EVI") {
    ro = array(0, c(d, d, g))
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro1 = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      
      ro[, , i] = ro1
    }
    
    Bk = matrix(0, nrow = d, ncol = g)
    for (i in 1:g)
      Bk[, i] = diag(ro[, , i])
    Omeg.hat = apply(Bk, 2, prod) ^ (1 / d)
    Bk  = sweep(Bk, 2, 1 / Omeg.hat, FUN = "*")
    Omeg.hat = sum(Omeg.hat) / n
    
    for (i in 1:g) {
      para[[i]]$sigma = Omeg.hat * diag(Bk[, i], d)
    }
  }
  
  if (type == "VVI") {
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      para[[i]]$sigma =  diag(diag(ro)) / ni[i]
    }
  }
  
  if (type == "EEE") {
    ro = matrix(c(numeric(d)), d, d)
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = ro + t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
    }
    for (i in 1:g) {
      para[[i]]$sigma = ro / n
    }
  }
  
  if (type == "EEV") {
    A = 0
    D = list()
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      eig.k = eigen(ro)
      Omeg = diag(eig.k$val)
      A = A + Omeg
      D[[i]] = eig.k$vec
    }
    Lam.A = A / n
    for (i in 1:g) {
      para[[i]]$sigma = D[[i]] %*%  Lam.A %*% t(D[[i]])
    }
  }
  
  if (type == "VEV") {
    ro = array(0, c(d, d, g))
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro1 = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      
      ro[, , i] = ro1
    }
    
    In.iter = 0
    D = list()
    A = diag(d)
    lam.K = numeric(g)
    for (i in 1:g)
      D[[i]] = diag(d)
    while (In.iter < max.in.iter) {
      In.iter = In.iter + 1
      AA = 0
      for (i in 1:g) {
        eig.k = eigen(ro[, , i])
        Omeg = diag(eig.k$val)
        D[[i]] = eig.k$vec
        lam.K[i] = tr(ro[, , i] %*% D[[i]] %*%
                        matrixcalc::matrix.inverse(A) %*% t(D[[i]])) / d /
          ni[i]
        AA = AA + Omeg / lam.K[i]
      }
      A = AA / (det(AA)) ^ (1 / d)
    }
    for (i in 1:g) {
      para[[i]]$sigma = lam.K[[i]] * D[[i]] %*%  A %*% t(D[[i]])
    }
  }
  
  if (type == "EVE") {
    D  = diag(d)
    Wk = array(0, c(d, d, g))
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      
      Wk[, , i] = ro
    }
    
    Ak = apply(Wk, 3, function(z, D) {
      diag(t(D) %*% z %*% (D))
    }, D = D)
    Ak = apply(Ak, 2, function(z) {
      z / prod(z) ^ (1 / length(z))
    })
    D = newD(
      D = D,
      Wk = Wk,
      Ak = Ak,
      G = g,
      d = d
    )
    In.iter = 0
    while (In.iter < max.in.iter) {
      D = newD(
        D = D,
        Wk = Wk,
        Ak = Ak,
        G = g,
        d = d
      )
      Ak = apply(Wk, 3, function(z, D) {
        diag(t(D) %*% z %*% (D))
      }, D = D)
      Ak = apply(Ak, 2, function(z) {
        z / prod(z) ^ (1 / length(z))
      })
      #conv = c(testval(Wk=Wk,Ak=Ak,D=D,G=G), conv[1] )
      In.iter = In.iter + 1
    }
    lam =  0
    for (i in 1:g)
      lam = lam  + sum(diag(D %*% diag(1 / Ak[, i]) %*% t(D) %*% Wk[, , i]))
    lam = lam / (n * d)
    
    for (i in 1:g) {
      para[[i]]$sigma = D %*% diag(lam * Ak[, i]) %*% t(D)
    }
    
  }
  
  if (type == "VVE") {
    D  = diag(d)
    Wk = array(0, c(d, d, g))
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      
      Wk[, , i] = ro
    }
    
    Ak = apply(Wk, 3, function(z, D) {
      diag(t(D) %*% z %*% (D))
    }, D = D)
    Ak = apply(Ak, 2, function(z) {
      z / prod(z) ^ (1 / length(z))
    })
    D = newD(
      D = D,
      Wk = Wk,
      Ak = Ak,
      G = g,
      d = d
    )
    In.iter = 0
    while (In.iter < max.in.iter) {
      D = newD(
        D = D,
        Wk = Wk,
        Ak = Ak,
        G = g,
        d = d
      )
      Ak = apply(Wk, 3, function(z, D) {
        diag(t(D) %*% z %*% (D))
      }, D = D)
      Ak = apply(Ak, 2, function(z) {
        z / prod(z) ^ (1 / length(z))
      })
      #conv = c(testval(Wk=Wk,Ak=Ak,D=D,G=G), conv[1] )
      In.iter = In.iter + 1
    }
    lam = numeric(g)
    for (i in 1:g)
      lam[i] = sum(diag(D %*% diag(1 / Ak[, i]) %*% t(D) %*% Wk[, , i])) / d / ni[i]
    for (i in 1:g) {
      para[[i]]$sigma = D %*% diag(lam[i] * Ak[, i]) %*% t(D)
    }
  }
  
  if (type == "VEE") {
    ro = array(0, c(d, d, g))
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro1 = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      
      ro[, , i] = ro1
    }
    
    In.iter = 0
    
    C = 0
    for (i in 1:g) {
      C = C + ro[, , i]
    }
    C = C / (det(C)) ^ (1 / d)
    while (In.iter < max.in.iter) {
      In.iter = In.iter + 1
      Lam.k = numeric(g)
      CC = 0
      for (i in 1:g) {
        Lam.k[i] = tr(ro[, , i] %*% matrixcalc::matrix.inverse(C)) / d / ni[i]
        CC = CC + ro[, , i] / Lam.k[i]
      }
      C = CC / (det(CC)) ^ (1 / d)
    }
    for (i in 1:g) {
      para[[i]]$sigma = C * Lam.k[i]
    }
  }
  
  if (type == "EVV") {
    Lam = 0
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - t(eij) %*% S3zlam - t(S3zlam) %*% eij
      
      Lam = Lam + (det(ro)) ^ (1 / d)
      para[[i]]$sigma = ro / (det(ro)) ^ (1 / d)
    }
    for (i in 1:g) {
      para[[i]]$sigma = para[[i]]$sigma * Lam / n
    }
  }
  
  if (type == "VVV") {
    for (i in 1:g) {
      nu = para$nu[i]
      S1z = S(Bhat[, i], 1, 1, 1 / nu) * zhat[, i]
      S2z = S(Bhat[, i], w2hat[, i], w2tild[, i], 1) * zhat[, i]
      S3z = S(Bhat[, i], what[, i], wtild[, i], 1 / sqrt(nu)) * zhat[, i]
      eij = sweep(Y, 2, para[[i]]$mu)
      eS1z = eij * S1z
      S3zlam = S3z %*% t(para[[i]]$delta)
      ro = t(eS1z) %*% eij + sum(S2z) * para[[i]]$delta %*% t(para[[i]]$delta) - 
        t(eij) %*% S3zlam - t(S3zlam) %*% eij
      para[[i]]$sigma = ro / ni[i]
    }
  }
  return(para)
}



my.Outer = function(para, n, m, logli.new, zhat, Bhat,
                    true.clus = c(), end, begin){
  
  BIC = -2 *  logli.new + m * log(n) 
  #HBIC = -2 *  logli.new + (m_0 + g - 1) * log(n) + m_1 * sum (log(n * pi))
  
  ent = sum(zhat * log(zhat))
  ICL = BIC + 2 * ent
  #HICL = HBIC + 2 * ent
  
  AIC = 2 * (-logli.new + m)
  CLC = 2 * (-logli.new + ent)
  EDC = -2 * logli.new + 0.2 * m * sqrt(n)
  
  #post.cate = apply(zhat, 1, which.max)
  post.cate = matrix(apply(zhat, 1, order), length(para$pi))[length(para$pi),]
  di = table(post.cate)
  
  Outlier.post = numeric(n)
  for (j in 1:n) {
    Ind = post.cate[j]
    Outlier.post[j] = ifelse(Bhat[j, Ind] >= 0.5, "good", "bad")
  }
  
  if (is.numeric(true.clus) == T) {
    km.clus = post.cate
    tab <- table(true.clus, km.clus)
    MCR = error.rate(true.clus, km.clus)
    ARI = rand.index(true.clus, km.clus)[2]
    
    variation.information = mcclust::vi.dist(true.clus, km.clus)
    #agree.proportion = MixSim :: ClassProp(true.clus, km.clus)
    Mallows  = dendextend::FM_index_R(true.clus, km.clus)
    #jacc_ind = clusteval::cluster_similarity(true.clus, km.clus)
    
    out = list(
      para.len = m,
      logli = logli.new,
      para = para,
      ni = di,
      post.cate = post.cate,
      outlier = Outlier.post,
      AIC = AIC,
      BIC = BIC,
      ICL = ICL,
      EDC = EDC ,
      CLC = CLC,
      #HICL = HICL,
      #HBIC = HBIC,
      MCR = MCR,
      Mallows.Index = Mallows,
      variation.information = variation.information,
      Cross.Tab = tab,
      #Rand.index = RI[1],
      adjusted.Rand.Index = ARI,
      #jacc_ind = jacc_ind,
      time = end - begin
    )
  }
  else{
    out = list(
      para.len = m,
      logli = logli.new,
      para = para,
      ni = di,
      post.cate = post.cate,
      outlier = Outlier.post,
      AIC = AIC,
      BIC = BIC,
      ICL = ICL,
      EDC = EDC ,
      CLC = CLC,
      #HICL = HICL,
      #HBIC = HBIC,
      time = end - begin
    )
  }
  return(out)
}



error.rate= function (clust1, clust2) {
  clust1 <- unclass(as.ordered(clust1))
  clust2 <- unclass(as.ordered(clust2))
  if ((n = length(clust1)) != length(clust2)) {
    warning("error: length not equal")
    return
  }
  if ((g = length(table(clust1))) != length(table(clust2))) {
    warning("the number of clusters are not equal")
    return
  }
  permute <- function(a) {
    n <- length(a)
    if (n == 1) 
      f <- a
    else {
      nm <- gamma(n)
      f <- array(0, c(n, n * nm))
      j <- 1
      for (i in a) {
        f[1, (j - 1) * nm + 1:nm] <- i
        f[-1, (j - 1) * nm + 1:nm] <- permute(setdiff(a, 
                                                      i))
        j <- j + 1
      }
    }
    f
  }
  id <- 1:n
  cmb <- permute(1:g)
  nperm <- ncol(cmb)
  rate <- rep(0, nperm)
  for (i in 1:nperm) {
    tmp <- rep(0, g)
    tc <- rep(0, n)
    for (j in 1:g) tc[clust2 == j] = cmb[j, i]
    for (j in 1:g) {
      tmp1 <- 0
      for (k in (1:g)[-j]) tmp1 <- tmp1 + length(intersect(id[clust1 == 
                                                                j], id[tc == k]))
      tmp[j] <- tmp1
    }
    rate[i] <- sum(tmp)/n
  }
  min(rate)
}



rand.index = function (LabelA, LabelB) {
  u <- unclass(as.ordered(LabelA))
  v <- unclass(as.ordered(LabelB))
  if ((N <- length(u)) != length(v)) 
    stop("Label A and B does not match!")
  row <- max(u)
  col <- max(v)
  nvect <- array(0, c(row, col))
  for (i in 1:row) {
    for (j in 1:col) {
      nvect[i, j] <- sum(u == i & v == j)
    }
  }
  SumsA <- rowSums(nvect)
  SumsB <- colSums(nvect)
  a = 0
  for (i in 1:row) a = a + choose(SumsA[i], 2)
  b = 0
  for (j in 1:col) b = b + choose(SumsB[j], 2)
  c <- a * b/choose(N, 2)
  d = 0
  for (i in 1:row) {
    for (j in 1:col) {
      d = d + choose(nvect[i, j], 2)
    }
  }
  arj <- (d - c)/((a + b)/2 - c)
  a = d
  b = 0
  for (l in 1:row) {
    for (i in 1:(col - 1)) {
      for (j in (i + 1):col) {
        b = b + nvect[l, i] * nvect[l, j]
      }
    }
  }
  c = 0
  for (l in 1:col) {
    for (i in 1:(row - 1)) {
      for (j in (i + 1):row) {
        c = c + nvect[i, l] * nvect[j, l]
      }
    }
  }
  rad = (choose(N, 2) - b - c)/choose(N, 2)
  ind <- c(rad, arj)
  names(ind) <- c("Rand Index (RI)", "Adjusted Rand Index (ARI)")
  ind
}
