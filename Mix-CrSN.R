#rm(list=ls(all=TRUE))
MIX.CrSN.EM = function (Y,
                        para = NULL,
                        Family = c("rSN", "ErSN"),
                        tol = 1e-5,
                        type = c(
                          "EII",
                          "VII",
                          "EEI",
                          "VEI",
                          "EVI",
                          "VVI",
                          "EEE",
                          "VEE",
                          "EVE",
                          "VVE",
                          "EEV",
                          "VEV",
                          "EVV",
                          "VVV"
                        ),
                        Stp.rule = c("Log.like", "Atiken"),
                        max.in.iter = 200,
                        max.iter = 100,
                        per = 100,
                        print = T,
                        class = c(),
                        lower = -100,
                        upper = 100)
{
  begin = proc.time()[1]
  
  TT = c(
    "EII",
    "VII",
    "EEI",
    "VEI",
    "EVI",
    "VVI",
    "EEE",
    "VEE",
    "EVE",
    "VVE",
    "EEV",
    "VEV",
    "EVV",
    "VVV"
  )
  if (all(type != TT))
    stop("Sigma structure is not correctly defined.")
  Y = as.matrix(Y)
  n = dim(Y)[1]
  d = dim(Y)[2]
  g = length(para$pi)
  
  d.rSN = function(Y,
                   mu,
                   Sigma,
                   lambda,
                   tau,
                   Family = c("rSN", "ErSN"),
                   log = F) {
    Y = as.matrix(Y, ncol = ncol(Y))
    if (Family == "rSN")
      tau = 0
    Omega = Sigma + lambda %*% t(lambda)
    xmu = sweep(Y, 2, mu)
    inv.Omega = matrixcalc::matrix.inverse(Omega)
    A = (c(lambda %*% inv.Omega %*% t(xmu)) + tau) /
      sqrt(1 - mahalanobis(lambda, 0, Omega))
    
    PDF = pnorm(A, log.p = TRUE) +
      mvtnorm::dmvnorm(Y,
                       mean = mu,
                       sigma = Omega,
                       log = TRUE) -
      pnorm(tau, log.p = TRUE)
    
    ifelse(log == T, return(PDF), return(exp(PDF)))
  }
  
  d.CrSN = function(Y,
                    para = NULL,
                    tau = NULL,
                    nu = NULL,
                    gam = NULL,
                    Family = c("rSN", "ErSN"),
                    log = F) {
    mu = para$mu
    Sigma = para$sigma
    lambda = para$delta
    
    Part1 = d.rSN(
      Y,
      mu = mu,
      Sigma = nu * Sigma,
      lambda = sqrt(nu) * lambda,
      tau = tau,
      Family = Family,
      log = T
    )
    Part2 = d.rSN(
      Y,
      mu = mu,
      Sigma = Sigma,
      lambda = lambda,
      tau = tau,
      Family = Family,
      log = T
    )
    Out = Part2 + log(gam + (1 - gam) * exp(Part1 - Part2))
    
    ifelse(log == T, return(Out), return(exp(Out)))
  }
  
  Log.like.CrSN <-
    function(Y,
             para = NULL,
             Family = c("rSN", "ErSN")) {
      fx = matrix(0, nrow = nrow(Y), ncol = length(para$pi))
      if (Family == "rSN")
        para$tau = rep(0, length(para$pi))
      
      for (k in 1:length(para$pi)) {
        fx[, k] = d.CrSN(
          Y,
          para = para[[k]],
          tau = para$tau[k],
          nu = para$nu[k],
          gam = para$gam[k],
          Family = Family,
          log = F
        )
      }
      
      val = apply(fx, 1, weighted.sum, wt = para$pi)
      #val[val == 0] <- .Machine$double.xmin
      val = sum(log(val))
      return(val)
    }
  
  m_0 = nSigma.par(Sig.Struc = type, p = d, G = g)
  
  if (Family == "rSN") {
    para$tau = rep(0, g)
    m = g * 2 * d + m_0 + 3 * g - 1
  }
  else{
    m = g * 2 * d + m_0 + 4 * g - 1
  }
  
  lk = logli.old = Log.like.CrSN(Y, para, Family = Family)
  iter = 0
  if (isTRUE(print)) {
    cat(
      "Fitting mixture of contaminated",
      Family,
      "distributions",
      " with G =",
      g,
      "components and",
      type,
      "structure of covariance matrix",
      "\n"
    )
    cat("Initial Log-likelihood = ", lk, "\n")
  }
  #Para.old = para
  repeat {
    iter = iter + 1
    ### E step
    zhat = Bhat = what = w2hat = wtild = w2tild = matrix(0, nrow = nrow(Y), ncol = g)
    
    for (k in 1:g) {
      PP = para[[k]]
      mu = PP$mu
      Sigma = PP$sigma
      lambda = PP$delta
      nu = para$nu[k]
      tau = para$tau[k]
      gam = para$gam[k]
      
      ZZ =  para$pi[k] * d.CrSN(Y,
                                para = PP,
                                tau = tau,
                                nu = nu,
                                gam = gam,
                                Family = Family,
                                log = F)
      #ZZ[ZZ == 0] <- .Machine$double.xmin
      zhat[, k] = ZZ
      
      Part1 = d.rSN(Y,
                    mu,
                    nu * Sigma,
                    sqrt(nu) * lambda,
                    tau,
                    Family = Family,
                    log = T)
      Part2 = d.rSN(Y,
                    mu,
                    Sigma,
                    lambda,
                    tau,
                    Family = Family,
                    log = T)
      Out = gam + (1 - gam) * exp(Part1 - Part2)
      Bhat[, k] = gam / Out
      
      # conditional expectations of good points
      O.inv = matrixcalc::matrix.inverse(Sigma + lambda %*% t(lambda))
      delta = sqrt(1 - mahalanobis(lambda, 0, O.inv, inverted = T))
      xmu = sweep(Y, 2, mu)
      xi = c(lambda %*% O.inv %*% t(xmu))
      alpha = (xi +  tau) / delta
      HF = exp(dnorm(alpha, log = TRUE) - pnorm(alpha, log.p = TRUE))
      what[, k] = xi + delta * HF
      w2hat[, k] = xi ^ 2 + delta ^ 2 - delta ^ 2 * alpha * HF + 2 * xi * delta * HF
      
      # conditional expectations of bad points
      O.inv = matrixcalc::matrix.inverse(nu * Sigma + nu * lambda %*% t(lambda))
      delta = sqrt(1 - mahalanobis(sqrt(nu) * lambda, 0, O.inv, inverted = T))
      xmu = sweep(Y, 2, mu)
      xi = c(sqrt(nu) * lambda %*% O.inv %*% t(xmu))
      alpha = (xi +  tau) / delta
      HF = exp(dnorm(alpha, log = TRUE) - pnorm(alpha, log.p = TRUE))
      wtild[, k] = xi + delta * HF
      w2tild[, k] = xi ^ 2 + delta ^ 2 - delta ^ 2 * alpha * HF + 2 * xi * delta * HF
    }
    
    for (k in 1:n) 
      if (all(zhat[k,] == 0))
        zhat[k,] = .Machine$double.xmin
    nij = rowSums(zhat)
    zhat = zhat / nij
    
    ni = colSums(zhat)
    para$pi = ni / n
    para$gam = colSums(zhat * Bhat) / ni
    for (k in 1:g) para$gam[k] = max(para$gam[k], 0.51)
      
    ### M---step
    para = M.step(Y, para, zhat, Bhat,
                  what, w2hat, wtild, w2tild,
                  type = type, g = g, max.in.iter = max.in.iter)
    
    for (k in 1:g) {
      f.nu = function(nuu) {
        -sum(
          zhat[, k] * d.CrSN(
            Y, para = para[[k]],
            tau = para$tau[k], nu = nuu,
            gam = para$gam[k],
            Family = Family,
            log = T
          )
        )
      }
      para$nu[k] = nlminb(para$nu[k], f.nu, lower = 1.001, upper = 200)$par
      #para$nu[k] = optim(para$nu[k], f.nu, method = "L-BFGS-B", lower = 1.001, upper = 20)$par
    }
    
    RA = rank(para$pi)
    paraO = para
    for(k in 1:g){
      Sel = which(RA == k)
      para[[k]] = paraO[[Sel]]
      para$pi[k] = paraO$pi[Sel]
      para$nu[k] = paraO$nu[Sel]
      para$gam[k] = paraO$gam[Sel]
    }
    
    if (Family == "rSN")
      logli.new = Log.like.CrSN(Y, para, Family = Family)
    
    
    if (Family == "ErSN") {
      for (k in 1:g) {
        PP = para[[k]]
        mu = PP$mu
        Sigma = PP$sigma
        lambda = PP$delta
        nu = para$nu[k]
        tau = para$tau[k]
        gam = para$gam[k]
        
        ZZ =  para$pi[k] * d.CrSN(Y,
                                  para = PP,
                                  tau = tau,
                                  nu = nu,
                                  gam = gam,
                                  Family = Family,
                                  log = F)
        #ZZ[ZZ == 0] <- .Machine$double.xmin
        zhat[, k] = ZZ
      }
      for (k in 1:n) 
        if (all(zhat[k,] == 0))
          zhat[k,] = .Machine$double.xmin
      nij = rowSums(zhat)
      zhat = zhat / nij
      
      for (k in 1:g) {
        f.tau = function(tau) {
          -sum(
            zhat[, k] * d.CrSN(
              Y, para = para[[k]],
              tau, nu = para$nu[k],
              gam = para$gam[k],
              Family = "ErSN",
              log = T
            )
          )
        }
        para$tau[k] = nlminb(para$tau[k], f.tau, lower = lower, upper = upper)$par
      }
      logli.new = Log.like.CrSN(Y, para, Family = Family)
    }
    
    if(is.infinite(logli.new)) logli.new = logli.old
    # Aikten's method
    lk = c(lk, logli.new)
    if (Stp.rule == "Log.like")
      epsilon = (logli.new - logli.old) / abs(logli.old)
    else {
      epsilon = Stop.rule(lk)
    }
    
    diff = logli.new - logli.old
    if (iter %% per == 0 | is.na(diff))
    {
      if (isTRUE(print)) {
        cat(
          'iter =',
          iter,
          '\t logli =',
          logli.new,
          '\t diff =',
          diff,
          "\t",
          Stp.rule,
          "'s diff =",
          epsilon,
          '\n'
        )
        cat(paste(rep("-", 60), sep = "", collapse = ""), "\n")
      }
    }
    if (is.na(diff)) {
      logli.new = logli.old
      break
    }
    if (is.na(epsilon)) {
      logli.new = logli.old
      break
    }
    if (iter > 1)
      if (epsilon < tol | iter == max.iter)
        break
    logli.old = logli.new
  }
  if (!is.na(epsilon)) {
    for (k in 1:g) {
      PP = para[[k]]
      mu = PP$mu
      Sigma = PP$sigma
      lambda = PP$delta
      nu = para$nu[k]
      tau = para$tau[k]
      gam = para$gam[k]
      
      ZZ =  para$pi[k] * d.CrSN(Y,
                                para = PP,
                                tau = tau,
                                nu = nu,
                                gam = gam,
                                Family = Family,
                                log = F)
      #ZZ[ZZ == 0] <- .Machine$double.xmin
      zhat[, k] = ZZ
      
      Part1 = d.rSN(Y,
                    mu,
                    nu * Sigma,
                    sqrt(nu) * lambda,
                    tau,
                    Family = Family,
                    log = T)
      Part2 = d.rSN(Y, mu, Sigma, lambda, tau, Family = Family, log = T)
      Out = gam + (1 - gam) * exp(Part1 - Part2)
      Bhat[, k] = gam / Out
    }
    zhat = zhat / rowSums(zhat)
  }
  
  Conveg = (epsilon < tol)
  end = proc.time()[1]
  
  if (isTRUE(print)) {
    cat('Convergence ', Conveg, '\n')
    cat(paste(rep("-", 60), sep = "", collapse = ""), "\n")
    cat(
      'Fitting mixture of contaminated',
      Family,
      ' distributions takes',
      end - begin,
      'seconds\n'
    )
  }
  
  if (is.numeric(class) == T) {
    out = my.Outer(para, n, m, logli.new, zhat, Bhat, 
                   true.clus = class, end, begin)
  }
  else{
    out = my.Outer(para, n, m, logli.new, zhat, Bhat,
                   true.clus = c(), end, begin)
  }
  return(out)
}

