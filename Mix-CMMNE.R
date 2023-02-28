

MIX.CMMNE.EM = function (Y,
                         para,
                         null.tau = T,
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
                         per,
                         print = T,
                         class = c())
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
  
  num.sig = nSigma.par(Sig.Struc = type, p = d, G = g)
  
  if (null.tau) {
    para$tau = rep(NULL, g)
    m = g * 2 * d + num.sig + 3 * g - 1
  }
  else{
    m = g * 2 * d + num.sig + 4 * g - 1
  }
  
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
  
  Log.like.CMMNE <- function(Y, para = NULL) {
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
    val = sum(log(val))
    return(val)
  }
  
  if(null.tau) para$tau = rep(1, g)
  
  lk = logli.old = Log.like.CMMNE(Y, para)
  iter = 0
  if (isTRUE(print)) {
    cat(
      "Fitting mixture of CMMNE distributions with G =",
      g,
      "components and",
      type,
      "structure of covariance matrix",
      "\n"
    )
    cat("Initial Log-likelihood = ", lk, "\n")
  }
  
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
      
      ZZ =  para$pi[k] * d.CMMNE(Y,
                                 para = PP,
                                 tau = tau,
                                 nu = nu,
                                 gam = gam,
                                 log = F)
      #ZZ[ZZ == 0] <- .Machine$double.xmin
      zhat[, k] = ZZ
      
      Part1 = d.MMNE(Y, mu, nu * Sigma, sqrt(nu) * lambda, tau, log = T)
      Part2 = d.MMNE(Y, mu, Sigma, lambda, tau, log = T)
      Out = gam + (1 - gam) * exp(Part1 - Part2)
      Bhat[, k] = gam / Out
      
      
      # conditional expectations of good points
      inv.sig = matrixcalc::matrix.inverse(Sigma)
      delta = sqrt(mahalanobis(lambda, 0, Sigma))
      xmu = sweep(Y, 2, mu)
      A = (c(lambda %*% inv.sig %*% t(xmu)) - 1) / delta
      if (null.tau == F)
        A = (c(lambda %*% inv.sig %*% t(xmu)) - tau) / delta
      
      HF = exp(dnorm(A, log = TRUE) - pnorm(A, log.p = TRUE))
      what[, k] = A / delta + HF / delta
      w2hat[, k] = A * what[, k] / delta + 1 / delta ^ 2
      
      # conditional expectations of bad points
      inv.sig = matrixcalc::matrix.inverse(nu * Sigma)
      delta = sqrt(mahalanobis(sqrt(nu) * lambda, 0, nu * Sigma))
      xmu = sweep(Y, 2, mu)
      A = (c(sqrt(nu) * lambda %*% inv.sig %*% t(xmu)) - 1) / delta
      if (null.tau == F)
        A = (c(sqrt(nu) * lambda %*% inv.sig %*% t(xmu)) - tau) / delta
      
      HF = exp(dnorm(A, log = TRUE) - pnorm(A, log.p = TRUE))
      wtild[, k] = A / delta + HF / delta
      w2tild[, k] = A / delta * wtild[, k] + 1 / delta ^ 2
    }
    
    for (k in 1:n)
      if (all(zhat[k,] == 0))
        zhat[k,] = .Machine$double.xmin
    
    nij = rowSums(zhat)
    zhat = zhat / nij
    
    ni = colSums(zhat)
    para$pi = ni / n
    para$gam = colSums(zhat * Bhat) / ni
    para$gam = ifelse(para$gam < 0.51, 0.51, para$gam)
    
    if (null.tau == F)
      para$tau = ni / colSums(zhat * S(Bhat, what, wtild, 1))
    
    ### M---step
    para = M.step(Y, para, zhat, Bhat,
                  what, w2hat, wtild, w2tild,
                  type = type, g = g, max.in.iter = max.in.iter)
    
    for (k in 1:g) {
      f.nu = function(nu) {
        -sum(
          zhat[, k] * d.CMMNE(
            Y, para = para[[k]],
            tau = para$tau[k], nu = nu,
            gam = para$gam[k],
            log = T
          )
        )
      }
      para$nu[k] = nlminb(para$nu[k], f.nu, lower = 1.001, upper = 20)$par
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
    
    logli.new = Log.like.CMMNE(Y, para)
    if(is.infinite(logli.new)) logli.new = logli.old
    #logli.new
    
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
    
    ZZ =  para$pi[k] * d.CMMNE(Y, para = PP, tau = tau, nu = nu, gam = gam, log = F)
    #ZZ[ZZ == 0] <- .Machine$double.xmin
    zhat[, k] = ZZ
    
    Part1 = d.MMNE(Y, mu, nu * Sigma, sqrt(nu) * lambda, tau, log = T)
    Part2 = d.MMNE(Y, mu, Sigma, lambda, tau, log = T)
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
    cat('Fitting Mixture of CMMNE distributions takes',
        end - begin,
        'seconds\n\n')
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



