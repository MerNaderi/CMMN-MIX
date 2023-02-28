
t.slice.gen = function(n1, n2, mu, sigma, lambda, nu, pi){
  
  # n1 number of good data
  # n2 number of noises
  #Mu and lambda are G * p
  #sigma is p * p * G array
  Class = numeric(n1)
  Point.Type = rep("Good", n1)
  p = dim(mu)[2]
  x = matrix(NA, nrow = n1, ncol = p)
  
  for(i in 1:n1){
    Z = rmultinom(1, 1, PI)
    j = Class[i] = which.max(Z)
    if( j == 2){
      Condi = F
      while(Condi == F){
        u = rgamma(1, nu[j]/2, nu[j]/2)
        z = mvtnorm :: rmvnorm(1, mean = rep(0, p), sigma = sigma[, , j])
        XX = mu[j, ] + z / sqrt(u)
        Condi = ((XX[2] + 0.8 * XX[1] + 5) >= 0)
      } 
    }
    if(j == 1){
      Condi = F
      while(Condi == F){
        u = rgamma(1, nu[j]/2, nu[j]/2)
        z = mvtnorm :: rmvnorm(1, mean = rep(0, p), sigma = sigma[, , j])
        XX = mu[j, ] + z / sqrt(u)
        Condi = ((XX[2] + 0.8 * XX[1] + 12) >= 0)
      } 
    }
    if(j == 3){
      Condi = F
      while(Condi == F){
        u = rgamma(1, nu[j]/2, nu[j]/2)
        z = mvtnorm :: rmvnorm(1, mean = rep(0, p), sigma = sigma[, , j])
        XX = mu[j, ] + z / sqrt(u)
        Condi = ((XX[2] + 0.8 * XX[1] - 8) >= 0)
      } 
    }
    x[i, ] = XX
  }
  
  x = rbind(x, cbind(runif(n2, -10, 10), runif(n2, -10, 10)))
  Point.Type = c(Point.Type, rep("bad", n2))
  
  out = list(X = x, Class = Class, Type = Point.Type)
  
  return(out)
}


d = 2
G = 3
sigma = array(0, dim = c(d, d, G))
sigma[,,1] = matrix( c(1, 0.7, 0.7, 1), ncol=d, byrow = T)
sigma[,,2] = matrix( c(2.33, 0.5, 0.5, 2), ncol=d, byrow = T)
sigma[,,3] = matrix( c(1.6 ,  0.58,  0.58, 2), ncol=d, byrow = T)

mu = matrix(0, nrow = G, ncol = d)
mu[1,] = matrix(c(-3,-9) , byrow = T, nrow = 2)
mu[2,] = matrix(c(-5,0) , byrow = T, nrow = 2)
mu[3,] = matrix(c(7,2)  , byrow = T, nrow = 2)

nu = c(2, 6, 15)
PI = c(0.25, 0.35, 0.4)

n1 = 500
n2 = 25

data = t.slice.gen(n1, n2, mu, sigma, lambda, nu, PI)
Y = data$X
Type = c(rep(0, n1),rep(3, n2))

plot(Y[,1], Y[,2], ylab = expression(Y[2]),
     xlab = expression(Y[1]), main = "", pch = 21 + Type, 
     bg = c("deeppink3", "lemonchiffon3", "deepskyblue3", "lightpink3")[unclass(c(data$Class, rep(4, n2)))], 
     cex.axis = 1.1, cex = 1.2, cex.main = 1)





d = 2
G = 3
sigma = array(0, dim = c(d, d, G))
sigma[,,1] = matrix( c(1, 0.7, 0.7, 1), ncol=d, byrow = T)
sigma[,,2] = matrix( c(2.33, 0.5, 0.5, 2), ncol=d, byrow = T)
sigma[,,3] = matrix( c(1.6 ,  0.58,  0.58, 2), ncol=d, byrow = T)

mu = matrix(0, nrow = G, ncol = d)
mu[1,] = matrix(c(-3,-2.5) , byrow = T, nrow = 2)
mu[2,] = matrix(c(-2, 2) , byrow = T, nrow = 2)
mu[3,] = matrix(c(1, -1)  , byrow = T, nrow = 2)


t.slice.gen = function(n1, n2, mu, sigma, lambda, nu, pi){
  
  # n1 number of good data
  # n2 number of noises
  #Mu and lambda are G * p
  #sigma is p * p * G array
  Class = numeric(n1)
  Point.Type = rep("Good", n1)
  p = dim(mu)[2]
  x = matrix(NA, nrow = n1, ncol = p)
  
  for(i in 1:n1){
    Z = rmultinom(1, 1, PI)
    j = Class[i] = which.max(Z)
    if( j == 2){
      Condi = F
      while(Condi == F){
        u = rgamma(1, nu[j]/2, nu[j]/2)
        z = mvtnorm :: rmvnorm(1, mean = rep(0, p), sigma = sigma[, , j])
        XX = mu[j, ] + z / sqrt(u)
        Condi = ((XX[2] + 0.8 * XX[1] + 1) >= 0)
      } 
    }
    if(j == 1){
      Condi = F
      while(Condi == F){
        u = rgamma(1, nu[j]/2, nu[j]/2)
        z = mvtnorm :: rmvnorm(1, mean = rep(0, p), sigma = sigma[, , j])
        XX = mu[j, ] + z / sqrt(u)
        Condi = ((XX[2] + 0.8 * XX[1] + 9) >= 0)
      } 
    }
    if(j == 3){
      Condi = F
      while(Condi == F){
        u = rgamma(1, nu[j]/2, nu[j]/2)
        z = mvtnorm :: rmvnorm(1, mean = rep(0, p), sigma = sigma[, , j])
        XX = mu[j, ] + z / sqrt(u)
        Condi = ((XX[2] + 0.8 * XX[1] - 2) >= 0)
      } 
    }
    x[i, ] = XX
  }
  
  x = rbind(x, cbind(runif(n2, -10, 10), runif(n2, -10, 10)))
  Point.Type = c(Point.Type, rep("bad", n2))
  
  out = list(X = x, Class = Class, Type = Point.Type)
  
  return(out)
}



nu = c(2, 6, 15)
PI = c(0.25, 0.35, 0.4)
data = t.slice.gen(n1, n2, mu, sigma, lambda, nu, PI)
Y = data$X
Type = c(rep(0, n1),rep(3, n2))

plot(Y[,1], Y[,2], ylab = expression(Y[2]),
     xlab = expression(Y[1]), main = "", pch = 21 + Type, 
     bg = c("deeppink3", "lemonchiffon3", "deepskyblue3", "lightpink3")[unclass(c(data$Class, rep(4, n2)))], 
     cex.axis = 1.1, cex = 1.2, cex.main = 1)


