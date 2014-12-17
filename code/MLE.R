
## Useful function for computing Maximum Likelihood

## Computes the log-likelihood for degrees of membership and param
ll <- function(G, theta, X) {
  
  N <- dim(G)[1]
  P <- length(theta)
  
  high_v <- function(p, n) {
    log(G[n, 1] * theta[[p]]$low[X[n, p] + 1] +
          G[n, 2] * theta[[p]]$high[X[n, p] + 1])
  } # Compute each n, p contribution
  
  sum(sapply(1:P, high_v, 1:N)) # Sum all contributions
  
}

## Find MLE using coordinate ascent approach
gomMLE <- function(X, G0, theta0, epsilon=1e-5, verbose=T) {
  
  N <- dim(G0)[1]
  P <- length(theta0)
  
  # Current parameter and log-likelihood
  g.l <- G0[, 1]
  theta <- theta0
  maxlik <- ll(G0, theta0, X)
  maxlik.minus_one <- -Inf
  count <- 0 # Iteration count
  
  # Iteration number
  it <- 0
  
  while (abs(maxlik - maxlik.minus_one) > epsilon) {
    
    it <- it + 1
    print(it)
    
    # Update previous maxlik and count
    maxlik.minus_one <- maxlik
    count <- count + 1
    
    # Optim on g compontents
    
    for (n in 1:N) {
      #print(n)
      opt.g <- optim(g.l[n], function(g) {
        if (n == 1) {
          g.l_col <-c(g, g.l[2:N])
        } else if (n == N) {
          g.l_col <-c(g.l[1:(N-1)], g)
        } else {
          g.l_col <-c(g.l[1:(n-1)], g, g.l[(n+1):N])
        }
        
        ll(matrix(data=c(g.l_col, 1 - g.l_col), nrow=N, ncol=2), theta, X)
      },
      method="L-BFGS-B", control=list(fnscale=-1), lower=1e-7, upper=1 - 1e-7)
      
      g.l[n] <- opt.g$par
      maxlik <- opt.g$value
    }
    if (verbose) {
      print(paste("optim on g, iteration", count))
      print(paste("max likelihood:", maxlik))
    }
    
    # Optim on theta.L compontents
    for (p in 1:P) {
      #print(p)
      thetaL <- theta[[p]]$low # Current thetas low
      xi.L0 <- log.scaling(thetaL)
      opt.L <- optim(xi.L0, function(xi.L) {
        theta.cur <- theta
        theta.cur[[p]]$low <- exp.scaling(xi.L)
        
        ll(matrix(data=c(g.l, 1 - g.l), nrow=N, ncol=2), theta.cur, X)
      },
      method="L-BFGS-B", control=list(fnscale=-1))
      
      xi.Lopt <- opt.L$par
      maxlik <- opt.L$value
      
      theta[[p]]$low <- exp.scaling(xi.Lopt)
    }
    if (verbose) {
      print(paste("optim on theta_l, iteration", count))
      print(paste("max likelihood:", maxlik))
    }
    
    # Optim on theta.H components
    for (p in 1:P) {
      #print(p)
      thetaH <- theta[[p]]$high # Current thetas high
      xi.H0 <- log.scaling(thetaH)
      opt.H <- optim(xi.H0, function(xi.H) {
        theta.cur <- theta
        theta.cur[[p]]$high <- exp.scaling(xi.H)
        
        ll(matrix(data=c(g.l, 1 - g.l), nrow=N, ncol=2), theta.cur, X)
      },
      method="L-BFGS-B", control=list(fnscale=-1), lower=-100, upper=100)
      
      xi.Hopt <- opt.H$par
      maxlik <- opt.H$value
      
      theta[[p]]$high <- exp.scaling(xi.Hopt)
    }
    if (verbose) {
      print(paste("optim on theta_h, iteration", count))
      print(paste("max likelihood:", maxlik))
    }
    
  }
  
  list("G.hat"=matrix(c(g.l, 1 - g.l), ncol=2, byrow=F), "theta.hat"=theta)
}

## Helper scaling functions, to obtain a constrained opt problem
exp.scaling <- function(x) {
  
  exp(x) / (sum(exp(x)))
}

log.scaling <- function(u) {
  log(u) - sum(log(u))
}