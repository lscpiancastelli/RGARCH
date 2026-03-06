#' Simulate R-GARCH data
#'
#'
#' @param N Length of the time series
#' @param n Number of items in the ranking
#' @param phi0 Intercept term in the mean equation
#' @param phi Coefficients for the lagged distance terms in the mean equation
#' @param alpha Coefficients for the lagged mean terms in the mean equation
#' @param init Initial rankings for the first max(p,q)+1 time points (optional)
#' @param burn_in Number of initial time points to discard (default: 100)
#' @return A list containing the simulated rankings, spread trajectory and mean trajectory
#' @export
rRGARCH = function(N, n, phi0, phi = NULL, alpha = NULL, init, burn_in = 100){

  if(is.null(phi)) phi = numeric(0)
  if(is.null(alpha)) alpha = numeric(0)

  p = length(phi)
  q = length(alpha)

  mu = theta = numeric(N + burn_in)
  Pi = matrix(NA, nrow = N + burn_in, ncol = n)

  lag_max = max(p, q)

  mu[1:(lag_max + 1)] = phi0

  if(missing(init)){
    Pi[1:(lag_max + 1), ] = t(replicate(lag_max + 1, sample(1:n, size = n, replace = FALSE)))
  }else{
    Pi[1:(lag_max + 1), ] = init
  }

  for(t in (lag_max + 2):(N + burn_in)){

    phi_part = 0
    if(p > 0){
      phi_part = matrix(phi, nrow = 1) %*%
        matrix(sapply(1:p, dist_lag, t, Pi), ncol = 1)
    }

    alpha_part = 0
    if(q > 0){
      alpha_part = matrix(alpha, nrow = 1) %*%
        matrix(mu[(t-1):(t-q)], ncol = 1)
    }

    mu[t] = phi0 + phi_part + alpha_part

    theta[t] = get_theta(mu[t], n)

    Pi[t, ] = BayesMallows::sample_mallows(Pi[t-1, ], alpha0 = theta[t] * n,
                             n_samples = 1, metric = "kendall")
  }

  Pi = Pi[(burn_in + 1):(N + burn_in), ]

  return(list(
    data = Pi,
    theta = theta[(burn_in + 1):(N + burn_in)],
    mu = mu[(burn_in + 1):(N + burn_in)]
  ))
}


lagonedist_t = function(t, Pi){

  p =1
  pi_m1 = Pi[t-p,]
  pi_lagged = Pi[t-p-1,]

  d = Rankcluster::distKendall(pi_m1, pi_lagged, type = 'ranking')
  return(d)
}

#' Compute the lag-1 Kendall distance for each time point in the ranking time series
#' @param Pi A matrix of rankings, rows correspond to the time points
#' @return A vector of lag-1 Kendall distances starting from t=3
#' @export
lagonedist = function(Pi){
  return(sapply(3:nrow(Pi), lagonedist_t, Pi = Pi))
}

dist_lag = function(p, t, Pi){

  pi_m1 = Pi[t-p,]
  pi_lagged = Pi[t-p-1,]

  d = Rankcluster::distKendall(pi_m1, pi_lagged, type = 'ranking')
  return(d)
}

f_aux2 = function(theta, mu, n){
  return( (Mean_Kendall(theta, n)-mu)^2 )
}


get_theta = function(mu, n){
  op = optim(0.1, f_aux2, mu = mu, n=n, method ='Brent', lower = 0, upper =15, control = list(abstol = 10^(-30)))
  return(op$par)
}


Mean_Kendall  = function(theta, n){

  t1 = (n*exp(-theta))/(1-exp(-theta))
  t2 = sapply(1:n, function(x, theta){ (x*exp(-x*theta))/(1-exp(-x*theta))}, theta )

  return(t1 - sum(t2))

}

mgf_j = function(j, k, theta){
  num = 1 - exp(-theta*(k-j+1))
  den = (k-j+1)*(1-exp(-theta))
  return(num/den)
}

kendall_psi = function(n, theta){
  return(factorial(n)*prod(mgf_j(seq(1,n-1,1), n, theta)))
}


#' Fits the R-GARCH model to ranking time series data
#'
#' @param Pi A matrix of rankings, rows correspond to the time points
#' @param p The order of the autoregressive component (default: 0)
#' @param q The order of the moving average component (default: 0)
#' @param inits Optional initial parameter values for optimization. If not provided, default starting values will be used based on the specified p and q.
#' @return A list containing the simulated rankings, spread trajectory and mean trajectory
#' @export
fitRGARCH = function(Pi, p = 0, q = 0, inits){

  if(p == 0 && q == 0){

    if(missing(inits)){
      return(fit_null(Pi))
    } else{
      return(fit_null(Pi, inits))
    }

  } else if(p > 0 && q == 0){

    if(missing(inits)){
      return(fit_AR(Pi, p))
    } else{
      return(fit_AR(Pi, p, inits))
    }

  } else if(p == 0 && q > 0){

    if(missing(inits)){
      return(fit_MA(Pi, q))
    } else{
      return(fit_MA(Pi, q, inits))
    }

  } else{

    if(missing(inits)){
      return(fit_AR_MA(Pi, p, q))
    } else{
      return(fit_AR_MA(Pi, p, q, inits))
    }

  }

}

fit_null = function(Pi, inits){

  if(missing(inits)){
    par =1
  } else{
    par = inits
  }

  op = optim(par=par, fn=loglik_ingarch_rank_null,
             control = list(fnscale =-1, factr = 10^(-50)),
             Pi= Pi, method = 'BFGS',

             hessian = TRUE)
  op$p = 0
  op$q=0
  op$Pi = Pi
  return(op)

}

loglik_ingarch_rank_null = function(pars, Pi){

  p=0
  q=0
  phi0 = pars

  N = nrow(Pi)
  n = ncol(Pi)

  mu = theta = numeric(N)
  mu[1:(max(p,q)+1)] = phi0
  loglik = numeric(N)

  for(t in (max(p,q)+2):(N)){

    mu[t] = phi0

    theta[t] = get_theta(mu[t], n)
    loglik[t] =-theta[t]*Rankcluster::distKendall(Pi[t, ], Pi[t-1,], type = 'ranking') - log(kendall_psi(n, theta[t]))


  }

  return(sum(loglik))
}


fit_AR = function(Pi, p, inits){

  if(missing(inits)){
    dists = lagonedist(Pi)
    par =c(mean(dists), rep(0.1, p))
  } else{
    par = inits
  }

  op = optim(par=par, fn=loglik_ingarch_rank_AR,
             control = list(fnscale =-1, factr = 10^(-50)),
             p=p,  Pi= Pi, method = 'L-BFGS-B',
             lower = c(10^(-50), 10^(-50)), upper = c(200, .9999),
             hessian = TRUE)
  op$p =p
  op$q=0
  op$Pi = Pi
  return(op)

}


loglik_ingarch_rank_AR = function(pars, p, Pi){

  q=0
  phi0 = pars[1]
  phi = pars[2:(p+1)]

  N = nrow(Pi)
  n = ncol(Pi)

  mu = theta = numeric(N)
  mu[1:(max(p,q)+1)] = phi0
  loglik = numeric(N)

  for(t in (max(p,q)+2):(N)){

    mu[t] = phi0 + matrix(phi, nrow=1)%*%matrix(sapply(1:p, dist_lag, t, Pi),ncol=1)


    theta[t] = get_theta(mu[t], n)
    loglik[t] =-theta[t]*Rankcluster::distKendall(Pi[t, ], Pi[t-1,], type = 'ranking') - log(kendall_psi(n, theta[t]))


  }

  return(sum(loglik))
}



fit_MA = function(Pi, q, inits ){

  if(missing(inits)){
    dists = lagonedist(Pi)
    par =c(mean(dists), rep(0.1, q))
  } else{
    par = inits
  }


  op = optim(par = par, fn=loglik_ingarch_rank_MA,
             control = list(fnscale =-1, factr = 10^(-30)),
             Pi= Pi, q=q, method = 'L-BFGS-B',
             lower = c(0.001, rep(.001, q)), upper = c(100, rep(.99, q)), hessian = FALSE)
  op$q =q
  op$p =0
  op$Pi = Pi
  return(op)

}



loglik_ingarch_rank_MA = function(pars, q, Pi){

  p=0
  phi0 = pars[1]
  alpha = pars[2:length(pars)]

  N = nrow(Pi)
  n = ncol(Pi)

  mu = theta = numeric(N)
  mu[1:(max(p,q)+1)] = phi0
  loglik = numeric(N)

  for(t in (max(p,q)+2):(N)){

    mu[t] = phi0 + matrix(alpha, nrow=1)%*%matrix(mu[(t-1):(t-q)])

    theta[t] = get_theta(mu[t], n)
    loglik[t] =-theta[t]*Rankcluster::distKendall(Pi[t, ], Pi[t-1,], type = 'ranking') - log(kendall_psi(n, theta[t]))


  }

  return(sum(loglik))
}

fit_AR_MA = function(Pi, p, q, inits ){

  if(missing(inits)){
    dists = lagonedist(Pi)
    par =c(mean(dists), rep(0.1, p), rep(0.1, q))
  } else{
    par = inits
  }

  op = optim(par = par, fn=loglik_ingarch_rank_AR_MA,
             control = list(fnscale =-1),
             p=p,  Pi= Pi, q=q, method = 'L-BFGS',
             lower = c(0.000000001, rep(.000000001, p+q)), upper = c(100, rep(.99, p+q)), hessian = TRUE)
  op$p =p; op$q =q; op$Pi = Pi
  return(op)

}


loglik_ingarch_rank_AR_MA = function(pars, p, q, Pi){

  phi0 = pars[1]
  phi = pars[2:(p+1)]
  alpha = pars[(p+2):length(pars)]

  N = nrow(Pi)
  n = ncol(Pi)

  mu = theta = numeric(N)
  mu[1:(max(p,q)+1)] = phi0
  loglik = numeric(N)

  for(t in (max(p,q)+2):(N)){

    mu[t] = phi0 + matrix(phi, nrow=1)%*%matrix(sapply(1:p, dist_lag, t, Pi),ncol=1)+
      matrix(alpha, nrow=1)%*%matrix(mu[(t-1):(t-q)])


    theta[t] = get_theta(mu[t], n)
    loglik[t] =-theta[t]*Rankcluster::distKendall(Pi[t, ], Pi[t-1,], type = 'ranking') - log(kendall_psi(n, theta[t]))


  }

  return(sum(loglik))
}


fittedRGARCH = function(fit){

  pars = fit$par
  p = fit$p
  q = fit$q
  Pi = fit$Pi

  dists = lagonedist(Pi)

  phi0 = pars[1]

  phi = if(p > 0) pars[2:(p + 1)] else numeric(0)

  alpha = if(q > 0) pars[(p + 2):(p + q + 1)] else numeric(0)

  N = nrow(Pi)
  n = ncol(Pi)

  mu = theta = numeric(N)
  mu[1:(max(p,q)+1)] = dists[1:(max(p,q)+1)]

  for(t in (max(p,q)+2):(N)){

    mu[t] = phi0 + matrix(phi, nrow=1)%*%matrix(sapply(1:p, dist_lag, t, Pi),ncol=1)+
      matrix(alpha, nrow=1)%*%matrix(mu[(t-1):(t-q)])


  }

  return(mu)
}


