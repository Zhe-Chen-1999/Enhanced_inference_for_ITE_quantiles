################################################################################
################# Helper functions for ITE quantiles inference #################
################################################################################
library(RIQITE)
library(extraDistr)
library(dplyr)

rank_score <- function( n, method.list = list(name = "Wilcoxon") ){
  
  if ( method.list$name == "DIM" ) {
    stop( "Can't calculate rank scores for DIM" )
  }
  
  if(method.list$name == "Wilcoxon"){
    score = c(1:n)
    score = score/max(score)
    return(score)
  }
  
  if(method.list$name == "Stephenson"){
    score = choose( c(1:n) - 1, method.list$s - 1 )
    score = score/max(score)
    return(score)
  }
  
}

#### (1-alpha) simultaneous confidence/prediction intervals (lower bound) for all effects quantiles among treated, control or all units ####
ci_lower_quantile_exp <- function( Z, Y, treat.method.list = list(name = "Stephenson", s = 6), control.method.list = list(name = "Stephenson", s = 6),
                                   score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL, alpha=0.05,  set = "treat", alpha.ratio.treat = 0.5,  tol = 10^(-3) ){
  n = length(Z)
  
  if(set == "treat"){
    ci.treat = RIQITE::ci_quantile(Z = Z, Y = Y, k.vec = (n-sum(Z)+1):n, alternative = "greater", method.list = treat.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm, alpha = alpha, tol = tol, switch = FALSE)
    ci.treat = as.numeric( ci.treat$lower )
    return(ci.treat)
  }
  
  if(set == "control"){
    if(is.null(stat.null)){
      stat.null.control = NULL
    }else{
      if(is.null(score)){
        score = rank_score(n, method.list = control.method.list) #### update later
      }
      stat.null.control = sum(score) - stat.null
    }
    
    if(is.null(Z.perm)){
      Z.perm.control = NULL
    }else{
      Z.perm.control = 1 - Z.perm
    }
    ci.control = RIQITE::ci_quantile(Z = 1-Z, Y = -Y, k.vec = (n-sum(1-Z)+1):n, alternative = "greater", method.list = control.method.list, score = score, stat.null = stat.null.control, nperm = nperm, Z.perm = Z.perm.control, alpha = alpha, tol = tol, switch = FALSE)
    ci.control = as.numeric( ci.control$lower )
    return(ci.control)
  }
  
  if(set == "all"){
    
    ci.treat = RIQITE::ci_quantile(Z = Z, Y = Y, k.vec = (n-sum(Z)+1):n, alternative = "greater", method.list = treat.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm, alpha = alpha*alpha.ratio.treat, tol = tol, switch = FALSE)
    ci.treat = as.numeric( ci.treat$lower )
    
    if(is.null(stat.null)){
      stat.null.control = NULL
    }else{
      if(is.null(score)){
        score = rank_score(n, method.list = control.method.list) #### update later
      }
      stat.null.control = sum(score) - stat.null
    }
    
    if(is.null(Z.perm)){
      Z.perm.control = NULL
    }else{
      Z.perm.control = 1 - Z.perm
    }
    ci.control = RIQITE::ci_quantile(Z = 1-Z, Y = -Y, k.vec = (n-sum(1-Z)+1):n, alternative = "greater", method.list = control.method.list, score = score, stat.null = stat.null.control, nperm = nperm, Z.perm = Z.perm.control, alpha = alpha*(1-alpha.ratio.treat), tol = tol, switch = FALSE)
    ci.control = as.numeric( ci.control$lower )
    
    ci.all = sort( c(ci.treat, ci.control) )
    
    cat("ci.treat:", ci.treat, "\n")
    cat("ci.control:", ci.control, "\n")
    
    return( ci.all)
  }
}

########################## Finite Population Inference #############################
#### correction term based on multivariate hypergeometric distribution ####

correct_hypergeom <- function(N = 100, n = 50, K.vec = c(60, 80), kk.vec = c(30, 40), ndraw = 10^4, hg.draw = NULL){
  if( length(K.vec) > 1 ){
    if(is.null(hg.draw)){
      hg.draw = rmvhyper(ndraw, diff( c(0, K.vec, N) ), n)
      hg.draw = hg.draw[, -1, drop=FALSE]
    }
    H.sum = apply( hg.draw[, ncol(hg.draw):1, drop=FALSE], 1, cumsum )
    H.sum = matrix(H.sum, nrow = ncol(hg.draw) )
    H.sum = H.sum[nrow(H.sum):1, , drop=FALSE]
    prob = mean( apply( H.sum - (n - kk.vec) > 0, 2, max ) )
    return(prob)
  }
  
  if(length(K.vec) == 1){
    prob = 1 - phyper(n-kk.vec, N-K.vec, K.vec, n)
    return(prob)
  }
}

## calculate the correction term Delta using the proposed choice of kk.vec (k_prime)
# Here the input "alpha" is the value of gamma * alpha
threshold_correct_hypergeom <- function( N = 100, n = 50, K.vec = c(60, 80), alpha = 0.05, ndraw = 10^4, hg.draw = NULL, tol = 10^(-2) ){
  
  if(length(K.vec) > 1){
    if(is.null(hg.draw)){
      hg.draw = rmvhyper(ndraw, diff( c(0, K.vec, N) ), n)
      hg.draw = hg.draw[, -1, drop=FALSE] 
    }
    
    kappa = NULL
    
    kappa.lower = 1/length(K.vec)
    kappa.upper = 1
    prob.lower = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = n-qhyper(1-kappa.lower*alpha, N - K.vec, K.vec, n), hg.draw = hg.draw )
    prob.upper = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = n-qhyper(1-kappa.upper*alpha, N - K.vec, K.vec, n), hg.draw = hg.draw )
    
    if(prob.upper <= alpha){
      kappa = kappa.upper
    }
    
    if(prob.lower > alpha){
      kappa = kappa.lower
    }
    
    while(is.null(kappa)){
      kappa.mid = (kappa.upper+kappa.lower)/2
      prob.mid = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = n-qhyper(1-kappa.mid*alpha, N - K.vec, K.vec, n), hg.draw = hg.draw )
      
      if(prob.mid <= alpha){
        kappa.lower = kappa.mid
      }else{
        kappa.upper = kappa.mid
      }
      
      if( kappa.upper - kappa.lower <= tol ){
        kappa = kappa.lower 
      }
    }
    
    kk.vec = n-qhyper(1-kappa*alpha, N - K.vec, K.vec, n)
    prob = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = kk.vec, hg.draw = hg.draw )
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }
  
  if(length(K.vec) == 1){
    kappa = 1
    temp = qhyper(1-alpha, N - K.vec, K.vec, n)
    kk.vec = n-qhyper(1-alpha, N - K.vec, K.vec, n)
    prob = correct_hypergeom( N = N, n = n, K.vec = K.vec, kk.vec = kk.vec)
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }
}


### generalize ####

ci_lower_quantile_gen <- function( Z, Y, N = 2*length(Z), K.vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N), 
                                   alpha=0.05, gamma = 0.5, ndraw = 10^4, treat.method.list = list(name = "Stephenson", s = 6),
                                   control.method.list = list(name = "Stephenson", s = 6), score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,  set = "all", alpha.ratio.treat = 0.5,  tol = 10^(-3) ){
  
  if(set == "all"){
    n.star = length(Z)
  }
  
  if(set == "control"){
    n.star = sum(1-Z)
  }
  
  if(set == "treat"){
    n.star = sum(Z)
  }
  
  if( N > n.star ){
    step1 = threshold_correct_hypergeom(N = N, n = n.star, K.vec = K.vec, alpha = gamma * alpha, ndraw = ndraw)
  }
  
  if(N == n.star){
    step1 = list(prob = 0, kk.vec = K.vec)
  }
  
  step2 = ci_lower_quantile_exp( Z, Y, treat.method.list = treat.method.list, control.method.list = control.method.list, 
                                 score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm, 
                                 alpha= alpha - step1$prob,  set = set, alpha.ratio.treat = alpha.ratio.treat,  tol = tol )
  step2 = c(-Inf, step2) # add lower confidence limit "-Inf" for case kk.vec = 0
  
  ci = step2[step1$kk.vec+1] # index +1 for case kk.vec = 0
  
  K.vec = sort(K.vec)
  conf.int = data.frame(k = K.vec, k_prime = step1$kk.vec, lower = ci, upper = Inf)
  return(conf.int)
}

# specify individual/simultaneous inference
ci_lower_quantile_generalize <- function(Z, Y, N, k_vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N), 
                                         alpha=0.05, gamma = 0.5, ndraw = 10^4, 
                                         treat.method.list = list(name = "Stephenson", s = 6), control.method.list = list(name = "Stephenson", s = 6), score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,  
                                         set = "all", alpha.ratio.treat = 0.5, tol = 10^(-3), simul = TRUE){
  
  k_vec = sort(k_vec)
  if (simul) {
    res = ci_lower_quantile_gen( Z = Z, Y = Y, N = N, K.vec = k_vec, 
                                 alpha = alpha, gamma = gamma, ndraw = ndraw, 
                                 treat.method.list = treat.method.list, control.method.list = control.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,
                                 set = set, alpha.ratio.treat = alpha.ratio.treat,  tol = tol )
  }
  else {
    res = NULL
    for (k in k_vec){
      res = rbind(res, 
                  ci_lower_quantile_gen( Z = Z, Y = Y, N = N, K.vec = c(k), 
                                         alpha = alpha, gamma = gamma, ndraw = ndraw, 
                                         treat.method.list = treat.method.list, control.method.list = control.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,
                                         set = set, alpha.ratio.treat = alpha.ratio.treat,  tol = tol))
    }
  }
  return(res)
}

################################################################################
#################  Methods to be compared for finite population ################
################################################################################
# M0: Method in Caughey et al. (2023) 
# M1: Method combining treated and control using Bonferroni correction (apply M0 twice)
# M2: Method using Berger and Boos (1994)'s approach

## The standard output for all methods: 
## J (number of testing quantiles) rows and three columns. Three columns are k, lower, and upper
################################################################################
# M0: Method in Caughey et al. (2023) 
method_caughey <- function(Z, Y, k_vec, method.list = list(name = "Wilcoxon"), nperm = 10^4,  alpha = 0.05){
  return(RIQITE::ci_quantile(Z = Z, Y = Y, k.vec = k_vec,
                             alternative = "greater", 
                             method.list = method.list,
                             nperm = nperm,  alpha = alpha, switch = FALSE))
}
################################################################################
# M1: Method combining treated and control using Bonferroni correction (apply M0 twice)
method_combine <- function(Z, Y, N, k_vec, 
                           control.method.list = list(name = "Stephenson", s = 6), 
                           treat.method.list = list(name = "Stephenson", s = 6),
                           simul = TRUE, stat.null = NULL, score = NULL, Z.perm = NULL, 
                           alpha=0.05, gamma = 0.5, alpha.ratio.treat = 0.5, 
                           ndraw = 10^5, nperm = 10^4, tol = 10^(-3)){
  
  return(ci_lower_quantile_generalize(Z = Z, Y = Y, N = N, k_vec = k_vec, simul = simul,
                                      alpha = alpha, gamma = gamma, ndraw = ndraw, control.method.list = control.method.list,
                                      treat.method.list = treat.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,
                                      set = "all", alpha.ratio.treat = alpha.ratio.treat, tol = tol))
}

## example

# method_combine( Z, Y, N = length(Z), k_vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N),simul = TRUE) 
################################################################################
# M2: Method using Berger and Boos (1994)'s approach
method_bergerboos <- function(Z, Y, N, k_vec, 
                              treat.method.list = list(name = "Stephenson", s = 6), 
                              control.method.list = list(name = "Stephenson", s = 6),
                              simul = TRUE, stat.null = NULL, score = NULL, Z.perm = NULL, 
                              alpha=0.05, gamma = 0.5, alpha.ratio.treat = 0.5, 
                              ndraw = 10^5, nperm = 10^4, tol = 10^(-3)){
  
  ci.gen1 = ci_lower_quantile_generalize(Z, Y, N, k_vec, treat.method.list = treat.method.list, 
                                         stat.null = stat.null, set = "treat", simul = simul,
                                         alpha=alpha/2, gamma = gamma, alpha.ratio.treat = alpha.ratio.treat, 
                                         score = score, Z.perm = Z.perm,
                                         ndraw = ndraw, nperm = nperm, tol = tol)
  
  ci.gen2 = ci_lower_quantile_generalize(Z, Y, N, k_vec, control.method.list = control.method.list, 
                                         stat.null = stat.null, set = "control", simul = simul,
                                         alpha=alpha/2, gamma = gamma, alpha.ratio.treat = alpha.ratio.treat, 
                                         score = score, Z.perm = Z.perm,
                                         ndraw = ndraw, nperm = nperm, tol = tol)
  return(pmax( ci.gen1, ci.gen2 ))
}

## example

# method_bergerboos( Z, Y, N = length(Z), k_vec = ceiling(c(0.6, 0.7, 0.8, 0.9)*N), simul = FALSE,
#                   method.list = list(name = "Stephenson", s = 6), stat.null = stat.null)

###################################################################################


########################## Super Population Inference #############################
#### correction term based on multinomial distribution ####
correct_multinomial <- function( n = 50, Beta.vec = c(0.6, 0.8), kk.vec = c(30, 40), ndraw = 10^4, mn.draw = NULL ){
  if( length(Beta.vec) > 1 ){
    if(is.null(mn.draw)){
      mn.draw = rmultinom(ndraw, size = n, prob = diff( c(0, Beta.vec, 1) ))
      mn.draw = mn.draw[-1, , drop=FALSE]
    }
    N.sum = apply( mn.draw[nrow(mn.draw):1, , drop=FALSE], 2, cumsum )
    N.sum = matrix(N.sum, nrow = nrow(mn.draw) )
    N.sum = N.sum[nrow(N.sum):1, , drop=FALSE]
    prob = mean( apply( N.sum - (n - kk.vec) > 0, 2, max ) )
    return(prob)
  }
  
  if(length(Beta.vec) == 1){
    prob = 1 - pbinom(n-kk.vec, size=n, prob=1-Beta.vec) 
    return(prob)
  }
}

## calculate the correction term Delta using the proposed choice of kk.vec (k_prime)
threshold_correct_multinomial <- function(n = 50, Beta.vec = c(0.6, 0.8), alpha = 0.05, ndraw = 10^4, mn.draw = NULL, tol = 10^(-3) ){
  
  if(length(Beta.vec) > 1){
    
    if(is.null(mn.draw)){
      mn.draw = rmultinom(ndraw, size = n, prob = diff( c(0, Beta.vec, 1) ))
      mn.draw = mn.draw[-1, , drop=FALSE]
    }
    
    kappa = NULL
    
    kappa.lower = 1/length(Beta.vec)
    kappa.upper = 1
    prob.lower = correct_multinomial(n = n, Beta.vec = Beta.vec, kk.vec = n-qbinom(1-kappa.lower*alpha, n, 1-Beta.vec), mn.draw = mn.draw )
    prob.upper = correct_multinomial(n = n, Beta.vec = Beta.vec, kk.vec = n-qbinom(1-kappa.upper*alpha, n, 1-Beta.vec), mn.draw = mn.draw )
    
    if(prob.upper <= alpha){
      kappa = kappa.upper
    }
    
    if(prob.lower > alpha){
      kappa = kappa.lower
    }
    
    while(is.null(kappa)){
      kappa.mid = (kappa.upper+kappa.lower)/2
      prob.mid = correct_multinomial(n = n, Beta.vec = Beta.vec, kk.vec = n-qbinom(1-kappa.mid*alpha, n, 1-Beta.vec), mn.draw = mn.draw  )
      
      if(prob.mid <= alpha){
        kappa.lower = kappa.mid
      }else{
        kappa.upper = kappa.mid
      }
      
      if( kappa.upper - kappa.lower <= tol ){
        kappa = kappa.lower 
      }
    }
    
    kk.vec = n-qbinom(1-kappa*alpha, n, 1-Beta.vec)
    prob = correct_multinomial(n = n, Beta.vec = Beta.vec, kk.vec = kk.vec, mn.draw = mn.draw )
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }
  
  if(length(Beta.vec) == 1){
    kappa = 1
    kk.vec = n-qbinom(1-alpha, n, 1-Beta.vec)
    prob = correct_multinomial(n = n, Beta.vec = Beta.vec, kk.vec = kk.vec)
    return(list(kappa = kappa, prob = prob, kk.vec = kk.vec))
  }
}

### generalize ####
ci_lower_quantile_gen_SP <- function( Z, Y, Beta.vec = c(0.5, 0.6, 0.7, 0.8, 0.9), 
                                      alpha=0.05, gamma = 0.5, ndraw = 10^5, treat.method.list = list(name = "Stephenson", s = 6),
                                      control.method.list = list(name = "Stephenson", s = 6), score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,  set = "all", alpha.ratio.treat = 0.5,  tol = 10^(-3) ){
  
  if(set == "all"){
    n.star = length(Z)
  }
  
  if(set == "control"){
    n.star = sum(1-Z)
  }
  
  if(set == "treat"){
    n.star = sum(Z)
  }
  
  step1 = threshold_correct_multinomial(n = n.star, Beta.vec = Beta.vec, alpha =  gamma * alpha, ndraw = ndraw)
  
  step2 = ci_lower_quantile_exp( Z, Y, treat.method.list = treat.method.list, control.method.list = control.method.list, 
                                 score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm, 
                                 alpha= alpha - step1$prob,  set = set, alpha.ratio.treat = alpha.ratio.treat,  tol = tol )
  step2 = c(-Inf, step2) 
  
  ci = step2[step1$kk.vec+1]
  
  Beta.vec = sort(Beta.vec)
  conf.int = data.frame(beta = Beta.vec, lower = ci, upper = Inf)
  return(conf.int)
}

# specify individual/simultaneous inference
ci_lower_quantile_generalize_SP <- function( Z, Y, beta_vec = c(0.5, 0.6, 0.7, 0.8, 0.9), 
                                             alpha=0.05, gamma = 0.5, ndraw = 10^4, treat.method.list = list(name = "Stephenson", s = 6),
                                             control.method.list = list(name = "Stephenson", s = 6), score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,  set = "all", alpha.ratio.treat = 0.5,  tol = 10^(-3), simul = TRUE ){
  
  beta_vec = sort(beta_vec)
  if (simul) {
    res = ci_lower_quantile_gen_SP( Z = Z, Y = Y, Beta.vec = beta_vec, 
                                    alpha = alpha, gamma = gamma, ndraw = ndraw, 
                                    treat.method.list = treat.method.list, control.method.list = control.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,
                                    set = set, alpha.ratio.treat = alpha.ratio.treat,  tol = tol )
  }
  else {
    res = NULL
    for (beta in beta_vec){
      res = rbind(res, 
                  ci_lower_quantile_gen_SP( Z = Z, Y = Y, Beta.vec = c(beta), 
                                            alpha = alpha, gamma = gamma, ndraw = ndraw, 
                                            treat.method.list = treat.method.list, control.method.list = control.method.list, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,
                                            set = set, alpha.ratio.treat = alpha.ratio.treat,  tol = tol))
    }
  }
  return(res)
}

############### Inferential Method for Super Population (M2: using Berger and Boos (1994)'s approach) ######################################
### combine results using treated (* = 1) and control units (* = 0) separately ###
SP_method_combine_star1_and_star0 <- function(Z, Y, beta_vec = c(0.5, 0.6, 0.7, 0.8, 0.9), 
                                              alpha=0.05, gamma = 0.5, ndraw = 10^5, treat.method.list = list(name = "Stephenson", s = 6),
                                              control.method.list = list(name = "Stephenson", s = 6), score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,  
                                              alpha.ratio.treat = 0.5,  tol = 10^(-3), simul = TRUE ){
  ## using treated units (* = 1)
  ci.gen1 = ci_lower_quantile_generalize_SP(Z, Y, beta_vec, treat.method.list = treat.method.list, control.method.list = control.method.list,
                                            stat.null = stat.null, set = "treat", simul = simul,
                                            alpha=alpha/2, gamma = gamma, alpha.ratio.treat = alpha.ratio.treat, 
                                            score = score, Z.perm = Z.perm,
                                            ndraw = ndraw, nperm = nperm, tol = tol)
  
  ## using control units (* = 0)
  ci.gen2 = ci_lower_quantile_generalize_SP(Z, Y, beta_vec, treat.method.list = treat.method.list, control.method.list = control.method.list,
                                            stat.null = stat.null, set = "control", simul = simul,
                                            alpha=alpha/2, gamma = gamma, alpha.ratio.treat = alpha.ratio.treat, 
                                            score = score, Z.perm = Z.perm,
                                            ndraw = ndraw, nperm = nperm, tol = tol)
  return(pmax( ci.gen1, ci.gen2 ))
}

#### using all experimental units (* = empty) ####
SP_method_star_empty <- function(Z, Y, beta_vec = c(0.5, 0.6, 0.7, 0.8, 0.9), 
                                 alpha=0.05, gamma = 0.5, ndraw = 10^5, treat.method.list = list(name = "Stephenson", s = 6),
                                 control.method.list = list(name = "Stephenson", s = 6), score = NULL, stat.null = NULL, nperm = 10^4, Z.perm = NULL,  
                                 alpha.ratio.treat = 0.5,  tol = 10^(-3), simul = TRUE ){
  
  ci.empty <- ci_lower_quantile_generalize_SP(Z = Z, Y = Y, beta_vec = beta_vec, treat.method.list = treat.method.list, control.method.list = control.method.list,
                                              stat.null = stat.null, set = "all", simul = simul,
                                              alpha=alpha, gamma = gamma, alpha.ratio.treat = alpha.ratio.treat, 
                                              score = score, Z.perm = Z.perm,
                                              ndraw = ndraw, nperm = nperm, tol = tol)
  return(ci.empty)
}

#################### Helper functions for simulation ###########################
### Data generating process ###

# Y0 ~ N(0, rho^2), Y1 ~ N(0, 1-rho^2)
# The total variance of Y1 and Y0 is 1
# So that the results for different rho^2 may be more comparable

simulate_dt_mixture <- function(N, n, k, rho, mu = 2){
  
  # Generate Y(0) and Y(1) for whole population
  Y_0 = rnorm(N, mean = 0, sd = rho)
  Y_1 = rnorm(N, mean = mu, sd = sqrt(1-rho^2))
  Tau = Y_1 - Y_0
  
  # Sample experimental units
  sample_index = sample(1:N, n)
  y_0 = Y_0[sample_index]
  y_1 = Y_1[sample_index]
  tau = y_1 - y_0
  
  # n_t treated units and n_c control units
  n_t = n * k
  n_c = n - n_t
  
  Z = c(rep(1, n_t), rep(0, n_c))
  Z = sample(Z, n)
  Y_obs = Z * y_1 + (1-Z)*y_0
  
  dt_N = data.frame(Y_0 = Y_0, Y_1 = Y_1, Tau = Tau)
  dt_n = data.frame(Z = Z, y_0 = y_0, y_1 = y_1, Y_obs = Y_obs, tau = tau)
  dt_n = slice(dt_n, sample(1:n()))
  
  return(list(dt_N, dt_n))
}

## example
# dt = simulate_dt_mixture(100, 50, 0.5, 0.2)
# dt[1]: population data Nx3
# dt[2]: sample data nx5

### Function to compare different inferential methods for ITE quantiles ###                                                                       
compare <- function(dt, k_vec, simul = FALSE, alpha=0.05){
  
  N = dim(dt[[1]])[1]
  
  Z = dt[[2]]$Z
  Y_obs = dt[[2]]$Y_obs
  
  # J is length of k_vec
  J = length(k_vec)
  
  # Method M0 is the same for simul = TRUE or FALSE
  m0_s6 = method_caughey(Z, Y_obs, k_vec = k_vec, method.list = list(name = "Stephenson", s = 6), alpha=alpha)
  
  # Method M1 depends on simul if N!=n
  m1_s6 = method_combine(Z, Y_obs, N = N,  k_vec = k_vec, simul = simul, alpha=alpha)
  
  # Method M2 depends on simul
  m2_s6 = method_bergerboos(Z, Y_obs, N = N, k_vec = k_vec,  simul = simul, alpha=alpha)
 
  res = rbind(m0_s6, m1_s6, m2_s6)
  
  res$method = c(rep('M0-S6', J),
                 rep('M1-S6', J),
                 rep('M2-S6', J)))
  return(res)
}
