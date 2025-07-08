# This script evaluates the effectiveness of professional development in Section 6.2

source('helper_functions.R')

# Load and randomly shuffle the teachers dataset
data( electric_teachers )
n = nrow( electric_teachers )
dat = electric_teachers[sample(n), ]

treat = dat$TxAny
outcome = dat$gain
n1 = sum(treat)
n0 = n - n1

# Function: Plot confidence intervals for effect quantiles
plot_quantile_lower <- function(ci.limit, k_start = NULL, main = NULL, 
                                x_custom = FALSE, x_custom_range = c(-20, 20),
                                quantiles = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3),
                                line = 3, fontsize = 1.5, numbersize = 1.2){
  n = length(ci.limit)
  ticks = ceiling(quantiles*n)
  if (is.null(k_start)) {
    k_start = n - sum(!is.nan(ci.limit) & !is.infinite(ci.limit))
  }
  ylim = c(k_start-1, n+1)
  if(x_custom)
    xlim = x_custom_range
  else
    xlim = range(ci.limit[ci.limit > -Inf]) * 1.1
  par(mar=c(4,4,1,5)+0.1)
  plot(NA, ylab = "k", xlab = expression("lower" ~ "confidence" ~ "limit" ~ "for" ~ tau[(k)]), ylim = ylim, xlim = xlim, yaxt = "n",
       main = main, cex.lab=fontsize, cex.axis = numbersize)
  for (k in 1:length(ci.limit)) {
    lines(c(max(ci.limit[k], min(xlim)-10*diff(xlim)), max(ci.limit) + 10*diff(xlim)), rep(k, 2), col = "grey")
  }
  points(ci.limit, c(1:length(ci.limit)), pch = 20)
  abline(v = 0, lty = 2)
  axis(side = 2, at = ticks, cex.axis=numbersize)
  axis(side = 4, at = ticks, labels = paste0(floor(ticks*100/n), '%'), cex=fontsize )
  mtext("quantile", side = 4, line = line, cex=fontsize)
}


######################################################################
################## Randomization Inference in CRE ####################
######################################################################

# Method in Caughey et al.(2023)
m0 = method_caughey(Z=treat, Y=outcome, k_vec = 1:n, method.list = list(name = "Stephenson", s = 6), nperm = 10^4,alpha = 0.1)
pdf(file="m0.pdf")
plot_quantile_lower(m0$lower, k_start = 65, main = NULL)
dev.off()

# The improved method outlined in Section 3.2
m1 = method_combine(Z=treat, Y=outcome, N = n, alpha = 0.1,
                    control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                    nperm = 10^4, k_vec = 1:n, simul = FALSE)
pdf(file="m1.pdf")
plot_quantile_lower(m1$lower, k_start = 65, main = NULL)
dev.off()

# The improved method outlined in Section 3.3
m2_pointwise = method_bergerboos(Z=treat, Y=outcome, N = n, alpha = 0.1,
                                 control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                                 nperm = 10^4, k_vec = 1:n, simul = FALSE)

m2_simul = method_bergerboos(Z=treat, Y=outcome, N = n, k_vec = ceiling(seq(0.5, 1, 0.1)*n), alpha = 0.1,
                             control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                             nperm = 10^4, simul = TRUE)

ci = rep(-Inf, n)

for (j in 2: length(m2_simul$k)) {
  ci[(m2_simul$k[j-1]):(m2_simul$k[j]-1)] = m2_simul$lower[j-1]
}

ci[m2_simul$k[length(m2_simul$k)]:n]=m2_simul$lower[length(m2_simul$k)]

pdf(file="m2_with_simul.pdf")
plot_quantile_lower(m2_pointwise$lower, k_start = 65, main = NULL)
points(ci, y = 1:n, col = "red")

for (i in 1:dim(m2_simul)[1]) {
  lines(c(m2_simul$lower[i],20), rep(m2_simul$k[i], 2), col = 'red')
}

dev.off()

################################################################################
####### Randomization Inference in Sampling-based Randomized Experiments #######
################################################################################

################# Inference for Finite Population #################
fp_all_res = NULL

for(N in c(1000, 5000, 10000)){
  for (simul in c(TRUE, FALSE)) {
    res = NULL
    
    ## star = empty
    ci.empty =  method_combine(Z=treat, Y=outcome, N = N, alpha=0.1, ndraw = 10^6,
                               control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6),
                               simul = simul,  k_vec = ceiling(c(0.5, 0.6, 0.7, 0.8, 0.9)*N), nperm = 10^4)
    ci.empty$star = "empty"
    res = rbind(res, ci.empty) 
    
    ## combine star = 1 and star = 0
    ci.combine = method_bergerboos(Z=treat, Y=outcome, N = N, alpha = 0.1, ndraw = 10^6,
                                   control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                                   nperm = 10^4, k_vec = ceiling(c(0.5, 0.6, 0.7, 0.8, 0.9)*N), simul = simul)
    ci.combine$star = "combine.star1.star0"
    res = rbind(res, ci.combine)
    
    res$inference = ifelse(simul==TRUE, "Simultaneous", "Pointwise")
    res$N = N
    fp_all_res = rbind(fp_all_res, res)
  }
}

fp_all_res$beta = fp_all_res$k/fp_all_res$N

############## Inference for Super Population #################
sp_all_res = NULL

for (simul in c(TRUE, FALSE)) {
  res = NULL
  
  ## star = empty
  ci.empty =  SP_method_star_empty(Z=treat, Y=outcome, beta_vec = c(0.5, 0.6, 0.7, 0.8, 0.9), 
                                   alpha=0.1, ndraw = 10^6, treat.method.list = list(name = "Stephenson", s = 6),
                                   control.method.list = list(name = "Stephenson", s = 6), nperm = 10^4, 
                                   simul = simul )
  ci.empty$star = "empty"
  res = rbind(res, ci.empty) 
  
  ## combine star = 1 and star = 0
  ci.combine = SP_method_combine_star1_and_star0(Z=treat, Y=outcome, beta_vec = c(0.5, 0.6, 0.7, 0.8, 0.9), 
                                                 alpha=0.1, ndraw = 10^6, treat.method.list = list(name = "Stephenson", s = 6),
                                                 control.method.list = list(name = "Stephenson", s = 6), nperm = 10^4, 
                                                 simul = simul )
  ci.combine$star = "combine.star1.star0"
  res = rbind(res, ci.combine)
  
  res$inference = ifelse(simul==TRUE, "Simultaneous", "Pointwise")
  sp_all_res = rbind(sp_all_res, res)
}

sp_all_res$N = Inf


################################################################################
########### Two-sided confidence intervals for the education experiment ########
################################################################################

# Function for drawing two-sided confidence intervals for the effect quantiles
plot_two_sided_CIs <- function(lb, ub, k_start = NULL, main = NULL, 
                               x_custom = FALSE, x_custom_range = c(-20, 20),
                               quantiles = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
                               line = 3, fontsize = 1.5, numbersize = 1.2){
  n = length(lb)
  ticks = ceiling(quantiles*n)
  
  if (is.null(k_start)){
    k_start = 1
  }
  
  ylim = c(k_start-1, n+1)
  
  if(x_custom)
    xlim = x_custom_range
  else
    xlim = range(c(lb[lb > -Inf], ub[ub<Inf])) * 1.1
  
  par(mar=c(4,4,1,5)+0.1)
  
  plot(NA, ylab = "k", xlab = expression("Two-sided" ~ "confidence" ~ "interval" ~ "for" ~ tau[(k)]), ylim = ylim, xlim = xlim, yaxt = "n",
       main = main, cex.lab=fontsize, cex.axis = numbersize)
  
  for (k in 1:length(lb)) {
    lines(c(max(lb[k], min(xlim)-10*diff(xlim)), min(ub[k], max(xlim) + 10*diff(xlim))), rep(k, 2), col = "grey")
  }
  
  points(lb, c(1:length(lb)), pch = 20)
  points(ub, c(1:length(ub)), pch = 20)
  abline(v = 0, lty = 2)
  axis(side = 2, at = ticks, cex.axis=numbersize)
  axis(side = 4, at = ticks, labels = paste0(floor(ticks*100/n), '%'), cex=fontsize )
  mtext("quantile", side = 4, line = line, cex=fontsize)
}

# Method in Caughey et al.(2023)
m0_lb = method_caughey(Z=treat, Y=outcome, k_vec = 1:n, method.list = list(name = "Stephenson", s = 6), nperm = 10^4,alpha = 0.1/2)
m0_ub = method_caughey(Z=1-treat, Y=outcome, k_vec = 1:n, method.list = list(name = "Stephenson", s = 6), nperm = 10^4,alpha = 0.1/2)

pdf(file="m0_two_sided_CIs.pdf")
plot_two_sided_CIs(m0_lb$lower, -m0_ub$lower[n:1], k_start = 1, main = NULL)
dev.off()


# Method outlined in Section 3.2
m1_lb = method_combine(Z=treat, Y=outcome, N = n, alpha = 0.1/2,
                       control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                       nperm = 10^4, k_vec = 1:n, simul = FALSE)


m1_ub = method_combine(Z=1-treat, Y=outcome, N = n, alpha = 0.1/2,
                       control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                       nperm = 10^4, k_vec = 1:n, simul = FALSE)

lb = m1_lb$lower
ub = -m1_ub$lower
ub = ub[233:1]
cbind(lb,ub)

pdf(file="m1_two_sided_CIs.pdf")
plot_two_sided_CIs(lb, ub, k_start = 1, main = NULL)
dev.off()

# Method outlined in Section 3.3
m2_pointwise_lb = method_bergerboos(Z=treat, Y=outcome, N = n, alpha = 0.1/2,
                                    control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                                    nperm = 10^4, k_vec = 1:n, simul = FALSE)

m2_pointwise_ub = method_bergerboos(Z=1-treat, Y=outcome, N = n, alpha = 0.1/2,
                                    control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                                    nperm = 10^4, k_vec = 1:n, simul = FALSE)


m2_simul_lb = method_bergerboos(Z=treat, Y=outcome, N = n, k_vec = ceiling(seq(0.1, 1, 0.1)*n), alpha = 0.1/2,
                                control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                                nperm = 10^4, simul = TRUE)

m2_simul_ub = method_bergerboos(Z=1-treat, Y=outcome, N = n, k_vec = n+1-ceiling(seq(0.1, 1, 0.1)*n), alpha = 0.1/2,
                                control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), 
                                nperm = 10^4, simul = TRUE)

cbind(m2_simul_lb$lower, -m2_simul_ub$lower[length(m2_simul_ub$lower):1])


pdf(file="m2_two_sided_CIs.pdf")
plot_two_sided_CIs(m2_pointwise_lb$lower, -m2_pointwise_ub$lower[n:1], k_start = 1, main = NULL)

xlim = range(c(m2_pointwise_lb$lower[m2_pointwise_lb$lower > -Inf], -m2_pointwise_ub$lower[-m2_pointwise_ub$lower<Inf])) * 1.1

for (i in 1:dim(m2_simul_lb)[1]) {
  lines(c(max(m2_simul_lb$lower[i], min(xlim)-10*diff(xlim)), 
          min(-m2_simul_ub$lower[length(m2_simul_ub$lower):1][i], max(xlim) + 10*diff(xlim))), 
        rep(m2_simul_lb$k[i], 2), col = 'red')
}

dev.off()


################################################################################
################################ Method Illustration ###########################
################################################################################

k = ceiling(233*0.75) # 75% quantile: k = 175

# Rescale the outcome
outcome[treat==1] = mean(outcome[treat==1]) + scale(outcome[treat==1]) * sd(outcome[treat==1]) * 3

## Method in Caughey et al.(2023) with original labels
m0_original = method_caughey(Z=treat, Y=outcome, k_vec = 1:n, method.list = list(name = "Stephenson", s = 6), nperm = 10^5, alpha = 0.1)$lower

# Lower bound and corresponding quantile among treated
m0_original[k]
k - n0
sum(m0_original>0)/n
sum(m0_original>0)/n1

## Method in Caughey et al.(2023) with switched labels
m0_switched = method_caughey(Z=1-treat, Y=-outcome, k_vec = 1:n, method.list = list(name = "Stephenson", s = 6), nperm = 10^5, alpha = 0.1)$lower

# Lower bound and corresponding quantile among control
m0_switched[k] 
k - n1
sum(m0_switched>0)/n
sum(m0_switched>0)/n0

## Improved method outlined in Section 3.3
ci.treat.BB = ci_lower_quantile_generalize(Z=treat, Y=outcome, N=n, k_vec = c(1:n), treat.method.list = list(name = "Stephenson", s = 6), set = "treat", simul = FALSE, alpha=0.1, gamma = 0.5, alpha.ratio.treat = 0.5, tol = 10^(-3))
ci.treat.BB[k,]
sum(ci.treat.BB$lower>0)/n

ci.treat.BB2 = ci_lower_quantile_generalize(Z=treat, Y=outcome, N=n, k_vec = c(1:n), treat.method.list = list(name = "Stephenson", s = 6), set = "treat", simul = FALSE, alpha=0.05, gamma = 0.5, alpha.ratio.treat = 0.5, tol = 10^(-3))
ci.control.BB2 = ci_lower_quantile_generalize(Z=treat, Y=outcome, N=n, k_vec = c(1:n), treat.method.list = list(name = "Stephenson", s = 6), set = "treat", simul = FALSE, alpha=0.05, gamma = 0.5, alpha.ratio.treat = 0.5, tol = 10^(-3))
sum(ci.treat.BB2$lower>0)/n
sum(ci.control.BB2$lower>0)/n

alpha.target = 0.1
correct.treat = threshold_correct_hypergeom(N = n, n = n1, K.vec = k, alpha = 0.5 * alpha.target, ndraw = 10^5)
correct.treat$kk.vec
m0_original[n0+correct.treat$kk.vec]
method_caughey(Z=treat, Y=outcome, k_vec = 1:n, method.list = list(name = "Stephenson", s = 6), nperm = 10^5, alpha = alpha.target-correct.treat$prob)$lower[n0+correct.treat$kk.vec]

correct.control = threshold_correct_hypergeom(N = n, n = n0, K.vec = k, alpha = 0.5 * alpha.target, ndraw = 10^5)
correct.control$kk.vec
m0_switched[n1+correct.control$kk.vec]
method_caughey(Z=1-treat, Y=-outcome, k_vec = 1:n, method.list = list(name = "Stephenson", s = 6), nperm = 10^5, alpha =alpha.target - correct.control$prob)$lower[n1+correct.control$kk.vec]

## Improved method outlined in Section 3.2
ci.treat = method_caughey(Z=treat, Y=outcome, k_vec = 1:n, method.list = list(name = "Stephenson", s = 6), nperm = 10^5, alpha = 0.1/2)$lower[(n0+1):n]
ci.control = method_caughey(Z=1-treat, Y=-outcome, k_vec = 1:n, method.list = list(name = "Stephenson", s = 6), nperm = 10^5, alpha = 0.1/2)$lower[(n1+1):n]
ci.cb = sort(c(ci.treat, ci.control))

quantile(ci.treat, probs = c(0:10)/10)
quantile(ci.control, probs = c(0:10)/10)

sum(ci.cb > 0)/n
ci.cb[k]
ci.treat[k-n0]
ci.control[k-n1]

cbind(ci.treat, m0_original[n0+c(1:n1)])

j.vec = c(max(1, k+1-(n0+1)) : min(n1+1, k+1-1))
cbind(j.vec, k+1-j.vec)

ci.combine = apply(cbind(ci.treat[j.vec], ci.control[k+1-j.vec]), 1, min, na.rm = TRUE)
j.vec[which(ci.combine == max(ci.combine))]

ci.treat[j.vec[20]]
ci.control[k+1-j.vec[20]]

sort(c(ci.treat, ci.control))[170]

ci.search = matrix(-Inf, nrow = n-k+1+1, ncol = 2)

ci.search[2:(n-k+2), 1] = ci.treat[ n1+1-c(1:(n-k+1)) ]
ci.search[1:(n-k+1), 2] = ci.control[ n0+1-rev(c(1:(n-k+1))) ]

sum(ci.treat>=10)
sum(ci.control>=10)

cbind( n1+2 - c(1:(n-k+2)), ci.search, n0+2 - rev(c(1:(n-k+2))) )

which.max(apply(ci.search, 1, min))

method_combine(Z=treat, Y=outcome, N = n, alpha = 0.1,
               control.method.list = list(name = "Stephenson", s = 6), treat.method.list = list(name = "Stephenson", s = 6), nperm = 10^4, k_vec = 1:n, simul = TRUE)

