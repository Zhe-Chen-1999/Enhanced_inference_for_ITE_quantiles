# This script analyzes the effect of smoking on the blood cadmium level in Section 6.3 

# Load and randomly shuffle the cadmium dataset
data("cadmium")

# Randomly permute indices
set.seed(1)
ind.perm = sample(c(1:n))

Y = cadmium$cadmium[ind.perm]
Z = cadmium$z[ind.perm]
block = as.factor(cadmium$mset)[ind.perm]

n = length(Z)
n1 = sum(Z)
n0 = n - n1
method.list.all = list()
method.list.all[[1]] = list(name = "Wilcoxon")

# Function: Plot confidence intervals for effect quantiles
plot_quantile_lower <- function(ci.limit, k_start = NULL, main = NULL, 
                                quantiles = c(1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3), 
                                xlab_margin =5, left_ylab_margin = 3, right_ylab_margin = 3.5, 
                                margin=c(5.2,6,3,6), fontsize = 3, numbersize = 3){
  n = length(ci.limit)
  ticks = ceiling(quantiles*n)
  
  if (is.null(k_start)) {
    k_start = n - sum(!is.nan(ci.limit) & !is.infinite(ci.limit))
  }
  
  ylim = c(k_start-1, n+1)
  xlim = range(ci.limit[ci.limit > -Inf]) * 1.1
  par(mar=margin+0.1)
  plot(NA, ylab = "", xlab="", ylim = ylim, xlim = xlim, yaxt = "n",
       main = main, cex.lab=fontsize, cex.axis = numbersize, cex = 2)
  title(xlab = expression("lower" ~ "confidence" ~ "limit" ~ "for" ~ tau[(k)]),line = xlab_margin, cex.lab=fontsize)
  title(ylab = expression(k),line = left_ylab_margin, cex.lab=fontsize)
  
  for (k in 1:length(ci.limit)) {
    lines(c(max(ci.limit[k], min(xlim)-10*diff(xlim)), max(ci.limit) + 10*diff(xlim)), rep(k, 2), col = "grey")
  }
  
  points(ci.limit, c(1:length(ci.limit)), pch = 20)
  abline(v = 0, lty = 2)
  axis(side = 1, cex.axis=numbersize)
  axis(side = 2, at = ticks, cex.axis=numbersize)
  axis(side = 4, at = ticks, labels = paste0(floor(ticks*100/n), '%'), cex=fontsize, cex.axis=numbersize)
  mtext("quantile", side = 4, line = right_ylab_margin, cex=fontsize)
}

###################### Randomization Inference ########################
library(QIoT)

## Compute 90% simultaneous confidence intervals using Su & Li (2023)
LB_treat = ci_quantile_scre(Z,Y,block,method.list.all=method.list.all, opt.method = "Greedy", switch = FALSE, null.max=10^4, confidence = 0.9)$LB

pdf(paste0("smoke_all_rand.pdf"),
    width = 5.5,
    height = 6,
    pointsize = 6)

par(mar=c(5, 5, 2, 5), 
    xaxs = "i",
    yaxs = "i",
    cex.axis = 2,
    cex.lab = 2)

plot_quantile_lower(LB_treat,margin=c(5.8,5.1,3,5),
                    quantiles = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3),
                    fontsize = 3, numbersize = 2.5, k_start = 812)
dev.off()

## Apply Su and Li (2023)'s approach to data with switched treatment labels and changed outcome signs 
LB_treat_switch = ci_quantile_scre(Z,Y,block,method.list.all=method.list.all, opt.method = "Greedy", switch = TRUE, null.max=10^4, confidence = 0.9)$LB

pdf(paste0("smoke_all_rand_switch.pdf"),
    width = 5.5,
    height = 6,
    pointsize = 6)

par(mar=c(5, 5, 2, 5), 
    xaxs = "i",
    yaxs = "i",
    cex.axis = 3,
    cex.lab  = 3)

plot_quantile_lower(LB_treat_switch, margin=c(5.8,5.1,3,5), 
                    quantiles = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3),
                    fontsize = 3, numbersize = 2.5, k_start = 812)
dev.off()

## Compute 90% CIs from our improved method: combine 95% prediction intervals for treated and control units 
ci.treat = ci_quantile_scre(Z,Y,block,quantiles= (n0+1):n, method.list.all=method.list.all, opt.method = "Greedy", switch = FALSE, null.max=10^4, confidence = 0.95)$LB
ci.control = ci_quantile_scre(Z,Y,block,quantiles= (n1+1):n, method.list.all=method.list.all, opt.method = "Greedy", switch = TRUE, null.max=10^4, confidence = 0.95)$LB
ci.all = sort( c(ci.treat, ci.control) )

pdf(paste0("smoke_improved_CIs.pdf"),
    width = 5.5,
    height = 6,
    pointsize = 6)

par(mar=c(5, 5, 2, 5), 
    xaxs = "i",
    yaxs = "i",
    cex.axis = 3,
    cex.lab  = 3)

plot_quantile_lower(ci.all, quantiles = c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3), margin=c(5.8,5.1,3,5),
                    fontsize = 3, numbersize = 2.5, k_start = 812)
dev.off()

######################## Sensitivity Analysis ##############################

# Find the largest Gamma values s.t. the 90% prediction intervals for the 55%, 60%, 70%, 80%, 90% and 100% quantiles of ITEs among treated do not cover zero
quant.seq = seq(from = 0.6, to = 1, by = 0.1)
gammas = rep(NA, length(quant.seq))
enum = 1
gam = 1

for(quant in quant.seq){
  k = ceiling(quant * n1) + n0
  while(1){
    pval = pval_quantile_sen(Z, Y, block, k, c=0, gam=gam, method.list.all=method.list.all, ties = "upper", switch = FALSE)$upper
    if(pval > 0.1){
      gammas[enum] = max(gam - 0.1, 1)
      enum = enum + 1
      break
    }
    gam = gam + 0.1
  }
  print(quant)
}

print(gammas)

pval_quantile_sen(Z, Y, block, k=ceiling(0.55*n1) + n0, c=0, gam=1, method.list.all=method.list.all, ties = "upper", switch = FALSE)$upper

gamma.seq = c(1.0, gammas)


## 90% simultaneous CIs for effect quantiles of all units under various bounds on confounding strength using Su & Li(2023)'s approach
## with switched treatment labels and changed outcome signs

LB.all.original = matrix(NA, nrow = n, ncol = length(gamma.seq))

for(kk in 1:length(gamma.seq)){
  gam = gamma.seq[kk]
  LB.all.original[, kk] = ci_quantile_sen(Z, Y, block, gam=gam, method.list.all=method.list.all, switch = TRUE, confidence = 0.9)$LB
  print(kk)
}

pdf(paste0("smoke_sen_original_switch.pdf"),
    width = 4,
    height = 4,
    pointsize = 4)

par(mar=c(4, 5, 2, 5), 
    xaxs = "i",
    yaxs = "i",
    cex.axis = 2,
    cex.lab  = 2)

plot(LB.all.original[, 1], c(1:n), ylim=c(780, n+5), type = "p", col = 1, pch = 20, yaxt = "n", ylab = "k", xlim = c(-1.5, 1.5), xlab = expression("lower" ~ "confidence" ~ "limit" ~ "for" ~ tau[(k)]))

for(kk in 2:ncol(LB.all.original)){
  points( LB.all.original[, kk], c(1:n), pch = 20, col = kk )
}

abline(v=0, lty = 2)
ticks = ceiling(c(0.55,quant.seq)*n)
axis(side = 2, at = ticks, cex.axis=2) 
axis(side = 4, at = ticks, labels = paste0(floor(ticks*100/n), '%'), cex=2, cex.axis=2)
mtext("quantile", side = 4, line = 3, cex=2)
legend("bottomright", col = 1:7, pch = 1, cex = 2, 
       legend = c(expression(Gamma == 1.0),
                  expression(Gamma == 1.3),
                  expression(Gamma == 2.2),
                  expression(Gamma == 4.0),
                  expression(Gamma == 8.3),
                  expression(Gamma == 38.4)))

dev.off()

## 90% simultaneous CIs for effect quantiles of all units under various bounds on confounding strength using our improved approach
LB.all.improved = matrix(NA, nrow = n, ncol = length(gamma.seq))

for(kk in 1:length(gamma.seq)){
  gam = gamma.seq[kk]
  ci.treat = ci_quantile_sen(Z, Y, block, gam=gam, quantiles= (n0+1):n,
                  method.list.all=method.list.all,switch = FALSE, confidence = 0.95)$LB
  ci.control = ci_quantile_sen(Z, Y, block, gam=gam, quantiles= (n1+1):n,
                               method.list.all=method.list.all,switch = TRUE, confidence = 0.95)$LB
  LB.all.improved[, kk] = sort( c(ci.treat, ci.control) )
  print(kk)
}

pdf(paste0("smoke_sen_improved.pdf"),
    width = 4,
    height = 4,
    pointsize = 4)

par(mar=c(4, 5, 2, 5), 
    xaxs = "i",
    yaxs = "i",
    cex.axis = 2,
    cex.lab  = 2)

plot(LB.all.improved[, 1], c(1:n), ylim=c(780, n+5), type = "p", col = 1, pch = 20, yaxt = "n", ylab = "k", xlim = c(-1.5, 1.5), xlab = expression("lower" ~ "confidence" ~ "limit" ~ "for" ~ tau[(k)]))

for(kk in 2:ncol(LB.all.improved)){
  points( LB.all.improved[, kk], c(1:n), pch = 20, col = kk )
}

abline(v=0, lty = 2)
ticks = ceiling(c(0.55,quant.seq)*n)
axis(side = 2, at = ticks, cex.axis=2)
axis(side = 4, at = ticks, labels = paste0(floor(ticks*100/n), '%'), cex=2, cex.axis=2)
mtext("quantile", side = 4, line = 3, cex=2)
legend("bottomright", col = 1:7, pch = 1, cex = 2, 
       legend = c(expression(Gamma == 1.0),
                  expression(Gamma == 1.3),
                  expression(Gamma == 2.2),
                  expression(Gamma == 4.0),
                  expression(Gamma == 8.3),
                  expression(Gamma == 38.4)))

dev.off()

## 90% simultaneous prediction intervals for all effect quantiles among treated units
LB.all = matrix(NA, nrow = n1, ncol = length(gamma.seq))

for(kk in 1:length(gamma.seq)){
  gam = gamma.seq[kk]
  LB.all[, kk] = ci_quantile_sen(Z, Y, block, gam=gam, method.list.all=method.list.all, switch = FALSE, confidence = 0.9)$LB[(n0+1):n]
  print(kk)
}

# Gamma = 2.2: 30.7% smokers has a positive effect
sum(LB.all[,3]>=0)/nrow(LB.all)

pdf(paste0("smoke_sen_treat.pdf"),
    width = 4,
    height = 4,
    pointsize = 4)

par(mar=c(4, 5, 2, 5), 
    xaxs = "i",
    yaxs = "i",
    cex.axis = 2,
    cex.lab  = 2)

plot(LB.all[, 1], c(1:n1), ylim=c(260, n1+5), type = "p", col = 1, pch = 20, yaxt = "n", ylab = "k", xlim = c(-1.5, 1.5), 
     xlab = expression("lower" ~ "confidence" ~ "limit" ~ "for" ~ tau[(k)]))

for(kk in 2:ncol(LB.all)){
  points( LB.all[, kk], c(1:n1), pch = 20, col = kk )
}

abline(v=0, lty = 2)
ticks = ceiling(c(0.55,quant.seq)*n1)
axis(side = 2, at = ticks, cex.axis=2) 
axis(side = 4, at = ticks, labels = paste0(floor(ticks*100/n1), '%'), cex=2, cex.axis=2)
mtext("quantile", side = 4, line = 3, cex=2)
legend("bottomright", col = 1:7, pch = 1, cex = 2, 
       legend = c(expression(Gamma == 1.0),
                  expression(Gamma == 1.3),
                  expression(Gamma == 2.2),
                  expression(Gamma == 4.0),
                  expression(Gamma == 8.3),
                  expression(Gamma == 38.4)))
dev.off()
