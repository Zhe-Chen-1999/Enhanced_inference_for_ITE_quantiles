# Simulation studies in Appendix A3 of the supplement
# This script runs simulations to evaluate the performance of the proposed methods under different settings.

source('helper_functions.R')
library(dplyr)
library(forcats)
library(ggplot2)
library(gtools)
library(tidyr)
library(xtable)

#################################################################################################
################### A3.1. Comparison between the two improved methods ###########################
#################################################################################################

## Run simulations 
run_simu <- function(N, n, k, rho_square, mu = 2, MC = 1000){
  
  rho = sqrt(rho_square)
  all_res_simul = NULL
  all_res_point = NULL
  
  k_vec = c(ceiling(0.1*N), ceiling(0.2*N), ceiling(0.3*N), ceiling(0.4*N), 
            ceiling(0.5*N), ceiling(0.6*N), ceiling(0.7*N), ceiling(0.8*N), ceiling(0.9*N), N)
  
  for (i in 1:MC){
    cat(i, '\n')
    dt = simulate_dt_mixture(N = N, n = n, k = k, rho =rho, mu = mu)
    
    # Simultaneous inference
    res_temp_simul = compare(dt, k_vec = k_vec, simul = TRUE, alpha=0.1)
    res_temp_simul$simultaneous = rep(TRUE, dim(res_temp_simul)[1])
    all_res_simul = rbind(all_res_simul, res_temp_simul)
    
    # Pointwise inference
    res_temp_point = compare(dt, k_vec = k_vec, simul = FALSE, alpha=0.1)
    res_temp_point$simultaneous = rep(FALSE, dim(res_temp_point)[1])
    all_res_point = rbind(all_res_point, res_temp_point)
  }
  
  res = rbind(all_res_simul, all_res_point)
  res$rho_square = rho_square
  res$tau = rep(sort(dt[[1]]$Tau)[k_vec], 2*length(unique(res$method)))
  return(res)
}

res = NULL
for (N in c(100)) {
  for (rho_square in c(0.1, 0.9)) {
    cat(N, rho_square, '\n')
    r = run_simu(N, n = 100, k = 0.5, rho_square, mu = 2, MC = 500)
    res = rbind(res, r)
  }
}

write.csv(res, file = './merged_results_for_500_simulations.csv', row.names = FALSE)

############## measurement I: Comparing the non-informativeness of the confidence intervals #################################### 

all_res = read.csv('merged_results_for_500_simulations.csv')

info_measure <- all_res %>%
  group_by(simultaneous, method, k) %>%
  summarise(lower.limit.Inf = mean(lower == -Inf), .groups = 'drop') %>%
  pivot_wider(
    names_from = k,
    values_from = lower.limit.Inf,
    names_prefix = "",
    names_glue = "{k}%"
  )

xtable(info_measure)

############# measurement II: Comparing the median of the lower confidence bounds ############################################## 

median0 <- all_res %>%
  filter(N == 100) %>%
  mutate(Inference = ifelse(simultaneous==TRUE, "Simultaneous", "Pointwise"))%>%
  group_by(N, n, rho_square, method, k, Inference) %>%
  summarise(median = median(lower)) 

N = 100
k_vec = c(ceiling(0.9*N), ceiling(0.8*N), ceiling(0.7*N), ceiling(0.6*N), ceiling(0.5*N))
median_plot = median0 %>%
  mutate(median_new = ifelse(median == -Inf, -1, median)) %>% 
  filter(k %in% k_vec, method %in% c('M0-S6','M1-S6','M2-S6')) %>%
  mutate(lab = ifelse(k == k_vec[1], '90 Percentile', 
                  ifelse(k == k_vec[2], '80 Percentile', 
                         ifelse(k == k_vec[3], '70 Percentile', 
                                ifelse(k == k_vec[4], '60 Percentile', 
                                       ifelse(k == k_vec[5], '50 Percentile', k)))))) %>%
  mutate(method = ifelse(method == 'M2-S6', 'M2',
                         ifelse(method == 'M1-S6', "M1",
                                ifelse(method == "M0-S6", "M0", method))))%>%
  ggplot(aes(x=rho_square, y = median_new, color=c(method), linetype=Inference)) +
  scale_linetype_manual(values=c( "dotted", "solid"))+
  geom_line()+ 
  geom_point()+
  facet_wrap(~lab, nrow = 1) + 
  theme_bw(base_size = 22)+
  scale_color_viridis_d()+
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), limits = c(0.1, 0.9))+
  scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2), labels = c(-Inf, -0.5, 0, 0.5, 1, 1.5, 2))+
  xlab(expression(rho^2))+
  ylab(expression(atop("Median of 90%", "lower confidence limits")))+
  labs(color = "Method")+ 
  theme(legend.position="top", axis.title.y.left = element_text(size = 18), panel.spacing.x = unit(8, "mm"))

ggsave("FP_method_comparison_N100_n100_s6_90percentCI_Median_1row.pdf", plot = median_plot, 
       scale = 1, width = 14, height = 5, 
       units = "in")

#############################################################################################################################
########### A3.2. Examine the robustness of the second improved method with respect to the hyper-parameter gamma ############
#############################################################################################################################

## Function for inferring effect quantiles using the second method with different gamma values
BB.gamma <- function(dt, k_vec, simul = FALSE, alpha=0.05, gamma = 0.5){
  N = dim(dt[[1]])[1]
  Z = dt[[2]]$Z
  Y_obs = dt[[2]]$Y_obs
  J = length(k_vec)
  m2_s6 = method_bergerboos(Z, Y_obs, N = N, k_vec = k_vec, method.list = list(name = "Stephenson", s = 6), simul = simul, alpha=alpha, gamma = gamma, ndraw = 10^6)
  res = m2_s6
  res$method = c(rep('M2-S6', J))
  return(res)
}

## Run simulations and save results
run_simu <- function(N, n, k, rho_square, mu = 2, MC = 1){
  
  rho = sqrt(rho_square)
  all_res_simul = NULL
  all_res_point = NULL
  
  k_vec = c(ceiling(0.3*N), ceiling(0.4*N), ceiling(0.5*N), ceiling(0.6*N), ceiling(0.7*N), ceiling(0.8*N), ceiling(0.9*N), N)
  
  for (i in 1:MC){
    dt = simulate_dt_mixture(N = N, n = n, k = k, rho = rho, mu = mu)
    
    for (gamma in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) {
      
      # Simultaneous
      res_temp_simul = BB.gamma(dt, k_vec = k_vec, simul = TRUE, alpha=0.1, gamma = gamma)
      res_temp_simul$simultaneous = rep(TRUE, dim(res_temp_simul)[1])
      res_temp_simul$gamma = gamma
      all_res_simul = rbind(all_res_simul, res_temp_simul)
      
      # Pointwise
      res_temp_point = BB.gamma(dt, k_vec = k_vec, simul = FALSE, alpha=0.1, gamma = gamma)
      res_temp_point$simultaneous = rep(FALSE, dim(res_temp_point)[1])
      res_temp_point$gamma = gamma
      all_res_point = rbind(all_res_point, res_temp_point)
    }
  }
  res = rbind(all_res_simul, all_res_point)
  res$N = N
  res$n = n
  res$p = k
  res$rho_square = rho_square
  res$mu = mu
  res$tau = rep(sort(dt[[1]]$Tau)[k_vec], 2*length(unique(res$method)))
  return(res)
}

all_res = NULL
for (N in c(100, 200, 500)) {
  for (rho_square in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
    cat(N, rho_square, '\n')
    r = run_simu(N, n = 100, k = 0.5, rho_square, mu = 2, MC = 500)
    all_res = rbind(all_res, r)
  }
}

### Measurement: median of the 90% lower confidence limits across 500 simulations
median <- all_res %>%
  filter(N == 100) %>%
  mutate(inference = ifelse(simultaneous==TRUE, "Simultaneous", "Pointwise"), gamma = as.factor(gamma))%>%
  group_by(rho_square, inference, k, gamma) %>%
  summarise(median = median(lower)) 

k_vec = c(ceiling(0.9*N), ceiling(0.8*N), ceiling(0.7*N), ceiling(0.6*N), ceiling(0.5*N))

## Median plot
median_plot = median %>%
  mutate(median_new = ifelse(median == -Inf, -1, median)) %>% 
  filter(k %in% k_vec) %>%
  mutate(lab =  ifelse(k == k_vec[1], '90 Percentile', 
                       ifelse(k == k_vec[2], '80 Percentile', 
                              ifelse(k == k_vec[3], '70 Percentile', 
                                     ifelse(k == k_vec[4], '60 Percentile', 
                                            ifelse(k == k_vec[5], '50 Percentile', 
                                                   ifelse(k == k_vec[6], '40 Percentile', 
                                                          ifelse(k == k_vec[7], '30 Percentile', k)))))))) %>%
  ggplot(aes(x=rho_square, y = median_new, color=gamma)) +
  geom_line()+ 
  geom_point()+
  facet_grid(inference ~ lab) +
  theme_bw(base_size = 25)+
  scale_color_viridis_d()+
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9), limits = c(0.1, 0.9))+
  scale_y_continuous(breaks = c(-1, -0.5 ,  0 , 0.5,  1,  1.5,  2), labels = c(-Inf, -0.5 ,  0 ,  0.5,  1, 1.5, 2))+
  xlab(expression(rho^2))+
  ylab("Median of 90% lower confidence limits")+ 
  labs(color = expression(gamma))+
  theme(legend.position = c(.85, .2), panel.border = element_rect(colour = "black", fill=NA))+ 
  theme(legend.position="top")+
  theme(panel.spacing.x = unit(8, "mm"), axis.title.y.left = element_text(size = 18))