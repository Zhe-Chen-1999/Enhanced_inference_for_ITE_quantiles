# This script performs inference for the total effect among treated units on a binary outcome, using the electric teacher dataset as an example.

library(RIQITE)
library(dplyr)

# Load and randomly shuffle the teachers dataset
data( electric_teachers )
n = nrow( electric_teachers ) 
dat = electric_teachers[sample(n), ]

df <- data.frame(
  ID = 1:n,
  Z = dat$TxAny, # 1 = treated, 0 = control
  Y = ifelse(dat$gain >= 10, 1, 0) # Dichotomize outcome
)

n1 = sum(df$Z) 
n0 = n - n1

################# Rigdon and Hudgens (2015)'s Method ###########################
 
AE.t.upper.CI = function(data,level = 0.1) {
  n = sum(data)
  m = sum(data[1, ])
  f = double()
  A1 = -data[1, 2]:data[1, 1]
  for (i in (1:length(A1))) {
    j = data
    j[1, 1] = j[1, 1] - A1[i]
    j[1, 2] = j[1, 2] + A1[i]
    f[i] = fisher.test(j)$p.value
  }
  lower1 = min(A1[f >= (level)])
  upper1 = max(A1[f >= (level)])
  output.all = list(lower=lower1,upper=upper1)
  return(output.all)
}

table = matrix(c(sum(df$Z*(df$Y)),sum(df$Z*(1-df$Y)),sum((1-df$Z)*(df$Y)),sum((1-df$Z)*(1-df$Y))), 2, 2, byrow=TRUE)
hudgens_res = AE.t.upper.CI(table) 

######################## Our method #######################

# Compute the minimum value of the test statistic under the null H_{c}: \sum_{Z_i*tau_i} = c 
min_stat <- function(df, c) {
  
  df <- df %>%
    mutate(Y_0_t0 = ifelse(Z == 0, Y, 0), # hypothesized Y(0) for ranking
           Y_0_t1 = ifelse(Z == 0, Y, 1),
           Rank_t0 = NA, Rank_t1 = NA, Rank_increase = NA) # Initialize rank column
  
  
  # Step 1: Compute Rank increase for treated units
  for (i in which(df$Z == 1)) {
    # Keep only the selected treated unit and all control units, preserving original order
    subset_indices <- sort(c(i, which(df$Z == 0)))  # Maintain original order
    subset_df <- df[subset_indices, ] 
    
    # Compute rank within this subset
    df$Rank_t0[i] <- rank(subset_df$Y_0_t0, ties.method = "first")[which(subset_indices == i)]
    df$Rank_t1[i] <- rank(subset_df$Y_0_t1, ties.method = "first")[which(subset_indices == i)]
    df$Rank_increase[i] <- df$Rank_t1[i] - df$Rank_t0[i]
  }
  
  # Step 2: Select (\sum_{Z_i*Y_i(0)}) treated units with the smallest Rank_increase
  num_ones_t = sum(df$Z*df$Y)-c 
  
  treated_smallest <- df %>%
    filter(Z == 1) %>%
    arrange(Rank_increase) %>%
    slice_head(n = num_ones_t) %>% 
    pull(ID)  # Extract IDs of selected units
  
  # Step 3: Assign Y(0) = 1 for selected treated units
  df <- df %>%
    mutate(Y0_minstat = ifelse(Z == 0, Y, ifelse(ID %in% treated_smallest, 1, 0)))
  
  stat.min = sum( rank(df$Y0_minstat, ties.method = "first")[df$Z==1] )
  return(stat.min)
}

# Score function of the ranks
rank_score <- function(n){
  score = c(1:n)
  return(score)
}

# Compute p-value for testing the null H_{c}: \sum_{Z_i*tau_i} = c
pval_Hc <- function(df, c,
                    score = NULL, stat.null = NULL,
                    nperm = 10^5, Z.perm = NULL){
  n = nrow(df)
  m = sum(df$Z)
  
  # get score if score is null #
  if(is.null(score)){
    score = rank_score( n )
  }
  
  # emp null dist #
  if(is.null(stat.null)){
    stat.null = null_dist(n, m, score = score, nperm = nperm, Z.perm = Z.perm )
  }
  
  # min stat value under H_{c} #
  stat.min = min_stat(df, c)
  
  # p-value #
  pval = mean( stat.null >= stat.min )
  
  return(pval)
}

# Compute the (1-alpha) lower confidence limit for the total effect among treated units
upper_CI <- function(alpha, df,  
                     score = NULL, stat.null = NULL, nperm = 10^5, Z.perm = NULL){
  
  n = nrow(df)
  df$ID = 1:n
  m = sum(df$Z)
  
  # binary search of the threshold c where p-val = alpha
  c.max = m-sum(df$Z*(1-df$Y))
  c.min = -m+sum(df$Z*df$Y)
  
  thres = NULL
  
  # emp null dist #
  if(is.null(stat.null)){
    stat.null = null_dist(n, m, score = score, nperm = nperm, Z.perm = Z.perm )
  }
  
  pval.upper = pval_Hc(df, c.max, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm)
  if(pval.upper <= alpha){
    thres = c.max
  }
  
  pval.lower = pval_Hc(df, c.min, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm)
  if(pval.lower > alpha){
    thres = c.min
  }
  
  while(is.null(thres)){
    c.mid = floor((c.max+c.min)/2)
    pval.mid = pval_Hc(df, c.mid, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm)
    
    if(pval.mid <= alpha){
      c.min = c.mid
    }else{
      c.max = c.mid
    }
    
    if( c.max - c.min <= 1 ){
      thres = c.max 
    }
  }
  
  pval.thres.max = pval_Hc(df, c.max, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm)
  pval.thres.min = pval_Hc(df, c.min, score = score, stat.null = stat.null, nperm = nperm, Z.perm = Z.perm)
  
  if(pval.thres.max >= 0.05 & pval.thres.min < 0.05){
    thres = c.max
  }else if(pval.thres.min >= 0.05){
    thres = c.min
  }else{
    thres = NA
  }
  
  return(c(ci.lower = thres, ci.upper = m-sum(df$Z*(1-df$Y))))
}

# Reshuffle the order of units 100 times and conduct inference based on each of the reshuffled data
res = NULL
for (i in 1:100) {
  cat(i, "\n")
  dat = df[sample(n), ]
  lower.limit = upper_CI(alpha = 0.05, dat)[1]
  dat$Y = 1 - dat$Y
  upper.limit = -upper_CI(alpha = 0.05, dat)[1]
  ci = c(lower.limit, upper.limit)
  res = rbind(res, ci)
}

result_vector <- na.omit(res)

# Summary statistics for the length of the prediction intervals
summary(result_vector[,2]-result_vector[,1])

# Histogram for the lower prediction bounds of the 90% prediction intervals
lower_result_vector <- na.omit(res[,1])
summary(lower_result_vector)

hist(lower_result_vector, 
     breaks = 30, 
     xlim = c(min(lower_result_vector)-10, max(c(lower_result_vector, hudgens_res$lower))+10),
     col = "lightblue", 
     main = NULL, 
     xlab = NULL, 
     ylab = "Frequency",
     border = "white",
     probability = TRUE)  

# Add density curve
lines(density(lower_result_vector), col = "blue", lwd = 2)

# Add a vertical red line for the threshold
abline(v = hudgens_res$lower, col = "red", lwd = 4, lty = 2)

legend("topright", 
       legend = c(paste("Rigdon and Hudgens's approach (2015)"), "Median of rank-based approach"), 
       col = c("red", "orange"), 
       lwd = 2, 
       lty = c(2, 2), 
       bg = "white")

# Histogram for the upper prediction bounds of the 90% prediction intervals
upper_result_vector <- na.omit(res[,2])
hist(upper_result_vector, 
     breaks = 30, 
     xlim = c(min(upper_result_vector)-10, max(c(upper_result_vector, hudgens_res$upper))+10),
     col = "lightblue", 
     main = NULL, 
     xlab = NULL, 
     ylab = "Frequency",
     border = "white",
     probability = TRUE) 

# Add density curve
lines(density(upper_result_vector), col = "blue", lwd = 2)

# Add a vertical red line for the threshold
abline(v = hudgens_res$upper, col = "red", lwd = 4, lty = 2)

legend("topright", 
       legend = c(paste("Rigdon and Hudgens's approach (2015)"), "Median of rank-based approach"), 
       col = c("red", "orange"), 
       lwd = 2, 
       lty = c(2, 2), 
       bg = "white")
