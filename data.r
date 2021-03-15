countn = c(15,1,
           6,1,2,
           6,6,
           7,2,3,2,
           16,9,3,3,1,
           57,38,17,2,2,
           119,81,45,6,1,1,
           173,118,57,16,3,1,
           136,103,50,13,6,1,1,
           54,51,32,5,1,1,
           13,15,12,3,1,
           4,3,1,
           1,1)

record = c(1,0,1,1,
           2,0,2,1,2,2,
           3,0,3,1,
           4,0,4,1,4,2,4,4,
           5,0,5,1,5,2,5,3,5,4,
           6,0,6,1,6,2,6,3,6,4,
           7,0,7,1,7,2,7,3,7,4,7,7,
           8,0,8,1,8,2,8,3,8,4,8,8,
           9,0,9,1,9,2,9,3,9,4,9,5,9,6,
           10,0,10,1,10,2,10,3,10,4,10,9,
           11,0,11,1,11,2,11,3,11,4,
           12,1,12,2,12,3,
           13,2,13,7)

record_m = matrix(record, ncol = 2, byrow = T)
data = matrix(ncol = 2)
for (i in 1:length(countn)) {
  data = rbind(data, matrix(rep(record_m[i,], countn[i]), ncol = 2, byrow = T))
}
data = cbind(data[-1,], Litter = 1:1328) %>% as.data.frame()
colnames(data)[1:2] = c("Nimplants","Ndead")

# MCEM
set.seed(232)
litter = 1328
mu = 1
sigma = runif(1)
out = matrix(NA, ncol = 1, nrow = litter)

for (i in 1:100) {
  # E-step
  ## Metropolis-Hastings algorithm
  alpha_initial = rnorm(litter, mean=0, sd=sigma)
  for (k in 1:1500) { # set a relatively large number for convergence
    alpha_update = rnorm(litter, mean=0, sd=sigma)
    
    for (j in 1:litter) {
      p_new = expit(mu+alpha_update[j])^data[j,2]*
        (1-expit(mu+alpha_update[j]))^(data[j,1]-data[j,2])
      p_old = expit(mu+alpha_initial[j])^data[j,2]*
        (1-expit(mu+alpha_initial[j]))^(data[j,1]-data[j,2])
      alpha = min(1, p_new/p_old)
      
      u = runif(1)
      if (u < alpha) {alpha_initial[j] = alpha_update[j]}
    }
    out = cbind(out, alpha_initial)
  }
  out = out[,-(1:1000)] # drop the non-converged values
  
  fr = function(x) {
    -mean(apply(out,2, function(i) {sum(x*data[,2]+data[,2]*i) -
        sum(data[,1]*log(1+exp(x+i)))
    }))}
  
  # M-step
  final_mu = NULL
  final_sigma = NULL
    final = optim(mu, fr, lower = -5, upper = 5, method = "L-BFGS-B")
    mu = final$par
    sigma = sqrt(mean(out^2))
    
    print(paste(i, mu, sigma))
    final_mu[i] = mu
    final_sigma[i] = sigma
}

final = cbind(final_mu, final_sigma)

# save(final,file="/Users/xuchenghuiyun/Desktop/final.Rdata")
# load(file="/Users/xuchenghuiyun/Desktop/final.Rdata")
