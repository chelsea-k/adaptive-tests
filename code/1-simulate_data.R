## Script for generating item responses and corresponding risk class 
##
## The item responses and risk outcome class are generated using the 
## Gaussian copula factor model. First, a MVN random variable is drawn
## according to a Gaussian factor model. Then the copula is achieved
## using the inversion method to get the proper marginals. 

# Directory for saving simulated data
save_dir = "simulated_data/"
dir.create(save_dir, recursive = TRUE)

# Set parameters for simulation
ntrain = 1000
ntest = 1000
n = ntrain + ntest
k = 2  # number of factors 
p = 5  # number of items 
n_resp = sample(2:5, p, replace=TRUE) # number of responses per item
n_resp = c(n_resp, 7) # last covariate will represent Age (7 levels: 12-18)
prop_at_risk = 0.15   # proportion in the "at-risk" group

# Get cutpoints for marginals (we assume uniform marginals here)
cuts = apply(as.matrix(n_resp), MARGIN = 1, function(x) {qnorm(c(0:x)/x)})

# Generate the factor loadings matrix Lambda. The entries of lambda
# are chosen to create a nice distribution of correlations between
# the items. We also add some 0's in there, being sure that we 
# don't have all zeros in the same row.

# We also draw factor loadings for the response y, which we treat as 
# another variable in the Gaussian copula factor model.

Lambda = matrix(sample(c(-5:-1,1:5)/sqrt(k),(p+1)*k,replace=TRUE),p+1,k)
zeros = matrix(sample(1:(p+1), floor((p+1)/10)*k), ncol=k)
for (j in 1:k){
  Lambda[zeros[,j],j] = 0
}

Lambda_y = matrix(1/sqrt(k),1,k)

# Draw factor scores Eta, compute Z = Lambda*Eta + Epsilon, scale Z by sd
Eta = matrix(rnorm(k*n), nrow=k, ncol=n)
Z = Lambda %*% Eta + rnorm((p+1)*n)
u = 1/(1+rowSums(Lambda^2))
Z = diag(sqrt(u))%*%Z # z_ij = z_ij / sqrt [1 + sum_{h=1}^k lambda_{jh}^2]

# Get X item responses based on cutoffs from marginals
X = t(apply(matrix(1:(p+1)), 1, function(j) cut(Z[j,], cuts[[j]],labels=FALSE)))

# Shift Age covariate to be from 12-18, not 1-7
X[p+1,] = X[p+1,] + 11 

# Compute Z_y and y
Z_y = Lambda_y %*% Eta + rnorm(n)
Z_y = Z_y/sqrt(1+sum(Lambda_y^2))
mu_Zy = qnorm(prop_at_risk)
Z_y = Z_y + mu_Zy 
y = 1*(Z_y>0)

# Store the final simulated data 
sim_data = cbind(y=t(y), data.frame(t(X)))
colnames(sim_data)[p+2] = "Age"
write.csv(sim_data, paste0(save_dir, "item_response_data_all.csv"), row.names=F)
write.csv(sim_data[1:ntrain,], paste0(save_dir, "item_response_data_train.csv"), row.names=F)
write.csv(sim_data[(ntrain+1):n,], paste0(save_dir, "item_response_data_test.csv"), row.names=F)

# Uncomment below to save more info about the generated data for other use
# copula.sim = list(sim_data = sim_data, Z=Z, Z_y = Z_y, factor.scores=Eta,
#                   Lambda=Lambda, Lambda_y = Lambda_y,
#                   mu_Zy = mu_Zy, k=k, p=p,
#                   ntrain=ntrain, ntest=ntest,
#                   prop_at_risk = prop_at_risk,
#                   cuts = cuts)
# save(copula.sim, file=paste0(save_dir, "copula.sim.txt"), ascii = TRUE)


