################################################################################################
#Simulate Datasets
################################################################################################

rm(list=ls())
set.seed(10032)
#define a function to draw samples from a Dirichlet distribution
rDirichlet <- function(alpha_vec){
  num <- length(alpha_vec)
  temp <- rgamma(num, shape = alpha_vec, rate = 1)
  return(temp / sum(temp))
}

n <- 50     #number of samples
n1 <- 50    #number of controls
n2 <- 50    #number of cases
m <- 600   #number of CpG sites
K <- 3       #underlying cell type number

#methylation profiles
mu <- matrix(rbeta(m*K,3,6),m,K)

#number of covariates
p <- 1

#covariates / phenotype
# X <- rbind(c(rep(0, n1),rep(1, n2)), runif(n, min=20, max=50))
# X <- matrix(c(rep(0, n1),rep(1, n2)),p,n)
X <- matrix(rnorm(p*n),p,n)
# X <- t(scale(t(X)))

#set risk-CpG sites under each cell type for each phenotype
# beta <- array(0, dim=c(m,K,p))
beta <- array(rnorm(m*K*p), dim=c(m,K,p))

#control vs case
# beta[1:300, 1, 1] <- 0.1   #cell type 1
# beta[301:600, 2, 1] <- 0.1 #cell type 2
# beta[201:500, 3, 1] <- 0.1 #cell type 3

# #age
# beta[601:900, 1, 2] <- 0.008  #cell type 1
# beta[601:800, 2, 2] <- 0.008  #cell type 2
# beta[701:900, 3, 2] <- 0.008 #cell type 3

#generate the cellular compositions
P <- sapply(1:n, function(i){
  if(X[1,i]==0){ #if control
    rDirichlet(c(4,2, 5))
  }else{
    rDirichlet(c(2,4, 1))
  }
})

#generate the observed methylation profiles
Ometh <- NULL
for(i in 1:n){
  utmp <- t(sapply(1:m, function(j){
    tmp1 <- colSums(X[ ,i] * t(beta[j, , ]))
    rnorm(K,mean=mu[j, ]+tmp1,sd=0.1)
  }))
  tmp2 <- colSums(P[ ,i] * t(utmp))
  Ometh <- cbind(Ometh, tmp2 + rnorm(m, sd = 0.1))
}

Ometh[Ometh > 1] <- 1
Ometh[Ometh < 0] <- 0



################################################################################################
#The Application of HIRE
################################################################################################

#return list by HIRE
ret_list <- HIRE_sq(Ometh, X, num_celltype=K,tol=1e-3,num_iter = 1000)

ret_list_full <- HIRE_sq_full(Ometh, X, num_celltype=K,tol=1e-3,num_iter = 500,P_t = ret_list$P_t,mu_t = ret_list$mu_t,sig_sqTiss_t = ret_list$sig_sqTiss_t,sig_sqErr_t = ret_list$sig_sqErr_t)

#case vs control
#the risk pattern for the first 1000 CpG sites
riskCpGpattern(ret_list$pvalues[1:1000, c(2,1,3)],
               main_title="Detected association pattern with disease status", hc_row_ind = FALSE)
#c(2,1,3) was used because of the label switching

#age
#the risk pattern for the first 1000 CpG sites
riskCpGpattern(ret_list$pvalues[1:1000, K+c(2,1,3)],
               main_title="Detected association pattern with age", hc_row_ind = FALSE)
