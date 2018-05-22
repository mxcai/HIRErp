rm(list=ls())
dir <- "/Users/cmx/Desktop/HIRE/"
# dir <- "/home/share/mingxuan/HIRE/"
X <- read.table(paste(dir,"Real_Data/RA/combined_covariates_no_name.txt",sep=""))
X <- data.matrix(X)
Ometh <- read.table(paste(dir,"Real_Data/RA/GSE42861_methylation_matr_after_quality_control_batch_effects_corrected_10K.txt",sep=""))
Ometh <- data.matrix(Ometh)


K <- 6

library(HIREewas)
system.time(fit_RA <- HIRE(Ometh, X, num_celltype=K))

library(HIRErp)
system.time(fit_RA_NULL <- HIRE_sq(Ometh, X, num_celltype=K,tol=1e-5,num_iter = 20000))
# system.time(fit_RA_full <- HIRE_sq_full(Ometh, X, num_celltype=K,tol=1e-3,num_iter = 500,P_t = fit_RA_NULL$P_t,mu_t = fit_RA_NULL$mu_t,sig_sqTiss_t = fit_RA_NULL$sig_sqTiss_t,sig_sqErr_t = fit_RA_NULL$sig_sqErr_t))


pval_HIRE <- cal_p(fit_RA$P_t, Ometh, X, cal_F = F)
pval_rp_null <- cal_p(fit_RA_NULL$P_t, Ometh, X, cal_F = F)

library(qqman)
png(file="/Users/cmx/Desktop/HIRE/pval_t2.png",width=6000,height=3000,res=600)
par(mfcol=c(2,6))
for(k in 1:K){
  qq(pval_HIRE$pvalues_t[,k])
  # qq(fit_RA$pvalues[,k])
  qq(pval_rp_null$pvalues_t[,k])
}
dev.off()

# png(file="/Users/cmx/Desktop/HIRE/pval_F1.png",width=5000,height=3000,res=600)
# par(mfcol=c(2,6))
# for(ell in 1:p){
#   qq(pval_HIRE$pvalues_F[,ell])
#   qq(pval_Null$pvalues_F[,ell])
# }
# dev.off()

