
//========================================================================================================================
//Define multiple dimensional arrays (Memory is not continuous here)
//========================================================================================================================

double**** make4Darray(int a1, int a2, int a3, int a4){
	double ****tmp;
	tmp = (double ****)malloc(a1*sizeof(double ***));
	for(int i=0; i<a1; i++){
		tmp[i] = (double ***) malloc(a2*sizeof(double **));
		for(int j=0; j<a2; j++){
			tmp[i][j] = (double **) malloc(a3*sizeof(double*));
			for(int k=0; k<a3; k++){
				tmp[i][j][k] = (double *)malloc(a4*sizeof(double));
			}
		}
	}
	return tmp;
}

void delet4Darray(double ****tmp, int a1, int a2, int a3, int a4){
	for(int i=0; i<a1; i++){
		for(int j=0; j<a2; j++){
			for(int k=0; k<a3; k++){
				free(tmp[i][j][k]);
			}
			free(tmp[i][j]);
		}
		free(tmp[i]);
	}
	free(tmp);
}

double*** make3Darray(int a1, int a2, int a3){
	double ***tmp;
	tmp = (double ***)malloc(a1*sizeof(double **));
	for(int i=0; i<a1; i++){
		tmp[i] = (double **) malloc(a2*sizeof(double *));
		for(int j=0; j<a2; j++){
			tmp[i][j] = (double *) malloc(a3*sizeof(double));
		}
	}
	return tmp;
}

void delet3Darray(double ***tmp, int a1, int a2, int a3){
	for(int i=0; i<a1; i++){
		for(int j=0; j<a2; j++){
			free(tmp[i][j]);
		}
		free(tmp[i]);
	}
	free(tmp);
}

double** make2Darray(int a1, int a2){
	double **tmp;
	tmp = (double **)malloc(a1*sizeof(double *));
	for(int i=0; i<a1; i++){
		tmp[i] = (double *) malloc(a2*sizeof(double));
	}
	return tmp;
}

void delet2Darray(double **tmp, int a1, int a2){
	for(int i=0; i<a1; i++){
		free(tmp[i]);
	}
	free(tmp);
}


//========================================================================================================================
//calculate the inverse of a matrix (the code is from http://www.cs.rochester.edu/u/brown/Crypto/assts/projects/adj.html)
//========================================================================================================================

/*Recursive definition of determinate using expansion by minors*/
double Determinant(double **a,int n){
   int i,j,j1,j2;
   double det = 0;
   double **m=NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = (double **) malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = (double *) malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}

/*Find the cofactor matrix of a square matrix*/

void CoFactor(double **a, int n, double **b){
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = (double **) malloc((n-1)*sizeof(double *));
   for (i=0;i<n-1;i++)
     c[i] = (double *) malloc((n-1)*sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}

/*Transpose of a square matrix, do it in place*/
void Transpose(double **a,int n){
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = tmp;
      }
   }
}

/*calculate the inverse*/
void inverse(double **a, int n, double **a_inv){
	double det;
	double  **cofac_a;
	cofac_a = (double **)malloc(n*sizeof(double*));
	for(int i =0; i <n;i++){
		cofac_a[i] = (double *)malloc(n*sizeof(double));
	}
	CoFactor(a, n, cofac_a);
	Transpose(cofac_a, n); //turn the cofacotor matrix into the adjoint matrix
	det = Determinant(a, n);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			a_inv[i][j] = cofac_a[i][j] / det;
		}
	}

	for(int i=0; i<n; i++){
		free(cofac_a[i]);
	}
	free(cofac_a);
}

//========================================================================================================================
//quadratic programming
//========================================================================================================================
void quadprog(double **Dmat, double *dvec, double *xvec, int K){ //K is the dimension of xvec
	if(K==1){
		xvec[0] = 1;
	}else{
		double **Dmat_inv, *x_star, s1, s2;
		double Dmat_inv_sum=0, x_star_sum=0, *Dmat_inv_rowsum;
		int num_negatives=0;
		double *ind_negative;
		Dmat_inv = (double **)malloc(K*sizeof(double*));
		x_star = (double *)malloc(K*sizeof(double));
		Dmat_inv_rowsum = (double *)malloc(K*sizeof(double));
		ind_negative = (double *)malloc(K*sizeof(double));
		for(int k=0; k<K; k++){
			Dmat_inv[k] = (double *)malloc(K*sizeof(double));
		}
		inverse(Dmat, K, Dmat_inv);
		for(int k=0; k<K; k++){
			s1 = 0;
			s2 = 0;
			for(int k1=0; k1<K; k1++){
				s1 += Dmat_inv[k][k1]*dvec[k1];
				Dmat_inv_sum += Dmat_inv[k][k1];
				s2 += Dmat_inv[k][k1];
			}
			x_star[k] = s1;
			Dmat_inv_rowsum[k] = s2;
			x_star_sum += s1;
		}
		for(int k=0; k<K; k++){
			xvec[k] = x_star[k] + (1-x_star_sum)/Dmat_inv_sum*Dmat_inv_rowsum[k];
			if(xvec[k]<0){
				num_negatives++;
				ind_negative[k] = 1;
			}else{
				ind_negative[k] = 0;
			}
		}
		free(x_star);
		free(Dmat_inv_rowsum);
		for(int k=0; k<K; k++){
			free(Dmat_inv[k]);
		}
		free(Dmat_inv);

		if(num_negatives == 0){
			free(ind_negative);
		}else{
			int Knew = K-num_negatives, i, j;
			double ** Dmat_new, *dvec_new, *xvec_sub;
			Dmat_new = (double **)malloc((Knew)*sizeof(double*));
			for(int k=0; k<Knew; k++){
				Dmat_new[k] = (double *) malloc((Knew)*sizeof(double));
			}
			dvec_new = (double *) malloc((Knew)*sizeof(double));
			xvec_sub = (double *) malloc((Knew)*sizeof(double));
			i = 0;

			for(int k1=0; k1<K; k1++){
				if(ind_negative[k1]==0){
					dvec_new[i] = dvec[k1];
					j = 0;
					for(int k2=0; k2<K; k2++){
						if(ind_negative[k2]==0){
							Dmat_new[i][j] = Dmat[k1][k2];
							j++;
						}
					}
					i++;
				}else{
					xvec[k1] = 0;
				}
			}

			quadprog(Dmat_new, dvec_new, xvec_sub, Knew);
			i=0;
			for(int k=0; k<K; k++){
				if(ind_negative[k]==0){
					xvec[k] = xvec_sub[i];
					i++;
				}
			}
			free(dvec_new);
			free(xvec_sub);
			for(int k=0; k<Knew;k++){
				free(Dmat_new[k]);
			}
			free(Dmat_new);
		}

	}

}


//the objective function when updating P_t
double val2(double **P_t, double *sig_sqErr_t,  int num_celltypes,
					int num_CpGsites,  double **Ometh, double ***E_Sigma, double ***E_mu, int i){
	int K,m;
	K = num_celltypes;
	m = num_CpGsites;

	double s1,s2,sum=0;

	for(int j=0; j<m; j++){
			s1=0;
			s2=0;
			for(int k=0; k<K; k++){
				for(int k2=0; k2<K; k2++){
					s1+= P_t[k][i]*E_Sigma[i][k][k2]*P_t[k2][i];
				}
				s2+= E_mu[i][j][k]*P_t[k][i];
			}
			s2 = Ometh[j][i] - s2;
			s2 = s2*s2;
			sum += (s1+s2)/(2*(*sig_sqErr_t));
	}

	return sum;
}

//========================================================================================================================
//generate random numbers
//========================================================================================================================

//generate one random number between 0 and 1 with 0 and 1 exclusive
double r_unif(){
	double temp;
  srand(1003211);
	do{
		temp = (double) rand() / (double) RAND_MAX;
	}while(temp>=1 || temp <=0);
	return temp;
}

//generate one noraml random number via Box-Muller method
double rnormal(double mu, double sd){
	double temp[2];
	temp[0] = r_unif();
	temp[1] = r_unif();
	double standardnormal = sqrt(-2 * log(temp[0])) * cos(2 * M_PI * temp[1]);
	return standardnormal * sd + mu;
}

//generate one gamma random number
double rgamma(double shape, double rate){
	double scale = 1.0 / rate;
	int shape_int = floor(shape);
	double s = 0;
	for(int i = 0; i < shape_int; i++){
		s = s - log(r_unif());
	}

	if(shape_int == shape){
		return scale * s;
	}else{
		double U, V, W, xi, eta;
		double delta = shape - shape_int;
		do{
			U = r_unif();
			V = r_unif();
			W = r_unif();
			if(U <= exp(1) / (exp(1) + delta)){
				xi = pow(V, 1 / delta);
				eta = W * pow(xi, delta - 1);
			}else{
				xi = 1 - log(V);
				eta = W * exp(-xi);
			}
		}while(eta > pow(xi, delta-1)*exp(-xi));
		return scale * ( xi + s);
	}

}

//generate one inverse gamma random number
double inv_rgamma(double shape, double rate){
	double temp = rgamma(shape, rate);
	return 1.0 / temp;
}

//generate a beta random number
double rbeta(double shape1, double shape2){
	double s1, s2;
	s1 = rgamma(shape1, 1);
	s2 = rgamma(shape2, 1);
	return s1 / (s1 + s2);
}

/* generate n Dirichlet random vectors with parameters being rep(2,K)
 K >= 2 and keep them in the P_t
 n is the column number of P_t
 K is the row number of P_t*/

/* void rdirichlet(double** P_t, int K, int n){
	double s;
	for(int i =0; i <n; i++){
		s = 0;
		for(int k = 0; k < K; k++){
			P_t[k][i] = rgamma(2.0, 1);
			s += P_t[k][i];
		}

		for(int k = 0; k < K; k++){
			P_t[k][i] /= s;
		}
	}
}*/

//========================================================================================================================
//functions in the algorithm
//========================================================================================================================

double normal_density_log(double x, double mean, double sd){
	return -0.5*log(2*M_PI) - log(sd) - pow(x-mean, 2) / (2*sd*sd);
}

// no use without lasso penalty
double positive_part(double x){
	if(x >= 0){
		return x;
	}else{
		return 0;
	}
}

double absolute(double x){
	if(x >= 0){
		return x;
	}else{
		return -x;
	}
}

double sign(double x){
	if(x > 0){
		return 1.0;
	}else if(x == 0){
		return 0.0;
	}else{
		return -1.0;
	}
}

//initialize sig_sqErr_t, sig_sqTiss_t, and beta_t
// P_t and mu_t have been initialized in the main function
// note that P_t was initialized using Houseman's quadratic programming method with a little modification

void initialize(double ***beta_t, double *sig_sqErr_t,
					double *sig_sqTiss_t, int num_celltypes, int num_samples,
					int num_CpGsites, int num_covariates){
	int K,m, p;
	K = num_celltypes;
	m = num_CpGsites;
	p = num_covariates;

	//initialize sig_sqTiss_t
	*sig_sqTiss_t = 0.01;

	for(int j=0; j<m; j++){
		for(int k=0; k<K; k++){
			//initialize beta_t
			for(int ell=0; ell<p; ell++){
				beta_t[j][k][ell] = 0;
			}
		}
	}
	//initialize sig_sqErr_t
	*sig_sqErr_t = 0.01;
}

//calculate the complete data likelihood
double complete_log_likelihood_CMX(double **P_t, double **mu_t, double ***beta_t, double *sig_sqErr_t,
                                   double **sig_sqTiss_t, int num_celltypes, int num_samples,
                                   int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled,
                                   double ***E_Sigma, double ***E_mu){
  int K,n, m, p;
  K = num_celltypes;
  n = num_samples;
  m = num_CpGsites;
  p = num_covariates;

  double llh = 0;
  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      double s1 = 0;
      for(int k1=0; k1<K; k1++){
        for(int k2=0; k2<K; k2++){
          s1 += P_t[k1][i]*E_Sigma[i][k1][k2]*P_t[k2][i];
        }
      }

      double s2 = 0;
      for(int k=0; k<K; k++){
        s2 += E_mu[i][j][k]*P_t[k][i];
      }

      double s3 = 0;
      for(int k=0; k<K; k++){
        double s4 = 0;
        for(int ell=0; ell<p; ell++){
          s4 += beta_t[j][k][ell]*X_scaled[ell][i];
        }
        s3 += -0.5*log(2*M_PI) - 0.5*log(sig_sqTiss_t[j][k]) - (E_Sigma[i][k][k]+pow(E_mu[i][j][k]-mu_t[j][k]-s4,2))/(2*sig_sqTiss_t[j][k]);
      }


      llh += -0.5*log(2*M_PI) - 0.5*log(sig_sqErr_t[j]) - (s1+pow(Ometh[j][i]-s2,2))/(2*sig_sqErr_t[j]) + s3;
    }
  }
  return llh;
}

//conduct cholesky decomposition
void cholesky(int N, double **M, double *diag)
{
  double *A = (double *) malloc(N*N*sizeof(double));

  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      A[j*N+i] = 1.0*M[i][j];
    }
  }

  for(int j=0;j<N;j++)
    diag[j] = A[N*j+j];
  for(int j=0;j<N;j++)
  {
    for(int k=0;k<j;k++)
      diag[j] -= A[N*k+j]*A[N*k+j];
    diag[j] = sqrt(diag[j]);
    for(int i=j+1;i<N;i++)
    {
      for(int k=0;k<j;k++)
        A[N*j+i] -= A[N*k+i]*A[N*k+j];
      A[N*j+i]/=diag[j];
    }
  }
  free(A);
}


//calculate the Q function
double Q_function(double **P_t, double **mu_t, double ***beta_t, double *sig_sqErr_t,
            double **sig_sqTiss_t, int num_celltypes, int num_samples,
            int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled,
            double ***E_Sigma, double ***E_mu){
  int K,n, m, p;
  K = num_celltypes;
  n = num_samples;
  m = num_CpGsites;
  p = num_covariates;
  double sum, s1, s2,s3;

  sum = 0;
  for(int j=0; j<m; j++){
    for(int i=0; i<n; i++){
      s1=0;
      for(int k=0; k<K;k++){
        s1 += P_t[k][i]*E_mu[i][j][k];
      }

      s3 = 0;
      for(int k=0; k<K; k++){
        s2 = 0;
        for(int ell=0; ell<p; ell++){
          s2 += beta_t[j][k][ell]*X_scaled[ell][i];
        }
        s3 += normal_density_log(E_mu[i][j][k],mu_t[j][k]+s2,sqrt(sig_sqTiss_t[j][k])) - E_Sigma[i][k][k]/(2*sig_sqTiss_t[j][k]);
      }

      s2 = 0;
      for(int k1 =0; k1<K; k1++){
        for(int k2=0; k2<K; k2++){
          s2+=P_t[k1][i] * E_Sigma[i][k1][k2] * P_t[k2][i];
        }
      }

      sum += normal_density_log(Ometh[j][i], s1, sqrt(sig_sqErr_t[j])) - s2/(2*sig_sqErr_t[j]) + s3;
    }
  }
  return sum;
}

//calculate the lower bound
double ELBO(double **P_t, double **mu_t, double ***beta_t, double *sig_sqErr_t,
                               double *sig_sqTiss_t, int num_celltypes, int num_samples,
                               int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled,
                               double ***E_Sigma, double ***E_mu){
  int K,n, m, p;
  K = num_celltypes;
  n = num_samples;
  m = num_CpGsites;
  p = num_covariates;
  double sum, s1, s2,s3, s4 ;
  double *tmp1 = (double *)malloc(K*sizeof(double));
  sum = 0;
  for(int j=0; j<m; j++){
    for(int i=0; i<n; i++){
      s1=0;
      for(int k=0; k<K;k++){
        s1 += P_t[k][i]*E_mu[i][j][k];
      }

      s3 = 0;
      for(int k=0; k<K; k++){
        s2 = 0;
        for(int ell=0; ell<p; ell++){
          s2 += beta_t[j][k][ell]*X_scaled[ell][i];
        }
        s3 += normal_density_log(E_mu[i][j][k],mu_t[j][k]+s2,sqrt(*sig_sqTiss_t)) - E_Sigma[i][k][k]/(2*(*sig_sqTiss_t));
      }

      s2 = 0;
      for(int k1 =0; k1<K; k1++){
        for(int k2=0; k2<K; k2++){
          s2+=P_t[k1][i] * E_Sigma[i][k1][k2] * P_t[k2][i];
        }
      }

      cholesky(K,E_Sigma[i],tmp1);
      s4 = 0;
      for(int k=0; k<K; k++){
        s4 += 2*log(tmp1[k]);
      }

      sum += normal_density_log(Ometh[j][i], s1, sqrt(*sig_sqErr_t)) - s2/(2*(*sig_sqErr_t)) + s3 + 0.5*s4 + 0.5*K*(1+log(2*M_PI));
    }
  }
  return sum;
  free(tmp1);
}

//calculate the observed data likelihood by CMX
double observed_log_likelihood_CMX(double **P_t, double **mu_t, double ***beta_t, double *sig_sqErr_t,
                                   double **sig_sqTiss_t, int num_celltypes, int num_samples,
                                   int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled){
  int K,n, m, p;
  K = num_celltypes;
  n = num_samples;
  m = num_CpGsites;
  p = num_covariates;

  double llh=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      double s1 = 0;
      double s2 = 0;
      for(int k=0; k<K ;k++){
        s1 += P_t[k][i]*sig_sqTiss_t[j][k]*P_t[k][i];
        double s3 = 0;
        for(int ell=0; ell<p; ell++){
          s3 += mu_t[j][k] + beta_t[j][k][ell]*X_scaled[ell][i];
        }
        s2 += P_t[k][i]*s3;
      }


      llh += -0.5*log(2*M_PI) - 0.5*log(s1+sig_sqErr_t[j]) - 0.5*pow(Ometh[j][i]-s2,2)/(s1+sig_sqErr_t[j]);
    }
  }
  return llh;
}

//calculate the observed data likelihood
double observed_log_likelihood(double **P_t, double **mu_t, double ***beta_t, double *sig_sqErr_t,
                               double *sig_sqTiss_t, int num_celltypes, int num_samples,
                               int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled){
  int K,n, m, p;
  K = num_celltypes;
  n = num_samples;
  m = num_CpGsites;
  p = num_covariates;
  double sum, s1, s2,s3 ;
  sum = 0;
  for(int j=0; j<m; j++){
    for(int i=0; i<n; i++){
      s1=0;
      s2=0;
      s3=0;
      for(int k=0; k<K;k++){
        s1 += P_t[k][i]*mu_t[j][k];
      }

      for(int k=0; k<K; k++){
        for(int ell=0; ell<p; ell++){
          s2 += P_t[k][i]*beta_t[j][k][ell]*X_scaled[ell][i];
        }
      }

      for(int k=0; k<K; k++){
        s3 += P_t[k][i]*P_t[k][i]*(*sig_sqTiss_t);
      }
      sum += normal_density_log(Ometh[j][i], s1+s2, sqrt(s3 + (*sig_sqErr_t)));
    }
  }
  return sum;
}

//conduct E step
void Estep(double **P_t, double **mu_t, double ***beta_t, double *sig_sqErr_t,
					double *sig_sqTiss_t, int num_celltypes, int num_samples,
					int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled,
					double ***E_Sigma, double ***E_mu ){
	int K,n, m, p;
	K = num_celltypes;
	n = num_samples;
	m = num_CpGsites;
	p = num_covariates;
	double g, g_tmp, s1, s2;
	double *tmp1 = (double *)malloc(K*sizeof(double));
	double *tmp2 = (double *)malloc(K*sizeof(double));

	// double lp0 = ELBO(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled, E_Sigma, E_mu);

	for(int i=0; i<n; i++){
	  s1 = 0;
	  for(int k=0; k<K; k++){
	    s1 += P_t[k][i]*P_t[k][i];
	  }
	  g = s1 * (*sig_sqTiss_t) / (*sig_sqErr_t);
	  for(int k1=0; k1<K; k1++){
	    for(int k2=k1; k2<K; k2++){
	      s2 = 1.0/(1.0+g) * P_t[k1][i] * P_t[k2][i] * (*sig_sqTiss_t) * (*sig_sqTiss_t);
	      g_tmp = s2 / (*sig_sqErr_t);
	      if(k2!=k1){
	        E_Sigma[i][k1][k2] = (-g_tmp);
	        E_Sigma[i][k2][k1] = E_Sigma[i][k1][k2];
	      }else{
	        E_Sigma[i][k1][k2] = (*sig_sqTiss_t) - g_tmp;
	      }
	    }
	  }
	}
	// double lp1 = ELBO(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled, E_Sigma, E_mu);
	// if(lp1<lp0){Rprintf("After E_Sigma: %lf\n", lp1-lp0);}

	for(int i=0; i<n; i++){
	  for(int j=0; j<m; j++){

			for(int k=0; k<K; k++){
				tmp1[k] = Ometh[j][i]*P_t[k][i]/(*sig_sqErr_t);
				s1 = 0;
				for(int ell=0; ell <p; ell++){
					s1+= beta_t[j][k][ell] * X_scaled[ell][i];
				}
				tmp2[k] =  (mu_t[j][k] + s1)/(*sig_sqTiss_t);
			}

			for(int k1=0; k1<K; k1++){
				s1 = 0;
				for(int k2=0; k2 <K; k2++){
					s1 += (tmp1[k2] + tmp2[k2])*E_Sigma[i][k2][k1];
				}
				E_mu[i][j][k1] = s1;
			}
		}
	}
	free(tmp1);
	free(tmp2);

	// double lp2 = ELBO(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled, E_Sigma, E_mu);
	// if(lp2<lp1){Rprintf("After E_mu: %lf\n", lp2-lp1);}

}



// conduct M step
void Mstep(double **P_t, double **mu_t, double ***beta_t, double *sig_sqErr_t,
					double *sig_sqTiss_t, int num_celltypes, int num_samples,
					int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled,
					double ***E_Sigma, double ***E_mu){
	int K,n, m, p;
	K = num_celltypes;
	n = num_samples;
	m = num_CpGsites;
	p = num_covariates;
	double s1, s2, y_tmp, sum_x_sq;
	//double **X_ext = make2Darray(p+1, n);
	//double **XX = make2Darray(p+1,p+1);
	//double **XX_inv = make2Darray(p+1,p+1);

	// double tp0 = observed_log_likelihood(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled);


	//==============================
	//update sig_sqTiss_t
	//==============================
	s1 = 0;
	for(int j=0; j<m; j++){
	  for(int k=0;k<K; k++){
	    for(int i=0; i<n; i++){
	      s2 = 0;
	      for(int ell=0; ell<p; ell++){
	        s2 += beta_t[j][k][ell]*X_scaled[ell][i];
	      }
	      s2 = pow(E_mu[i][j][k] - mu_t[j][k] - s2,2);
	      s1 += (s2 + E_Sigma[i][k][k]);
	    }
	  }
	}
	*sig_sqTiss_t = s1 / (n*m*K);


	// double tp3 = observed_log_likelihood(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled);
	// if(tp3<tp0){Rprintf("After sig_sqTiss: %lf\n", tp3-tp0);}

	//=================================
	// use coordinate descent
	//=================================

	//use coordinate descent algorithm to update mu_t and beta_t

	for(int j=0; j<m; j++){
		for(int k=0; k<K; k++){
			s1 = 0;
			for(int i=0; i<n; i++){
				s2 = 0;
				for(int ell=0; ell<p; ell++){
					s2 += beta_t[j][k][ell]*X_scaled[ell][i];
				}
				s1 += E_mu[i][j][k] - s2;
			}

			mu_t[j][k] = s1 / n;
			// if(mu_t[j][k]>1){
			// 	mu_t[j][k] = 1;
			// }else if(mu_t[j][k]<0){
			// 	mu_t[j][k] = 0;
			// }

			// for(int ell=0; ell<p; ell++){
			// 	s1 = 0;
			// 	sum_x_sq=0;
			// 	for(int i=0; i<n; i++){
			// 		s2 = 0;
			// 		for(int ell2=0; ell2<p; ell2++){
			// 			if(ell2 == ell){
			// 				continue;
			// 			}
			// 			s2+=beta_t[j][k][ell2]*X_scaled[ell2][i];
			// 		}
			// 		y_tmp = E_mu[i][j][k] - mu_t[j][k] -s2;
			// 		s1 += y_tmp * X_scaled[ell][i];
			// 		sum_x_sq += (X_scaled[ell][i] * X_scaled[ell][i]);
			// 	}
			//
			// 	beta_t[j][k][ell] = s1 / (sum_x_sq);
			// }
		}
	}

	// double tp1 = observed_log_likelihood(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled);
	// if(tp1<tp3){Rprintf("After mu, beta: %lf\n", tp1-tp3);}




	//==============================
	//update sig_sqErr_t
	//==============================
	double sum;
	sum = 0;
	for(int j=0; j<m; j++){
	  for(int i=0; i<n; i++){
	    s1 = 0;
	    for(int k1 =0; k1<K; k1++){
	      for(int k2=0; k2<K; k2++){
	        s1+=P_t[k1][i] * E_Sigma[i][k1][k2] * P_t[k2][i];
	      }
	    }
	    s2 = 0;
	    for(int k=0; k<K; k++){
	      s2 += E_mu[i][j][k]*P_t[k][i];
	    }
	    s2 = pow(Ometh[j][i] - s2,2);
	    sum += (s1+s2);
	  }
	}
	*sig_sqErr_t = sum/(n*m);

	// double tp2 = observed_log_likelihood(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled);
	// if(tp2<tp1){Rprintf("After sig_sqErr: %lf\n", tp2-tp1);}

	//==============================
	//update P_t
	//==============================
	double *dvec, **Dmat, **Dmat_inv, *xvec, *xvec0, tar2a, tar2b;
	xvec = (double *)malloc(K*sizeof(double));
	xvec0 = (double *)malloc(K*sizeof(double));
	dvec = (double *)malloc(K*sizeof(double));
	Dmat = (double **)malloc(K*sizeof(double*));
	Dmat_inv = (double **)malloc(K*sizeof(double*));
	for(int k=0; k<K; k++){
		Dmat[k] =  (double *)malloc(K*sizeof(double));
		Dmat_inv[k] =  (double *)malloc(K*sizeof(double));
	}
	for(int i=0; i<n; i++){

		for(int k=0; k <K; k++){
			for(int k2=k; k2<K; k2++){
				s1 = 0;
				for(int j=0; j<m; j++){
					s1 += (E_Sigma[i][k][k2] + E_mu[i][j][k]*E_mu[i][j][k2]) / (*sig_sqErr_t);
				}
				Dmat[k][k2] = s1;
				if(k2 > k){
					Dmat[k2][k] = Dmat[k][k2];
				}
			}

			s2 = 0;
			for(int j=0; j<m; j++){
				s2 += Ometh[j][i]*E_mu[i][j][k]/(*sig_sqErr_t);
			}
			dvec[k] = s2;
		}

		for(int k=0; k<K; k++){
			xvec0[k] = P_t[k][i];
		}

		tar2a = val2(P_t, sig_sqErr_t, K, m, Ometh, E_Sigma, E_mu, i);

		quadprog(Dmat, dvec, xvec, K);

		for(int k=0; k<K; k++){
			P_t[k][i] = xvec[k];
		}
		tar2b = val2(P_t, sig_sqErr_t, K, m, Ometh, E_Sigma, E_mu, i);

		if(tar2b > tar2a){//avoid offset effects when estimates are stable
			for(int k=0; k<K; k++){
				P_t[k][i] = xvec0[k];
			}
		}
	}

	// double tp4 = observed_log_likelihood(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled);
	// if(tp4<tp3){Rprintf("After P_t: %lf\n", tp4-tp3);}

	free(xvec);
	free(xvec0);
	free(dvec);
	for(int k=0; k<K; k++){
		free(Dmat[k]);
		free(Dmat_inv[k]);
	}
	free(Dmat);
	free(Dmat_inv);

}




//calculate the BIC to select an appropriate K
double BIC(double **P_t, double **mu_t, double ***beta_t, double *sig_sqErr_t,
					double **sig_sqTiss_t, int num_celltypes, int num_samples,
					int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled){
	int K,n, m, p;
	K = num_celltypes;
	n = num_samples;
	m = num_CpGsites;
	p = num_covariates;
	double sum, s1, s2,s3 ;
	sum = 0;
	for(int j=0; j<m; j++){
		for(int i=0; i<n; i++){
			s1=0;
			s2=0;
			s3=0;
			for(int k=0; k<K;k++){
				s1 += P_t[k][i]*mu_t[j][k];
			}

			for(int k=0; k<K; k++){
				for(int ell=0; ell<p; ell++){
					s2 += P_t[k][i]*beta_t[j][k][ell]*X_scaled[ell][i];
				}
			}

			for(int k=0; k<K; k++){
					s3 += P_t[k][i]*P_t[k][i]*sig_sqTiss_t[j][k];
			}
			sum += normal_density_log(Ometh[j][i], s1+s2, sqrt(s3 + sig_sqErr_t[j]));
		}
	}

	sum = (-2.0)*sum + log(n)*((K-1)*n + m*(1+2*K+p*K));
	return sum;
}

//carry out the EM algorithm
void EmEwas(double **P_t, double **mu_t, double ***beta_t, double *sig_sqErr_t,
					double *sig_sqTiss_t, int num_celltypes, int num_samples,
					int num_CpGsites, int num_covariates, double **Ometh, double **X_scaled,
					double *BICpointer, double tol, int num_iteration, double ***E_Sigma, double ***E_mu){ //tol is the tolerance deciding when the algorithm stops

	int K,n, m, p, t=1;
	K = num_celltypes;
	n = num_samples;
	m = num_CpGsites;
	p = num_covariates;
	// double ****E_Sigma;
	// E_Sigma = make4Darray(n, m, K, K);
	// double ***E_mu;
	// E_mu = make3Darray(n,m,K);

	double ploglike0=-1e10, ploglike1=1;
	initialize(beta_t, sig_sqErr_t, sig_sqTiss_t, K, n,m, p);

	// double tp0 = -1e10;
	while( (absolute(ploglike1-ploglike0) / absolute(ploglike0) > tol) && (t<=num_iteration)){ //relative tolerance
		ploglike0 = ploglike1;
		Estep(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled, E_Sigma, E_mu);


		// double diff_loglike_ELBO = ploglike0-ELBO(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled, E_Sigma, E_mu);

		Mstep(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled, E_Sigma, E_mu);

		ploglike1 = observed_log_likelihood(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled);
		// double tp1 = ELBO(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled, E_Sigma, E_mu);
		// Rprintf("Iteration: %d\t observed log-like: %lf\t ELBO: %lf\t log-like difference: %lf\t ELBO difference: %lf\t (log-like-ELBO at E-step): %e\n", t, ploglike1, tp1, ploglike1-ploglike0, tp1-tp0, diff_loglike_ELBO);
		Rprintf("Iteration: %d\t observed-data log likelihood: %lf\n", t, ploglike1);
		// tp0 = tp1;
		if(ploglike1 < ploglike0){
		printf("\t\t\t\t\t warning!!!\n");
		}
		t++;
	}

	// *BICpointer = BIC(P_t,mu_t,beta_t,sig_sqErr_t,sig_sqTiss_t, K,n,m,p, Ometh, X_scaled);
	// delet4Darray(E_Sigma, n, m, K, K);
	// delet3Darray(E_mu, n, m, K);
}


