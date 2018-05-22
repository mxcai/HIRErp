// functions for SQUAREM
#include <stdbool.h>
#include "HIRE.h"

double crossprod(double *a, double *b, double len){
  double c=0;
  for(int i=0; i < len; i++){
    c += a[i] * b[i];
  }
  return c;
}
/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////          Null model          ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
void makeVec(double *para, double **P_t, double **mu_t, double *sig_sqErr_t, double *sig_sqTiss_t,
             int m, int n, int K){
  for(int k=0; k<K; k++){
    for(int i=0; i<n; i++){
      para[(k*n)+i] = P_t[k][i];
    }
  }
  for(int j=0; j<m; j++){
    for(int k=0; k<K; k++){
      para[(K*n)+(j*K)+k] = mu_t[j][k];
    }
  }
  para[K*n + m*K] = *sig_sqErr_t;
  para[K*n + m*K + 1] = *sig_sqTiss_t;

  // Rprintf("para[0]: %lf\t;P_t[0][0]: %lf\t;para[1]: %lf\t;P_t[0][1]: %lf\t;para[2]: %lf\t;P_t[0][2]: %lf\n",para[0],P_t[0][0],para[n],P_t[1][0],para[K*n-1],P_t[K-1][n-1]);
  // Rprintf("%lf,%lf\n",sig_sqErr_t, sig_sqTiss_t);
}

void deVec(double *para, double **P_t, double **mu_t, double *sig_sqErr_t, double *sig_sqTiss_t,
             int m, int n, int K){
  for(int k=0; k<K; k++){
    for(int i=0; i<n; i++){
      P_t[k][i] = para[(k*n)+i];
    }
  }
  for(int j=0; j<m; j++){
    for(int k=0; k<K; k++){
      mu_t[j][k] = para[(K*n)+(j*K)+k];
    }
  }
  *sig_sqErr_t = para[K*n + m*K];
  *sig_sqTiss_t = para[K*n + m*K + 1];

  // Rprintf("para[0]: %lf\t;P_t[0][0]: %lf\t;para[1]: %lf\t;P_t[0][1]: %lf\t;para[2]: %lf\t;P_t[0][2]: %lf\n",para[0],P_t[0][0],para[n],P_t[1][0],para[K*n-1],P_t[K-1][n-1]);
  // Rprintf("%lf,%lf\n",sig_sqErr_t, sig_sqTiss_t);
}

double loglike(double *para, int K, int m, int n, double **Ometh){
  double **P_t,  **mu_t;
  P_t = make2Darray(K, n);
  mu_t = make2Darray(m, K);
  double sig_sqErr_t, sig_sqTiss_t;
  double *psig_sqTiss_t;
  double *psig_sqErr_t;
  psig_sqTiss_t = &sig_sqTiss_t;
  psig_sqErr_t = &sig_sqErr_t;

  deVec(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,m,n,K);

  double sum, s1, s2;
  sum = 0;
  for(int j=0; j<m; j++){
    for(int i=0; i<n; i++){
      s1=0;
      s2=0;
      for(int k=0; k<K;k++){
        s1 += P_t[k][i]*mu_t[j][k];
      }
      for(int k=0; k<K; k++){
        s2 += P_t[k][i]*P_t[k][i]*sig_sqTiss_t;
      }
      sum += normal_density_log(Ometh[j][i], s1, sqrt(s2 + sig_sqErr_t));
    }
  }
  delet2Darray(P_t, K, n);
  delet2Darray(mu_t, m, K);
  return sum;
}

void EMstep(double *para, int K, int m, int n, double **Ometh, double *para_new){ //tol is the tolerance deciding when the algorithm stops

  double **P_t,  **mu_t;
  P_t = make2Darray(K, n);
  mu_t = make2Darray(m, K);
  double sig_sqErr_t, sig_sqTiss_t;
  double *psig_sqTiss_t;
  double *psig_sqErr_t;
  psig_sqTiss_t = &sig_sqTiss_t;
  psig_sqErr_t = &sig_sqErr_t;

  double ***E_Sigma;
  E_Sigma = make3Darray(n, K, K);
  double ***E_mu;
  E_mu = make3Darray(n,m,K);

  deVec(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,m,n,K);

  // E-step
  double g, g_tmp, s1, s2;
  double *tmp1 = (double *)malloc(K*sizeof(double));
  double *tmp2 = (double *)malloc(K*sizeof(double));

  for(int i=0; i<n; i++){
    s1 = 0;
    for(int k=0; k<K; k++){
      s1 += P_t[k][i]*P_t[k][i];
    }
    g = s1 * sig_sqTiss_t / sig_sqErr_t;
    for(int k1=0; k1<K; k1++){
      for(int k2=k1; k2<K; k2++){
        s2 = 1.0/(1.0+g) * P_t[k1][i] * P_t[k2][i] * sig_sqTiss_t * sig_sqTiss_t;
        g_tmp = s2 / sig_sqErr_t;
        if(k2!=k1){
          E_Sigma[i][k1][k2] = (-g_tmp);
          E_Sigma[i][k2][k1] = E_Sigma[i][k1][k2];
        }else{
          E_Sigma[i][k1][k2] = sig_sqTiss_t - g_tmp;
        }
      }
    }
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){

      for(int k=0; k<K; k++){
        tmp1[k] = Ometh[j][i]*P_t[k][i]/sig_sqErr_t;
        tmp2[k] =  (mu_t[j][k])/sig_sqTiss_t;
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

  // M-step
  //==============================
  //update sig_sqTiss_t
  //==============================
  s1 = 0;
  for(int j=0; j<m; j++){
    for(int k=0;k<K; k++){
      for(int i=0; i<n; i++){
        s2 = pow(E_mu[i][j][k] - mu_t[j][k],2);
        s1 += (s2 + E_Sigma[i][k][k]);
      }
    }
  }
  sig_sqTiss_t = s1 / (n*m*K);

  //=================================
  // use coordinate descent
  //=================================

  //use coordinate descent algorithm to update mu_t and beta_t

  for(int j=0; j<m; j++){
    for(int k=0; k<K; k++){
      s1 = 0;
      for(int i=0; i<n; i++){
        s1 += E_mu[i][j][k];
      }

      mu_t[j][k] = s1 / n;
    }
  }

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
  sig_sqErr_t = sum/(n*m);

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
          s1 += (E_Sigma[i][k][k2] + E_mu[i][j][k]*E_mu[i][j][k2]) / sig_sqErr_t;
        }
        Dmat[k][k2] = s1;
        if(k2 > k){
          Dmat[k2][k] = Dmat[k][k2];
        }
      }

      s2 = 0;
      for(int j=0; j<m; j++){
        s2 += Ometh[j][i]*E_mu[i][j][k]/sig_sqErr_t;
      }
      dvec[k] = s2;
    }

    for(int k=0; k<K; k++){
      xvec0[k] = P_t[k][i];
    }

    tar2a = val2(P_t, psig_sqErr_t, K, m, Ometh, E_Sigma, E_mu, i);

    quadprog(Dmat, dvec, xvec, K);

    for(int k=0; k<K; k++){
      P_t[k][i] = xvec[k];
    }
    tar2b = val2(P_t, psig_sqErr_t, K, m, Ometh, E_Sigma, E_mu, i);

    if(tar2b > tar2a){//avoid offset effects when estimates are stable
      for(int k=0; k<K; k++){
        P_t[k][i] = xvec0[k];
      }
    }
  }
  free(xvec);
  free(xvec0);
  free(dvec);
  for(int k=0; k<K; k++){
    free(Dmat[k]);
    free(Dmat_inv[k]);
  }
  free(Dmat);
  free(Dmat_inv);

  makeVec(para_new,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,m,n,K);


  delet2Darray(P_t, K, n);
  delet2Darray(mu_t, m, K);
  delet3Darray(E_Sigma, n, K, K);
  delet3Darray(E_mu, n, m, K);
}

void fixptfn(double *para, double **Ometh, int K, int m, int n, double tol, int num_iteration){
  int L = K*n + m*K + 2;
  int feval = 0;
  double *p1, *p2, *q1, *q2, *pnew;
  double sr2, sq2, sv2, srv, alpha, lold, lnew;
  bool extrap;
  double stepmax =1;
  double mstep = 4;
  double dec = 1;
  p1 = (double *) malloc(L*sizeof(double));
  p2 = (double *) malloc(L*sizeof(double));
  q1 = (double *) malloc(L*sizeof(double));
  q2 = (double *) malloc(L*sizeof(double));
  pnew = (double *) malloc(L*sizeof(double));

  lold = loglike(para, K, m, n, Ometh);
  int level = 1;

  double check0 = lold;
  double check1 = lold;

  for(int iter=0; iter<num_iteration; iter++){
    extrap = true;
    EMstep(para, K, m, n, Ometh, p1);
    feval += 1;
    sr2 = 0;
    for(int j=0; j<L; j++){
      q1[j] = p1[j] - para[j];
      sr2 += pow(q1[j],2);
    }
    // if(sqrt(sr2) < tol){
    //   break;
    // }

    EMstep(p1, K, m, n, Ometh, p2);
    feval += 1;
    sq2 = 0;
    for(int j=0; j<L; j++){
      q2[j] = p2[j] - p1[j];
      sq2 += pow(q2[j],2);
    }
    // if(sqrt(sq2) < tol){
    //   break;
    // }

    sv2 = 0;
    srv = 0;
    for(int j=0; j<L; j++){
      sv2 += pow(q2[j]-q1[j],2);
      srv += q1[j]*(q2[j]-q1[j]);
    }
    alpha= sqrt(sr2/sv2);
    alpha = 1>(stepmax<alpha?stepmax:alpha)?1:(stepmax<alpha?stepmax:alpha);

    for(int j=0; j<L; j++){
      pnew[j] = para[j] + 2*alpha*q1[j] + pow(alpha,2)*(q2[j]-q1[j]);
    }
    if(absolute(alpha-1)>0.01){
      EMstep(pnew, K, m, n, Ometh, pnew);  //stabiliztion step
      feval += 1;
    }

    lnew = loglike(pnew, K, m, n, Ometh);
    level += 1;
    if(lnew < (lold-dec)){
      for(int j=0; j<L; j++){
        pnew[j] = p2[j];
      }
      lnew = loglike(pnew, K, m, n, Ometh);
      level += 1;
      if(alpha == stepmax){
        stepmax = 1>(stepmax/mstep)?1:(stepmax/mstep);
      }
      alpha = 1;
      extrap = false;
      // Rprintf("Extrapolation: %d\t Step length: %lf\n", extrap, alpha);
    }

    if(alpha == stepmax){
      stepmax = mstep * stepmax;
    }

    for(int j=0; j<L; j++){
      para[j] = pnew[j];
    }
    lold = lnew;

    Rprintf("iter: %d\t Log-likelihood: %lf\t Extrapolation: %d\t Step length: %lf\n", iter, lnew, extrap, alpha);

    if(iter%10 == 0){
      check1 = lnew;
      if(absolute(check1 - check0) < tol*check0){
        break;
      }
      check0 = lnew;
    }
  }
  free(p1);
  free(p2);
  free(q1);
  free(q2);
  free(pnew);
}




/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////          full model          ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
void makeVec_full(double *para, double **P_t, double **mu_t, double *sig_sqErr_t, double *sig_sqTiss_t, double ***beta_t,
             int m, int n, int K, int j, int ell){
  for(int k=0; k<K; k++){
    for(int i=0; i<n; i++){
      para[(k*n)+i] = P_t[k][i];
    }
  }
  for(int j=0; j<m; j++){
    for(int k=0; k<K; k++){
      para[(K*n)+(j*K)+k] = mu_t[j][k];
    }
  }

  for(int k=0; k<K; k++){
    para[(K*n)+(m*K) + k] = beta_t[j][k][ell];
  }

  para[K*n + m*K + K] = *sig_sqErr_t;
  para[K*n + m*K + K + 1] = *sig_sqTiss_t;
}
void makeVec_full1(double *para, double **P_t, double **mu_t, double *sig_sqErr_t, double *sig_sqTiss_t, double *beta_t,
                  int m, int n, int K){
  for(int k=0; k<K; k++){
    for(int i=0; i<n; i++){
      para[(k*n)+i] = P_t[k][i];
    }
  }
  for(int j=0; j<m; j++){
    for(int k=0; k<K; k++){
      para[(K*n)+(j*K)+k] = mu_t[j][k];
    }
  }

  for(int k=0; k<K; k++){
    para[(K*n)+(m*K) + k] = beta_t[k];
  }

  para[K*n + m*K + K] = *sig_sqErr_t;
  para[K*n + m*K + K + 1] = *sig_sqTiss_t;
}

void deVec_full(double *para, double **P_t, double **mu_t, double *sig_sqErr_t, double *sig_sqTiss_t, double ***beta_t,
           int m, int n, int K, int j, int ell){
  for(int k=0; k<K; k++){
    for(int i=0; i<n; i++){
      P_t[k][i] = para[(k*n)+i];
    }
  }
  for(int j=0; j<m; j++){
    for(int k=0; k<K; k++){
      mu_t[j][k] = para[(K*n)+(j*K)+k];
    }
  }


  for(int k=0; k<K; k++){
    beta_t[j][k][ell] = para[(K*n)+(m*K) + k];
  }

  *sig_sqErr_t = para[K*n + m*K + K];
  *sig_sqTiss_t = para[K*n + m*K + K + 1];

}
void deVec_full1(double *para, double **P_t, double **mu_t, double *sig_sqErr_t, double *sig_sqTiss_t, double *beta_t,
                int m, int n, int K){
  for(int k=0; k<K; k++){
    for(int i=0; i<n; i++){
      P_t[k][i] = para[(k*n)+i];
    }
  }
  for(int j=0; j<m; j++){
    for(int k=0; k<K; k++){
      mu_t[j][k] = para[(K*n)+(j*K)+k];
    }
  }


  for(int k=0; k<K; k++){
    beta_t[k] = para[(K*n)+(m*K) + k];
  }

  *sig_sqErr_t = para[K*n + m*K + K];
  *sig_sqTiss_t = para[K*n + m*K + K + 1];

}

double loglike_full(double *para, int K, int m, int n, int p, double **Ometh, double **X, int J, int ELL){
  double **P_t,  **mu_t, *beta_t;
  P_t = make2Darray(K, n);
  mu_t = make2Darray(m, K);
  beta_t = (double *)malloc(K*sizeof(double));
  double sig_sqErr_t, sig_sqTiss_t;
  double *psig_sqTiss_t;
  double *psig_sqErr_t;
  psig_sqTiss_t = &sig_sqTiss_t;
  psig_sqErr_t = &sig_sqErr_t;

  deVec_full1(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,beta_t,m,n,K);

  double sum, s1, s2, s3;
  sum = 0;
  for(int j=0; j<m; j++){
    for(int i=0; i<n; i++){
      s1 = 0;
      s2 = 0;
      s3 = 0;
      for(int k=0; k<K;k++){
        s1 += P_t[k][i]*mu_t[j][k];
      }

      if(j==J){
        for(int k=0; k<K; k++){
          s2 += P_t[k][i]*beta_t[k]*X[ELL][i];
        }
      }

      for(int k=0; k<K; k++){
        s3 += P_t[k][i]*P_t[k][i]*sig_sqTiss_t;
      }
      sum += normal_density_log(Ometh[j][i], s1+s2, sqrt(s3 + sig_sqErr_t));
    }
  }
  delet2Darray(P_t, K, n);
  delet2Darray(mu_t, m, K);
  free(beta_t);

  return sum;
}

void EMstep_full(double *para, int K, int m, int n, int p, double **Ometh, double **X, double *para_new, int J, int ELL){ //tol is the tolerance deciding when the algorithm stops

  double **P_t,  **mu_t, *beta_t;
  P_t = make2Darray(K, n);
  mu_t = make2Darray(m, K);
  beta_t = (double *)malloc(K*sizeof(double));
  double sig_sqErr_t, sig_sqTiss_t;
  double *psig_sqTiss_t;
  double *psig_sqErr_t;
  psig_sqTiss_t = &sig_sqTiss_t;
  psig_sqErr_t = &sig_sqErr_t;

  double ***E_Sigma;
  E_Sigma = make3Darray(n, K, K);
  double ***E_mu;
  E_mu = make3Darray(n,m,K);

  deVec_full1(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,beta_t,m,n,K);

  // E-step
  double g, g_tmp, s1, s2, y_tmp, sum_x_sq;
  double *tmp1 = (double *)malloc(K*sizeof(double));
  double *tmp2 = (double *)malloc(K*sizeof(double));

  for(int i=0; i<n; i++){
    s1 = 0;
    for(int k=0; k<K; k++){
      s1 += P_t[k][i]*P_t[k][i];
    }
    g = s1 * sig_sqTiss_t / sig_sqErr_t;
    for(int k1=0; k1<K; k1++){
      for(int k2=k1; k2<K; k2++){
        s2 = 1.0/(1.0+g) * P_t[k1][i] * P_t[k2][i] * sig_sqTiss_t * sig_sqTiss_t;
        g_tmp = s2 / sig_sqErr_t;
        if(k2!=k1){
          E_Sigma[i][k1][k2] = (-g_tmp);
          E_Sigma[i][k2][k1] = E_Sigma[i][k1][k2];
        }else{
          E_Sigma[i][k1][k2] = sig_sqTiss_t - g_tmp;
        }
      }
    }
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){

      for(int k=0; k<K; k++){
        tmp1[k] = Ometh[j][i]*P_t[k][i]/sig_sqErr_t;
        s1 = 0;
        if(j==J){
            s1 = beta_t[k]*X[ELL][i];
        }
        tmp2[k] =  (mu_t[j][k] + s1)/sig_sqTiss_t;
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

  // M-step
  //==============================
  //update sig_sqTiss_t
  //==============================
  s1 = 0;
  for(int j=0; j<m; j++){
    for(int k=0;k<K; k++){
      for(int i=0; i<n; i++){
        s2 = 0;
        if(j==J){
          s2 = beta_t[k]*X[ELL][i];
        }
        s2 = pow(E_mu[i][j][k] - mu_t[j][k] - s2,2);
        s1 += (s2 + E_Sigma[i][k][k]);
      }
    }
  }
  sig_sqTiss_t = s1 / (n*m*K);

  //=================================
  // use coordinate descent
  //=================================

  //use coordinate descent algorithm to update mu_t and beta_t

  for(int j=0; j<m; j++){
    for(int k=0; k<K; k++){
      s1 = 0;
      for(int i=0; i<n; i++){
        s2 = 0;
        if(j==J){
          s2 = beta_t[k]*X[ELL][i];
        }
        s1 += E_mu[i][j][k] - s2;
      }

      mu_t[j][k] = s1 / n;
    }
  }

  for(int k=0; k<K; k++){
    s1 = 0;
    sum_x_sq=0;
    for(int i=0; i<n; i++){
      y_tmp = E_mu[i][J][k] - mu_t[J][k];
      s1 += y_tmp * X[ELL][i];
      sum_x_sq += (X[ELL][i] * X[ELL][i]);
    }
    beta_t[k] = s1 / (sum_x_sq);
  }

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
  sig_sqErr_t = sum/(n*m);

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
          s1 += (E_Sigma[i][k][k2] + E_mu[i][j][k]*E_mu[i][j][k2]) / sig_sqErr_t;
        }
        Dmat[k][k2] = s1;
        if(k2 > k){
          Dmat[k2][k] = Dmat[k][k2];
        }
      }

      s2 = 0;
      for(int j=0; j<m; j++){
        s2 += Ometh[j][i]*E_mu[i][j][k]/sig_sqErr_t;
      }
      dvec[k] = s2;
    }

    for(int k=0; k<K; k++){
      xvec0[k] = P_t[k][i];
    }

    tar2a = val2(P_t, psig_sqErr_t, K, m, Ometh, E_Sigma, E_mu, i);

    quadprog(Dmat, dvec, xvec, K);

    for(int k=0; k<K; k++){
      P_t[k][i] = xvec[k];
    }
    tar2b = val2(P_t, psig_sqErr_t, K, m, Ometh, E_Sigma, E_mu, i);

    if(tar2b > tar2a){//avoid offset effects when estimates are stable
      for(int k=0; k<K; k++){
        P_t[k][i] = xvec0[k];
      }
    }
  }
  free(xvec);
  free(xvec0);
  free(dvec);
  for(int k=0; k<K; k++){
    free(Dmat[k]);
    free(Dmat_inv[k]);
  }
  free(Dmat);
  free(Dmat_inv);

  makeVec_full1(para_new,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,beta_t,m,n,K);

  delet2Darray(P_t, K, n);
  delet2Darray(mu_t, m, K);
  free(beta_t);
  delet3Darray(E_Sigma, n, K, K);
  delet3Darray(E_mu, n, m, K);
}

void fixptfn_full(double *para, double **Ometh, double **X, int K, int m, int n, int p, double tol,
                  int num_iteration, double **loglike1, int J, int ELL){
  int L = K*n + m*K + K + 2;
  int feval = 0;
  double *p1, *p2, *q1, *q2, *pnew;
  double sr2, sq2, sv2, srv, alpha, lold, lnew;
  bool extrap;
  double stepmax =1;
  double mstep = 4;
  double dec = 1;
  p1 = (double *) malloc(L*sizeof(double));
  p2 = (double *) malloc(L*sizeof(double));
  q1 = (double *) malloc(L*sizeof(double));
  q2 = (double *) malloc(L*sizeof(double));
  pnew = (double *) malloc(L*sizeof(double));

  lold = loglike_full(para, K, m, n, p, Ometh, X, J, ELL);
  int level = 1;

  double check0 = lold;
  double check1 = lold;

  for(int iter=0; iter<num_iteration; iter++){
    extrap = true;
    EMstep_full(para, K, m, n, p, Ometh, X, p1, J, ELL);
    // EMstep_full(para, K, m, n, p, Ometh, X, pnew, J, ELL);
    // lnew = loglike_full(pnew, K, m, n, p, Ometh, X, J, ELL);
    feval += 1;
    sr2 = 0;
    for(int j=0; j<L; j++){
      q1[j] = p1[j] - para[j];
      sr2 += pow(q1[j],2);
    }
    // if(sqrt(sr2) < tol){
    //   break;
    // }

    EMstep_full(p1, K, m, n, p, Ometh, X, p2, J, ELL);
    feval += 1;
    sq2 = 0;
    for(int j=0; j<L; j++){
      q2[j] = p2[j] - p1[j];
      sq2 += pow(q2[j],2);
    }
    // if(sqrt(sq2) < tol){
    //   break;
    // }

    sv2 = 0;
    srv = 0;
    for(int j=0; j<L; j++){
      sv2 += pow(q2[j]-q1[j],2);
      srv += q1[j]*(q2[j]-q1[j]);
    }
    alpha= sqrt(sr2/sv2);
    alpha = 1>(stepmax<alpha?stepmax:alpha)?1:(stepmax<alpha?stepmax:alpha);

    for(int j=0; j<L; j++){
      pnew[j] = para[j] + 2*alpha*q1[j] + pow(alpha,2)*(q2[j]-q1[j]);
    }
    if(absolute(alpha-1)>0.01){
      EMstep_full(pnew, K, m, n, p, Ometh, X, pnew, J, ELL);  //stabiliztion step
      feval += 1;
    }

    lnew = loglike_full(pnew, K, m, n, p, Ometh, X, J, ELL);
    level += 1;
    if(lnew < (lold-dec)){
      for(int j=0; j<L; j++){
        pnew[j] = p2[j];
      }
      lnew = loglike_full(pnew, K, m, n, p, Ometh, X, J, ELL);
      level += 1;
      if(alpha == stepmax){
        stepmax = 1>(stepmax/mstep)?1:(stepmax/mstep);
      }
      alpha = 1;
      extrap = false;
      // Rprintf("Extrapolation: %d\t Step length: %lf\n", extrap, alpha);
    }

    if(alpha == stepmax){
      stepmax = mstep * stepmax;
    }

    for(int j=0; j<L; j++){
      para[j] = pnew[j];
    }
    lold = lnew;

    Rprintf("Phenotype: %d\t CpGsite: %d\t Iteration: %d\t observed-data log likelihood: %lf\t Extrapolation: %d\t Step length:%lf\n", ELL+1, J+1, iter, lnew, extrap, alpha);

    if(iter%10 == 0){
      check1 = lnew;
      if((check1 - check0) < tol){
        break;
      }
      check0 = lnew;
    }
  }
  loglike1[J][ELL] = lnew;
  free(p1);
  free(p2);
  free(q1);
  free(q2);
  free(pnew);
}

