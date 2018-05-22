#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include "HIRE.h"
#include "HIRE_sq.h"

SEXP getListElement (SEXP list, char *str) {

	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;

	for(i = 0; i < length(list); i++) {
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0){
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}

	return elmt;
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////   ordinary EM   /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
SEXP EmEwasRcallC(SEXP args) {
	int nProtected = 0;
	double **P_t,  **mu_t, ***beta_t;
	double sig_sqErr_t, sig_sqTiss_t;
	double **Ometh, **X;
	double BIC=-1;
	double *BICpointer;
	BICpointer = &BIC;
	int K, n, m, p;
	double tol;
	int num_iteration;

	SEXP list, P_init, mu_init, beta_init, beta_init_dim, \
			sig_sqTiss_init, sig_sqErr_init, Ometh_r, Ometh_r_dim, X_r, X_r_dim, tol_r, num_iter;

	//input variables from R
	list = args;
	P_init = getListElement(list, "P_init");
	//P_init_dim = getListElement(list, "P_init_dim");
	mu_init = getListElement(list, "mu_init");
	//mu_init_dim = getListElement(list, "mu_init_dim");
	beta_init = getListElement(list, "beta_init");
	beta_init_dim = getListElement(list, "beta_init_dim");
	Ometh_r = getListElement(list, "Ometh_r");
	Ometh_r_dim = getListElement(list, "Ometh_r_dim");
	X_r = getListElement(list, "X_r");
	X_r_dim = getListElement(list, "X_r_dim");
	sig_sqTiss_init = getListElement(list, "sig_sqTiss_init");
	sig_sqErr_init = getListElement(list, "sig_sqErr_init");
	tol_r = getListElement(list, "tol_r");
	num_iter = getListElement(list, "num_iter");

	m = INTEGER(Ometh_r_dim)[0];  // CpG site number
	n = INTEGER(Ometh_r_dim)[1];  // sample number
	p = INTEGER(X_r_dim)[0];      // covariate number
	K = INTEGER(beta_init_dim)[1];// cell number
	tol = REAL(tol_r)[0]; 		  // tolerance
	num_iteration = INTEGER(num_iter)[0]; // maximum iteration number

	//declare and initialize variables in C
	Ometh = make2Darray(m, n);
	X = make2Darray(p, n);
	P_t = make2Darray(K, n);
	mu_t = make2Darray(m, K);
	beta_t = make3Darray(m, K, p);
	sig_sqTiss_t = REAL(sig_sqTiss_init)[0];
	sig_sqErr_t = REAL(sig_sqErr_init)[0];
	double *psig_sqTiss_t;
	double *psig_sqErr_t;
	psig_sqTiss_t = &sig_sqTiss_t;
	psig_sqErr_t = &sig_sqErr_t;


	int i, j, k, ell;

	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			Ometh[j][i] = REAL(Ometh_r)[j+m*i];
		}

		for(ell=0; ell<p; ell++){
			X[ell][i] = REAL(X_r)[ell+p*i];
		}

		for(k=0; k<K; k++){
			P_t[k][i] = REAL(P_init)[k+K*i];
		}
	}

	for(j=0; j<m; j++){
		for(k=0; k<K; k++){
			mu_t[j][k] = REAL(mu_init)[j+m*k];
		}
	}


	double ***E_Sigma;
	E_Sigma = make3Darray(n, K, K);
	double ***E_mu;
	E_mu = make3Darray(n,m,K);

	EmEwas(P_t, mu_t, beta_t, psig_sqErr_t, psig_sqTiss_t, K, n,
					m, p, Ometh, X, BICpointer, tol, num_iteration,E_Sigma,E_mu);

	//return P_t
	SEXP ret_P_t, ret_P_t_dim;
	PROTECT(ret_P_t = allocVector(REALSXP, K*n));
	++nProtected;
	for(k=0; k < K; k++){
		for(i=0; i < n; i++){
			REAL(ret_P_t)[k + K*i] = P_t[k][i];
		}
	}
	PROTECT(ret_P_t_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_P_t_dim)[0] = K;
	INTEGER(ret_P_t_dim)[1] = n;
	setAttrib(ret_P_t, R_DimSymbol, ret_P_t_dim);

	//return mu_t
	SEXP ret_mu_t, ret_mu_t_dim;
	PROTECT(ret_mu_t = allocVector(REALSXP, m*K));
	++nProtected;
	for(j=0; j < m; j++){
		for(k=0; k < K; k++){
			REAL(ret_mu_t)[j + m*k] = mu_t[j][k];
		}
	}
	PROTECT(ret_mu_t_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_mu_t_dim)[0] = m;
	INTEGER(ret_mu_t_dim)[1] = K;
	setAttrib(ret_mu_t, R_DimSymbol, ret_mu_t_dim);

	//return beta_t
	SEXP ret_beta_t, ret_beta_t_dim;
	PROTECT(ret_beta_t = allocVector(REALSXP, m*K*p));
	++nProtected;
	for(j=0; j < m; j++){
		for(k=0; k < K; k++){
			for(ell=0; ell < p; ell++){
				REAL(ret_beta_t)[j + m*k + m*K*ell] = beta_t[j][k][ell];
			}
		}
	}
	PROTECT(ret_beta_t_dim = allocVector(INTSXP, 3));
	++nProtected;
	INTEGER(ret_beta_t_dim)[0] = m;
	INTEGER(ret_beta_t_dim)[1] = K;
	INTEGER(ret_beta_t_dim)[2] = p;
	setAttrib(ret_beta_t, R_DimSymbol, ret_beta_t_dim);

	//return E_Sigma
	SEXP ret_E_Sigma, ret_E_Sigma_dim;
	PROTECT(ret_E_Sigma = allocVector(REALSXP, n*K*K));
	++nProtected;
	for(i=0; i < n; i++){
	  for(int k1=0; k1 < K; k1++){
	    for(int k2=0; k2 < K; k2++){
	      REAL(ret_E_Sigma)[i + n*k1 + n*K*k2] = E_Sigma[i][k1][k2];
	    }
	  }
	}
	PROTECT(ret_E_Sigma_dim = allocVector(INTSXP, 3));
	++nProtected;
	INTEGER(ret_E_Sigma_dim)[0] = n;
	INTEGER(ret_E_Sigma_dim)[1] = K;
	INTEGER(ret_E_Sigma_dim)[2] = K;
	setAttrib(ret_E_Sigma, R_DimSymbol, ret_E_Sigma_dim);

	//return mu
	SEXP ret_E_mu, ret_E_mu_dim;
	PROTECT(ret_E_mu = allocVector(REALSXP, n*m*K));
	++nProtected;
	for(i=0; i < n; i++){
	  for(j=0; j < m; j++){
	    for(int k=0; k < K; k++){
	       REAL(ret_E_mu)[i + n*j + n*m*k] = E_mu[i][j][k];
	    }
	  }
	}
	PROTECT(ret_E_mu_dim = allocVector(INTSXP, 3));
	++nProtected;
	INTEGER(ret_E_mu_dim)[0] = n;
	INTEGER(ret_E_mu_dim)[1] = m;
	INTEGER(ret_E_mu_dim)[2] = K;
	setAttrib(ret_E_mu, R_DimSymbol, ret_E_mu_dim);

	//return sig_sqErr_t
	SEXP ret_sig_sqErr_t;
	PROTECT(ret_sig_sqErr_t = allocVector(REALSXP, 1));
	++nProtected;
	REAL(ret_sig_sqErr_t)[0] = sig_sqErr_t;

	//return sig_sqTiss_t
	SEXP ret_sig_sqTiss_t;
	PROTECT(ret_sig_sqTiss_t = allocVector(REALSXP, 1));
	++nProtected;
	REAL(ret_sig_sqTiss_t)[0] = sig_sqTiss_t;

	//return BIC
	SEXP ret_BIC;
	PROTECT(ret_BIC = allocVector(REALSXP, 1));
	++nProtected;
	REAL(ret_BIC)[0] = BIC;

	//return list
	SEXP retList;
	PROTECT(retList = allocVector(VECSXP, 8));
	++nProtected;

	SEXP names; // components names in the return list
	PROTECT(names = allocVector(STRSXP, 8));
	++nProtected;

	SET_STRING_ELT(names, 0, mkChar("P_t"));
	SET_STRING_ELT(names, 1, mkChar("mu_t"));
	SET_STRING_ELT(names, 2, mkChar("beta_t"));
	SET_STRING_ELT(names, 3, mkChar("sig_sqErr_t"));
	SET_STRING_ELT(names, 4, mkChar("sig_sqTiss_t"));
	SET_STRING_ELT(names, 5, mkChar("BIC"));
	SET_STRING_ELT(names, 6, mkChar("E_Sigma"));
	SET_STRING_ELT(names, 7, mkChar("E_mu"));

	//add elements to the return list
	SET_VECTOR_ELT(retList, 0, ret_P_t);
	SET_VECTOR_ELT(retList, 1, ret_mu_t);
	SET_VECTOR_ELT(retList, 2, ret_beta_t);
	SET_VECTOR_ELT(retList, 3, ret_sig_sqErr_t);
	SET_VECTOR_ELT(retList, 4, ret_sig_sqTiss_t);
	SET_VECTOR_ELT(retList, 5, ret_BIC);
	SET_VECTOR_ELT(retList, 6, ret_E_Sigma);
	SET_VECTOR_ELT(retList, 7, ret_E_mu);
	setAttrib(retList, R_NamesSymbol, names);

	UNPROTECT(nProtected);
	delet2Darray(P_t, K, n);
	delet2Darray(mu_t, m, K);
	delet2Darray(Ometh, m, n);
	delet2Darray(X, p, n);
	delet3Darray(beta_t, m, K, p);

	delet3Darray(E_Sigma, n, K, K);
	delet3Darray(E_mu, n, m, K);

	return retList;
}


////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////   SQUAREM   /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
SEXP EmEwasRcallC_sq(SEXP args) {
  int nProtected = 0;
  double **P_t,  **mu_t, ***beta_t;
  double sig_sqErr_t, sig_sqTiss_t;
  double **Ometh, **X;
  double BIC=-1;
  double *BICpointer;
  BICpointer = &BIC;
  int K, n, m, p;
  double tol;
  int num_iteration;

  SEXP list, P_init, mu_init, beta_init, beta_init_dim, \
  sig_sqTiss_init, sig_sqErr_init, Ometh_r, Ometh_r_dim, X_r, X_r_dim, tol_r, num_iter;

  //input variables from R
  list = args;
  P_init = getListElement(list, "P_init");
  //P_init_dim = getListElement(list, "P_init_dim");
  mu_init = getListElement(list, "mu_init");
  //mu_init_dim = getListElement(list, "mu_init_dim");
  beta_init = getListElement(list, "beta_init");
  beta_init_dim = getListElement(list, "beta_init_dim");
  Ometh_r = getListElement(list, "Ometh_r");
  Ometh_r_dim = getListElement(list, "Ometh_r_dim");
  X_r = getListElement(list, "X_r");
  X_r_dim = getListElement(list, "X_r_dim");
  sig_sqTiss_init = getListElement(list, "sig_sqTiss_init");
  sig_sqErr_init = getListElement(list, "sig_sqErr_init");
  tol_r = getListElement(list, "tol_r");
  num_iter = getListElement(list, "num_iter");

  m = INTEGER(Ometh_r_dim)[0];  // CpG site number
  n = INTEGER(Ometh_r_dim)[1];  // sample number
  p = INTEGER(X_r_dim)[0];      // covariate number
  K = INTEGER(beta_init_dim)[1];// cell number
  tol = REAL(tol_r)[0]; 		  // tolerance
  num_iteration = INTEGER(num_iter)[0]; // maximum iteration number

  //declare and initialize variables in C
  Ometh = make2Darray(m, n);
  X = make2Darray(p, n);
  P_t = make2Darray(K, n);
  mu_t = make2Darray(m, K);
  beta_t = make3Darray(m, K, p);
  sig_sqTiss_t = REAL(sig_sqTiss_init)[0];
  sig_sqErr_t = REAL(sig_sqErr_init)[0];
  double *psig_sqTiss_t;
  double *psig_sqErr_t;
  psig_sqTiss_t = &sig_sqTiss_t;
  psig_sqErr_t = &sig_sqErr_t;


  int i, j, k, ell;

  for(i=0; i<n; i++){
    for(j=0; j<m; j++){
      Ometh[j][i] = REAL(Ometh_r)[j+m*i];
    }

    for(ell=0; ell<p; ell++){
      X[ell][i] = REAL(X_r)[ell+p*i];
    }

    for(k=0; k<K; k++){
      P_t[k][i] = REAL(P_init)[k+K*i];
    }
  }

  for(j=0; j<m; j++){
    for(k=0; k<K; k++){
      mu_t[j][k] = REAL(mu_init)[j+m*k];
    }
  }

  int L = K*n + m*K + 2;
  double *para;
  para = (double *) malloc(L*sizeof(double));

  initialize(beta_t, psig_sqErr_t, psig_sqTiss_t, K, n,m, p);
  makeVec(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,m,n,K);

  fixptfn(para, Ometh, K, m, n, tol, num_iteration);

  deVec(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,m,n,K);

  //return P_t
  SEXP ret_P_t, ret_P_t_dim;
  PROTECT(ret_P_t = allocVector(REALSXP, K*n));
  ++nProtected;
  for(k=0; k < K; k++){
    for(i=0; i < n; i++){
      REAL(ret_P_t)[k + K*i] = P_t[k][i];
    }
  }
  PROTECT(ret_P_t_dim = allocVector(INTSXP, 2));
  ++nProtected;
  INTEGER(ret_P_t_dim)[0] = K;
  INTEGER(ret_P_t_dim)[1] = n;
  setAttrib(ret_P_t, R_DimSymbol, ret_P_t_dim);

  //return mu_t
  SEXP ret_mu_t, ret_mu_t_dim;
  PROTECT(ret_mu_t = allocVector(REALSXP, m*K));
  ++nProtected;
  for(j=0; j < m; j++){
    for(k=0; k < K; k++){
      REAL(ret_mu_t)[j + m*k] = mu_t[j][k];
    }
  }
  PROTECT(ret_mu_t_dim = allocVector(INTSXP, 2));
  ++nProtected;
  INTEGER(ret_mu_t_dim)[0] = m;
  INTEGER(ret_mu_t_dim)[1] = K;
  setAttrib(ret_mu_t, R_DimSymbol, ret_mu_t_dim);

  //return beta_t
  SEXP ret_beta_t, ret_beta_t_dim;
  PROTECT(ret_beta_t = allocVector(REALSXP, m*K*p));
  ++nProtected;
  for(j=0; j < m; j++){
    for(k=0; k < K; k++){
      for(ell=0; ell < p; ell++){
        REAL(ret_beta_t)[j + m*k + m*K*ell] = beta_t[j][k][ell];
      }
    }
  }
  PROTECT(ret_beta_t_dim = allocVector(INTSXP, 3));
  ++nProtected;
  INTEGER(ret_beta_t_dim)[0] = m;
  INTEGER(ret_beta_t_dim)[1] = K;
  INTEGER(ret_beta_t_dim)[2] = p;
  setAttrib(ret_beta_t, R_DimSymbol, ret_beta_t_dim);

  //return sig_sqErr_t
  SEXP ret_sig_sqErr_t;
  PROTECT(ret_sig_sqErr_t = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(ret_sig_sqErr_t)[0] = sig_sqErr_t;

  //return sig_sqTiss_t
  SEXP ret_sig_sqTiss_t;
  PROTECT(ret_sig_sqTiss_t = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(ret_sig_sqTiss_t)[0] = sig_sqTiss_t;

  //return BIC
  SEXP ret_BIC;
  PROTECT(ret_BIC = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(ret_BIC)[0] = BIC;

  //return para
  SEXP ret_para;
  PROTECT(ret_para = allocVector(REALSXP, L));
  ++nProtected;
  for(j=0; j < L; j++){
    REAL(ret_para)[j] = para[j];
  }

  //return list
  SEXP retList;
  PROTECT(retList = allocVector(VECSXP, 7));
  ++nProtected;

  SEXP names; // components names in the return list
  PROTECT(names = allocVector(STRSXP, 7));
  ++nProtected;

  SET_STRING_ELT(names, 0, mkChar("P_t"));
  SET_STRING_ELT(names, 1, mkChar("mu_t"));
  SET_STRING_ELT(names, 2, mkChar("beta_t"));
  SET_STRING_ELT(names, 3, mkChar("sig_sqErr_t"));
  SET_STRING_ELT(names, 4, mkChar("sig_sqTiss_t"));
  SET_STRING_ELT(names, 5, mkChar("BIC"));
  SET_STRING_ELT(names, 6, mkChar("para"));

  //add elements to the return list
  SET_VECTOR_ELT(retList, 0, ret_P_t);
  SET_VECTOR_ELT(retList, 1, ret_mu_t);
  SET_VECTOR_ELT(retList, 2, ret_beta_t);
  SET_VECTOR_ELT(retList, 3, ret_sig_sqErr_t);
  SET_VECTOR_ELT(retList, 4, ret_sig_sqTiss_t);
  SET_VECTOR_ELT(retList, 5, ret_BIC);
  SET_VECTOR_ELT(retList, 6, ret_para);
  setAttrib(retList, R_NamesSymbol, names);

  UNPROTECT(nProtected);
  delet2Darray(P_t, K, n);
  delet2Darray(mu_t, m, K);
  delet2Darray(Ometh, m, n);
  delet2Darray(X, p, n);
  delet3Darray(beta_t, m, K, p);

  free(para);

  return retList;
}



////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////   SQUAREM-full   ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
SEXP EmEwasRcallC_sq_full(SEXP args) {
  int nProtected = 0;
  double **P_t,  **mu_t, ***beta_t;
  double sig_sqErr_t, sig_sqTiss_t;
  double **Ometh, **X;
  double BIC=-1;
  double *BICpointer;
  BICpointer = &BIC;
  int K, n, m, p;
  double tol;
  int num_iteration;

  SEXP list, P_init, mu_init, beta_init, beta_init_dim, \
  sig_sqTiss_init, sig_sqErr_init, Ometh_r, Ometh_r_dim, X_r, X_r_dim, tol_r, num_iter;

  //input variables from R
  list = args;
  P_init = getListElement(list, "P_init");
  //P_init_dim = getListElement(list, "P_init_dim");
  mu_init = getListElement(list, "mu_init");
  //mu_init_dim = getListElement(list, "mu_init_dim");
  beta_init = getListElement(list, "beta_init");
  beta_init_dim = getListElement(list, "beta_init_dim");
  Ometh_r = getListElement(list, "Ometh_r");
  Ometh_r_dim = getListElement(list, "Ometh_r_dim");
  X_r = getListElement(list, "X_r");
  X_r_dim = getListElement(list, "X_r_dim");
  sig_sqTiss_init = getListElement(list, "sig_sqTiss_init");
  sig_sqErr_init = getListElement(list, "sig_sqErr_init");
  tol_r = getListElement(list, "tol_r");
  num_iter = getListElement(list, "num_iter");

  m = INTEGER(Ometh_r_dim)[0];  // CpG site number
  n = INTEGER(Ometh_r_dim)[1];  // sample number
  p = INTEGER(X_r_dim)[0];      // covariate number
  K = INTEGER(beta_init_dim)[1];// cell number
  tol = REAL(tol_r)[0]; 		  // tolerance
  num_iteration = INTEGER(num_iter)[0]; // maximum iteration number

  //declare and initialize variables in C
  Ometh = make2Darray(m, n);
  X = make2Darray(p, n);
  P_t = make2Darray(K, n);
  mu_t = make2Darray(m, K);
  beta_t = make3Darray(m, K, p);
  sig_sqTiss_t = REAL(sig_sqTiss_init)[0];
  sig_sqErr_t = REAL(sig_sqErr_init)[0];
  double *psig_sqTiss_t;
  double *psig_sqErr_t;
  psig_sqTiss_t = &sig_sqTiss_t;
  psig_sqErr_t = &sig_sqErr_t;


  int i, j, k, ell;

  for(i=0; i<n; i++){
    for(j=0; j<m; j++){
      Ometh[j][i] = REAL(Ometh_r)[j+m*i];
    }

    for(ell=0; ell<p; ell++){
      X[ell][i] = REAL(X_r)[ell+p*i];
    }

    for(k=0; k<K; k++){
      P_t[k][i] = REAL(P_init)[k+K*i];
    }
  }

  for(j=0; j<m; j++){
    for(k=0; k<K; k++){
      mu_t[j][k] = REAL(mu_init)[j+m*k];
    }
  }

  int L = K*n + m*K + K + 2;
  double *para;
  para = (double *) malloc(L*sizeof(double));

  //initialize(beta_t, psig_sqErr_t, psig_sqTiss_t, K, n,m, p);
  //makeVec(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,m,n,K);
  //fixptfn(para, Ometh, K, m, n, tol, num_iteration);
  //deVec(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t,m,n,K);

  double **loglike1;
  loglike1 = make2Darray(m,p);

  for(int j=0; j<m; j++){
    for(int ell=0; ell<p; ell++){
      makeVec_full(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t, beta_t,m,n,K,j,ell);
      fixptfn_full(para, Ometh, X, K, m, n, p, tol, num_iteration, loglike1, j, ell);
      deVec_full(para,P_t, mu_t, psig_sqErr_t, psig_sqTiss_t, beta_t,m,n,K,j,ell);
    }
  }




  //return P_t
  SEXP ret_P_t, ret_P_t_dim;
  PROTECT(ret_P_t = allocVector(REALSXP, K*n));
  ++nProtected;
  for(k=0; k < K; k++){
    for(i=0; i < n; i++){
      REAL(ret_P_t)[k + K*i] = P_t[k][i];
    }
  }
  PROTECT(ret_P_t_dim = allocVector(INTSXP, 2));
  ++nProtected;
  INTEGER(ret_P_t_dim)[0] = K;
  INTEGER(ret_P_t_dim)[1] = n;
  setAttrib(ret_P_t, R_DimSymbol, ret_P_t_dim);

  //return mu_t
  SEXP ret_mu_t, ret_mu_t_dim;
  PROTECT(ret_mu_t = allocVector(REALSXP, m*K));
  ++nProtected;
  for(j=0; j < m; j++){
    for(k=0; k < K; k++){
      REAL(ret_mu_t)[j + m*k] = mu_t[j][k];
    }
  }
  PROTECT(ret_mu_t_dim = allocVector(INTSXP, 2));
  ++nProtected;
  INTEGER(ret_mu_t_dim)[0] = m;
  INTEGER(ret_mu_t_dim)[1] = K;
  setAttrib(ret_mu_t, R_DimSymbol, ret_mu_t_dim);

  //return beta_t
  SEXP ret_beta_t, ret_beta_t_dim;
  PROTECT(ret_beta_t = allocVector(REALSXP, m*K*p));
  ++nProtected;
  for(j=0; j < m; j++){
    for(k=0; k < K; k++){
      for(ell=0; ell < p; ell++){
        REAL(ret_beta_t)[j + m*k + m*K*ell] = beta_t[j][k][ell];
      }
    }
  }
  PROTECT(ret_beta_t_dim = allocVector(INTSXP, 3));
  ++nProtected;
  INTEGER(ret_beta_t_dim)[0] = m;
  INTEGER(ret_beta_t_dim)[1] = K;
  INTEGER(ret_beta_t_dim)[2] = p;
  setAttrib(ret_beta_t, R_DimSymbol, ret_beta_t_dim);

  //return sig_sqErr_t
  SEXP ret_sig_sqErr_t;
  PROTECT(ret_sig_sqErr_t = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(ret_sig_sqErr_t)[0] = sig_sqErr_t;

  //return sig_sqTiss_t
  SEXP ret_sig_sqTiss_t;
  PROTECT(ret_sig_sqTiss_t = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(ret_sig_sqTiss_t)[0] = sig_sqTiss_t;

  //return BIC
  SEXP ret_BIC;
  PROTECT(ret_BIC = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(ret_BIC)[0] = BIC;

  //return loglike1
  SEXP ret_loglike1, ret_loglike1_dim;
  PROTECT(ret_loglike1 = allocVector(REALSXP, m*p));
  ++nProtected;
  for(j=0; j < m; j++){
    for(ell=0; ell < p; ell++){
      REAL(ret_loglike1)[j + m*ell] = loglike1[j][ell];
    }
  }
  PROTECT(ret_loglike1_dim = allocVector(INTSXP, 2));
  ++nProtected;
  INTEGER(ret_loglike1_dim)[0] = m;
  INTEGER(ret_loglike1_dim)[1] = p;
  setAttrib(ret_loglike1, R_DimSymbol, ret_loglike1_dim);

  //return para
  SEXP ret_para;
  PROTECT(ret_para = allocVector(REALSXP, L));
  ++nProtected;
  for(j=0; j < L; j++){
    REAL(ret_para)[j] = para[j];
  }

  //return list
  SEXP retList;
  PROTECT(retList = allocVector(VECSXP, 8));
  ++nProtected;

  SEXP names; // components names in the return list
  PROTECT(names = allocVector(STRSXP, 8));
  ++nProtected;

  SET_STRING_ELT(names, 0, mkChar("P_t"));
  SET_STRING_ELT(names, 1, mkChar("mu_t"));
  SET_STRING_ELT(names, 2, mkChar("beta_t"));
  SET_STRING_ELT(names, 3, mkChar("sig_sqErr_t"));
  SET_STRING_ELT(names, 4, mkChar("sig_sqTiss_t"));
  SET_STRING_ELT(names, 5, mkChar("BIC"));
  SET_STRING_ELT(names, 6, mkChar("loglike1"));
  SET_STRING_ELT(names, 7, mkChar("para"));

  //add elements to the return list
  SET_VECTOR_ELT(retList, 0, ret_P_t);
  SET_VECTOR_ELT(retList, 1, ret_mu_t);
  SET_VECTOR_ELT(retList, 2, ret_beta_t);
  SET_VECTOR_ELT(retList, 3, ret_sig_sqErr_t);
  SET_VECTOR_ELT(retList, 4, ret_sig_sqTiss_t);
  SET_VECTOR_ELT(retList, 5, ret_BIC);
  SET_VECTOR_ELT(retList, 6, ret_loglike1);
  SET_VECTOR_ELT(retList, 7, ret_para);
  setAttrib(retList, R_NamesSymbol, names);

  UNPROTECT(nProtected);
  delet2Darray(P_t, K, n);
  delet2Darray(mu_t, m, K);
  delet2Darray(Ometh, m, n);
  delet2Darray(X, p, n);
  delet3Darray(beta_t, m, K, p);

  free(para);
  delet2Darray(loglike1,m,p);
  return retList;
}
