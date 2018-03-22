#include "matrix.h"
#include "linalg.h"
#include "matrixext.h"
#include "lasvdgp.h"
#include "lasvdgpms.h"
#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

void predlasvdgpCall(lasvdGP* lasvdgp, double* pmean, double* ps2mode, double* ps2mean,
		     double* rvmode, double* rvmean)
{
  int i, n0, tlen, nbas;
  double **resid, **coeff;
  double *cmean, *cs2, *cdf, *bassq, ress2mode, ress2mean;
  GPsep **gpseps;

  assert(pmean);
  assert(ps2mode);
  assert(ps2mean);
  
  gpseps = lasvdgp->gpseps;
  n0 = lasvdgp -> n0;
  tlen = lasvdgp -> tlen;
  nbas = lasvdgp -> nbas;
  coeff = new_zero_matrix(nbas,n0);
  for(i=0; i < nbas; ++i)
    linalg_daxpy(n0,lasvdgp->reds[i], gpseps[i]->Z,1,coeff[i],1);
  resid = new_p_submatrix_rows(lasvdgp->feaidx,lasvdgp->resp, n0, tlen, 0);
  linalg_dgemm(CblasNoTrans,CblasTrans,tlen,n0,nbas,-1.0,&(lasvdgp->basis),tlen,
	       coeff,n0,1.0,resid,tlen);
  /* Y-USV^T */
  ress2mode = var_vector(*resid,(double)(n0*tlen+2), n0*tlen);
  ress2mean = var_vector(*resid,(double)(n0*tlen-2), n0*tlen);
  *rvmode = ress2mode;
  *rvmean = ress2mean;
  cmean = new_vector(nbas);
  cs2 = new_vector(nbas);
  cdf = new_vector(nbas);
  for(i=0; i<nbas; ++i)
    predGPsep_lite(gpseps[i], 1, &(lasvdgp->xpred), cmean+i, cs2+i, cdf+i,NULL);
  prod_vector(cmean,lasvdgp->reds, nbas);
  prod_vector(cs2,lasvdgp->reds, nbas);
  prod_vector(cs2,lasvdgp->reds, nbas);
  linalg_dgemv(CblasNoTrans,tlen,nbas,1.0,&(lasvdgp->basis),tlen,cmean,1,0.0,pmean,1);
  bassq = new_sq_vector(lasvdgp->basis,tlen*nbas);
  linalg_dgemv(CblasNoTrans,tlen,nbas,1.0,&bassq,tlen,cs2,1,0.0,ps2mode,1);
  dupv(ps2mean,ps2mode,tlen);
  sum_vector_scalar(ps2mode,ress2mode,tlen);
  sum_vector_scalar(ps2mean,ress2mean,tlen);
  delete_matrix(coeff);
  delete_matrix(resid);
  free(cmean);
  free(cs2);
  free(cdf);
  free(bassq);
}

SEXP lasvdgp_Call(SEXP X0_, SEXP design_, SEXP resp_, SEXP nn_,
		  SEXP n0_, SEXP nfea_, SEXP nsvd_, SEXP nadd_,
		  SEXP frac_, SEXP gstart_, SEXP resvdThres_,
		  SEXP every_, SEXP numstarts_, SEXP maxit_, SEXP verb_)
{
  SEXP ans, rglst;
  /* rglist list of range parameters */
  int i, j, M, N, m, tlen, nn, n0, nfea, nsvd;
  int nadd, resvdThres, every, maxit, verb, nbas, isMode, numstarts;
  double frac, gstart;
  double **X0, **design, **resp;
  double **pmean, **ps2mode, **ps2mean;
  double **range, *xpred, *ress2mode, *ress2mean;
  lasvdGP *lasvdgp;

  assert(isMatrix(X0_));
  assert(isMatrix(design_));
  assert(isMatrix(resp_));
  M = ncols(X0_);
  N = ncols(design_);
  m = nrows(design_);
  tlen = nrows(resp_);
  nn = INTEGER(nn_)[0];
  n0 = INTEGER(n0_)[0];
  nfea = INTEGER(nfea_)[0];
  nsvd = INTEGER(nsvd_)[0];
  nadd = INTEGER(nadd_)[0];
  resvdThres = INTEGER(resvdThres_)[0];
  every = INTEGER(every_)[0];
  numstarts = INTEGER(numstarts_)[0];
  maxit = INTEGER(maxit_)[0];
  verb = INTEGER(verb_)[0];
  frac = REAL(frac_)[0];
  gstart = REAL(gstart_)[0];
  X0 = new_matrix_bones(REAL(X0_),M,m);
  design = new_matrix_bones(REAL(design_),N,m);
  resp = new_matrix_bones(REAL(resp_),N,tlen);
  PROTECT(ans = allocVector(VECSXP,6));
  SET_VECTOR_ELT(ans,0,allocMatrix(REALSXP, tlen, M));
  SET_VECTOR_ELT(ans,1,allocMatrix(REALSXP, tlen, M));
  SET_VECTOR_ELT(ans,2,allocMatrix(REALSXP, tlen, M));
  SET_VECTOR_ELT(ans,3,allocVector(REALSXP, M));
  SET_VECTOR_ELT(ans,4,allocVector(REALSXP, M));
  pmean = new_matrix_bones(REAL(VECTOR_ELT(ans,0)), M, tlen);
  ps2mode = new_matrix_bones(REAL(VECTOR_ELT(ans,1)), M, tlen);
  ps2mean = new_matrix_bones(REAL(VECTOR_ELT(ans,2)), M, tlen);
  ress2mode = REAL(VECTOR_ELT(ans,3));
  ress2mean = REAL(VECTOR_ELT(ans,4));
  PROTECT(rglst = allocVector(VECSXP,M));
  for(i = 0; i < M; ++i)
  {
    xpred = X0[i];
    lasvdgp = newlasvdGP(xpred, design, resp, N, m, tlen, nn, n0,
			 nfea, nsvd, nadd, frac, gstart);
    jmlelasvdGPms(lasvdgp, numstarts, maxit, verb);
    iterlasvdGPms(lasvdgp, resvdThres, every, numstarts, maxit, verb);
    predlasvdgpCall(lasvdgp, pmean[i], ps2mode[i], ps2mean[i], ress2mode+i, ress2mean+i);
    nbas = lasvdgp->nbas;
    SET_VECTOR_ELT(rglst, i, allocMatrix(REALSXP,m,nbas));
    range = new_matrix_bones(REAL(VECTOR_ELT(rglst, i)), nbas, m);
    for(j = 0; j < nbas; ++j)
      dupv(range[j], lasvdgp->gpseps[j]->d, m);
    deletelasvdGP(lasvdgp);
    free(range);
  }
  SET_VECTOR_ELT(ans,5,rglst);
  free(X0);
  free(design);
  free(resp);
  free(pmean);
  free(ps2mode);
  free(ps2mean);
  UNPROTECT(2);

  return ans;
}
