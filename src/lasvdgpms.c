#include <assert.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <Rmath.h>
#include "lasvdgp.h"
#include "matrix.h"
#include "linalg.h"
#include "matrixext.h"
#include "linalgext.h"
#include "lhsdefines.h"

static const double dab1 = 1.5;
static const double numdab2 = 3.907364;
static const double quanp = 0.1;
static const double sqreps = 1.490116119384766E-8;
static double gab[2] = {0.0, 0.0};

static void genStarts(int numstarts, int dim, double logdmin, double logdmax,
		      double loggmin, double loggmax, double** dstarts, double* gstarts)
{
  int i, j, din, dup;
  int **results;
  double tval, eps, dnumstarts;
  double ldrange, lgrange;
  din = dim + 1;
  dup = 1;
  results = new_imatrix(numstarts,din);
  maximinLHS_C(&numstarts, &din, &dup, *results);
  dnumstarts = (double) numstarts;
  ldrange = logdmax - logdmin;
  lgrange = loggmax - loggmin;
  for(i = 0; i < numstarts; ++i)
  {
    for(j = 0; j < dim; ++j)
    {
      eps = runif(0.0, 1.0);		/* look at Rmath.h */
      tval = (double)(results[i][j]-1);
      tval += eps;
      tval /= dnumstarts;
      dstarts[i][j] = exp(logdmin+tval*ldrange);
    }
    eps = runif(0.0, 1.0);
    tval = (double)(results[i][dim]-1);
    tval += eps;
    tval /= dnumstarts;
    gstarts[i] = exp(loggmin+tval*lgrange);
  }
  delete_imatrix(results);
}
void jmlelasvdGPms(lasvdGP *lasvdgp, unsigned int numstarts,
		   unsigned int maxit, unsigned int verb)
{
  double *dmin, *dmax;
  double dab[2], grange[2]={sqreps,lasvdgp->gstart};
  double dstart, ddmin, ddmax, dab2;
  double **dstarts, *gstarts;
  double *optd, optg, optllik, llik;
  int i, j, dim, dits, gits, dconv;
  getDs(lasvdgp->gpseps[0]->X,lasvdgp->n0,lasvdgp->m, &dstart, &ddmin, &ddmax,
	&dab2);
  dmin = new_const_vector(ddmin,lasvdgp->m);
  dmax = new_const_vector(ddmax,lasvdgp->m);
  dab[0] = dab1;
  dab[1] = dab2;
  /* get the start points */
  dim = lasvdgp-> m;
  dstarts = new_matrix(numstarts,dim);
  gstarts = new_vector(numstarts);
  genStarts(numstarts, dim, log(ddmin), log(ddmax),
	    log(grange[0]), log(grange[1]),
	    dstarts, gstarts);
  optd = new_vector(dim);
  for(i=0; i<lasvdgp->nbas; ++i)
  {
    optllik = -DBL_MAX;
    for(j = 0; j < numstarts; ++j)
    {
      newparamsGPsep(lasvdgp->gpseps[i], dstarts[j], gstarts[j]);
      jmleGPsep(lasvdgp->gpseps[i], maxit, dmin, dmax,
		grange, dab, gab, verb, &dits,
		&gits, &dconv, 1); /* fromR not sure */
      llik = llikGPsep(lasvdgp->gpseps[i], dab, gab);
      if(llik > optllik)
      {
	optllik = llik;
	dupv(optd,lasvdgp->gpseps[i]->d,dim);
	optg = lasvdgp->gpseps[i]->g;
      }
    }
    newparamsGPsep(lasvdgp->gpseps[i], optd, optg);
  }
  lasvdgp->hasfitted = 1;
  free(dmin);
  free(dmax);
  free(gstarts);
  free(optd);
  delete_matrix(dstarts);
}

void iterlasvdGPms(lasvdGP* lasvdgp, unsigned int resvdThres,
		   unsigned int every, unsigned int numstarts,
		   unsigned int maxit, unsigned int verb)
{
  int i, n0, nn, niter, nadd, nrem;
  nn = lasvdgp -> nn;
  n0 = lasvdgp -> n0;
  nadd = lasvdgp -> nadd;
  niter = ceil_divide(nn-n0,nadd);
  for(i = 1; i <= niter; ++i)
  {
    n0 = lasvdgp->n0;
    nrem = nn - n0;
    nadd = lasvdgp-> nadd;
    nadd = nadd<nrem ? nadd : nrem;
    lasvdgp -> nadd = nadd;
    selectNewPoints(lasvdgp);
    if(lasvdgp -> nappsvd >= resvdThres)
    {
      renewlasvdGP(lasvdgp);
      jmlelasvdGPms(lasvdgp, numstarts, maxit,verb);
      continue;
    }
    if(i % every == 0)
      jmlelasvdGP(lasvdgp, maxit, verb);
  }
  /* finishing off */
  if(lasvdgp->nappsvd>0)
  {
    renewlasvdGP(lasvdgp);
    jmlelasvdGPms(lasvdgp, numstarts, maxit, verb);
    return;
  }
  if(lasvdgp->hasfitted == 0)
    jmlelasvdGP(lasvdgp, maxit, verb);
}
void lasvdGPms_worker(double** X0, double **design, double **resp,
		      unsigned int M, unsigned int N, unsigned int m,
		      unsigned int tlen, unsigned int nn, unsigned int n0,
		      unsigned int nfea, unsigned int nsvd, unsigned int nadd,
		      double frac, double gstart, unsigned int resvdThres,
		      unsigned int every, unsigned int numstarts,
		      unsigned int maxit, unsigned int verb,
		      double **pmean, double **ps2)
{
  int i;
  double *xpred;
  lasvdGP *lasvdgp;
  for(i = 0; i < M; ++i)
  {
    xpred = X0[i];
    lasvdgp = newlasvdGP(xpred, design, resp, N, m, tlen, nn, n0,
			 nfea, nsvd, nadd, frac, gstart);
    jmlelasvdGPms(lasvdgp, numstarts, maxit, verb);
    iterlasvdGPms(lasvdgp, resvdThres, every, numstarts, maxit, verb);
    predlasvdGP(lasvdgp, pmean[i], ps2[i]);
    deletelasvdGP(lasvdgp);
  }
}

void lasvdGPms_R(double *X0_, double *design_, double *resp_, int* M_,
		 int *N_, int *m_, int *tlen_, int *nn_, int *n0_,
		 int *nfea_, int* nsvd_, int *nadd_, double *frac_,
		 double *gstart_, int *resvdThres_, int *every_,
		 int *numstarts_, int *maxit_, int *verb_,
		 double *pmean_, double *ps2_)
{
  double **X0, **design, **resp;
  double **pmean, **ps2;
  X0 = new_matrix_bones(X0_,*M_, *m_);
  design = new_matrix_bones(design_,*N_,*m_);
  resp = new_matrix_bones(resp_,*N_, *tlen_);
  pmean = new_matrix_bones(pmean_,*M_,*tlen_);
  ps2 = new_matrix_bones(ps2_,*M_,*tlen_);
  lasvdGPms_worker(X0,design,resp,*M_, *N_, *m_, *tlen_, *nn_, *n0_,
		   *nfea_, *nsvd_, *nadd_, *frac_, *gstart_, *resvdThres_,
		   *every_, *numstarts_, *maxit_, *verb_, pmean, ps2);
  free(X0);
  free(design);
  free(resp);
  free(pmean);
  free(ps2);
}
