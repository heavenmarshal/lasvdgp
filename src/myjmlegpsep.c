#include "matrix.h"
#include "gp_sep.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "lbfgs.h"
#include "rhelp.h"
#include <stdio.h>
#define SDEPS 1.490116e-08

struct mycallinfo_sep {
  GPsep *gpsep;
  double *dab;
  double *gab;
  double dmin;
  double dmax;
  int its;  /* updated but not used since lbfgsb counts fmin and gr evals */
  int verb;
};
static double trans(double x, double ymin, double ymax)
{
  double y;
  y = atan(x)/PI+0.5;
  y = ymin+(ymax-ymin)*y;
  return y;
}
static double invtrans(double y, double ymin, double ymax)
{
  double x;
  x = (y-ymin)/(ymax-ymin)-0.5;
  x = tan(PI*x);
  return x;
}
static double dtrans(double x, double ymin, double ymax)
{
  double dy;
  dy = (ymax-ymin)/PI;
  dy /= (1.0+x*x);
  return dy;
}
static lbfgsfloatval_t fnllik_sep(void* oinfo, const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *df, const int n,
				  const lbfgsfloatval_t step)
{
  double llik, *p;
  int psame, k, m;
  struct mycallinfo_sep* info = (struct mycallinfo_sep*) oinfo;
  m = info->gpsep->m;
  assert(n == m);
  p = (double*) malloc(n*sizeof(double));
  for(k=0; k<n; ++k) p[k] = trans(x[k],info->dmin,info->dmax);

  psame = 1;
  for(k=0; k<n; k++)
    if(p[k] != info->gpsep->d[k]) { psame = 0; break; }
  if(!psame) {
    (info->its)++;
    newparamsGPsep(info->gpsep, p, info->gpsep->g);
  }
  llik = llikGPsep(info->gpsep, info->dab, info->gab);

  dllikGPsep(info->gpsep, info->dab, df);
  for(k=0; k<n; k++) df[k] = 0.0-df[k]*dtrans(x[k],info->dmin,info->dmax);
  if(info->verb > 0) {
    MYprintf(MYstdout, "fmin it=%d, d=(%g", info->its, info->gpsep->d[0]);
    for(k=1; k<m; k++) MYprintf(MYstdout, " %g", info->gpsep->d[k]);
    if(n == m) MYprintf(MYstdout, "), llik=%g\n", llik);
    MYprintf(MYstdout, "dd=(%g", df[0]);
    for(k=1; k<n; k++) MYprintf(MYstdout, " %g", df[k]);
    MYprintf(MYstdout, ")\n");
  }
  free(p);
  return 0.0-llik;
}
static int fnllik_sep_trace(void* info, const lbfgsfloatval_t *x,
			    const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
			    const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm,
			    const lbfgsfloatval_t step, int n, int k, int ls)
{
    int i;
    MYprintf(MYstdout, "in lbfgs iteration %d, x=(%g",k, x[0]);
    for(i=1; i<n; ++i) MYprintf(MYstdout," %g",x[i]);
    MYprintf(MYstdout,")\n");
    MYprintf(MYstdout, "gradient=(%g",g[0]);
    for(i=1; i<n; ++i) MYprintf(MYstdout," %g",g[i]);
    MYprintf(MYstdout,")\n");
    MYprintf(MYstdout,"nllik=%g\n",fx);
    MYprintf(MYstdout,"norm of x = %g\n",xnorm);
    MYprintf(MYstdout,"norm of gradient = %g\n",gnorm);
    MYprintf(MYstdout,"step = %g\n",step);
    return 0;
}
void mymleGPsep(GPsep* gpsep, double dmin, double dmax, double *ab,
		const unsigned int maxit, int verb, double *p, int *its,
		char *msg, int *conv)
{
  double rmse;
  int k;
  double *dold;
  /* create structure for Brent_fmin */
  struct mycallinfo_sep info;
  info.gpsep = gpsep;
  info.dab = ab;
  info.gab = NULL;
  info.dmin = dmin;
  info.dmax = dmax;
  info.its = 0;
  info.verb = verb-6;

  lbfgs_parameter_t bfgsparams;

  lbfgs_parameter_init(&bfgsparams);
  bfgsparams.max_iterations = maxit;
  bfgsparams.delta = 1.0e-7;
  bfgsparams.linesearch=LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
  /* copy the starting value */
  for(k=0; k<gpsep->m; ++k) p[k]=invtrans(gpsep->d[k],dmin,dmax);

  dold = new_dup_vector(gpsep->d, gpsep->m);

  if(verb > 0) {
    MYprintf(MYstdout, "(d=[%g", gpsep->d[0]);
    for(k=1; k<gpsep->m; k++) MYprintf(MYstdout, ",%g", gpsep->d[k]);
    MYprintf(MYstdout, "], llik=%g) ", llikGPsep(gpsep, ab, NULL));
  }

  /* set ifail argument and verb/trace arguments */
  lbfgs_progress_t itertrace = (verb>0)? fnllik_sep_trace: NULL;

  *conv = lbfgs(gpsep->m, p, NULL, fnllik_sep, itertrace, (void *) &info, &bfgsparams);
  *its = info.its;
  /* check if parameters in p are new */
  rmse = 0.0;
  for(k=0; k<gpsep->m; k++) p[k] = trans(p[k],dmin,dmax);
  for(k=0; k<gpsep->m; k++) rmse += sq(p[k] - gpsep->d[k]);
  if(sqrt(rmse/k) > SDEPS) MYprintf(MYstderr,"stored d not same as d-hat\n");
  rmse = 0.0;
  for(k=0; k<gpsep->m; k++) rmse += sq(p[k] - dold[k]);
  if(sqrt(rmse/k) < SDEPS) {
    sprintf(msg, "lbfgs initialized at minima");
    *conv = 0;
    *its = 0;
  }

  /* print progress */
  if(verb > 0) {
    MYprintf(MYstdout, "-> %d lbfgsb its -> (d=[%g", *its, gpsep->d[0]);
    for(k=1; k<gpsep->m; k++) MYprintf(MYstdout, ",%g", gpsep->d[k]);
    MYprintf(MYstdout, "], llik=%g)\n", llikGPsep(gpsep, ab, NULL));
  }

  /* clean up */
  free(dold);
}
void myjmleGPsep(GPsep *gpsep, int maxit, double dmin, double dmax,
		 double *grange, double *dab, double *gab, int verb,
		 int *dits, int *gits, int *dconv)
{
  unsigned int i;
  int dit, git;
  char msg[60];
  double *d;

  /* sanity checks */
  assert(gab && dab);
  assert(dmin && dmax && grange);

  /* auxillary space for d-parameter values(s) */
  d = new_vector(gpsep->m);

  /* loop over coordinate-wise iterations */
  *dits = *gits = 0;
  for(i=0; i<100; i++) {
    mymleGPsep(gpsep, dmin, dmax, dab, maxit, verb, d, &dit, msg, dconv);
    *dits += dit;
    mleGPsep_nug(gpsep, grange[0], grange[1], gab, verb, &git);
    *gits += git;
    if((git <= 2 && (dit <= gpsep->m+1 && *dconv == 0)) || *dconv > 1) break;
  }
  if(i == 100 && verb > 0) warning("max outer its (N=100) reached");

  /* clean up */
  free(d);
}
