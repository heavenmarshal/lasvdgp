#include "matrix.h"
#include "gp_sep.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "lbfgsb.h"
#include "rhelp.h"
#include <stdio.h>
#define SDEPS 1.490116e-08

struct mycallinfo_sep {
  GPsep *gpsep;
  double *dab;
  double *gab;
  int its;   /* updated but not used since lbfgsb counts fmin and gr evals */
  int verb;
};

static double nllik_sep(int n, double *p, struct mycallinfo_sep *info)
{
  double llik;
  int psame, k, m;

  /* sanity check */
  m = info->gpsep->m;
  assert(n == m);

  /* check if parameters in p are new */
  psame = 1;
  for(k=0; k<n; k++) {
    if(p[k] != info->gpsep->d[k]) { psame = 0; break; }
  }

  /* update GP with new parameters */
  if(!psame) {
    (info->its)++;
    newparamsGPsep(info->gpsep, p, info->gpsep->g);
  }

  /* evaluate likelihood with potentially new paramterization */
  llik = llikGPsep(info->gpsep, info->dab, info->gab);

  /* progress meter */
  if(info->verb > 0) {
    MYprintf(MYstdout, "fmin it=%d, d=(%g", info->its, info->gpsep->d[0]);
    for(k=1; k<m; k++) MYprintf(MYstdout, " %g", info->gpsep->d[k]);
    if(n == m) MYprintf(MYstdout, "), llik=%g\n", llik);
    else MYprintf(MYstdout, "), g=%g, llik=%g\n", info->gpsep->g, llik);
  }

  /* done */
  return 0.0-llik;
}

static void ndllik_sep(int n, double *p, double *df, struct mycallinfo_sep *info)
{
  int dsame, k;

  /* sanity check */
  assert(n == info->gpsep->m);

  /* check if parameters in p are new */
  dsame = 1;
  for(k=0; k<n; k++) if(p[k] != info->gpsep->d[k]) { dsame = 0; break; }

  /* update GP with new parameters */
  if(!dsame) {
    (info->its)++;
    newparamsGPsep(info->gpsep, p, info->gpsep->g);
  }

  /* evaluate likelihood with potentially new paramterization */
  dllikGPsep(info->gpsep, info->dab, df);

  /* negate values */
  for(k=0; k<n; k++) df[k] = 0.0-df[k];

  /* progress meter */
  if(info->verb > 1) {
    MYprintf(MYstdout, "grad it=%d, d=(%g", info->its, info->gpsep->d[0]);
    for(k=1; k<n; k++) MYprintf(MYstdout, " %g", info->gpsep->d[k]);
    MYprintf(MYstdout, "), dd=(%g", df[0]);
    for(k=1; k<n; k++) MYprintf(MYstdout, " %g", df[k]);
    MYprintf(MYstdout, ")\n");
  }
}
void mymleGPsep(GPsep* gpsep, double* dmin, double *dmax, double *ab,
		const unsigned int maxit, int verb, double *p, int *its,
		char *msg, int *conv)
{
  double rmse;
  int k, lbfgs_verb;
  double *dold;

  /* create structure for Brent_fmin */
  struct mycallinfo_sep info;
  info.gpsep = gpsep;
  info.dab = ab;
  info.gab = NULL;
  info.its = 0;
  info.verb = verb-6;

  /* copy the starting value */
  dupv(p, gpsep->d, gpsep->m);
  dold = new_dup_vector(gpsep->d, gpsep->m);

  if(verb > 0) {
    MYprintf(MYstdout, "(d=[%g", gpsep->d[0]);
    for(k=1; k<gpsep->m; k++) MYprintf(MYstdout, ",%g", gpsep->d[k]);
    MYprintf(MYstdout, "], llik=%g) ", llikGPsep(gpsep, ab, NULL));
  }

  /* set ifail argument and verb/trace arguments */
  *conv = 0;
  if(verb <= 1) lbfgs_verb = 0;
  else lbfgs_verb = verb - 1;

  /* call the C-routine behind R's optim function with method = "L-BFGS-B" */
  lbfgsb_C(gpsep->m, p, dmin, dmax, (lbfgsb_fmin) nllik_sep,
	   (lbfgsb_fgrad) ndllik_sep, conv, (void*) &info,
	   SDEPS, its, maxit, msg, lbfgs_verb);

  /* check if parameters in p are new */
  rmse = 0.0;
  for(k=0; k<gpsep->m; k++) rmse += sq(p[k] - gpsep->d[k]);
  if(sqrt(rmse/k) > SDEPS) MYprintf(MYstderr,"stored d not same as d-hat\n");
  rmse = 0.0;
  for(k=0; k<gpsep->m; k++) rmse += sq(p[k] - dold[k]);
  if(sqrt(rmse/k) < SDEPS) {
    sprintf(msg, "lbfgs initialized at minima");
    *conv = 0;
    its[0] = its[1] = 0;
  }

  /* print progress */
  if(verb > 0) {
    MYprintf(MYstdout, "-> %d lbfgsb its -> (d=[%g", its[1], gpsep->d[0]);
    for(k=1; k<gpsep->m; k++) MYprintf(MYstdout, ",%g", gpsep->d[k]);
    MYprintf(MYstdout, "], llik=%g)\n", llikGPsep(gpsep, ab, NULL));
  }

  /* clean up */
  free(dold);
}

void myjmleGPsep(GPsep *gpsep, int maxit, double *dmin, double *dmax,
		 double *grange, double *dab, double *gab, int verb,
		 int *dits, int *gits, int *dconv)
  {
    unsigned int i;
    int dit[2], git;
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
      mymleGPsep(gpsep, dmin, dmax, dab, maxit, verb, d, dit, msg, dconv);
      if(dit[1] > dit[0]) dit[0] = dit[1];
      *dits += dit[0];
      mleGPsep_nug(gpsep, grange[0], grange[1], gab, verb, &git);
      *gits += git;
      if((git <= 2 && (dit[0] <= gpsep->m+1 && *dconv == 0)) || *dconv > 1) break;
    }
    if(i == 100 && verb > 0) MYprintf(stderr,"max outer its (N=100) reached\n");

    /* clean up */
    free(d);
  }
