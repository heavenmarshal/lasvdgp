#ifndef __MYJMLEGPSEP_H__
#define __MYJMLEGPSEP_H__
#include "gp_sep.h"
void myjmleGPsep(GPsep *gpsep, int maxit, double *dmin, double *dmax,
		 double *grange, double *dab, double *gab, int verb,
		 int *dits, int *gits, int *dconv);
#endif
