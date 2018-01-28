#include <stdlib.h>
#include <R.h>

const char jobz = 'S';

int linalg_dgesdd(double **X, unsigned int nrow, unsigned int ncol,
		   double *s, double *u, double **vt)
{
  int info = 0, lwork = -1;
  int nsv = nrow<ncol? nrow : ncol;
  int *iwork = (int *) malloc(8*(size_t)(nsv)*sizeof(int));
  double tmp, *work;

  F77_CALL(dgesdd)(&jobz,&nrow,&ncol,*X,&nrow,s,u,&nrow,
		   *vt,&nsv,&tmp,&lwork,iwork, &info);
  if(info != 0)
    error("error code %d from Lapack routine '%s'", info, "dgesdd");
  lwork = (int) tmp;

  work = (double*) malloc(lwork * sizeof(double));

  F77_CALL(dgesdd)(&jobz,&nrow,&ncol,*X,&nrow,s,u,&nrow,
		   *vt,&nsv,work,&lwork,iwork,&info);
  free(work);
  free(iwork);
  return info;
  /* if(info != 0) */
  /*   error("error code %d from Lapack routine '%s'", info, "dgesdd"); */
  
}
