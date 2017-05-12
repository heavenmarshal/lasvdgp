#include "matrix.h"
#include "covar.h"
#include <math.h>
#include <assert.h>
void vector_minmax(double* vin, int len, double *min, double *max)
{
  int i;
  double tmax = vin[0];
  double tmin = vin[0];
  for(i = 0; i < len; ++i)
  {
    if(vin[i] > tmax)
    {
      tmax = vin[i];
      continue;
    }
    if(vin[i] < tmin)
      tmin = vin[i];
  }
  *min = tmin;
  *max = tmax;
}

int *nearest_indices(const unsigned int m, const unsigned int nref, double **Xref,
		     const unsigned int n, double **X, int* segs, const int nsegs)

{
  int i, close, start;
  int *oD;
  double **D;

  /* calculate distances to reference location(s), and so-order X & Z */
  D = new_matrix(nref, n);
  distance(Xref, nref, X, n, m, D);
  if(nref > 1) min_of_columns(*D, D, nref, n);

  close = segs[0];
  /* partition based on "close"st */

  oD = iseq(0, n-1);
  if(n > close) quick_select_index(*D, oD, n, close);

  /* now partition based on start */
  for(i = 1; i < nsegs; ++i)
  {
    close = segs[i-1];
    start = segs[i];
    quick_select_index(*D, oD, close, start);
  }
  delete_matrix(D);
  return(oD);
}

int fracvlen(double *v, double frac, unsigned int len)
{
  int i;
  double sum, psum;
  sum = sumv(v,len);
  psum = 0;
  for(i = 0; i < len; ++i)
  {
    psum += v[i];
    if(psum / sum > frac)
      break;
  }
  return i+1;
}
void sub_p_matrix_rows_col(double* vec, int* p, double **mat,
			   unsigned int col, unsigned int lenp)
{
  int i;
  assert(vec); assert(vec); assert(mat);
  for(i=0; i<lenp; ++i)
    vec[i] = mat[p[i]][col];
}
double* new_const_vector(double scalar, unsigned int n)
{
  int i;
  double *vec;
  vec = new_vector(n);
  for(i=0; i<n; ++i)
    vec[i] = scalar;
  
  return vec;
}
void sum_vector_scalar(double *v, double scalar, unsigned int n)
{
  int i;
  for(i=0; i<n; ++i)
    v[i] += scalar;
}
/* find the idx in vi that equal to val, if fail return -1 */
int find_int(int *vi, int val, unsigned int n)
{
  int i;
  for(i=0; i<n; ++i)
    if(vi[i] == val)
      return i;
  return -1;
}
void prod_vector(double *v1, double *v2, unsigned int n)
{
  int i;
  for(i = 0; i < n; ++i)
    v1[i] *= v2[i];
}
void divid_vector(double *v1, double *v2, unsigned int n)
{
  int i;
  for(i = 0; i < n; ++i)
    v1[i] /= v2[i];
}
double var_vector(double *v, double dividor, unsigned int n)
{
  int i;
  double sx = 0.0, sx2 = 0.0, var;
  for(i=0; i < n; ++i)
  {
    sx += v[i];
    sx2 += v[i]*v[i];
  }
  var = sx2 - sx*sx/(double)n;
  var /= dividor;
  return var;
}

double* new_sq_vector(double *v, unsigned int n)
{
  int i;
  double *res = new_vector(n);
  for(i=0; i<n; ++i)
    res[i] = v[i] * v[i];
  return res;
}

int ceil_divide(int n1, int n2)
{
  double dn1 = (double) n1;
  double dn2 = (double) n2;
  int res = (int) ceil(dn1/dn2);
  return res;
}
/* do not use square root for distance */
void distance_sym_vec(double **X, int n, int m, double *dist)
{
  int i, j, k, idx;
  double tmp;
  for(i = 0, idx = 0; i<n; ++i)
    for(j = i+1; j<n; ++j, ++idx)
    {
      tmp = 0.0;
      for(k = 0; k < m; ++k)
	tmp += sq(X[i][k] - X[j][k]);
      dist[idx] = tmp;
    }
}

int remove_nonpos(double *v, int n)
{
  int head = 0, tail;
  double tmp;
  for(tail = n-1; v[tail]<=0; --tail);
  while(head <= tail)
  {
    if(v[head] <= 0.0)
    {
      tmp = v[tail];
      v[tail] = v[head];
      v[head] = tmp;
      for(;v[tail]<=0; --tail);
    }
    head++;
  }
  return head;
}
void get_col(double *v, double **M, int col, int nrow)
{
  int i;
  for(i=0; i<nrow; ++i)
    v[i] = M[i][col];
}
