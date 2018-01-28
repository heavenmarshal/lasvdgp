#ifndef __MATRIXEXT_H__
#define __MATRIXEXT_H__
void vector_minmax(double*, int, double*, double*);
int *nearest_indices(const unsigned int, const unsigned int, double**,
		     const unsigned int, double**, int*, const int);
int fracvlen(double*, double, unsigned int);
void sub_p_matrix_rows_col(double*, int*, double**,
			   unsigned int, unsigned int);
double* new_const_vector(double, unsigned int);
void sum_vector_scalar(double*, double, unsigned int);
int find_int(int*, int, unsigned int);
void prod_vector(double *, double *, unsigned int);
void divid_vector(double *, double *, unsigned int);
double var_vector(double *, double , unsigned int);
double* new_sq_vector(double *, unsigned int);
int ceil_divide(int, int);
void distance_sym_vec(double**, int, int, double*);
int remove_nonpos(double*, int);
void get_col(double*, double **, int, int);
#endif
