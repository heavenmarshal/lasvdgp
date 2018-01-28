#include"CommonDefines.hpp"
#include"CRandom.hpp"
#include<time.h>
extern "C"{
    
#include"matrix.h"
    void genmaximinLHS(unsigned int n, unsigned int k, double **maximinlhs)
    {

	bclib::CRandomStandardUniform runif;
	bclib::matrix<int> intMat = bclib::matrix<int>(n,k);
	unsigned int col, row;
	unsigned int seed1 = (int) time(NULL);
	double eps, dn;
	runif.setSeed(seed1,seed1+1024);
	lhslib::maximinLHS(n,k,1,intMat,runif);
	dn = static_cast<double>(n);

	for(col = 0; col < k; ++col)
	    for(row = 0; row < n; ++row)
	    {
		eps = runif.getNextRandom();
		maximinlhs[row][col] = static_cast<double>(intMat(row, col) - 1) + eps;
		maximinlhs[row][col] /= dn;
	    }
    }
}
