#ifndef __LASVDGPMS_H__
#define __LASVDGPMS_H__
void jmlelasvdGPms(lasvdGP *lasvdgp, unsigned int numstarts,
		   unsigned int maxit, unsigned int verb);

int iterlasvdGPms(lasvdGP* lasvdgp, unsigned int resvdThres,
		  unsigned int every, unsigned int numstarts,
		  unsigned int maxit, unsigned int verb);

#endif 
