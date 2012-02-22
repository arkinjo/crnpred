#ifndef __chkaccu_h_
#define __chkaccu_h_

#include "sov.h"

/* $Id: chkaccu.h,v 1.2 2006/02/20 08:48:06 akinjo Exp $ */

extern void check_accuracy_ss(ans_t *aans, char psec[], 
			      double *q3, double *sov3);

extern void check_accuracy_cn(pssm_t *apssm, ans_t *aans, double pred[], 
			      double *cor, double *deva);

extern void check_accuracy_rwco(ans_t *aans, double pred[], 
				double *cor, double *deva);

extern double cnave2[]; /* average contact number for each residue type */
#endif /* __chkaccu_h_ */
