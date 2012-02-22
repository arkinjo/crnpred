/* 	$Id: chkaccu.c,v 1.2 2006/02/20 08:30:52 akinjo Exp $	 */

#ifndef lint
static char vcid[] = "$Id: chkaccu.c,v 1.2 2006/02/20 08:30:52 akinjo Exp $";
#endif /* lint */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blast.h"


/* average contact number for each residue type */
double cnave2[] = {
  25.430 , /* A */
  21.038 , /* R */
  20.093 , /* N */
  18.594 , /* D */
  29.647 , /* C */
  20.206 , /* Q */
  18.008 , /* E */
  22.505 , /* G */
  23.572 , /* H */
  29.469 , /* I */
  28.173 , /* L */
  18.452 , /* K */
  26.466 , /* M */
  28.057 , /* F */
  20.350 , /* P */
  21.420 , /* S */
  22.747 , /* T */
  26.913 , /* W */
  26.627 , /* Y */
  28.656 , /* V */
};

void check_accuracy_cn(pssm_t *apssm, ans_t *aans, double pred[], 
		       double *cor, double *deva)
{
  int i,len;
  double ave=0.0,sd=0.0;
  double c=0.0,d=0.0;
  int ind, be0,be1,nn[2][2];
  double q2;

  len = aans->len;
  nn[0][0] = nn[1][1] = nn[1][0] = nn[0][1] = 0;
  for(i = 0; i < len; i++) {
    ave += pred[i];
    ind = strchr(amino1,apssm->seq[i]) - amino1;

    if(ind>=0) {
      be0 = (cnave2[ind] < aans->nvec[i]*log(len)) ? 0 : 1;
      be1 = (cnave2[ind] < pred[i]*log(len)) ? 0 : 1;
      nn[be0][be1]++;
    }
  }
  ave /= (double)len;
  q2 = 100.0*(nn[0][0] + nn[1][1])/(double)len;
  fprintf(stderr, "Q2: %15.7f %5d %5d\n", q2, nn[0][0]+nn[1][1], len);
  fprintf(stderr, "Matthews correlation: %5d %5d %5d %5d (11,00,10,01)\n",
	  nn[1][1],nn[0][0],nn[1][0], nn[0][1]);

  for(i = 0; i < len; i++) {
    sd += (pred[i]-ave)*(pred[i]-ave);
    c += (aans->nvec[i]-aans->ave_cn)*(pred[i]-ave);
    d += (aans->nvec[i] - pred[i])*(aans->nvec[i] - pred[i]);
  }
  sd = sqrt(sd/(double)aans->len);
  *cor = c/(aans->len*sd*aans->sd_cn);
  *deva = sqrt(d/(double)aans->len)/aans->sd_cn;
  return;
}

void check_accuracy_rwco(ans_t *aans, double pred[], double *cor, double *deva)
{
  int i;
  double ave=0.0,sd=0.0;
  double c=0.0,d=0.0;

  for(i = 0; i < aans->len; i++) {
    ave += pred[i];
  }
  ave /= (double)aans->len;

  for(i = 0; i < aans->len; i++) {
    sd += (pred[i]-ave)*(pred[i]-ave);
    c += (aans->ovec[i]-aans->ave_co)*(pred[i]-ave);
    d += (aans->ovec[i] - pred[i])*(aans->ovec[i] - pred[i]);
  }
  sd = sqrt(sd/(double)aans->len);
  *cor = c/(aans->len*sd*aans->sd_co);
  *deva = sqrt(d/(double)aans->len)/aans->sd_co;
  return;
}

