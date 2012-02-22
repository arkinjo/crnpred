#include <stdio.h>
#include "blast.h"

/* 	$Id: secdp.c,v 1.2 2006/02/20 08:35:37 akinjo Exp $	 */

#ifndef lint
static char vcid[] = "$Id: secdp.c,v 1.2 2006/02/20 08:35:37 akinjo Exp $";
#endif /* lint */

static const double scagap = 1.0;
static  char *sec = "HEC";

static double gap[3][3] = 
/* { */
/*   {       1.0113036 ,      -5.8465444 ,      -1.6370586 , }, */
/*   {      -4.2604538 ,       1.3119079 ,      -0.9302283 , }, */
/*   {      -1.6658889 ,      -0.9093533 ,       0.6393491 , }, */
/* }; */
{
  {      -1.1949517 ,      -8.4615910 ,      -3.5438586 , },
  {      -6.8755005 ,      -1.7119301 ,      -3.2458195 , },
  {      -3.5726889 ,      -3.2249445 ,      -0.9679954 , },
};



double mymax3(double st[], int *ind) {
  if(st[2] >= st[1] && st[2] >= st[0]) {
    *ind = 2;
    return st[2];
  }
  else if(st[1] >= st[2] && st[1] >= st[0]) {
    *ind = 1;
    return st[1];
  }

  *ind = 0;
  return st[0];
}

void secdp(int len, double y[][3], char psec[])
{
  double t[MAXSEQ][3];
  int bt[MAXSEQ][3];
  double st[3];


  int i,j,ind;
  double tmax;

  for(i=0; i < len; i++) {
    for(j=0; j<3; j++) {
      t[i][j] = 0.0;
      bt[i][j] = 0;
    }
  }

  for(j=0;j<3;j++) t[0][j] = y[0][j];

  for(i = 1; i<len; i++) {
    for(j=0; j<3; j++) {
      st[0] = t[i-1][0] + gap[0][j]*scagap;
      st[1] = t[i-1][1] + gap[1][j]*scagap;
      st[2] = t[i-1][2] + gap[2][j]*scagap;
      tmax = mymax3(st, &ind);
      t[i][j] = tmax + y[i][j];
      bt[i][j] = ind;
    }
  }

/*   for(i=0; i< len; i++) { */
/*     printf ("%5d : ", i+1); */
/*     for(j=0; j<3; j++) { */
/*       printf("%8.3f %5d %8.3f: ", t[i][j],bt[i][j], y[i][j]); */
/*     } */
/*     printf("\n"); */
/*   } */

  tmax = mymax3(t[len-1], &ind);
  psec[len-1] = sec[ind];
  for(i=len-1; i> 0; i--) {
    psec[i-1] = sec[bt[i][ind]];
    ind = bt[i][ind];
  }

  return;
}  

void secmax(int len, double y[][3], char psec[])
{
  int i,ind;
  double tmax;

  for(i=0; i< len; i++) {
    tmax = mymax3(y[i],&ind);
    psec[i] = sec[ind];
  }
  return;
}
