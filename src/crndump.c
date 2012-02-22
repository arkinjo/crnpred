/* 	$Id: crndump.c,v 1.2 2006/02/20 08:33:09 akinjo Exp $	 */

#ifndef lint
static char vcid[] = "$Id: crndump.c,v 1.2 2006/02/20 08:33:09 akinjo Exp $";
#endif /* lint */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "blast.h"
#include "xesn.h"
#include "eprintf.h"

int main(int argc, char *argv[])
{
  int i;
  pssm_t apssm;
  FILE *op;
  char *dump_file;

  double uu[MAXSEQ][NI_ESN];  
  double xx[MAXSEQ][NDIM];
  double hx[MAXSEQ][NDIM];
  double win0[NDIM][NI_ESN], *win0p[NDIM];
  wmatsp_t w0[NWELEM];
  int nw0;

  if(argc<5) {
    eprintf("Usage: train wmat win pssm dump.out\n");
  }

  fprintf(stderr,"# K1,N,L= %d,%d,%d\n", NI_ESN,NDIM, 1);
  for(i= 0; i < NDIM; i++) {
    win0p[i] = &win0[i][0];
  }

  nw0 = read_wmatsp(argv[1], w0);
  read_wmat(argv[2], NDIM, NI_ESN, win0p);
  read_pssm(argv[3], &apssm);
  dump_file = estrdup(argv[4]);

  for(i=0; i < apssm.len; i++) {
    pssm2iunit(WINDOW_ESN, i, &apssm, uu[i]);
  }
  init_esn_state(apssm.len, win0, uu, hx,xx);
  get_stationary_state(apssm.len, nw0, w0, hx, xx);

  op = fopen(dump_file,"w");
  write_xx(dump_file, apssm.len, xx);

#ifdef DEBUG
  for(i=1; i <= apssm.len; i++) {
    int j;
    for(j=0; j< NDIM; j++) {
      printf(" %15.7e", xx[i][j]);
    }
    printf("\n");
  }
#endif
  return 0;
}
