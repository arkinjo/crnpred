#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blast.h"
#include "xesn.h"
#include "eprintf.h"
#include "chkaccu.h"

#include "xpredm.h"
/* 	$Id: xpredm.c,v 1.4 2006/03/13 08:47:44 akinjo Exp $	 */

#ifndef lint
static char vcid[] = "$Id: xpredm.c,v 1.4 2006/03/13 08:47:44 akinjo Exp $";
#endif /* lint */


static void big_pred(int len, double wet[][5],
		     double psec[][MAXSEQ][3], double pcn[][MAXSEQ],
		     double prwco[][MAXSEQ],
		     double qsec[][3], double qcn[], double qrwco[])
{
  double uu[NI_HIGH];
  int i,j;
  
  for(i=0; i< len; i++) {
    qsec[i][0] = qsec[i][1] = qsec[i][2] = qcn[i] = qrwco[i] = 0.0;
    make_high_input(i, len, NWET, psec, pcn, prwco, uu);
    for(j=0; j<NI_HIGH; j++) {
      qsec[i][0] += wet[j][0]*uu[j];
      qsec[i][1] += wet[j][1]*uu[j];
      qsec[i][2] += wet[j][2]*uu[j];
    }
    /* simple summation is better for CN and RWCO */
    for(j=0; j<NWET; j++) {
      qcn[i] += pcn[j][i];
      qrwco[i] += prwco[j][i];
    }
  }
  
  return;
}

int main(int argc, char *argv[])
{
  int i,j,n;
  pssm_t apssm;
  wlist3_t wetl[MAXWFILES];
  int nwet;
  char csec[MAXSEQ],scn2[MAXSEQ];
  double psec[MAXSEQ][3];
  double pcn[MAXSEQ], prwco[MAXSEQ];
  double qsec[NWET][MAXSEQ][3];
  double qcn[NWET][MAXSEQ], qrwco[NWET][MAXSEQ];

  static double hx[MAXSEQ][NDIM];
  static double xx[MAXSEQ][NDIM];
  double uu[MAXSEQ][NI_ESN];
  static double win0[NDIM][NI_ESN], *win0p[NDIM];
  static wmatsp_t w0[NWELEM];
  int nw0;
  static double big_wet[NI_HIGH][5];

  if(argc<2) {
    eprintf("Usage: xpredm pssm\n");
  }

  read_pssm(argv[1],&apssm);
  apssm.fname = estrdup(argv[1]);
  read_big_wet(big_wet);

  nwet = read_wlist3(wetl);
  if(nwet != NWET) {
    eprintf("number of weights files incorrect: expected %5d, got %5d\n",
	    NWET, nwet);
  }
  for(i= 0; i < NDIM; i++) {
    win0p[i] = &win0[i][0];
  }

  init_xpredm(&apssm, csec, psec, pcn, scn2, prwco);
  for(n = 0; n < NWET; n++) {
    for(i=0; i < apssm.len; i++) {
      qsec[n][i][0] = qsec[n][i][1] = qsec[n][i][2] = 0.0;
      qcn[n][i] = qrwco[n][i] = 0.0;
    }
  }

  for(i=0; i < apssm.len; i++) {
    pssm2iunit(WINDOW_ESN, i, &apssm, uu[i]);
  }

  for(n=0; n<nwet; n++) {
    nw0 = read_wmatsp(wetl[n].fw0, w0);
    read_wmat(wetl[n].fwin, NDIM, NI_ESN, win0p);
    init_esn_state(apssm.len, win0, uu, hx, xx);
    get_stationary_state(apssm.len, nw0, w0, hx, xx);

    xsspred1(&apssm, xx, wetl[n].fwout_ss, qsec[n]);
    xcnpred1(&apssm, xx, wetl[n].fwout_cn, qcn[n]);
    xrwcopred1(&apssm, xx, wetl[n].fwout_rwco, qrwco[n]);
  }

  big_pred(apssm.len, big_wet, qsec, qcn, qrwco, psec, pcn, prwco);

  finalize_xpredm(nwet, &apssm, csec, psec, pcn, scn2, prwco);

#ifndef GTOP
  print_result(&apssm, csec, psec, pcn, scn2, prwco);
#else
  print_result_gtop(&apssm, csec, psec, pcn, scn2, prwco);
#endif

  free(apssm.fname);
  return 0;
}
