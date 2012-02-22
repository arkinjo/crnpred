#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blast.h"
#include "xesn.h"
#include "eprintf.h"
#include "chkaccu.h"

#include "xpredm.h"

/* 	$Id: lpredm.c,v 1.2 2006/02/20 08:33:52 akinjo Exp $	 */

#ifndef lint
static char vcid[] = "$Id: lpredm.c,v 1.2 2006/02/20 08:33:52 akinjo Exp $";
#endif /* lint */

/** Linear prediction. */
int main(int argc, char *argv[])
{
  int i,n;
  pssm_t apssm;
  wlist3_t wetl[MAXWFILES];
  int nwet;
  char psec[MAXSEQ],scn2[MAXSEQ];
  double wsec[MAXSEQ][3];
  double pcn[MAXSEQ], prwco[MAXSEQ];

  double uu[MAXSEQ][NI_ESN];

  if(argc<2) {
    eprintf("Usage: xpredm pssm\n");
  }

  read_pssm(argv[1],&apssm);
  apssm.fname = estrdup(argv[1]);

  nwet = read_wlist3_lin(wetl);

  init_xpredm(&apssm, psec, wsec, pcn, scn2, prwco);

  for(i=0; i < apssm.len; i++) {
    pssm2iunit(WINDOW_ESN, i, &apssm, uu[i]);
  }

  for(n=0; n<nwet; n++) {
    lsspred1(&apssm, wetl[n].fwout_ss, wsec);
    lcnpred1(&apssm, wetl[n].fwout_cn, pcn);
    lrwcopred1(&apssm, wetl[n].fwout_rwco, prwco);
  }

  finalize_xpredm(nwet, &apssm, psec, wsec, pcn, scn2, prwco);

#ifndef GTOP
  print_result(&apssm, psec, wsec, pcn, scn2, prwco);
#else
  print_result_gtop(&apssm, psec, wsec, pcn, scn2, prwco);
#endif

  free(apssm.fname);
  return 0;
}
