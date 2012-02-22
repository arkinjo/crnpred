#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "blast.h"
#include "xesn.h"
#include "eprintf.h"
#include "chkaccu.h"

#include "xpredm.h"
/* 	$Id: xpredmsub.c,v 1.5 2006/03/27 02:49:36 akinjo Exp $	 */

#ifndef lint
static char vcid[] = "$Id: xpredmsub.c,v 1.5 2006/03/27 02:49:36 akinjo Exp $";
#endif /* lint */



static char WMATS[12];
static char *WMATS_LIN = "WMATS_LIN";
static char *BIGWET = "WMAT_ENS";

static char *append_env(char *base, char *name)
{
  char dest[500];

  dest[0] = '\0';
  if(base == NULL) 
    strcat(dest,".");
  else
    strcat(dest, base);

  strcat(dest,"/");
  strcat(dest,name);

  return estrdup(dest);
}

int read_wlist3(wlist3_t wl[])
{
  int n;
  char line[1001];
  char fw0[100],fwi[100],fswo[100], fcwo[100], frwo[100];
  char *crnpath;
  FILE *ic;
  char *fname;

  crnpath = getenv("CRNPRED_DIR");
  fprintf(stderr, "CRNPRED_DIR= %s\n", crnpath);
  sprintf(WMATS,"WMATS%d", NDIM);
  fname = append_env(crnpath, WMATS);

  if((ic = fopen(fname,"r")) == NULL) {
    eprintf("read_wlist: cannot open file: %s\n", fname);
  }

  n=0;
  while(fgets(line,1000,ic) != NULL) {
    if(line[0] == '#') continue;
    sscanf(line,"%s %s %s %s %s", fw0,fwi,fswo,fcwo, frwo);

    wl[n].fw0 = append_env(crnpath, fw0);
    wl[n].fwin = append_env(crnpath, fwi);
    wl[n].fwout_ss = append_env(crnpath, fswo);
    wl[n].fwout_cn = append_env(crnpath, fcwo);
    wl[n].fwout_rwco = append_env(crnpath, frwo);
    n++;
    if(n>MAXWFILES) {
      eprintf("read_wlist: too many weight files! (> %d)\n",MAXWFILES);
    }
  }

  fclose(ic);
  /*
  {
    int i;
    for(i=0; i< n; i++) {
      fprintf(stderr, "read_wlist: %s %s %s %s %s\n",
	      wl[i].fw0, wl[i].fwin,
	      wl[i].fwout_ss, wl[i].fwout_cn, wl[i].fwout_rwco);
    }
  }
  */

  return n;
}

void read_big_wet(double wet[][5])
{
  double *wop[NI_HIGH];
  char *fname, *crnpath;
  int i;

  crnpath = getenv("CRNPRED_DIR");
  fname = append_env(crnpath, BIGWET);
  for(i=0; i<NI_HIGH; i++) wop[i] = wet[i];
  read_wmat(fname, NI_HIGH, 5, wop);
  /* wet[i][3] (for CN), wet[i][4] (for RWCO) are not used. */

  return;
}

int read_wlist3_lin(wlist3_t wl[])
{
  int n;
  char line[1001];
  char fswo[100], fcwo[100], frwo[100];
  char *crnpath;
  FILE *ic;
  char *fname;

  crnpath = getenv("CRNPRED_DIR");
  fprintf(stderr, "CRNPRED_DIR= %s\n", crnpath);
  fname = append_env(crnpath, WMATS_LIN);

  if((ic = fopen(fname,"r")) == NULL) {
    eprintf("read_wlist: cannot open file: %s\n", fname);
  }

  n=0;
  while(fgets(line,1000,ic) != NULL) {
    if(line[0] == '#') continue;
    sscanf(line,"%s %s %s", fswo,fcwo, frwo);

    wl[n].fwout_ss = append_env(crnpath, fswo);
    wl[n].fwout_cn = append_env(crnpath, fcwo);
    wl[n].fwout_rwco = append_env(crnpath, frwo);
    n++;
    if(n>MAXWFILES) {
      eprintf("read_wlist: too many weight files! (> %d)\n",MAXWFILES);
    }
  }

  fclose(ic);

  return n;
}

void print_result(pssm_t *apssm, char psec[], double wsec[][3],
		  double pcn[], char scn2[], double prwco[])
{
  int i,j;


  printf(">prediction for: %s\n\n", apssm->fname);

  for(i = 0; i < apssm->len/60; i++) {
    printf("%-10s", "#");
    for(j=0; j< 60; j++) putchar(((j+1)%10) ? ' ' : '*');
    putchar('\n');

    printf("%-10s", "AA:");
    for(j=0; j< 60; j++) putchar(apssm->seq[i*60 + j]);
    putchar('\n');

    printf("%-10s", "SS:");
    for(j=0; j< 60; j++) putchar(psec[i*60 + j]);
    putchar('\n');

    printf("%-10s", "CN:");
    for(j=0; j< 60; j++) putchar(scn2[i*60 + j]);
    putchar('\n');
  }

  if(apssm->len % 60 != 0) {
    printf("%-10s", "#");
    for(i = 60*(apssm->len/60); i< apssm->len; i++) 
      putchar(((i+1)%10) ? ' ' : '*');
    putchar('\n');

    printf("%-10s", "AA:");
    for(i = 60*(apssm->len/60); i< apssm->len; i++) putchar(apssm->seq[i]);
    putchar('\n');

    printf("%-10s", "SS:");
    for(i = 60*(apssm->len/60); i< apssm->len; i++) putchar(psec[i]);
    putchar('\n');

    printf("%-10s", "CN:");
    for(i = 60*(apssm->len/60); i< apssm->len; i++) putchar(scn2[i]);
    putchar('\n');
  }
  printf("//\n\n");

  printf(">#   AA : SS P_H P_E P_C : CN     : RWCO\n");
  for(i = 0; i < apssm->len; i++) {
    double conf;
    double z,c0,c1,c2;
    const double beta = 3.0;

    c0 = exp(beta*(wsec[i][0]-wsec[i][2]));
    c1 = exp(beta*(wsec[i][1]-wsec[i][2]));
    c2 = 1.0;
    z = c0 + c1 + c2;
    c0 *= 100.0/z;
    c1 *= 100.0/z;
    c2 *= 100.0/z;
    switch (psec[i]) {
    case 'H':
      conf=c0;
      break;
    case 'E':
      conf=c1;
      break;
    default:
      conf=c2;
    }

    printf("%5d %c : %c  %3.0f %3.0f %3.0f : %c %4.0f : %4.0f\n", 
	   i+1, apssm->seq[i], psec[i], c0, c1, c2, scn2[i], pcn[i], prwco[i]);
    /*
    fprintf(stderr,"### %8.3f %8.3f %8.3f\n", wsec[i][0], wsec[i][1], wsec[i][2]);
    */
  }

}

void print_result_gtop(pssm_t *apssm, char psec[], double wsec[][3],
		       double pcn[], char scn2[], double prwco[])
{
  int i;

  printf("%-15s","CRNPRED_SS");
  for(i = 0; i < apssm->len; i++) 
    (psec[i] == 'C') ? putchar('c') : putchar(psec[i]);
  putchar('\n');

  printf("%-15s","CRNPRED_CN");
  for(i = 0; i < apssm->len; i++) 
    putchar(scn2[i]);
  putchar('\n');

  return;
}

void init_xpredm(pssm_t *apssm, char psec[], double wsec[][3],
		 double pcn[], char scn2[], double prwco[])
{
  int i;

  for(i=0; i< apssm->len; i++) {
    psec[i] = scn2[i] = '.';
    wsec[i][0] = wsec[i][1] = wsec[i][2] = 0.0;
    pcn[i] = prwco[i] = 0.0;
  }

  return;
}

void xcnpred1(pssm_t *apssm, double xx[][NDIM], char *fwout, double pred[])
{
  int i,j;
  double wout[NPAR_N];
  double uu0[NI_UNITS_N];
  double dall[NPAR_N];
  double eta;

  read_wmat1(fwout,NPAR_N,wout);
  for(i=0; i<apssm->len; i++) {
    pssm2iunit(HALF_WINDOW_N, i, apssm, uu0);
    append_uxi(NI_UNITS_N, uu0, xx[i+1], dall);
    eta = 0.0;
    for(j=0; j<NPAR_N; j++) {
      eta += wout[j]*dall[j];
    }
    pred[i] += eta;
  }
  return;
}

void xrwcopred1(pssm_t *apssm, double xx[][NDIM], char *fwout, double pred[])
{
  int i,j;
  double wout[NPAR_O];
  double uu0[NI_UNITS_O];
  double dall[NPAR_O];
  double eta;

  read_wmat1(fwout,NPAR_O,wout);
  for(i=0; i<apssm->len; i++) {
    pssm2iunit(HALF_WINDOW_O, i, apssm, uu0);
    append_uxi(NI_UNITS_O, uu0, xx[i+1], dall);
    eta = 0.0;
    for(j=0; j<NPAR_O; j++) {
      eta += wout[j]*dall[j];
    }
    pred[i] += eta;
  }
  return;
}

void xsspred1(pssm_t *apssm, double xx[][NDIM], char *fwout, double pred[][3])
{
  int i,j;
  double wout[NPAR_S][3],*wop[NPAR_S];
  double uu0[NI_UNITS_S];
  double dall[NPAR_S];
  double eta0,eta1,eta2;
  double mu[3],z;
  double tpred[MAXSEQ][3];
  char psec[MAXSEQ];
  int ind;

  for(i=0; i<NPAR_S; i++) wop[i] = wout[i];
  read_wmat(fwout,NPAR_S,3,wop);

  for(i=0; i<apssm->len; i++) {
    pssm2iunit(HALF_WINDOW_S, i, apssm, uu0);
    append_uxi(NI_UNITS_S, uu0, xx[i+1], dall);
    eta0 = eta1 = eta2 = 0.0;
    for(j=0; j<NPAR_S; j++) {
      eta0 += wout[j][0]*dall[j];
      eta1 += wout[j][1]*dall[j];
      eta2 += wout[j][2]*dall[j];
    }
    mu[0] = exp(eta0);
    mu[1] = exp(eta1);
    mu[2] = exp(eta2);
    z = mu[0] + mu[1] + mu[2];
    tpred[i][0] = eta0;
    tpred[i][1] = eta1;
    tpred[i][2] = eta2;
  }
  secmax(apssm->len, tpred, psec);
  for(i=0; i<apssm->len; i++) {
    if (psec[i]=='H') ind=0;
    else if(psec[i] == 'E') ind=1;
    else ind = 2;
    pred[i][ind] += 1.0;
    /*
    pred[i][(ind+1)%3] -= 1.0;
    pred[i][(ind+2)%3] -= 1.0;
    */
  }

  return;
}

void finalize_xpredm(int nwet, pssm_t *apssm, char psec[], double wsec[][3],
		     double pcn[], char scn2[], double prwco[])
{
  int i, ind;
  
  for(i=0; i< apssm->len; i++) {
    pcn[i] *= log(apssm->len) /(double)nwet;
    ind = strchr(amino1, apssm->seq[i]) - amino1;
    scn2[i] = (cnave2[ind] < pcn[i]) ? 'E' : 'B';

    prwco[i] *= (double) (apssm->len) /(double)nwet;
    /*
    wsec[i][0] /= (double)nwet;
    wsec[i][1] /= (double)nwet;
    wsec[i][2] /= (double)nwet;
    */
  }

  smooth_spred(apssm->len, wsec);
  secmax(apssm->len, wsec, psec);

  return;
}

void smooth_spred(int len, double pred[][3])
{
  int i,j;
  double tpred[MAXSEQ][3];

  for(i=0; i<len; i++) {
    for(j=0;j<3;j++)
      tpred[i][j] = 0.0;
  }

  for(j=0;j<3;j++) {
    tpred[0][j] = 0.5*(pred[0][j]+pred[1][j]);
    tpred[len-1][j] = 0.5*(pred[len-1][j]+pred[len-2][j]);
  }
  for(i=1;i<len-1;i++) {
    for(j=0;j<3;j++)
      tpred[i][j] = 
	0.25*(pred[i-1][j]+pred[i+1][j]) + 0.5*pred[i][j];
  }

  for(i=0; i<len; i++) {
    for(j=0;j<3;j++)
      pred[i][j] = tpred[i][j];
  }
			       
  return;
}

void lcnpred1(pssm_t *apssm, char *fwout, double pred[])
{
  int i,j;
  double wout[NI_UNITS_N];
  double uu0[NI_UNITS_N];
  double eta;

  read_wmat1(fwout,NI_UNITS_N,wout);
  for(i=0; i<apssm->len; i++) {
    pssm2iunit(HALF_WINDOW_N, i, apssm, uu0);
    eta = 0.0;
    for(j=0; j<NI_UNITS_N; j++) {
      eta += wout[j]*uu0[j];
    }
    pred[i] += eta;
  }
  return;
}

void lrwcopred1(pssm_t *apssm, char *fwout, double pred[])
{
  int i,j;
  double wout[NI_UNITS_O];
  double uu0[NI_UNITS_O];
  double eta;

  read_wmat1(fwout,NI_UNITS_O,wout);
  for(i=0; i<apssm->len; i++) {
    pssm2iunit(HALF_WINDOW_O, i, apssm, uu0);
    eta = 0.0;
    for(j=0; j<NI_UNITS_O; j++) {
      eta += wout[j]*uu0[j];
    }
    pred[i] += eta;
  }
  return;
}

void lsspred1(pssm_t *apssm, char *fwout, double pred[][3])
{
  int i,j;
  double wout[NI_UNITS_S][3],*wop[NI_UNITS_S];
  double uu0[NI_UNITS_S];
  double eta0,eta1,eta2;
  double mu[3],z;
  double tpred[MAXSEQ][3];
  char psec[MAXSEQ];
  int ind;

  for(i=0; i<NI_UNITS_S; i++) wop[i] = wout[i];

  read_wmat(fwout,NI_UNITS_S,3,wop);

  for(i=0; i<apssm->len; i++) {
    pssm2iunit(HALF_WINDOW_S, i, apssm, uu0);
    eta0 = eta1 = eta2 = 0.0;
    for(j=0; j<NI_UNITS_S; j++) {
      eta0 += wout[j][0]*uu0[j];
      eta1 += wout[j][1]*uu0[j];
      eta2 += wout[j][2]*uu0[j];
    }
    mu[0] = exp(eta0);
    mu[1] = exp(eta1);
    mu[2] = exp(eta2);
    z = mu[0] + mu[1] + mu[2];
    tpred[i][0] = eta0;
    tpred[i][1] = eta1;
    tpred[i][2] = eta2;
  }
  secmax(apssm->len, tpred, psec);
  for(i=0; i<apssm->len; i++) {
    if (psec[i]=='H') ind=0;
    else if(psec[i] == 'E') ind=1;
    else ind = 2;
    pred[i][ind] += 1.0;
    /*
    pred[i][(ind+1)%3] -= 1.0;
    pred[i][(ind+2)%3] -= 1.0;
    */
  }

  return;
}

void make_high_input(int ind, int len, int nwet,
		     double psec[][MAXSEQ][3], double pcn[][MAXSEQ], 
		     double prwco[][MAXSEQ], double uu[])
{
  int i,j,k;
  int is,ie,n;

  is = ind - WINDOW_HIGH;
  ie = ind + WINDOW_HIGH;
  for(n = 0; n < NI_HIGH; n++) uu[n] = 0.0;

  n = 0;
  for(i=is; i <= ie; i++) {
    if(i<0 || i>= len) {
      for(j=0; j<nwet*5; j++) {
	uu[n++] = 0.0;
      }
      uu[n++] = 1.0;
    }
    else {
      for(k = 0; k < nwet; k++) {
	uu[n++] = psec[k][i][0];
	uu[n++] = psec[k][i][1];
	uu[n++] = psec[k][i][2];
	uu[n++] = pcn[k][i];
	uu[n++] = prwco[k][i];
      }
      uu[n++] = (i==ind) ? 1.0 : 0.0;
    }
  }

  return;
}
