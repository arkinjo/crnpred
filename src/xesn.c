#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "xesn.h"
#include "eprintf.h"
/* 	$Id: xesn.c,v 1.2 2006/02/20 08:38:12 akinjo Exp $	 */

#ifndef lint
static char vcid[] = "$Id: xesn.c,v 1.2 2006/02/20 08:38:12 akinjo Exp $";
#endif /* lint */

#define MAXITER 1000
#define WIDENT 16

#define PERIODIC_BOUNDARY 1
#undef PERIODIC_BOUNDARY  /* comment out if you want periodic boundary condition */

#ifdef LINEAR
#define NET_FUNC(x) (x)
#else
#define NET_FUNC(x) (tanh(x))
#endif

int read_wlist(char *fname, wlist_t wl[])
{
  int n;
  char line[1001];
  char fw0[200],fwi[200],fwo[200];

  FILE *ic;

  if((ic = fopen(fname,"r")) == NULL) {
    eprintf("read_wlist: cannot open file: %s\n", fname);
  }

  n=0;
  while(fgets(line,1000,ic) != NULL) {
    if(line[0] == '#') continue;
    sscanf(line,"%s %s %s", fw0,fwi,fwo);
    wl[n].fw0 = estrdup(fw0);
    wl[n].fwin = estrdup(fwi);
    wl[n].fwout = estrdup(fwo);
    n++;
    if(n>MAXWFILES) {
      eprintf("read_wlist: too many weight files! (> %d)\n",MAXWFILES);
    }
  }

  fclose(ic);

  return n;
}

void read_wmat(char *fname, int nrow, int ncol, double *wmat[])
{
  int i,j;
  int line_len = WIDENT*2000;
  char line[WIDENT*2000+1];
  char subline[20];
  FILE *ic;

  for(i=0; i < nrow; i++)
    for(j=0; j < ncol; j++)
      wmat[i][j] = 0.0;

  if((ic = fopen(fname,"r")) == NULL) {
    eprintf("read_wmat: cannot open file: %s\n", fname);
  }

  i = 0;
  while(fgets(line,line_len,ic) != NULL) {
    for(j=0; j < ncol; j++) {
      strncpy(subline,line+WIDENT*j,WIDENT);
      wmat[i][j] = atof(subline);
      if(j>ncol) {
	eprintf("read_wmat: matrix col. size exceeded > %5d\n", ncol);
      }
    }
    i++;
    if(i>nrow) {
      eprintf("read_wmat: matrix row size exceeded > %5d\n", nrow);
    }
  }

  fclose(ic);
}

#undef WIDENT
void read_wmat1(char *fname, int nrow, double wmat[])
{
  int i;
  int line_len = 200;
  char line[201];
  FILE *ic;

  for(i=0; i < nrow; i++)
    wmat[i] = 0.0;

  if((ic = fopen(fname,"r")) == NULL) {
    eprintf("read_wmat: cannot open file: %s\n", fname);
  }

  i = 0;
  while(fgets(line,line_len,ic) != NULL) {
    wmat[i] = atof(line);
    i++;
    if(i>nrow) {
      eprintf("read_wmat: matrix row size exceeded > %5d\n", nrow);
    }
  }

  fclose(ic);
}

/* read sparse wmat */
int read_wmatsp(char *fname, wmatsp_t wmat[])
{
  int i,j,n;
  double x;
  char line[101];
  FILE *ic;

  for(n=0; n < NWELEM; n++) {
    wmat[n].i = wmat[n].j = -1;
    wmat[n].x = 0.0;
  }

  if((ic = fopen(fname,"r")) == NULL) {
    eprintf("read_wmatsp: cannot open file: %s\n", fname);
  }

  n = 0;
  while(fgets(line,100,ic) != NULL) {
    sscanf(line,"%d %d %lf", &i, &j, &x);
    wmat[n].i = i;
    wmat[n].j = j;
    wmat[n].x = x;
    n++;
    if(n>NWELEM) {
      eprintf("read_wmatsp: matrix row size exceeded > %5d\n", NWELEM);
    }
  }

  fclose(ic);

  return n;
}

void init_esn_state(int len, double win0[][NI_ESN], double uu[][NI_ESN], 
		    double hx[][NDIM], double xx[][NDIM])
{
  int i,j,k;
  double waf;

  for(i=0; i <= len+1; i++) {
    for(j=0; j < NDIM; j++) {
      xx[i][j] = 0.0;
      hx[i][j] = 0.0;
    }
  }

  for(i=1; i <= len; i++) {
    for(j=0; j < NDIM; j++) {
      waf = 0.0;
      for(k = 0; k < NI_ESN; k++) {
	waf += win0[j][k]*uu[i-1][k];
      }
      hx[i][j] = waf;
      xx[i][j] = NET_FUNC(waf);
    }
  }

#ifdef PERIODIC_BOUNDARY
  fprintf(stderr, "periodic boundary condition imposed.\n");
  for(j=0; j < NDIM; j++) {
    xx[0][j] = xx[len][j];
    xx[len+1][j] = xx[1][j];
    hx[0][j] = hx[len][j];
    hx[len+1][j] = hx[1][j];
  }
#endif
  return;
}

void get_stationary_state(int len, int nw0, wmatsp_t w0[], 
			  double hx[][NDIM], double xx[][NDIM])
{
  int iter,i,j, k, n;
  double rmse,dev;
  static double wa[NDIM];
  static double xxt[MAXSEQ][NDIM];

  for(iter = 0; iter < MAXITER; iter++) {
    for(i = 1; i <= len; i++) {
      for(j=0; j < NDIM; j++) {
	wa[j] = hx[i][j];
      }
      for(n=0; n < nw0; n++) {
	j = w0[n].i;
	k = w0[n].j;
	wa[j] += w0[n].x*BETA*(xxt[i-1][k]+xx[i+1][k]);
      }
      for(j=0; j < NDIM; j++) {
	xxt[i][j] = xx[i][j] + OMEGA_SOR*(NET_FUNC(wa[j]) - xx[i][j]);
      }
    }
#ifdef PERIODIC_BOUNDARY
    for(j=0; j < NDIM; j++) {
      xxt[0][j] = xxt[len][j];
      xxt[len+1][j] = xxt[1][j];
    }
#endif
    for(i = len; i >= 1; i--) {
      for(j=0; j < NDIM; j++) {
	wa[j] = hx[i][j];
      }
      for(n=0; n < nw0; n++) {
	j = w0[n].i;
	k = w0[n].j;
	wa[j] += w0[n].x*BETA*(xxt[i-1][k]+xx[i+1][k]);
      }
      for(j=0; j < NDIM; j++) {
	xx[i][j] = xxt[i][j] + OMEGA_SOR*(NET_FUNC(wa[j]) - xxt[i][j]);
      }
    }
#ifdef PERIODIC_BOUNDARY
    for(j=0; j < NDIM; j++) {
      xx[0][j] = xx[len][j];
      xx[len+1][j] = xx[1][j];
    }
#endif
    /* check convergence */
    rmse = 0.0;
    for(i = 1; i <= len; i++) {
      for(j=0; j < NDIM; j++) {
	dev = xxt[i][j] - xx[i][j];
	rmse += dev*dev;
      }
    }
    rmse = sqrt(rmse/(len*NDIM));
    if(rmse < CONV_LIM) break;
#ifdef DEBUG
    fprintf(stderr,"get_stationary_state: %5d : rmse = %15.7e\n", iter, rmse);
#endif
  }

  fprintf(stderr,"get_stationary_state: %5d : rmse = %15.7e\n", iter, rmse);
  if(iter == MAXITER) {
    weprintf("get_stationary_state: did not converge!\n");
  }

  return;
}


/* this is for testing. very slow. */
void get_stationary_state_jacobi(int len, int nw0, wmatsp_t w0[], 
				 double hx[][NDIM], double xx[][NDIM])
{
  int iter,i,j, k,n;
  double rmse,dev;
  static double wa[NDIM];
  static double xxt[MAXSEQ][NDIM];

  for(iter = 0; iter < MAXITER; iter++) {
    for(i = 0; i <= len+1; i++) {
      for(j=0; j < NDIM; j++) {
	xxt[i][j] = xx[i][j];
      }
    }
    for(i = 1; i <= len; i++) {
      for(j=0; j < NDIM; j++) {
	wa[j] = hx[i][j];
      }
      for(n=0; n < nw0; n++) {
	j = w0[n].i;
	k = w0[n].j;
	wa[j] += w0[n].x*BETA*(xxt[i-1][j]+xx[i+1][j]);
      }
      for(j=0; j < NDIM; j++) {
	xx[i][j] = NET_FUNC(wa[j]);
      }
    }
#ifdef PERIODIC_BOUNDARY
    for(j=0; j < NDIM; j++) {
      xx[0][j] = xx[len][j];
      xx[len+1][j] = xx[1][j];
    }
#endif

    /* check convergence */
    rmse = 0.0;
    for(i = 1; i <= len; i++) {
      for(j=0; j < NDIM; j++) {
	dev = xxt[i][j] - xx[i][j];
	rmse += dev*dev;
      }
    }
    rmse = sqrt(rmse/(len*NDIM));
    if(rmse < CONV_LIM) break;
#ifdef DEBUG
    fprintf(stderr,"get_stationary_state_jacobi: %5d: rmse = %15.7e\n", iter, rmse);
#endif
  }

  fprintf(stderr,"get_stationary_state_jacobi: %5d: rmse = %15.7e\n", iter, rmse);
  if(iter == MAXITER) {
    weprintf("get_stationary_state_jacobi: did not converge!\n");
  }

  return;
}


void write_xx(char *fname, int len, double xx[][NDIM])
{
  FILE *oc;
  const size_t nsize = (size_t) (len * NDIM);

  if((oc = fopen(fname,"w")) == NULL) {
    eprintf("write_xx: cannot open file: %s\n", fname);
  }

  if(fwrite(xx[1], sizeof(double), nsize, oc) < nsize) {
    eprintf("write_xx: error occurred during writing: %s\n", fname);
  }
  fclose(oc);

  return;
}

void read_xx(char *fname, int len, double xx[][NDIM])
{
  FILE *ic;
  const size_t nsize = (size_t) (len * NDIM);

  if((ic = fopen(fname,"r")) == NULL) {
    eprintf("read_xx: cannot open file: %s\n", fname);
  }

  if(fread(xx[1], sizeof(double), nsize, ic)< nsize) {
    eprintf("read_xx: error occurred during reading: %s\n", fname);
  }
  fclose(ic);

  return;
}

void append_uxi(int niu, double uui[], double xxi[], double dalli[])
{
  int i,ind;
  
  ind = 0;
  for(i = 0; i < niu; i++) {
    dalli[ind++] = uui[i];
  }
  for(i = 0; i < NDIM; i++) {
    dalli[ind++] = xxi[i];
  }

  return;
}
