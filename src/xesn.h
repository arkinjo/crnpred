#ifndef __xesn_h_
#define __xesn_h_

#include "blast.h"

#ifndef NDIM
#define NDIM 100
#endif

#ifndef BETA
#define BETA 0.5
#endif

/* coefficient for successive over-relaxation */
#ifndef OMEGA_SOR
#define OMEGA_SOR 1.4
#endif

#ifndef WINDOW_ESN
#define WINDOW_ESN 0
#endif

#define NI_ESN ((2*WINDOW_ESN+1)*21)

#define CONV_LIM 1.e-3   /* convergence limit */

#ifndef MAXBLOCK
#define MAXBLOCK 20
#endif

#define NWELEM (100*5000)   /* max # of elements of wmat */

#define MAXWFILES 300

#define NPAR_N (NDIM+NI_UNITS_N)
#define NPAR_O (NDIM+NI_UNITS_O)
#define NPAR_S (NDIM+NI_UNITS_S)

typedef struct _wmatsp_t {
  int i,j;
  double x;
} wmatsp_t;

typedef struct _wlist_t {
  char *fw0;
  char *fwin;
  char *fwout;
} wlist_t;

int read_wlist(char *file, wlist_t []);
void read_wmat(char *file, int nrow, int ncol, double *wmat[]);
void read_wmat1(char *file, int nrow, double wmat[]);
int read_wmatsp(char *file, wmatsp_t wmat[]);

void write_xx(char *file, int len, double xx[][NDIM]);
void read_xx(char *file, int len, double xx[][NDIM]);

void append_uxi(int niu, double uui[], double xxi[], double dalli[]);

void init_esn_state(int len, double win0[][NI_ESN], double uu[][NI_ESN], 
		    double hx[][NDIM], double xx[][NDIM]);

void get_stationary_state(int len, int nw0, wmatsp_t w0[], 
			  double hx[][NDIM], double xx[][NDIM]);

void get_stationary_state_jacobi(int len, int nw0, wmatsp_t w0[], 
				 double hx[][NDIM], double xx[][NDIM]);

#endif /* __xesn_h_ */
