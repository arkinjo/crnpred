#ifndef __xpredm_h_
#define __xpredm_h_

#define NWET 20
#define WINDOW_HIGH 3
#define NI_HIGH ((NWET*5+1)*(2*WINDOW_HIGH+1))

typedef struct _wlist3_t {
  char *fw0;    /* file name of random unitary matrix */
  char *fwin;   /* file name of input random matrix */
  char *fwout_ss; /* file name of output matrix, SS */
  char *fwout_cn; /* file name of output matrix, CN */
  char *fwout_rwco; /* file name of output matrix, RWCO */
} wlist3_t;

int read_wlist3(wlist3_t wl[]);
int read_wlist3_lin(wlist3_t wl[]);
void read_big_wet(double wet[][5]);

void print_result(pssm_t *apssm, char psec[], double wsec[][3],
		  double pcn[], char scn2[], double prwco[]);

void print_result_gtop(pssm_t *apssm, char psec[], double wsec[][3],
		  double pcn[], char scn2[], double prwco[]);

void init_xpredm(pssm_t *apssm, char psec[], double wsec[][3],
		 double pcn[], char scn2[], double prwco[]);

void xsspred1(pssm_t *apssm, double xx[][NDIM], char *fwout, double pred[][3]);
void xcnpred1(pssm_t *apssm, double xx[][NDIM], char *fwout, double pred[]);
void xrwcopred1(pssm_t *apssm, double xx[][NDIM], char *fwout, double pred[]);

void lsspred1(pssm_t *apssm, char *fwout, double pred[][3]);
void lcnpred1(pssm_t *apssm, char *fwout, double pred[]);
void lrwcopred1(pssm_t *apssm, char *fwout, double pred[]);

void smooth_spred(int len, double pred[][3]);
void finalize_xpredm(int nwet, pssm_t *apssm, char psec[], double wsec[][3],
		     double pcn[], char scn2[], double prwco[]);

void make_high_input(int ind, int len, int nwet,
		     double psec[][MAXSEQ][3], double pcn[][MAXSEQ], 
		     double prwco[][MAXSEQ], double uu[]);
#endif
