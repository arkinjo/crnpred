#ifndef __blast_h_
#define __blast_h_

/* 	$Id: blast.h,v 1.2 2006/02/20 08:29:20 akinjo Exp $	 */

/* data structure for PSI-BLAST profile (PSSM) */

/* maximum sequence length */
#define MAXSEQ 1000
/* secondary structure */
#ifndef HALF_WINDOW_S
#define HALF_WINDOW_S 9
#endif

/* contact number */
#ifndef HALF_WINDOW_N
#define HALF_WINDOW_N 9
#endif

/* residue-wise contact order */
#ifndef HALF_WINDOW_O
#define HALF_WINDOW_O 26
#endif

#if NI_LEVEL==3
# define NI_UNITS_S ((2*HALF_WINDOW_S+1)*21 + 40)
# define NI_UNITS_N ((2*HALF_WINDOW_N+1)*21 + 40)
# define NI_UNITS_O ((2*HALF_WINDOW_O+1)*21 + 40)
#elif (NI_LEVEL==1 || NI_LEVEL==2)
# define NI_UNITS_S ((2*HALF_WINDOW_S+1)*21 + 20)
# define NI_UNITS_N ((2*HALF_WINDOW_N+1)*21 + 20)
# define NI_UNITS_O ((2*HALF_WINDOW_O+1)*21 + 20)
#else
# define NI_UNITS_S ((2*HALF_WINDOW_S+1)*21)
# define NI_UNITS_N ((2*HALF_WINDOW_N+1)*21)
# define NI_UNITS_O ((2*HALF_WINDOW_O+1)*21)
#endif

#define NO_UNITS_S (3)  
#define NO_UNITS_N (1)  
#define NO_UNITS_O (1)  

typedef struct _pssm {
  int len;
  char *fname;
  char seq[MAXSEQ];
  double prof[MAXSEQ][21];
  double compo[20]; /* AA composition */

  /* site-dependent composition weighted by seq. sep.*/
  double wlcompo[MAXSEQ][20]; 
} pssm_t;

typedef struct _ans {
  int len;
  int nh, ne, nc;
  double ave_cn,sd_cn;
  double ave_co,sd_co;
  char ss[MAXSEQ];
  double svec[MAXSEQ][NO_UNITS_S]; 
  double nvec[MAXSEQ];
  double ovec[MAXSEQ];
} ans_t;


#define MAXFILES 720

typedef struct _mydata {
  int nent;
  long ndat;
  pssm_t pssm[MAXFILES];
  ans_t ans[MAXFILES];
  char *xxfile[MAXFILES];
} mydata_t;

extern const char amino1[]; /* one-letter amino acid code */

extern void read_pssm(const char *file, pssm_t *prof);
extern void read_ans(const char *file, ans_t *ans);
extern void pssm2iunit(const int window, const int ind, const pssm_t *apssm, 
                       double iunit[]);
extern int read_file_list(char list[], mydata_t *adata);

double mymax3(double st[], int *ind);
extern void secdp(int len, double y[][3], char psec[]);
extern void secmax(int len, double y[][3], char psec[]);
#endif /* __blast_h_ */
