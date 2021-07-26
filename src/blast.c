/* 	$Id: blast.c,v 1.2 2006/02/20 08:31:42 akinjo Exp $	 */

#ifndef lint
static char vcid[] = "$Id: blast.c,v 1.2 2006/02/20 08:31:42 akinjo Exp $";
#endif /* lint */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>


#include "eprintf.h"
#include "blast.h"

const char amino1[] = "ARNDCQEGHILKMFPSTWYV";

enum {MAXBUF = 201, MAXLINE=300};

/* undef this if PSI score matrix is to be used instead of AA frequency.*/
/*#define FREQPROF 0*/

#ifdef FREQPROF 
/* reading __frequency__ profile */
#define WIDENT 4
#define OFFENT 70
#define SCAENT 0.01
#else
/* reading psi __score__ profile */
#define WIDENT 3
#define OFFENT 9
#define SCAENT 1.0
#endif

static void get_compo(pssm_t *x)
{
  int i,ind;

  for(i=0; i<20; i++) x->compo[i] = 0.0;
  for(i=0; i<x->len; i++) {
    ind = strchr(amino1, x->seq[i]) - amino1;
    if(ind>=0 && ind < 20) {
      x->compo[ind] += 1.0;
    } else {
      x->compo[0] += 1.0;
    }
  }

  for(i=0; i<20; i++) x->compo[i] /= x->len;

  return;
}

static void get_wlcompo(pssm_t *x)
{
  int i,j,ind;
  double lsep, lsep0;

  for(i=0; i<x->len; i++) 
    for(j=0; j<20; j++) 
      x->wlcompo[i][j] = 0.0;

  for(i=0; i<x->len; i++) {
    lsep0 = 0.0;
    for(j=0; j<x->len; j++) {
      lsep = abs(i-j);
      if (lsep < 2) continue;
      lsep = log(lsep);
      lsep0 += lsep;
      ind = strchr(amino1, x->seq[j]) - amino1;
      if(ind>=0) {
	x->wlcompo[i][ind] += lsep;
      }
    }
    for(j=0; j<20; j++) 
      if(x->compo[j]>0.0) {
	x->wlcompo[i][j] = (x->wlcompo[i][j]/(lsep0 * x->compo[j]));
      }/*
      else {
	x->wlcompo[i][j] = 1.0;
	}*/
  }


  return;
}


/* read PSI-Blast profile */
void read_pssm(const char *fname, pssm_t *x)
{
  FILE *ic;
  char line[MAXBUF]; 
  char subline[MAXBUF];
  int i,j;
  char aa;

  for(i = 0; i < MAXSEQ; i++) {
    for(j = 0; j < 21; j++) {
      x->prof[i][j] = 0.0;
    }
  }
      
  if((ic = fopen(fname,"r")) == NULL) {
    eprintf("read_pssm: cannot open profile :%s:\n", fname);
  }

  i=0;
  while(fgets(line,MAXBUF,ic)!=NULL) {
    if(strlen(line)< 150) continue;
    if(!isdigit(line[4])) continue;

    aa = x->seq[i] = line[6];

    for(j = 0; j < 20; j++) {
      strncpy(subline,line+OFFENT+j*WIDENT,WIDENT);
      x->prof[i][j] = (SCAENT*atof(subline)); 
    }
    i++;
    if(i>= MAXSEQ) {
      eprintf("read_pssm: maximum sequence length (%d) exceeded in %s\n",
	      MAXSEQ, fname);
    }
  }

  fclose(ic);

  x->len = i;

  get_compo(x);
  get_wlcompo(x);

  return;
}

#undef WIDENT 
#undef OFFENT 

void read_ans(const char *fname, ans_t *ans)
{
  FILE *ic;
  char line[MAXBUF]; 
  int i;
  int ires, len;
  char aa,sec;
  int nh,ne,nc;
  double acn,scn,aco,sco;
  double cn,rwco;

  if((ic = fopen(fname,"r")) == NULL) {
    eprintf("read_ans: cannot open profile :%s:\n", fname);
  }

  ans->len = 0;
  for(i = 0; i < MAXSEQ; i++) {
    ans->svec[i][0] = ans->svec[i][1] = ans->svec[i][2] = 0.0;
  }

  i=0;
  while(fgets(line,MAXBUF,ic)!=NULL) {
    if(line[0] == '#') {
      sscanf(line, "# %d %d %d %d %lf %lf %lf %lf", 
	     &len, &nh, &ne, &nc, &acn,&scn,&aco,&sco);
      ans->nh = nh;
      ans->ne = ne;
      ans->nc = nc;
      ans->ave_cn = acn;
      ans->sd_cn = scn;
      ans->ave_co = aco;
      ans->sd_co = sco;
      continue;
    }
    sscanf(line, " %d %c %c %lf %lf", &ires, &aa, &sec,&cn,&rwco);
    ans->ss[i] = sec;
    if(sec == 'H') {
      ans->svec[i][0] =  1.0;
    }
    else if(sec == 'E') {
      ans->svec[i][1] =  1.0;
    }
    else {
      ans->svec[i][2] =  1.0;
    }
    ans->nvec[i] = cn;
    ans->ovec[i] = rwco;
    
    i++;
    if(i>= MAXSEQ) {
      eprintf("read_ans: maximum sequence length (%d) exceeded in %s\n",
	      MAXSEQ, fname);
    }
  }

  if(i != len) {
    eprintf("read_ans: sequence length inconsistent: %d %d\n", len, i);
  }
  
  ans->len = i;

  fclose(ic);

  return;
}

void pssm2iunit(const int window, const int ind, const pssm_t *apssm, 
                double iunit[])
{
  int is, ie, i, j, n, nunits;
  static const double term[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
				0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
				1.0};

  is = ind - window;
  ie = ind + window;
  nunits = (2*window + 1)*21 + 1;

  for(n=0; n <nunits; n++) iunit[n] = 0.0;

  n=0;
  for(i = is; i <= ie; i++) {
    if(i<0 ||i >= apssm->len) {
      for(j = 0; j < 21; j++) {
	iunit[n++] = term[j];
      }
    }
    else {
      for(j=0; j < 20; j++) {
	iunit[n++] = apssm->prof[i][j];
      }
      iunit[n++] = (i==ind) ? 1.0 : 0.0;
    }
  }

#if (NI_LEVEL==1 || NI_LEVEL==3)
  /* amino acid composition */
  for(j=0; j<20; j++) {
    iunit[n++] = apssm->compo[j];
  }
#endif
#if (NI_LEVEL==2 || NI_LEVEL==3)
  for(j=0; j<20; j++) {
    iunit[n++] = apssm->wlcompo[ind][j];
  }
#endif

  return;
}

int read_file_list(char list[], mydata_t *adata)
{
  FILE *ic;
  char line[MAXLINE];
  char f1[MAXLINE],f2[MAXLINE], f3[MAXLINE];
  int n=0;
  long ndata=0;
 
  if((ic = fopen(list,"r")) == NULL) {
    eprintf("read_file_list: cannot open file: %s\n", list);
  }

  while(fgets(line,MAXLINE,ic) != NULL) {
    if(line[0] == '#') continue;
    sscanf(line,"%s %s %s", f1,f2,f3);

    read_ans(f1, adata->ans + n);
    read_pssm(f2, adata->pssm + n);
    adata->pssm[n].fname = estrdup(f2);
    adata->xxfile[n] = estrdup(f3);

    ndata += (long) adata->pssm[n].len;
    n++;
    if(n>= MAXFILES) {
      eprintf("read_file_list: maximum number of files exceeded: MAXFILES=%d\n",
              MAXFILES);
    }
  }

  fclose(ic);
  adata->nent = n;
  adata->ndat = ndata;

  return n;
}

