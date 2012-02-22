/*-----------------------------------------------------------
/
/   Program:      sov.c
/
/   Secondary structure prediction accuracy evaluation
/
/   SOV (Segment OVerlap) measure 
/
/   Copyright by Adam Zemla (11/16/1996)
/   Email: adamz@llnl.gov  
/
/------------------------------------------------------------   
/
/   Compile:      cc sov.c -o sov -lm
/
/------------------------------------------------------------*/
/*
 * modified by Akira Kinjo (2005).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sov.h"

/* 	$Id: sov.c,v 1.2 2006/02/20 08:37:17 akinjo Exp $	 */

#ifndef lint
static char vcid[] = "$Id: sov.c,v 1.2 2006/02/20 08:37:17 akinjo Exp $";
#endif /* lint */

/*-------------------------------------------------------------
/
/  default_parameters - default parameters for SOV program
/
/------------------------------------------------------------*/
void sov_default_parameters(sov_t *pdata)
{
  pdata->input=0;
  pdata->order=0;
  pdata->sov_method=1;
  pdata->sov_delta=1.0;
  pdata->sov_delta_s=0.5;
  pdata->sov_out=0;

  return;
}

/*-----------------------------------------------------------
/
/    sov - evaluate SSp by the Segment OVerlap quantity (SOV)
/                Input: secondary structure segments
/
/------------------------------------------------------------*/
double evalsov(int n_aa, char sss1[MAXRES], char sss2[MAXRES], sov_t *pdata)
{
  int i, k, length1, length2, beg_s1, end_s1, beg_s2, end_s2;
  int j1, j2, k1, k2, minov, maxov, d, d1, d2, n, multiple;
  char s1, s2, sse[3];
  double out; 
  double s, x;

  sse[0]='#';
  sse[1]='#';
  sse[2]='#';

  if(pdata->sov_what==0) {
    sse[0]='H';
    sse[1]='E';
    sse[2]='C';
  }
  if(pdata->sov_what==1) {
    sse[0]='H';
    sse[1]='H';
    sse[2]='H';
  }
  if(pdata->sov_what==2) {
    sse[0]='E';
    sse[1]='E';
    sse[2]='E';
  }
  if(pdata->sov_what==3) {
    sse[0]='C';
    sse[1]='C';
    sse[2]='C';
  }
  n=0;
  for(i=0;i<n_aa;i++) {
    s1=sss1[i];
    if(s1==sse[0] || s1==sse[1] || s1==sse[2]) {
      n++;
    }
  }
  out=0.0;
  s=0.0;
  length1=0;
  length2=0;
  i=0;
  while(i<n_aa) {
    beg_s1=i;
    s1=sss1[i];
    while(sss1[i]==s1 && i<n_aa) {
      i++;
    }
    end_s1=i-1;
    length1=end_s1-beg_s1+1;
    multiple=0;
    k=0;
    while(k<n_aa) {
      beg_s2=k;
      s2=sss2[k];
      while(sss2[k]==s2 && k<n_aa) {
        k++;
      }
      end_s2=k-1;
      length2=end_s2-beg_s2+1;
      if(s1==sse[0] || s1==sse[1] || s1==sse[2]) {
        if(s1==s2 && end_s2>=beg_s1 && beg_s2<=end_s1) {
          if(multiple>0 && pdata->sov_method==1) {
            n=n+length1;
          }
          multiple++;
          if(beg_s1>beg_s2) {
            j1=beg_s1;
            j2=beg_s2;
          }
          else {
            j1=beg_s2;
            j2=beg_s1;
          }
          if(end_s1<end_s2) {
            k1=end_s1;
            k2=end_s2;
          }
          else {
            k1=end_s2;
            k2=end_s1;
          }
          minov=k1-j1+1;
          maxov=k2-j2+1;
          d1=floor(length1*pdata->sov_delta_s);
          d2=floor(length2*pdata->sov_delta_s);
          if(d1>d2) d=d2;
          if(d1<=d2 || pdata->sov_method==0) d=d1;
          if(d>minov) {
            d=minov;
          }
          if(d>maxov-minov) {
            d=maxov-minov;
          }
          x=pdata->sov_delta*d;
          x=(minov+x)*length1;
          if(maxov>0) {
            s=s+x/maxov;
          }
          else {
            printf("\n ERROR! minov = %-4d maxov = %-4d length = %-4d d = %-4d   %4d %4d  %4d %4d",
                     minov,maxov,length1,d,beg_s1+1,end_s1+1,beg_s2+1,end_s2+1);
          }
          if(pdata->sov_out==2) {
            printf("\n TEST: minov = %-4d maxov = %-4d length = %-4d d = %-4d   %4d %4d  %4d %4d",
                   minov,maxov,length1,d,beg_s1+1,end_s1+1,beg_s2+1,end_s2+1);
          }
        }
      }
    }
  }
  if(pdata->sov_out==2) {
    printf("\n TEST: Number of considered residues = %d",n);
  }
  if(n>0) {
    out=s/n;
  }
  else {
    out=1.0;
  }
  return out*100.0;
}

/*-----------------------------------------------------------
/
/    Q3 - evaluate SSp by the residues predicted correctly (Q3)
/                Input: secondary structure segments
/
/------------------------------------------------------------*/
double evalq3(int n_aa, char sss1[MAXRES], char sss2[MAXRES], sov_t *pdata)
{
  int i, n;
  double out;
  char s, sse[3];

  sse[0]='#';
  sse[1]='#';
  sse[2]='#';

  if(pdata->q3_what==0) {
    sse[0]='H';
    sse[1]='E';
    sse[2]='C';
  }
  if(pdata->q3_what==1) {
    sse[0]='H';
    sse[1]='H';
    sse[2]='H';
  }
  if(pdata->q3_what==2) {
    sse[0]='E';
    sse[1]='E';
    sse[2]='E';
  }
  if(pdata->q3_what==3) {
    sse[0]='C';
    sse[1]='C';
    sse[2]='C';
  }

  n=0;
  out=0.0;
  for(i=0;i<n_aa;i++) {
    s=sss1[i];
    if(s==sse[0] || s==sse[1] || s==sse[2]) {
      n++;
      if(sss1[i]==sss2[i]) {
        out=out + 1.0;
      }
    }
  }
  return out;
}
      
/* Matthews' correlation */
void evalmc(int n_aa, char typ, char sss1[MAXRES], char sss2[MAXRES], 
    sov_t *pdata, int *ntt, int *nff, int *ntf,int *nft)
{
  int i;
  char s1,s2;

  *ntt=*ntf=*nft=*nff = 0;
  for(i=0;i<n_aa;i++) {
    s1=sss1[i];
    s2=sss2[i];
    if(s1 == typ && s2 == typ) {
      (*ntt)++;
    }
    else if(s1 == typ && s2 != typ) {
      (*ntf)++;
    }
    else if(s1 != typ && s2 == typ) {
      (*nft)++;
    }
    else if(s1 != typ && s2 != typ) {
      (*nff)++;
    }
  }

  return;

}

void check_accuracy_ss(ans_t *aans, char psec[], double *q3, double *sov3)
{
  sov_t pdata;
  int ntt,nff,ntf,nft;
  int nh,ne,nc,i,ph,pe,pc;
  double qs[3];
  
  sov_default_parameters(&pdata);

  pdata.q3_what = 0;
  *q3 = 100.0*evalq3(aans->len,aans->ss,psec,&pdata)/(double)aans->len;
  pdata.q3_what = 1;
  qs[0] = evalq3(aans->len,aans->ss,psec,&pdata);
  pdata.q3_what = 2;
  qs[1] = evalq3(aans->len,aans->ss,psec,&pdata);
  pdata.q3_what = 3;
  qs[2] = evalq3(aans->len,aans->ss,psec,&pdata);
  
  nh = ne = nc = 0;
  ph = pe = pc = 0;
  for(i = 0; i < aans->len; i++) {
    if (aans->ss[i] == 'H') nh++;
    else if (aans->ss[i] == 'E') ne++;
    else if (aans->ss[i] == 'C') nc++;

    if (psec[i] == 'H') ph++;
    else if (psec[i] == 'E') pe++;
    else if (psec[i] == 'C') pc++;

  }
  fprintf(stderr,"Q1: %5.0f %5d %5d %5.0f %5d %5d %5.0f %5d %5d (HEC)\n", 
	  qs[0], nh, ph, qs[1], ne,pe, qs[2],nc,pc);

  pdata.sov_what = 0;
  *sov3 = evalsov(aans->len,aans->ss,psec,&pdata);
  pdata.sov_what = 1;
  qs[0] = evalsov(aans->len,aans->ss,psec,&pdata);
  pdata.sov_what = 2;
  qs[1] = evalsov(aans->len,aans->ss,psec,&pdata);
  pdata.sov_what = 3;
  qs[2] = evalsov(aans->len,aans->ss,psec,&pdata);
  fprintf(stderr,"SOV1: %5.0f %5.0f %5.0f (HEC)\n", qs[0], qs[1],qs[2]);

  evalmc(aans->len,'H', aans->ss, psec, &pdata,&ntt,&nff,&ntf,&nft);
  fprintf(stderr,"MC_H: %5d %5d %5d %5d (11,00,10,01)\n", ntt, nff, ntf, nft);

  evalmc(aans->len,'E', aans->ss, psec, &pdata,&ntt,&nff,&ntf,&nft);
  fprintf(stderr,"MC_E: %5d %5d %5d %5d (11,00,10,01)\n", ntt, nff, ntf, nft);

  evalmc(aans->len,'C', aans->ss, psec, &pdata,&ntt,&nff,&ntf,&nft);
  fprintf(stderr,"MC_C: %5d %5d %5d %5d (11,00,10,01)\n", ntt, nff, ntf, nft);

  return;
}
