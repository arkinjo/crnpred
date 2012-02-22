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
/   modified by A.R. Kinjo (21/02/2005)
/
/------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "blast.h"

#define MAXRES   MAXSEQ

typedef struct {
  int     input;
  int     order;
  int     q3_what;
  int     sov_what;
  int     sov_method;
  double   sov_delta;
  double   sov_delta_s;
  int     sov_out;
} sov_t;

