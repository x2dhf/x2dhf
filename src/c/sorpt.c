
// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2023  Jacek Kobus 

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "sorpt.h"

void sorpt_ ( void *args )
{
  struct sor_t  *params = (struct sor_t *)args;

  int isym          = params->isym;
  int *indx         = params->indx;
  int *indxex       = params->indxex;  
  int *indx6a       = params->indx6a;
  int *indx6b       = params->indx6b;
  int *indx7        = params->indx7; 

  //double *b         = params->b;
  double *b         = &c_supplptr_.cw_suppl[c_i4b_.i4barr[1]-1];
  //double *d         = params->d;
  double *d         = &c_supplptr_.cw_suppl[c_i4b_.i4barr[2]-1];
  double *e         = params->e;
  //double *e         = &c_supplptr_.cw_suppl[c_i4b_.i4barr[3]-1];
  double *g         = params->g;
  //double *g         = &c_supplptr_.cw_suppl[c_i4b_.i4barr[9]-1];
  double *uext      = params->uext;  

  int ngrd6a=c_interface_17_.ngrd6ac;
  int ngrd6b=c_interface_17_.ngrd6bc;
  int ngrd7=c_interface_17_.ngrd7c;

  int nnu1=c_interface_17_.nnu1c;
  int nnu2=c_interface_17_.nnu2c;
  int nnu3=c_interface_17_.nnu3c;
  int nnu4=c_interface_17_.nnu4c;
  int nnu5=c_interface_17_.nnu5c;
  int isstart=c_interface_17_.isstartc;
  int isstop=c_interface_17_.isstopc;
  int maxsor2=c_interface_17_.maxsor2c;

  double *dmu1=c_interface_40_.dmu1c;
  double *dmu2=c_interface_41_.dmu2c;
  double *dnu1=c_interface_42_.dnu1c;
  double *dnu2=c_interface_43_.dnu2c;
  double *exeven=c_interface_44_.exevenc;
  double omega=c_interface_45_.omegac;
  double omega1=c_interface_45_.omega1c;  
  
  double ddnu1,ddnu2,ddmu1,ddmu2,diffop;
  int i,miter;
  
  for ( miter=1; miter<=maxsor2; miter++) {
    for ( i = isstart; i <=isstop; i++ ) {
      int kext = indxex[i-1] - 1;
      int k  = indx[i-1] - 1;
      
      ddmu2 = dmu2[0] * ( uext[kext-nnu4] + uext[kext+nnu4] ) + 
	dmu2[1] * ( uext[kext-nnu3] + uext[kext+nnu3] ) +   
	dmu2[2] * ( uext[kext-nnu2] + uext[kext+nnu2] ) + 
	dmu2[3] * ( uext[kext-nnu1] + uext[kext+nnu1] );
      
      ddmu1 = dmu1[0] * ( uext[kext-nnu4] - uext[kext+nnu4] ) + 
	dmu1[1] * ( uext[kext-nnu3] - uext[kext+nnu3] ) +   
	dmu1[2] * ( uext[kext-nnu2] - uext[kext+nnu2] ) +   
	dmu1[3] * ( uext[kext-nnu1] - uext[kext+nnu1] );
      
      ddnu2 = dnu2[0] * ( uext[kext-  4] + uext[kext+  4] ) + 
	dnu2[1] * ( uext[kext-  3] + uext[kext+  3] ) + 
	dnu2[2] * ( uext[kext-  2] + uext[kext+  2] ) + 
	dnu2[3] * ( uext[kext-  1] + uext[kext+  1] );
      
      ddnu1 = dnu1[0] * ( uext[kext-  4] - uext[kext+  4] ) + 
	dnu1[1] * ( uext[kext-  3] - uext[kext+  3] ) +   
	dnu1[2] * ( uext[kext-  2] - uext[kext+  2] ) +   
	dnu1[3] * ( uext[kext-  1] - uext[kext+  1] );

      diffop     = ddmu2 + b[k]*ddmu1 + ddnu2 + d[k]*ddnu1;
      uext[kext] = omega * (g[k]-diffop)/e[k] + omega1 * uext[kext];
    }

    if (isym == 1) {
      for ( i=1; i<=ngrd7; i++) {
	int ij=indx7[i-1]-1;
	uext[ij]=exeven[0]*uext[ij+nnu1]+exeven[1]*uext[ij+nnu2]+ 
	  exeven[2]*uext[ij+nnu3]+exeven[3]*uext[ij+nnu4]+exeven[4]*uext[ij+nnu5];
      }
      
      for ( i=1; i<=ngrd6a; i++) {      
	int ij=indx6a[i-1]-1;
	uext[ij]=exeven[0]*uext[ij+1]+exeven[1]*uext[ij+2]+
	  exeven[2]*uext[ij+3]+exeven[3]*uext[ij+4]+exeven[4]*uext[ij+5];
      }

      for ( i=1; i<=ngrd6b; i++) {            
	int ij=indx6b[i-1]-1;
	uext[ij]=exeven[0]*uext[ij-1]+exeven[1]*uext[ij-2]+
	  exeven[2]*uext[ij-3]+exeven[3]*uext[ij-4]+exeven[4]*uext[ij-5];
      } 
    } else {
      for ( i=1; i<=ngrd7; i++) {
	int ij=indx7[i-1]-1;
	uext[ij]=0.0;
      }

      for ( i=1; i<=ngrd6a; i++) {
	int ij=indx6a[i-1]-1;
	uext[ij]=0.0;
      }

      for ( i=1; i<=ngrd6b; i++) {
	int ij=indx6b[i-1]-1;
	uext[ij]=0.0;
      }
    }
  }
}

