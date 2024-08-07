// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2023  Jacek Kobus 


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "sorpt.h"

extern void sorpt_ ( void *args );
extern void putinc_ (int * nnu, int * mxnmu, int * isymetry, double fun[], double work[]);
extern void putoutc_ (int * nnu, int * mxnmu, double fun[], double work[]);

void relax_single_pot_pthread_ ( void *args )  
{
  struct exchsor_t *paramsIN = (struct exchsor_t *)args;
  
  int threadsNum = paramsIN->threadsNum;
  int threadNum  = paramsIN->threadNum;

  double deltam2   = paramsIN->deltam2;     
  double *ef     = paramsIN->ef;      
  double *excp   = paramsIN->excp;        
  double *f3     = paramsIN->f3;
  double *gf     = paramsIN->gf;
  double *psi    = paramsIN->psi;  
  
  int ib1        = paramsIN->ib1;
  int ib2        = paramsIN->ib2;  
  int ibexc      = paramsIN->ibexc;  
  int isym       = paramsIN->isym; 
  //int ipc            = paramsIN->ipc; 

#ifdef TRACE  
  //printf("TRACE: relax_single_pot_pthread/Thread number %ld\n", pthread_self());
printf("TRACE:  relax_single_pot_pthread/Thread number %d\n", threadNum);
#endif
  
  struct sor_t paramsOUT[threadsNum];

  int nnu=c_interface_17_.nnu1c-8;
  int mxnmu=c_interface_17_.mxnmuc;
  int mxsize=c_interface_17_.mxsizec;
  int maxsor1=c_interface_17_.maxsor1c;

  paramsOUT[threadsNum].threadsNum=threadsNum;  
  paramsOUT[threadNum].threadNum=threadNum;

  //printf("nnu mxnmu mxsize maxsor1 %5d %5d %5d %5d \n",nnu,mxnmu,mxsize, maxsor1);

  paramsOUT[threadNum].indx=paramsIN->indx;
  paramsOUT[threadNum].indxex=paramsIN->indxex;
  paramsOUT[threadNum].indx6a=paramsIN->indx6a;
  paramsOUT[threadNum].indx6b=paramsIN->indx6b;
  paramsOUT[threadNum].indx7=paramsIN->indx7;
  paramsOUT[threadNum].b=paramsIN->b;
  paramsOUT[threadNum].d=paramsIN->d;
  paramsOUT[threadNum].isym=isym;

  //printf ("threadNum,nnu mxnmu isym %5d %5d %5d %5d %5d  \n", threadNum,nnu, mxnmu, isym, mxsize);
  //printf ("threadNum, ib1,ib2, %5d %5d %5d  \n", threadNum, ib1, ib2);
  //printf ("nnu mxnmu isym  \n");

  double  *lhs = (double*) malloc(mxsize * sizeof(double));
  double  *rhs = (double*) malloc(mxsize * sizeof(double));
  double  *wk2 = (double*) malloc((nnu+8)*(mxnmu+8) * sizeof(double));
  int i;

  //printf("threadsNum threadNum mxsize= %d %d %d \n",threadsNum,threadNum,mxsize);
  //printf("threadNum,nnu,mxnmu ibexc isym %8d %8d %8d %8d %8d \n",threadNum,nnu,mxnmu,ibexc,isym);

  //printf("deltam2 isym %12.4e %5d %5d \n",deltam2,ibexc,isym);
  // prepare right-hand side of the Poisson equation
  for (i=0; i<mxsize; i++) {
    rhs[i]=psi[ib1+i-1]*psi[ib2+i-1]*gf[i];
  }
  paramsOUT[threadNum].g=rhs;

  // prepare left-hand side of the Poisson equation include the diagonal part of
  // the differentiation operator in lhs
  if (deltam2==0.0) {
    for (i=0; i<mxsize; i++) {
      lhs[i]=f3[i];
    }
  } else {
    for (i=0; i<mxsize; i++) {
      lhs[i]=f3[i] + ef[i]*deltam2;
    }
  }
  paramsOUT[threadNum].e=lhs;
  
  for (i=1; i<=maxsor1; i++) {
#ifdef TRACE2
    printf ("TRACE2: relax_single_pot_pthread/i threadNum isym ibexc %5d %5d %5d %5d \n",i, threadNum, isym, ibexc);
#endif

    putinc_ ( &nnu, &mxnmu, &isym, &excp[ibexc-1], wk2);
    paramsOUT[threadNum].uext=wk2;

    sorpt_ ( (void *) &paramsOUT[threadNum] );
    
    putoutc_ ( &nnu, &mxnmu, &excp[ibexc-1], wk2);
  }

  free(lhs);
  free(rhs);
  free(wk2);
}  

