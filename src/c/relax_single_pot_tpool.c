
// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2023  Jacek Kobus 

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <stdbool.h>
#include "sorpt.h"

extern struct exchsor_t params[max_threads4pots];

extern pthread_barrier_t barrier4pots;
extern int worker1_threads_done;
extern int nthreads;
extern pthread_mutex_t mutex;
extern void sorpt_ ( void *args );
extern void putinc_ (int * nnu, int * mxnmu, int * isymetry, double fun[], double work[]);
extern void putoutc_ (int * nnu, int * mxnmu, double fun[], double work[]);
extern bool quit;

extern pthread_cond_t exchsor_start[max_threads4pots];
extern pthread_cond_t exchsor_stop[max_threads4pots];
extern pthread_cond_t exchsor_ready;

extern int exchsor_threads_done;
extern int exchsor_threads_waiting;

void relax_single_pot_tpool_ ( void * args )
{
  // struct exchsor_t *paramsIN = (struct exchsor_t *)args;
  //struct exchsor_t *params = (struct exchsor_t *) arg;  

  int threadNum= *(int *) args;
  
  //threadNum=params[threadNum].threadNum;
  //int threadsNum=params[threadNum].threadsNum;
  int thread=threadNum;
  //int iorb=c_interface_46_.iorbc;

  //pthread_t threads[threadsNum];
  struct sor_t paramsOUT;

  int nnu=c_interface_17_.nnu1c-8;
  int mxnmu=c_interface_17_.mxnmuc;
  int mxsize=c_interface_17_.mxsizec;
  int maxsor1=c_interface_17_.maxsor1c;
  
  //printf("nnu mxnmu mxsize maxsor1 %5d %5d %5d %5d \n",nnu,mxnmu,mxsize, maxsor1);

  //printf ("threadNum,nnu mxnmu isym %5d %5d %5d %5d %5d  \n", threadNum,nnu, mxnmu, isym, mxsize);
  //printf ("threadNum, ib1,ib2, %5d %5d %5d  \n", threadNum, ib1, ib2);
  //printf ("nnu mxnmu isym  \n");

  double  *lhs = (double*) malloc(mxsize * sizeof(double));
  double  *rhs = (double*) malloc(mxsize * sizeof(double));
  double  *wk2 = (double*) malloc((nnu+8)*(mxnmu+8) * sizeof(double));
  int i;

  //printf("threadsNum threadNum mxsize= %d %d %d \n",threadsNum,threadNum,mxsize);
  //printf("threadNum,nnu,mxnmu ibexc isym %8d %8d %8d %8d %8d \n",threadNum,nnu,mxnmu,ibexc,isym);

  //printf("idel2 isym %12.4e %5d %5d \n",idel2,ibexc,isym);
  //exit(2);
  // prepare right-hand side of the Poisson equation
  int nexchpot;

  while (1) {

    int iorbc=c_interface_46_.iorbc;
    nexchpot=c_interface_8_.nexchpotsc[c_interface_46_.iorbc-1];
    int threadsNum=c_interface_46_.nthreadsc;
    //int threadsNum=nexchpot;
    
#ifdef TRACE1
    printf("TRACE1: relax_single_pot_tpool/iorbc threadNum: %2d %2d \n",iorbc,threadNum);    
    printf("TRACE1: relax_single_pot_tpool/iorbc threadsNum nthreads: %2d %2d %2d\n",iorbc,threadsNum, nthreads);    
    printf("TRACE1: relax_single_pot_tpool/iorbc threadsNum nthreads: %2d %2d %2d\n",
	   iorbc,c_interface_46_.nthreadsc, nthreads);
#endif	   
    pthread_mutex_lock(&mutex);
    exchsor_threads_waiting++;
    if (exchsor_threads_waiting==c_interface_46_.nthreadsc) { 
      pthread_cond_signal(&exchsor_ready); 
    }
#ifdef TRACE2
    printf("TRACE2: relax_single_pot_tpool/thread=%1d \t wait exchsor_start\n",thread);
#endif
    pthread_cond_wait(&exchsor_start[thread],&mutex); 
#ifdef TRACE3    
    printf("TRACE3: relax_single_pot_tpool/thread=%1d %s\n", thread, quit ? "true" : "false");
#endif
    if (quit == true) {
      pthread_mutex_unlock(&mutex);
      free(lhs);
      free(rhs);
      free(wk2);
      break;
    }
    pthread_mutex_unlock(&mutex);
 
    threadsNum=c_interface_8_.nexchpotsc[c_interface_46_.iorbc-1];
    
#ifdef TRACE2
    printf("TRACE2: relax_single_pot_tpool/iorb threadNum threadsNum nthreads: %2d %2d %2d %2d \n",
	   iorbc, thread, c_interface_46_.iorbc,threadNum,threadsNum, nthreads);

    printf("TRACE2: relax_single_pot_tpool/iorb thread threadsNum: %2d %2d %2d\n",
	   c_interface_46_.iorbc,thread,c_interface_8_.nexchpotsc[c_interface_46_.iorbc-1]);
#endif

    paramsOUT.indx=&c_sorptr_.cw_sor[c_iadex_.iadnorc-1];        
    paramsOUT.indx=params[thread].indx;
    //paramsOUT.indxex=params[thread].indxex;
    paramsOUT.indxex=&c_sorptr_.cw_sor[c_iadex_.iadextc-1];
    //paramsOUT.indx6a=params[thread].indx6a;
    paramsOUT.indx6a=&c_sorptr_.cw_sor[c_iadex_.iadex1c-1];
    //paramsOUT.indx6b=params[thread].indx6b;
    paramsOUT.indx6b=&c_sorptr_.cw_sor[c_iadex_.iadex2c-1];
    //paramsOUT.indx7=params[thread].indx7;
    paramsOUT.indx7=&c_sorptr_.cw_sor[c_iadex_.iadex3c-1];
    //paramsOUT.b=params[thread].b;
    //paramsOUT.d=params[thread].d;
    paramsOUT.b=&c_supplptr_.cw_suppl[c_i4b_.i4barr[1]-1];
    //paramsOUT.d=&c_supplptr_.cw_suppl[c_i4b_.i4barr[2]-1];

    //int isym  = c_interface_6_.isymsc[thread][iorbc-1];
    int isym  = c_interface_6_.isymsc[thread][c_interface_46_.iorbc-1];
    paramsOUT.isym=c_interface_6_.isymsc[thread][c_interface_46_.iorbc-1];

    int in1 = c_interface_1_.ins1[thread][c_interface_46_.iorbc-1];
    int in2 = c_interface_2_.ins2[thread][c_interface_46_.iorbc-1];
    int *i1b=c_interface_7_.i1bc;
    
    int ib1=i1b[in1-1];
    int ib2=i1b[in2-1];
    
    int ibexc = c_interface_3_.ibexcpc[thread][c_interface_46_.iorbc-1];
    int idel  = c_interface_4_.deltam4potc[thread][c_interface_46_.iorbc-1];

    double idel2=(double) idel*idel;
    
    for (i=0; i<mxsize; i++) {
      //    //pthread_attr_init(&attr);                                                                                               
      //printf ("exchsorsingle1: thread i %1d %5d\n",thread,i);
      //rhs[i]=params[thread].psi[ib1+i-1]*params[thread].psi[ib2+i-1]*params[thread].gf[i];
      rhs[i]=c_orbptr_.cw_orb[ib1+i-1]*c_orbptr_.cw_orb[ib2+i-1]*params[thread].gf[i];
    }

    //printf("exchsorptdrv1: orbptr: %f %2f \n",c_orbptr_.orbptr[0],c_orbptr_.orbptr[9]);
    //printf("exchsorptdrv1: orbptr: %f  \n",c_orbptr_.cw_orb[ib1+10]);
    paramsOUT.g=rhs;
    
    // prepare left-hand side of the Poisson equation include the diagonal part of
    // the differentiation operator in lhs
    if (idel2==0.0) {
      for (i=0; i<mxsize; i++) {
	lhs[i]=params[thread].f3[i];
      }
    } else {
      for (i=0; i<mxsize; i++) {
	lhs[i]=params[thread].f3[i] + params[thread].ef[i]*idel2;
      }
    }
    paramsOUT.e=lhs;
    
    for (i=1; i<=maxsor1; i++) {
#ifdef TRACE2
      printf ("TRACE2: relax_single_pot_tpool/iorbc, thread nnu mxnmu isym ibexc %5d %5d %5d %5d %5d %5d \n",
	      c_interface_46_.iorbc, thread, nnu, mxnmu, isym, c_interface_3_.ibexcpc[thread][c_interface_46_.iorbc-1]);
#endif
	      putinc_ ( &nnu, &mxnmu, &isym, &params[thread].excp[c_interface_3_.ibexcpc[thread][c_interface_46_.iorbc-1]-1], wk2);
      paramsOUT.uext=wk2;
      sorpt_ ( (void *) &paramsOUT );
      putoutc_ ( &nnu, &mxnmu, &params[thread].excp[c_interface_3_.ibexcpc[thread][c_interface_46_.iorbc-1]-1], wk2);
    }
#ifdef TRACE
    printf ("TRACE:  relax_single_pot_tpool/sorpt done for thread: %5d \n",threadNum);
#endif
    
    pthread_mutex_lock(&mutex);    
    exchsor_threads_done++;


#ifdef TRACE    
    printf("TRACE:  relax_single_pot_tpool/exchsor_threads_done %1d %1d \n",thread,exchsor_threads_done);
#endif

    if (exchsor_threads_done==threadsNum) {
#ifdef TRACE
      printf("TRACE:  relax_single_pot_tpool/thread=%1d exchsor_threads_done= %1d \n",thread,exchsor_threads_done);	
#endif

      exchsor_threads_done=0;
      exchsor_threads_waiting=0;
      for (i=0; i<c_interface_46_.nthreadsc; i++) {
	if (i != thread) {
#ifdef TRACE2
	  printf("TRACE2: relax_single_pot_tpool/thread signal work_stop to thread %1d %1d \n",thread,i);
#endif
	  pthread_cond_signal(&exchsor_stop[i]);
	}
      }
      pthread_mutex_unlock(&mutex);    
      continue;
    }
#ifdef TRACE2
    printf("TRACE2: relax_single_pot_tpool/wait on work-stop %1d \n",threadNum);
#endif
    pthread_cond_wait(&exchsor_stop[thread],&mutex);
    pthread_mutex_unlock(&mutex);    
  }
  
  return;
}  

