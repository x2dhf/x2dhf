
// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2023  Jacek Kobus 


#include <stdlib.h>
#include <stdio.h>
#include <sched.h>
#include <unistd.h>
#include <time.h>
#include <stdalign.h> 
#include <stdbool.h>
#include <pthread.h>
#include "sorpt.h"

extern pthread_barrier_t barrier4mcsor;
extern pthread_mutex_t mutex4mcsor;
extern pthread_cond_t work_stop;
extern pthread_cond_t work_start;
extern pthread_cond_t main_ready;
extern bool quit4mcsor;
extern int worker_threads_done;
extern int worker_threads;
//extern struct sor_t params4mcsor;
extern int nthreads;

#ifdef MUTEX
extern pthread_mutex_t mutex4mcsor;
#endif

extern int worker_threads_waiting;

void mcsor_single_colour_tpool_ (void *args) {
    
  struct sor_t *params4mcsor = (struct sor_t *) args;
  
  int threadsNum    = params4mcsor->threadsNum;  
  int threadNum     = params4mcsor->threadNum;
  int thread        = threadNum;
  int isym          = params4mcsor->isym;
  int *indxex       = &c_sorptr_.cw_sor[c_iadex_.iadextc-1];
  int *indx         = &c_sorptr_.cw_sor[c_iadex_.iadnorc-1];    
  int *indx6a       = &c_sorptr_.cw_sor[c_iadex_.iadex1c-1];
  int *indx6b       = &c_sorptr_.cw_sor[c_iadex_.iadex2c-1];
  int *indx7        = &c_sorptr_.cw_sor[c_iadex_.iadex3c-1];  

  double *b         = &c_supplptr_.cw_suppl[c_i4b_.i4barr[0]-1];  
  double *d         = &c_supplptr_.cw_suppl[c_i4b_.i4barr[2]-1];         

  double *e         = params4mcsor->e;
  double *g         = params4mcsor->g;
  double *uext      = params4mcsor->uext;          

  int *nstart=c_interface_9_.nstartc;
  int *nstop=c_interface_10_.nstopc;  
  int *nstart6a=c_interface_11_.nstart6ac;
  int *nstop6a=c_interface_12_.nstop6ac;  
  int *nstart6b=c_interface_13_.nstart6bc;
  int *nstop6b=c_interface_14_.nstop6bc;  
  int *nstart7=c_interface_15_.nstart7c;
  int *nstop7=c_interface_16_.nstop7c;  

  int nnu1=c_interface_17_.nnu1c;
  int nnu2=c_interface_17_.nnu2c;
  int nnu3=c_interface_17_.nnu3c;
  int nnu4=c_interface_17_.nnu4c;
  int nnu5=c_interface_17_.nnu5c;
  int maxsor2=c_interface_17_.maxsor2c;

  double *dmu1=c_interface_40_.dmu1c;
  double *dmu2=c_interface_41_.dmu2c;
  double *dnu1=c_interface_42_.dnu1c;
  double *dnu2=c_interface_43_.dnu2c;
  double *exeven=c_interface_44_.exevenc;
  double omega=c_interface_45_.omegac;
  double omega1=c_interface_45_.omega1c;  
  
  int i,ij,k,kext,miter;
  int colour,start,stop;
  double ddnu1,ddnu2,ddmu1,ddmu2,diffop;

/* #ifdef TRACE3 */
/*   printf("TRACE3: mcsorpt: 1 threadNum threadsNum %5d %5d \n",thread,threadsNum);   */
/*   printf("TRACE3: mcsorpt: 1 threadNum thread %5d %5d \n",threadNum,thread); */
/*   printf("TRACE3: mcsorpt: 1 threadNum indx     %5d %15p \n",threadNum,indx);     */
/*   printf("TRACE3: mcsorpt: 1 threadNum indxex   %5d %15p \n",threadNum,indxex); */
/*   printf("TRACE3: mcsorpt: 1 threadNum dmu1[01]               \n"); */
/*   printf("TRACE3: mcsorpt: 1 threadNum dmu1[01] %5d %15e %15e \n",threadNum,dmu1[0],dmu1[1]); */
/* #endif */
  //int thread=threadNum;

  while (1) {
    if (threadsNum>1) {
      //printf("mcsorpt: thread=%1d worker_threads_done= %1d \n",thread,worker_threads_done);	
      pthread_mutex_lock(&mutex4mcsor);
      worker_threads_waiting++;
      if (worker_threads_waiting==threadsNum) { 
	pthread_cond_signal(&main_ready); 
      }
#ifdef TRACE2
      printf("TRACE2: mcsorpt/thread=%1d \t wait work_start\n",thread);
#endif
      pthread_cond_wait(&work_start,&mutex4mcsor); 
#ifdef TRACE2    
      printf("TRACE2: mcsorpt/thread=%1d %s\n", thread, quit4mcsor ? "true" : "false");
#endif
      if (quit4mcsor == true) {
	pthread_mutex_unlock(&mutex4mcsor);
	break;
      }
      pthread_mutex_unlock(&mutex4mcsor);
    }

#ifdef TRACE2    
    printf("TRACE2: mcsorpt/thread=%1d  false & continue 1\n", thread);
#endif


#ifdef TRACE3
    //printf("TRACE3: mcsorpt: 1 threadNum %5d  \n",threadNum);
    //printf("TRACE3: mcsorpt: 1 threadNum indx   %5d %15p \n",threadNum,indx);  
    //printf("TRACE3: mcsorpt: 1 threadNum indx[100] %5d %5d \n",threadNum,indx[100]);
    //printf("TRACE3: mcsorpt: 1 threadNum g[01] %5d %15e %15e \n",threadNum,g[0],g[1]);
#endif
    
#ifdef TRACE3    
    //printf("TRACE3: mcsorpt: 2 threadNum nnu1    %5d %5d \n",threadNum,nnu1);
    //printf("TRACE3: mcsorpt: 2 threadNum nnu5    %5d %5d \n",threadNum,nnu5);
    //printf("TRACE3: mcsorpt: 2 threadNum maxsor2 %5d %5d \n",threadNum,maxsor2);
    //printf("TRACE3: mcsorpt: 2 threadNum omega   %5d %15e
    //\n",threadNum,omega);
    //printf("TRACE3: 1 mcsorpt: threadNum maxsor2 %5d %5d \n",threadNum,maxsor2);
    /* printf("TRACE3: mcsorpt: 2 threadNum indx     %5d %15p \n",threadNum,params4mcsor->indx);     */
    /* printf("TRACE3: mcsorpt: 2 threadNum indxex   %5d %15p \n",threadNum,params4mcsor->indxex); */
    /* printf("TRACE3: mcsorpt: 2 b[10,100] d[10,100] \n"); */
    /* printf("TRACE3: mcsorpt: 2 threadNum b[10,100]%5d %15e %15e \n",threadNum,b[0],b[1]); */
    /* printf("TRACE3: mcsorpt: 2 threadNum d[10,100] %5d %15e %15e\n",threadNum, d[10],d[100]);     */
    /* printf("TRACE3: mcsorpt: 2 e[10,100] indx[10,100] \n"); */
    /* printf("TRACE3: mcsorpt: 2 threadNum e[10,100]%5d %15e %15e \n",threadNum,e[0],e[1]); */
    /* printf("TRACE3: mcsorpt: 2 threadNum indx[10,100] %5d %5d %5d\n",threadNum, indx[10],indx[100]);     */
    
#endif
    for ( miter=1; miter<=maxsor2; miter++) {

#ifdef TRACE2    
      printf("TRACE2: mcsorpt/thread=%1d  miter=%2d \n", thread, miter);
#endif
      
      if ( threadsNum > 1 ) {
	pthread_barrier_wait(&barrier4mcsor);
      }

      
      colour=0;
      start=nstart[threadNum*5+colour];
      stop= nstop[threadNum*5+colour];
      
/* #ifdef TRACE3     */
/*       printf("TRACE3: mcsorpt: threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop); */
/* #endif */
      
      for ( i = start; i <=stop; i++ ) {
	kext = indxex[i-1] - 1;
	k  = indx[i-1] - 1;
	
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
#ifdef MUTEX
	pthread_mutex_lock(&mutex4mcsor);
#endif
	uext[kext] = omega * (g[k]-diffop)/e[k] + omega1 * uext[kext];
#ifdef MUTEX
	pthread_mutex_unlock(&mutex4mcsor);
#endif
	
	//if (k % 1000 == 0) {
	//printf (" ZZZ %3d %8d%16.6e%16.6e%16.6e%16.6e%16.6e\n", threadNum,k+1,b[k],d[k],g[k],e[k],uext[kext]);
	//printf (" ZZZ %8d%16.6e%16.6e%16.6e%16.6e%16.6e\n", k+1,b[k],d[k],g[k],e[k],uext[kext]);
	//printf (" ZZZ %5d %12.4e  %12.4e %12.4e %12.4e %12.4e %12.4e \n", k+1,b[k],d[k],g[k],e[k],diffop,uext[kext]);
	//printf (" XXX %8d%16.6e%16.6e%16.6e%16.6e\n", k+1,ddmu2,ddmu1,ddnu2,ddnu1);
	//}
      }
      if ( threadsNum > 1 ) {

#ifdef TRACE3    
	printf("TRACE3: mcsorpt/before barrier/colour threadNum        0 %5d \n",threadNum);
#endif

	pthread_barrier_wait(&barrier4mcsor);      
#ifdef TRACE3    
	printf("TRACE3: mcsorpt/barrier passed/colour threadNum        0 %5d \n",threadNum);
#endif
      }
      
      
      colour=1;
      start=nstart[threadNum*5+colour];
      stop=nstop[threadNum*5+colour];
#ifdef TRACE3    
      printf("TRACE3: mcsorpt/threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop);
#endif
      
      for ( i = start; i <=stop; i++ ) {
	kext = indxex[i-1] - 1;
	k  = indx[i-1] - 1;
	ddmu2 =
	  dmu2[0] * ( uext[kext-nnu4] + uext[kext+nnu4] ) + 
	  dmu2[1] * ( uext[kext-nnu3] + uext[kext+nnu3] ) +   
	  dmu2[2] * ( uext[kext-nnu2] + uext[kext+nnu2] ) + 
	  dmu2[3] * ( uext[kext-nnu1] + uext[kext+nnu1] );
	
	ddmu1 =
	  dmu1[0] * ( uext[kext-nnu4] - uext[kext+nnu4] ) + 
	  dmu1[1] * ( uext[kext-nnu3] - uext[kext+nnu3] ) +   
	  dmu1[2] * ( uext[kext-nnu2] - uext[kext+nnu2] ) +   
	  dmu1[3] * ( uext[kext-nnu1] - uext[kext+nnu1] );
	
	ddnu2 =
	  dnu2[0] * ( uext[kext-  4] + uext[kext+  4] ) + 
	  dnu2[1] * ( uext[kext-  3] + uext[kext+  3] ) + 
	  dnu2[2] * ( uext[kext-  2] + uext[kext+  2] ) + 
	  dnu2[3] * ( uext[kext-  1] + uext[kext+  1] );
	
	ddnu1 =
	  dnu1[0] * ( uext[kext-  4] - uext[kext+  4] ) + 
	  dnu1[1] * ( uext[kext-  3] - uext[kext+  3] ) +   
	  dnu1[2] * ( uext[kext-  2] - uext[kext+  2] ) +   
	  dnu1[3] * ( uext[kext-  1] - uext[kext+  1] );
	
	diffop     = ddmu2 + b[k]*ddmu1 + ddnu2 + d[k]*ddnu1;
#ifdef MUTEX
	pthread_mutex_lock(&mutex4mcsor);
#endif
	uext[kext] = omega * (g[k]-diffop)/e[k] + omega1 * uext[kext];
#ifdef MUTEX
	pthread_mutex_unlock(&mutex4mcsor);
#endif
      }
      
      if ( threadsNum > 1 ) {

#ifdef TRACE3    
	printf("TRACE3: mcsorpt/before barrier/colour threadNum        1 %5d \n",threadNum);
#endif

	pthread_barrier_wait(&barrier4mcsor);
#ifdef TRACE3    
	printf("TRACE3: mcsorpt/barrier passed/colour threadNum        1 %5d \n",threadNum);
#endif
      }
      
      colour=2;
      start=nstart[threadNum*5+colour];
      stop=nstop[threadNum*5+colour];
      
#ifdef TRACE3    
      printf("TRACE3: mcsorpt/threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop);
#endif
      
      for ( i = start; i <=stop; i++ ) {
	kext = indxex[i-1] - 1;
	k  = indx[i-1] - 1;
	ddmu2 =
	  dmu2[0] * ( uext[kext-nnu4] + uext[kext+nnu4] ) + 
	  dmu2[1] * ( uext[kext-nnu3] + uext[kext+nnu3] ) +   
	  dmu2[2] * ( uext[kext-nnu2] + uext[kext+nnu2] ) + 
	  dmu2[3] * ( uext[kext-nnu1] + uext[kext+nnu1] );
	
	ddmu1 =
	  dmu1[0] * ( uext[kext-nnu4] - uext[kext+nnu4] ) + 
	  dmu1[1] * ( uext[kext-nnu3] - uext[kext+nnu3] ) +   
	  dmu1[2] * ( uext[kext-nnu2] - uext[kext+nnu2] ) +   
	  dmu1[3] * ( uext[kext-nnu1] - uext[kext+nnu1] );
	
	ddnu2 =
	  dnu2[0] * ( uext[kext-  4] + uext[kext+  4] ) + 
	  dnu2[1] * ( uext[kext-  3] + uext[kext+  3] ) + 
	  dnu2[2] * ( uext[kext-  2] + uext[kext+  2] ) + 
	  dnu2[3] * ( uext[kext-  1] + uext[kext+  1] );
	
	ddnu1 =
	  dnu1[0] * ( uext[kext-  4] - uext[kext+  4] ) + 
	  dnu1[1] * ( uext[kext-  3] - uext[kext+  3] ) +   
	  dnu1[2] * ( uext[kext-  2] - uext[kext+  2] ) +   
	  dnu1[3] * ( uext[kext-  1] - uext[kext+  1] );
	
	diffop     = ddmu2 + b[k]*ddmu1 + ddnu2 + d[k]*ddnu1;                  
#ifdef MUTEX
	pthread_mutex_lock(&mutex4mcsor);
#endif
	uext[kext] = omega * (g[k]-diffop)/e[k] + omega1 * uext[kext];
#ifdef MUTEX
	pthread_mutex_unlock(&mutex4mcsor);
#endif
      }

      
      if ( threadsNum > 1 ) {

#ifdef TRACE3
	printf("TRACE3: mcsorpt/before barrier/colour threadNum        2 %5d \n",threadNum);
#endif

	pthread_barrier_wait(&barrier4mcsor);
#ifdef TRACE3    
	printf("TRACE3: mcsorpt/barrier passed/colour threadNum        2 %5d \n",threadNum);
#endif
      }
      
      colour=3;
      start=nstart[threadNum*5+colour];
      stop=nstop[threadNum*5+colour];
      
#ifdef TRACE3    
      printf("TRACE3: mcsorpt/threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop);
#endif
      
      for ( i = start; i <=stop; i++ ) {
	kext = indxex[i-1] - 1;
	k  = indx[i-1] - 1;
	ddmu2 =
	  dmu2[0] * ( uext[kext-nnu4] + uext[kext+nnu4] ) + 
	  dmu2[1] * ( uext[kext-nnu3] + uext[kext+nnu3] ) +   
	  dmu2[2] * ( uext[kext-nnu2] + uext[kext+nnu2] ) + 
	  dmu2[3] * ( uext[kext-nnu1] + uext[kext+nnu1] );
	
	ddmu1 =
	  dmu1[0] * ( uext[kext-nnu4] - uext[kext+nnu4] ) + 
	  dmu1[1] * ( uext[kext-nnu3] - uext[kext+nnu3] ) +   
	  dmu1[2] * ( uext[kext-nnu2] - uext[kext+nnu2] ) +   
	  dmu1[3] * ( uext[kext-nnu1] - uext[kext+nnu1] );
	
	ddnu2 =
	  dnu2[0] * ( uext[kext-  4] + uext[kext+  4] ) + 
	  dnu2[1] * ( uext[kext-  3] + uext[kext+  3] ) + 
	  dnu2[2] * ( uext[kext-  2] + uext[kext+  2] ) + 
	  dnu2[3] * ( uext[kext-  1] + uext[kext+  1] );
	
	ddnu1 =
	  dnu1[0] * ( uext[kext-  4] - uext[kext+  4] ) + 
	  dnu1[1] * ( uext[kext-  3] - uext[kext+  3] ) +   
	  dnu1[2] * ( uext[kext-  2] - uext[kext+  2] ) +   
	  dnu1[3] * ( uext[kext-  1] - uext[kext+  1] );
	
	diffop     = ddmu2 + b[k]*ddmu1 + ddnu2 + d[k]*ddnu1;                        
#ifdef MUTEX
	pthread_mutex_lock(&mutex4mcsor);
#endif
	uext[kext] = omega * (g[k]-diffop)/e[k] + omega1 * uext[kext];
#ifdef MUTEX
	pthread_mutex_unlock(&mutex4mcsor);
#endif
      }
      
      if ( threadsNum > 1 ) {
#ifdef TRACE3
	printf("TRACE3: mcsorpt/before barrier/colour threadNum        3 %5d \n",threadNum);
#endif

	pthread_barrier_wait(&barrier4mcsor);
#ifdef TRACE3    
	printf("TRACE3: mcsorpt/barrier passed/colour threadNum        3 %5d \n",threadNum);
#endif
      }
      
      colour=4;
      start=nstart[threadNum*5+colour];
      stop=nstop[threadNum*5+colour];
      
#ifdef TRACE3    
      printf("TRACE3: mcsorpt/threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop);
#endif
      
      for ( i = start; i <= stop; i++ ) {
	kext = indxex[i-1] - 1;
	k  = indx[i-1] - 1;
	ddmu2 =
	  dmu2[0] * ( uext[kext-nnu4] + uext[kext+nnu4] ) + 
	  dmu2[1] * ( uext[kext-nnu3] + uext[kext+nnu3] ) +   
	  dmu2[2] * ( uext[kext-nnu2] + uext[kext+nnu2] ) + 
	  dmu2[3] * ( uext[kext-nnu1] + uext[kext+nnu1] );
	
	ddmu1 =
	  dmu1[0] * ( uext[kext-nnu4] - uext[kext+nnu4] ) + 
	  dmu1[1] * ( uext[kext-nnu3] - uext[kext+nnu3] ) +   
	  dmu1[2] * ( uext[kext-nnu2] - uext[kext+nnu2] ) +   
	  dmu1[3] * ( uext[kext-nnu1] - uext[kext+nnu1] );
      
	ddnu2 =
	  dnu2[0] * ( uext[kext-  4] + uext[kext+  4] ) + 
	  dnu2[1] * ( uext[kext-  3] + uext[kext+  3] ) + 
	  dnu2[2] * ( uext[kext-  2] + uext[kext+  2] ) + 
	  dnu2[3] * ( uext[kext-  1] + uext[kext+  1] );
	
	ddnu1 =
	  dnu1[0] * ( uext[kext-  4] - uext[kext+  4] ) + 
	  dnu1[1] * ( uext[kext-  3] - uext[kext+  3] ) +   
	  dnu1[2] * ( uext[kext-  2] - uext[kext+  2] ) +   
	  dnu1[3] * ( uext[kext-  1] - uext[kext+  1] );
	
	diffop     = ddmu2 + b[k]*ddmu1 + ddnu2 + d[k]*ddnu1;
#ifdef MUTEX
	pthread_mutex_lock(&mutex4mcsor);
#endif
	uext[kext] = omega * (g[k]-diffop)/e[k] + omega1 * uext[kext];
#ifdef MUTEX
	pthread_mutex_unlock(&mutex4mcsor);
#endif
      }
      
      if ( threadsNum > 1 ) {
#ifdef TRACE3
	printf("TRACE3: mcsorpt/before barrier/colour threadNum        4 %5d \n",threadNum);
#endif

	pthread_barrier_wait(&barrier4mcsor);
#ifdef TRACE3    
	printf("TRACE3: mcsorpt/barrier passed/colour threadNum        4 %5d \n",threadNum);
#endif
      }
      
      // isym is orbital dependent and is being changed in the calling routine
      if (isym == 1) {
	start=nstart7[threadNum];
	stop=nstop7[threadNum];
#ifdef MUTEX
	pthread_mutex_lock(&mutex4mcsor);      
#endif
	for (i=start; i<=stop; i++) {
	  ij=indx7[i-1]-1;
	  uext[ij]=exeven[0]*uext[ij+nnu1]+exeven[1]*uext[ij+nnu2]+ 
	    exeven[2]*uext[ij+nnu3]+exeven[3]*uext[ij+nnu4]+exeven[4]*uext[ij+nnu5];
	}
#ifdef MUTEX
	pthread_mutex_unlock(&mutex4mcsor);
#endif
	
	if ( threadsNum > 1 ) {
#ifdef TRACE3    
	  printf("TRACE3: mcsorpt/barrier 6-a threadNum %2d \n",threadNum);
#endif
	  pthread_barrier_wait(&barrier4mcsor);
	}
	
	start=nstart6a[threadNum];
	stop=nstop6a[threadNum];
#ifdef MUTEX
	pthread_mutex_lock(&mutex4mcsor);
#endif
	for (i=start; i<=stop; i++) {      
	  ij=indx6a[i-1]-1;
	  uext[ij]=exeven[0]*uext[ij+1]+exeven[1]*uext[ij+2]+
	    exeven[2]*uext[ij+3]+exeven[3]*uext[ij+4]+exeven[4]*uext[ij+5];
	}
#ifdef MUTEX
	pthread_mutex_unlock(&mutex4mcsor);
#endif
      
	if ( threadsNum > 1 ) {
#ifdef TRACE3    
	  printf("TRACE3: mcsorpt/barrier 6-b threadNum %2d \n",threadNum);
#endif
	  pthread_barrier_wait(&barrier4mcsor);
	}
	
	start=nstart6b[threadNum];
	stop=nstop6b[threadNum];
#ifdef MUTEX
	pthread_mutex_lock(&mutex4mcsor);
#endif
	for (i=start; i<=stop; i++) {            
	  ij=indx6b[i-1]-1;
	  uext[ij]=exeven[0]*uext[ij-1]+exeven[1]*uext[ij-2]+
	    exeven[2]*uext[ij-3]+exeven[3]*uext[ij-4]+exeven[4]*uext[ij-5];
	}
#ifdef MUTEX
	pthread_mutex_unlock(&mutex4mcsor);
#endif
	if ( threadsNum > 1 ) {
#ifdef TRACE3    
	  printf("TRACE3: mcsorpt/barrier 6-c threadNum %2d \n",threadNum);
#endif
	  pthread_barrier_wait(&barrier4mcsor);
	}
	
      } else {
	
	start=nstart7[threadNum];
	stop=nstop7[threadNum];
	for (i=start; i<=stop; i++) {
	  ij=indx7[i-1]-1;
	  uext[ij]=0.0;
	}
	
	start=nstart6a[threadNum];
	stop=nstop6a[threadNum];
	for (i=start; i<=stop; i++) {
	  ij=indx6a[i-1]-1;
	  uext[ij]=0.0;
	}
	
	start=nstart6b[threadNum];
	stop=nstop6b[threadNum];
	for (i=start; i<=stop; i++) {
	  ij=indx6b[i-1]-1;
	  uext[ij]=0.0;
	}
      }
    }

    if (threadsNum>1) {

#ifdef TRACE3   
	printf("TRACE3: mcsorpt/end of macsor %1d \n",threadNum);
#endif

      pthread_mutex_lock(&mutex4mcsor);    
      worker_threads_done++;
#ifdef TRACE2    
      printf("TRACE2: mcsorpt/worker_threads_done++ %1d %1d \n",threadNum,worker_threads_done);
#endif

      if (worker_threads_done==threadsNum) {
#ifdef TRACE2
	printf("TRACE2: mcsorpt/thread=%1d worker_threads_done= %1d \n",thread,worker_threads_done);	
	printf("TRACE2: mcsorpt/signal work_stop  %1d \n",threadNum);
#endif
	worker_threads_done=0;
	worker_threads_waiting=0;
	pthread_cond_broadcast(&work_stop);
	pthread_mutex_unlock(&mutex4mcsor);    
	continue;
      }
#ifdef TRACE2
      printf("TRACE2: mcsorpt/wait on work-stop %1d \n",threadNum);
#endif
      pthread_cond_wait(&work_stop,&mutex4mcsor);
      pthread_mutex_unlock(&mutex4mcsor);    
    }
  }
}

