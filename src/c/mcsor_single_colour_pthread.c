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

#ifdef MUTEX
extern pthread_mutex_t mutex4mcsor;
#endif

void mcsor_single_colour_pthread_ (void *args) {                             

  struct sor_t *params4mcsor = (struct sor_t *) args;
  
  int threadsNum    = params4mcsor->threadsNum;  
  int threadNum     = params4mcsor->threadNum;
  int isym          = params4mcsor->isym;        
  int *indx         = params4mcsor->indx;
  int *indxex       = params4mcsor->indxex;  
  int *indx6a       = params4mcsor->indx6a;
  int *indx6b       = params4mcsor->indx6b;
  int *indx7        = params4mcsor->indx7; 
  
  double *b         = params4mcsor->b;
  double *d         = params4mcsor->d;          
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

#ifdef TRACE2
  printf("TRACE2: mcsor_single_colour_pthread/threadNum maxsor2 %5d %5d \n",threadNum,maxsor2);
#endif

  for ( miter=1; miter<=maxsor2; miter++) {
    colour=0;
    start=nstart[threadNum*5+colour];
    stop= nstop[threadNum*5+colour];

#ifdef TRACE2    
    printf("TRACE2: mcsor_single_colour_pthread/threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop);
#endif
    
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
      pthread_barrier_wait(&barrier4mcsor);      
#ifdef TRACE3    
      printf("TRACE3: mcsor_single_colour_pthread/barrier passed/colour threadNum        0 %5d \n",threadNum);
#endif
    }

    
    colour=1;
    start=nstart[threadNum*5+colour];
    stop=nstop[threadNum*5+colour];
#ifdef TRACE3    
    printf("TRACE2: mcsor_single_colour_pthread/threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop);
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
      pthread_barrier_wait(&barrier4mcsor);
#ifdef TRACE3    
    printf("TRACE3: mcsor_single_colour_pthread/barrier passed/colour threadNum        1 %5d \n",threadNum);
#endif
    }

    colour=2;
    start=nstart[threadNum*5+colour];
    stop=nstop[threadNum*5+colour];

#ifdef TRACE3    
    printf("TRACE2: mcsor_single_colour_pthread/threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop);
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
      pthread_barrier_wait(&barrier4mcsor);
#ifdef TRACE3    
    printf("TRACE3: mcsor_single_colour_pthread/barrier passed/colour threadNum        2 %5d \n",threadNum);
#endif
    }

    colour=3;
    start=nstart[threadNum*5+colour];
    stop=nstop[threadNum*5+colour];

#ifdef TRACE3    
    printf("TRACE2: mcsor_single_colour_pthread/threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop);
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
      pthread_barrier_wait(&barrier4mcsor);
#ifdef TRACE3    
    printf("TRACE3: mcsor_single_colour_pthread/barrier passed/colour threadNum        3 %5d \n",threadNum);
#endif
    }

    colour=4;
    start=nstart[threadNum*5+colour];
    stop=nstop[threadNum*5+colour];

#ifdef TRACE3    
    printf("TRACE2: mcsor_single_colour_pthread/threadNum colour start stop %5d %5d %5d %5d \n",threadNum,colour,start,stop);
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
      pthread_barrier_wait(&barrier4mcsor);
#ifdef TRACE3    
    printf("TRACE3: mcsor_single_colour_pthread/barrier passed/colour threadNum        4 %5d \n",threadNum);
#endif
    }
    
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
	printf("TRACE3: mcsor_single_colour_pthread/barrier 6-a threadNum %2d \n",threadNum);
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
	printf("TRACE3: mcsor_single_colour_pthread/barrier 6-b threadNum %2d \n",threadNum);
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
	printf("TRACE3: mcsor_single_colour_pthread/barrier 6-c threadNum %2d \n",threadNum);
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
      
      if ( threadsNum > 1 ) {
	pthread_barrier_wait(&barrier4mcsor);
      }
      
      start=nstart6a[threadNum];
      stop=nstop6a[threadNum];
      for (i=start; i<=stop; i++) {
	ij=indx6a[i-1]-1;
	uext[ij]=0.0;
      }
      
      if ( threadsNum > 1 ) {
	pthread_barrier_wait(&barrier4mcsor);
      }
      
      start=nstart6b[threadNum];
      stop=nstop6b[threadNum];
      for (i=start; i<=stop; i++) {
	ij=indx6b[i-1]-1;
	uext[ij]=0.0;
      }
      if ( threadsNum > 1 ) {
	pthread_barrier_wait(&barrier4mcsor);
      }
    }
  }
}
  
