// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2023  Jacek Kobus 

   /*
     mcsorptdrv is a driver routine that controls creating p-threads needed
     to perform SOR relaxations using a multi-threaded version of the MCSOR
     algorithm. 
   */

#include <stdlib.h>
#include <stdio.h>
#include <sched.h>
#include <unistd.h>
#include <stdalign.h> 
#include <stdbool.h>
#include <time.h>
#include <pthread.h>
#include "sorpt.h"

//struct sor_t params4sor[max_num_threads_mcsor];

pthread_barrier_t barrier4mcsor;


#ifdef MUTEX
pthread_mutex_t mutex4mcsor = PTHREAD_MUTEX_INITIALIZER;
#endif


extern void mcsor_single_colour_pthread_ (void *args);

void mcsor_pthread_ (int *isym, int *nthreads, double *wk2, double *lhs, double *rhs)
{
  int thread;
  int threadsNum = *nthreads;
  pthread_t threads[threadsNum];
  struct sor_t params[threadsNum];

  unsigned int count=threadsNum;
  pthread_barrier_init(&barrier4mcsor, NULL, count);
  for (thread = 0; thread < threadsNum; thread++) {
    params[thread].threadsNum = threadsNum;
    params[thread].threadNum = thread;
    //params[thread].isym=*isym;
    params[thread].isym= *(int *) isym;
    //params[thread].indx=iadnor;
    //params[thread].indxex=iadext;
    //params[thread].indx6a=indx6a;
    //params[thread].indx6b=indx6b;
    //params[thread].indx7=indx7;
    //params[thread].b=b;
    //params[thread].d=d;
    params[thread].e=lhs;
    params[thread].g=rhs;
    params[thread].uext=wk2;

    params[thread].indxex=&c_sorptr_.cw_sor[c_iadex_.iadextc-1];
    params[thread].indx=&c_sorptr_.cw_sor[c_iadex_.iadnorc-1];    
    params[thread].indx6a=&c_sorptr_.cw_sor[c_iadex_.iadex1c-1];
    params[thread].indx6b=&c_sorptr_.cw_sor[c_iadex_.iadex2c-1];
    params[thread].indx7=&c_sorptr_.cw_sor[c_iadex_.iadex3c-1];
    params[thread].b=&c_supplptr_.cw_suppl[c_i4b_.i4barr[0]-1];
    params[thread].d=&c_supplptr_.cw_suppl[c_i4b_.i4barr[2]-1];
  }

  if ( threadsNum == 1 ) {
    mcsor_single_colour_pthread_ ( (void *)&params[0] );
  }
  else {
#ifdef TRACE
      printf("TRACE:  mcsor_pthread/barrier_count = %d\n",count);
#endif
      
      for (thread = 0; thread < threadsNum; thread++) {
	//for (thread = 0; thread < 1; thread++) {      
	pthread_create(&threads[thread], NULL, (void *) mcsor_single_colour_pthread_, (void *) &params[thread]);
#ifdef TRACE2
      printf("TRACE2: mcsor_pthread/pthread_created= %d\n",thread);
#endif
    }

      for (thread = 0; thread < threadsNum; thread++) {
	int status;
#ifdef TRACE3
	printf("TRACE3: mcsor_pthread/pthread_join= %d\n",thread);
#endif
	pthread_join(threads[thread], (void **) &status);
      }
  }
  pthread_barrier_destroy(&barrier4mcsor);            
}


