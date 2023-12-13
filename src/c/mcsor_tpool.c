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
#include "tpool.h"  
//#include "ISO_Fortran_binding.h"

extern void mcsor_single_colour_tpool_ (void *args);
pthread_mutex_t mutex4mcsor;
struct sor_t params4mcsor[max_threads4mcsor];

pthread_barrier_t barrier4mcsor;
pthread_mutex_t mcsor = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t work_start = PTHREAD_COND_INITIALIZER;
pthread_cond_t work_stop  = PTHREAD_COND_INITIALIZER;
pthread_cond_t main_ready = PTHREAD_COND_INITIALIZER;

int worker_threads_done;
int worker_threads_waiting;
bool quit4mcsor;
int nthreads4mcsor;

pthread_t threads4mcsor[max_threads4mcsor];
int work_threads4mcsor;
int nthreads_mcsor;

void tpoolStart4mcsor_ ( int *num_threads4mcsor)
{
  //pthread_attr_t attr;
  int thread;
 
  int threadsNum=*(int *)num_threads4mcsor;
  nthreads_mcsor=threadsNum;
  nthreads4mcsor=threadsNum;
  if (threadsNum==1) {
    return;
  }

  int numofcpus = sysconf(_SC_NPROCESSORS_ONLN); // Get the number of logical CPUs.
  bool adjusted=false;  
  if ( (int) *num_threads4mcsor > numofcpus ) {  
    *num_threads4mcsor=numofcpus;
    adjusted=true;
  }

  if ( (int) *num_threads4mcsor > max_threads4mcsor ) {
    printf ("%20d threads cannot be created; "
	    "increase max_threads4mcsor (see sortpt.h)\n",*num_threads4mcsor);
    return;
  }

  if ( nthreads_mcsor > 1 ) {
  
    unsigned int count=nthreads_mcsor;
    pthread_barrier_init(&barrier4mcsor, NULL, count);

    pthread_mutex_init(&mutex4mcsor,NULL);
    pthread_cond_init(&work_start,NULL);
    pthread_cond_init(&work_stop,NULL);
    pthread_cond_init(&main_ready,NULL);

#ifdef TRACE
    printf("TRACE: tpoolstart4mcsor: barrier4mcsor %p \n",&barrier4mcsor); 
    //printf("TRACE: mcsorptdrv2:      count         %2d\n",count);
#endif

    //pthread_attr_init(&attr);
    //pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_mutex_lock(&mutex4mcsor);
    worker_threads_done=0;
    worker_threads_waiting=0;
    pthread_mutex_unlock(&mutex4mcsor);
    
    for (thread=0; thread<nthreads_mcsor; thread++) {
      threadIDs[thread]=thread;
      pthread_create( &threads4mcsor[thread], NULL, (void *)mcsor_single_colour_tpool_,  &params4mcsor[thread]);
    }

    printf ("      TPOOL: %2d MCSORPT threads created \n",*num_threads4mcsor);
    if ( adjusted ) {  
      printf ("%65s\n","number of threads adjusted to available logical CPUs");
    }
  }
  return;
}

void tpoolStop4mcsor_ (int *num_threads_mcsor )
{
  int thread;
  nthreads_mcsor=*(int *)num_threads_mcsor;

  if ( nthreads_mcsor == 1 ) {
    return;
  }
       
#ifdef TRACE2
    printf("end: all work done; winding up ... \n");  
    printf("end: work_threads4mcsor=%1d\n",nthreads_mcsor); 
    printf("end: signaling threads to return \n");  
#endif   

    pthread_mutex_lock(&mutex4mcsor);  
    quit4mcsor=true;                // If it's a round after job_rounds
                                    // have been already performed, tell worker threads to quit.
    pthread_cond_broadcast(&work_start);
    pthread_mutex_unlock(&mutex4mcsor);

    //printf("nthreads_mcsor: %1d\n",nthreads_mcsor);
    for (thread=0; thread<nthreads_mcsor; thread++) {
      //printf("join: %1d\n",thread);
      pthread_join(threads4mcsor[thread], NULL);
    };
    printf ("      TPOOL: %2d MCSORPT threads destroyed\n", (int) *num_threads_mcsor); 

    //pthread_barrier_destroy(&barrier4mcsor);
#ifdef TRACE2
    printf ("   TPOOL: barrier4mcsor destroyed\n");
#endif
}

void mcsor_tpool_ ( int *isym,  int *nthreads, double *wk2, double *lhs, double *rhs)
{
  static bool threads_created4mcsor=false;
  int thread;
  int i,max;
  int threadsNum = *nthreads;
  nthreads_mcsor= *nthreads;
  
  if (! threads_created4mcsor) {
    for (thread=0; thread < threadsNum; thread++) {
      params4mcsor[thread].threadsNum = threadsNum;
      params4mcsor[thread].threadNum = thread;
      //params4mcsor[thread].isym=isym;
      params4mcsor[thread].isym= *(int *) isym;
      //params4mcsor[thread].indx=iadnor;
      //params4mcsor[thread].indxex=iadext;
      //params4mcsor[thread].indx6a=indx6a;
      //params4mcsor[thread].indx6b=indx6b;
      //params4mcsor[thread].indx7=indx7;
      //params4mcsor[thread].b=b;
      //params4mcsor[thread].d=d;
      params4mcsor[thread].e=lhs;
      params4mcsor[thread].g=rhs;
      params4mcsor[thread].uext=wk2;

      //      params4mcsor[thread].indxex=&c_sorptr_.cw_sor[c_iadex_.iadextc-1];
      //params4mcsor[thread].indx=&c_sorptr_.cw_sor[c_iadex_.iadnorc-1];    
      //params4mcsor[thread].indx6a=&c_sorptr_.cw_sor[c_iadex_.iadex1c-1];
      //params4mcsor[thread].indx6b=&c_sorptr_.cw_sor[c_iadex_.iadex2c-1];
      //params4mcsor[thread].indx7=&c_sorptr_.cw_sor[c_iadex_.iadex3c-1];
      //params4mcsor[thread].b=&c_supplptr_.cw_suppl[c_i4b_.i4barr[0]-1];
      //params4mcsor[thread].d=&c_supplptr_.cw_suppl[c_i4b_.i4barr[2]-1];
    }
    
    tpoolStart4mcsor_ ( &threadsNum);
    threads_created4mcsor=true;
  }
  if ( threadsNum == 1 ) {
    mcsor_single_colour_tpool_ ( ((void *) &params4mcsor[0]));
  } else {

    for (thread=0; thread < threadsNum; thread++) {
      pthread_mutex_lock(&mutex4mcsor);
      worker_threads_done=0;
      pthread_mutex_unlock(&mutex4mcsor);
#ifdef TRACE
      printf("TRACE/mcsor_tpool: dispatch work for thread %2d \n",thread);  
#endif
      
      pthread_mutex_lock(&mutex4mcsor);
      if (worker_threads_done==0) {
#ifdef TRACE2
	printf("dispatch: broadcast new_work_ready \n");
#endif
	pthread_cond_broadcast(&work_start);
      }
      pthread_mutex_unlock(&mutex4mcsor);
      
      pthread_mutex_lock(&mutex4mcsor);
#ifdef TRACE2
      printf("mcsor_tpool: wait main_ready \n");  
#endif
      pthread_cond_wait(&main_ready,&mutex4mcsor);
      pthread_mutex_unlock(&mutex4mcsor);
    }
  }
}



