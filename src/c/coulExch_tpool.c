// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2023  Jacek Kobus 

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <time.h>
#include <stdalign.h> 
#include <stdbool.h>
#include "sorpt.h"  
   //#include "tpool.h"


int threadIDs[max_threads4pots];

struct exchsor_t params[max_threads4pots];

//typedef struct exchsor_t* params[max_threads4pots];

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;  

pthread_cond_t exchsor_start[max_threads4pots];
pthread_cond_t exchsor_stop[max_threads4pots];
pthread_cond_t exchsor_ready = PTHREAD_COND_INITIALIZER;

pthread_t threads[max_threads4pots];
int work_threads;
int exchsor_threads_done;
int exchsor_threads_waiting;
bool quit;
int thread,nthreads;

extern void relax_single_pot_tpool_( void *arg );

void tpoolStart4pots ( int *num_threads4pots )
{
  //  int thread;
  int threadsNum=*(int *)num_threads4pots;

  int numofcpus = sysconf(_SC_NPROCESSORS_ONLN); // Get the number of logical CPUs.
  bool adjusted=false;  
  if ( (int) *num_threads4pots > numofcpus ) {  
    *num_threads4pots=numofcpus;
    adjusted=true;
  }

  threadsNum=*(int *)num_threads4pots;
  
  if ( (int) *num_threads4pots > max_threads4pots ) {
    printf ("%20d threads cannot be created; "
	    "increase max_threads4pots (see sortpt.h)\n",*num_threads4pots);
    return;
  }

  
  pthread_mutex_init(&mutex,NULL);

  pthread_cond_init(&exchsor_ready,NULL);
  
  if ( nthreads >= 1 ) {
    pthread_mutex_lock(&mutex);
    exchsor_threads_done=0;
    pthread_mutex_unlock(&mutex);

    pthread_cond_init(&exchsor_ready,NULL);    
    for (thread=0; thread<nthreads; thread++) {
      threadIDs[thread]=thread;
      pthread_cond_init(&exchsor_start[thread],NULL);
      pthread_cond_init(&exchsor_stop[thread],NULL);
#ifdef TRACE1      
      printf ("TRACE1/tpoolStart4pots: created thread %2d %2d\n", thread,params[thread].threadNum);
      printf ("TRACE1/tpoolStart4pots: nthreads %2d \n", nthreads);
#endif
      pthread_create( &threads[thread], NULL, (void *) relax_single_pot_tpool_, &params[thread].threadNum);
    }

    printf ("      TPOOL: %2d SORPT threads created \n",threadsNum);

    if ( adjusted ) {  
      printf ("%65s\n","number of threads adjusted to available logical CPUs");
    }
  }
  return;
}

void tpoolStop4pots_ (int *num_threads_pots )
{
  int nthreads,thread;
  nthreads=*(int *)num_threads_pots;

  pthread_mutex_lock(&mutex);  
  quit=true;

  for (thread=0; thread<nthreads; thread++) {
    pthread_cond_signal(&exchsor_start[thread]);
  }
  pthread_mutex_unlock(&mutex);    
  for (thread=0; thread<nthreads; thread++) {
    pthread_join(threads[thread], NULL); 
  };

#ifdef TRACE2  
  printf ("\n ... %2d thread(s) destroyed\n\n", (int) *num_threads_pots); 
#endif
}

void coulExch_tpool_ () 
{
  static bool threads_created=false;
  
  int thread;
  int ib1, ib2, ibexc, idel, isym, in1, in2, nexchpot;

  int *i1b=c_interface_7_.i1bc;      
  int *nexchpots=c_interface_8_.nexchpotsc;    
  int maxpots=c_interface_8_.maxpotsc;
  
  int iorb=c_interface_46_.iorbc;

  nexchpot=nexchpots[iorb-1];
  int threadsNum=nexchpot;
  threadsNum=c_interface_46_.nthreadsc;
  nthreads=maxpots;
#ifdef TRACE
  printf("TRACE/coulExch_tpool: iorb threadsNum, nthreads: %2d %2d %2d \n",iorb,threadsNum, nthreads);
#endif
 
  for (thread = 0; thread < maxpots; thread++) {
    params[thread].threadNum = thread;
    params[thread].indxex=&c_sorptr_.cw_sor[c_iadex_.iadextc-1];
    params[thread].indx=  &c_sorptr_.cw_sor[c_iadex_.iadnorc-1];    
    params[thread].indx6a=&c_sorptr_.cw_sor[c_iadex_.iadex1c-1];
    //params[thread].indx6b=c_iadex_.indx6b;
    params[thread].indx6b=&c_sorptr_.cw_sor[c_iadex_.iadex2c-1];
    //params[thread].indx7=c_iadex_.indx7;
    params[thread].indx7= &c_sorptr_.cw_sor[c_iadex_.iadex3c-1];
    //params[thread].b=b;
    params[thread].b= &c_supplptr_.cw_suppl[c_i4b_.i4barr[1]-1];
    //params[thread].d=d;
    params[thread].d= &c_supplptr_.cw_suppl[c_i4b_.i4barr[2]-1];
    params[thread].ef=&c_supplptr_.cw_suppl[c_i4b_.i4barr[3]-1];
    //params[thread].ef=ef;
    params[thread].excp=c_exchptr_.cw_exch;
    //params[thread].excp=excp;
    params[thread].f3=&c_supplptr_.cw_suppl[c_i4b_.i4barr[7]-1];
    params[thread].gf=&c_supplptr_.cw_suppl[c_i4b_.i4barr[9]-1];
    //params[thread].f3=f3;
    //params[thread].gf=gf;
    params[thread].psi=c_orbptr_.cw_orb;
    //params[thread].psi=psi;
  }  

  for (thread = 0; thread < nexchpots[iorb-1]; thread++) {
    in1 = c_interface_1_.ins1[thread][iorb-1];
    in2 = c_interface_2_.ins2[thread][iorb-1];     
    
    ib1=i1b[in1-1];
    ib2=i1b[in2-1];
    
    ibexc = c_interface_3_.ibexcpc[thread][iorb-1];
    idel  = c_interface_4_.deltam4potc[thread][iorb-1];
    //ipc   = c_interface_5_.ipcsc[thread][iorb-1];
    isym  = c_interface_6_.isymsc[thread][iorb-1];
    
    //printf ("%CM 8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",iorb, i+1, in1, in2, ibexc, idel,isym,ib1, ib2); 
    
    params[thread].threadNum=thread;            
    params[thread].threadsNum=threadsNum;            
    params[thread].in1=in1;
    params[thread].in2=in2;
    params[thread].ib1=ib1;
    params[thread].ib2=ib2;   
    params[thread].ibexc=ibexc;
    params[thread].isym=isym; 
    params[thread].deltam2=(double) idel*idel;        
  }
 
  pthread_mutex_lock(&mutex);
  exchsor_threads_done=0;
  pthread_mutex_unlock(&mutex);

  if (! threads_created) {
    tpoolStart4pots(&nthreads);
    threads_created=true;
    usleep(500000);
  }
  
  for ( thread=0; thread<nexchpots[iorb-1]; thread++ ) {
    pthread_mutex_lock(&mutex);
    pthread_cond_signal(&exchsor_start[thread]);
    pthread_mutex_unlock(&mutex);
#ifdef TRACE1
    printf("TRACE1: coulExch_tpool/dispatch work for thread %2d\n",thread);  
#endif
    //}
  }
  
  pthread_mutex_lock(&mutex);
  pthread_cond_wait(&exchsor_ready,&mutex);
  pthread_mutex_unlock(&mutex);

  return;
}



