
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

extern void relax_single_pot_pthread_ ( void *args );

// symbols e and g must be renamed since they are used in sor4exch 

void coulExch_pthread_ ()
{

  int ib1, ib2, ibexc, idel, isym, in1, in2, nexchpot, thread,  threadsNum;
  int *nexchpots=c_interface_8_.nexchpotsc;    

  int iorb=c_interface_46_.iorbc;
  
  nexchpot=nexchpots[iorb-1];
  threadsNum=nexchpot;
  pthread_t threads[threadsNum];
  struct exchsor_t params[threadsNum];

  for (thread = 0; thread < threadsNum; thread++) {

    params[thread].threadsNum = threadsNum;
    params[thread].threadNum = thread;

    //params[thread].indx=iadnor;
    //params[thread].indxex=iadext;
    //params[thread].indx6a=indx6a;
    //params[thread].indx6b=indx6b;
    //params[thread].indx7=indx7;
    params[thread].indxex=&c_sorptr_.cw_sor[c_iadex_.iadextc-1];
    params[thread].indx=  &c_sorptr_.cw_sor[c_iadex_.iadnorc-1];    
    params[thread].indx6a=&c_sorptr_.cw_sor[c_iadex_.iadex1c-1];
    params[thread].indx6b=&c_sorptr_.cw_sor[c_iadex_.iadex2c-1];
    params[thread].indx7= &c_sorptr_.cw_sor[c_iadex_.iadex3c-1];

    //params[thread].b=b;
    //params[thread].d=d;
    //params[thread].ef=ef;
    //params[thread].excp=excp;
    //params[thread].f3=f3;
    //params[thread].gf=gf;
    //params[thread].psi=psi;

    params[thread].b=&c_supplptr_.cw_suppl[c_i4b_.i4barr[1]-1];
    params[thread].d=&c_supplptr_.cw_suppl[c_i4b_.i4barr[2]-1];
    params[thread].ef=&c_supplptr_.cw_suppl[c_i4b_.i4barr[3]-1];
    params[thread].excp=c_exchptr_.cw_exch;
    params[thread].f3=&c_supplptr_.cw_suppl[c_i4b_.i4barr[7]-1];
    params[thread].gf=&c_supplptr_.cw_suppl[c_i4b_.i4barr[9]-1];
    params[thread].psi=c_orbptr_.cw_orb;

  }  

 if ( threadsNum == 1 ) {
   thread=0;
   nexchpot=0;
   in1 = c_interface_1_.ins1[nexchpot][iorb-1];
   in2 = c_interface_2_.ins2[nexchpot][iorb-1];     
   
   ib1=c_interface_7_.i1bc[in1-1];
   ib2=c_interface_7_.i1bc[in2-1];      
   
   ibexc = c_interface_3_.ibexcpc[nexchpot][iorb-1];
   idel  = c_interface_4_.deltam4potc[nexchpot][iorb-1];
   isym  = c_interface_6_.isymsc[nexchpot][iorb-1];
   
   params[thread].in1=in1;
   params[thread].in2=in2;
   params[thread].ib1=ib1;
   params[thread].ib2=ib2;   
   params[thread].ibexc=ibexc;
   params[thread].deltam2=(double) idel*idel;
   params[thread].isym=isym;
   params[thread].threadNum=thread;
   relax_single_pot_pthread_ ((void *) &params[thread]);
 }
 else
   {
     for ( thread=0; thread<nexchpots[iorb-1]; thread++ ) {

       in1 = c_interface_1_.ins1[thread][iorb-1];
       in2 = c_interface_2_.ins2[thread][iorb-1];     

       ib1=c_interface_7_.i1bc[in1-1];
       ib2=c_interface_7_.i1bc[in2-1];      
       
       ibexc = c_interface_3_.ibexcpc[thread][iorb-1];
       idel  = c_interface_4_.deltam4potc[thread][iorb-1];
       //ipc   = c_interface_5_.ipcsc[thread][iorb-1];
       isym  = c_interface_6_.isymsc[thread][iorb-1];

       params[thread].in1=in1;
       params[thread].in2=in2;
       params[thread].ib1=ib1;
       params[thread].ib2=ib2;   
       params[thread].ibexc=ibexc;
       params[thread].isym=isym; 
       params[thread].threadNum=thread;
       params[thread].deltam2=(double) idel*idel;        
     }

     for (thread = 0; thread < threadsNum; thread++) {
	pthread_create (&threads[thread], NULL, (void *) relax_single_pot_pthread_,((void *) &params[thread]));
     }

     for (thread = 0; thread < threadsNum; thread++) {
       int status;
       pthread_join(threads[thread], (void **) &status);
     }
   }
}


