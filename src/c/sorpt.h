// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2023  Jacek Kobus 

#ifndef PTHREADADD
#define PTHREADADD 

// maxorb must be equal to maxorb parameter (see params.f90)
// x2dhfctl make sure that this is the case
#define maxorb 36
#define max_threads4pots 45
#define max_threads4mcsor 16

struct exchsor_t {
  int threadsNum;
  int threadNum;
  int *indx;
  int *indxex;
  int *indx6a;
  int *indx6b;
  int *indx7;  

  double *b;
  double *d;
  double *uext;
  double *psi;
  double *excp;
  double *ef;  
  double *f3;  
  double *gf;
  double deltam2;
  
  int ib1;
  int ib2;   
  int ibexc; 
  int in1;
  int in2;  
  int isym;
  //int *iorb;
};

struct sor_t {
  //pthread_barrier_t barrier4mcsor;
  int threadsNum;
  int threadNum;
  int isym;

  int *indx;
  int *indxex;
  int *indx6a;
  int *indx6b;
  int *indx7;  

  double *b;
  double *d;
  double *e;
  double *g;
  double *uext;
};

//extern struct commonBlock_0_t c_interface_0_;


struct commonBlock_1_t {
    int ins1[2*maxorb][maxorb];
};
extern struct commonBlock_1_t c_interface_1_;

struct commonBlock_2_t {
    int ins2[2*maxorb][maxorb];
};
extern struct commonBlock_2_t c_interface_2_;

struct commonBlock_3_t {
    int ibexcpc[2*maxorb][maxorb];
};
extern struct commonBlock_3_t c_interface_3_;

struct commonBlock_4_t {
    int deltam4potc[2*maxorb][maxorb];
};
extern struct commonBlock_4_t c_interface_4_;

extern struct commonBlock_5_t c_interface_5_;

struct commonBlock_6_t {
    int isymsc[2*maxorb][maxorb];
};
extern struct commonBlock_6_t c_interface_6_;

struct commonBlock_7_t {
    int i1bc[maxorb];
};
extern struct commonBlock_7_t c_interface_7_;

struct commonBlock_8_t {
  int nexchpotsc[maxorb];
  int maxpotsc; 
};
extern struct commonBlock_8_t c_interface_8_;

struct commonBlock_9_t {
    int nstartc[5*max_threads4mcsor];
};
extern struct commonBlock_9_t c_interface_9_;

struct commonBlock_10_t {
    int nstopc[5*max_threads4mcsor];
};
extern struct commonBlock_10_t c_interface_10_;

struct commonBlock_11_t {
    int nstart6ac[5*max_threads4mcsor];
};
extern struct commonBlock_11_t c_interface_11_;

struct commonBlock_12_t {
    int nstop6ac[5*max_threads4mcsor];
};
extern struct commonBlock_12_t c_interface_12_;

struct commonBlock_13_t {
    int nstart6bc[5*max_threads4mcsor];
};
extern struct commonBlock_13_t c_interface_13_;

struct commonBlock_14_t {
    int nstop6bc[5*max_threads4mcsor];
};
extern struct commonBlock_14_t c_interface_14_;

struct commonBlock_15_t {
    int nstart7c[5*max_threads4mcsor];
};
extern struct commonBlock_15_t c_interface_15_;

struct commonBlock_16_t {
    int nstop7c[5*max_threads4mcsor];
};
extern struct commonBlock_16_t c_interface_16_;

struct commonBlock_17_t {
  int mxnmuc;
  int mxsizec;
  int ngrd6ac;
  int ngrd6bc;
  int ngrd7c;
  int nnu1c;
  int nnu2c;
  int nnu3c;
  int nnu4c;
  int nnu5c;
  int isstartc;
  int isstopc;
  int maxsor1c;
  int maxsor2c;
};
extern struct commonBlock_17_t c_interface_17_;



struct commonBlock_18_t {
  int isstartc;
  int isstopc;
  int maxsor1c;
  int maxsor2c;
};

struct commonBlock_20_t {
  int msleep;
};
extern struct commonBlock_20_t c_interface_20_;

struct commonBlock_40_t {
  double dmu1c[4];
};
extern struct commonBlock_40_t c_interface_40_;

struct commonBlock_41_t {
  double dmu2c[4];
};
extern struct commonBlock_41_t c_interface_41_;

struct commonBlock_42_t {
  double dnu1c[4];
};
extern struct commonBlock_42_t c_interface_42_;

struct commonBlock_43_t {
  double dnu2c[4];
};
extern struct commonBlock_43_t c_interface_43_;

struct commonBlock_44_t {
  double exevenc[5];
};
extern struct commonBlock_44_t c_interface_44_;

struct commonBlock_45_t {
  double omegac;
  double omega1c;
};
extern struct commonBlock_45_t c_interface_45_;


struct commonBlock_46_t {
  int iorbc;
  int isymc;
  int nthreadsc;
};
extern struct commonBlock_46_t c_interface_46_;


struct commonBlock_47_t {
  void *cw_sor4p;
};
extern struct commonBlock_47_t c_interface_47_;


struct cb_1_t {
  double *cw_orb;
};
extern struct cb_1_t c_orbptr_;

struct cb_2_t {
  double *cw_exch;
};
extern struct cb_2_t c_exchptr_;

struct cb_3_t {
  double *cw_suppl;
};
extern struct cb_3_t c_supplptr_;

struct cb_5_t {
  int *cw_sor;
};
extern struct cb_5_t c_sorptr_;

struct cb_10_t {
  int i4barr[20];
};
extern struct cb_10_t c_i4b_;

struct cb_11_t {
  int iadex1c;
  int iadex2c;
  int iadex3c;
  int iadextc;
  int iadnorc;
};
extern struct cb_11_t c_iadex_;

#endif

