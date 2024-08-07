
// SPDX-License-Identifier: GPL-2.0-or-later
// Copyright (C) 2023  Jacek Kobus 

#ifndef __TPOOL4T_H__
#define __TPOOL4T_H__

#include <stdbool.h>
#include <stddef.h>

struct tpool;
typedef struct tpool tpool_t;

typedef void (*thread_func_t)(void *arg);

tpool_t *tpool_create(int num);
void tpool_destroy(tpool_t *tm);

bool tpool_add_work(tpool_t *tm, thread_func_t func, void *arg);
void tpool_wait(tpool_t *tm);

int tpoolstart_ (int *num_threads_coulexch);
void tpoolstop_();
void tpoolstop4mcsor_ (int *num_threads4mcsor);

typedef struct tpool_work tpool_work_t;

extern struct exchsor_t params[max_threads4pots];
extern struct sor_t params4sor[max_threads4mcsor];

extern void mcsor_single_colour_tpool_ (void *args);

struct tpool {
  tpool_work_t    *work_first;
  tpool_work_t    *work_last;
  pthread_mutex_t  work_mutex;
  pthread_cond_t   work_cond;
  pthread_cond_t   working_cond;
  int              working_cnt;
  int              thread_cnt;
  bool             stop;
};

struct tpool_work {
  thread_func_t     func;
  void              *arg;
  struct tpool_work *next;
};


extern tpool_t *tm;
extern int threadsNum;

extern int threadIDs[max_threads4pots];

extern pthread_barrier_t barrier;

struct tpool;
typedef struct tpool tpool_t;
typedef void (*thread_func_t)(void *arg);

//bool tpool_add_work(tpool_t *tm, thread_func_t func, void *arg);
void tpool_wait(tpool_t *tm);

#endif

