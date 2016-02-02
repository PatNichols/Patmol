#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_
#include <time.h>
#include <sys/time.h>

#ifdef _USE_RT_

typedef struct { 
double acc; 
struct timespec ts;
struct timespec tf;
} stopwatch_t;

#define stopwatch_stop(t)\
  {	clock_gettime(CLOCK_MONOTONIC,(t)->tf);\
    (t)->acc += (double)((t)->tf.tv_sec - (t)->ts.tv_sec) + 1.e-9 * ( (t)->tf.tv_nsec - (t)->ts.tv_nsec);\
  }

#define stopwatch_start(t) {\
    clock_gettime(CLOCK_MONOTONIC,(t)->ts);\
  }

#else

typedef struct { 
double acc; 
struct timeval ts;
struct timeval tf;
} stopwatch_t;

#define stopwatch_stop(t)\
  {	gettimeofday(&((t)->tf),0x0);\
    (t)->acc += (double)((t)->tf.tv_sec - (t)->ts.tv_sec) + 1.e-6 * ( (t)->tf.tv_usec - (t)->ts.tv_usec);\
  }

#define stopwatch_start(t) { gettimeofday(&((t)->ts),0x0);\
  }

#endif

#define stopwatch_clear(t) (t)->acc=0.
#define stopwatch_elapsed_time(t) (t)->acc

#endif
