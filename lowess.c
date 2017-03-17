// Taken from: https://github.com/blindglobe/common-lisp-stat/blob/master/lib/lowess.c

/*
  Translated from RATFOR lowess code of W. S. Cleveland as obtained from NETLIB
*/

#include "math.h"
#include <stdlib.h>
#include <stdio.h>

#define FALSE 0
#define TRUE 1

long
min (long x, long y)
{
  return((x < y) ? x : y);
}

long
max (long x, long y)
{
  return((x > y) ? x : y);
}

static double pow2(double x) { return(x * x); }
static double pow3(double x) { return(x * x * x); }
/* static */ double fmax(double x, double y) { return (x > y ? x : y); }

static int compar(const void *aa, const void *bb)
{
  const double* a = (double*)aa;
  const double* b = (double*)bb;

  if (*a < *b) return(-1);
  else if (*a > *b) return(1);
  else return(0);
}

static void
lowest(double *x, double *y, size_t n, double xs, double *ys, long nleft, long nright,
       double *w, int userw, double *rw, int *ok)
{
  double range, h, h1, h9, a, b, c, r;
  long j, nrt;

  range = x[n - 1] - x[0];
  h = fmax(xs - x[nleft], x[nright] - xs);
  h9 = .999 * h;
  h1 = .001 * h;

  /* compute weights (pick up all ties on right) */
  a = 0.0;        /* sum of weights */
  for(j = nleft; j < n; j++) {
    w[j]=0.0;
    r = fabs(x[j] - xs);
    if (r <= h9) {    /* small enough for non-zero weight */
      if (r > h1) w[j] = pow3(1.0-pow3(r/h));
      else w[j] = 1.0;
      if (userw) w[j] = rw[j] * w[j];
      a += w[j];
    }
    else if (x[j] > xs) break;  /* get out at first zero wt on right */
  }
  nrt = j - 1;  /* rightmost pt (may be greater than nright because of ties) */
  if (a <= 0.0) *ok = FALSE;
  else { /* weighted least squares */
    *ok = TRUE;

    /* make sum of w[j] == 1 */
    for (j = nleft; j <= nrt; j++) w[j] = w[j] / a;

    if (h > 0.0) {     /* use linear fit */

      /* find weighted center of x values */
      for (j = nleft, a = 0.0; j <= nrt; j++) a += w[j] * x[j];

      b = xs - a;
      for (j = nleft, c = 0.0; j <= nrt; j++) 
	c += w[j] * (x[j] - a) * (x[j] - a);

      if(sqrt(c) > .001 * range) {
	/* points are spread out enough to compute slope */
	b = b/c;
	for (j = nleft; j <= nrt; j++) 
	  w[j] = w[j] * (1.0 + b*(x[j] - a));
      }
    }
    for (j = nleft, *ys = 0.0; j <= nrt; j++) *ys += w[j] * y[j];
  }
}

static void
sort(double *x, size_t n)
{
  extern void qsort();

  qsort(x, n, sizeof(double), compar);
}



int
lowess(double *x, double *y, size_t n,
       double f, size_t nsteps,
       double delta, double *ys, double *rw, double *res)
{
  int iter, ok;
  long i, j, last, m1, m2, nleft, nright, ns;
  double d1, d2, denom, alpha, cut, cmad, c9, c1, r;
  
  if (n < 2) { ys[0] = y[0]; return(1); }
  ns = max(min((long) (f * n), n), 2);  /* at least two, at most n points */
  for(iter = 1; iter <= nsteps + 1; iter++){      /* robustness iterations */
    nleft = 0; nright = ns - 1;
    last = -1;        /* index of prev estimated point */
    i = 0;   /* index of current point */
    do {
      while(nright < n - 1){
	/* move nleft, nright to right if radius decreases */
	d1 = x[i] - x[nleft];
	d2 = x[nright + 1] - x[i];
	/* if d1 <= d2 with x[nright+1] == x[nright], lowest fixes */
	if (d1 <= d2) break;
	/* radius will not decrease by move right */
	nleft++;
	nright++;
      }
      lowest(x, y,
	     n, x[i],
	     &ys[i],
	     nleft, nright,
	     res, (iter > 1), rw, &ok);
      /* fitted value at x[i] */
      if (! ok) ys[i] = y[i];
      /* all weights zero - copy over value (all rw==0) */
      if (last < i - 1) { /* skipped points -- interpolate */
	denom = x[i] - x[last];    /* non-zero - proof? */
	for(j = last + 1; j < i; j = j + 1){
	  alpha = (x[j] - x[last]) / denom;
	  ys[j] = alpha * ys[i] + (1.0 - alpha) * ys[last];
	}
      }
      last = i;        /* last point actually estimated */
      cut = x[last] + delta;     /* x coord of close points */
      for(i=last + 1; i < n; i++) {     /* find close points */
	if (x[i] > cut) break;     /* i one beyond last pt within cut */
	if(x[i] == x[last]) {      /* exact match in x */
	  ys[i] = ys[last];
	  last = i;
	}
      }
      i = max(last + 1,i - 1);
      /* back 1 point so interpolation within delta, but always go forward */
    } while(last < n - 1);
    for (i = 0; i < n; i++)      /* residuals */
      res[i] = y[i] - ys[i];
    if (iter > nsteps) break; /* compute robustness weights except last time */
    for (i = 0; i < n; i++) 
      rw[i] = fabs(res[i]);
    sort(rw,n);
    m1 = 1 + n / 2; m2 = n - m1 + 1;
    cmad = 3.0 * (rw[m1] + rw[m2]);      /* 6 median abs resid */
    c9 = .999 * cmad; c1 = .001 * cmad;
    for (i = 0; i < n; i++) {
      r = fabs(res[i]);
      if(r <= c1) rw[i] = 1.0;      /* near 0, avoid underflow */
      else if(r > c9) rw[i] = 0.0;  /* near 1, avoid underflow */
      else rw[i] = pow2(1.0 - pow2(r / cmad));
    }
  }
  return(0);
}

// To test this file on its own, set the #if 0 to #if 1, and compile and run with
// 
// gcc -Wall -std=gnu99 -lm lowess.c -o lowess && ./lowess | python
//
//
#if 0
int main()
{
    const long n = 200;

    double *x = (double *) malloc(sizeof(double)*n);
    double *y = (double *) malloc(sizeof(double)*n);

    for(long i = 0; i < n; i++) {
        x[i] = i*2*3.14/n;
        y[i] = sin(x[i]) + (1.0*rand()/(RAND_MAX + 1.0));
    }

    const double f = 0.25;
    const size_t nsteps = 3;

    const double delta = 0.3; // fixme what is this really?

    double *ys  = (double *) malloc(sizeof(double)*n);
    double *rw  = (double *) malloc(sizeof(double)*n); // what are these?
    double *res = (double *) malloc(sizeof(double)*n); // what are these?

    lowess(x, y, n,
           f, nsteps,
           delta, ys, rw, res);

    printf("x = [");
    for(long i = 0; i < n; i++) {
        printf("%f, ", x[i]);
    }
    printf("]\n");

    printf("y = [");
    for(long i = 0; i < n; i++) {
        printf("%f, ", y[i]);
    }
    printf("]\n");

    printf("ys = [");
    for(long i = 0; i < n; i++) {
        printf("%f, ", ys[i]);
    }
    printf("]\n");

    printf("import pylab as plt\n");
    printf("plt.plot(x, y, label='xy')\n");
    printf("plt.plot(x, ys, label='xys')\n");
    printf("plt.legend()\n");
    printf("plt.show()\n");

    free(x);
    free(y);
    free(ys);
    free(rw);
    free(res);

    return 0;
}
#endif
