#include "catalog.h"
/*
   checks of the event weighting in the Michael inversion: an event
   with weight 2 equals the event entered twice, scaling all weights
   does not change the result, zero weight equals omitting the event,
   and the normal equation and design matrix solvers agree. returns
   the number of failures
*/
#define NQ 12

static double maxdiff(double *a, double *b)
{
  double d = 0;
  int k;
  for(k=0;k < 6;k++)
    if(fabs(a[k]-b[k]) > d)
      d = fabs(a[k]-b[k]);
  return d;
}

static int check(char *name, double *a, double *b, double tol)
{
  double d = maxdiff(a,b);
  printf("%-55s max diff %.2e %s\n",name,d,(d <= tol)?("PASS"):("FAIL"));
  return (d <= tol)?(0):(1);
}

int main(void)
{
  long int seed = -11;
  int i, j, nobs, fail = 0;
  double ang[6*(NQ+1)], ang2[6*(NQ+1)], w[NQ+1], w1[NQ+1], w2[NQ+1];
  double s1[6], s2[6], st, dp, rk, *sl, *am;

  for(i=0;i < NQ;i++){		/* random planes (auxiliary plane not needed) */
    st = ran2(&seed)*2*M_PI;
    dp = 0.2 + ran2(&seed)*1.3;
    rk = (ran2(&seed)-0.5)*2*M_PI;
    ang[6*i]   = st;   ang[6*i+1] = dp;     ang[6*i+2] = rk;
    ang[6*i+3] = st+1; ang[6*i+4] = dp*0.9; ang[6*i+5] = rk+0.5;
    w[i] = 0.2 + ran2(&seed);
  }
  /* 1: weight 2 on event 0 equals event 0 entered twice */
  for(i=0;i < NQ;i++)
    w1[i] = 1.0;
  w1[0] = 2.0;
  solve_stress_michael_specified_plane(NQ,ang,w1,s1,BC_FALSE);
  memcpy(ang2,ang,6*NQ*sizeof(double));
  memcpy(ang2+6*NQ,ang,6*sizeof(double));
  for(i=0;i <= NQ;i++)
    w2[i] = 1.0;
  solve_stress_michael_specified_plane(NQ+1,ang2,w2,s2,BC_FALSE);
  fail += check("weight 2 equals duplicated event",s1,s2,1e-12);

  /* 2: scaling all weights does not change the solution */
  for(i=0;i < NQ;i++)
    w2[i] = 3.7*w[i];
  solve_stress_michael_specified_plane(NQ,ang,w,s1,BC_FALSE);
  solve_stress_michael_specified_plane(NQ,ang,w2,s2,BC_FALSE);
  fail += check("all weights times 3.7",s1,s2,1e-12);

  /* 3: normal equations and design matrix solver agree (random weights) */
  nobs = 0;
  sl = (double *)malloc(sizeof(double)*3);
  am = (double *)malloc(sizeof(double)*15);
  for(i=0;i < NQ;i++)
    michael_assign_to_matrix(ang+6*i,&nobs,&sl,&am);
  michael_solve_lsq(5,3,nobs,am,sl,w,s2);
  fail += check("normal equations vs design matrix, random weights",s1,s2,1e-10);
  free(sl);free(am);

  /* 4: zero weight equals leaving the event out */
  for(i=0;i < NQ;i++)
    w1[i] = 1.0;
  w1[3] = 0.0;
  solve_stress_michael_specified_plane(NQ,ang,w1,s1,BC_FALSE);
  for(i=j=0;i < NQ;i++)
    if(i != 3){
      memcpy(ang2+6*j,ang+6*i,6*sizeof(double));
      w2[j] = 1.0;
      j++;
    }
  solve_stress_michael_specified_plane(NQ-1,ang2,w2,s2,BC_FALSE);
  fail += check("zero weight equals omitted event",s1,s2,1e-12);

  printf("%d failure(s)\n",fail);
  return fail;
}
