#include "catalog.h"




/* 
   
   this is the part of the code which depends on GMT psmeca utilities
   
   these used to be linked, now included here. proceed with caution!



*/

/* 
   convert a moment tensor to double couple and compute the best fit fault planes
   
   return in radians

*/
void tensor2fpangle(BC_CPREC *m, BC_CPREC *strike1, BC_CPREC *dip1, BC_CPREC *rake1,
		    BC_CPREC *strike2, BC_CPREC *dip2, BC_CPREC *rake2)
{
  /* meca stuff */
  struct M_TENSOR mt;		/* GMT meca convention */
  struct AXIS t_axis, p_axis,n_axis;
  struct nodal_plane np1,np2;
  mt.expo=1;
  mt.f[0] = m[BC_RR];mt.f[1] = m[BC_TT];mt.f[2] = m[BC_PP];
  mt.f[3] = m[BC_RT];mt.f[4] = m[BC_RP];mt.f[5] = m[BC_TP];

  GMT_momten2axe(mt,&t_axis,&n_axis,&p_axis);
  axe2dc(t_axis, p_axis, &np1, &np2); /* axes to double couple nodal planes */

  *strike1 = np1.str;*dip1 = np1.dip;*rake1 = np1.rake;
  *strike2 = np2.str;*dip2 = np2.dip;*rake2 = np2.rake;
}

/* 

   below modified from GMT 4.5.18 - their copyright 

*/
/***************************************************************************************/
void GMT_momten2axe(struct M_TENSOR mt,struct AXIS *T,struct AXIS *N,struct AXIS *P) {
  /* This version uses GMT_jacobi and does not suffer from the convert_matrix bug */
  int j, nrots;
  int np = 3;
  double *a, *d, *b, *z, *v;
  double az[3], pl[3];

  /* should check, but these are tiny */
  a = (double *) malloc ( (size_t)(np*np)* sizeof(double));
  d = (double *) malloc ( (size_t)(np) * sizeof(double));
  b = (double *) malloc ( (size_t)(np) * sizeof(double));
  z = (double *) malloc ( (size_t)(np) * sizeof(double));
  v = (double *) malloc ( (size_t)(np*np) * sizeof(double));

  a[0]=mt.f[0];	a[1]=mt.f[3];	a[2]=mt.f[4];
  a[3]=mt.f[3];	a[4]=mt.f[1];	a[5]=mt.f[5];
  a[6]=mt.f[4];	a[7]=mt.f[5];	a[8]=mt.f[2];

  if (GMT_jacobi (a, &np, &np, d, v, b, z, &nrots))
    fprintf(stderr," GMT_momten2axe: Eigenvalue routine failed to converge in 50 sweeps.\n");

  for (j = 0; j < np; j++) {
    pl[j] = asin(-v[j*np]);
    az[j] = atan2(v[j*np+2], -v[j*np+1]);
    if (pl[j] <= 0.) {
      pl[j] = -pl[j];
      az[j] += M_PI;
    }
    if (az[j] < 0)
      az[j] += TWO_PI;
    else if (az[j] > TWO_PI)
      az[j] -= TWO_PI;
  }
  T->val = d[0];	T->e = mt.expo;	T->str = az[0]; T->dip = pl[0];
  N->val = d[1];	N->e = mt.expo; N->str = az[1]; N->dip = pl[1];
  P->val = d[2];	P->e = mt.expo; P->str = az[2]; P->dip = pl[2];

  free(a);free(d);free(b);free(z);free(v);
}


/***************************************************************************************/
/*
  Calculate double couple from principal axes.
  Angles are in degrees.

  Genevieve Patau, 16 juin 1997
*/

void axe2dc(struct AXIS T,struct AXIS P,
	    struct nodal_plane *NP1,struct nodal_plane *NP2)

{

  double pp, dp, pt, dt;
  double p1, d1, p2, d2;
  double cdp, sdp, cdt, sdt;
  double cpt, spt, cpp, spp;
  double amz, amy, amx;
  double im;

  pp = P.str ; dp = P.dip;
  pt = T.str ; dt = T.dip;

  sincos (dp, &sdp, &cdp);
  sincos (dt, &sdt, &cdt);
  sincos (pt, &spt, &cpt);
  sincos (pp, &spp, &cpp);

  cpt *= cdt; spt *= cdt;
  cpp *= cdp; spp *= cdp;

  amz = sdt + sdp; amx = spt + spp; amy = cpt + cpp;
  d1 = atan2(sqrt(amx*amx + amy*amy), amz);
  p1 = atan2(amy, -amx);
  if (d1 > M_PI_2) {
    d1 = M_PI - d1;
    p1 += M_PI;
    if (p1 > TWO_PI)
      p1 -= TWO_PI;
  }
  if (p1 < 0.)
    p1 += TWO_PI;

  amz = sdt - sdp; amx = spt - spp; amy = cpt - cpp;
  d2 = atan2(sqrt(amx*amx + amy*amy), amz);
  p2 = atan2(amy, -amx);
  if (d2 > M_PI_2) {
    d2 = M_PI - d2;
    p2 += M_PI;
    if (p2 > TWO_PI)
      p2 -= TWO_PI;
  }
  if (p2 < 0.)
    p2 += TWO_PI;

  NP1->dip = d1; NP1->str = p1;
  NP2->dip = d2; NP2->str = p2;

  im = 1;
  if (dp > dt) im = -1;
  NP1->rake = computed_rake2(NP2->str,NP2->dip,NP1->str,NP1->dip,im);
  NP2->rake = computed_rake2(NP1->str,NP1->dip,NP2->str,NP2->dip,im);
}
/********************************************************************/
double computed_rake2(double str1,double dip1,double str2,double dip2,double fault)

/*
  Compute rake in the second nodal plane when strike and dip
  for first and second nodal plane are given with a double
  characterizing the fault :
  +1. inverse fault
  -1. normal fault.
  Angles are in radians
*/

/* Genevieve Patau */

{
  double rake2, sinrake2;
  double sd, cd, ss, cs;

  sincos (str1 - str2, &ss, &cs);

  sd = sin(dip1);        cd = cos(dip2);
  if (fabs(dip2 - HALF_PI) < BC_EPSIL)
    sinrake2 = fault * cd;
  else
    sinrake2 = -fault * sd * cs / cd;

  rake2 = atan2(sinrake2, - fault * sd * ss);

  return(rake2);
}

#define MAX_SWEEPS 50

int	GMT_jacobi (double *a, int *n, int *m, double *d, double *v, double *b, double *z, int *nrots) {
  /*
   *
   * Find eigenvalues & eigenvectors of a square symmetric matrix by Jacobi's
   * method.  Given A, find V and D such that A = V * D * V-transpose, with
   * V an orthogonal matrix and D a diagonal matrix.  The eigenvalues of A
   * are on diag(D), and the j-th column of V is the eigenvector corresponding
   * to the j-th diagonal element of D.  Returns 0 if OK, -1 if it fails to
   * converge in MAX_SWEEPS.
   *
   * a is sent as a square symmetric matrix, of size n, and row dimension m.
   * Only the diagonal and super-diagonal elements of a will be used, so the
   * sub-diagonal elements could be used to preserve a, or could have been
   * destroyed by an earlier attempt to form the Cholesky decomposition of a.
   * On return, the super-diagonal elements are destroyed.  The diagonal and
   * sub-diagonal elements are unchanged.
   * d is returned as an n-vector containing the eigenvalues of a, sorted 
   * so that d[i] >= d[j] when i < j.  d = diag(D).
   * v is returned as an n by n matrix, V, with row dimension m, and the
   * columns of v are the eigenvectors corresponding to the values in d.
   * b is an n-vector of workspace, used to keep a copy of the diagonal
   * elements which is updated only after a full sweep.  
   * z is an n-vector of workspace, used to accumulate the updates to
   * the diagonal values of a during each sweep.  This reduces round-
   * off problems.
   * nrots is the number of rotations performed.  Bounds on round-off error
   * can be estimated from this if desired.
   *
   * Numerical Details:
   * The basic algorithms is in many textbooks.  The idea is to make an
   * infinite series (which turns out to be at quadratically convergent)
   * of steps, in each of which A_new = P-transpose * A_old * P, where P is
   * a plane-rotation matrix in the p,q plane, through an angle chosen to
   * zero A_new(p,q) and A_new(q,p).  The sum of the diagonal elements
   * of A is unchanged by these operations, but the sum of squares of 
   * diagonal elements of a is increased by 2 * |A_old(p,q)| at each step.
   * Although later steps make non-zero again the previously zeroed entries,
   * the sum of squares of diagonal elements increases with each rotation,
   * while the sum of squares of off-diagonals keeps decreasing, so that
   * eventually A_new is diagonal to machine precision.  This should 
   * happen in a few (3 to 7) sweeps.
   *
   * If only the eigenvalues are wanted then there are faster methods, but
   * if all eigenvalues and eigenvectors are needed, then this method is
   * only somewhat slower than the fastest method (Householder tri-
   * diagonalization followed by symmetric QR iterations), and this method
   * is numerically extremely stable.
   *
   * C G J Jacobi ("Ueber ein leichtes Vefahren, die in der Theorie der 
   * Saekularstoerungen vorkommenden Gelichungen numerisch aufzuloesen",
   * Crelle's Journal, v. 30, pp. 51--94, 1846) originally searched the
   * entire (half) matrix for the largest |A(p,q)| to select each step.
   * When the method was developed for machine computation (R T Gregory,
   * "Computing eigenvalues and eigenvectors of a symmetric matrix on
   * the ILLIAC", Math. Tab. and other Aids to Comp., v. 7, pp. 215--220,
   * 1953) it was done with a series of "sweeps" through the upper triangle,
   * visiting all p,q in turn.  Later, D A Pope and C Tompkins ("Maximizing
   * functions of rotations - experiments concerning speed of diagonalization
   * of symmetric matrices using Jacobi's method", J Assoc. Comput. Mach.
   * v. 4, pp. 459--466, 1957) introduced a variant that skips small
   * elements on the first few sweeps.  The algorithm here was given by 
   * Heinz Rutishauser (1918--1970) and published in Numer. Math. v. 9, 
   * pp 1--10, 1966, and in Linear Algebra (the Handbook for Automatic 
   * Computation, v. II), by James Hardy Wilkinson and C. Reinsch (Springer-
   * Verlag, 1971).  It also appears in Numerical Recipes.
   *
   * This algorithm takes care to avoid round-off error in several ways.
   * First, although there are four values of theta in (-pi, pi] that
   * would zero A(p,q), there is only one with magnitude <= pi/4.
   * This one is used.  This is most stable, and also has the effect
   * that, if A_old(p,p) >= A_old(q,q) then A_new(p,p) > A_new(q,q).
   * Two copies of the diagonal elements are maintained in d[] and b[].
   * d[] is updated immediately in each rotation, and each new rotation
   * is computed based on d[], so that each rotation gets the benefit
   * of the previous ones.  However, z[] is also used to accumulate
   * the sum of all the changes in the diagonal elements during one sweep, 
   * and z[] is used to update b[] after each sweep.  Then b is copied
   * to d.  In this way, at the end of each sweep, d is reset to avoid
   * accumulating round-off.
   *
   * This routine determines whether y is small compared to x by testing
   * if (fabs(y) + fabs(x) == fabs(x) ).  It is assumed that the
   * underflow which may occur here is nevertheless going to allow this
   * expression to be evaluated as TRUE or FALSE and execution to 
   * continue.  If the run environment doesn't allow this, the routine
   * won't work properly.
   *
   * programmer:	W. H. F. Smith, 7 June, 1991.
   * Revised:	PW: 12-MAR-1998 for GMT 3.1
   * Revision by WHF Smith, March 03, 2000, to speed up loop indexes.
   */
  int	p, q, pp, pq, mp1, pm, qm, nsweeps, j, jm, i, k;
  double	sum, threshold, g, h, t, theta, c, s, tau;


  /* Begin by initializing v, b, d, and z.  v = identity matrix,
     b = d = diag(a), and z = 0:  */
	
  memset ((void *)v, 0, (size_t)((*m)*(*n)*sizeof(double)) );
  memset ((void *)z, 0, (size_t)((*n)*sizeof(double)) );
	
  mp1 = (*m) + 1;
	
  for (p = 0, pp = 0; p < (*n); p++, pp+=mp1) {
    v[pp] = 1.0;
    b[p] = a[pp];
    d[p] = b[p];
  }

  /* End of initializations.  Set counters and begin:  */

  (*nrots) = 0;
  nsweeps = 0;

  while (nsweeps < MAX_SWEEPS) {
	
    /* Sum off-diagonal elements of upper triangle.  */
    sum = 0.0;
    for (q = 1, qm = (*m); q < (*n); q++, qm += (*m) ) {
      for (p = 0, pq = qm; p < q; p++, pq++) {
	sum += fabs(a[pq]);
      }
    }
		
    /* Exit this loop (converged) when sum == 0.0  */
    if (sum == 0.0) break;


    /* If (nsweeps < 3) do only bigger elements;  else all  */
    threshold =  (nsweeps < 3) ? 0.2 * sum / ( (*n) * (*n) ) : 0.0;

    /* Now sweep whole upper triangle doing Givens rotations:  */
		
    for (q = 1, qm = (*m); q < (*n); q++, qm += (*m) ) {
      for (p = 0, pm = 0, pq = qm; p < q; p++, pm += (*m), pq++) {
	/* In 3/2000 I swapped order of these loops,
	   to allow simple incrementing of pq  */
			
	if (a[pq] == 0.0) continue;	/* New 3/2000  */
			
	g = 100.0 * fabs(a[pq]);
				
	/* After four sweeps, if g is small relative
	   to a(p,p) and a(q,q), skip the 
	   rotation and set a(p,q) to zero.  */

	if ( (nsweeps > 3) && ( (fabs(d[p])+g) == fabs(d[p]) ) && ( (fabs(d[q])+g) == fabs(d[q]) ) ) {
	  a[pq] = 0.0;
	}
	else if (fabs(a[pq]) > threshold) {

	  h = d[q] - d[p];
					
	  if (h == 0.0) {
	    t = 1.0;	/* This if block is new 3/2000  */
	  }
	  else if ( (fabs(h)+g) ==  fabs(h) ) {
	    t = a[pq] / h;
	  }
	  else {
	    theta = 0.5 * h / a[pq];
	    t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta) );
	    if (theta < 0.0) t = -t;
	  }

	  c = 1.0 / sqrt(1.0 + t*t);
	  s = t * c;
	  tau = s / (1.0 + c);
					
	  h = t * a[pq];
	  z[p] -= h;
	  z[q] += h;
	  d[p] -= h;
	  d[q] += h;
	  a[pq] = 0.0;

	  for (j = 0; j < p; j++) {
	    g = a[j + pm];
	    h = a[j + qm];
	    a[j + pm] = g - s * (h + g * tau);
	    a[j + qm] = h + s * (g - h * tau);
	  }
	  for (j = p+1, jm = (*m)*(p+1); j < q; j++, jm += (*m) ) {
	    g = a[p + jm];
	    h = a[j + qm];
	    a[p + jm] = g - s * (h + g * tau);
	    a[j + qm] = h + s * (g - h * tau);
	  }
	  for (j = q+1, jm = (*m)*(q+1); j < (*n); j++, jm += (*m) ) {
	    g = a[p + jm];
	    h = a[q + jm];
	    a[p + jm] = g - s * (h + g * tau);
	    a[q + jm] = h + s * (g - h * tau);
	  }

	  for (j = 0; j < (*n); j++) {
	    g = v[j + pm];
	    h = v[j + qm];
	    v[j + pm] = g - s * (h + g * tau);
	    v[j + qm] = h + s * (g - h * tau);
	  }

	  (*nrots)++;
	}
      }
    }
		
    /* End of one sweep of the upper triangle.  */
		
    nsweeps++;

    for (p = 0; p < (*n); p++) {
      b[p] += z[p];	/* Update the b copy of diagonal  */
      d[p] = b[p];	/* Replace d with b to reduce round-off error  */
      z[p] = 0.0;	/* Clear z.  */
    }
  }

  /* Get here via break when converged, or when nsweeps == MAX_SWEEPS.
     Sort eigenvalues by insertion:  */

  for (i = 0; i < (*n)-1; i++) {
    k = i;
    g = d[i];
    for (j = i+1; j < (*n); j++) {  /* Find max location  */
      if (d[j] >= g) {
	k = j;
	g = d[j];
      }
    }
    if (k != i) {  /*  Need to swap value and vector  */
      d[k] = d[i];
      d[i] = g;
      p = i * (*m);
      q = k * (*m);
      for (j = 0; j < (*n); j++) {
	g = v[j + p];
	v[j + p] = v[j + q];
	v[j + q] = g;
      }
    }
  }

  /* Return 0 if converged; else print warning and return -1:  */

  if (nsweeps == MAX_SWEEPS) {
    fprintf (stderr, "GMT_jacobi:  Failed to converge in %i sweeps\n", nsweeps);
    return(-1);
  }
  return(0);
}


