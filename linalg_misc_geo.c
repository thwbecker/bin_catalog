#include "catalog.h"

void ranger(BC_CPREC *z)
/* makes z in 0 to 360 */
{
  while(*z >= 360)
    *z-= 360;
  while(*z < 0)
    *z += 360;
}

BC_CPREC max_x_from_int_vector(BC_CPREC *x, int *y, int n)
{
  int i;
  int max = -1e8;
  BC_CPREC xmax;
  xmax = y[0];
  for(i=0;i<n;i++){
    if(y[i] > max){
      max = y[i];
      xmax = x[i];
    }
  }
  return xmax;
}
/* tin, tou[6] symmetric storage  */
void rotate_vec6(BC_CPREC *tin, BC_CPREC *tout, BC_CPREC alpha, BC_CPREC beta, BC_CPREC gamma)
{
  BC_CPREC r[3][3],xin[3][3],xout[3][3];
  get_gen_rot(r,alpha,beta,gamma); /* get rotation matrix */
  sixsymtomat(tin,xin);		   /* convert to 3x3 */
  rotate_mat(xin,xout,r);	   /* rotate 3x3 */
  mattosixsym(xout,tout);	   /* convert to [6] */
}
/* convert symmetric [6] storage to [3][3] matrix */
void sixsymtomat(BC_CPREC *in6,BC_CPREC out[3][3])
{
  out[BC_R][BC_R] = in6[BC_RR];
  out[BC_R][BC_THETA]= out[BC_THETA][BC_R] = in6[BC_RT];
  out[BC_R][BC_PHI] = out[BC_PHI][BC_R] = in6[BC_RP];
  out[BC_THETA][BC_THETA] = in6[BC_TT];
  out[BC_THETA][BC_PHI] = out[BC_PHI][BC_THETA] = in6[BC_TP];
  out[BC_PHI][BC_PHI] = in6[BC_PP];
}
void mattosixsym(BC_CPREC in[3][3], BC_CPREC *out6)
{
  out6[BC_RR] = in[BC_R][BC_R] ;
  out6[BC_RT] = in[BC_R][BC_THETA];
  out6[BC_RP] = in[BC_R][BC_PHI];
  out6[BC_TT] = in[BC_THETA][BC_THETA];
  out6[BC_TP] = in[BC_THETA][BC_PHI];
  out6[BC_PP] = in[BC_PHI][BC_PHI];
}

/*
  
  obtain a general rotation matrix with angles alpha, beta, and gamma
  (given in degrees) as defined in Dahlen and Tromp, p. 921
  
*/
void get_gen_rot(BC_CPREC r[3][3],BC_CPREC alpha,
		 BC_CPREC beta, BC_CPREC gamma)
{
  BC_CPREC ralpha,sin_alpha,cos_alpha;
  BC_CPREC rbeta,sin_beta,cos_beta;
  BC_CPREC rgamma,sin_gamma,cos_gamma;

  ralpha=BC_DEG2RAD(alpha);
  sin_alpha = sin(ralpha);
  cos_alpha = cos(ralpha);

  rbeta=BC_DEG2RAD(beta);
  sin_beta = sin(rbeta);
  cos_beta = cos(rbeta);

  rgamma=BC_DEG2RAD(gamma);
  sin_gamma = sin(rgamma);
  cos_gamma = cos(rgamma);
  
  r[BC_R][BC_R] = cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma; 
  r[BC_R][BC_THETA] = sin_alpha*cos_beta*cos_gamma + cos_alpha*sin_gamma;
  r[BC_R][BC_PHI] = -sin_beta*cos_gamma;
  r[BC_THETA][BC_R] = -cos_alpha*cos_beta*sin_gamma -sin_alpha*cos_gamma;
  r[BC_THETA][BC_THETA] = -sin_alpha*cos_beta*sin_gamma +cos_alpha*cos_gamma;
  r[BC_THETA][BC_PHI] = sin_beta*sin_gamma;
  r[BC_PHI][BC_R] = cos_alpha*sin_beta;
  r[BC_PHI][BC_THETA] = sin_alpha*sin_beta;
  r[BC_PHI][BC_PHI] = cos_beta;
}

/* 
   rotate a 3x3 tensor using a general rotation matrix r 
   xout = r . xin . r^T
*/
void rotate_mat(BC_CPREC xin[3][3],BC_CPREC xout[3][3],
		BC_CPREC r[3][3])
{
  BC_CPREC tmp[3][3];
  int i,j,k;
  // calculate xin . r^T
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      tmp[i][j] = 0.0;
      for(k=0;k<3;k++)
	tmp[i][j] += xin[i][k] * r[j][k];
    }
  // calculate r . tmp
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      xout[i][j] = 0.0;
      for(k=0;k<3;k++)
	xout[i][j] += r[i][k] * tmp[k][j];
    }
}
void remove_trace(BC_CPREC *a)
{
  BC_CPREC t3;
  t3 = (a[BC_RR]  + a[BC_TT] + a[BC_PP])/3.0;
  a[BC_RR] -= t3;
  a[BC_PP] -= t3;
  a[BC_TT] -= t3;
}
BC_CPREC distance_cart(BC_CPREC x1,BC_CPREC y1,BC_CPREC x2, BC_CPREC y2)
{
  BC_CPREC tmp1,tmp2;
  tmp1 = x1-x2;
  tmp1*=tmp1;
  tmp2 = y1-y2;
  tmp2*=tmp2;
  return sqrt(tmp1+tmp2);
}


/* 

   compute distance on sphere, input in degree, return in km 

   cosines for lat1 and lat2 are precomputed

*/
BC_CPREC distance_geo(BC_CPREC lon1,BC_CPREC lat1,
		   BC_CPREC lon2,BC_CPREC lat2,
		   BC_CPREC coslat1, BC_CPREC coslat2)
{
  BC_CPREC tmp1,tmp2,tmp3;

  lon1/=BC_PIF;
  lat1/=BC_PIF;
  lon2/=BC_PIF;
  lat2/=BC_PIF;
  tmp1 = sin((lat1 - lat2)/2.0);
  tmp1 = tmp1 * tmp1;
  tmp2 = sin((lon1 - lon2)/2.0);
  tmp2 = tmp2 * tmp2;
  tmp2 *= coslat1;
  tmp2 *= coslat2;
  tmp3 = sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3)*BC_RADIUS;
}
/* compute the distance between quake i of catalog one and quake j (or
   bin j) of catalog two, depending on cartesian coordinates or not */

BC_CPREC distance(struct cat *c1,struct cat *c2,int i,int j)
{
  if(c1->is_xy){
    return distance_cart(c1->quake[i].lon,c1->quake[i].lat,
			 c2->quake[j].lon,c2->quake[j].lat);
    
  }else{
    return distance_geo(c1->quake[i].lon,c1->quake[i].lat,
			c2->quake[j].lon,c2->quake[j].lat,
			c1->quake[i].coslat,c2->quake[j].coslat);
    
  }
}
void get_index_vector(int **ind, int n, int random,long *seed)
{
  int i;
  *ind = (int *)realloc(*ind,n * sizeof(int));
  if(!ind){BC_MEMERROR("index");}
  if(random > 0 ){
    for(i=0;i<n;i++)
      *(*ind+i) = (int)(BC_RGEN(seed) * n);
  }else{
    for(i=0;i < n;i++)
      *(*ind+i) = i;
  }
}
void normalize_tens6(BC_CPREC *m6)	/* normalize a tensor given in 0...5
				   format */
{
  BC_CPREC norm;
  int i;
  norm = tensor6_norm(m6);
  for(i=0;i<6;i++)
    m6[i] /= norm;
  /* test */
  //fprintf(stderr,"old norm: %12.5e new norm: %12.5e\n",norm,tensor6_norm(m6));
}

/* tensor in [6] notation */
BC_CPREC tensor6_norm(BC_CPREC *m6)	/* plug in upper triangle */
{
  BC_CPREC t[3][3], ret;
  int i,j;
  tens6to3by3(m6,t);		/* convert to 3x3 */
  ret = 0;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      ret += t[i][j] * t[i][j];
  return sqrt(ret);
}

void tens6to3by3(BC_CPREC *m, BC_CPREC t[3][3])
{
  t[BC_R][BC_R] = m[BC_RR];
  t[BC_THETA][BC_R] = t[BC_R][BC_THETA] = m[BC_RT];
  t[BC_THETA][BC_THETA] = m[BC_TT];
  t[BC_PHI][BC_THETA] = t[BC_THETA][BC_PHI] = m[BC_TP];
  t[BC_PHI][BC_PHI] = m[BC_PP];
  t[BC_R][BC_PHI] = t[BC_PHI][BC_R] = m[BC_RP];
}
/* 
   quick and dirty way to compute standard deviationq 

*/
BC_CPREC std_quick(int n, BC_CPREC sum, BC_CPREC sum2)
{
  double nf;
  nf = (double)n;
  return sqrt((nf * (double)sum2 - (double)sum*(double)sum) / (nf*(nf-1)));
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 5e-15
#define RNMX (1.0-EPS)
BC_CPREC ran2(long int *idum)
{
  int j;
  long k;
  static int ntabp7 = NTAB + 7;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  BC_CPREC temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=ntabp7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) 
    *idum += IM1;
  k=idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) 
    idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) 
    iy += IMM1;
  if ((temp=AM*iy) > RNMX) 
    return RNMX;
  else 
    return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
