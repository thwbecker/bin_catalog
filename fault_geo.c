#include "catalog.h"
/* 

   fault, geometry, and generic earthquake stuff routines
   
 */
void stridip(BC_CPREC n,BC_CPREC e,BC_CPREC u,BC_CPREC *strike,BC_CPREC *dip)
/* finds the strike and dip of a plane given its normal */
/* vector, output is in degrees north of east and then  */
/* uses a right hand rule for the dip of the plane */
{
  BC_CPREC x;
  if(u <0.) {
    n= -n;
    e= -e;
    u= -u;
  }
  *strike=atan2(e,n)*BC_PIF;
  *strike= *strike-90.;
  if(*strike < 0.)
    *strike+= 360.;
  if(*strike > 360.)
    *strike-= 360.;
  x=sqrt(n*n+e*e);   /* x is the horizontal magnitude */
  *dip=atan2(x,u)*BC_PIF;
}
void ranger(BC_CPREC *z)
/* makes z in 0 to 360 */
{
  while(*z >= 360)
    *z-= 360;
  while(*z < 0)
    *z += 360;
}

/* code from Andy Michael, as the SLICK package */

/* switcher routine 

   find alternate fault plane from given strike(dir, CW from N), dip
   and rake
*/
void find_alt_plane(BC_CPREC ddir1,BC_CPREC dip1,BC_CPREC rake1,
		    BC_CPREC *ddir2,BC_CPREC *dip2,BC_CPREC *rake2)
{
  BC_CPREC z,z2,z3,s1,s2,s3,sin_z,cos_z,sin_z2,cos_z2,sin_z3,cos_z3;
  BC_CPREC  n1,n2,h1,h2;
  
  z=ddir1/BC_PIF;
  if(dip1==90)
    dip1=89.99999;
  z2=dip1/BC_PIF;
  z3=rake1/BC_PIF;

  sin_z  = sin(z); cos_z  = cos(z);
  sin_z2 = sin(z2);cos_z2 = cos(z2);
  sin_z3 = sin(z3);cos_z3 = cos(z3);
  
  /* slick vector in plane 1 */
  s1= -cos_z3*cos_z-sin_z3*sin_z*cos_z2;
  s2=  cos_z3*sin_z-sin_z3*cos_z*cos_z2;
  s3= sin_z3*sin_z2;

  n1=sin_z*sin_z2;  /* normal vector to plane 1 */
  n2=cos_z*sin_z2;

  h1= -s2; /* strike vector of plane 2 */
  h2= s1;
  /* note h3=0 always so we leave it out */
  stridip(s2,s1,s3,&z,&z2);
  z+= 90.;
  *ddir2=z;
  ranger(ddir2);
  *dip2=z2;
  z= h1*n1 + h2*n2;
  z/= sqrt(h1*h1 + h2*h2);
  z=acos(z);
  if(s3>=0)
    *rake2= z*BC_PIF;
  else
    *rake2= -z*BC_PIF;
}

/*
 
convert event in Aki convention (angles in deg) into moment, magnitude
is Ml

*/

void aki2mom(BC_CPREC strike, BC_CPREC dip, BC_CPREC rake, BC_CPREC mag, 
	     BC_CPREC *m, BC_CPREC *m0)
{
  BC_CPREC phi,delta,gamma,phi2,delta2,sin_2delta,cos_gamma,cos_2delta;
  BC_CPREC sin_phi,cos2_phi,sin2_phi,sin_delta,sin_gamma,cos_delta,
    cos_phi,cos_2phi,sin_2phi;
  // angles
  phi =   BC_DEG2RAD(strike);		        // strike
  delta = BC_DEG2RAD(dip);			// dip
  gamma = BC_DEG2RAD(rake);			// rake

  /* use moment */
  *m0 = mag2mom(mag);

  //
  // trigonometry
  phi2 = phi * 2.0;
  delta2 = delta * 2.0;

  sin_phi   = sin(phi);  
  cos_phi   = cos(phi);  
  sin_2phi  = sin(phi2);  
  cos_2phi  = cos(phi2);
  sin2_phi = sin_phi*sin_phi;
  cos2_phi = cos_phi*cos_phi;
  sin_delta = sin(delta);
  cos_delta = cos(delta);
  sin_2delta= sin(delta2);
  cos_2delta= cos(delta2);
  sin_gamma = sin(gamma);
  cos_gamma = cos(gamma);
  //
  // moment components, x, y, z, = East, North, Up = P, -T, R
  //
  m[BC_RR] =  (sin_2delta * sin_gamma); /*                                                     UU */
  m[BC_RT] = -(cos_delta  * cos_gamma * cos_phi  +     cos_2delta * sin_gamma * sin_phi); /*  -NU */
  m[BC_RP] =  (cos_delta  * cos_gamma * sin_phi  -     cos_2delta * sin_gamma * cos_phi); /*   EU */
  m[BC_TT] = -(sin_delta  * cos_gamma * sin_2phi +     sin_2delta * sin_gamma * sin2_phi); /*  NN */
  m[BC_TP] = -(sin_delta  * cos_gamma * cos_2phi + 0.5*sin_2delta * sin_gamma * sin_2phi); /* -EN */
  m[BC_PP] =  (sin_delta  * cos_gamma * sin_2phi -     sin_2delta * sin_gamma * cos2_phi); /*  EE */
}
/*  */
BC_CPREC mean_hor_strain(BC_CPREC *a)
{
  BC_CPREC s11,s12,s22,r,x1,x2,fms,sms;
  s11 =  a[BC_PP];
  s12 = -a[BC_TP];
  s22 =  a[BC_TT];
  x1 = (s11 + s22)/2.0;
  x2 = (s11 - s22)/2.0;
  r = x2 * x2 + s12 * s12;
  if(r > 0.0){
    r = sqrt(r);
    fms = x1 + r;
    sms = x1 - r;
  }else{
    fms = sms = x1;
  }
  return ((fms+sms)/2.0);
}
/* 

moment conversions
 from magnitude to moment 
*/
BC_CPREC mag2mom(BC_CPREC mag)
{				     /* Hanks and Kanamori (1979) */
  return pow(10.0,3./2.*mag + 9.1); /* M0 in Nm */
}
BC_CPREC mom2mag(BC_CPREC mom)
{
  return ((2./3.)*(log10(mom) - 9.1)); /* in Nm */
}

BC_CPREC mag2pot(BC_CPREC mag)
{				/* from Ben Zion? can't remember,
				   nonsense anyway */
				/* scalar potency */
  if(mag < 3.5){		// lower scaling
    return pow(10.0,1.45 * mag - 5.69);
  }else{ //
    return pow(10.0,1.08 * mag - 4.87);
  } 
}

/* input is a [6] tensor */
BC_CPREC scalar_mom(BC_CPREC *m)	/* compute scalar moment from tensor */
{
  return tensor6_norm(m) * M_SQRT1_2; /* norm/sqrt(2) */
}
