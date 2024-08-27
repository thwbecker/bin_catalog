#include "catalog.h"
/*

modified from STRESSINVERSE https://www.ig.cas.cz/en/stress-inverse/

Vavryčuk, V., 2014. Iterative joint inversion for stress and fault
orientations from focal mechanisms, Geophysical Journal International,
199, 69-77, doi: 10.1093/gji/ggu224.


%*************************************************************************%
%                                                                         %
%   function SLIP_DEVIATION                                               %
%                                                                         %
%   callculating the deviation between the theoretical and observed slip  %
%                                                                         %
%   input: fault normal n, slip direction u, stress tensor                %
%   output: slip_deviation_1 - fault is identified by strike and dip      %
%           slip_deviation_2 - fault is the second nodal plane            %
%                                                                         %
%*************************************************************************%



returns average dot product

*/
void calc_misfits_from_single_angle_set(BC_CPREC *stress,BC_CPREC *angles, int nquakes, BC_CPREC *sdev)
{
  BC_CPREC m_smat[3][3],lsdev[2];
  int iquake,iquake6;
  my6stress2m3x3(stress,m_smat);
  /* compute misfit */
  sdev[0]=sdev[1]=0.0;
  for(iquake=iquake6=0;iquake < nquakes;iquake++,iquake6+=6){
    slip_deviation_mmat_single(m_smat,(angles+iquake6),(lsdev),(lsdev+1));
    //fprintf(stderr,"%g %g\n",lsdev[0],lsdev[1]);
    sdev[0]+=lsdev[0];
    sdev[1]+=lsdev[1];
  }
  sdev[0]/=(BC_CPREC)nquakes;
  sdev[1]/=(BC_CPREC)nquakes;
}


/* 
   
   my stress vector format
   angles[3]: strike dip rake in radian
   output is in radian as well

 */

void slip_deviation_svec_single(BC_CPREC *svec,BC_CPREC *angles,
				BC_CPREC *slip_deviation_1,BC_CPREC *slip_deviation_2)
{
  BC_CPREC mmat[3][3];
  my6stress2m3x3(svec,mmat);	/* convert my stress vector to a 3x3
				   matrix in Michael format */
  slip_deviation_mmat_single(mmat,angles,slip_deviation_1,slip_deviation_2);
}
/* 
   
   stress tensor is [3][3] on input with all entries filled
   Michael format, ENU
   strike dip rake in radian

   returns dot product


 */

/* strike, dip, rake in radians */
void slip_deviation_mmat_single(BC_CPREC tau[3][3],BC_CPREC *angles,
				BC_CPREC *slip_deviation_1,BC_CPREC *slip_deviation_2)
{
  BC_CPREC u[3],n[3],u0[3],n0[3];
  BC_CPREC cos_rake,sin_rake,sin_dip,cos_dip,sin_strike,cos_strike;
  sincos(angles[0],&sin_strike,&cos_strike);  
  sincos(angles[1],&sin_dip,&cos_dip);
  sincos(angles[2],&sin_rake,&cos_rake);
  

  
  //--------------------------------------------------------------------------
  //  fault normals and slip directions
  //--------------------------------------------------------------------------
  u0[0] =  cos_rake*cos_strike + cos_dip*sin_rake*sin_strike;
  u0[1] =  cos_rake*sin_strike - cos_dip*sin_rake*cos_strike;
  u0[2] = -sin_rake*sin_dip;
   
  n0[0] = -sin_dip*sin_strike;
  n0[1] =  sin_dip*cos_strike;
  n0[2] = -cos_dip;

  //--------------------------------------------------------------------------
  // calculation of slip_deviation_1
  //--------------------------------------------------------------------------
  n[0] = n0[0]; n[1] = n0[1]; n[2] = n0[2];
  u[0] = u0[0]; u[1] = u0[1]; u[2] = u0[2];
  /* first plane */
  *slip_deviation_1 = slip_deviation_dotp(n,u,tau);
  
  //--------------------------------------------------------------------------
  // calculation of slip_deviation_2
  //--------------------------------------------------------------------------
  n[0] = u0[0]; n[1] = u0[1]; n[2] = u0[2];
  u[0] = n0[0]; u[1] = n0[1]; u[2] = n0[2];
  
  if (n[2]>0){
    n[0] = -n[0];
    u[0] = -u[0];
  }; // vertical component is always negative!
  *slip_deviation_2 = slip_deviation_dotp(n,u,tau);

}
/* actually compute the dot product */

BC_CPREC slip_deviation_dotp(BC_CPREC *n,BC_CPREC *u,BC_CPREC tau[3][3])
{
  BC_CPREC tau_normal,tau_normal_square,tau_shear_square,
    tau_shear,tau_total,tau_total_square;
  BC_CPREC traction[3],shear_traction[3],tmp[3];
#ifdef DEBUG
  BC_CPREC unity1,unity2;
#endif
  //--------------------------------------------------------------------------
  // shear and normal stresses 
  //--------------------------------------------------------------------------
  tau_normal =
    tau[0][0]*n[0]*n[0] + tau[0][1]*n[0]*n[1] + tau[0][2]*n[0]*n[2] +
    tau[1][0]*n[1]*n[0] + tau[1][1]*n[1]*n[1] + tau[1][2]*n[1]*n[2] +
    tau[2][0]*n[2]*n[0] + tau[2][1]*n[2]*n[1] + tau[2][2]*n[2]*n[2];
  
  tau_normal_square = tau_normal * tau_normal;
  
  tmp[0] = tau[0][0]*n[0] + tau[0][1]*n[1] + tau[0][2]*n[2];
  tmp[1] = tau[1][0]*n[0] + tau[1][1]*n[1] + tau[1][2]*n[2];
  tmp[2] = tau[2][0]*n[0] + tau[2][1]*n[1] + tau[2][2]*n[2];
  tau_total_square   = tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2];
  
  tau_shear_square   = tau_total_square - tau_normal_square;
  
  tau_shear  = sqrt(tau_shear_square);
  tau_total  = sqrt(tau_total_square);
  
  //--------------------------------------------------------------------------
  // projection of stress into the fault plane
  //--------------------------------------------------------------------------
  // traction
  traction[0] = tau[0][0]*n[0] + tau[0][1]*n[1] + tau[0][2]*n[2];
  traction[1] = tau[1][0]*n[0] + tau[1][1]*n[1] + tau[1][2]*n[2];
  traction[2] = tau[2][0]*n[0] + tau[2][1]*n[1] + tau[2][2]*n[2];
  
  // projection of the traction into the fault plane
  shear_traction[0] = (traction[0] - tau_normal*n[0])/tau_shear;
  shear_traction[1] = (traction[1] - tau_normal*n[1])/tau_shear;
  shear_traction[2] = (traction[2] - tau_normal*n[2])/tau_shear;

#ifdef DEBUG
  // checking whether the calculations are correct
  unity1 = sqrt( shear_traction[0]*shear_traction[0] + shear_traction[1]*shear_traction[1] + shear_traction[2]*shear_traction[2]);
  unity2 = sqrt( u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  if((fabs(unity1-1)>1e-8)||(fabs(unity2-1)>1e-8)){
    fprintf(stderr,"|s| %g |u| %g\n",unity1,unity2);
  }
#endif
  
  // deviation between the slip and stress direction in radian
  //return acos(shear_traction[0]*u[0] + shear_traction[1]*u[1] + shear_traction[2]*u[2]); /*  */
  /* dot product */

  /* slip the sign because of stress sign convention? */
  return -(shear_traction[0]*u[0] + shear_traction[1]*u[1] + shear_traction[2]*u[2]);
}
