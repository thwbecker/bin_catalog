#include "catalog.h"
/*

  stress_inversion_vavrycuk.c

  Alternative Vavrycuk-type iterative stress inversion that follows the
  *exact* strategy of the MATLAB STRESSINVERSE package (Vavrycuk, 2014,
  GJI 199, 69-77), as opposed to the swap-until-stable optimizer in
  stability_criterion.c / stress_inversion.c.

  Differences from optimize_angles_via_instability():

  1. fixed iteration count per friction (no swap-until-converged, no
     bail-with-best). Batch reselection: at each iteration the more
     unstable of the two original nodal planes is chosen for every
     event simultaneously, then the Michael tensor is re-solved.

  2. initial guess is the average of N_realizations randomized-plane
     Michael tensors (each max-abs-eigenvalue normalized), matching
     linear_stress_inversion_Michael.m + the averaging loop in
     stress_inversion.m.

  3. the friction scan carries tau across friction values (tau is not
     reset between frictions), and the optimum friction is the one with
     maximum mean instability, followed by a final pass. This mirrors
     stress_inversion.m line for line.

  4. the instability / plane-selection is a direct port of
     stability_criterion.m, with sigma1 = SMALLEST eigenvalue (ascending
     sort, sigma_vector_1 = vector(:,j(1))) and shape_ratio =
     (sigma1-sigma2)/(sigma1-sigma3). The Michael tensor is first
     converted from the spherical storage used here (RR,RT,RP,TT,TP,PP)
     into the MATLAB Cartesian frame, then eigendecomposed in that
     frame, so the MATLAB fault normal and the eigenvectors are
     consistent. This differs from stability_criterion_eig() here, which
     labels sigma1 = largest eigenvalue and computes the normal in the
     Michael (dip-azimuth) frame; the two are not interchangeable for
     the iteration, so the selection is reimplemented rather than reused.

  Verification (10 SoCal mechanisms in becker_subset_angles.dat):
  reproduces the MATLAB STRESSINVERSE result to four digits, at fixed
  mu = 0.6 and for the full friction scan. Tensor (max|eig|-normalized,
  RR,RT,RP,TT,TP,PP): -0.4013 0.3881 0.3984 0.3092 -0.6770 0.0921;
  R = 0.5893; sigma1 az/pl 131.5/42.9; friction_opt = 0.500;
  mean instability 0.99316 at mu = 0.6, 0.99489 at the optimum.

  Reused, unmodified, from the existing sources:
     solve_stress_michael_specified_plane()  (Michael leasq, zero-trace)
     calc_eigensystem_vec6()                 (symmetric eigensystem)
     max_ev_normalize_tens6()                (max|eig| normalization)
     swap_angles(), find_alt_plane(), assign_quake_angles(), ran2()

  (c) 2026, written to match V. Vavrycuk's MATLAB code; reuses
  T. Becker's bin_catalog primitives. See README / COPYRIGHT.

*/

/*
  shared instability kernel, MATLAB (Vavrycuk reference) convention.

  vavrycuk_eigen: convert the Michael tensor from the spherical storage
  used here (RR,RT,RP,TT,TP,PP) to the MATLAB Cartesian 3x3
  [[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]], packed for calc_eigensystem_vec6
  as [xx,xy,xz,yy,yz,zz] (verified relation: xx=TT, yy=PP, zz=RR,
  xy=-TP, xz=RT, yz=-RP), and eigendecompose ascending so sigma[0] is
  the smallest eigenvalue (sigma1 in the .m file). Returns the three
  eigenvectors (v1=min .. v3=max) and sf = 1 - 2R.

  Working in this frame is what makes the MATLAB fault normal and the
  eigenvectors mutually consistent. This is the single definition of
  instability used by both the iterative solver and the averaging
  routine, so the two never drift apart.
*/
void vavrycuk_eigen(const BC_CPREC *stress,
		  BC_CPREC *v1, BC_CPREC *v2, BC_CPREC *v3,
		  BC_CPREC *sf)
{
  BC_CPREC m6[6], sigma[3], svec[9];
  m6[0] =  stress[BC_TT];   /* xx */
  m6[1] = -stress[BC_TP];   /* xy */
  m6[2] =  stress[BC_RT];   /* xz */
  m6[3] =  stress[BC_PP];   /* yy */
  m6[4] = -stress[BC_RP];   /* yz */
  m6[5] =  stress[BC_RR];   /* zz */
  calc_eigensystem_vec6(m6, sigma, svec, BC_TRUE, BC_FALSE); /* ascending */
  v1[0] = svec[0]; v1[1] = svec[1]; v1[2] = svec[2]; /* sigma1 (min) */
  v2[0] = svec[3]; v2[1] = svec[4]; v2[2] = svec[5]; /* sigma2       */
  v3[0] = svec[6]; v3[1] = svec[7]; v3[2] = svec[8]; /* sigma3 (max) */
  *sf = 1.0 - 2.0 * (sigma[0] - sigma[1]) / (sigma[0] - sigma[2]);
}

/*
  Vavrycuk instability of the two nodal planes of one event (6 angles,
  radians), given the eigensystem from vavrycuk_eigen, sf = 1-2R, friction
  mu, and ff = mu + sqrt(1+mu*mu). Result in inst[0] (plane 1), inst[1]
  (plane 2). sigma1 = smallest eigenvalue, MATLAB Cartesian normal.
*/
void vavrycuk_plane_inst(const BC_CPREC *v1, const BC_CPREC *v2,
		       const BC_CPREC *v3, BC_CPREC sf,
		       BC_CPREC mu, BC_CPREC ff,
		       const BC_CPREC *a6, BC_CPREC *inst)
{
  BC_CPREC ss, cs, sd, cd, nx, ny, nz, p1, p2, p3, p1s, p2s, p3s, tn, ts, tmp;
  int p, off;
  for (p = 0; p < 2; p++) {            /* p=0 -> plane 1, p=1 -> plane 2 */
    off = p * 3;
    sincos(a6[off],     &ss, &cs);
    sincos(a6[off + 1], &sd, &cd);
    /* MATLAB Cartesian fault normal (strike used directly, no +pi/2) */
    nx = -sd * ss;
    ny =  sd * cs;
    nz = -cd;
    p1 = nx * v1[0] + ny * v1[1] + nz * v1[2];
    p2 = nx * v2[0] + ny * v2[1] + nz * v2[2];
    p3 = nx * v3[0] + ny * v3[1] + nz * v3[2];
    p1s = p1 * p1; p2s = p2 * p2; p3s = p3 * p3;
    tn  = p1s + sf * p2s - p3s;                    /* normalized normal */
    tmp = p1s + sf * sf * p2s + p3s - tn * tn;
    ts  = (tmp > 0.0) ? sqrt(tmp) : 0.0;           /* normalized shear  */
    inst[p] = (ts - mu * (tn - 1.0)) / ff;
  }
}

/*
  pick, for each event, the more unstable of its two original nodal
  planes under stress tensor "stress" and friction "mu", writing the
  selected planes (selected plane first in each 6-block) into sel.
  returns the mean of the per-event chosen instabilities in *mean_inst.
*/
void vavrycuk_select_planes(int n, BC_CPREC *angles, BC_CPREC mu,
			  BC_CPREC *stress, BC_CPREC *sel,
			  BC_CPREC *mean_inst)
{
  BC_CPREC v1[3], v2[3], v3[3], sf, ff, inst[2], acc = 0.0;
  int j, j6;

  vavrycuk_eigen(stress, v1, v2, v3, &sf);
  ff = mu + sqrt(1.0 + mu * mu);

  for (j = j6 = 0; j < n; j++, j6 += 6) {
    vavrycuk_plane_inst(v1, v2, v3, sf, mu, ff, (angles + j6), inst);
    memcpy(sel + j6, angles + j6, 6 * sizeof(BC_CPREC));
    if (inst[1] > inst[0]) {           /* keep the more unstable plane first */
      swap_angles(sel + j6);
      acc += inst[1];
    } else {
      acc += inst[0];
    }
  }
  *mean_inst = acc / (BC_CPREC)n;
}

/*
  mean Vavrycuk instability of a stress tensor over a set of events,
  MATLAB convention, consistent with stress_inversion_vavrycuk (replaces
  calc_average_instability, which used the opposite sigma convention).

    ainst[0] = mean over events of the MORE unstable plane (the fault
               the stress predicts), i.e. the meaningful quality metric
    ainst[1] = mean over events of the LESS unstable plane

  angles are the two nodal planes per event (order does not matter, the
  max is taken). weights are accepted for API symmetry but, as in the
  MATLAB reference, the mean is unweighted.
*/
void vavrycuk_average_instability(int n, BC_CPREC *angles, BC_CPREC *weights,
                                BC_CPREC mu, BC_CPREC *stress, BC_CPREC *ainst)
{
  BC_CPREC v1[3], v2[3], v3[3], sf, ff, inst[2], hi = 0.0, lo = 0.0;
  int j, j6;
  vavrycuk_eigen(stress, v1, v2, v3, &sf);
  ff = mu + sqrt(1.0 + mu * mu);
  for (j = j6 = 0; j < n; j++, j6 += 6) {
    vavrycuk_plane_inst(v1, v2, v3, sf, mu, ff, (angles + j6), inst);
    if (inst[1] > inst[0]) { hi += inst[1]; lo += inst[0]; }
    else                   { hi += inst[0]; lo += inst[1]; }
  }
  ainst[0] = hi / (BC_CPREC)n;
  ainst[1] = lo / (BC_CPREC)n;
}

/*
  MATLAB-style iterative joint stress / fault inversion.

  inputs:
    n, angles (6 per event, radians), weights
    fmin, fmax, finc  friction scan (set fmin == fmax for fixed mu)
    n_iter            iterations per friction   (MATLAB N_iterations,   6)
    n_real            realizations for init     (MATLAB N_realizations, 10)
    seed              RNG seed (ran2 convention; pass a negative long)
    norm_type         normalization of the RETURNED tensor only:
                      BC_STRESS_NORM_EV     -> max abs eigenvalue
                      BC_STRESS_NORM_TENSOR -> tensor (Frobenius) norm
                      (the iteration itself is scale invariant, so this
                       affects the returned scale only, not R, the
                       instability, the resolved planes, or convergence)
  outputs:
    stress    (6) tensor normalized per norm_type, R,theta,phi order
    shape_ratio   (sigma1-sigma2)/(sigma1-sigma3), ascending convention
    fopt          optimum friction
    minst         mean instability at the optimum
    sel_out   (6n) resolved planes, selected plane first (may be NULL)

  note: tau is NOT normalized between iterations. vavrycuk_select_planes
  and the Michael re-solve depend on tau only through scale invariant
  quantities (eigenvector directions, shape ratio) and the selected
  planes, so the per-iteration normalization that the MATLAB code does
  is redundant and is skipped here to avoid an extra eigensolve per
  step. The single normalization of the returned tensor is done once.
*/
void stress_inversion_vavrycuk(int n, BC_CPREC *angles, BC_CPREC *weights,
                             BC_CPREC fmin, BC_CPREC fmax, BC_CPREC finc,
                             int n_iter, int n_real, long int *seed,
                             BC_CPREC *stress, BC_CPREC *shape_ratio,
                             BC_CPREC *fopt, BC_CPREC *minst,
                             BC_CPREC *sel_out, int norm_type)
{
  BC_CPREC *sel, *rnd, raw[6], tau[6], acc[6], tnorm[6], sg[3], sv[9];
  BC_CPREC mean_inst, best_mean, fopt_l, mu;
  size_t asize = 6 * sizeof(BC_CPREC) * n;
  int it, r, j, j6, k;

  sel = (BC_CPREC *)malloc(asize);
  rnd = (BC_CPREC *)malloc(asize);
  if ((!sel) || (!rnd)) BC_MEMERROR("stress_inversion_vavrycuk");

  /* -------- initial guess: averaged randomized-plane Michael --------
     each realization is max|eig| normalized before averaging, matching
     linear_stress_inversion_Michael.m. this is only n_real calls and is
     not on the hot path; the converged result is initialization
     independent, so the choice here does not affect the answer. */
  for (k = 0; k < 6; k++) acc[k] = 0.0;
  for (r = 0; r < n_real; r++) {
    for (j = j6 = 0; j < n; j++, j6 += 6) {
      memcpy(rnd + j6, angles + j6, 6 * sizeof(BC_CPREC));
      if (BC_RGEN(seed) < 0.5)        /* random nodal plane to the front */
        swap_angles(rnd + j6);
    }
    solve_stress_michael_specified_plane(n, rnd, weights, raw, BC_FALSE);
    max_ev_normalize_tens6(raw, tnorm);
    for (k = 0; k < 6; k++)
      acc[k] += tnorm[k];
  }
  memcpy(tau, acc, 6 * sizeof(BC_CPREC)); /* tau0 (scale irrelevant to iter) */

  /* -------- friction scan, tau carried across (as MATLAB) ----------- */
  best_mean = -1e30; fopt_l = fmin;
  for (mu = fmin; mu <= fmax + 1e-9; mu += finc) {
    for (it = 0; it < n_iter; it++) {
      vavrycuk_select_planes(n, angles, mu, tau, sel, &mean_inst);
      solve_stress_michael_specified_plane(n, sel, weights, raw, BC_FALSE);
      memcpy(tau, raw, 6 * sizeof(BC_CPREC)); /* no norm: selection is scale invariant */
    }
    if (mean_inst > best_mean) {
      best_mean = mean_inst; fopt_l = mu;
    }
    if (finc <= 0.0)
      break;             /* guard against fmin==fmax, finc 0 */
  }

  /* -------- final pass at optimum friction (tau carries over) -------- */
  for (it = 0; it < n_iter; it++) {
    vavrycuk_select_planes(n, angles, fopt_l, tau, sel, &mean_inst);
    solve_stress_michael_specified_plane(n, sel, weights, raw, BC_FALSE);
    memcpy(tau, raw, 6 * sizeof(BC_CPREC));
  }
  /* resolve planes and mean instability consistent with the FINAL tensor
     (the loop above leaves mean_inst one solve behind tau) */
  vavrycuk_select_planes(n, angles, fopt_l, tau, sel, &mean_inst);

  /* -------- normalize the returned tensor once, per request -------- */
  if (norm_type == BC_STRESS_NORM_TENSOR)
    normalize_tens6(tau, tau);
  else
    max_ev_normalize_tens6(tau, tau);    /* BC_STRESS_NORM_EV (default) */

  /* -------- outputs -------- */
  memcpy(stress, tau, 6 * sizeof(BC_CPREC));
  calc_eigensystem_vec6(tau, sg, sv, BC_FALSE, BC_FALSE); /* ascending */
  *shape_ratio = (sg[0] - sg[1]) / (sg[0] - sg[2]);
  *fopt  = fopt_l;
  *minst = mean_inst;
  if (sel_out)
    memcpy(sel_out, sel, asize);

  free(sel); free(rnd);
}

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
   
   stress tensor is [3][3] on input with all entries filled
   Michael format, ENU
   strike dip rake in radian

   returns dot product

   
*/

/* strike, dip, rake in radians */
void slip_deviation_mmat_single(BC_CPREC tau[3][3],BC_CPREC *angles,
				BC_CPREC *slip_dotp_1,BC_CPREC *slip_dotp_2)
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
  // calculation of slip_dotp_1
  //--------------------------------------------------------------------------
  n[0] = n0[0]; n[1] = n0[1]; n[2] = n0[2];
  u[0] = u0[0]; u[1] = u0[1]; u[2] = u0[2];
  /* first plane */
  *slip_dotp_1 = slip_dev_dotp(n,u,tau);
  
  //--------------------------------------------------------------------------
  // calculation of slip_dotp_2
  //--------------------------------------------------------------------------
  n[0] = u0[0]; n[1] = u0[1]; n[2] = u0[2];
  u[0] = n0[0]; u[1] = n0[1]; u[2] = n0[2];
  
  if (n[2]>0){
    n[0] = -n[0];
    u[0] = -u[0];
  }; // vertical component is always negative!
  *slip_dotp_2 = slip_dev_dotp(n,u,tau);

}
/* 

   actually compute the dot product  - the close to unity, the better
   
*/

BC_CPREC slip_dev_dotp(BC_CPREC *n,BC_CPREC *u,BC_CPREC tau[3][3])
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
  
  // dotp between the slip and stress direction in radian
  //return acos(shear_traction[0]*u[0] + shear_traction[1]*u[1] + shear_traction[2]*u[2]); /*  */
  /* dot product */


  return (shear_traction[0]*u[0] + shear_traction[1]*u[1] + shear_traction[2]*u[2]);
}
