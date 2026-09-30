#include "catalog.h"
/*

  stress_inversion_vavrycuk.c

  Alternative Vavrycuk-type iterative stress inversion that follows the
  strategy of the MATLAB STRESSINVERSE package (Vavrycuk, 2014,
  GJI 199, 69-77), as opposed to the swap-until-stable optimizer in
  stability_criterion.c / stress_inversion.c.

  Differences from optimize_angles_via_instability():

  1. at most a fixed iteration count per friction (no swap-until-converged, no
     bail-with-best; the iterations stop early only when the plane
     selection repeats, which cannot change the result). Batch reselection: at each iteration the more
     unstable of the two original nodal planes is chosen for every
     event simultaneously, then the Michael tensor is re-solved.

  2. initial guess is the average of N_realizations randomized-plane
     Michael tensors (each max-abs-eigenvalue normalized), matching
     linear_stress_inversion_Michael.m + the averaging loop in
     stress_inversion.m. Different from the MATLAB code, the scan is
     also run from each individual realization, and the result with the
     largest mean instability is kept (see stress_inversion_vavrycuk).

  3. the friction scan carries tau across friction values (tau is not
     reset between frictions), and the optimum friction is the one with
     maximum mean instability, followed by a final pass at the optimum
     starting from the tensor found there.

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
/*
  eigenvalues (ascending) and unit eigenvectors of a symmetric 3x3
  matrix by cyclic Jacobi rotations. a is overwritten. vec[k][i] is
  component i of the eigenvector of val[k]. converges to machine
  precision in a few sweeps and is several times faster than the
  general EISPACK route for this size, which matters since the
  Vavrycuk iteration needs one eigensystem per plane selection.
*/
static void vavrycuk_eig3(BC_CPREC a[3][3], BC_CPREC val[3], BC_CPREC vec[3][3])
{
  BC_CPREC v[3][3] = {{1,0,0},{0,1,0},{0,0,1}}, off, theta, t, c, sn, tau, tmp, d[3];
  int sweep, p, q, r, i, j, k;
  static const int pq[3][2] = {{0,1},{0,2},{1,2}};
  for (sweep = 0; sweep < 50; sweep++) {
    off = a[0][1]*a[0][1] + a[0][2]*a[0][2] + a[1][2]*a[1][2];
    if (off < 1e-30 * (a[0][0]*a[0][0] + a[1][1]*a[1][1] + a[2][2]*a[2][2]) || off == 0.0)
      break;
    for (k = 0; k < 3; k++) {
      p = pq[k][0]; q = pq[k][1];
      if (a[p][q] == 0.0) continue;
      theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
      t = ((theta >= 0) ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta*theta + 1.0));
      c = 1.0 / sqrt(t*t + 1.0); sn = t * c; tau = sn / (1.0 + c);
      tmp = a[p][q];
      a[p][p] -= t * tmp;
      a[q][q] += t * tmp;
      a[p][q] = a[q][p] = 0.0;
      for (r = 0; r < 3; r++) {
        if ((r == p) || (r == q)) continue;
        BC_CPREC arp = a[r][p], arq = a[r][q];
        a[r][p] = a[p][r] = arp - sn * (arq + tau * arp);
        a[r][q] = a[q][r] = arq + sn * (arp - tau * arq);
      }
      for (r = 0; r < 3; r++) {
        BC_CPREC vrp = v[r][p], vrq = v[r][q];
        v[r][p] = vrp - sn * (vrq + tau * vrp);
        v[r][q] = vrq + sn * (vrp - tau * vrq);
      }
    }
  }
  for (i = 0; i < 3; i++) d[i] = a[i][i];
  /* sort ascending */
  int idx[3] = {0, 1, 2};
  for (i = 0; i < 2; i++)
    for (j = i + 1; j < 3; j++)
      if (d[idx[j]] < d[idx[i]]) { k = idx[i]; idx[i] = idx[j]; idx[j] = k; }
  for (k = 0; k < 3; k++) {
    val[k] = d[idx[k]];
    for (i = 0; i < 3; i++) vec[k][i] = v[i][idx[k]];
  }
}

void vavrycuk_eigen(const BC_CPREC *stress,
		  BC_CPREC *v1, BC_CPREC *v2, BC_CPREC *v3,
		  BC_CPREC *sf)
{
  BC_CPREC a[3][3], sigma[3], vec[3][3];
  int i;
  /* MATLAB Cartesian frame: xx=TT, xy=-TP, xz=RT, yy=PP, yz=-RP, zz=RR */
  a[0][0] =  stress[BC_TT];
  a[0][1] = a[1][0] = -stress[BC_TP];
  a[0][2] = a[2][0] =  stress[BC_RT];
  a[1][1] =  stress[BC_PP];
  a[1][2] = a[2][1] = -stress[BC_RP];
  a[2][2] =  stress[BC_RR];
  vavrycuk_eig3(a, sigma, vec); /* ascending */
  for (i = 0; i < 3; i++) {
    v1[i] = vec[0][i];            /* sigma1 (min) */
    v2[i] = vec[1][i];            /* sigma2       */
    v3[i] = vec[2][i];            /* sigma3 (max) */
  }
  *sf = 1.0 - 2.0 * (sigma[0] - sigma[1]) / (sigma[0] - sigma[2]);
}

/* for testing: vavrycuk_eigen with the general EISPACK based routine */
void vavrycuk_eigen_eispack(const BC_CPREC *stress,
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
  normals of both nodal planes of each event in the MATLAB Cartesian
  frame used by vavrycuk_plane_inst: nrm[6*j + 3*p + k] for event j,
  plane p, component k
*/
static void vavrycuk_normals(int n, const BC_CPREC *angles, BC_CPREC *nrm)
{
  BC_CPREC ss, cs, sd, cd;
  int j, p, off;
  for (j = 0; j < n; j++) {
    for (p = 0; p < 2; p++) {
      off = 6 * j + 3 * p;
      sincos(angles[off],     &ss, &cs);
      sincos(angles[off + 1], &sd, &cd);
      nrm[off]     = -sd * ss;
      nrm[off + 1] =  sd * cs;
      nrm[off + 2] = -cd;
    }
  }
}

/* instability of one plane from its precomputed normal, as vavrycuk_plane_inst */
static inline BC_CPREC vavrycuk_inst_from_normal(const BC_CPREC *nv,
						 const BC_CPREC *v1, const BC_CPREC *v2,
						 const BC_CPREC *v3, BC_CPREC sf,
						 BC_CPREC mu, BC_CPREC ff)
{
  BC_CPREC p1, p2, p3, p1s, p2s, p3s, tn, ts, tmp;
  p1 = nv[0] * v1[0] + nv[1] * v1[1] + nv[2] * v1[2];
  p2 = nv[0] * v2[0] + nv[1] * v2[1] + nv[2] * v2[2];
  p3 = nv[0] * v3[0] + nv[1] * v3[1] + nv[2] * v3[2];
  p1s = p1 * p1; p2s = p2 * p2; p3s = p3 * p3;
  tn  = p1s + sf * p2s - p3s;
  tmp = p1s + sf * sf * p2s + p3s - tn * tn;
  ts  = (tmp > 0.0) ? sqrt(tmp) : 0.0;
  return (ts - mu * (tn - 1.0)) / ff;
}

/*
  plane selection with precomputed normals: out[j] = 1 if the second
  plane of event j is more unstable under stress/mu, else 0. returns the
  mean instability of the chosen planes, and in *changed whether out
  differs from ref
*/
static BC_CPREC vavrycuk_choose(int n, const BC_CPREC *nrm, BC_CPREC mu,
				const BC_CPREC *stress, const unsigned char *ref,
				unsigned char *out, BC_BOOLEAN *changed)
{
  BC_CPREC v1[3], v2[3], v3[3], sf, ff, i0, i1, acc = 0.0;
  unsigned char c;
  int j;
  vavrycuk_eigen(stress, v1, v2, v3, &sf);
  ff = mu + sqrt(1.0 + mu * mu);
  *changed = BC_FALSE;
  for (j = 0; j < n; j++) {
    i0 = vavrycuk_inst_from_normal(nrm + 6 * j,     v1, v2, v3, sf, mu, ff);
    i1 = vavrycuk_inst_from_normal(nrm + 6 * j + 3, v1, v2, v3, sf, mu, ff);
    c = (i1 > i0) ? 1 : 0;
    acc += (c) ? i1 : i0;
    out[j] = c;
    if (c != ref[j]) *changed = BC_TRUE;
  }
  return acc / (BC_CPREC)n;
}

/* Michael solution for the chosen planes from precomputed normal equations */
static void vavrycuk_solve_choice(int n, const BC_CPREC *ne,
				  const unsigned char *choice, BC_CPREC *stress)
{
  BC_CPREC sum[BC_MICHAEL_NNE];
  const BC_CPREC *nep;
  int j, k;
  for (k = 0; k < BC_MICHAEL_NNE; k++) sum[k] = 0.0;
  for (j = 0; j < n; j++) {
    nep = ne + (2 * j + choice[j]) * BC_MICHAEL_NNE;
    for (k = 0; k < BC_MICHAEL_NNE; k++) sum[k] += nep[k];
  }
  michael_normal_eq_solve(sum, stress);
}

/*
  up to n_iter iterations of plane selection and Michael re-solve at
  friction mu, starting from tau. choice holds the planes that tau was
  solved from (2 = none, for the initial guess) and is updated with tau.
  the loop stops early when the selection repeats, since the re-solve
  would then reproduce tau exactly. returns the mean instability of the
  selection under the returned tau, which is left in sel (scratch array
  of length n).
*/
static BC_CPREC vavrycuk_iterate(int n, const BC_CPREC *nrm, const BC_CPREC *ne,
				 BC_CPREC mu, int n_iter, BC_CPREC *tau,
				 unsigned char *choice, unsigned char *sel)
{
  BC_BOOLEAN changed;
  BC_CPREC mean_inst;
  int it;
  for (it = 0; it < n_iter; it++) {
    mean_inst = vavrycuk_choose(n, nrm, mu, tau, choice, sel, &changed);
    if (!changed)             /* fixed point: tau already solves this selection */
      return mean_inst;
    memcpy(choice, sel, n);
    vavrycuk_solve_choice(n, ne, choice, tau);
  }
  return vavrycuk_choose(n, nrm, mu, tau, choice, sel, &changed);
}

/*
  friction scan from one starting tensor: tau is carried from one
  friction to the next with up to n_iter iterations at each, the
  optimum is the friction with the largest mean instability, followed
  by a final pass at the optimum starting from the tensor found there.
  on return, tau is the final tensor, sel the plane selection under it,
  *fopt the optimum friction. returns the mean instability. choice,
  cbest are scratch arrays of length n.
*/
static BC_CPREC vavrycuk_scan(int n, const BC_CPREC *nrm, const BC_CPREC *ne,
			      BC_CPREC fmin, BC_CPREC finc, int nfric, int n_iter,
			      BC_CPREC *tau, BC_CPREC *fopt,
			      unsigned char *choice, unsigned char *cbest,
			      unsigned char *sel)
{
  BC_CPREC mean_inst, best_mean = -1e30, tbest[6], mu;
  int ifric;
  memset(choice, 2, n);       /* tau is not the solution of any selection */
  memcpy(tbest, tau, 6 * sizeof(BC_CPREC));
  memcpy(cbest, choice, n);
  *fopt = fmin;
  for (ifric = 0; ifric < nfric; ifric++) {
    mu = fmin + ifric * finc;
    mean_inst = vavrycuk_iterate(n, nrm, ne, mu, n_iter, tau, choice, sel);
    if (mean_inst > best_mean) {
      best_mean = mean_inst; *fopt = mu;
      memcpy(tbest, tau, 6 * sizeof(BC_CPREC));
      memcpy(cbest, choice, n);
    }
  }
  memcpy(tau, tbest, 6 * sizeof(BC_CPREC));
  memcpy(choice, cbest, n);
  return vavrycuk_iterate(n, nrm, ne, *fopt, n_iter, tau, choice, sel);
}

/*
  alternative candidate (the approach of the previous version of this
  code): converge a reference tensor at the middle of the friction
  range, pick the friction that maximizes the mean instability under
  that fixed tensor, then iterate at that friction. same interface and
  outputs as vavrycuk_scan.
*/
static BC_CPREC vavrycuk_scan_ref(int n, const BC_CPREC *nrm, const BC_CPREC *ne,
				  BC_CPREC fmin, BC_CPREC finc, int nfric, int n_iter,
				  BC_CPREC *tau, BC_CPREC *fopt,
				  unsigned char *choice, unsigned char *sel)
{
  BC_CPREC v1[3], v2[3], v3[3], sf, ff, mu, i0, i1, acc, best = -1e30;
  int ifric, j;
  memset(choice, 2, n);
  vavrycuk_iterate(n, nrm, ne, fmin + 0.5 * (nfric - 1) * finc, n_iter, tau, choice, sel);
  vavrycuk_eigen(tau, v1, v2, v3, &sf);
  *fopt = fmin;
  for (ifric = 0; ifric < nfric; ifric++) {
    mu = fmin + ifric * finc;
    ff = mu + sqrt(1.0 + mu * mu);
    for (j = 0, acc = 0.0; j < n; j++) {
      i0 = vavrycuk_inst_from_normal(nrm + 6 * j,     v1, v2, v3, sf, mu, ff);
      i1 = vavrycuk_inst_from_normal(nrm + 6 * j + 3, v1, v2, v3, sf, mu, ff);
      acc += (i1 > i0) ? i1 : i0;
    }
    acc /= (BC_CPREC)n;
    if (acc > best) { best = acc; *fopt = mu; }
  }
  return vavrycuk_iterate(n, nrm, ne, *fopt, n_iter, tau, choice, sel);
}

/*
  MATLAB-style iterative joint stress / fault inversion.

  inputs:
    n, angles (6 per event, radians), weights
    fmin, fmax, finc  friction scan (set fmin == fmax for fixed mu)
    n_iter            iterations per friction   (MATLAB N_iterations,   6)
    n_real            random plane realizations for the starting tensors
                      (MATLAB N_realizations, 10)
    seed              RNG seed (ran2 convention; pass a negative long)
    norm_type         normalization of the RETURNED tensor only:
                      BC_STRESS_NORM_EV     -> max abs eigenvalue
                      BC_STRESS_NORM_TENSOR -> tensor (Frobenius) norm
  outputs:
    stress    (6) tensor normalized per norm_type, R,theta,phi order
    shape_ratio   (sigma1-sigma2)/(sigma1-sigma3), ascending convention
    fopt          optimum friction
    minst         mean instability at the optimum
    sel_out   (6n) resolved planes, selected plane first (may be NULL)

  starting tensors: as in the MATLAB code, the average of n_real
  randomized-plane Michael tensors (each max|eig| normalized). the
  iteration can end in different local maxima of the mean instability
  depending on the start, mainly for small numbers of events. if
  n_real > 1, the scan (see vavrycuk_scan) is therefore also run from
  each of the n_real individual realizations, and, from the averaged
  start, the reference-tensor variant of the previous version of this
  code (vavrycuk_scan_ref) is added as a candidate. the result with the
  largest mean instability is returned (the averaged start with the
  full scan wins ties). the result is thus never worse, in terms of
  mean instability, than either the single averaged start or the
  previous version.

  implementation notes: normal equations and fault normals of both
  planes are computed once per call, so each iteration costs O(n)
  without trigonometry. iterations at a friction stop when the plane
  selection repeats (the result would not change). tau is not
  normalized between iterations; the selection depends on tau only
  through scale invariant quantities.
*/
void stress_inversion_vavrycuk(int n, BC_CPREC *angles, BC_CPREC *weights,
                             BC_CPREC fmin, BC_CPREC fmax, BC_CPREC finc,
                             int n_iter, int n_real, long int *seed,
                             BC_CPREC *stress, BC_CPREC *shape_ratio,
                             BC_CPREC *fopt, BC_CPREC *minst,
                             BC_CPREC *sel_out, int norm_type)
{
  BC_CPREC *ne, *nrm, *starts, raw[6], tau[6], sg[3], sv[9];
  BC_CPREC best_tau[6], mean_inst, best_mean = -1e30, fo, best_fo = fmin;
  unsigned char *choice, *cbest, *sel, *best_sel;
  int r, j, j6, k, nfric, nstart, is;

  if (n_real < 1) n_real = 1;
  nstart = (n_real > 1) ? (n_real + 1) : 1;
  ne       = (BC_CPREC *)malloc(sizeof(BC_CPREC) * BC_MICHAEL_NNE * 2 * n);
  nrm      = (BC_CPREC *)malloc(sizeof(BC_CPREC) * 6 * n);
  starts   = (BC_CPREC *)calloc(6 * (n_real + 1), sizeof(BC_CPREC));
  choice   = (unsigned char *)malloc(n);
  cbest    = (unsigned char *)malloc(n);
  sel      = (unsigned char *)malloc(n);
  best_sel = (unsigned char *)malloc(n);
  if ((!ne) || (!nrm) || (!starts) || (!choice) || (!cbest) || (!sel) || (!best_sel))
    BC_MEMERROR("stress_inversion_vavrycuk");
  michael_setup_normal_eq(n, angles, weights, ne);
  vavrycuk_normals(n, angles, nrm);

  /* -------- starting tensors: starts[0] average, starts[1..] realizations
     same random sequence as before (plane 2 used if the draw is < 0.5) */
  for (r = 0; r < n_real; r++) {
    for (j = 0; j < n; j++)
      choice[j] = (BC_RGEN(seed) < 0.5) ? 1 : 0;
    vavrycuk_solve_choice(n, ne, choice, raw);
    max_ev_normalize_tens6(raw, starts + 6 * (r + 1));
    for (k = 0; k < 6; k++)
      starts[k] += starts[6 * (r + 1) + k];
  }
  if ((finc > 0.0) && (fmax > fmin + 1e-9))
    nfric = (int)((fmax - fmin) / finc + 1e-9) + 1;
  else
    nfric = 1;
  for (is = 0; is < nstart + ((nstart > 1) ? 1 : 0); is++) {
    if (is < nstart) {
      memcpy(tau, starts + 6 * is, 6 * sizeof(BC_CPREC));
      mean_inst = vavrycuk_scan(n, nrm, ne, fmin, finc, nfric, n_iter, tau, &fo,
				choice, cbest, sel);
    } else {			/* reference tensor variant from the averaged start */
      memcpy(tau, starts, 6 * sizeof(BC_CPREC));
      mean_inst = vavrycuk_scan_ref(n, nrm, ne, fmin, finc, nfric, n_iter, tau, &fo,
				    choice, sel);
    }
    if ((is == 0) || (mean_inst > best_mean)) { /* first start always sets the result */
      best_mean = mean_inst; best_fo = fo;
      memcpy(best_tau, tau, 6 * sizeof(BC_CPREC));
      memcpy(best_sel, sel, n);
    }
  }
  /* -------- normalize the returned tensor once, per request -------- */
  if (norm_type == BC_STRESS_NORM_TENSOR)
    normalize_tens6(best_tau, best_tau);
  else
    max_ev_normalize_tens6(best_tau, best_tau);    /* BC_STRESS_NORM_EV (default) */

  /* -------- outputs -------- */
  memcpy(stress, best_tau, 6 * sizeof(BC_CPREC));
  calc_eigensystem_vec6(best_tau, sg, sv, BC_FALSE, BC_FALSE); /* ascending */
  *shape_ratio = (sg[0] - sg[1]) / (sg[0] - sg[2]);
  *fopt  = best_fo;
  *minst = best_mean;
  if (sel_out) {
    for (j = j6 = 0; j < n; j++, j6 += 6) {
      memcpy(sel_out + j6, angles + j6, 6 * sizeof(BC_CPREC));
      if (best_sel[j] == 1)        /* selection under the final tau */
	swap_angles(sel_out + j6);
    }
  }
  free(ne); free(nrm); free(starts); free(choice); free(cbest); free(sel); free(best_sel);
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

/* ascending comparator for BC_CPREC, used by the bootstrap below */
static int vavrycuk_cmp_dbl(const void *a, const void *b)
{
  BC_CPREC da = *(const BC_CPREC *)a, db = *(const BC_CPREC *)b;
  return (da < db) ? -1 : ((da > db) ? 1 : 0);
}

/*
  bootstrap error estimate for the best fit friction.

  the friction that maximizes mean instability is the least resolved
  output of the inversion: the mean instability versus friction curve is
  typically very flat near its maximum, so a single best friction can be
  far more precise looking than the data support. to attach an
  uncertainty, resample the n events with replacement n_boot times,
  rerun the friction scan on each resample, and summarize the spread of
  the resulting optima. this captures both sources of friction
  uncertainty at once: the number of events in the bin, and the flatness
  of the peak, which makes the argmax wander between resamples.

  inputs:
    n, angles (6 per event, radians), weights
    fmin, fmax, finc   friction grid. a coarse finc (e.g. 0.02) is fine
                       and recommended: the bootstrap spread dominates
                       the grid resolution and a coarse grid keeps the
                       cost down (each resample is a full scan).
    n_iter, n_real     as in stress_inversion_vavrycuk
    n_boot             number of bootstrap resamples (e.g. 200). if < 1,
                       only the point estimate is returned with zero error.
    seed               RNG seed (ran2 convention, pass a negative long)
  outputs:
    fbest    best friction from the full, unresampled data (point estimate)
    fmean    mean of the bootstrap optima
    fstd     standard deviation of the bootstrap optima
    f16,f84  16th and 84th percentile bounds, a roughly 1 sigma interval
             that does not assume the distribution is symmetric

  note: the returned interval is quantized to the friction grid, and
  individual resamples can be degenerate (events drawn repeatedly); both
  wash out for n_boot of order 100 or more and do not bias the spread.
*/
void vavrycuk_friction_error(int n, BC_CPREC *angles, BC_CPREC *weights,
                             BC_CPREC fmin, BC_CPREC fmax, BC_CPREC finc,
                             int n_iter, int n_real, int n_boot, long int *seed,
                             BC_CPREC *fbest, BC_CPREC *fmean, BC_CPREC *fstd,
                             BC_CPREC *f16, BC_CPREC *f84)
{
  BC_CPREC *rangles, *rweights, *fb, stress[6], shape_ratio, minst, fopt, s1, s2;
  size_t asize = 6 * sizeof(BC_CPREC) * n;
  int b, j, j6, idx, i16, i84, n_real_boot, nvalid = 0;

  /* point estimate from the full data */
  stress_inversion_vavrycuk(n, angles, weights, fmin, fmax, finc,
                            n_iter, n_real, seed, stress, &shape_ratio,
                            fbest, &minst, NULL, BC_STRESS_NORM_EV);

  if (n_boot < 1) {                 /* point estimate only */
    *fmean = *fbest; *fstd = 0.0; *f16 = *fbest; *f84 = *fbest;
    return;
  }
  rangles  = (BC_CPREC *)malloc(asize);
  rweights = (BC_CPREC *)malloc(sizeof(BC_CPREC) * n);
  fb       = (BC_CPREC *)malloc(sizeof(BC_CPREC) * n_boot);
  if ((!rangles) || (!rweights) || (!fb)) BC_MEMERROR("vavrycuk_friction_error");

  /* the resamples use a single random realization as the start, i.e.
     one scan per resample instead of n_real + 1. the result can depend
     on the start for small bins, which adds some scatter to the
     bootstrap optima; in tests with 5 and 10 realizations the mean and
     standard deviation of the optima changed by amounts comparable to
     the bootstrap noise, at 5 to 10 times the cost. the point estimate
     above keeps the full n_real. */
  n_real_boot = 1;
  s1 = s2 = 0.0;
  for (b = 0; b < n_boot; b++) {
    /* draw n events with replacement */
    for (j = j6 = 0; j < n; j++, j6 += 6) {
      idx = (int)(BC_RGEN(seed) * (BC_CPREC)n);
      if (idx >= n) idx = n - 1;     /* guard the measure-zero endpoint */
      memcpy(rangles + j6, angles + idx * 6, 6 * sizeof(BC_CPREC));
      rweights[j] = weights[idx];
    }
    stress_inversion_vavrycuk(n, rangles, rweights, fmin, fmax, finc,
                              n_iter, n_real_boot, seed, stress, &shape_ratio,
                              &fopt, &minst, NULL, BC_STRESS_NORM_EV);
    /* resamples that draw too few distinct events can give a singular
       system and a non-finite solution, those are skipped */
    if (!finite(minst))
      continue;
    fb[nvalid++] = fopt;
    s1 += fopt;
    s2 += fopt * fopt;
  }
  if (nvalid == 0) {
    *fmean = *fstd = *f16 = *f84 = NAN;
    free(rangles); free(rweights); free(fb);
    return;
  }
  *fmean = s1 / (BC_CPREC)nvalid;
  *fstd  = sqrt(fabs(s2 / (BC_CPREC)nvalid - (*fmean) * (*fmean)));

  /* percentile bounds from the sorted bootstrap distribution */
  qsort(fb, (size_t)nvalid, sizeof(BC_CPREC), vavrycuk_cmp_dbl);
  i16 = (int)(0.16 * (BC_CPREC)nvalid);
  i84 = (int)(0.84 * (BC_CPREC)nvalid);
  if (i84 >= nvalid) i84 = nvalid - 1;
  *f16 = fb[i16];
  *f84 = fb[i84];

  free(rangles); free(rweights); free(fb);
}
