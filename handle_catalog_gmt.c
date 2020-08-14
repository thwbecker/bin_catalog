#include "catalog.h"
/* 
   
   this is the part of the code which depends on GMT psmeca utilities
   this is regrettably for GMT 4.5.18
   
*/
/* for fault planes */
#include "gmt.h"
#include "meca.h"		/* GTM meca */
/* from meca stuff */
void GMT_momten2axe(struct M_TENSOR ,struct AXIS *,struct AXIS *,struct AXIS *);
void axe2dc(struct AXIS ,struct AXIS ,struct nodal_plane *,struct nodal_plane *);
/* 
   convert a moment tensor to double couple and compute the best fit fault planes
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
