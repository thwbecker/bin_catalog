#include "catalog.h"
/* 
   calculate gutenberg richter following Marzocci 2013

   m[nm] magnitudes, with dm spacing (e.g. 0.1)
   mmin: magnitude threshold (e.g. 2)

   returns b and error thereof, sb
*/
void calc_b_value_marzocci(BC_CPREC *m, int nm, BC_CPREC dm, BC_CPREC mmin,
			   BC_CPREC *b, BC_CPREC *sb)
{
  // M_LN10 = 2.30258509299405
  BC_CPREC mu,sum,std,p; 
  int i,nm_ac;
  for(nm_ac=0,mu=0.,i=0;i<nm;i++){ /* mean for magnitudes above
				   completeness */
    if(m[i] >= mmin){
      mu += m[i];
      nm_ac++;
    }
  }
  mu /= (BC_CPREC)nm_ac;
  
  if(dm > 0){			/* for magnitudes with dm spacing */
    /* Tinti and Mulargia (1987) approach */
    /* this doesn't quite work for sb, why? */
    p = 1. + dm/(mu-mmin); 	/* eq. 3.10 */
    *b = 1/(M_LN10*dm) * log(p);		    /* b estimate, eq. 3.9 */
    *sb = -(1-p)/(M_LN10*dm*sqrt((BC_CPREC)nm_ac*p)); /* error */
  }else{
    for(std=0.,i=0;i < nm_ac;i++){
      sum = m[i]-mu;
      std += sum*sum;
    }
    std /= ((BC_CPREC)nm_ac)*(((BC_CPREC)nm_ac)-1.);
    std = sqrt(std);
    *b = 1/(M_LN10*(mu-mmin));
    *sb = 2.3*(*b)*(*b)*std;
  }
}
/* 
   calculate gutenberg richter 

   m[nm] magnitudes, with dm spacing (e.g. 0.1)
   mmin: magnitude threshold (e.g. 2)

   returns b and error thereof, sb
*/
void calc_b_value_thomas(BC_CPREC *m, int nm, BC_CPREC dm, BC_CPREC mmin,
			 BC_CPREC *b, BC_CPREC *sb)
{
  // M_LN10 = 2.30258509299405
  BC_CPREC mu,sum,sum2,std; 
  int i,nm_ac;
  for(nm_ac=0.,sum=sum2=0.,i=0;i<nm;i++){	/* mean for magnitudes above
					   completeness */
    if(m[i] >= mmin){
      sum += m[i];
      sum2 += m[i]*m[i];
      nm_ac++;
    }
  }
  mu = sum/(BC_CPREC)nm_ac;		/* mean magnitude */
  std = std_quick(nm_ac, sum,sum2);
  /* Utsu (1966), Bender (1983) dm corrected equation, eq. 3.1 of
     Marzocci */
  *b =  ( 1. / ( mu - (mmin - dm/2)) )/M_LN10; /* b value (1/ln(10) = log10(e)) */
  /* standard deviation */
  *sb = 2.3 * std/sqrt((BC_CPREC)nm_ac) * (*b)*(*b); /* Shi and Bolt (1982), eq. 2.5 */
}
/* 
   m[nm] magnitudes, mmin: magnitude threshold

   simple max likelihood from Aki, this is not a good idea

   dm unused!
*/
void  calc_b_value_ml(BC_CPREC *m, int nm, BC_CPREC dm, BC_CPREC mmin,BC_CPREC *b, BC_CPREC *sb)
{
  // M_LN10 = 2.30258509299405
  BC_CPREC mu; 
  int i,nm_ac;
  for(nm_ac=0,mu=0.,i=0;i<nm;i++){	/* mean for magnitudes above
					   completeness */
    if(m[i] > mmin){
      mu += m[i];
      nm_ac++;
    }
  }
  mu /= (BC_CPREC)nm_ac;
  
  *b = 1/(mu-mmin)/M_LN10;	/* eq. 2.3 from Marzocci */
  *sb = *b/sqrt((BC_CPREC)nm_ac);	/* eq. 2.4  */
}
