#include <stdio.h>
#include <math.h>
#include "catalog.h"


  
void michael_leasq(BC_CPREC *a,int m,int n,BC_CPREC *x,BC_CPREC *b,
		   BC_CPREC *a2,BC_CPREC *c,BC_CPREC *psis) /* finds the least
						      squares solution
						      of ax=b */

/* a: the coefficients matrix with n rows and m columns */
/* x: the solution vector of length m */
/* b: the constant vector of length n */
 /* a2: a square matrix of size m for internal use */
/* c: vector of length m for internal use */
 /* m, n: see above */
 /*psis: pointer to the variance */

/* steps 1 a2= a transpose a */
/*       2 c= a transpose b */
/*       3 solve a2x=c by gaussian elimination */

{
  michael_atransa(a,m,n,a2);
  michael_atransb(a,m,n,b,c);
  michael_gaus(a2,m,x,c);
  michael_sigsq(a,m,n,x,b,psis);
}

#define MS_SUB(I,J) (J+I*m)

void michael_gaus(BC_CPREC *a,int m,BC_CPREC *x,BC_CPREC *b) /* solves ax=b for x by gaussian elimination */
{
  int i,i2,i3;
  BC_CPREC d,hold,fact;
  /* take care of special cases */
  if(m < 2){
    x[1]=0.;
    if(m == 1)
      x[1] = b[1]/a[1];
    return;
  }

  for(i=0;i<m;++i){     /* loop for each pivot */
    for(i2=i+1;i2<m;++i2){   /* loop for each row below a pivot */
      /* see if element below pivot is 0 */
      if(a[MS_SUB(i2,i)]==0.)
	continue;
      
      /* if element below pivot > pivot flop rows  */
      if(dabs(a[MS_SUB(i2,i)]) > dabs(a[MS_SUB(i,i)])){
	hold=b[i];
	b[i]=b[i2];
	b[i2]=hold;
	for(i3=i;i3<m;++i3){
	  hold=a[MS_SUB(i,i3)];
	  a[MS_SUB(i,i3)]=a[MS_SUB(i2,i3)];
	  a[MS_SUB(i2,i3)]=hold;
	}
      }
      
      /* do the elimination */ 
      fact=a[MS_SUB(i2,i)]/a[MS_SUB(i,i)];
      a[MS_SUB(i2,i)]=0.;
      for(i3=i+1;i3<m;++i3)
	a[MS_SUB(i2,i3)]=a[MS_SUB(i2,i3)]-fact*a[MS_SUB(i,i3)];
      b[i2]=b[i2]-fact*b[i];
    }
  }
  /* solve the equations */
  x[m-1]=b[m-1]/a[m*m-1];
  for(i=m-2;i> -1;--i){
    d=b[i];
    for(i2=i+1;i2<m;++i2)
      d=d-x[i2]*a[MS_SUB(i,i2)];
    x[i] = d/a[MS_SUB(i,i)];
  }
}


BC_CPREC dabs(BC_CPREC cc)
{  
  BC_CPREC dd;
  dd=cc;  
  if(dd<0)
    dd= -dd; 
  return (dd); 
}



void michael_atransa(BC_CPREC *a,int m,int n,BC_CPREC *b) /* computes b=a transpose*a  */

  
{
  int i,j,k;
  BC_CPREC bb;
  
  for(i=0;i<m;++i){
    for(j=i;j<m;++j){
      bb=0;
      for(k=0;k<n;++k) bb=bb+a[MS_SUB(k,i)]*a[MS_SUB(k,j)];
      b[MS_SUB(i,j)]=bb;
      b[MS_SUB(j,i)]=bb;
    }
  }
}



void michael_atransb(BC_CPREC *a,int m,int n,BC_CPREC *b,BC_CPREC *c) /* computes c= atranspose * b */
{   
  int i,i2;

  for(i=0;i<m;++i){ 
    c[i]=0.;
    for(i2=0;i2<n;++i2)
      c[i]=c[i]+a[MS_SUB(i2,i)]*b[i2];
  }

}

void michael_sigsq(BC_CPREC *a,int m,int n,BC_CPREC *x,BC_CPREC *b,BC_CPREC *psis) /* computes the variance of a */
/* single observation */
/* matrix of n rows and m columns */
/* data vector length n*/
/* solution vector length m */
/* where to put answer */
/* see above */
{ 
  int i,j;    /* loop variables */
  BC_CPREC y,z,z2;     /* sum variables */
  z=0;
  for(i=0;i<n;++i){  /* loop over rows */
    y=0;
    for(j=0;j<m;++j)
      y+= a[MS_SUB(i,j)]*x[j];
    z+= (b[i]-y)*(b[i]-y);
  }
  if(n!=m){
    z2=(BC_CPREC)n-m;
    z=z/z2;
  }
  *psis=z;
}
