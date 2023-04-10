#include "eigen.h"
/*
  
  given real, symmetric a 3 x 3 matrix A where only the RR, RT, RP, TT, TP, and PP
  components have to be filled (the corresponding other ones will be overwritten)

  since this is for symmteric matrices, doesn't matter if we have C or
  FORTRAN storage

  calculate the eigenvalues eval[3]
  and the corresponding eigenvectors, evec is 3x3

  val will have the eigenvalues in ascending order, 
  e1 = val[2],e2 = val[1],e3 = val[0]
  or 
  e1 = val[E1],e2 = val[E2],e3 = val[E3]

  with e1 > e2 > e3

  vec will be each normalized to unity 

  r, theta, phi
  
  vec[6,7,8] is the vector that corresponds to the highest  eigenvalue,     val[2]
  vec[3,4,5] is the vector that corresponds to the intermediate eigenvalue, val[1]
  vec[0,1,2] is the vector that corresponds to the smallest eigenvalue,     val[0]


  WARNING:

  if icalc_vectors is not set to TRUE, will only calculate the eigenvalues

  uses EISPACK routine

  we have also defined macros E1, E2, E3 that refer to the indices
  2, 1, and 0, respectively (see above)

*/
void calc_eigensystem_sym_9(COMP_PRECISION *a,COMP_PRECISION *eval,
			    COMP_PRECISION *evec,
			    short int icalc_vectors)
{
  static int n=3;// dimension, don't use a define since we want to pass n to 'rs'
  int i,ierr,matz,j;
  COMP_PRECISION fv1[3],fv2[3],loca[9];
  //
  // assign the symmetric values of the matrix
  // this is unnecessary, do it for safety
  // (will, however, overwrite the presumed symmetric entries)
  //
  a[TR]=a[RT];a[PR]=a[RP];a[PT]=a[TP];
  memcpy(loca,a,(size_t)9*sizeof(COMP_PRECISION));
  //
  // EISPACK matz flag, 0: only eigenvalues, !=0: values + vectors
  matz=(icalc_vectors)?(1):(0);
  // call the appropriate symmetric matrix eigensystem routine
  // from EISPACK
  SROUT(&n,&n,loca,eval,&matz,evec,fv1,fv2,&ierr);
  if(ierr){
    fprintf(stderr,"calc_eigensystem_sym: runtime error %i in routine\n",ierr);
    exit(-1);
  }
  for(i=0;i < 3;i++){
    /*  */
    if(!finite(eval[i]))
      fprintf(stderr,"calc_eigensystem_sym: WARNING: eigenvalue %i not finite\n",i);
    if(icalc_vectors){
      /* make sure eigenvectors are pointing in standard direction */
      if(evec[i*3+0] < 0)
	for(j=0;j<3;j++){
	  evec[i*3+j] = -evec[i*3+j];
	}
    }
  }
}

void calc_eigensystem_sym_3x3(COMP_PRECISION a[3][3],COMP_PRECISION *eval,
			      COMP_PRECISION *evec,
			      short int ivec)
{
  COMP_PRECISION b[9];
  b[XX] = a[0][0];
  b[XY] = a[0][1];
  b[XZ] = a[0][2];
  b[YY] = a[1][1];
  b[YZ] = a[1][2];
  b[ZZ] = a[2][2];
  calc_eigensystem_sym_9(b,eval,evec,ivec);
}





#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 5000

void indexx(int n,COMP_PRECISION *arr,int *indx)
{
  int i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  COMP_PRECISION a;
  if((istack=(int *)malloc(sizeof(int)*NSTACK))==NULL){
    fprintf(stderr,"indexx: memerrror\n");exit(-1);
  }
  for (j=1;j<=n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j];
	a=arr[indxt];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	SWAP(indx[l+1],indx[ir])
	  }
      if (arr[indx[l]] > arr[indx[ir]]) {
	SWAP(indx[l],indx[ir])
	  }
      if (arr[indx[l+1]] > arr[indx[l]]) {
	SWAP(indx[l+1],indx[l])
	  }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	SWAP(indx[i],indx[j])
	  }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK) {
	fprintf(stderr,"increase NSTACK in indexx\n");exit(-1);}
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free(istack);
}
#undef M
#undef NSTACK
#undef SWAP


