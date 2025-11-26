/*


  calculate eigenvalues and eigenvectors for general or symmetricx matrices

  a11 a12 a13 ... a1m
  a21 a22 ...
  ...
  an1 ...         anm
  
  if 3D_SYMMETRIC is set, will expect a11 a12 a13 a22 a23 a33 only


*/
#include "eigen.h"

int main(int argc, char **argv)
{
  int n,m,i,j,ierr,matz,*iv1,*sortindex,mn,readcnt,twon,nexec,offset;
  int is_nan;
  COMP_PRECISION *a,*eigenr,*eigeni,*z,*fv1,*fv2,*sortarr,norm,tmp,sum;
#ifndef THREED_SYMMETRIC
  if(argc!=3){
    fprintf(stderr,"%s m n\nreads m x n matrix A from stdin\n",
	    argv[0]);
#ifdef ONLY_VALUES
    fprintf(stderr,"and writes eigenvalues (Re Im) to stdout\n");
#else
    fprintf(stderr,"and writes eigenvalues (Re Im) and eigenvectors to stdout\n");
#endif
    exit(-1);
  }
  sscanf(argv[1],"%i",&m);
  sscanf(argv[2],"%i",&n);
  fprintf(stderr,"%s: expecting %i by %i matrix \n",argv[0],m,n);
#else
  if(argc != 1){
    fprintf(stderr,"%s\nreads 3-D symmetric matrix (11,12,13,22,23,33) from stdin\n",
	    argv[0]);
#ifdef ONLY_VALUES
    fprintf(stderr,"and writes eigenvalues to stdout\n");
#else
    fprintf(stderr,"and writes eigenvalues and eigenvectors to stdout\n");
    fprintf(stderr,"format:\n\neval_3 evec_3^1 evec_3^2 evec_3^3 eval_2 evec_2^1 evec_2^2 evec_2^3 eval_1 evec_1^1 evec_1^2 evec_1^3\n");
    fprintf(stderr,"\neval_1 > eval_2 > eval_3\n\n");
#endif
    exit(-1);
  }
  m=n=3;
#endif
  mn = n * m;
  twon = n*2;
  // only eigenvalues
#ifdef ONLY_VALUES
  matz=0;
#else
  matz=1;
#endif
  if((a=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*mn))==NULL){
    fprintf(stderr,"%s: memerror with matrix\n",argv[0]);
    exit(-1);
  }
  if((z=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*mn+10))==NULL){
    fprintf(stderr,"%s: memerror with matrix\n",argv[0]);
    exit(-1);
  }
  if((sortindex=(int *)malloc(sizeof(int)*n+5))==NULL){
    fprintf(stderr,"%s: memerror with vector\n",argv[0]);
    exit(-1);
  }
  if((eigenr=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n))==NULL){
    fprintf(stderr,"%s: memerror with vector\n",argv[0]);
    exit(-1);
  }
  if((sortarr=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n))==NULL){
    fprintf(stderr,"%s: memerror with vector\n",argv[0]);
    exit(-1);
  }
  if((fv1=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n))==NULL){
    fprintf(stderr,"%s: memerror with vector\n",argv[0]);
    exit(-1);
  }
#ifdef THREED_SYMMETRIC
  if((fv2=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n))==NULL){
    fprintf(stderr,"%s: memerror with vector\n",argv[0]);
    exit(-1);
  }
#else
  if((eigeni=(COMP_PRECISION *)malloc(n*sizeof(COMP_PRECISION)))==NULL){
    fprintf(stderr,"%s: memerror with vector\n",argv[0]);
    exit(-1);
  }
  if((iv1=(int *)malloc(sizeof(int)*n))==NULL){
    fprintf(stderr,"%s: memerror with vector\n",argv[0]);
    exit(-1);
  }
#endif
  nexec=0;
  do{
#ifdef THREED_SYMMETRIC
    readcnt  = fscanf(stdin,FLT_FORMAT,(a       ));// a11
    readcnt += fscanf(stdin,FLT_FORMAT,(a+n     ));// a12
    readcnt += fscanf(stdin,FLT_FORMAT,(a+twon  ));// a13
    readcnt += fscanf(stdin,FLT_FORMAT,(a+n   +1));// a22
    readcnt += fscanf(stdin,FLT_FORMAT,(a+twon+1));// a23
    readcnt += fscanf(stdin,FLT_FORMAT,(a+twon+2));// a33
    readcnt += 3;
    a[1]=a[n];// a21
    a[2]=a[twon];// a31
    a[2+n]=a[twon+1];//a32
#else
    for(readcnt=i=0;i<m;i++)
      for(j=0;j<n;j++)
	readcnt += fscanf(stdin,FLT_FORMAT,(a+i*n+j));
#endif
    /* check for NaN */
    is_nan=0;
    for(i=0;i < mn;i++)
      if(!finite(a[i]))
	is_nan=1;

    if(readcnt == mn){
      nexec++;
      if(!is_nan){
#ifdef THREED_SYMMETRIC
	SROUT(&m, &n, a, eigenr, &matz, z, fv1, fv2, &ierr);
#else
	GROUT(&m, &n, a, eigenr, eigeni, &matz, z, iv1, fv1, &ierr);
#endif
      }else{
	ierr=0;
#ifdef THREED_SYMMETRIC
	for(i=0;i<n;i++){
	  eigenr[i] = nan("");
	}
#else
	for(i=0;i<n;i++){
	  eigenr[i] = nan("");
	  eigeni[i] = nan("");
	}
#endif	
	for(i=0;i < mn;i++)
	  z[i] = nan("");
      }
      
      if(ierr){
	fprintf(stderr,"%s: runtime error %i in rg routine\n",argv[0],ierr);
	exit(-1);
      }else{
#ifdef THREED_SYMMETRIC
	  // output is already in ascending order
	  for(i=0;i<n;i++)
	    sortindex[i]=i;
#else
	  //
	  // sort the eigenvalues in ascending order for the general 
	  // routines
	  //
	  for(j=1,i=0;i<n;i++,j++){
	    sortarr[i]= eigenr[i]*eigenr[i] + eigeni[i]*eigeni[i];
	    sortindex[i]=j;
	  }
	  indexx(n,sortarr-1,sortindex-1);
	  for(i=0;i<n;i++)
	    sortindex[i]--;
#endif
#ifdef ONLY_VALUES
	for(i=0;i<n;i++){
#ifdef THREED_SYMMETRIC
	  // print eigenvalue
	  fprintf(stdout,"%12g  ",eigenr[sortindex[i]]);
#else
	  // print real and imaginary part of eigenvalues
	  fprintf(stdout,"%12g %12g  ",eigenr[sortindex[i]],
		  eigeni[sortindex[i]]);
#endif
	}
	printf("\n");
#else
	for(i=0;i<n;i++){
#ifdef THREED_SYMMETRIC
	  fprintf(stdout,"%12g  ",eigenr[sortindex[i]]);
	  if(i>0){
	    if(eigenr[sortindex[i]] < eigenr[sortindex[i-1]]){
	      fprintf(stderr,"%s: sort error: i: %i e(i): %g e(i-1): %g\n",
		      argv[0],i,eigenr[sortindex[i]],
		      eigenr[sortindex[i-1]]);
	      exit(-1);
	    }
	  }
#else
	  fprintf(stdout,"%12g %12g  ",eigenr[sortindex[i]],
		  eigeni[sortindex[i]]);
#endif
#ifdef NORMALIZE
	  //
	  // normalize eigenvectors and have them point in the 
	  // positive direction
	  //
	  norm = 0.0;
	  for(sum = 0.0,j=0;j < m;j++){
	    offset = sortindex[i]*m+j; /* if I didn't set this like
					  that, the intel compiler
					  would fuck it up with -tpp7
					  -xW options
					*/
	    /* get the norm of the eigenvector */
	    norm  += z[offset] * z[offset];
	    sum +=  z[offset];
	  }
	  norm = sqrt(norm);
	 
	  if(sum < 0)
	    norm *= -1.0;
#else
	  // normalize eigenvectors with last element
	  norm = z[sortindex[i]*m+m-1];
#endif
	  for(j=0;j<m;j++){
	    fprintf(stdout,"%12g  ",z[sortindex[i]*m+j]/norm);
	  }
	}
	printf("\n");
#endif
      }
    }
  }while(readcnt == mn);
  if(nexec == 0)
    fprintf(stderr,"%s: read error\n",argv[0]);
  return 0;
}

