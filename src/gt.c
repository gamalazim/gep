
// cc -c gt.c

/********************************************
 * File gt.c, contains general subroutines,
 * to use them in other files, include gt.h.
 ********************************************/

#include "gt.h"
#include "debug.h"

// cp -r /usr/local/lib64/R/include ./

void Error(char *error_text) {
  extern void memCleanUp();

  memCleanUp(); // Do Memory clean up - analyze.c/memCleanUp() --
  PROBLEM "%s", error_text RECOVER(NULL_ENTRY);
}

// This function counts records in a file and rewinds the file
long records(FILE *fpd) {
  char line[LINELEN];
  long  count=0, size=sizeof line;

  rewind(fpd);
  fgets(line, size, fpd); count++; /* count a line */

  if(strlen(line) >= LINELEN-1) Error("records(): ---  Increase LINELEN  ---");

  while (fgets(line, size, fpd)) count++;
  rewind(fpd);

  return(count);
}

/*------------------------------------------------------------------
  max_int() receives a one dim'l array, t, along with its length n.
  It returns an int value, the max value found in the array
  ------------------------------------------------------------------*/
int max_int(int *t, int n) {
  int i, Max=0;
  for(i=1; i<= n; i++)
    {
      if(Max < t[i])
	Max = t[i];
    }
  return (Max);
}

/*-----------------------------------------------------------------------
  min_int() receives a one dim'l array, t[1:n], along with its length n.
  It returns an int value, the minimum value found in the array
  -----------------------------------------------------------------------*/
int min_int(int *t, int n)
{
  int i, Min=t[1];
  for(i=1; i<= n; i++)
    {
      if(Min > t[i])
	Min = t[i];
    }
  return (Min);
}

/* Returns the complementary error function erc(x) with error less than
   1.2 x 10^-7 */
double erfcc(double x) {
  float t, z, ans;

  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
      t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
      t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}

/* if called, returns the time of the system */
struct tm *Time() {
  struct tm *systime;
  time_t t;

  t = time(NULL);
  systime = localtime(&t);

  return systime;
}


/* if called, returns the time and the date of the system on two different lines */
void timedate() {
  struct tm *systime;
  time_t t;

  t = time(NULL);
  systime = localtime(&t);

  printf("Time: %.2d:%.2d:%.2d\n", systime->tm_hour, systime->tm_min, systime->tm_sec);
  printf("Date: %.2d/%.2d/%.4d\n\n", systime->tm_mon+1, systime->tm_mday, 1900+systime->tm_year);

}

/* same as timedate but it prints out to a file */
void ftimedate(FILE *fp) {
  struct tm *systime;
  time_t t;

  t = time(NULL);
  systime = localtime(&t);

  fprintf(fp,"Time: %.2d:%.2d:%.2d\n", systime->tm_hour, systime->tm_min, systime->tm_sec);
  fprintf(fp,"Date: %.2d/%.2d/%.2d\n", systime->tm_mon+1, systime->tm_mday, systime->tm_year);

}

/* swap(&x[c], &y[d]) is how to call swap to interchange
    a[c] <==> a[d].
    You can use it with variables too,i.e., swap(&x, &y) */
void swap(int *x, int *y) {
  int temp;
  if(x == y) return;
  temp = *x;
  *x = *y;
  *y = temp;
}

/* Same as swap() but works for doubles */
void dswap(double *x, double *y) {
  double temp;
  temp = *x;
  *x = *y;
  *y = temp;
}

 /*
  * raise x to the nth power --
  * better use the math.h built-in 'double pow(double, double)' as:
  *    pow(x, (double) n) that returns a double
 */
double DPower(double x, int n) {
  int i; double y=x;
  if(n == 0) return 1;
  for(i=2; i<=n; i++)
    y *= x;
  return y;
}

/* sort int array (A) of size (n) in ascending order if opt =1 or descending if opt =2 */
void iSort(int *A, int n, int opt) {
  int i; int j;
  if(opt == 1) {
    for(i=1; i<n; i++)
      for(j=i; j<=n; j++)
	if(A[j] < A[i]) swap(&A[i], &A[j]);
  }
  else if(opt ==2) {
    for(i=1; i<n; i++)
      for(j=i; j<=n; j++)
	if(A[j] > A[i]) swap(&A[i], &A[j]);
  }
  else Error("Unknown Option in iSort\n");
}

// Multiply elements of vector elements[] from 'from' to 'to'
// This function returns at least 1 even if to < from ... needed for some algorithms
int iVectorMultiply(int *elements, int from, int to) {
  int product=1;
  while(from <= to)
    product *= elements[from++];

  return product;
}


////////////////////////////////////////////////
// 200803
// Normal distribution utility functions:
// Was normal.h
////////////////////////////////////////////////
/*=================================================================================
  File: NORMAL.h
  --------------
       THIS FILE CONTAINS THE TOOLS USUALLY USED WITH THE NORMAL DIST. FOR:
       --------------------------------------------------------------------
         1. COMPUTING THE NORMAL DENSITY AT ANY POINT --> ndf()
         2. COMPUTING THE CUMULATIVE DENSITY UP TIL AN UPPER LIMIT --> normal_cd()
         3. COMPUTING THE INVERSE OF 2, i.e. GET THE UPPER LIMIT IF YOU GIVE THE
            CUMULATIVE DENSITY --> cdfa()
===================================================================================
*/

double ndf(double x)
{
  return(0.3989422804014327*exp(-square(x)/2.0));
}

double normal_cd(double z)
/*
      =======================================================
      This function receives a double value corresponding
      to the upper limit of the integration of the normal
      density. It returns the cumulative density from -Inf.
      to the upper limit denoted by uraw!

      *
      *  |  *
      *    |    *
      *     |     *
      *     P|     .*
      *<-------|---->.  *
      *     CUM DENS   .    *
      *       (OUTPUT)   .      *
      ---------------|-----|----------
      -Inf          0.0    z
      (input)

      PS:  "z"  could be positive, negative, or zero!
      ======================================================
*/
{
  double area;

  if (z == 0.0)
    area = 0.5;

  else if (z < 0.0)
    area = 0.5*erfcc(-z/SQRT2);

  else if (z > 0.0)
    area = 0.5 + 0.5*(1.0 - erfcc(z/SQRT2));

  else
    Error("Illegal Upper Limit");

  return(area);

}

/*******************************************************/
         /**************************/

double cdfa(double area)
/*
      ====================================================================
      This routine is used to get the upper limit of integration after
      receiving the area under the normal density from -Inf to the some
      upper limit !!
      ====================================================================
*/

{
  double t,zp,  c0=2.515517, c1=.802853, c2=.010328, d1=1.432788, d2=.189269, d3=.001308;

  if(area == 0.5) return(0.0);

  else if(area < 0.5) {
    area = 0.5 + fabs(area-0.5);
    t = sqrt(-2*log(1-area));
    zp = t - (c0 + c1*t + c2*square(t))/(1 + d1*t + d2*square(t) + d3*cube(t));
    return(-zp);
  }

  else if (area > 0.5){
    t=sqrt(-2*log(1-area));
    zp = t - (c0 + c1*t + c2*square(t))/(1 + d1*t + d2*square(t) + d3*cube(t));
    return(zp);}
}
/*******************************************************/
         /**************************/

///////////////////////////////////
// 200803
// Matrix Set-Up Tool Functions
// Was matrix.h
///////////////////////////////////
/*************************************************************
 * File: MATRIX.h, contains subroutines to allocate and free *
 * memory for int, float, or double arrays!                  *
 *************************************************************/


int *ivector(int nl, int nh) {
  int *v;

  v=(int *) malloc((size_t) (nh-nl+1+NR_END)*sizeof(int));
  if (!v) Error("allocation failure in ivector()");
  return v-nl+NR_END;
}

void free_ivector(int *v, int nl, int nh) {
  free((char*) (v+nl-NR_END));
}

// allocate an int matrix with subscript range m[nrl..nrh][ncl..nch]
int **imatrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  int **m;

  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) Error("allocation failure 1 in imatrix()");
  m += NR_END;
  m -= nrl;

  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) Error("allocation failure 2 in imatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh; i++) m[i]=m[i-1]+ncol;

  return m;
}

// free an int matrix allocated by imatrix()
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch) {
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

// **dmatrix() allocates a double matrix
double **dmatrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  double **m;

  m = (double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if(!m) Error("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;

  m[nrl]= (double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) Error("allocation failure 2 in dmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

  return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch) {
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

double *dvector(long nl, long nh) {
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) Error("allocation failure in dvector()");
  return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh) {
  free((FREE_ARG) (v+nl-NR_END));
}

/* You need to pass the matrix a, that you want to get a column or part of
   it in the vector vec, vec, the particular col of a, start and end rows */
void GetVector(float **a, float *vec, long col, long start, long end) {
  long i, veclen;
  veclen = end-start+1;
  for(i=1; i<=veclen; i++)
    vec[i] = a[i+start-1][col];
}


// point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch]
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
                  long newrl, long newcl) {
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  float **m;

  m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
  if (!m) Error("allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;

  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

  return m;
}

// allocate a char matrix with subscript range m[nrl..nrh][ncl..nch]
char **cmatrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
  char **m;

  m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(char*)));
  if (!m) Error("allocation failure 1 in cmatrix()");
  m += NR_END;
  m -= nrl;

  m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(char)));
  if (!m[nrl]) Error("allocation failure 2 in cmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh; i++) m[i]=m[i-1]+ncol;

  return m;
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch) {
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

// For obtaining the k extremes
// Return a vector or ordered indexes 1 to n
void OrderIndex(int *B, int N, int *idx) {
  int i, j;
  for(i=1; i<=N; i++) idx[i]=i;
  for(i=1; i<N; i++) {
    for(j=i; j<=N; j++)
      if(B[idx[i]] > B[idx[j]])
        swap(&idx[i], &idx[j]);
  }
}

void IndexSort(int *arr, int n, int *indx) {
  int i, indxt, ir=n, itemp, j, k, l=1;
  int jstack=0, *istack;
  int a;
  istack=ivector(1, NSTACK);
  for(j=1; j<=n; j++) indx[j]=j;
  for(;;) {
    if(ir-1 < MM) {
      for(j=l+1; j<=ir; j++) {
        indxt=indx[j];
        a=arr[indxt];
        for(i=j-1; i>=1; i--) {
          if(arr[indx[i]] <= a) break;
          indx[i+1]=indx[i];
        }
        indx[i+1] = indxt;
      }
      if(jstack == 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    } else {
      k = (l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if(arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP(indx[l+1], indx[ir]);
      }
      if(arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l], indx[ir]);
      }
      if(arr[indx[l+1]] > arr[indx[l]]) {
        SWAP(indx[l+1], indx[l]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for(;;) {
        do i++; while(arr[indx[i]] < a);
        do j--; while(arr[indx[j]] > a);
        if(j < i) break;
        SWAP(indx[i], indx[j]);
      }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if(jstack > NSTACK) Error("NSTACK too small in IndexSort.");
      if(ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir = j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l = i;
      }
    }
  }
  free_ivector(istack, 1, NSTACK);
}
