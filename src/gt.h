#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <stdbool.h> // where bool, true, and false are defined 10/10/2005
#include <ctype.h>
#include <R.h>

#define NR_END 1
#define FREE_ARG char*
//#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
//#define FMAX(a,b) (maxarg1=(a), maxarg2=(b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
//#define IMIN(a,b) (iminarg1=(a), iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define square(x) x*x        /* Compute the SQUARE of x */
#define cube(x) x*x*x        /* Compute the CUBE of x */
#define SQRT2 sqrt(2) /* 1.414213562 */
#define LINELEN 10000

#define SWAP(a, b) itemp=(a);(a)=(b);(b)=itemp;
#define MM 7
#define NSTACK 80000

#define Max(a,b) ((a)<(b)?(b):(a))
#define Min(a,b) ((a)>(b)?(b):(a))

typedef char string[80];
struct tm *Time();

static float maxarg1, maxarg2;
static int iminarg1, iminarg2;
float pythag(float, float);
static float sqrarg;
double erfcc(double);
void Error(char *);
long records(FILE *);
int max_int(int *, int);
int min_int(int *, int);
double erfcc(double);
void timedate();
void ftimedate(FILE *);
void swap(int *, int *);
void dswap(double *, double *);
double DPower(double, int);
void iSort(int *, int, int);
int iVectorMultiply(int *, int, int);

// Normal Distribution Tool Functions
// ==================================
double ndf(double);
double normal_cd(double);
double cdfa(double);


// Matrix Allocation Tool Functions
// ==================================
int *ivector(int, int);
int *Ivector(char *, int, int);
void free_ivector(int *, int, int);
float **matrix(int, int, int, int);
void free_matrix(float **, int, int, int, int);
float *vector(long, long);
void free_vector(float *, long, long);
int **imatrix(long, long, long, long);
void free_imatrix(int **, long, long, long, long);
double **dmatrix(long, long, long, long);
double **Dmatrix(char *, long, long, long, long);
void free_dmatrix(double **, long, long, long, long);
double *dvector(long, long);
double *Dvector(char *, long, long);
void free_dvector(double *, long, long);
void GetVector(float **, float *, long, long, long);
float **submatrix(float **, long, long, long, long, long, long);
char **cmatrix(long, long, long, long);
void free_cmatrix(char **, long, long, long, long);
void OrderIndex(int *, int, int *);
void IndexSort(int *, int, int *);



