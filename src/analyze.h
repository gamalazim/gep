/*==========
  analyze.h
  ========== */

#define SolStart 1 /* -not finished- 1: for initial guess. 2: for using existing solutions on disk */
#define D_min .0000001
#define SOLEPS 0.0001 /* This will end iteration at convergence to the fifth decimal */
#define DBUG 0
#define eps 1.0e-20
#define SOL 1  /* if(SOL) print out solutions in DenseVectorSolve() */

// Constants to be defined in light of each other
#define MaxNumFactors 200
#define MaxNumLoci 10
#define MaxNumCovFact 10 // Factors with covariance matrices - other than the pedigree factor
#define _StructNum 23 // This should = 3 + MaxNumCovFact + MaxNumLoci

#define CHECKSOL 1
#define Rows2Build 10
#define StopFactor 5  // Factor Used For Calculating Stoppoing Criterion

/// Memory Options
/// ##############
int available; // space in fspak
int max_size; // nonzeros of coeff or cov inv matrices
int AllocMem; // available memory in MB's for building a section or submatrix in C


int solit;
int TM;
int FI; // Fit Intercept
int NCATS, nxi;
int tmrnd;
int RSM;    /*  Relationship Matrix for animals  */
int GRM;     /*  for markers                      */
int *DCM;   /*  Diagonal Cov Matrix              */

int ORDER; /* set it to one to skip ordering */

double vccc, solcc, oldsolcc; /* variance components and solutions conv. criterion for TM */
double lambda; /* Global Variance Ratio (p[0]=1)/p[1] */
double *rrr;
double RELF;
int SM;
int Cov; /* Read from model file as either 1 to indicate 1 covariate, or 0 */

/* int variables */
int NumRecords,NumFactors,NumFixed,NumRandom, NumLoci, NumFR,
SumLevels, SumFixedLevels, SumRandomLevels,yField,
SetUpStage, nind, NOP, StructNum;

/* double variables */
double TSS, yPy, LEV, UEV, spr;
/* double h2 = .3; - Already in hdecl - heritability (total add var) */

/* int 1-dim'l arrays */
int Levels[MaxNumFactors], CumLevels[MaxNumFactors],Field[MaxNumFactors];

/* double 1-dim'l arrays */
double *b, *oldb, *R, *y, *rowval;
double **Q, **L, *vR, *ee, *vee; /* TM vectors and matrices */

/* int 2-dim'l arrays */
int **FactorsArray;

// double 2-dim'l arrays
double **cvrt; // Array of covariate terms to be fit

int CLast;

/* Stuff for the linked lists of C */
int *c_col, *c_link; double *c_val;

int LSCLast, *lsc_col, *lsc_link; double *lsc_val;   /* Stuff for the linked lists of LSC */
int *ilink, *icol, IGLast;

struct sparse {
  int *ia, *ja, n, nze;
  double *a;
} StructSp[_StructNum]; /* 3 sruct's for C, LSC, Ainv, 10 cov mat's, 10 qtl loci */

char modelEq[MaxNumFactors][128];

// Function Def's in analyze.c:
// ============================
     void Out4OtherProgs();
     void RSolutions(char *, char **, char **);
     void DenseVectorSolve(double *);
     void MPM(int, int);
     void MPQ(int, int, int, int);
     void DoCov();
     void InsertAList(int, double *, int, int);
     void Hock2CList(int, int, double, int);
     void CList(int, int, double);
     void CVec(double *, int, int, double);
     int StDex(int, int, int);
     double GetVal(int, int);
     void AddVariance(double *);
     void GetLSC();
     void SetUp(double *);
     void KeepLSC();
     void sform(int, int);
     void cFree(int);
     void SPrint(int, char *);
     void OutList();
     void OutOneList(int);
     void MakeSpC(double *);
     void SortCList(int);
     void Struct2Default(int);
     int GetRank(int *, int);
     void AddInvs(double *);
     void FillME(char *);
     void MMParms(double *, int, int, int, int *, int *, int, double *, double *);
     void Finish();
     void SpC2(double *);
     void TM_Finish();
     void memCleanUp();
     int TakahashiInverse(int);
     void fitIntercept();
     void getMemOptions(int *);
     void writeSparseCoef(int);

     extern double fspak_();
