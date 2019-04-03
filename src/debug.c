/*
 * This file is for things written for debugging reasons
 */

#include "gt.h"
#include "analyze.h"
#include "debug.h"

void OutModel(

	int NumFixed,
	int SumFixedLevels,
	int nCov,
	int NumRandom,
	int NumFactors,
	int *Levels,
	int *CumLevels,
        int Cov) {

  int i;
  FILE *fp;

  if((fp=fopen("ModelFile","w"))==NULL)
    Error("Debug <OutModel>: Can't open ModelFile as: w");

  fprintf(fp, "Number of Fixed Factors = %d\n", NumFixed);
  fprintf(fp, "Sum of Fixed Levels = %d\n", SumFixedLevels);
  fprintf(fp, "Number of Covariates = %d\n", nCov);
  fprintf(fp, "Number of Random Factors = %d\n", NumRandom);

  fprintf(fp,"\n\nNumFactors = %d\n",NumFactors);
  for(i=1; i<=NumFactors; i++)
    fprintf(fp,"Levels[%d] = %d      CumLevels[%d] = %d\n", i, Levels[i], i, CumLevels[i]);

  if(Cov) {
    fprintf(fp, "\nExtra covariate equation(s) = %d\n", Cov);
    fprintf(fp, "Total No. of Eq'ns = SumLevels = %d\n", SumLevels);
  }

  fprintf(fp, "\n\n--- E N D ---\n\n");
  fclose(fp);
}

void OutWb(double *Wb) {
  int i;

  FILE *fp;

  if((fp=fopen("WbFile","w"))==NULL)
    Error("Debug <OutWb>: Can't open WbFile");

  for(i=1; i<=NumRecords; i++)
    fprintf(fp,"%16.8f\n", Wb[i]);

  fclose(fp);
}

void OutRHS() {
  int i;
  FILE *fp;

  if((fp=fopen("RHS","w"))==NULL)
    Error("Debug <OutRHS>: Can't open RHSFile");

  for(i=1; i<=SumLevels; i++)
    fprintf(fp,"%12.d %16.8f\n", i, R[i]);

  fclose(fp);
}

// Function to read RHS from disc and use it instead of RHS based on data
// Utilized in the development version for machine learnin purposes
void OverrideRHS() {
  int i;
  FILE *fp;

  if((fp=fopen("rhs","r"))==NULL)
    Error("Debug <OverrideRHS>: Can't open RHSFile");

  for(i=1; i<=SumLevels; i++)
    fscanf(fp,"%lf", &R[i]);

  fclose(fp);
}


 /* Moved from analyze.c */
void PrintAList(int s, int x, FILE *fp) {
  int i;
  /* printf("The %d Structure:\n",s); */
  for(i=StructSp[s].ia[x-1]; i<StructSp[s].ia[x]; i++)
    fprintf(fp,"%d %d %g\n", x, StructSp[s].ja[i-1], StructSp[s].a[i-1]);
}


void OutData(char *dataFileName) {
  int i, j;
  FILE *fpd;

  if((fpd=fopen(dataFileName,"w"))==NULL)
    Error("Debug <OutData>: Can't Open file to output data");

  for(i=1; i<=NumRecords; i++) {
    for(j=1; j<=NumFactors; j++)
      fprintf(fpd, "%d ", FactorsArray[i][j]);

    if(Cov) {
      for(j=1; j<=Cov; j++)
	fprintf(fpd, "%f ", cvrt[i][j]);
    }
    fprintf(fpd,"  %f\n", y[i]);
  }
  fclose(fpd);
}

void CoefMat() {
  int i;
  FILE *fp;

  if((fp=fopen("SC","w"))==NULL)
    Error("Can't Open file for LNL SC -- Coef Mat, Sparse");

  fprintf(fp, "SPARSE FULL-STORED [Full-Rank ord, %d] [nze, %d]\n", StructSp[0].n, StructSp[0].nze);
  for(i=1; i<=SumLevels; i++) PrintAList(0, i, fp);
  fclose(fp);
}

void OutSolVector() {
  int i;
  FILE *fp;

  if((fp=fopen("solVector","w"))==NULL)
    Error("Debug <OutSolVector>: Can't open solVector");

  for(i=1; i<=SumLevels; i++)
    fprintf(fp,"%12.d %16.8f\n", i, b[i-1]);

  fclose(fp);
}







