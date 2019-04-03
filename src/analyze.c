
#include "gt.h"
#include "analyze.h"
#include "debug.h"

void getMemOptions(int *memOptions) {
  available = 300000000;
  max_size  = 30000000;
  AllocMem  = 895;

  if(memOptions[0] != 0 && memOptions[0] < 100000000) available = 100000000;
  if(memOptions[0] >= 100000000) available = memOptions[0];

  if(memOptions[1] != 0 && memOptions[1] < 15000000) max_size = 15000000;
  if(memOptions[1] >= 15000000) max_size = memOptions[1];

  if(memOptions[2] != 0 && memOptions[2] < 512) AllocMem = 512;
  if(memOptions[2] >= 512) AllocMem = memOptions[2];
}


void TM_Finish() {

  if(Q) free_dmatrix(Q, 1, nxi, 1, nxi);
  if(L) free_dmatrix(L, 1, NumRecords, 1, nxi);
  if(vR) free_dvector(vR, 1, NumRecords);
  if(ee) free_dvector(ee, 1, NumRecords);
  if(vee) free_dvector(vee, 1, nxi);
  if(oldb) free_dvector(oldb, 0, SumLevels);

  Q = L = NULL;
  vR = ee = vee = oldb = NULL;
}

// Memory Cleaning After a call to gt.c/Error() --
void memCleanUp() {
  Finish();
  TM_Finish();

  // Free Data Containers --
  if(FactorsArray) free_imatrix(FactorsArray,1,NumRecords,1,NumFR+2*GRM);
  if(y) free_dvector(y,1,NumRecords);
  if(cvrt) free_dmatrix(cvrt, 1, NumRecords, 1, Cov);

  FactorsArray = NULL;
  cvrt = NULL;
  y = NULL;

}

int TakahashiInverse(int w) {
  int flag, needed, *is;
  int iout=1, ioor=2, fnode[1], bnode[1], fx[1], fnseg, bnseg, feqb, irank;
  int opt, avail = available;

  is = ivector(0, avail-1);

  /* RESTART */
  opt = 70;
  fspak_(&opt, &StructSp[w].n, StructSp[w].ia, StructSp[w].ja, StructSp[w].a,b, &flag,
	  &iout, &ioor, &avail, &needed, is, fnode, bnode,  &fnseg, &bnseg, fx, &feqb, &irank);

  /* SYMBOLIC FACTORIZATION */
  opt = 20;
  fspak_(&opt, &StructSp[w].n, StructSp[w].ia, StructSp[w].ja, StructSp[w].a,b, &flag,
	  &iout, &ioor, &avail, &needed, is, fnode, bnode, &fnseg, &bnseg, fx, &feqb, &irank);

  /* NUMERICAL FACTORIZATION */
  opt = 40;
  fspak_(&opt, &StructSp[w].n, StructSp[w].ia, StructSp[w].ja, StructSp[w].a,b, &flag,
	  &iout, &ioor, &avail, &needed, is, fnode, bnode,&fnseg, &bnseg, fx, &feqb, &irank);

  /* -----------------------------------------------------------
      Get the inverse using the Takahashi (1971/or 1973?) method
      and put the NZE in ia, ja, a - yes you loose the original
      matrix!. Notice that fspak computes only needed elements of
      the inverse, ie, those that are not zero in the lower-right
      corner of the C matrix for example.
     ----------------------------------------------------------- */
  opt = 60;
  fspak_(&opt, &StructSp[w].n, StructSp[w].ia, StructSp[w].ja, StructSp[w].a,b, &flag,
	  &iout, &ioor, &avail, &needed, is, fnode, bnode,&fnseg, &bnseg, fx, &feqb, &irank);

  StructSp[w].nze = StructSp[w].ia[StructSp[w].n]-1;

  if(is) { free_ivector(is, 0, avail-1); is = NULL; }
  return flag;
}


void OutTakSparseInverse () {
  int i;
  FILE *fpscinv;

  /* OPEN FILE FOR SPARSE C-INVERSE */
  if((fpscinv=fopen("mixed.sci","w"))==NULL)
    Error("Can't Open file for sparse-stored coefficient matrix inverse");

  fprintf(fpscinv, "SPARSE FULL-STORED [Full-Rank ord, %d] [nze, %d]\n", StructSp[0].n, StructSp[0].nze);
  for(i=1; i<=SumLevels; i++) PrintAList(0, i, fpscinv);

  fclose(fpscinv);
}

void Out4OtherProgs() {
  int i, j, c, k;
  FILE *fps, *fpscinv;

  if((fps=fopen("mixed.sol","w"))==NULL)
    Error("Can't open file for model solutions");

  /* OPEN FILE FOR SPARSE C-INVERSE */
  if((fpscinv=fopen("mixed.sci","w"))==NULL)
    Error("Can't Open file for sparse-stored coefficient matrix inverse");

  fprintf(fpscinv, "SPARSE FULL-STORED [Full-Rank ord, %d] [nze, %d]\n", StructSp[0].n, StructSp[0].nze);
  for(i=1; i<=SumLevels; i++) PrintAList(0, i, fpscinv);

  if(DBUG) OutData("Data.frank");

  if(TM) {
    c = 1+nxi;
    if(nxi)
      for(i=0; i<nxi; i++)
        fprintf(fps,"%d %lf\n", i+1, b[i]); /* thresholds */

    fprintf(fps,"%d %lf\n", 1, b[nxi]); /* Overall mean */
    fprintf(fps,"%d %lf\n", 1, 0.); /* 1st eq deleted */
    for(i=1; i<=Levels[1]; i++)
      fprintf(fps,"%d %lf\n", i+1, b[nxi+i]);
  }
  else {
    for(i=1; i<=Levels[1]; i++) /* 1st factor is complete no mu with LM */
      fprintf(fps,"%d %lf\n", i, b[i-1]);
  }

  for(i=2; i<=NumFixed; i++) {
    fprintf(fps,"%d %lf\n", 1, 0.); /* 1st eq'n deleted */
    for(j=1; j<=Levels[i]; j++) {
      k = c + CumLevels[i-1] + j - 1;
      fprintf(fps,"%d %lf\n", j+1, b[k]);
    }
  }

  for(i=1; i<=NumRandom; i++) {
    for(j=1; j<=Levels[NumFixed+i]; j++) {
      k = c + CumLevels[NumFixed+i-1] + j - 1;
      fprintf(fps,"%d %lf\n", j, b[k]);
    }
  }
  if(Cov) {
    k = SumLevels; /* To make next statement short */
    fprintf(fps,"%d %lf\n", 1, b[k-1]);
  }

  fclose(fps); fclose(fpscinv);
}

void RSolutions(char *solFile, char **ff, char **funCall) {
  int i, j, k, c=0, flag;
  double *mse;
  char varname[20], num[6];
  FILE *fpb;

  if((fpb=fopen(solFile,"w"))==NULL)
    Error("analyze.c/RSolutions(): Can't open a file for solutions");

  // if(DBUG) CoefMat(); // different format
  // writeSparseCoef(0); // write coef matrix to disc - file "LHS"
  // OutRHS(); // write RHS vector to disc - file "RHS"

  flag = TakahashiInverse(0); /* Invert StructSp[0] using Takahashi method */
  if(flag) Error("analyze.c/RSolutions(): TAKAHASHI INVERSE FAILED");

  /* MSE ... */
  mse = dvector(1, SumLevels);
  for(i=1; i<=SumLevels; i++) {
    k = StDex(i, i, 0);
    mse[i] = StructSp[0].a[k];
  }

  fprintf(fpb, "Formula:   %s\n", funCall[0]);
  fprintf(fpb, "Data:      %s\n", funCall[1]);
  fprintf(fpb, "sol.file:  %s\n", solFile);
  fprintf(fpb, "TM:        %s\n", funCall[2]);


  if(FI) { // if(TM) FI = 1
    c = 1;
    if(nxi) { c += nxi; fprintf(fpb, "\n%20s %20s %20s\n", "Thresholds",  "Estimate", "MSE"); }
    for(i=0; i<nxi; i++)
      fprintf(fpb,"%20d %20.6f %20.8f\n",
              i+1, b[i], mse[i+1]); /* thresholds - boundary pts */

    fprintf(fpb, "\n%20s %20s %20s\n", "Intercept", "Estimate", "MSE");
    fprintf(fpb,"%20d %20.8f %20.8f\n",
            1, b[nxi], mse[nxi+1]); /* Overall mean */

    if(NumFixed >= 1) {
      fprintf(fpb, "\n%20s\n%20s %20s %20s\n", ff[0], "Fixed_Fact_1", "Estimate", "MSE");
      fprintf(fpb,"%20d %20.6f %20.8f\n",
	      1, 0., 0.); /* 1st eq deleted */
      for(i=1; i<=Levels[1]; i++)
	fprintf(fpb,"%20d %20.6f %20.8f\n",
		i+1, b[nxi+i], mse[nxi+i+1]);
    }
  }
  else {
    if(NumFixed >= 1)  {
      fprintf(fpb, "\n%20s\n%20s %20s %20s\n", ff[0], "Fixed_Fact_1", "Estimate", "MSE");
      for(i=1; i<=Levels[1]; i++) /* 1st factor is complete no mu with LM */
	fprintf(fpb,"%20d %20.6f %20.8f\n", i, b[i-1], mse[i]);
    }
  }
  for(i=2; i<=NumFixed; i++) {
    strcpy(varname, "Fixed_Fact_"); sprintf(num, "%d", i); strcat(varname, num);
    fprintf(fpb, "\n%20s\n%20s %20s %20s\n", ff[i-1], varname, "Estimate", "MSE");
    fprintf(fpb,"%20d %20.6f %20.8f\n",
            1, 0., 0.); /* 1st eq'n deleted */
    for(j=1; j<=Levels[i]; j++) {
      k = c + CumLevels[i-1] + j - 1;
      fprintf(fpb,"%20d %20.6f %20.8f\n",
              j+1, b[k], mse[k+1]);
    }
  }

  for(i=1; i<=NumRandom; i++) {
    strcpy(varname, "Random_Fact_"); sprintf(num, "%d", i); strcat(varname, num);
    fprintf(fpb, "\n%20s\n%20s %20s %20s %20s\n", ff[NumFixed+i-1], varname, "Estimate", "MSE", "Reliability");
    for(j=1; j<=Levels[NumFixed+i]; j++) {
      k = c + CumLevels[NumFixed+i-1] + j - 1;
      fprintf(fpb,"%20d %20.6f %20.8f %20.8f\n",
              j, b[k], mse[k+1], (1 - mse[k+1]/rrr[i]) );
    }
  }

  for(i=1; i<=NumLoci; i++) {
    strcpy(varname, "Random_Locus_"); sprintf(num, "%d", i); strcat(varname, num);
    if(i == 1) fprintf(fpb, "\n%20s\n%20s %20s %20s %20s\n", ff[NumFR], varname, "Estimate", "MSE", "Reliability");
    if(i > 1) fprintf(fpb, "\n%20s %20s %20s %20s\n", varname, "Estimate", "MSE", "Reliability");
    for(j=1; j<=Levels[NumFR+i]; j++) {
      k = c + CumLevels[NumFR+i-1] + j - 1;
      fprintf(fpb,"%20d %20.6f %20.8f %20.8f\n",
              j, b[k], mse[k+1], (1 - mse[k+1]/rrr[NumRandom+i]) );
    }
  }

  for(k = (SumLevels-Cov+1); k<=SumLevels; k++) {
    j = (k - SumLevels + Cov);
    strcpy(varname, "Covariate_"); sprintf(num, "%d", j); strcat(varname, num);

    if(NumLoci) fprintf(fpb,"\n%20s\n%20s %20s %20s\n", ff[NumFR + j], varname, "Estimate", "MSE");
    else fprintf(fpb,"\n%20s\n%20s %20s %20s\n", ff[NumFR + j - 1], varname, "Estimate", "MSE");

    fprintf(fpb,"%20d %20.6f %20.8f\n",
            1, b[k-1],mse[k]);
  }

  // Out4OtherProgs();
  // OutTakSparseInverse(); // Takahashi sparse inverse 3/12/2009
  fclose(fpb);
  free_dvector(mse, 1, SumLevels);

} //RSolutions() with file name from R

void DenseVectorSolve(double *chromosome) {
  int i;

  int flag, needed, *is;
  int iout=1, ioor=2, fnode[1], bnode[1], fx[1], fnseg,
  bnseg, feqb, irank;
  int opt, avail = available;

  is = ivector(0, avail-1);

  /*-------------
  find ordering
  ------------- */

  /* ORDER = 1; */   /* Don't Order */
  if(!ORDER) { /* ordering once */
    ORDER++;
    printf("Ordering ...\n");
    opt = 10;
    fspak_(&opt, &StructSp[0].n, StructSp[0].ia, StructSp[0].ja,
           StructSp[0].a,
            b, &flag, &iout, &ioor, &avail, &needed, is, fnode, bnode,
            &fnseg, &bnseg, fx, &feqb, &irank);
    printf("Ordering   D O N E ...\n");
    /* printf("SumLev = %d\n", SumLevels); */
  }

  /*--------------------------------------------
  use the ordering vectors stored in unit ioor
  --------------------------------------------*/
  else {
    opt = 70;
    fspak_(&opt, &StructSp[0].n, StructSp[0].ia, StructSp[0].ja, StructSp[0].a,
            b, &flag, &iout, &ioor, &avail, &needed, is, fnode, bnode,
            &fnseg, &bnseg, fx, &feqb, &irank);
  }
  /*----------------------
  symbolic factorization
  ----------------------*/
  opt = 20;
  fspak_(&opt, &StructSp[0].n, StructSp[0].ia, StructSp[0].ja, StructSp[0].a,
          b, &flag, &iout, &ioor, &avail, &needed, is, fnode, bnode,
          &fnseg, &bnseg, fx, &feqb, &irank);

  /*-----------------------
  numerical factorization
  -----------------------*/
  opt = 40;
  fspak_(&opt, &StructSp[0].n, StructSp[0].ia, StructSp[0].ja, StructSp[0].a,
          b, &flag, &iout, &ioor, &avail, &needed, is, fnode, bnode,
          &fnseg, &bnseg, fx, &feqb, &irank);

  /* -----------------------------
  Solve - output solutions to b
  ----------------------------- */
  /* This is the cheapest step */
  for(i=1; i<=SumLevels; i++) b[i-1] = R[i];
  opt = 50;
  fspak_(&opt, &StructSp[0].n, StructSp[0].ia, StructSp[0].ja, StructSp[0].a,
          b, &flag, &iout, &ioor, &avail, &needed, is, fnode, bnode,
          &fnseg, &bnseg, fx, &feqb, &irank);

  if(CHECKSOL) {
    double *b2; b2 = dvector(0, SumLevels-1);
    for(i=1; i<=SumLevels; i++) b2[i-1] = R[i];
    opt = 56;
    fspak_(&opt, &StructSp[0].n, StructSp[0].ia, StructSp[0].ja, StructSp[0].a,
            b, &flag, &iout, &ioor, &avail, &needed, is, fnode, bnode,
            &fnseg, &bnseg, b2, &feqb, &irank);
    // printf("Mean(C-RHS) = %g\n", b2[0]);
    if(flag) printf("WARNING: solutions are not numerically correct in current iteration\n");
    free_dvector(b2, 0, SumLevels-1); b2 = NULL;
  }

  if(is) { free_ivector(is, 0, avail-1); is = NULL; }
}

void MPM(int f1, int f2) {
  int i,j, k, reposi=0, reposj=0, offset, sdx;
  double val = 1.0;
  double **mat2alloc;

  offset = 0; if(FI) offset++; if(TM) offset += nxi;
  // if(TM) tm_offset = nxi+1; /* threshold and the overall mean equations */

  if(f1 == f2) { /* do the diagonal elements */
    reposi = offset + CumLevels[f1-1];
    for(k=1; k<=NumRecords; k++) {
      if(FactorsArray[k][f1]) {
        i = reposi + FactorsArray[k][f1];
        if(TM) { R[i] += ee[k]; c_val[i] += vR[k]; }
        else { R[i] += y[k]; c_val[i]++; }
      }
    }
  }
 /*
  * In the following I used ( (8*X/1024) * (Y/1024) ) instead of 8*X*Y/(1024*1024),
  * where X and Y are Levels[f1] and Levels[f2]. This is because 8*X*Y overflows
  * (grows as a huge number that cannot be represented).
  */
  else if( ((8*Levels[f1]/1024) * (Levels[f2]/1024)) < AllocMem) {
    reposi = offset + CumLevels[f1-1];
    reposj = offset + CumLevels[f2-1];

    mat2alloc = dmatrix(reposi+1, reposi+Levels[f1], reposj+1, reposj+Levels[f2]);

    for(i=reposi+1; i<=reposi+Levels[f1]; i++)
      for(j=reposj+1; j<=reposj+Levels[f2]; j++) mat2alloc[i][j]=0;

    for(k=1; k<=NumRecords; k++) {
      if(TM) val = vR[k];
      if(FactorsArray[k][f1] && FactorsArray[k][f2]) {
        mat2alloc[reposi+FactorsArray[k][f1]][reposj+FactorsArray[k][f2]] += val;
	/*
        i = reposi + FactorsArray[k][f1];
        j = reposj + FactorsArray[k][f2];
        if(TM) val = vR[k];
        CList(i, j, val);
    */
      }
    }
    for(i=reposi+1; i<=reposi+Levels[f1]; i++)
      InsertAList(i, mat2alloc[i], reposj+1, reposj+Levels[f2]);
    free_dmatrix(mat2alloc, reposi+1, reposi+Levels[f1], reposj+1, reposj+Levels[f2]);
    mat2alloc = NULL;
  }

  else {
    int sectn = Rows2Build;
    int Nsectn;

    while (8*sectn*Levels[f2]/(1024*1024) < AllocMem) sectn += Rows2Build;
    if(sectn > Rows2Build) sectn -= Rows2Build; else Error("AllocMem is too small ...");
    Nsectn = (int) ceil(Levels[f1]/sectn);

    /* printf("For Factors %d X %d, Allocating sectioned matrices of size %d MB and dim %d X %d\n",
           f1, f2, (int) rint( (8*sectn*Levels[f2]/(1024*1024)) ), sectn, Levels[f2]);
    printf("Where Levels[%d] = %d  and  Levels[%d] = %d\n", f1, Levels[f1], f2, Levels[f2]); */

    if(Levels[f1] > Levels[f2]) {
      fprintf(stderr,"\n\n  .......................\n");
      fprintf(stderr,"\n  ...  W A R N I N G  ...\n");
      fprintf(stderr,"\n  .......................\n");
      fprintf(stderr,
              "This program runs much faster if you rank both fixed and random factors\n");
      fprintf(stderr,
              "so that factors with less levels appear first.\n");
      fprintf(stderr,
              "The program hocks row-linked lists and longer lists minimize uses of the deadly\n");
      fprintf(stderr,
              "time consuming InsertAList(). Interrupt and reorder factors in data file!\n\n");
    }

    /* printf("Nsctn = %d\n", Nsectn);
      Do rows 1 to Nsectn*sectn */
    for(sdx=0; sdx < Nsectn; sdx++) {

      reposi = offset + CumLevels[f1-1];
      reposj = offset + CumLevels[f2-1];

      mat2alloc = dmatrix(reposi+(sdx*sectn)+1, reposi+(sdx+1)*sectn, reposj+1, reposj+Levels[f2]);
      for(i=reposi+(sdx*sectn)+1; i<=reposi+(sdx+1)*sectn; i++)
        for(j=reposj+1; j<=reposj+Levels[f2]; j++) mat2alloc[i][j]=0;

      for(k=1; k<=NumRecords; k++) {
        if(TM) val = vR[k];
        if(FactorsArray[k][f1] > (sdx*sectn) && FactorsArray[k][f1] <= (sdx+1)*sectn
           && FactorsArray[k][f2]) {
          mat2alloc[reposi+FactorsArray[k][f1]][reposj+FactorsArray[k][f2]] += val;
	  /*
          i = reposi + FactorsArray[k][f1];
          j = reposj + FactorsArray[k][f2];
          if(TM) val = vR[k];
          CList(i, j, val);
      */
           }
      }
      for(i=reposi+(sdx*sectn)+1; i<=reposi+(sdx+1)*sectn; i++)
        InsertAList(i, mat2alloc[i], reposj+1, reposj+Levels[f2]);

      free_dmatrix(mat2alloc, reposi+(sdx*sectn)+1, reposi+(sdx+1)*sectn, reposj+1, reposj+Levels[f2]);
      mat2alloc = NULL;
    }

    /*  Do rows Nsectn*sectn +1 to Levels[f1] 'less that one sectn' */
    //printf("Finishing rows %d  to  %d\n", Nsectn*sectn+1, Levels[f1]);
    reposi = offset + CumLevels[f1-1];
    reposj = offset + CumLevels[f2-1];

    mat2alloc = dmatrix(reposi+(Nsectn*sectn)+1, reposi+Levels[f1], reposj+1, reposj+Levels[f2]);
    for(i=reposi+(Nsectn*sectn)+1; i<=reposi+Levels[f1]; i++)
      for(j=reposj+1; j<=reposj+Levels[f2]; j++) mat2alloc[i][j]=0;

    for(k=1; k<=NumRecords; k++) {
      if(TM) val = vR[k];
      if(FactorsArray[k][f1] > (Nsectn*sectn) && FactorsArray[k][f2]) {
        mat2alloc[reposi+FactorsArray[k][f1]][reposj+FactorsArray[k][f2]] += val;
	  /*
        i = reposi + FactorsArray[k][f1];
        j = reposj + FactorsArray[k][f2];
        if(TM) val = vR[k];
        CList(i, j, val);
      */
      }
    }
    for(i=reposi+(Nsectn*sectn)+1; i<=reposi+Levels[f1]; i++)
      InsertAList(i, mat2alloc[i], reposj+1, reposj+Levels[f2]);

    free_dmatrix(mat2alloc, reposi+(Nsectn*sectn)+1, reposi+Levels[f1], reposj+1, reposj+Levels[f2]);
    mat2alloc = NULL;
  }
}



void MPQ(int f1, int f2, int fL1, int fL2) {
  int i,j, k, reposi=0, reposj=0, offset;
  double val = 1.0;

  offset = 0; if(FI) offset++; if(TM) offset += nxi;

  if(fL1 == fL2 && f1 == f2) { /* do the diagonal elements */
    reposi = offset + CumLevels[fL1-1];
    for(k=1; k<=NumRecords; k++) {
      if(FactorsArray[k][f1]) {
        i = reposi + FactorsArray[k][f1];
        if(TM) { R[i] += ee[k]; c_val[i] += vR[k]; }
        else { R[i] += y[k]; c_val[i]++; }
      }
    }
  }

  else {
    reposi = offset + CumLevels[fL1-1];
    reposj = offset + CumLevels[fL2-1];

    for(k=1; k<=NumRecords; k++) {
      if(FactorsArray[k][f1] && FactorsArray[k][f2]) {
        i = reposi + FactorsArray[k][f1];
        j = reposj + FactorsArray[k][f2];
        if(TM) val = vR[k];
        CList(i, j, val);
      }
    }
  }
}



void DoCov() {

  /* Adjust the LS matrix C to the covariates by adding end columns for
  [X Z]'cvrt - also do the R[Last] = cvrt'y */

  int i,j, J, k, cc, reposi, reposj;

  J = 0; if(FI) J++; if(TM) J += nxi;
  //if(TM) J = nxi + 1; else J = 0;

  for(i=1; i<=NumFR; i++) {
    for(k=1; k<=NumRecords; k++) {
      if(FactorsArray[k][i]) /* if !deleted 1st eq'n */  {
        reposi = J + CumLevels[i-1] + FactorsArray[k][i];
        if(TM) {
          for(cc=1; cc<=Cov; cc++) {
            reposj = J + CumLevels[NumFactors] + cc;
            CList(reposi, reposj, (vR[k]*cvrt[k][cc]));
          }
        }
        else {
          for(cc=1; cc<=Cov; cc++) {
            reposj = J + CumLevels[NumFactors] + cc;
            CList(reposi, reposj, cvrt[k][cc]);
          }
        }
      }
    }
  }

  /// qtl' x cvrt
  for(i=NumFR+1; i<=NumFactors; i++) {
    for(k=1; k<=NumRecords; k++) {
      reposi = J + CumLevels[i-1] + FactorsArray[k][NumFR+1];
      if(TM) {
	for(cc=1; cc<=Cov; cc++) {
	  reposj = J + CumLevels[NumFactors] + cc;
	  CList(reposi, reposj, (vR[k]*cvrt[k][cc]));
	  CList(reposi+1, reposj, (vR[k]*cvrt[k][cc]));
	}
      }
      else {
	for(cc=1; cc<=Cov; cc++) {
	  reposj = J + CumLevels[NumFactors] + cc;
	  CList(reposi, reposj, cvrt[k][cc]);
	  CList(reposi+1, reposj, cvrt[k][cc]);
	}
      }
    }
  }

  // cvrt by cvrt submatrix
  reposi = reposj = J + CumLevels[NumFactors];
  // --already set to 0 in SetUp-- for(i=1; i<=Cov; i++) R[reposi+i] = 0.;

  for(k=1; k<=NumRecords; k++) {
    if(TM) {
      for(i=1; i<=Cov; i++) {
        for(j=i; j<=Cov; j++)
          CList(reposi+i, reposj+j, cvrt[k][i]*vR[k]*cvrt[k][j]);
        R[reposi+i] += (cvrt[k][i]*ee[k]);
      }
    }
    else {
      for(i=1; i<=Cov; i++) {
        for(j=i; j<=Cov; j++)
          CList(reposi+i, reposj+j, cvrt[k][i]*cvrt[k][j]);
        R[reposi+i] += (cvrt[k][i]*y[k]);
      }
    }
  }

}


void InsertAList(int rindex, double *row, int from, int to) {
  int i , x=CLast;
  int l = c_link[rindex];

  if(l == 0) c_link[rindex] = CLast;
  else {
    while(1) {
      if(c_link[l] == 0) { c_link[l] = CLast; break; }
      l = c_link[l];
    }
  }

  for(i=from; i<=to; i++) {
    if(row[i] != 0) {
      c_link[CLast] = CLast+1;
      c_col[CLast] = i; c_val[CLast] = row[i];
      CLast++;
    }
  }
  if(CLast > x) c_link[CLast-1] = 0; /* seal the list */
  else c_link[rindex] = 0; // just in case no elements were inserted

  if(CLast > max_size) Error("max_size exceeded");

}

void Hock2CList(int x, int y, double val, int o) {
  /* o: is an optn =1 if you want to insert a value at only col y
  or = 2 if you want to insert the value at col y and y+1 */

  int i;
  while(c_link[x])   /* Roll to end of current list */
    x = c_link[x];

  for(i=0; i<o; i++) {
    c_link[x]=CLast; c_link[CLast]=0;
    c_col[CLast]=y+i; c_val[CLast] = val;
    x = CLast;
    CLast++;
  }

}

// CList adds 'val' to existing elements or creates a new element
void CList(int x, int y, double val) {
  if(val == 0.) return; /* don't add a 0. neither creat a zero elem't
  ^^^^^^^^^^^^^^^^^^^ */
  while(1) {
    if(c_col[x]==y) { c_val[x] += val; return; }
    if(c_col[x] > y) {
      c_link[CLast]=c_link[x]; c_link[x]=CLast;
      c_col[CLast]=c_col[x]; c_val[CLast]=c_val[x];
      c_col[x]=y; c_val[x]=val;
      CLast++; return;
    }

    if(c_link[x] == 0) {
      c_link[x]=CLast; c_link[CLast]=0;
      c_col[CLast]=y; c_val[CLast] = val;
      CLast++; return;
    }
    x = c_link[x];
  }
  return;
}


void CVec(double *CVector, int i, int j, double v) {
  int k = i - 1;
  int m = (k*k - k)/2;
  int re = i - j + 1;
  CVector[m + re] += v;
}


int StDex(int x, int y, int opt) {
  /* index of (x,y) in StructSp[opt] */
  int i;
  for(i=StructSp[opt].ia[x-1]; i<StructSp[opt].ia[x]; i++) {
    if(StructSp[opt].ja[i-1] == y) return i-1;
  }
  printf("WARNING: Element (%d, %d) does not exist in StructSp[%d]\n",
         x, y, opt);
  return 0;
}


double GetVal(int mn, int mx) {
  if(mn==0 || mx==0) return 0.;
  while(mn > 0) {
    if(mx == c_col[mn]) return c_val[mn];
    mn = c_link[mn]; }
    return 0.;
}


void AddVariance(double *p) {

  GetLSC(); /* Get a copy of LSC in c_??? */

  SpC2(p); /* Add Cov inv matrices */

  sform(0, SumLevels); /* convert to ija format */

}

void GetLSC() {
  int i;
  CLast = LSCLast+1;

  for(i=1; i<=LSCLast; i++) {
    c_val[i] = lsc_val[i];
    c_col[i] = lsc_col[i];
    c_link[i] = lsc_link[i];
  }
}

// The Overall Mean Equation: 1'[1 X Z cvrt W] = 1'y
// For Linear Model
// --------------------------------------------------
void fitIntercept() {
  int i, j, m, cc;
  double *vec2alloc;

  vec2alloc = dvector(1, SumLevels);
  for(i=1; i<=SumLevels; i++) vec2alloc[i]=0;

  // 1. diagonal mu x mu = [ 1'1 ]
  c_val[1] =  NumRecords;
  for(i=1; i<=NumRecords; i++) {
    // 2. [1' X Z]
    for(j=1; j<=NumFR; j++)
      if(FactorsArray[i][j]) {
	m = 1 + CumLevels[j-1]+FactorsArray[i][j];
	vec2alloc[m]++;
      }

    // 3. [ 1'qtl ]
    if(GRM) {
      for(j=1; j<=NumLoci; j++) {
        m = 1 + CumLevels[NumFR+j-1]+FactorsArray[i][NumFR+1];
        vec2alloc[m]++;
        vec2alloc[m+1]++;
      }
    }

    // 4. [1' cvrt[][]]
    if(Cov)
      for(cc=1; cc<=Cov; cc++) {
	m = 1+CumLevels[NumFactors]+cc;
	vec2alloc[m] += (cvrt[i][cc]);
      }

    // 5. RHS 1'y
    R[1] += y[i];
  }
  InsertAList(1, vec2alloc, 2, SumLevels);  /* Insert Mu Eqn except the diagonal element */
  free_dvector(vec2alloc, 1, SumLevels);
  vec2alloc = NULL;
}

void SetUp(double *p)
{
  int i,j;
  double temp;

  for(i=1; i<=SumLevels; i++) R[i] = 0.0;
  for(i=1; i<=SumLevels; i++) { c_val[i] = c_link[i]=0; c_col[i]= i; }
  CLast = SumLevels+1;

  /* ----------------------------------------------------------
  Build the Upper triangle of the left hand side matrix (C)
  except the parts corresponding to random QTLs and cvrts.
  Do R[] too.
  ---------------------------------------------------------- */

  // Fit intercept for linear model
  if(FI) fitIntercept();
  for(i=1; i<=NumFR; i++)
    for(j=i; j<=NumFR; j++)
      MPM(i,j);

  if(GRM) {
    // Factors x QTLs
    for(i=1; i<=NumFR; i++)
      for(j=1; j<=NumLoci; j++) {
        MPQ(i, NumFR+1, i, NumFR+j);
        MPQ(i, NumFR+2, i, NumFR+j);
      }
    // QTLs x QTLs
    for(i=1; i<=NumLoci; i++)
      for(j=i; j<=NumLoci; j++) {
        MPQ(NumFR+1, NumFR+1, NumFR+i, NumFR+j);
	MPQ(NumFR+1, NumFR+2, NumFR+i, NumFR+j);
	if(j > i) MPQ(NumFR+2, NumFR+1, NumFR+i, NumFR+j);
	MPQ(NumFR+2, NumFR+2, NumFR+i, NumFR+j);
      }
  }

/*-----------------------------------------------------------
  Do the covariate part to be the last 'Cov' equations
  The equations are set in a model that is NOT reprametrized.
  -----------------------------------------------------------*/
  if(Cov) DoCov();

  TSS = 0;
  for(i=1;i<=NumRecords;i++)
    TSS += square(y[i]);

  if(FactorsArray) free_imatrix(FactorsArray,1,NumRecords,1,NumFR+2*GRM);
  if(y) free_dvector(y,1,NumRecords);
  if(Cov) free_dmatrix(cvrt, 1, NumRecords, 1, Cov);
  FactorsArray = NULL;
  y = NULL;

  /*---------------------------------------------------
  Sparse Struct for the Least Squares C before adding
  variance matrices.
  ---------------------------------------------------*/
  // sform(1, SumLevels);
  // writeSparseCoef(1);

  // KeepLSC(); /* Keep a copy of LSC -- needed only if VC estimation */

  if(NumRandom | NumLoci) SpC2(p);  /* Add var inv matrices to the LSC linked lists */

  if(0) { /* expensive and useless */
    for(i=SumFixedLevels+1; i<=SumLevels; i++)
      SortCList(i);
    /* OutList(); */
  }

  sform(0, SumLevels);
  //SPrint(0, "C.s.2");
  //Struct2Default(0);
  //Error("Fake ... ");
}



void KeepLSC() {

  int i;
  LSCLast = CLast-1;
  lsc_val = dvector(1, LSCLast);
  lsc_col = ivector(1, LSCLast);
  lsc_link = ivector(1, LSCLast);

  for(i=1; i<=LSCLast; i++) {
    lsc_val[i] = c_val[i];
    lsc_col[i] = c_col[i];
    lsc_link[i] = c_link[i];
  }
}

void sform(int opt, int ord) {
  /* -------------------------------------------
  Transforms the C linked lists into ija form
  -------------------------------------------   */
  int i, j, k = CLast-1;

  k = 2*(k-ord)+ord; /* in case of full stored sparse inverse - i.e., both UT and LT stored*/

  StructSp[opt].n = ord;
  StructSp[opt].a = dvector(0,k-1);
  StructSp[opt].ia = ivector(0, ord);
  StructSp[opt].ja = ivector(0,k-1);

  /*  ------------------------
  work  the upper triangle
  ------------------------  */
  k=0;
  for(i=0; i<ord; i++) {
    StructSp[opt].ia[i] = k+1;
    j = i+1;
    while(j > 0) {
      if(c_val[j] != 0.) {
        StructSp[opt].a[k] = c_val[j];
        StructSp[opt].ja[k] = c_col[j];
        k++;
      }
      j = c_link[j];

    }
  }

  StructSp[opt].ia[ord] = k+1;
  StructSp[opt].nze = k;
}

void cFree(int opt) {
  if(StructSp[opt].a) { free_dvector(StructSp[opt].a, 0, StructSp[opt].nze - 1); StructSp[opt].a = NULL; }
  if(StructSp[opt].ja) { free_ivector(StructSp[opt].ja, 0, StructSp[opt].nze - 1); StructSp[opt].ja = NULL; }
  if(StructSp[opt].ia) { free_ivector(StructSp[opt].ia, 0, StructSp[opt].n); StructSp[opt].ia = NULL; }
}

void SPrint (int opt, char *name) {
  /* Print Structure Number 'opt' ! */
  int i;
  FILE *fp;
  if( (fp= fopen(name, "w"))==NULL) {
    printf("SPrint(): %s\n", name);
    Error("can't open file"); }

  printf("\nprinted struct %d into %s\n", opt, name);

    fprintf(fp,"SPARSE FORM, ORDER %d, NZE %d\n", StructSp[opt].n, StructSp[opt].nze);
    for(i=0;i<StructSp[opt].nze;i++) {
      if(i<StructSp[opt].n+1)
        fprintf(fp,"%3.d %3.d %3.d     %g\n",i+1, StructSp[opt].ia[i], StructSp[opt].ja[i], StructSp[opt].a[i]);
      else
        fprintf(fp,"%3.d     %3.d     %g\n",i+1, StructSp[opt].ja[i], StructSp[opt].a[i]);
    }
    fclose(fp);
}

void OutList() {
  /* This is to output c_* list to disk - FILE 'clist' */
  int i;

  FILE *fp;

  if( (fp= fopen("clist", "w"))==NULL)
    Error("can't open clist");
  for(i=1; i<CLast; i++)
    fprintf(fp, "%6d %3d  %3d  %f\n", i, c_col[i], c_link[i], c_val[i]);
  fclose(fp);
}

void OutOneList(int ldx) {
  /* This is to output one UT row (ldx) from c_* list to disk - FILE 'c1list' */
  int l = c_link[ldx]; /* Element next to the diagonal */

  FILE *fp;

  if( (fp= fopen("c1list", "w"))==NULL)
    Error("can't open c1list");

  /* Diagonal Element */
  fprintf(fp, "%d %g\n", c_col[ldx], c_val[ldx]);

  while(l > 0) {
    fprintf(fp, "%d %g\n", c_col[l], c_val[l]);
    l = c_link[l];
  }

  fclose(fp);
}

void MakeSpC(double *p) {
  int i, j, jj, k, l, i1, i2, nze, ord;
  double *dumvec, EV=p[0], alphaP = EV/p[1],
    alphaQ = EV/p[NumRandom];

  ord = StructSp[1].n; /* same order as LS C */
  StructSp[0].a = dvector(0,max_size);
  StructSp[0].ia = ivector(0, ord);
  StructSp[0].ja = ivector(0,max_size);

  StructSp[0].n = ord;
  for(i=0; i<=SumFixedLevels; i++)
    StructSp[0].ia[i] = StructSp[1].ia[i];
  nze = StructSp[1].ia[SumFixedLevels]-1;
  for(i=0; i<nze; i++) {
    StructSp[0].ja[i] = StructSp[1].ja[i];
    StructSp[0].a[i] = StructSp[1].a[i]; }

  k = SumFixedLevels+1; l = SumFixedLevels+nind; nze++;
  dumvec = dvector(k,SumLevels);
  for(i=k; i<=l; i++) {
    for(j=i; j<=SumLevels; j++) dumvec[j]=0;
    i1 = StructSp[1].ia[i-1]; i2 = StructSp[1].ia[i]-1;
    for(j=i1; j<=i2; j++) dumvec[StructSp[1].ja[j-1]]=StructSp[1].a[j-1] ;

    i1 = StructSp[2].ia[i-k]; i2 =StructSp[2].ia[i-k+1]-1;
    for(j=i1; j<=i2; j++) dumvec[k + StructSp[2].ja[j-1]-1]+=alphaP*StructSp[2].a[j-1];

    for(j=i; j<=SumLevels; j++)
      if(dumvec[j] != 0) {
        StructSp[0].ja[nze-1]=j;
        StructSp[0].a[nze-1]=dumvec[j];
        nze++; }

    StructSp[0].ia[i] = nze;
  }

  free_dvector(dumvec,k,SumLevels);
  dumvec = NULL;

  /* We could get rid of the dumvec for the intermediary
    factors. However, I have no brains to do it now */

  if(NumRandom > 2) {
    for(jj=2; jj<NumRandom; jj++) {
      k = CumLevels[NumFixed+jj-1]+1; l = CumLevels[NumFixed+jj];
      dumvec = dvector(k,SumLevels);
      for(i=k; i<=l; i++) {
        for(j=i; j<=SumLevels; j++) dumvec[j]=0;
        i1 = StructSp[1].ia[i-1]; i2 = StructSp[1].ia[i]-1;
        for(j=i1; j<=i2; j++) dumvec[StructSp[1].ja[j-1]]=StructSp[1].a[j-1]+EV/p[jj];

        for(j=i; j<=SumLevels; j++)
          if(dumvec[j] != 0) {
            StructSp[0].ja[nze-1]=j;
            StructSp[0].a[nze-1]=dumvec[j];
            nze++; }

        StructSp[0].ia[i] = nze;
      }
      free_dvector(dumvec,k,SumLevels);
      dumvec = NULL;
    }
  }

  k = CumLevels[NumFactors-1]+1; l = SumLevels;
  dumvec = dvector(k,l);
  for(i=k; i<=l; i++) {
    for(j=i; j<=l; j++) dumvec[j]=0;
    i1 = StructSp[1].ia[i-1]; i2 = StructSp[1].ia[i]-1;
    for(j=i1; j<=i2; j++) dumvec[StructSp[1].ja[j-1]]=StructSp[1].a[j-1];

    i1 = StructSp[3].ia[i-k]; i2 =StructSp[3].ia[i-k+1]-1;
    for(j=i1; j<=i2; j++) dumvec[k + StructSp[3].ja[j-1]-1]+=alphaQ*StructSp[3].a[j-1];

    for(j=i; j<=l; j++)
      if(dumvec[j] != 0) {
        StructSp[0].ja[nze-1]=j;
        StructSp[0].a[nze-1]=dumvec[j];
        nze++; }

    StructSp[0].ia[i] = nze;
  }
  StructSp[0].nze = nze-1;

  free_dvector(dumvec,k,SumLevels);
  dumvec = NULL;

}

void SortCList(int i) {
  int j;
  while(c_link[i] > 0) {
    j = c_link[i];
    while(j > 0) {
      if(c_col[i] > c_col[j]) {
        swap(&c_col[j], &c_col[i]);
        dswap(&c_val[j], &c_val[i]); }
        j = c_link[j]; }
        i = c_link[i]; }
}

void Struct2Default(int c) {
  int i, j, n = StructSp[c].n;
  double **matc;
  FILE *fp;

  if((fp=fopen("C.d", "w"))==NULL)
    Error("can't open file for Default C");

  matc = dmatrix(1, n, 1, n);

  for(i=1; i<=n; i++)
    for(j=1; j<=n; j++)
      matc[i][j] = 0.;

  for(i=0; i<n; i++)
    for(j=StructSp[c].ia[i]; j<StructSp[c].ia[i+1]; j++)
      matc[i+1][StructSp[c].ja[j-1]] = StructSp[c].a[j-1];

  for(i=1; i<=n; i++) {
    for(j=1; j<i; j++)
      fprintf(fp,"%lf ", matc[j][i]);
    for(j=i; j<=n; j++)
      fprintf(fp,"%lf ", matc[i][j]);
    fprintf(fp,"\n");
  }

  free_dmatrix(matc, 1, n, 1, n);
  matc = NULL;
  fclose(fp);
}

int GetRank(int *O, int id) {
  int i;
  for(i=1; i<=nind; i++)
    if(O[i] == id)
      return i;

  /* printf("I could not find %d\n", id); */
  return 0;
}


void AddInvs(double *p) {
  int i, j, jj, k, l, i1, i2, nze;
  double *dumvec, EV=p[0], alphaP = EV/p[1],
  alphaQ = EV/p[NumRandom];

  nze = StructSp[1].ia[SumFixedLevels]-1;

  k = SumFixedLevels+1; l = SumFixedLevels+nind; nze++;
  dumvec = dvector(k,SumLevels);
  for(i=k; i<=l; i++) {
    for(j=i; j<=SumLevels; j++) dumvec[j]=0;
    i1 = StructSp[1].ia[i-1]; i2 = StructSp[1].ia[i]-1;
    for(j=i1; j<=i2; j++) dumvec[StructSp[1].ja[j-1]]=StructSp[1].a[j-1] ;

    i1 = StructSp[2].ia[i-k]; i2 =StructSp[2].ia[i-k+1]-1;
    for(j=i1; j<=i2; j++) dumvec[k + StructSp[2].ja[j-1]-1]+=alphaP*StructSp[2].a[j-1];

    for(j=i; j<=SumLevels; j++)
      if(dumvec[j] != 0) {
      StructSp[0].ja[nze-1]=j;
      StructSp[0].a[nze-1]=dumvec[j];
      nze++; }

      StructSp[0].ia[i] = nze;
  }

  free_dvector(dumvec,k,SumLevels);
  dumvec = NULL;

  /* We could get rid of the dumvec for the intermediary
  factors. However, I have no brains to do it now */

  if(NumRandom > 2) {
    for(jj=2; jj<NumRandom; jj++) {
      k = CumLevels[NumFixed+jj-1]+1; l = CumLevels[NumFixed+jj];
      dumvec = dvector(k,SumLevels);
      for(i=k; i<=l; i++) {
        for(j=i; j<=SumLevels; j++) dumvec[j]=0;
        i1 = StructSp[1].ia[i-1]; i2 = StructSp[1].ia[i]-1;
        for(j=i1; j<=i2; j++) dumvec[StructSp[1].ja[j-1]]=StructSp[1].a[j-1]+EV/p[jj];

        for(j=i; j<=SumLevels; j++)
          if(dumvec[j] != 0) {
          StructSp[0].ja[nze-1]=j;
          StructSp[0].a[nze-1]=dumvec[j];
          nze++; }

          StructSp[0].ia[i] = nze;
      }
      free_dvector(dumvec,k,SumLevels);
      dumvec = NULL;
    }
  }

  k = CumLevels[NumFactors-1]+1; l = SumLevels;
  dumvec = dvector(k,l);
  for(i=k; i<=l; i++) {
    for(j=i; j<=l; j++) dumvec[j]=0;
    i1 = StructSp[1].ia[i-1]; i2 = StructSp[1].ia[i]-1;
    for(j=i1; j<=i2; j++) dumvec[StructSp[1].ja[j-1]]=StructSp[1].a[j-1];

    i1 = StructSp[3].ia[i-k]; i2 =StructSp[3].ia[i-k+1]-1;
    for(j=i1; j<=i2; j++) dumvec[k + StructSp[3].ja[j-1]-1]+=alphaQ*StructSp[3].a[j-1];

    for(j=i; j<=l; j++)
      if(dumvec[j] != 0) {
      StructSp[0].ja[nze-1]=j;
      StructSp[0].a[nze-1]=dumvec[j];
      nze++; }

      StructSp[0].ia[i] = nze;
  }
  StructSp[0].nze = nze-1;

  free_dvector(dumvec,k,SumLevels);
  dumvec = NULL;

}

// Model Terms - perhaps will remain in R
// Needed to write solutions on disc
void FillME(char *line) {
  int j=0, k=1, l;
  while(line[j] != '\0' && line[j] != '\n') {
    l = 0;
    while(line[j] != ' ' && line[j] != '\n') {
      modelEq[k][l] = line[j];
      l++; j++;
    }
    modelEq[k][l] = '\0';
    k++;
    j++;
  }
}

void MMParms(double *RDataVec,
             int RNC,
             int RN,
             int RFI,
             int *Rp,
             int *Roptn,
             int loci,
             double *qtlVar,
             double *randomVar)
{

  /*
  * THE FOLLOWING ARE THE PARAMETERS FOR THE MIXED MODEL
  * OR THINGS THAT ARE DONE ONLY ONCE
  */

  int i, j, k;
  char line[1024];

  RELF = 10.0;
  SM = 1;
  SumFixedLevels=0;
  SumRandomLevels=0;
  nind=0;
  CLast = 0;
  StructNum = _StructNum;
  NumLoci = loci;
  nxi = 0;

  c_val = dvector(1, max_size);
  c_col = ivector(1, max_size);
  c_link = ivector(1, max_size);
  R = dvector(1, max_size);
  b = dvector(0, max_size-1);

  for(i=0; i<=max_size-1; i++) b[i] = 0; // needed, if TM=TRUE

  NumRecords = RN; // records(fpd);

  NumFixed =Rp[0];
  NumRandom =Rp[1];
  NCATS = RNC;
  FI = RFI; // fit intercept, FI=0 || FI = 1

  if(NCATS > 1) { FI = 1; TM = 1; nxi = NCATS - 2; /* Threshold Models analysis */ }
  if(NCATS == 0) TM = 0; /* Linear Models analysis */
  if(NCATS == 1) Error("Can't analyze one category - use 0 for linear");
  if(Rp[2] == 0) RSM = 0; // Do Not Build Ainverse for animals

  GRM = Rp[3]; // 0 for no QTl or 1 for a random QTL component.

  NumFactors = NumFixed + NumRandom + NumLoci;
  NumFR = NumFixed + NumRandom; // to replace old NumFactors

  for(i=1; i<=Rp[0]; i++) Levels[i] = Rp[4+i];

  k = 1;
  if(Rp[2] == 1) { // if animal/sire factor, ranim() was used
    k = 2;
    Levels[NumFixed+1] = Rp[5+Rp[0]+Rp[1]]; // No. ind. in animal factor
    NumRandom++;
    NumFactors++;
    NumFR++;
  }

  DCM = ivector(1, NumRandom);
  if(Rp[2] == 1) DCM[1] = 0; // animal/sire

  // add levels of other/all random factors
  for(i=1; i<=Rp[1]; i++) {
    Levels[NumFixed+k] = Rp[4+NumFixed+i];
    if(Roptn[i-1] == 0) DCM[k] = 1; // no cov matrix was supplied ...
    else DCM[k] = 0;
    k++;
  }

  // add levels of QTL random factors thru the animal factor, which should be provided
  // if(NumLoci) is not needed, loop will not start with NumLoci = 0
  for(i=1; i<=NumLoci; i++) Levels[NumFR + i] = 2*Rp[5+Rp[0]+Rp[1]]; // = 2 * nind
  Cov = Levels[NumFactors+1] = Rp[4]; /* levels of covariate at the end */

  FactorsArray = imatrix(1, NumRecords, 1, NumFR+2*GRM);

  y = dvector(1, NumRecords);
  if(Cov) cvrt = dmatrix(1, NumRecords, 1, Cov);

  /// Read R Response Variable and R Factors */
  for(i=1; i<=NumRecords; i++) y[i] = RDataVec[i-1]; // copy respVar to y
  for(j=1; j<=NumFixed; j++)
    for(i=1; i<=NumRecords; i++)
      FactorsArray[i][j] = (int) RDataVec[j*NumRecords+i-1];

  /// 200803
  // Needed for VC estimation by simplex

  // Marked QTL's are not considered in NumRandom but are considered in NOP.
  NOP = NumRandom+1+NumLoci; // loci = 0 if Rp[3] = 0, ie, in case of no MQTL

  // starting values for the parameters
  rrr = dvector(0, (NOP-1));
  rrr[0] = 1;      /*  Residual */

  k = 1;
  for(j=NumFixed+1; j<=NumFR; j++) {
    for(i=1; i<=NumRecords; i++)
      FactorsArray[i][j] = (int) RDataVec[j*NumRecords+i-1];

    rrr[k] = randomVar[k-1];
    k++;
  }
  lambda = rrr[1];

  // This is not going to start if loci == 0, hence qtlVar[i] will not be accessed!
  for(i=1; i<=loci; i++) rrr[NumRandom+i] = qtlVar[i-1];

  /// Fill up 2 extra colums to needed to build QTL sections in C
  if(NumLoci) {
    for(i=1; i<=NumRecords; i++) {
      FactorsArray[i][NumFR+1] = 2*FactorsArray[i][NumFixed+1]-1;
      FactorsArray[i][NumFR+2] = 2*FactorsArray[i][NumFixed+1];
    }
  }


// Read Covariates
  k = NumRecords + NumFR*NumRecords -1;
  if(Cov) {
    for(j=1; j<=Cov; j++)
      for(i=1; i<=NumRecords; i++) {
      k++;
      cvrt[i][j] = RDataVec[ k ];
      }
  }

  // OutData("RDataVec");

  // Delete 1st equation from 2nd and more fixed factors
  if(FI) k = 1; else k = 2;
  for(i=k; i<=NumFixed; i++) {
    Levels[i]--;
    for(j=1; j<=NumRecords; j++)
      FactorsArray[j][i]--;
  }

  /* build an array of commulative levels */
  CumLevels[0] = 0;
  for(i=1; i<=NumFactors; i++)
    CumLevels[i] = CumLevels[i-1]+Levels[i];
  if(Cov) CumLevels[NumFactors+1] = CumLevels[NumFactors]+Cov;
  SumFixedLevels = CumLevels[NumFixed];
  SumRandomLevels = CumLevels[NumFactors]-SumFixedLevels; // includes QT Loci
  SumLevels = CumLevels[NumFactors]; // includes QT Loci
  if(Cov) SumLevels = CumLevels[NumFactors] + Cov;
  if(FI) SumLevels++;

  OutModel(NumFixed,SumFixedLevels,Cov,NumRandom,NumFactors,Levels,CumLevels, Cov);
  timedate();
}

void Finish() {
  int i;
  /*if(TM) { cFree(0); if(RSM) cFree(2); }
  else  { for(i=0; i<StructNum; i++) cFree(i); }*/


  for(i=0; i<StructNum; i++) cFree(i);

  if(c_val) free_dvector(c_val, 1, max_size);
  if(c_col) free_ivector(c_col, 1, max_size);
  if(c_link) free_ivector(c_link,1, max_size);
  if(R) free_dvector(R, 1, max_size);
  if(b) free_dvector(b, 0, max_size-1);
  if(DCM) free_ivector(DCM, 1, NumRandom);

  c_val = R = b = NULL;
  c_col = c_link = DCM = NULL;
}

void SpC2(double *p) {
  int i, j, jj, k, off=SumFixedLevels;
  int at, to;
  double l;
  double *rl;

  if(NumRandom == 0 && NumLoci == 0) return;

  rl = dvector(1, (NumRandom+NumLoci));
  for(i=1; i<=(NumRandom+NumLoci); i++) rl[i] = p[0]/p[i];

  if(FI) off++; if(TM) off += nxi;

  k = 1;
  if(RSM) { // RSM == 1 for sire model or 2 for animal model
    // nind = StructSp[2].n; -- already set in kinshipCi.c/kinship(()
    k=2; // latter on, Start at next/second random factor
    l = rl[1];
    for(i=0; i<nind; i++) {
      at = StructSp[2].ia[i];
      to = StructSp[2].ia[i+1];
      at--;  to--;  /* for arrays starting at 0 index */
      c_val[off+i+1] += l*StructSp[2].a[at];
      for(jj=at+1; jj<to; jj++)
        CList(off+i+1, off+StructSp[2].ja[jj], l*StructSp[2].a[jj]);
    }
  }

  /* ADD VAR RATIO/matrix OF ALL (if k=1) OR OTHER (if k=2) RANDOM FACTORS */
  for(j=k; j<=NumRandom; j++) {
    off = CumLevels[NumFixed+j-1]; if(FI) off++; if(TM) off += nxi;
    if(DCM[j]) {
      for(i=1; i<=Levels[NumFixed+j]; i++)
        c_val[off+i] += rl[j];
    }
    else {
      printf("Adding supplied cov mat - fvr = %f, j = %d - Levels[%d] = %d\n\n", rl[j], j, (NumFixed+j), Levels[NumFixed+j]);
      for(i=0; i<Levels[NumFixed+j]; i++) {
        at = StructSp[(1+j)].ia[i];
        to = StructSp[(1+j)].ia[i+1];
        at--;  to--;  /* for arrays starting at 0 index */
        c_val[off+i+1] += rl[j]*StructSp[(1+j)].a[at];
        for(jj=at+1; jj<to; jj++)
          CList(off+i+1, off+StructSp[(1+j)].ja[jj], rl[j]*StructSp[(j+1)].a[jj]);
      }
    }
  }

  // Do MQTL inverse matrices ...
  if(GRM) {
    for(j=1; j<=NumLoci; j++) {
      off = CumLevels[NumFR+j-1]; if(FI) off++; if(TM) off += nxi;
      for(i=0; i<Levels[NumFR+j]; i++) {
        at = StructSp[(2+NumRandom+j)].ia[i];
        to = StructSp[(2+NumRandom+j)].ia[i+1];
        at--;  to--;  /* for arrays starting at 0 index */
        c_val[off+i+1] += rl[NumRandom+j]*StructSp[(2+NumRandom+j)].a[at];
        for(jj=at+1; jj<to; jj++)
          CList(off+i+1, off+StructSp[(2+NumRandom+j)].ja[jj], rl[NumRandom+j]*StructSp[(2+NumRandom+j)].a[jj]);
      }
    }
  }

  free_dvector(rl, 1, (NumRandom+NumLoci));
  rl = NULL;
}
