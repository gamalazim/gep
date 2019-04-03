#include "gt.h"
#include "analyze.h"
#include "debug.h"

void OutTMatrices() {
  int i, j;
  FILE *fp;

  if((fp=fopen("TMatrices","w"))==NULL)
    Error("Can't open TMatrices File");

  fprintf(fp,"Q and vee ...\n");
  for(i=1; i<=nxi; i++) {
    fprintf(fp,"%3.d ", i);
    for(j=1; j<=nxi; j++)
      fprintf(fp,"%16.8g ", Q[i][j]);
    fprintf(fp,"      %16.8g\n", vee[i]);
  }

  fprintf(fp,"\nL, R, and ee ...\n");
  for(i=1; i<=NumRecords; i++) {
    fprintf(fp,"%3.d ", i);
    for(j=1; j<=nxi; j++)
      fprintf(fp,"%16.8g ", L[i][j]);
    fprintf(fp,"      %16.8g ", vR[i]);
    fprintf(fp,"      %16.8g\n", ee[i]);
  }

  fclose(fp);
}

void Initiate(int opt) {
  /*
   * To initiate the iterative algorithm, we need to specify initial values
   * for the solutions. For the boundary points, we use (Nj - N1)/N,
   * where Nj represents the total number of individuals in categories 1 to j,
   * N1 is the number of individuals in category 1, and N = NumRecords.
   * We use -N1/N for mu, and 0's for the rest of the solutions.
   */
  int i;

  if(opt==1) {
    int *catcount;
    catcount = ivector(1, NCATS);

    for(i=1; i<=NCATS; i++) catcount[i] = 0;
    for(i=1; i<=NumRecords; i++) catcount[((int) y[i])]++; /* count categories */
    for(i=1; i<NCATS; i++) catcount[i+1] += catcount[i];
    for(i=1; i<=SumLevels; i++) b[i-1] = 0;

    for(i=2; i<NCATS; i++)
      b[i-2] = ((double) (catcount[i]-catcount[1]))/((double) NumRecords);
    b[NCATS - 2] = -((double) catcount[1])/((double) NumRecords); /* initial mu */

    free_ivector(catcount, 1, NCATS);
    catcount = NULL;
  }
  else {
    FILE *fp;
    if((fp=fopen("SolStart","r"))==NULL)
      Error("Initiate(opt = 2): Can't open SolStart for reading solutions");
    for(i=1; i<=SumLevels; i++) fscanf(fp,"%lf", &b[i-1]);
    fclose(fp);
  }

  for(i=0; i<=SumLevels; i++) oldb[i] = b[i];
}

// Wb =  [X Z qtl cvrt]*b
// ======================
void MakeWb(double *Wb) {
  int i, j, k, cc;
  for(i=1; i<=NumRecords; i++) {
    Wb[i] = b[nxi]; /* mu -> the vector of 1's in X */
    for(j=1; j<=NumFR; j++) {
      if(FactorsArray[i][j]) {
	k=nxi+CumLevels[j-1]+FactorsArray[i][j]; /*ignored -1 to account for mu*/
	Wb[i] += b[k]; /* Add fixed and random solutions
			  to produce (xa + zu) */
      }
    }

    if(GRM) {
      for(j=NumFR+1; j<=NumFactors; j++) {
        k=nxi+CumLevels[j-1]+FactorsArray[i][NumFR+1]; /* ignored -1 to account for mu */
        Wb[i] += b[k];
        Wb[i] += b[k+1];
      }
    }

    if(Cov) for(cc=1; cc<=Cov; cc++) {
	k = nxi+CumLevels[NumFactors]+cc;
	Wb[i] += (cvrt[i][cc]*b[k]);
    }
  }
}

double phi(double *mu, double *xi, int i, int k, int flag) {
  /*
   * (flag = 1): return the NDF(i, k)
   * (flag = 2): return NDF(i, k-1) - NDF(i, k)
   * (flag = 3): return CDF(i, k) - CDF(i, k-1)
   *
   * k is the boundary point! if(k==0), it's -infinity;
   * if(k==1), it's 0; if(k==2), it's 2nd bp; ...;
   * if(k==NCATS), it's +infinity;
   */
  switch(flag) {
  case 1:
    if(k == 0 || k == NCATS) return 0;
    else return (ndf((xi[k]-mu[i])));

  case 2:
    if(k<=0) Error("< -infinity value when computing delta");
    else if(k == 1)
      return (-(ndf(-mu[i])));
    else if(k==NCATS) return (ndf(xi[k-1]-mu[i]));
    else return (ndf(xi[k-1]-mu[i]) - ndf(xi[k]-mu[i]));

  case 3:
    if(k<=0) Error("< -infinity value when computing DELTA");
    else if(k == 1)
      return (normal_cd(-mu[i]));
    else if(k==NCATS) return (1.0 - normal_cd(xi[k-1]-mu[i]));
    else return (normal_cd(xi[k]-mu[i]) - normal_cd(xi[k-1]-mu[i]));

  default:
    Error("unknown flag in phi()");
    return 0; // not needed but just to avoid compiler warning
  }
}

void TMatrices(double *mu, double *xi) {
  int i, j, k, yi;
  double *delta, *DELTA, *fi, *lv, *rv;

  /* Allocate local arrays */
  delta = dvector(1, NCATS);
  DELTA = dvector(1, NCATS);
  if(nxi) {
    fi = dvector(0, NCATS);
    lv = dvector(1, nxi);
    rv = dvector(1, nxi);
    for(i=1; i<=nxi; i++) lv[i] = rv[i] = Q[i][i] = Q[i][i+1] = 0.0;
  }

  for(i=1; i<=NumRecords; i++) {
    vR[i] = 0.0;
    for(j=1; j<=NCATS; j++) { /* make delta, DELTA, vR */
      delta[j] = phi(mu, xi, i, j, 2);
      DELTA[j] = phi(mu, xi, i, j, 3);
      if(DELTA[j] < D_min) DELTA[j] = D_min;
      vR[i] += square(delta[j])/DELTA[j];
    }
    yi = y[i];
    ee[i] = delta[yi]/DELTA[yi];
    if(nxi) {
      for(j=0; j<=NCATS; j++) fi[j] = phi(mu, xi, i, j, 1);
      for(k=2; k<NCATS; k++)
	L[i][k-1]= fi[k]*(delta[k]/DELTA[k]-delta[k+1]/DELTA[k+1]); /* L[i,] */

      if(yi>1 && yi<NCATS)
	lv[yi-1] += (fi[yi]/DELTA[yi]);
      if(yi>2)
	rv[yi-2] += (fi[yi-1]/DELTA[yi]);

      for(j=2; j<NCATS; j++) { /* make UT of tridiagonal Q */
	Q[j-1][j-1] += ( square(fi[j])*(1.0/DELTA[j]+1.0/DELTA[j+1]) );
	Q[j-1][j] -= ( fi[j]*fi[j+1]/DELTA[j+1] ); /* UT element */
      }
    }
  }
  for(i=1; i<=nxi; i++) vee[i] = lv[i] - rv[i];

  free_dvector(delta, 1, NCATS);
  free_dvector(DELTA, 1, NCATS);
  delta = DELTA = NULL;
  if(nxi) {
    free_dvector(fi, 0, NCATS);
    free_dvector(lv, 1, nxi);
    free_dvector(rv, 1, nxi);
    fi = lv = rv = NULL;
  }
}

void DInvBeta(double *p) {
  int i, k, kk, j, jj, at, to, off;
  int order;
  int c = SumFixedLevels;
  double li; /* lambda inverse */
  double *rli;

  if(NumRandom == 0 && NumLoci == 0) return;

  rli = dvector(1, (NumRandom+NumLoci));
  for(i=1; i<=(NumRandom+NumLoci); i++) rli[i] = p[0]/p[i];

  if(TM) c += (1 + nxi); /* add threshold and overall mean eqn's */

  /* Diagonal elements of Inv(D) */
  kk = 1;
  if(RSM) {
    kk = 2;
    li = rli[1];
    order = StructSp[2].n; /* Order of NRM or MGS matrix. order = nind, any could be used */
    for(i=1; i<=order; i++) {
      k = StructSp[2].ia[i-1]; /* UT row-list start */
      R[c+i] -= (li*StructSp[2].a[k-1]*b[c+i-1]);
      /* b is [0, ] and R is [1, ] - StructSp[2].a[0, nze-1]*/
    }
  }
  else {
    for(j=kk; j<=NumRandom; j++) {
      off = CumLevels[NumFixed+j-1];
      if(DCM[j]) { /* Diagonal elements of Inv(D)i, i=kk to NumRandom */
        for(i=1; i<=Levels[NumFixed+j]; i++) {
          R[off+i] -= (rli[j]*b[off+i-1]);
          /* b is [0, ] and R is [1, ] - StructSp[2].a[0, nze-1]*/
        }
      }
      else {  // RHS the whole D_inv not only its diagonal, when the inverse is supplied
        for(i=1; i<StructSp[(1+j)].n; i++) {
          k = StructSp[(1+j)].ia[i-1] + 1; /* Second to the diagonal element */
          to = StructSp[(1+j)].ia[i] - 1; /* last element of the (i-1)th list */
          for(jj=k; jj<=to; jj++) {
            R[off+i] -= (rli[j]*StructSp[(1+j)].a[jj-1]*b[off+StructSp[(1+j)].ja[jj-1]-1]);
            R[off+StructSp[(1+j)].ja[jj-1]] -= (rli[j]*StructSp[(1+j)].a[jj-1]*b[off+i-1]); /* LT */
          }
        }
      }
    }
  }

  /* Off diagonal elements and the LOWER TRIANGULAR elements */
  if(RSM) {
    for(i=1; i<order; i++) {
      k = StructSp[2].ia[i-1] + 1; /* Second to the diagonal element */
      to = StructSp[2].ia[i] - 1; /* last element of the (i-1)th list */
      for(j=k; j<=to; j++) {
	R[c+i] -= (li*StructSp[2].a[j-1]*b[c+StructSp[2].ja[j-1]-1]);
	R[c+StructSp[2].ja[j-1]] -= (li*StructSp[2].a[j-1]*b[c+i-1]); /* LT */
      }
    }
  }

  // Do MQTL inverse matrices ...
  if(GRM) {
    for(j=1; j<=NumLoci; j++) {
      off = CumLevels[NumFR+j-1]; if(TM) off += (1 + nxi);
      for(i=0; i<Levels[NumFR+j]; i++) {
        at = StructSp[(2+NumRandom+j)].ia[i];
        to = StructSp[(2+NumRandom+j)].ia[i+1];
        at--;  to--;  /* for arrays starting at 0 index */

        for(jj=at; jj<to; jj++) {
          R[off+i+1] -= rli[NumRandom+j]*StructSp[(2+NumRandom+j)].a[jj]*b[off + StructSp[(2+NumRandom+j)].ja[jj]-1];
          R[off + StructSp[(2+NumRandom+j)].ja[jj]] -= rli[NumRandom+j]*StructSp[(2+NumRandom+j)].a[jj]*b[off+i]; // no -1, i starts at 0
        }
      }
    }
  }

  if(rli) { free_dvector(rli, 1, (NumRandom+NumLoci)); rli = NULL; }
}

void TM_SetUp(double *p) {
  /*
   * Set up threshold model equations for the first time using initial
   * guess of the solutions if tmrnd == 1 and using the most recent set
   * of solutions if tmrnd > 1.
   */
  int i, j, k, m, ii, cc;
  double *Wb, *xi;
  double *vec2alloc; /* a vector to allocate (e.g., for the Mu equation) */

  vec2alloc = dvector(nxi+1, SumLevels); // SumLevels = nxi + 1 + NumFactors + Cov
  for(i=1; i<=SumLevels; i++) R[i] = 0.0;
  for(i=1; i<=SumLevels; i++) { c_val[i] = c_link[i]=0; c_col[i]= i; }
  CLast = SumLevels+1;

/* THRESHOLD MODEL EQUATIONS:

|<-- Coefficient Matrix --> |   |<---- RHS ----->|
_                         _     _              _
|*------------------------- |   |                |
|  *  Q  ||    L'[X  Z]    ||   |                |
|    *   ||     = L'W      ||   |      vee       |
|      * ||                ||   |                |
|        * ---------------- |   |----------------|
|          *                |   |                |
|            *    W'RW      |   | W'ee - Inv(D)B |
|              *            |   |                |
|                *          |   |                |
|                  *        |   |                |

*/

  /******************************
  * Build Q, L, vee, vR, and ee
  *****************************/
  xi = dvector(1, NCATS-1); xi[1] = 0; /* 1st bp = 0 */
  for(i=2; i<NCATS; i++) xi[i] = b[i-2]; /* Top nxi soln's are for bp's */
  Wb = dvector(1, NumRecords);
  MakeWb(Wb); /* =  [X Z qtl cvrt]*b  */
  TMatrices(Wb, xi);

  if(1/*DBUG*/) {
    //RSolutions("dbugSolFile");
    //OutWb(Wb);
    //OutTMatrices();
  }

  /*
  * Build UT of [Q  L'W] = [Q  L'[X Z cvrt]] = [Q  L'X  L'Z  L'cvrt]
  * Also build the RHS
  */
  if(nxi) { /* if > 2 categories */
    /* -- Q -- */
    for(i=1; i<=nxi; i++) {
      R[i] = vee[i];
      c_val[i] = Q[i][i];
      for(j=(i+1); j<=nxi; j++) Hock2CList(i, j, Q[i][j], 1);
    }

     /* -- L'W -- */
    for(i=1; i<=nxi; i++) {
      for(ii=(nxi+1); ii<=SumLevels; ii++) vec2alloc[ii]=0;
      for(j=1; j<=NumRecords; j++) {
	vec2alloc[nxi+1] += L[j][i]; /* Overall mean */
	if(Cov) {
	    for(cc=1; cc<=Cov; cc++) {
		m = nxi+1+CumLevels[NumFactors]+cc;
		vec2alloc[m] += (L[j][i]*cvrt[j][cc]); /* Covariates */
	    }
	}
	for(k=1; k<=NumFR; k++)
	  if(FactorsArray[j][k]) {
	    m = nxi+1+CumLevels[k-1]+FactorsArray[j][k];
	    vec2alloc[m] += L[j][i];
	  }

        if(GRM) {
          for(k=NumFR+1; k<=NumFactors; k++) {
            m = nxi+1+CumLevels[k-1]+FactorsArray[j][NumFR+1];
            vec2alloc[m] += L[j][i];
            vec2alloc[m+1] += L[j][i];
          }
        }

      }
      InsertAList(i, vec2alloc, (nxi+1), SumLevels);
    }
  } /* if(nxi) */

  /*
  * Build UT of [W'RW] = | X'RX  X'RZ |
  *                      | Z'RX  Z'RZ |
  * & Build  | X'ee |
  *          | Z'ee |
  * ---------------------------------------------------------- */

  // The Overall Mean Equation: 1'R[1 X Z qtl cvrt] = 1'ee
  // ------------------------------------------------------
  k=nxi+1;
  for(ii=(nxi+1); ii<=SumLevels; ii++) vec2alloc[ii]=0;

  for(i=1; i<=NumRecords; i++) {

      // Diagonal mu x mu [ 1' R 1 ]
    c_val[k] += vR[i];
      // vec2alloc[k] += vR[i];

      // [1' R W]
    for(j=1; j<=NumFR; j++)
      if(FactorsArray[i][j]) {
        m = k + CumLevels[j-1]+FactorsArray[i][j];
        vec2alloc[m] += vR[i];
      }

    // [1' R qtl]
    if(GRM)
      for(j=NumFR+1; j<=NumFactors; j++) {
        m = k+CumLevels[j-1]+FactorsArray[i][NumFR+1];
        vec2alloc[m] += vR[i];
        vec2alloc[m+1] += vR[i];
      }


      // [1' R cvrt]
    if(Cov)
      for(cc=1; cc<=Cov; cc++) {
        m = k+CumLevels[NumFactors]+cc;
        vec2alloc[m] += (vR[i]*cvrt[i][cc]);
      }

      // RHS 1'ee
    R[k] += ee[i];
  }

  InsertAList(k, vec2alloc, k+1, SumLevels);  /* Insert Mu Eqn except the diagonal element */
  if(vec2alloc) { free_dvector(vec2alloc, nxi+1, SumLevels); vec2alloc = NULL; }

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

  /* ----------------------------
  * Subtract D-1*B from the RHS */
  if(NumRandom | NumLoci) DInvBeta(p);

  /* --------------------------
    - Do the covariates part -
    -------------------------- */
  if(Cov) DoCov();
  if(NumRandom | NumLoci) SpC2(p);  /* Add var inv matrices to the LSC linked lists */

  cFree(0); /* free previous allocations of StructSp[ 0 ] -- cFree is conditional */
  sform(0, SumLevels);

  if(DBUG) {
    Struct2Default(0); /* OutPut C */
    OutRHS();
    OutSolVector();
  }

  if(Wb) { free_dvector(Wb, 1, NumRecords); Wb = NULL; }
  if(xi) { free_dvector(xi, 1, nxi); xi = NULL; }

  //Error("Fake error ...");
} /* TM_SetUp */


void TM_Prep() {

  SM = 1;
  vccc = 100.0; solcc = oldsolcc = 100.0;  /* No conversion yet */
  nxi = NCATS-2;

//  if(RSM) {
//    grm(); /* Create Lists */
//    mgs();
//    sform(2, nind); /* sparse form of the relationship matrix */
//  }

  /*
   * There are NCATS - 2 unknown boundary points (bp):
   * (-inf)---(0)---(1st bp)---(2nd bp)---(...)---(+inf)
   * We estimate those between 0 and +inf
   */
  SumLevels += (NCATS - 2);
  // SumLevels++;  already done in MMParms() /* for the overall mean */

  /* Allocate memory for threshold matrices */
  Q = dmatrix(1, nxi, 1, nxi);
  L = dmatrix(1, NumRecords, 1, nxi);
  vR = dvector(1, NumRecords);
  ee = dvector(1, NumRecords);
  vee = dvector(1, nxi);
  oldb = dvector(0, SumLevels);

  Initiate(SolStart); /* solutions' first guess
			 Two options (1/2): 1: for initial guess. 2: for using
			 existing solutions on disk */
}

int violate() {
  int i;

  /* if(0 > b[0]+oldb[0]) return 1; */

  for(i=0; i<(nxi-1); i++)
    if( (b[i]+oldb[i]) > (b[i+1]+oldb[i+1]) )
      return 1;

  return 0;
}

void shrink(double ShrinkFactor) {
  int i;

  for(i=0; i<SumLevels; i++)
    b[i] *= ShrinkFactor;
}

void SimpleStat(double *data, int from, int to, double *ss) {
  /*
   * Compute simple stats for elements 'from' to 'to' of the vector 'data'
   * Place computed statistics in ss
   * 5/22/02 ss is of size 2 and contains mean and variance
   */
  int i, n = (to-from+1);
  double sum1, sum2;

  sum1 = sum2 = 0.;
  for(i=from; i<=to; i++) {
    sum1 += data[i];
    sum2 += square(data[i]);
  }

  ss[0] = sum1/n;  /* Mean */
  ss[1] = (sum2 - (square(sum1)/n))/(n-1); /* Variance */
}


void Burst(double *ss) {
  int i, k=0;
  double sd = sqrt(ss[1]); /* Stand Dev. of solutions */

  for(i=0; i<SumLevels; i++)
    if(fabs(b[i]) > 10*sd) {
      b[i] = 0;
      k++;
    }
  if(k) printf("Warning ... %d Solutions became > 10*SD and were 0'ed\n", k);
}

void StopCritDebug() {

  int i, k, from , to;
  double num = 0, denum = 0;

  from = 1 + nxi + SumFixedLevels;
  to = CumLevels[NumFactors];
  num = 0, denum = 0;
  for(i=from; i<to; i++) {
    num += (b[i]*b[i]);
    denum += ((oldb[i]+b[i])*(oldb[i]+b[i]));
  }
  printf("ALL: %lf\n", (sqrt(num/denum)));

  for(k=(NumFixed+1); k<=NumFactors; k++) {
    from = CumLevels[k-1];
    to = CumLevels[k];

    num=denum=0.;
    for(i=from; i<to; i++) {
      num += (b[i]*b[i]);
      denum += ((oldb[i]+b[i])*(oldb[i]+b[i]));
    }
    printf("Factor(%d): %lf\n", k, (sqrt(num/denum)));
  }
}



double Get_SolStopCrit(char *opt) {

  /*
     b[0 : SumLevels-1]: difference between two subseq iteratives.
     oldb[0 : SumLevels-1]: old set of solutions
  */

  int i=0, from=0, to=0;
  double num = 0., denum = 0.;

  if(!strcmp(opt, "all")) {
    from = 1 + nxi + SumFixedLevels;
    to = CumLevels[NumFactors];
  }
  else if(!strcmp(opt, "one")) {
    printf("Factor Used in Calculating Stopping Criterion: %d\n", StopFactor);
    from = CumLevels[StopFactor-1];
    to = CumLevels[StopFactor];
  }

  else
    Error("Unknown Option to Get_SolStopCrit()");

  num = 0, denum = 0;
  for(i=from; i<to; i++) {
    num += (b[i]*b[i]);
    denum += ((oldb[i]+b[i])*(oldb[i]+b[i]));
  }
  return(sqrt(num/denum));
}

void TM_Solve(double *p) {
  int i;
  double sf= 0.1 /* shrink factor */;
  double big;
  int ibig;
  double ss[2]; /* maen and variance */

  for(i=0; i<SumLevels; i++)
    oldb[i] = b[i];

  DenseVectorSolve(p);

  SimpleStat(b, 0, SumLevels-1, ss); /* mean = ss[0]; var = ss[1] of the diff vector */

  big = fabs(b[0]); ibig = 0;
  for(i=1; i<SumLevels; i++)
    if(fabs(b[i]) > big) { big = fabs(b[i]); ibig = i; }

  ibig++;
  printf("FABS(Max(b): %d.  %g\n",ibig, big);
  printf("RHS[%d]:  %g\n", ibig, R[ibig]);
  printf("Mean:  %g -- Variance:  %g -- SD:  %g\n", ss[0], ss[1], sqrt(ss[1]));

  while(violate()) {
    printf("Shrinking Components of b\n");
    shrink(sf); /* shrink if 0 < t2 < ... < t(nxi+1) */
    sf /= 10;
  }

  Burst(ss);  /* delete bloated solutions */

  /* Compute solutions stopping criterion, Misztal et al. 1989 - JDS 72:1557 */
  oldsolcc = solcc;
  solcc = Get_SolStopCrit("all"); /* use random solutions to get stopping criterion */

  StopCritDebug();  // Debugging Tool ...

  shrink(.8); // shrink anyway ..
  for(i=0; i<SumLevels; i++)
    b[i] += oldb[i];

  printf("\n\n%d. SOLUTIONS STOPPING CRITERION: %g\n",solit, solcc);
  printf("A Formal Criterion = 10 to the negative p corresponds to p sig decimal degits\n");
  for(i=0; i<nxi; i++) printf("(t%d. %g    %g)\n",i+1, b[i], (b[i]-oldb[i]));
  printf("\n");
  printf("Overall Mean = %g\n", b[nxi]);
  printf("\n\n");
  printf("TM_Solve ... DONE ...\n"); timedate();
}
