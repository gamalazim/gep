
// Gamal Abdel-Azim 1/29/2008
// Kinships and their inverses ...

#include "kinshipSetUp.h"
#include "gt.h"
#include "debug.h"

// createLists() is a function that constructs linked lists
void createLists(int *sire, int *dam, int n) {
  int i, j, tn, l, minp, maxp;
  double chk;

  col = ivector(1,maxSz);
  links = ivector(1,maxSz);
  icol = ivector(1,ImaxSz);
  ilinks = ivector(1,ImaxSz);

  for(i=1; i<=n; i++) links[i] = ilinks[i] = col[i] = icol[i] = 0;
  GLast = IGLast = n+1;

  for(i=1; i<=n; i++)  /* add the sire by dam element */
    List(sire[i], dam[i]);

  for(i=n-1; i>0; i--) {
    chk = ((double) GLast)/((double) maxSz);
    if(chk > 0.94) Error("Increase maxSz");
    l = i;
    while(l > 0) { /* while not a terminal link */
      List(sire[i], col[l]);
      List(dam[i], col[l]);
      l = links[l]; }
  }

  /* --------------------------
      The inverse linked lists
     -------------------------- */
  for(i=1;i<=n;i++) {
    chk = ((double) IGLast)/((double) ImaxSz);
    if(chk > 0.94) Error("Increase ImaxSz");
    maxp=Max(sire[i],dam[i]); minp=Min(sire[i],dam[i]);
    IList(maxp,minp); IList(i,sire[i]); IList(i,dam[i]);
  }
}

// Free arrays allocated in createLists()
void freeLists() {

  if(col) free_ivector(col,1,maxSz);
  if(links) free_ivector(links, 1,maxSz);
  if(icol) free_ivector(icol,1,ImaxSz);
  if(ilinks) free_ivector(ilinks, 1,ImaxSz);

  col = links = icol = ilinks = NULL;
}

void fillB(int locus, int *sire, int *dam, int *gtypes, int **B, int n) {
  int i, r, c, from = (locus-1)*(2*n), to1 = ((locus-1)*(2*n) + n - 1), to2 = locus*(2*n)-1;

  i = from; r = 1; c = 1;
  while(i <= to1) {
    B[r][c] = gtypes[i];
    i++;
    r++;
  }

  i = to1 + 1; r = 1; c = 2;
  while(i <= to2) {
    B[r][c] = gtypes[i];
    i++;
    r++;
  }
}

// Inverse of Gametic Relationship Matrix -- LOCUS number 'locus'
void grm(int locus, int *sire, int *dam, int *gtypes, int n) {

  int i, j, tn, l, k=0, minp, maxp, **B;
  double **Q;

  f = dvector(0,n); f[0]=0.;
  Q = dmatrix(1,2,1,4);
  values = dmatrix(1, GLast, 1, 4);
  B = imatrix(1,n,1,2);

  for(i=1; i<GLast; i++) values[i][1]=values[i][2]=values[i][3]=values[i][4]=0.;
  fillB(locus, sire, dam, gtypes, B, n);

  /* ---------------------------------------
     start computations and inversion steps
     --------------------------------------- */
  for(i=1;i<=n;i++) {
    if(sire[i]==0 && dam[i]==0) {
      f[i]=0.;
      MakeQ(Q, i,sire[i],dam[i],B);
      Add2G(i, sire[i], dam[i], Q);
      k++;
    }
    else break;
  }

  /* for(i=1; i<=n; i++) SortIList(i); */

  iv3 = dmatrix(1,n,1,3);
  iv4 = dmatrix(1,IGLast,1,4);

  for(i=1;i<=n;i++) iv3[i][1]=iv3[i][2]=iv3[i][3]=0.;
  for(i=1;i<=k;i++) iv3[i][1]=iv3[i][3]=1.;
  for(i=1;i<=IGLast;i++) iv4[i][1]=iv4[i][2]=iv4[i][3]=iv4[i][4]=0.;

  for(i=k+1; i<=n; i++) {
    MakeQ(Q, i,sire[i],dam[i],B);
    Add2Inverse(i, sire[i], dam[i], Q);
    Add2G(i, sire[i], dam[i], Q);
  }

  /* PrintIList(n); */
  MakeUT(n);

  free_imatrix(B,1,n,1,2);
  free_dmatrix(iv3,1,n, 1,3);
  free_dmatrix(iv4,1,IGLast, 1,4);
  free_dvector(f,0,n);
  free_dmatrix(Q,1,2,1,4);
  free_dmatrix(values,1,GLast,1,4);

  B = NULL;
  iv3 = iv4 = Q = values = NULL;
  f = NULL;
}


void List(int x, int y) {
  int minp, maxp, j;
  if(x == 0 || y == 0) return;
  if(x == y) return;
  minp = Min(x, y); maxp = Max(x, y);
  if(col[maxp] == 0) { col[maxp] = minp; return; }
  while(1) {
    if(col[maxp]==minp) return;
    if(links[maxp] == 0) {
      links[maxp]=GLast; links[GLast]=0;
      col[GLast]=minp;
      GLast++; return; }
    maxp = links[maxp];
  }
  return;
}

void IList(int maxp, int minp) {
  int j;
  if(maxp == 0 || minp == 0) return;
  if(maxp == minp) return;
  if(icol[maxp] == 0) { icol[maxp] = minp; return; }
  while(1) {
    if(icol[maxp]==minp) return;

    if(icol[maxp] > minp) {
      ilinks[IGLast]=ilinks[maxp]; ilinks[maxp]=IGLast;
      icol[IGLast]=icol[maxp]; icol[maxp]=minp;
      IGLast++; return; }

    if(ilinks[maxp] == 0) {
      ilinks[maxp]=IGLast; ilinks[IGLast]=0;
      icol[IGLast]=minp;
      IGLast++; return; }

    maxp = ilinks[maxp];
  }
  return;
}


void SortIList(int i) {
  int j;
  while(ilinks[i] > 0) {
    j = ilinks[i];
    while(j > 0) {
      if(icol[i] > icol[j])
	swap(&icol[j], &icol[i]);
      j = ilinks[j]; }
    i = ilinks[i]; }
}

void PrintList(int n) {
  int i, j;
  FILE *fp;

  if( (fp= fopen("list", "w"))==NULL)
    Error("can't open a list file");

  for(i=1; i<GLast; i++) {
    fprintf(fp,"%3d  %3d  ", col[i], links[i]);
    for(j=1;j<=4;j++) fprintf(fp,"%8.5f ", values[i][j]);
    if(i<=n) fprintf(fp, "%8.5f   ", f[i]);
    fprintf(fp, "\n"); }

  fclose(fp);
}

void PrintIList(int n) {
  int i, j;
  FILE *fp;

  if( (fp= fopen("ilist", "w"))==NULL)
    Error("can't open an ilist file");

  for(i=1; i<IGLast; i++) {
    fprintf(fp,"%3d  %3d  ", icol[i], ilinks[i]);
    for(j=1;j<=4;j++) fprintf(fp,"%f ", iv4[i][j]);
    if(i<=n) for(j=1;j<=3;j++) fprintf(fp,"%f ", iv3[i][j]);
    fprintf(fp, "\n"); }
  fclose(fp);
  return;
}

void Add2Inverse(int in, int s, int d, double **Q) {
  int i, j, k, l, p,q,i1,i2,j1,j2,TO[4], dex;
  double **delta, **T, **P, **g, *v, c, rootm, komc;

  delta=dmatrix(1,2,1,6); T=dmatrix(1,2,1,2); P=dmatrix(1,4,1,4); g=dmatrix(1,2,1,2);

  v = IKinship(s,d);
  TO[0]=0; TO[1]=s; TO[2]=d; TO[3]=in;

  for(i=1; i<=2; i++) for(j=1; j<=2; j++) {
    T[i][j] = 0.;
    if(Q[1][1]+Q[1][2] != 0.) T[i][j] += (Q[1][i]*Q[2][2+j])/(Q[1][1]+Q[1][2]);
    if(Q[1][3]+Q[1][4] != 0.) T[i][j] += (Q[2][i]*Q[1][2+j])/(Q[1][3]+Q[1][4]); }

  f[in] = T[1][1]*v[1]+T[1][2]*v[2]+T[2][1]*v[3]+T[2][2]*v[4];

  /* for(i=1;i<=4;i++) P[i][i] = 1.; */
  if(s) P[1][1]=P[2][2]=1.; else P[1][1]=P[2][2]=0.;
  if(d) P[3][3]=P[4][4]=1.; else P[3][3]=P[4][4]=0.;

  P[1][2]=P[2][1]=f[s]; P[3][4]=P[4][3]=f[d]; P[1][3]=P[3][1]=v[1];
  P[1][4]=P[4][1]=v[2]; P[2][3]=P[3][2]=v[3]; P[2][4]=P[4][2]=v[4];
/*
    printf("in = %d\n",in);
    printf("%f %f %f %f\n%f %f %f %f\n",Q[1][1],Q[1][2],Q[1][3],Q[1][4],
	   Q[2][1],Q[2][2],Q[2][3],Q[2][4]);
*/
  T[1][1]=T[1][2]=T[2][1]=T[2][2]=0.;
  for(k=1;k<=2;k++) for(l=1;l<=2;l++) for(i=1;i<=4;i++) for(j=1;j<=4;j++)
    T[l][k] += Q[l][j]*P[j][i]*Q[k][i]; /* T = Q*Csd*t(Q) */

  g[1][1]=1.-T[1][1]; g[1][2]=f[in]-T[1][2]; g[2][1]=f[in]-T[2][1]; g[2][2]=1.-T[2][2];

  /* Check for singularity (usually with incomplete marker data) */
  if(g[1][1]*g[2][2]-g[1][2]*g[2][1] == 0.0) {
    /* printf("Add2Inverse: ZERO DETERMINANT was encountered with individual %d\n",in); */
    /* exit(1); */ return;}

  c = sqrt((g[2][2] - (g[2][1]*g[2][1]/g[1][1])));
  rootm = sqrt(g[1][1]);

  if(c != 0) {
    komc = g[2][1]/(g[1][1]*c);
    delta[1][5]=1./rootm; delta[1][6]=0.;
    delta[2][5]= -komc;
    delta[2][6]= 1./c;

    for(j=1;j<=4;j++) {
      delta[1][j]=-Q[1][j]/rootm;
      delta[2][j]=komc*Q[1][j]-Q[2][j]/c; }

    /* for(i=1;i<=6;i++) for(j=i;j<=6;j++) {
       c = delta[1][i]*delta[1][j] + delta[2][i]*delta[2][j]; */
    /* inverse[to[i]][to[j]] += c; */
    /* fprintf(fpi,"%d %d %f\n", to[i], to[j], c); */
    /* } */

    /* ----------------------------
       Working With limited memory
       ---------------------------- */
    for(i=1;i<=3;i++)
      for(j=1;j<=i;j++) {
	if(i == j && TO[i]) { /* if diagonal AND known */
	  q=2*i; p=q-1;
	  iv3[TO[i]][1] += delta[1][p]*delta[1][p]+delta[2][p]*delta[2][p];
	  iv3[TO[i]][3] += delta[1][q]*delta[1][q]+delta[2][q]*delta[2][q];
	  iv3[TO[i]][2] += delta[1][q]*delta[1][p]+delta[2][q]*delta[2][p]; }
	else {
	  i2=2*i; i1=i2-1; j2=2*j; j1=j2-1; p=2; q=3;
	  if(i == 2) {
	    if(s > d) { p=3; q=2; dex=GetDex(TO[j],TO[i]); }
	    else dex = GetDex(TO[i], TO[j]); }
	  else dex = GetDex(TO[i], TO[j]);

	  if(dex) {
	    iv4[dex][1] += delta[1][i1]*delta[1][j1]+delta[2][i1]*delta[2][j1];
	    iv4[dex][p] += delta[1][i1]*delta[1][j2]+delta[2][i1]*delta[2][j2];
	    iv4[dex][q] += delta[1][i2]*delta[1][j1]+delta[2][i2]*delta[2][j1];
	    iv4[dex][4] += delta[1][i2]*delta[1][j2]+delta[2][i2]*delta[2][j2];
	  }
	}
      }
  } /* if(c != 0) */

  else  /* if chol isn't working, direct-invert D (g) */
    UseWang(Q, g, TO);

}

void UseWang(double **Q, double **D, int *TO) {

  int i, j, k, dex, maxp, minp;
  double **QI2; /* [-Qi I2]' */
  double det = det = (D[1][1]*D[2][2]) - (D[1][2]*D[2][1]);
  double **W;

  W = dmatrix(1,6,1,6);
  QI2 = dmatrix(1,6,1,2);

  /* Invert D in D */
  dswap(&D[1][1], &D[2][2]);
  D[1][1] /= det; D[2][2] /= det;
  D[1][2] /= (-det); D[2][1] /= (-det);

  /* make [-Q  I2]' */
  for(j=1; j<=4; j++) for(k=1; k<=2; k++) QI2[j][k] = -Q[k][j];
  QI2[5][1] = QI2[6][2] = 1.0; QI2[5][2] = QI2[6][1] = 0.0;

  /* Make W */
  NMNp(QI2, D, 6, 2, W);

  /* Add diagonal blocks */
  for(i=1; i<=3; i++) {
    if(TO[i]) {
      iv3[TO[i]][1] += W[2*i-1][2*i-1]; /* diagonal elem't      */
      iv3[TO[i]][3] += W[2*i][2*i];     /* diagonal elem't      */
      iv3[TO[i]][2] += W[2*i][2*i-1];   /* off diagonal elem't  */
    }
  }

  /* Add R(i,s) and R(i,d) */
  for(i=1; i<=2; i++) {
    dex = GetDex(TO[3], TO[i]);
    if(dex) {
      iv4[dex][1] += W[5][2*i-1];
      iv4[dex][2] += W[5][2*i];
      iv4[dex][3] += W[6][2*i-1];
      iv4[dex][4] += W[6][2*i];
    }
  }

  /* Add R(s,d) or R(d,s) */
  if(TO[1] && TO[2]) {
    if(TO[1] > TO[2]) {  /* if(s > d) then add R(s,d) */
      dex = GetDex(TO[1], TO[2]);
      iv4[dex][1] += W[1][3];
      iv4[dex][2] += W[1][4];
      iv4[dex][3] += W[2][3];
      iv4[dex][4] += W[2][4];
    }
    else {
      dex = GetDex(TO[2], TO[1]);
      iv4[dex][1] += W[3][1];
      iv4[dex][2] += W[3][2];
      iv4[dex][3] += W[4][1];
      iv4[dex][4] += W[4][2];
    }
  }
  free_dmatrix(W, 1, 6, 1, 6);
  free_dmatrix(QI2, 1,6,1,2);
  W = QI2 = NULL;
}

void NMNp(double **N, double **Mm, int row, int col, double **product)
/*=================================================================
  a function to multiply N * Mm * N' -- where N is a rowxcol matrix,
  Mm is a colxcol square matrix and N' is the transpose of N
  =================================================================*/
{
  int i, j, k, l;

  /* set product = 0.0 */
  for(i=1; i<=row; i++)
    for(j=1; j<=row; j++)
      product[i][j] = 0.0;

  for(k=1; k<=row; k++)
    for(l=1;l<=row;l++)
      for(i=1; i<=col; i++)
        for(j=1; j<=col; j++)
          product[l][k] += N[l][j]*Mm[j][i]*N[k][i];
}

int GetDex(int maxp, int minp) {
  if(maxp==0 || minp==0) return 0;
  while(maxp > 0) {
    if(minp == icol[maxp]) return maxp;
    maxp = ilinks[maxp]; }
  return 0;
}

void Add2G(int i, int s, int d, double **Q) {
  int j = i, k, p;
  double *v;
  if(s==0 || d==0) return;
  while(j > 0) {
    p = s;
    for(k=0;k<2;k++) {
      v = IKinship(p, col[j]);
      values[j][1] += (Q[1][1+2*k]*v[1]+Q[1][2+2*k]*v[3]);
      values[j][2] += (Q[1][1+2*k]*v[2]+Q[1][2+2*k]*v[4]);
      values[j][3] += (Q[2][1+2*k]*v[1]+Q[2][2+2*k]*v[3]);
      values[j][4] += (Q[2][1+2*k]*v[2]+Q[2][2+2*k]*v[4]);
      p = d; }
    j = links[j]; }
}

double *IKinship(int x, int y) {
  int minp, maxp, i;
  double temp, *v; v = dvector(1, 4);

  if(x==0 || y==0) v[1]=v[2]=v[3]=v[4]=0.;
  else if(x==y) { v[1]=v[4]=1.; v[2]=v[3]=f[x]; }
  else { minp=Min(x,y); maxp=Max(x,y);
	 while(maxp > 0) {
	   if(minp == col[maxp]) {
	     for(i=1;i<=4;i++) v[i]=values[maxp][i]; break; }
	   maxp = links[maxp]; } }
  if(x<y) { temp=v[2]; v[2]=v[3]; v[3]=temp; }
  return v;
}

void MakeQ(double **Q, int i, int s, int d, int **B) {
  int l=2, c, x1, x2, k, j, g[7], o1, o2, c1, c2;
  double r = rr;
  for(k=1;k<=2;k++) for(j=1;j<=4;j++) Q[k][j] = 0.;

  if(s==0 || d==0) return;

  g[1]=B[i][1]; g[2]=B[i][2]; g[3]=B[s][1];
  g[4]=B[s][2]; g[5]=B[d][1]; g[6]=B[d][2];

  for(k=1;k<=2;k++) {
    x1 = x2 = 0;
    if(g[l]==g[5] || g[l]==g[6])
      for(j=3;j<=4;j++) if(g[k]==g[j]) x1++;
    if(g[l]==g[3] || g[l]==g[4])
      for(j=5;j<=6;j++) if(g[k]==g[j]) x2++;

    if((x1+x2) == 3){ x1 *= 2; x2*=2; }
    else { c=x1; x1+=x2; x2+=c; }

    o1=5; o2=6; c1=3; c2=4;
    c = x1; /* for j = 1, c will act in place of x1 */
    for(j=1; j<=2; j++) {
      if(g[l]==g[o1] || g[l]==g[o2])
        if(g[k]==g[c1] || g[k]==g[c2]) {
          if(g[c1]==g[c2]) Q[k][c1-2]=Q[k][c2-2]=1./c;
          else if(g[k]==g[c1]) { Q[k][c1-2]=(1-r)/c; Q[k][c2-2]=r/c; }
          else { Q[k][c1-2]=r/c; Q[k][c2-2]=(1-r)/c;} }
      swap(&o1, &c1); swap(&o2, &c2);
      c = x2; /* for j = 2, c will act in place of x2 */
    }
    l--;
  }
}

void PrintQ(double **Q) {
  int i, j;
  for(i=1;i<=2;i++) {for(j=1;j<=4;j++) printf("%f ",Q[i][j]); printf("\n");}
}

void PrintP(double **P) {
  int i, j;
  for(i=1;i<=4;i++) {for(j=1;j<=4;j++) printf("%f ",P[i][j]); printf("\n");}
}

// Transform the LT into an UT and evacuate in c_list
// Specific for grm inverse
void MakeUT(int n) {
  int h[3], p[3], *index0, i, k, l;

  index0 = ivector(1, 2*n);
  CLast = 2*n+1;
  /* Evacuate the diagonal blocks first */
  for(i=1; i<=n; i++) {
    c_link[2*i-1] = CLast; c_link[CLast] = 0.;
    c_col[2*i-1] = 2*i-1; c_col[CLast] = 2*i;
    c_val[2*i-1] = iv3[i][1]; c_val[CLast] = iv3[i][2];
    index0[2*i-1] = CLast; CLast++;

    c_link[2*i] = 0;
    c_col[2*i] = 2*i;
    c_val[2*i] = iv3[i][3];
    index0[2*i] = 2*i;
  }

  for(i=2; i<=n; i++) {
    l = i;
    while(l > 0) {
      if(icol[l] == 0) break;
      p[1] = 2*icol[l]-1; p[2] = 2*icol[l];
      h[1] = index0[p[1]]; h[2] = index0[p[2]];
      for(k=1; k<=2; k++) {
        c_link[h[k]] = CLast; c_link[CLast] = CLast+1; c_link[CLast+1] = 0;
        c_col[CLast] = 2*i-1; c_col[CLast+1] = 2*i;
        c_val[CLast] = iv4[l][0+k]; c_val[CLast+1] = iv4[l][2+k];
        index0[p[k]] = CLast+1; CLast += 2;
      }
      l = ilinks[l];
    }
  }
  free_ivector(index0, 1, 2*n);
  index0 = NULL;
}

/* Use Mewissen's method to compute D of LDL' */
void nrm(int *sire, int *dam, int n) {
  int i, j, k;
  double *kinMLD;

  kinMLD = dvector(1, n); // kinship matrix local diag
  kinMGD = dvector(1, n);
  kinMGOff = dvector(1, IGLast);

  GetD(kinMLD, sire, dam, n);

  /* Set A-1 to D-1 and prepare lists */
  for(i=1; i<=n; i++)  kinMGD[i] = 1./kinMLD[i];
  for(i=1; i<=IGLast; i++) kinMGOff[i] = 0.;

  /* Inverse computing (Henderson 76) */
  for(i=1; i<=n; i++) A_inverse(sire[i], dam[i], i, kinMGD[i]);

  UT(n); /* convert the LT diag/off into UT c_?? */

  free_dvector(kinMLD, 1, n);
  free_dvector(kinMGD, 1, n);
  free_dvector(kinMGOff, 1, IGLast);
  kinMLD = kinMGD = kinMGOff = NULL;
}

void A_inverse(int s, int d, int i, double b) {
  int dex;
  double p,q;

  if(s == 0) return;  /* s > d */
  if(d > 0) {  /* s must be greater than 0 too */
    p = .25*b;  q = -.5*b;
    kinMGD[s] += p; kinMGD[d] += p; /* Use the icol-ilink lists of G */
    dex = GetDex(s,d); kinMGOff[dex] += p;
    dex = GetDex(i,d); kinMGOff[dex] += q;
    dex = GetDex(i,s); kinMGOff[dex] += q;
    return; }

  /* the following is possible only if d==0 && s>0 */
  p = .25*b;  q = -.5*b;
  kinMGD[s] += p;
  dex = GetDex(i,s); kinMGOff[dex] += q;
}

void GetD(double *D, int *sire, int *dam, int n) {
  int i, j, k, is, id, ks, kd;
  double c, fi;
  int *ANC;
  double *Fvec, *Lvec;

  ANC = ivector(1,n);
  Fvec = dvector(0, n);  // Inb. Coeff's
  Lvec = dvector(0, n);

  for(i=0; i<=n; i++) {
    Lvec[i]=0.;
    Fvec[i]=0.;
  }
  for(i=1; i<=n; i++) {
    ANC[i]=0;
    D[i]=0.;
  }

  Fvec[0] = -1;
  for(i=1; i<=n; i++) {
    is = sire[i];
    id = dam[i];
    sire[i] = Max(is,id);
    dam[i] = Min(is,id);
    D[i] = .5 - .25*(Fvec[is]+Fvec[id]);
    if(is == 0 || id == 0)
      Fvec[i] = 0.;
    else if(sire[i-1]==sire[i] && dam[i-1]==dam[i])
      Fvec[i] = Fvec[i-1];
    else {
      fi = -1.; Lvec[i] = 1.; j = i;
      while(j != 0) {
	k = j; c = .5*Lvec[k];
	ks = sire[k]; kd = dam[k];
	if(ks > 0) {
	  while(ANC[k] > ks) k = ANC[k];
	  Lvec[ks] += c;
	  if(ks != ANC[k]) {
	    ANC[ks] = ANC[k];
	    ANC[k] = ks; }
	  if(kd > 0) {
	    while(ANC[k] > kd) k = ANC[k];
	    Lvec[kd] += c;
	    if(kd != ANC[k]) {
	      ANC[kd] = ANC[k];
	      ANC[k] = kd; }
	  }
	}
	fi += Lvec[j]*Lvec[j]*D[j];
	Lvec[j] = 0.; k = j; j = ANC[j]; ANC[k] = 0;
      }
      Fvec[i] = fi;
    }
  }

  free_ivector(ANC, 1, n);
  free_dvector(Fvec, 0, n);
  free_dvector(Lvec, 0, n);
  ANC = NULL;
  Fvec = Lvec = NULL;
}

// Convert to UT and evacuate into c_??
void UT(int n) {
  int i, l, *index0;
  index0 = ivector(1, n);
  /* Work the diagonal elements */
  for(i=1; i<=n; i++) {
    c_val[i] = kinMGD[i];
    c_link[i] = 0;
    c_col[i] = index0[i] = i; /* point to where lists end */
  }

  CLast = n+1;
  /* Work the off diags */
  for(i=2; i<=n; i++) {
    l = i;
    while(l > 0) {
      if(icol[l] == 0) break; /* empty row list */
      /* Hock kinMGOff[l] to row the icol[l] list of c_??
	 and assign to it col i */
      Hock(icol[l], i, kinMGOff[l], index0);
      l = ilinks[l];
    }
  }
  free_ivector(index0, 1, n);
  index0 = NULL;
}

// Hock (val) to the (row)th list of c_?? and
// point to the last element in this list in index0
void Hock(int row, int col, double val, int *index0) {
  c_link[ index0[row] ] = CLast; c_link[CLast] = 0;
  c_val[CLast] = val;  c_col[CLast] = col;
  index0[row] = CLast; /* the end has moved to CLast */
  CLast++;
}

// Never called
void Print(int n) {
  int i, j;
  FILE *fp;

  if( (fp= fopen("NRM_ilist", "w"))==NULL)
    Error("Print in SNRM2.h: can't open NRM_ilist file");

  for(i=1; i<IGLast; i++) {
    fprintf(fp,"%10d  %10d  %10.3f", icol[i], ilinks[i], kinMGOff[i]);
    if(i<=n) fprintf(fp,"  %10.3f\n", kinMGD[i]);
    else fprintf(fp,"\n");
  }
  printf("Printed out the list of A inverse to NRM_ilist\n");
  printf("%d animals were involved in a list of length %d\n",n,IGLast);
  fclose(fp);
  return;
}
/*
 * April 25th 2002
 * To build relationship inverse using list of sires of sires and
 * maternal grandsires of sires - Mrode's Linear Models book p 33.
 *
 * Objective: To be included in LNL program. Mainly used with threshold
 * model computations.
 *
 * Pedigree should be provided and composed of two columns one for sires
 * of sires and the second for maternal grand sires of sires. Recoding
 * is assumed, ie, the ith row contains sires and MGS of individual i
 *
 * First the LT is built and then converted to UT and put in the c_*
 * linked lists in accordance with similar practices and already-written
 * functions for NRM and GRM. Changing this although might
 * produce cleaner code, will not help in speed and will require
 * long time for changing mgs(), nrm(), and grm() functions.
 *
 * Inbreeding is NOT accounted for in computing the inverse!
 */
void MGSInv(int ind, int sr, int dm, double a, double b, double c) {
  int i, m, n;

  if(sr && dm) {
    m = Min(sr, dm); n = Max(sr, dm);
    kinMGD[ind] += (a);
    kinMGD[sr] += (a/4.0);
    kinMGD[dm] += (a/16.0);
    i = GetDex(n, m); kinMGOff[i] += (a/8.0);
    i = GetDex(ind, sr); kinMGOff[i] -= (a/2.0);
    i = GetDex(ind, dm); kinMGOff[i] -= (a/4.0);
  }

  else if(sr && !dm) {
    kinMGD[ind] += (b);
    kinMGD[sr] += (b/4.0);
    i = GetDex(ind, sr); kinMGOff[i] -= (b/2.0);
  }

  else if(dm && !sr) {
    kinMGD[ind] += (c);
    kinMGD[dm] += (c/16.0);
    i = GetDex(ind, dm); kinMGOff[i] -= (c/4.0);
  }
  else kinMGD[ind] += 1.0;
}

void mgs(int *sire, int *dam, int n) {
  int maxnze, i, j, k;
  double a, b, c;

  maxnze = 3*n; /* max number of the offdiagonal nze's */
  kinMGD = dvector(1, n);
  kinMGOff = dvector(1, maxnze);

  for(i=1; i<=n; i++) kinMGD[i] = 0;
  for(i=1; i<=IGLast; i++) kinMGOff[i] = 0;

  /* elements of D-1 if (a) both ss and mgs are known,
     (b) only ss is known, and (c) only mgs is known */
  a = 16.0/11.0; b = 4.0/3.0; c = 16.0/15.0;

  for(i=1; i<=n; i++)
    MGSInv(i, sire[i], dam[i], a, b, c);

  UT(n);

  free_dvector(kinMGD, 1, n);
  free_dvector(kinMGOff, 1, maxnze);
  kinMGD = kinMGOff = NULL;
}
















