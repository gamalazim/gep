
#include <stdio.h>
#include <R.h>
#include "analyze.h"

void Error(char *);

void writeSparse(char **fileName, int *ia, int *ja, double *a, int *order, int *nze){

  int i;
  FILE *fp;
  if( (fp= fopen(fileName[0], "w"))==NULL) {
    printf("writeSparse(): %s\n", fileName[0]);
    Error("can't open file"); }

  fprintf(fp,"SPARSE FORM, ORDER %d, NZE %d\n", (*order), (*nze));
  for(i=1; i<=(*nze); i++) {
      if(i <= (*order))
        fprintf(fp,"%3.d %3.d %3.d     %g\n",i, ia[i-1], ja[i-1], a[i-1]);
      else
        fprintf(fp,"%3.d     %3.d     %g\n",i, ja[i-1], a[i-1]);
    }
    fclose(fp);

}

void writeSparseCoef(int s) {

  int i;
  FILE *fp;
  if( (fp= fopen("LHS", "w"))==NULL) {
    printf("writeSparseCoef(): %s\n", "LHS");
    Error("can't open file"); }

  fprintf(fp,"SPARSE FORM, ORDER %d, NZE %d\n", (StructSp[s].n), (StructSp[s].nze));
  for(i=1; i<=(StructSp[s].nze); i++) {
    if(i <= (StructSp[s].n))
      fprintf(fp,"%d %d %d     %g\n",i, StructSp[s].ia[i-1], StructSp[s].ja[i-1], StructSp[s].a[i-1]);
    else
      fprintf(fp,"%d %d     %g\n",i, StructSp[s].ja[i-1], StructSp[s].a[i-1]);
  }
  fclose(fp);

}
