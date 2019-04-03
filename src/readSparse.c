
#include <stdio.h>
#include <R.h>

void Error(char *);

//void readSparse(char **fileName, int *ia, int *ja, double *a, int *order, int *nze) {
void readSparse(char **fileName, int *ia, int *ja, double *a, int *order, int *nze) {

  int i, j;
  char line[125];

  FILE *fp;
  if( (fp= fopen(fileName[0], "r"))==NULL) {
    printf("readSparse(): %s\n", fileName[0]);
    Error("can't open file"); }

  fgets(line, sizeof(line), fp);
//  sscanf(line, "%*[^0123456789]%d%*[^0123456789]%d", &local_order, &local_nze);

  for(i=1; i<=(*order); i++)
    fscanf(fp,"%d %d %d %lf", &j, &ia[i-1], &ja[i-1], &a[i-1]);

  for(i=((*order)+1); i<=(*nze); i++)
    fscanf(fp,"%d %d %lf\n", &j, &ja[i-1], &a[i-1]);

  fclose(fp);

}
