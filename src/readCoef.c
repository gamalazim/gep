#include <stdio.h>
#include <R.h>
#define maxNumVar 1024 // Limit num. of variables to read to 1024

void Error(char *);

void readCoef(char **fileName, int *pv, int *nCats, int *FI, double *Est, double *Err, double *Rel) {
  int i, j, x, ct;
  int numVar= 2 + pv[0] + pv[1] + pv[2] + pv[4]; // num of variables to read
  //int *allLevels = (int *) malloc(numVar * sizeof(int)); // calling ivector fails!!
  int allLevels[maxNumVar]; // avoid anytype of mem. allocation from this fun, it just crashes R
  char line[1024];

  FILE *fp;
  if( (fp= fopen(fileName[0], "r"))==NULL) {
    printf("readCoef(): %s\n", fileName[0]);
    Error("can't open file"); }

  // Construct a vector with levels at each model variable

  //if(allLevels == NULL) Error("Allocation Failure in allLevels");
  for(i=1; i<=numVar; i++) allLevels[i] = 0;

  // [2] Fill up the vector
  if((*nCats) > 2) allLevels[1] = ((*nCats) - 2);      // Thresholds
  if((*FI)) allLevels[2] = 1;                          // Intercept
  for(i=1; i<=pv[0]; i++) allLevels[2+i] = pv[4+i];    // Fixed Factors
                                                    // Animal Factor
  if(pv[2] == 1) { allLevels[2+pv[0]+1] = pv[4+pv[0]+pv[1]+1]; x = 2+pv[0]+1; }
  else x = 2+pv[0];
                                                    // Other Random Factors
  for(i=1; i<=pv[1]; i++) allLevels[x+i] = pv[4+pv[0]+i];

  allLevels[2+pv[0]+pv[1]+pv[2]+1] = pv[4];            // Covariates

  // Skip top 4 lines
  for(i=1; i<=4; i++) fgets(line, sizeof(line), fp);

  ct = 0; // set counter (ct) to 0 before reading solutions

  // Read Thresholds
  if(allLevels[1] > 0) {
    for(i=1; i<=2; i++) fgets(line, sizeof(line), fp); // Empty line and title line
    for(i=1; i<=allLevels[1]; i++) {
      fgets(line, sizeof(line), fp);
      sscanf(line,"%d %lf %lf", &j, &Est[ct], &Err[ct]);
      ct++;
    }
  }

  // Read Intercept
  if(allLevels[2] > 0) {
    for(i=1; i<=2; i++) fgets(line, sizeof(line), fp); // Empty line and title line
    fgets(line, sizeof(line), fp);
    sscanf(line,"%d %lf %lf", &j, &Est[ct], &Err[ct]);
    ct++;
  }

  for(i=3; i<=numVar; i++) {
    for(j=1; j<=3; j++) fgets(line, sizeof(line), fp); // 1 Empty line + 1 Title line + 1 Factor type line
    if(i <= (2+pv[0]) || i > (2+pv[0]+pv[1]+pv[2]) ) {
      for(j=1; j<=allLevels[i]; j++) {
	fgets(line, sizeof(line), fp);
	sscanf(line,"%d %lf %lf", &x, &Est[ct], &Err[ct]);
	ct++;
      }
    }
    else {
      for(j=1; j<=allLevels[i]; j++) {
	fgets(line, sizeof(line), fp);
        sscanf(line,"%d %lf %lf %lf", &x, &Est[ct], &Err[ct], &Rel[ct]);
	ct++;
      }
    }
  }

  //free(allLevels);
  fclose(fp);
}

