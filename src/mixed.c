
/* ===========
 *  main C prog
 * =========== */
#include "gt.h"
#include "analyze.h"
#include "TMAnalyze.h"
#include "kinshipCi.h"

void mixed(char **ff, char **funCall,
           double *RDataVec,
	         int *RNC,
	         int *RN,
           int *RFI,
	         int *Rparam,
	         int *Roptn,
	         int *Rsr,
	         int *Rdm,
	         double *Riaj,
	         int *fitLoci,
	         double *qtlVar,
	         int *gtypes,
           double *randomVar,
           char **solFile,
           int *maxit,
           int *orderOpt,
           int *memOptions
	   )
{

  getMemOptions(memOptions);

  /* Mixed Model Parameters - things that are done once  */
  MMParms(RDataVec, (*RNC), (*RN), (*RFI), Rparam, Roptn, (*fitLoci), qtlVar, randomVar);
  Rprintf("Completed MMParms and getting to kinship\n");
  kinship(Rparam, Roptn, Rsr, Rdm, Riaj, (*fitLoci), gtypes);
  // Rprintf("Completed kinship work\n");
  ORDER = 0; if((*orderOpt) == 1) ORDER = 1;
  tmrnd = 0;
  if(TM) {
    TM_Prep();

    solit = 0;
    do {
      solit++;
      Rprintf("Iteration No. %d\n", solit);
      TM_SetUp(rrr);
      TM_Solve(rrr);

      if(solcc <= SOLEPS) break;
      // if(solcc > oldsolcc) break;
      /* if(solit == 3) break; */ /* get out after 3 iterates */
      // cFree(0); /* free previous allocations of StructSp[ 0 ] */
    } while(solit < (*maxit));

    if(solcc > SOLEPS)
      printf("\n\nSolutions did not converge to the desired criterion (%lf) - use maxit = k, k > 100 \n\n\n", SOLEPS);
    if(solcc > oldsolcc)
      printf("\n\nSolutions Diverged - Exited Before Divergence ... \n\n");
    if(solcc <= SOLEPS)
      printf("\n\nSolutions Converged ... \n\n");

    RSolutions(solFile[0], ff, funCall);

    TM_Finish();
    printf("Done TM_Finish\n"); timedate();
  }

  else {
    Rprintf("Before SetUp\n");
    SetUp(rrr);
    Rprintf("Finished SetUp\n");
    DenseVectorSolve(rrr); printf("Done DenseVectorSolve ... \n");
    RSolutions(solFile[0], ff, funCall);
  }

  Finish(); /* free arrays and close files */
  printf("ANALYSIS COMPLETED AT:  "); timedate(); printf("\n");
}

















