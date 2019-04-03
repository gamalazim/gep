
/* =============================================================
 * 200803
 *  main C prog for relationship inverse setup - a C interface
 * Interface:
 * C: Other C functions
 * ============================================================= */

// cc -c this
// cc -c kinshipSetup.c - where kinshipSetup.h is also #included

#include "gt.h"
#include "analyze.h"
#include "kinshipSetUp.h"
#include "kinshipCi.h"
#include "debug.h"

  // length(buildOptn) = NumRandom + 2, last 2 are for animal and QTL.
  // Options of the first NumRandom are either 0 or 1, 0 or 2 for the first of the
  // last two (animal) and 0 or 3 for the second.
  // build options in kinship are buildOptn =

// For random() factors there are 2 options:
  // (0) No cov matrix inverse, nor pedigree
  // (1) -- Receive Inverse form R --

// For the ranim() factor there are three optins
  // (0) no animal factor to be built
  // (1) pedigree array is received with mgs's in the dam column
  // (2) pedigree array is received from R with dams in the dam column


// For rqtl() factors there are three optins
  // (0) no rqtl factors to be built
  // (3) Receive genotype info and fit a factor for each locus in a full model
  // (1) augment all qtl factors in a 'Minimized Animal Model Approach' - not yet implemented


// 'factor' is the (order - 1) of the factor whose cov. matrix starting pos is being
// searched for in the iaj vector. iajStartingPos is the 0 pos for the 1st factor (factor=0).
int getIajStartingPos(int factor, double *iaj, int *o) {
    int i, j, iajStartingPos = 0;
    for(i=0; i<factor; i++) {
      if(o[i] == 1) {
        j = iajStartingPos + Levels[NumFixed+2+i];
        iajStartingPos += ( Levels[NumFixed+2+i] + 2*iaj[j-1] );
      }
      else iajStartingPos += 3;
    }
    return iajStartingPos;
}


// 'factor' is the order of the factor whose cov. matrix is being extracted from the iaj vector.
void covMatrix(int factor, int iajStartingPos, double *iaj, int structNo) {
  int i, j, nze;
  int ord;

  if(RSM) ord = Levels[NumFixed+1+factor];
  else ord = Levels[NumFixed+factor];

  StructSp[structNo].n = ord;
  StructSp[structNo].ia = ivector(0, ord);
  nze = StructSp[structNo].nze = iaj[iajStartingPos + ord -1];
  StructSp[structNo].a = dvector(0, nze);
  StructSp[structNo].ja = ivector(0, nze);

  // Fill ia
  // =======
  StructSp[structNo].ia[0] = 1;
  for(i=1; i<=ord; i++) {
    j = iajStartingPos + i - 1;
    StructSp[structNo].ia[i] = iaj[j] + 1;
  }
  // Fill a and ja
  // =============
  for(i=0; i< nze; i++) {
    j = ord + iajStartingPos + i;
    StructSp[structNo].ja[i] = iaj[j];
    StructSp[structNo].a[i] = iaj[(j+nze)];
  }
}


void kinship(int *Rp, int *buildOptn, int *sire, int *dam, double *iaj, int loci, int *gtypes) {
  // sire/dam: brought from R - no need to allocate or free

  int i, k, startAt, COUNT;

  if(loci > MaxNumLoci) Error("Number of QT Loci exceeded limit");

  if(buildOptn[Rp[1]] == 1) { // dam = mgs
    RSM = 1;
    nind = Levels[NumFixed+1];
    createLists(sire, dam, nind);
    mgs(sire, dam, nind);
    sform(2, nind);
  }

  if(buildOptn[Rp[1]] == 2) { // dam = dam
    RSM = 2;
    nind = Rp[5+Rp[0]+Rp[1]];
    Rprintf("in buildOptn 2, just before CreateLists!\n");
    createLists(sire, dam, nind);
    Rprintf("CreateLists -- Done!\n");
    nrm(sire, dam, nind);
    Rprintf("nrm -- Done!\n");

    sform(2, nind);
    Rprintf("sform -- Done!\n");
    // SPrint(2, "nrm");
    // Struct2Default(2);
  }

  COUNT = 0;
  for(i=0; i<Rp[1]; i++) {
    if(buildOptn[i] == 1) { // if a cov mat inv is supplied
      COUNT++;
      if(COUNT > MaxNumCovFact) Error("Number of supplied cov matrices exceeded limit");

      if(RSM)
        k = 3 + i;
      else
        k = 2 + i;

      startAt = getIajStartingPos(i, iaj, buildOptn);
      covMatrix( (i+1), startAt, iaj, k );
      //SPrint(k, "sup.sp");
    }
  }

  if(buildOptn[ (Rp[1]+1) ] == 3) {
    GRM = 1;

    if(buildOptn[Rp[1]] != 2) { // just in case MQTLs are fit w/t a polygenic component, otherwise nind was set before.
      nind = Rp[5+Rp[0]+Rp[1]];
      createLists(sire, dam, nind);
    }

    for(i= 1; i<=loci; i++) {
      grm(i, sire, dam, gtypes, nind);
      k = 2 + NumRandom + i;
      sform(k, (nind*2) );
    }
    // Struct2Default(2 + NumRandom +  3);
  }
  freeLists();
}














