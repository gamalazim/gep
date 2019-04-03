
// 200802
// A C function to put upper tri received from R in ia/ja/a format
// warnings: the function recalculates nonzero elements and the count 
// must agree with the count in R. If needed a check must be done to stop w/error
// if nze > *NZE.

#include <R.h>
#include <stdio.h>

void makeSparse(double *data, int *N, int *NZE, int *ia, int *ja, double *a) {
    int i, j, nze=0;
    int ute = 0; // ((*N)*((*N) - 1))/2 + (*N) = all upper tri elements
	
  //printf("N = %d, NZE = %d\n\n", *N, *NZE);
    
    for(i=0; i<(*N); i++) {
      ia[i] = nze; 
 
      for(j=i; j<(*N); j++) {

        if(data[ute] != 0) {
          a[nze] = data[ute];
          ja[nze] = j+1;
          nze++;
          ia[i]++;
        }

        ute++;
      }
    }
}
