
#define rr .05 /* recombination fraction */
#define maxSz  150000000  /* size of required G */
#define ImaxSz 300000000   /* size of the UT of the inverse linked lists */

double **values, *f, **iv3, **iv4;
int *col, *links, *icol, *ilinks;
int GLast, IGLast;
double *kinMGD, *kinMGOff; /* kinMGD for diag elem'ts and kinMGOff for the off diags */

void createLists(int *, int *, int);
void freeLists();
void grm(int, int *, int *, int *, int);
void nrm(int *, int *, int);
void MGSInv(int, int, int, double, double, double);
void mgs(int *, int *, int);
void fillB(int, int *, int *, int *, int **, int);

void List(int, int);
void IList(int, int);
void SortIList(int);
void PrintList(int);
void PrintIList(int);
void Add2Inverse(int, int, int, double **);
void UseWang(double **, double **, int *);
void NMNp(double **, double **, int, int, double **);
int GetDex(int, int);
void Add2G(int, int, int, double **);
double *IKinship(int, int);
void MakeQ(double **, int, int, int, int **);
void PrintQ(double **);
void PrintP(double **);
void MakeUT(int);
void A_inverse(int, int, int, double);
void GetD(double *, int *, int *, int);
void UT(int);
void Hock(int, int, double, int *);
void Print(int);

extern int CLast;
extern int *c_col;
extern int *c_link;
extern double *c_val;
