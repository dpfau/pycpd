/*
* Andriy Myronenko
* Adapted for Python by David Pfau
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

void comp_p(double*, double*, double*, double*, double*, double*, double*, double*, int, int, int);
void comp_truncate(double*, double*, double*, double*, double*, double*, double*, double*, int, int, int, double*);
void comp_correspondence(double*, double*, double*, double*, double*, int, int, int);
