

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

void adjoint_eq_2_node(
    ITG *nk,
    ITG *nactdof,     /* (mi[2]+1) * (*nk) */
    ITG *nboun,
    ITG *nodeboun,    /* length *nboun, 1-based nodes */
    ITG *ndirboun,    /* length *nboun, values 1..mi[2] */
    ITG *mi,          /* mi[2] = #mech DOFs per node */
    double *lambda,   /* ladj, size (mi[2]+1)*(*nk) in Fortran layout */
    double *B_adj     /* adjoint solution vector (length = neq) */
)
{
    FORTRAN(resultsini_adjoint_linstatic_nompc,
            (nk, lambda, nactdof, B_adj,
             nodeboun, ndirboun, nboun, mi));
}
