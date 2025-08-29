/* linear C3D4 stress (CalculiX-style ordering) + SIMP
   Primitives kept:
     - constant-strain tet gradient (dN/dx) and volume
     - build isotropic D(E,nu) (Voigt: xx,yy,zz,xy,xz,yz; shear uses tensor strain)
     - B (6x12), eps = B u_e, sigma = D eps
     - von Mises and per-element stress output

   Conventions:
     - Connectivity is 0-based: conn[4] are node indices into co[3*Nnode], u[3*Nnode].
     - Voigt order = (xx, yy, zz, xy, xz, yz).
     - Strain is TENSOR form: ε_xy = 0.5*(du_x/dy + du_y/dx).
*/

#include <math.h>
#include <stddef.h>

/* Get C3D4 connectivity from CalculiX-style arrays.
   - i_elem_1b: element index i (1-based, as in Fortran)
   - ipkon: length >= ne, each entry 1-based pointer into kon
   - kon: global connectivity, node IDs are 1-based
   - out conn[4]: 0-based node indices for use in C
   Returns 0 on success, <0 on error.
*/
int ccx_get_conn_C3D4(const int *ipkon, const int *kon,
                      int i_elem_1b, int conn[4])
{
    if (i_elem_1b <= 0) return -1;
    /* Convert ipkon(i) from 1-based to 0-based index into kon[] */
    int idx0 = ipkon[i_elem_1b - 1];            /* still 1-based */
    if (idx0 <= 0) return -2;                   /* negative => inactive in CCX */
    idx0 -= 1;                                  /* now 0-based */

    /* kon holds 1-based node IDs. Make them 0-based for C. */
    conn[0] = kon[idx0 + 0] - 1;
    conn[1] = kon[idx0 + 1] - 1;
    conn[2] = kon[idx0 + 2] - 1;
    conn[3] = kon[idx0 + 3] - 1;

    /* Optional sanity: ensure all >=0 */
    for (int a=0; a<4; ++a) if (conn[a] < 0) return -3;
    return 0;
}


typedef struct 
{
    double dNdx[4][3];  /* gradients of N_a: [a=0..3][x,y,z] */
    double volume;      /* element volume */
} Tet4Geom;

/* determinant of a 3x3 */
static inline double det3(const double A[3][3])
{
    return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
         - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
         + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
}


/* compute constant dN/dx and volume for a linear tet */
static int tet4_geom(const double (*x)[3], Tet4Geom *tg)
{
    /* node positions: x[a][k], a=0..3, k=0..2 */
    double J[3][3] = {
        { x[1][0]-x[0][0], x[2][0]-x[0][0], x[3][0]-x[0][0] },
        { x[1][1]-x[0][1], x[2][1]-x[0][1], x[3][1]-x[0][1] },
        { x[1][2]-x[0][2], x[2][2]-x[0][2], x[3][2]-x[0][2] }
    };
    double detJ = det3(J);
    double sixV = detJ;         /* 6*Volume */
    double V    = sixV / 6.0;
    if (V <= 0.0) return -1;    /* inverted or degenerate */

    /* For linear tet: ∇N = J^{-T} ∇_ξ N, with reference grads:
       ∇_ξ N0 = [ -1, -1, -1 ]
       ∇_ξ N1 = [  1,  0,  0 ]
       ∇_ξ N2 = [  0,  1,  0 ]
       ∇_ξ N3 = [  0,  0,  1 ]  */
    /* Compute J^{-1} via adj/det */
    double invJ[3][3];
    invJ[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1]) / detJ;
    invJ[0][1] = -(J[0][1]*J[2][2]-J[0][2]*J[2][1]) / detJ;
    invJ[0][2] =  (J[0][1]*J[1][2]-J[0][2]*J[1][1]) / detJ;

    invJ[1][0] = -(J[1][0]*J[2][2]-J[1][2]*J[2][0]) / detJ;
    invJ[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0]) / detJ;
    invJ[1][2] = -(J[0][0]*J[1][2]-J[0][2]*J[1][0]) / detJ;

    invJ[2][0] =  (J[1][0]*J[2][1]-J[1][1]*J[2][0]) / detJ;
    invJ[2][1] = -(J[0][0]*J[2][1]-J[0][1]*J[2][0]) / detJ;
    invJ[2][2] =  (J[0][0]*J[1][1]-J[0][1]*J[1][0]) / detJ;

    /* ∇N in physical coords: ∇N = invJ^T * ∇_ξ N */
    const double dNr[4][3] = {
        {-1.0, -1.0, -1.0},
        { 1.0,  0.0,  0.0},
        { 0.0,  1.0,  0.0},
        { 0.0,  0.0,  1.0}
    };
    for(int a=0;a<4;++a){
        for(int k=0;k<3;++k){
            /* dNdx[a][k] = sum_r invJ[r][k] * dNr[a][r] = (invJ^T * dNr[a])_k */
            tg->dNdx[a][k] = invJ[0][k]*dNr[a][0] + invJ[1][k]*dNr[a][1] + invJ[2][k]*dNr[a][2];
        }
    }
    tg->volume = V;
    return 0;
}


/* isotropic 3D elasticity matrix, Voigt order (xx,yy,zz,xy,xz,yz), tensor shear */
static void iso_D(double E, double nu, double D[6][6])
{
    const double lam = (nu*E)/((1.0+nu)*(1.0-2.0*nu));
    const double G   = E/(2.0*(1.0+nu));
    /* zero */
    for(int i=0;i<6;++i) for(int j=0;j<6;++j) D[i][j]=0.0;
    /* normal blocks */
    D[0][0]=lam+2*G; D[0][1]=lam;     D[0][2]=lam;
    D[1][0]=lam;     D[1][1]=lam+2*G; D[1][2]=lam;
    D[2][0]=lam;     D[2][1]=lam;     D[2][2]=lam+2*G;
    /* shear (tensor) */
    D[3][3]=G; D[4][4]=G; D[5][5]=G;
}

/* gather element nodal displacement vector u_e (12) */
static void gather_ue(const double *u, const int *conn, double ue[12])
{
    for(int a=0;a<4;++a){
        int n = conn[a];
        ue[3*a+0] = u[3*n+0];
        ue[3*a+1] = u[3*n+1];
        ue[3*a+2] = u[3*n+2];
    }
}


/* B (6x12) for a linear tet from dNdx */

static void build_B(const Tet4Geom *tg, double B[6][12])
{
    /* zero */
    for(int i=0;i<6;++i) for(int j=0;j<12;++j) B[i][j]=0.0;
    for(int a=0;a<4;++a){
        const double dNx = tg->dNdx[a][0];
        const double dNy = tg->dNdx[a][1];
        const double dNz = tg->dNdx[a][2];
        const int c = 3*a;
        /* normal */
        B[0][c+0] = dNx;                 /* exx row */
        B[1][c+1] = dNy;                 /* eyy row */
        B[2][c+2] = dNz;                 /* ezz row */
        /* tensor shear: ε_xy, ε_xz, ε_yz */
        B[3][c+0] = 0.5*dNy; B[3][c+1] = 0.5*dNx;             /* ε_xy */
        B[4][c+0] = 0.5*dNz; B[4][c+2] = 0.5*dNx;             /* ε_xz */
        B[5][c+1] = 0.5*dNz; B[5][c+2] = 0.5*dNy;             /* ε_yz */
    }
}


static void matvec_6x6(const double A[6][6], const double x[6], double y[6])
{
    for(int i=0;i<6;++i)
    {
        double s=0.0;
        for(int j=0;j<6;++j) s += A[i][j]*x[j];
        y[i]=s;
    }
}


static void matvec_6x12(const double A[6][12], const double x[12], double y[6])
{
    for(int i=0;i<6;++i)
    {
        double s=0.0;
        for(int j=0;j<12;++j) s += A[i][j]*x[j];
        y[i]=s;
    }
}


static double von_mises_from_sigma(const double s[6])
{
    const double sxx=s[0], syy=s[1], szz=s[2], sxy=s[3], sxz=s[4], syz=s[5];
    const double t1 = (sxx - syy)*(sxx - syy)
                    + (syy - szz)*(syy - szz)
                    + (szz - sxx)*(szz - sxx);
    const double t2 = 6.0*(sxy*sxy + sxz*sxz + syz*syz);
    return sqrt(0.5*t1 + 0.5*t2);
}

/* SIMP modulus */
static inline double simp_E(double rho, double p, double E0, double Emin)
{
    if (rho<0.0) rho=0.0; if (rho>1.0) rho=1.0;
    return Emin + pow(rho, p)*(E0 - Emin);
}

/* ==== PUBLIC: compute stress & von Mises for one C3D4 element ==== */
int element_stress_C3D4_linear(
    const double *co,    /* [3*Nnode] node coords (x,y,z interleaved: co[3*n + k]) */
    const double *u,     /* [3*Nnode] nodal displacements */
    const int *conn,     /* [4] 0-based node indices for the element */
    double rhoPhys, double p, double E0, double Emin, double nu,
    double sigma[6],     /* out: (xx,yy,zz,xy,xz,yz) at the unique GP */
    double *sigma_vm,    /* out: von Mises */
    double *volume       /* out: element volume */)
    {
        /* gather element coords */
        double x[4][3];
        for(int a=0;a<4;++a)
        {
            int n=conn[a];
            x[a][0]=co[3*n+0]; x[a][1]=co[3*n+1]; x[a][2]=co[3*n+2];
        }

        Tet4Geom tg;
        if (tet4_geom(x, &tg)!=0) return -1; /* inverted/degenerate */

        /* SIMP elasticity, fixed nu */
        const double E = simp_E(rhoPhys, p, E0, Emin);

        /* build operators */
        double D[6][6]; iso_D(E, nu, D);
        double B[6][12]; build_B(&tg, B);

        /* strain = B * u_e */
        double ue[12]; gather_ue(u, conn, ue);
        double eps[6]; matvec_6x12(B, ue, eps);

        /* stress = D * strain */
        matvec_6x6(D, eps, sigma);

        /* Compute von-Misses */
        *sigma_vm = von_mises_from_sigma(sigma);

        if (volume) *volume = tg.volume;
        return 0;
    }