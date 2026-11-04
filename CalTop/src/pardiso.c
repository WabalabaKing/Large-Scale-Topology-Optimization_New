/*     CalculiX - A 3-dimensional finite element program                   */
/*              Copyright (C) 1998-2018 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#ifdef PARDISO

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <unistd.h>
#include "mkl_spblas.h"
#include "CalculiX.h"
#include "pardiso.h"
#include "bddc.h"

#ifdef PROFILING_ON
#include <TAU.h>
#endif

ITG *icolpardiso=NULL,*pointers=NULL,iparm[64];
long long pt[64];
double *aupardiso=NULL;
/* double dparm[64];  not used */
ITG nthread_mkl=0;
int g_useNK=0;   /* --NK: Newton-Krylov in nonlingeo (see end of file) */

/* factor kept by pardiso_factor until pardiso_cleanup */
static ITG factor_alive=0,alive_neq=0,alive_symmetryflag=0;

/* Newton-Krylov state (separate PARDISO handle, see end of file) */
static ITG nk_alive=0;
static void nk_solve(double *b,ITG *neq,ITG *symmetryflag,ITG *nrhs,double tol);
static void nk_release(void);

/* p-multigrid state (--MG true, see end of file) */
int g_useMG=0;
static ITG mg_alive=0;
static void mg_solve(double *b,ITG *neq,ITG *nrhs,double tol,double atol,
		     ITG adjoint);
static void mg_release(void);
static ITG mg_store(double *ad, double *au, double *adb, double *aub,
                double *sigma,ITG *icol, ITG *irow,
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3);
static double mg_tol(int newton);
static double mg_atol=0.;
static ITG mg_storing=0;
static void mg_options(void);
static ITG mg_directadj;

/* CALTOP_SOLVER=bddc adjoint state (pardiso_store_factor) */
static ITG bd_adj=0,bd_fb_mtype=-2,*bd_fb_p=NULL,*bd_fb_c=NULL;
static double *bd_fb_a=NULL,*bd_adj_ad=NULL,*bd_adj_au=NULL;
static ITG *bd_adj_jq=NULL,*bd_adj_irow=NULL;
static ITG nk_adopt(ITG *neq,ITG mtype,ITG *p,ITG *c,double *a);
static double *mg_guess;              /* converged displacement (pardiso_mg_guess) */
static ITG mg_guessn;
static void nk_numeric(ITG *neq);

/* PARDISO two-level factorization (iparm[23]=1): off (tested, no gain) */
static ITG pardiso_twolevel(void)
{
  return 0;
}
static ITG mg_mkl_mv;
#define NK_TOL_NEWTON 1.e-3  /* FGMRES rel. residual, Newton solves           */
#define NK_TOL_SOLVE  1.e-8  /* FGMRES rel. residual, other solves            */
#define NK_KRYLOV_MAXIT 30   /* FGMRES iterations before refactoring          */
#define NK_GIVEUP_MIN   10   /* earliest iteration to give up on FGMRES when  */
                             /*   the residual trend needs > MAXIT iterations */
#define NK_KRYLOV_M     15   /* restart length (memory ~ (2m+1)*neq doubles)  */
#define NK_ORDERING      3   /* iparm[1]: 2 METIS, 3 parallel METIS           */
#define MG_TOL_SOLVE_ 1.e-8  /* p-multigrid CG rel. residual, other solves    */
/* char envMKL[32];   moved to pardiso.h */
//char envMKL[32];

/*
 * ----------------------------------------------------------------------
 * Function: pardiso_factor
 * ----------------------------------------------------------------------
 * Purpose:
 *    Prepares and factors the global linear system matrix using the
 *    Intel MKL PARDISO sparse direct solver. This function converts
 *    the CalculiX internal sparse matrix format into the CSR format
 *    required by PARDISO, applies an optional shift (K - sigma*M),
 *    configures threading, initializes PARDISO control parameters,
 *    and performs the symbolic and numerical factorization (phase = 12).
 *
 * Arguments:
 *    ad       - [in]  Diagonal entries of stiffness matrix K.
 *    au       - [in]  Off-diagonal entries of stiffness matrix K
 *                     (format depends on symmetryflag and inputformat).
 *    adb      - [in]  Diagonal entries of auxiliary matrix (e.g., mass) M.
 *    aub      - [in]  Off-diagonal entries of M.
 *    sigma    - [in]  Scalar shift; if non-zero, factorizes K - sigma*M.
 *
 *    icol     - [in]  Number of off-diagonal entries in each column of K.
 *    irow     - [in]  Row indices for off-diagonal terms in K.
 *
 *    neq      - [in]  Number of equations (matrix dimension).
 *    nzs      - [in]  Number of off-diagonal nonzeros.
 *
 *    symmetryflag - [in]  0 = symmetric, 1 = unsymmetric matrix.
 *    inputformat   - [in]  Storage scheme used by CalculiX
 *                          (1 = split triangular, 3 = column-based unsymmetric).
 *
 *    jq       - [in]  Indices of upper triangular entries (inputformat = 1).
 *    nzs3     - [in]  Auxiliary offset for upper triangular storage.
 *
 * Output:
 *    No explicit return. Allocates and fills global arrays:
 *       pointers     - CSR row pointer array (1-based)
 *       icolpardiso  - CSR column indices
 *       aupardiso    - CSR matrix values
 *       pt[], iparm[] configured for later solve call
 *
 * Notes:
 *    - Sets MKL_NUM_THREADS based on CCX_NPROC_EQUATION_SOLVER,
 *      MKL_NUM_THREADS, or OMP_NUM_THREADS (first call only).
 *    - Performs PARDISO symbolic + numeric factorization (phase = 12).
 *    - Must be followed by pardiso_solve() and pardiso_cleanup().
 *
 * ----------------------------------------------------------------------
 */

static void pardiso_threads(void)
	{
  		char *env;
  		ITG nthread,nthread_v;

		/* set MKL_NUM_THREADS to min(CCX_NPROC_EQUATION_SOLVER,OMP_NUM_THREADS)
   			must be done once  */
  		if (nthread_mkl == 0) 
		{
    		nthread=1;
    		env=getenv("MKL_NUM_THREADS");
    		if(env) 
			{
      			nthread=atoi(env);
			}
    		else 
			{
      			env=getenv("OMP_NUM_THREADS");
      			if(env) 
				{
					nthread=atoi(env);
				}
    		}
    		env=getenv("CCX_NPROC_EQUATION_SOLVER");
    		if(env) 
			{
      			nthread_v=atoi(env);
      			if (nthread_v <= nthread) 
				{
					nthread=nthread_v;
				}
    		}
    		if (nthread < 1) 
			{
				nthread=1;
			}
    		sprintf(envMKL,"MKL_NUM_THREADS=%" ITGFORMAT "",nthread);  
    		putenv(envMKL);
    		nthread_mkl=nthread;
  		}
    
//  		printf(" number of threads =% d\n\n",nthread_mkl);
	}

/* CalculiX storage -> PARDISO CSR (pointers, icolpardiso, aupardiso);
   returns the PARDISO matrix type */
static ITG pardiso_csr(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
	{
  		ITG i,j,k,l,mtype,*irowpardiso=NULL,kflag,kstart,n,k2,ndim;

  		if(*symmetryflag==0)
		{

      		/* symmetric matrix; the subdiagonal entries are stored
         	column by column in au, the diagonal entries in ad;
         	pardiso needs the entries row per row */      

      		mtype=-2;

      		/* CalTop: rows filled in parallel; row i is the diagonal followed
      		   by column i of the strict lower triangle (same layout as the
      		   original reverse loop) */
      		ndim=*neq+*nzs;

      		NNEW(pointers,ITG,*neq+1);
      		NNEW(icolpardiso,ITG,ndim);
      		NNEW(aupardiso,double,ndim);

      		{
      		    ITG *l0=NULL;
      		    NNEW(l0,ITG,*neq+1);
      		    for(i=0;i<*neq;i++) l0[i+1]=l0[i]+icol[i];
      		    for(i=0;i<*neq;i++) pointers[i]=l0[i]+i+1;
      		    pointers[*neq]=ndim+1;
#pragma omp parallel for private(j,k,l)
      		    for(i=0;i<*neq;i++){
      			k=l0[i]+i;
      			icolpardiso[k]=i+1;
      			aupardiso[k]=(*sigma==0.) ? ad[i] : ad[i]-*sigma*adb[i];
      			for(j=0;j<icol[i];j++){
      			    l=l0[i]+j;
      			    icolpardiso[k+1+j]=irow[l];
      			    aupardiso[k+1+j]=(*sigma==0.) ? au[l]
      				: au[l]-*sigma*aub[l];
      			}
      		    }
      		    SFREE(l0);
      		}
  		}
		else
		{
      		mtype=11;

      		if(*inputformat==3)
			{

          		/* off-diagonal terms  are stored column per
            	 column from top to bottom in au;
             	diagonal terms are stored in ad  */

	  			ndim=*neq+*nzs;
	  			NNEW(pointers,ITG,*neq+1);
	  			NNEW(irowpardiso,ITG,ndim);	  
	  			NNEW(icolpardiso,ITG,ndim);
	  			NNEW(aupardiso,double,ndim);
	  
	  			k=0;
	  			k2=0;
	  			for(i=0;i<*neq;i++)
				{
	      			for(j=0;j<icol[i];j++)
					{
		  				if(au[k]>1.e-12||au[k]<-1.e-12)
						{
		      				icolpardiso[k2]=i+1;
		      				irowpardiso[k2]=irow[k];
		      				aupardiso[k2]=au[k];
		      				k2++;		  
		  				}
		  				k++;	      
	      			}	  
	  			}  
	  			/* diagonal terms */  
	  			for(i=0;i<*neq;i++)
				{
	      			icolpardiso[k2]=i+1;
	      			irowpardiso[k2]=i+1;
	      			aupardiso[k2]=ad[i];
	      			k2++;	  
	  			}
	  			ndim=k2;
	  
	  			/* pardiso needs the entries row per row; so sorting is
	     		needed */ 
	  
	  			kflag=2;
	  			FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
			    	&ndim,&kflag));
	  
	  			/* sorting each row */
	  
	  			k=0;
	  			pointers[0]=1;
	  			for(i=0;i<*neq;i++)
				{
	      			j=i+1;
	      			kstart=k;
	      			do
					{
		  				if(irowpardiso[k]!=j )
						{
		      				n=k-kstart;		  
		      				FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
				       		&n,&kflag));
		     	 			pointers[i+1]=k+1;
		      				break;  
		  				}
						else
						{
		      				if(k+1==ndim)
							{
			  					n=k-kstart+1;	  
			  					FORTRAN(isortid,(&icolpardiso[kstart],
                            	      &aupardiso[kstart],&n,&kflag));
			  					break;	       
		      				}
							else
							{
			  					k++;	       
		      				}  
		  				}
	      			}
					while(1);
	  			}
	  			pointers[*neq]=ndim+1;
	  			SFREE(irowpardiso);

      		}
			else if(*inputformat==1)
			{
	  
          		/* lower triangular matrix is stored column by column in
            	 au, followed by the upper triangular matrix row by row;
             		the diagonal terms are stored in ad */

          		/* reordering lower triangular matrix */

	  			ndim=*nzs;
	  			NNEW(pointers,ITG,*neq+1);
	  			NNEW(irowpardiso,ITG,ndim);
	  			NNEW(icolpardiso,ITG,ndim);
	  			NNEW(aupardiso,double,ndim);
				
	  			k=0;
	  			for(i=0;i<*neq;i++)
				{
	      			for(j=0;j<icol[i];j++)
					{
		  				icolpardiso[k]=i+1;
		  				irowpardiso[k]=irow[k];
		  				aupardiso[k]=au[k];
		  				k++;
	      			}
	  			}
	  
	  			/* pardiso needs the entries row per row; so sorting is
	     			needed */
	  
	  			kflag=2;
	  			FORTRAN(isortiid,(irowpardiso,icolpardiso,aupardiso,
	  				&ndim,&kflag));
	  
	  			/* sorting each row */
	  
	  			k=0;
	  			pointers[0]=1;
	  			if(ndim>0)
				{
	      			for(i=0;i<*neq;i++)
					{
		  				j=i+1;
		  				kstart=k;
		  				do
						{
	  	      				/* end of row reached */

		      				if(irowpardiso[k]!=j)
							{
			  					n=k-kstart;
			  					FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
					   				&n,&kflag));
			  					pointers[i+1]=k+1;
			  					break;
		      				}
							else
							{
		          				/* end of last row reached */
			  					if(k+1==ndim)
								{
			      					n=k-kstart+1;
			      					FORTRAN(isortid,(&icolpardiso[kstart],&aupardiso[kstart],
					       				&n,&kflag));
			      					break;
			  					}
								else
								{
			      					/* end of row not yet reached */
			      					k++;
			  					}
		      				
							}
		  				}
						while(1);
	      			}
	  			}
	  			pointers[*neq]=ndim+1;
	  			SFREE(irowpardiso);

	  			/* composing the matrix: lower triangle + diagonal + upper triangle */

	  			ndim=*neq+2**nzs;
	  			RENEW(icolpardiso,ITG,ndim);
	  			RENEW(aupardiso,double,ndim);
	  			k=ndim;

	  			for(i=*neq-1;i>=0;i--)
				{
	      			l=k+1;
	      			for(j=jq[i+1]-1;j>=jq[i];j--)
					{
		  				icolpardiso[--k]=irow[j-1];
		  				aupardiso[k]=au[j+*nzs3-1];
	      			}

	      			icolpardiso[--k]=i+1;
	      			aupardiso[k]=ad[i];

	      			for(j=pointers[i+1]-1;j>=pointers[i];j--)
					{
		  				icolpardiso[--k]=icolpardiso[j-1];
		  				aupardiso[k]=aupardiso[j-1];
	      			}

	      			pointers[i+1]=l;
	  			}
	  			pointers[0]=1;
      		}
  		}

  		return mtype;
	}

void pardiso_factor(double *ad, double *au, double *adb, double *aub, 
                double *sigma,ITG *icol, ITG *irow, 
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
	{

	//	#ifdef PROFILING_ON
    //    	TAU_PROFILE_TIMER(t_pardiso_factor,"PARDISO: Factorization","",TAU_USER);
	//		TAU_PROFILE_START(t_pardiso_factor);
    //	#endif

  		ITG i,maxfct=1,mnum=1,phase=12,nrhs=1,*perm=NULL,mtype,
      		msglvl=0,error=0;
  		double *b=NULL,*x=NULL;

  		iparm[0]=0;
  		pardiso_threads();

  		for(i=0;i<64;i++)
		{
			pt[i]=0;
		}

  		mtype=pardiso_csr(ad,au,adb,aub,sigma,icol,irow,neq,nzs,
  			  symmetryflag,inputformat,jq,nzs3);

  		FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
			   pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  		factor_alive=1;
  		alive_neq=*neq;
  		alive_symmetryflag=*symmetryflag;

	//	#ifdef PROFILING_ON
    //    	TAU_PROFILE_STOP(t_pardiso_factor);
    //	#endif
				   
  		return;
	}
/*
 * ----------------------------------------------------------------------
 * Function: pardiso_solve
 * ----------------------------------------------------------------------
 * Purpose:
 *    Solves the linear system A * x = b using the factorization
 *    produced earlier by pardiso_factor(). The solution overwrites b.
 *    This function invokes PARDISO with phase = 33 (solve only).
 *
 * Arguments:
 *    b            - [in/out]  Right-hand side vector(s). On return,
 *                              contains the computed solution(s).
 *    neq          - [in]      Number of equations.
 *    symmetryflag - [in]      0 = symmetric system, 1 = unsymmetric.
 *    nrhs         - [in]      Number of right-hand sides.
 *
 * Output:
 *    - Solution is written into b.
 *    - No matrix memory is modified; uses the existing PARDISO pt[] data.
 *
 * Notes:
 *    - Requires that pardiso_factor() has already been called.
 *    - Allocates a temporary solution vector x and copies results into b.
 *    - Thread count is inherited from the first call to pardiso_factor().
 *
 * ----------------------------------------------------------------------
 */

void pardiso_solve(double *b, ITG *neq,ITG *symmetryflag,ITG *nrhs)
{

//	#ifdef PROFILING_ON
//        TAU_PROFILE_TIMER(t_pardiso_solve,"PARDISO: Solve","",TAU_USER);
//		TAU_PROFILE_START(t_pardiso_solve);
//    #endif

  ITG maxfct=1,mnum=1,phase=33,*perm=NULL,mtype,
    	msglvl=0,i,error=0;
  double *x=NULL;

  /* CALTOP_SOLVER=bddc adjoints (pardiso_store_factor) */
  if(bd_adj){
    char *ev=getenv("CALTOP_BDDC_TOL_ADJ");
    double tol=(ev!=NULL&&atof(ev)>0.) ? atof(ev) : 1.e-5;
    ITG r;
    for(r=0;r<*nrhs;r++)        /* a solve that fails stops the run */
      bddc_solve(&b[(size_t)r**neq],*neq,tol,1,bd_adj_ad,bd_adj_au,
		 bd_adj_jq,bd_adj_irow,
		 (mg_guess!=NULL&&mg_guessn==*neq) ? mg_guess : NULL);
    return;
  }
  /* matrix stored by pardiso_store_factor / pardiso_factor_solve */
  if(!factor_alive&&mg_alive)
  {
    mg_solve(b,neq,nrhs,mg_tol(0),0.,1);
    return;
  }
  if(!factor_alive&&nk_alive)
  {
    nk_solve(b,neq,symmetryflag,nrhs,NK_TOL_SOLVE);
    return;
  }

  if(*symmetryflag==0)
  {
    mtype=-2;
  }
  else
  {
    mtype=11;
  }

  iparm[0]=0;

	/* pardiso_factor has been called befor, MKL_NUM_THREADS=nthread_mkl is set*/
  NNEW(x,double,*nrhs**neq);

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,nrhs,iparm,&msglvl,
                   b,x,&error));

  for(i=0;i<*nrhs**neq;i++){b[i]=x[i];}
  SFREE(x);

//	#ifdef PROFILING_ON
//    	TAU_PROFILE_STOP(t_pardiso_solve);
//	#endif

  return;
}

/*
 * ----------------------------------------------------------------------
 * Function: pardiso_cleanup
 * ----------------------------------------------------------------------
 * Purpose:
 *    Releases all PARDISO internal memory associated with the matrix
 *    factorization and frees the CSR arrays allocated during
 *    pardiso_factor(). Uses PARDISO phase = -1 (memory release).
 *
 * Arguments:
 *    neq          - [in]  Number of equations.
 *    symmetryflag - [in]  0 = symmetric matrix, 1 = unsymmetric.
 *
 * Output:
 *    No return value. Frees global arrays:
 *        icolpardiso
 *        aupardiso
 *        pointers
 *    and clears all PARDISO internal structures in pt[].
 *
 * Notes:
 *    - Must be called after all solves are complete.
 *    - Safe to call even in one-shot solve (used by pardiso_main()).
 *
 * ----------------------------------------------------------------------
 */


void pardiso_cleanup(ITG *neq,ITG *symmetryflag)
{

//	#ifdef PROFILING_ON
//        TAU_PROFILE_TIMER(t_pardiso_cleanup,
//                          "PARDISO: Cleanup",
//                          "",
//                          TAU_USER);

//        TAU_PROFILE_START(t_pardiso_cleanup);
//    #endif

  ITG maxfct=1,mnum=1,phase=-1,*perm=NULL,nrhs=1,mtype,
      msglvl=0,error=0;
  double *b=NULL,*x=NULL;

  if(*symmetryflag==0)
  {
    mtype=-2;
  }
  else
  {
    mtype=11;
  }

  FORTRAN(pardiso,(pt,&maxfct,&mnum,&mtype,&phase,neq,aupardiso,
		   pointers,icolpardiso,perm,&nrhs,iparm,&msglvl,
                   b,x,&error));

  SFREE(icolpardiso);
  SFREE(aupardiso);
  SFREE(pointers);

  factor_alive=0;

//    #ifdef PROFILING_ON
//        TAU_PROFILE_STOP(t_pardiso_cleanup);
//    #endif

  return;
}

/* 1 if a matrix is held for pardiso_solve (pardiso_factor or NK) */
ITG pardiso_factor_alive(void)
{
  return factor_alive||nk_alive||mg_alive||bd_adj;
}

/* release the held factorization(s), if any */
void pardiso_release(void)
{
  if(factor_alive) pardiso_cleanup(&alive_neq,&alive_symmetryflag);
  mg_release();
  nk_release();
  bd_adj=0;
  if(bd_fb_p!=NULL){SFREE(bd_fb_p);SFREE(bd_fb_c);SFREE(bd_fb_a);}
}

/*
 * ----------------------------------------------------------------------
 * Function: pardiso_main
 * ----------------------------------------------------------------------
 * Purpose:
 *    Convenience wrapper that performs the complete PARDISO workflow
 *    for a single solve:
 *       1) Assemble and factor the system matrix (pardiso_factor),
 *       2) Solve A * x = b for the given right-hand side(s)
 *          (pardiso_solve),
 *       3) Release all PARDISO and matrix-related memory
 *          (pardiso_cleanup).
 *
 * Arguments:
 *    ad       - [in]  Diagonal entries of stiffness matrix K.
 *    au       - [in]  Off-diagonal entries of stiffness matrix K.
 *    adb      - [in]  Diagonal entries of auxiliary matrix (e.g., mass) M.
 *    aub      - [in]  Off-diagonal entries of M.
 *    sigma    - [in]  Scalar shift; if non-zero, K - sigma*M is factored.
 *
 *    b        - [in/out]  Right-hand side(s). On return, contains the
 *                          solution vector(s).
 *
 *    icol     - [in]  Number of off-diagonal entries in each column of K.
 *    irow     - [in]  Row indices for off-diagonal terms.
 *
 *    neq      - [in]  Number of equations (matrix dimension).
 *    nzs      - [in]  Number of off-diagonal nonzeros.
 *
 *    symmetryflag - [in]  0 = symmetric system, 1 = unsymmetric system.
 *    inputformat   - [in]  CalculiX matrix storage format
 *                          (1 = split triangular, 3 = column-based).
 *
 *    jq       - [in]  Column pointer array for upper triangle
 *                     (used when inputformat = 1).
 *    nzs3     - [in]  Offset into au for upper triangular entries.
 *
 *    nrhs     - [in]  Number of right-hand sides.
 *
 * Output:
 *    - Overwrites b with the solution(s) of the linear system.
 *    - Internally allocates and frees all data structures for PARDISO.
 *
 * Notes:
 *    - Returns immediately if *neq == 0 (empty system).
 *    - Suitable for single load case / one-shot solves.
 *      For repeated solves with the same matrix and different RHS,
 *      call pardiso_factor(), then pardiso_solve() repeatedly,
 *      and finally pardiso_cleanup().
 *
 * ----------------------------------------------------------------------
 */
void pardiso_main(double *ad, double *au, double *adb, double *aub, 
         double *sigma,double *b, ITG *icol, ITG *irow, 
	ITG *neq, ITG *nzs,ITG *symmetryflag,ITG *inputformat,
	ITG *jq, ITG *nzs3,ITG *nrhs)
	{

	//	#ifdef PROFILING_ON
    //    	TAU_PROFILE_TIMER(t_pardiso_total,"PARDISO: Total","",TAU_USER);
    //    	TAU_PROFILE_START(t_pardiso_total);
    //	#endif

  		if(*neq==0)
		{
		//	#ifdef PROFILING_ON
      //  		TAU_PROFILE_STOP(t_pardiso_total);
    	//	#endif
			return;
		}
  		pardiso_factor(ad,au,adb,aub,sigma,icol,irow, 
			neq,nzs,symmetryflag,inputformat,jq,nzs3);

  		pardiso_solve(b,neq,symmetryflag,nrhs);

  		pardiso_cleanup(neq,symmetryflag);

	//	#ifdef PROFILING_ON
    //    	TAU_PROFILE_STOP(t_pardiso_total);
    //	#endif

  		return;
	}

/*
 * ----------------------------------------------------------------------
 * Newton-Krylov (--NK true, nonlingeo only)
 * ----------------------------------------------------------------------
 * Separate PARDISO handle; the stock functions above are not affected.
 *  - symbolic analysis (phase 11) is redone only when the CSR pattern
 *    changes (FNV-1a hash); numerical factorization is phase 22
 *  - a new matrix is solved by right-preconditioned FGMRES(m) with the
 *    last factorization as preconditioner (phase 33); refactor only if
 *    FGMRES does not converge in NK_KRYLOV_MAXIT iterations, or its
 *    residual trend shows it would not (checked from NK_GIVEUP_MIN on)
 *  - pardiso_factor_solve: Newton step (FGMRES tol NK_TOL_NEWTON)
 *    pardiso_store_factor: store and factor the converged K_T in this
 *    handle (symbolic analysis reused); pardiso_solve then back-solves
 * ----------------------------------------------------------------------
 */

static long long nk_pt[64];
static ITG nk_iparm[64],*nk_pointers=NULL,*nk_icol=NULL;
static double *nk_au=NULL;
static ITG nk_valid=0,nk_stale=0,nk_neq=-1,nk_mtype=0,nk_nnz=-1;
static unsigned long long nk_hash=0ULL;
static long n_symbolic=0,n_factor=0,n_kry_ok=0,n_kry_fail=0,n_kry_iter=0,
            n_backsolve=0;
static double t_symbolic=0.,t_factor=0.,t_krylov=0.,t_backsolve=0.;

/* single-precision factor of the converged K_T for the adjoints (iparm[27]=1,
   FGMRES in double precision on top): off; the code path is kept */
static ITG nk_single=0;
static float *nk_auf=NULL;


/* CALTOP_PARDISO_PERM_DIR=<folder>: PARDISO's fill-reducing ordering
   (most of the analysis time) depends only on the sparsity pattern, which
   is the same in every optimization iteration. The first run saves it
   there, keyed by a hash of the pattern; later runs pass it back
   (iparm[4]=1) instead of running METIS again. Results are unchanged. */
static const char *perm_dir(void)
{
  static int init=0;
  static const char *dir=NULL;
  if(!init){
    char *env=getenv("CALTOP_PARDISO_PERM_DIR");
    init=1;
    if(env!=NULL&&env[0]!='\0') dir=env;
  }
  return dir;
}
static void perm_name(char *buf,size_t len,const char *tag,
		      unsigned long long h,ITG n,ITG nnz)
{
  snprintf(buf,len,"%s/caltop_perm_%s_%d_%d_%016llx.bin",perm_dir(),tag,
	   (int)n,(int)nnz,h);
}
/* returns a permutation read from the folder, or NULL */
static ITG *perm_load(const char *tag,unsigned long long h,ITG n,ITG nnz)
{
  char name[4096];
  FILE *f;
  ITG *p=NULL,i;
  long long hdr[3];
  char *seen=NULL;
  if(perm_dir()==NULL) return NULL;
  perm_name(name,sizeof(name),tag,h,n,nnz);
  if((f=fopen(name,"rb"))==NULL) return NULL;
  if(fread(hdr,sizeof(long long),3,f)!=3||hdr[0]!=n||hdr[1]!=nnz||
     (unsigned long long)hdr[2]!=h){fclose(f);return NULL;}
  NNEW(p,ITG,n);
  if(fread(p,sizeof(ITG),n,f)!=(size_t)n){fclose(f);SFREE(p);return NULL;}
  fclose(f);
  NNEW(seen,char,n);                  /* must be a permutation of 1..n */
  for(i=0;i<n;i++){
    if(p[i]<1||p[i]>n||seen[p[i]-1]){SFREE(seen);SFREE(p);return NULL;}
    seen[p[i]-1]=1;
  }
  SFREE(seen);
  printf(" PARDISO: fill-reducing ordering (%s) reused from %s\n",tag,name);
  return p;
}
static void perm_save(const char *tag,unsigned long long h,ITG n,ITG nnz,
		      const ITG *p)
{
  char name[4096],tmp[4200];
  FILE *f;
  long long hdr[3];
  if(perm_dir()==NULL||p==NULL) return;
  perm_name(name,sizeof(name),tag,h,n,nnz);
  snprintf(tmp,sizeof(tmp),"%s.tmp%d",name,(int)getpid());
  if((f=fopen(tmp,"wb"))==NULL){
    printf(" *WARNING: cannot write %s\n",tmp);return;
  }
  hdr[0]=n;hdr[1]=nnz;hdr[2]=(long long)h;
  if(fwrite(hdr,sizeof(long long),3,f)!=3||
     fwrite(p,sizeof(ITG),n,f)!=(size_t)n){
    fclose(f);remove(tmp);return;
  }
  fclose(f);
  if(rename(tmp,name)==0)              /* atomic: parallel runs are safe */
    printf(" PARDISO: fill-reducing ordering (%s) saved to %s\n",tag,name);
  else remove(tmp);
}

static ITG *nk_perm=NULL;              /* ordering for the next analysis */

static void nk_call(ITG phase,ITG mtype,ITG *neq,ITG nrhs,double *b,
                    double *x,const char *where)
{
  ITG maxfct=1,mnum=1,*perm=(phase==11) ? nk_perm : NULL,msglvl=0,error=0;
  if(nk_single){                       /* float values and vectors */
    float *bf=NULL,*xf=NULL;
    ITG i,n=*neq*nrhs;
    if(phase==33&&b!=NULL){
      NNEW(bf,float,n);NNEW(xf,float,n);
#pragma omp parallel for
      for(i=0;i<n;i++) bf[i]=(float)b[i];
    }
    FORTRAN(pardiso,(nk_pt,&maxfct,&mnum,&mtype,&phase,neq,(double*)nk_auf,
                     nk_pointers,nk_icol,perm,&nrhs,nk_iparm,&msglvl,
                     (double*)bf,(double*)xf,&error));
    if(bf!=NULL){
#pragma omp parallel for
      for(i=0;i<n;i++) x[i]=(double)xf[i];
      SFREE(bf);SFREE(xf);
    }
  }else
  FORTRAN(pardiso,(nk_pt,&maxfct,&mnum,&mtype,&phase,neq,nk_au,
                   nk_pointers,nk_icol,perm,&nrhs,nk_iparm,&msglvl,
                   b,x,&error));
  if(error!=0){
    printf(" *ERROR in %s: pardiso phase %d returned error=%d\n",
           where,(int)phase,(int)error);
    FORTRAN(stop,());
  }
}

static void nk_set_iparm(ITG mtype)
{
  ITG i;
  for(i=0;i<64;i++) nk_iparm[i]=0;
  nk_iparm[0]=1;
  nk_iparm[1]=NK_ORDERING;
  if(mtype==-2){              /* MKL defaults, real symmetric indefinite */
    nk_iparm[9]=8;nk_iparm[10]=0;nk_iparm[12]=0;nk_iparm[20]=1;
    nk_iparm[23]=pardiso_twolevel();
    nk_iparm[27]=nk_single;           /* 1: single precision */
  }else{                      /* MKL defaults, real unsymmetric */
    nk_iparm[9]=13;nk_iparm[10]=1;nk_iparm[12]=1;
  }
  nk_iparm[17]=-1;            /* report nnz(L) */
}

static void nk_release(void)
{
  ITG neq=nk_neq;
  if(!nk_alive) return;
  printf(" PARDISO (NK): symbolic=%ld factorizations=%ld back-solves=%ld"
         " FGMRES ok=%ld (avg %.1f it) failed->refactor=%ld\n"
         "               time: symbolic %.1f s, factor %.1f s, FGMRES %.1f s,"
         " back-solve %.1f s\n",
         n_symbolic,n_factor,n_backsolve,n_kry_ok,
         n_kry_ok>0 ? (double)n_kry_iter/(double)n_kry_ok : 0.,n_kry_fail,
         t_symbolic,t_factor,t_krylov,t_backsolve);
  nk_call(-1,nk_mtype,&neq,1,NULL,NULL,"nk_release");
  SFREE(nk_pointers);SFREE(nk_icol);SFREE(nk_au);
  nk_pointers=NULL;nk_icol=NULL;nk_au=NULL;
  if(nk_auf!=NULL){SFREE(nk_auf);nk_auf=NULL;}
  nk_single=0;
  nk_alive=0;nk_valid=0;nk_stale=0;nk_neq=-1;nk_mtype=0;nk_nnz=-1;
  nk_hash=0ULL;
  n_symbolic=n_factor=n_kry_ok=n_kry_fail=n_kry_iter=n_backsolve=0;
  t_symbolic=t_factor=t_krylov=t_backsolve=0.;
}

static unsigned long long nk_pattern_hash(ITG neq,ITG nnz)
{
  unsigned long long h=1469598103934665603ULL;
  ITG i;
  for(i=0;i<=neq;i++){h^=(unsigned long long)nk_pointers[i];h*=1099511628211ULL;}
  for(i=0;i<nnz;i++){h^=(unsigned long long)nk_icol[i];h*=1099511628211ULL;}
  return h;
}

static void nk_numeric(ITG *neq)
{
  double t0=omp_get_wtime(),dt;
  nk_call(22,nk_mtype,neq,1,NULL,NULL,"pardiso factor (NK)");
  nk_valid=1;nk_stale=0;n_factor++;
  dt=omp_get_wtime()-t0;t_factor+=dt;
  printf(" PARDISO: factorization %.2f s%s\n",dt,
	 nk_single ? " (single precision, CALTOP_MG_ADJ_SINGLE)" : "");
}

/* take over a CSR matrix (1-based; arrays owned by the NK handle from
   now on); symbolic analysis only if the pattern changed; returns 1 if
   it was redone (no usable factor) */
/* analysis (phase 11), with a saved fill-reducing ordering if there is
   one (CALTOP_PARDISO_PERM_DIR) */
static void nk_analysis(ITG *neq)
{
  ITG nnz=nk_pointers[*neq]-1,save=0;
  unsigned long long h=nk_pattern_hash(*neq,nnz)^(unsigned long long)nk_mtype;
  nk_perm=NULL;
  if(perm_dir()!=NULL){
    nk_perm=perm_load("fine",h,*neq,nnz);
    if(nk_perm!=NULL) nk_iparm[4]=1;
    else{NNEW(nk_perm,ITG,*neq);nk_iparm[4]=2;save=1;}
  }
  nk_call(11,nk_mtype,neq,1,NULL,NULL,"pardiso symbolic (NK)");
  if(save) perm_save("fine",h,*neq,nnz,nk_perm);
  nk_iparm[4]=0;
  if(nk_perm!=NULL){SFREE(nk_perm);nk_perm=NULL;}
}
static ITG nk_adopt(ITG *neq,ITG mtype,ITG *p,ITG *c,double *a)
{
  ITG i,nnz;
  double t0;
  unsigned long long h;

  if(nk_pointers!=NULL){SFREE(nk_pointers);SFREE(nk_icol);SFREE(nk_au);}
  nk_pointers=p;nk_icol=c;nk_au=a;

  nnz=nk_pointers[*neq]-1;
  if(nk_auf!=NULL){SFREE(nk_auf);nk_auf=NULL;}
  if(mtype!=-2) nk_single=0;
  if(nk_single){
    NNEW(nk_auf,float,nnz);
#pragma omp parallel for
    for(i=0;i<nnz;i++) nk_auf[i]=(float)nk_au[i];
  }
  h=nk_pattern_hash(*neq,nnz);
  if(nk_alive&&nk_neq==*neq&&nk_mtype==mtype&&nk_nnz==nnz&&nk_hash==h&&
     nk_iparm[27]==nk_single)
    return 0;

  if(nk_alive){
    ITG neq0=nk_neq;
    nk_call(-1,nk_mtype,&neq0,1,NULL,NULL,"nk_store");
  }
  for(i=0;i<64;i++) nk_pt[i]=0;
  nk_set_iparm(mtype);
  nk_mtype=mtype;
  t0=omp_get_wtime();
  nk_analysis(neq);
  t_symbolic+=omp_get_wtime()-t0;n_symbolic++;
  nk_alive=1;nk_valid=0;nk_stale=0;
  nk_neq=*neq;nk_nnz=nnz;nk_hash=h;
  return 1;
}

/* store a new matrix in the NK handle; symbolic analysis only if the
   pattern changed; returns 1 if it was redone (no usable factor) */
static ITG nk_store(double *ad, double *au, double *adb, double *aub,
                double *sigma,ITG *icol, ITG *irow,
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
{
  ITG mtype,*sp=pointers,*si=icolpardiso,*np,*nc;
  double *sa=aupardiso,*na;

  if(factor_alive) pardiso_cleanup(&alive_neq,&alive_symmetryflag);
  pardiso_threads();

  /* pardiso_csr fills the stock arrays; move them to the NK handle */
  mtype=pardiso_csr(ad,au,adb,aub,sigma,icol,irow,neq,nzs,
		    symmetryflag,inputformat,jq,nzs3);
  np=pointers;nc=icolpardiso;na=aupardiso;
  pointers=sp;icolpardiso=si;aupardiso=sa;
  return nk_adopt(neq,mtype,np,nc,na);
}

static double vdot(ITG n,const double *x,const double *y)
{
  double s=0.;
  ITG i;
#pragma omp parallel for reduction(+:s)
  for(i=0;i<n;i++) s+=x[i]*y[i];
  return s;
}

/* y = A x with the stored CSR (upper triangle for mtype -2) */
static void nk_matvec(ITG n,const double *x,double *y)
{
  sparse_matrix_t A;
  struct matrix_descr d;
  sparse_status_t st;
  st=mkl_sparse_d_create_csr(&A,SPARSE_INDEX_BASE_ONE,n,n,nk_pointers,
                             nk_pointers+1,nk_icol,nk_au);
  if(st!=SPARSE_STATUS_SUCCESS){
    printf(" *ERROR in nk_matvec: mkl_sparse_d_create_csr=%d\n",(int)st);
    FORTRAN(stop,());
  }
  if(nk_mtype==-2){
    d.type=SPARSE_MATRIX_TYPE_SYMMETRIC;
    d.mode=SPARSE_FILL_MODE_UPPER;
    d.diag=SPARSE_DIAG_NON_UNIT;
  }else{
    d.type=SPARSE_MATRIX_TYPE_GENERAL;
  }
  st=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,1.0,A,d,x,0.0,y);
  mkl_sparse_destroy(A);
  if(st!=SPARSE_STATUS_SUCCESS){
    printf(" *ERROR in nk_matvec: mkl_sparse_d_mv=%d\n",(int)st);
    FORTRAN(stop,());
  }
}

/*
 * Right-preconditioned flexible GMRES(m) for A x = b, x0 = 0, with the
 * existing factor as M (flexible: MKL may apply iterative refinement in
 * phase 33). Returns the iteration count on ||b-Ax|| <= 10^-L ||b||
 * (true residual), or -1.
 */
static ITG nk_fgmres(ITG *neq,const double *b,double *x,double tol)
{
  ITG n=*neq,m=NK_KRYLOV_M,maxit=NK_KRYLOV_MAXIT,it=0,i,j,k,jl,giveup=0;
  double bnorm,beta,rnorm,tmp,*V=NULL,*Z=NULL,
         *H=NULL,*cs=NULL,*sn=NULL,*g=NULL,*y=NULL,*w=NULL;

  bnorm=sqrt(vdot(n,b,b));
  for(i=0;i<n;i++) x[i]=0.;
  if(bnorm==0.) return 1;

  NNEW(V,double,(size_t)(m+1)*n);
  NNEW(Z,double,(size_t)m*n);
  NNEW(H,double,(m+1)*m);
  NNEW(cs,double,m);NNEW(sn,double,m);NNEW(g,double,m+1);NNEW(y,double,m);
  NNEW(w,double,n);

  memcpy(V,b,sizeof(double)*n);
  beta=bnorm;

  while(1){
    for(i=0;i<n;i++) V[i]/=beta;
    for(i=0;i<=m;i++) g[i]=0.;
    g[0]=beta;
    jl=-1;

    for(j=0;j<m&&it<maxit;j++){
      double *vj=&V[(size_t)j*n],*zj=&Z[(size_t)j*n],*vj1=&V[(size_t)(j+1)*n];
      nk_call(33,nk_mtype,neq,1,vj,zj,"pardiso preconditioner (NK)");
      nk_matvec(n,zj,w);
      for(k=0;k<=j;k++){                      /* modified Gram-Schmidt */
        double *vk=&V[(size_t)k*n],hk=vdot(n,w,vk);
        H[k+j*(m+1)]=hk;
#pragma omp parallel for
        for(i=0;i<n;i++) w[i]-=hk*vk[i];
      }
      tmp=sqrt(vdot(n,w,w));
      H[j+1+j*(m+1)]=tmp;
      if(tmp>0.){
#pragma omp parallel for
        for(i=0;i<n;i++) vj1[i]=w[i]/tmp;
      }
      for(k=0;k<j;k++){                       /* previous Givens */
        double h0=H[k+j*(m+1)],h1=H[k+1+j*(m+1)];
        H[k+j*(m+1)]  = cs[k]*h0+sn[k]*h1;
        H[k+1+j*(m+1)]=-sn[k]*h0+cs[k]*h1;
      }
      {                                       /* new Givens */
        double h0=H[j+j*(m+1)],h1=H[j+1+j*(m+1)],r=sqrt(h0*h0+h1*h1);
        if(r==0.){cs[j]=1.;sn[j]=0.;}else{cs[j]=h0/r;sn[j]=h1/r;}
        H[j+j*(m+1)]=r;H[j+1+j*(m+1)]=0.;
        g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];
      }
      it++;jl=j;
      if(fabs(g[j+1])<=tol*bnorm||tmp==0.) break;
      /* give up early: at the rate reached so far the tolerance needs
         more than maxit iterations (or the residual stagnates) */
      if(it>=NK_GIVEUP_MIN){
        double red=fabs(g[j+1])/bnorm;
        if(red>=1.||(double)it*log(tol)/log(red)>(double)maxit){
          giveup=1;break;
        }
      }
    }

    for(k=jl;k>=0;k--){                       /* x += Z y,  H y = g */
      double sk=g[k];
      for(i=k+1;i<=jl;i++) sk-=H[k+i*(m+1)]*y[i];
      y[k]=(H[k+k*(m+1)]!=0.) ? sk/H[k+k*(m+1)] : 0.;
    }
    for(k=0;k<=jl;k++){
      double *zk=&Z[(size_t)k*n],yk=y[k];
#pragma omp parallel for
      for(i=0;i<n;i++) x[i]+=yk*zk[i];
    }

    nk_matvec(n,x,w);                         /* true residual */
#pragma omp parallel for
    for(i=0;i<n;i++) V[i]=b[i]-w[i];
    rnorm=sqrt(vdot(n,V,V));
    if(rnorm<=tol*bnorm) break;
    if(giveup||it>=maxit||rnorm!=rnorm){it=-1;break;}
    beta=rnorm;                               /* restart */
  }

  SFREE(V);SFREE(Z);SFREE(H);SFREE(cs);SFREE(sn);SFREE(g);SFREE(y);SFREE(w);
  return it;
}

/* solve with the stored matrix: back-solve if the factor is current,
   FGMRES on the stale factor otherwise (refactor if FGMRES fails) */
/* single-precision factor not good enough: analysis and factorization
   again in double precision */
static void nk_to_double(ITG *neq)
{
  ITG i,neq0=nk_neq;
  double t0;
  nk_call(-1,nk_mtype,&neq0,1,NULL,NULL,"nk_to_double");
  nk_single=0;
  if(nk_auf!=NULL){SFREE(nk_auf);nk_auf=NULL;}
  for(i=0;i<64;i++) nk_pt[i]=0;
  nk_set_iparm(nk_mtype);
  t0=omp_get_wtime();
  nk_analysis(neq);
  t_symbolic+=omp_get_wtime()-t0;n_symbolic++;
  nk_numeric(neq);
}
static void nk_solve(double *b,ITG *neq,ITG *symmetryflag,ITG *nrhs,double tol)
{
  ITG i,r,it,mtype=(*symmetryflag==0) ? -2 : 11;
  double *x=NULL,t0,dt;

  if(!nk_alive||!nk_valid||mtype!=nk_mtype||*neq!=nk_neq){
    printf(" *ERROR in pardiso_solve: no matching NK matrix stored\n");
    FORTRAN(stop,());
  }
  NNEW(x,double,*nrhs**neq);
  for(r=0;r<*nrhs;r++){
    double *br=&b[(size_t)r**neq],*xr=&x[(size_t)r**neq];
    if(nk_single){
      t0=omp_get_wtime();
      it=nk_fgmres(neq,br,xr,tol);
      dt=omp_get_wtime()-t0;t_krylov+=dt;
      if(it>0){
        n_kry_ok++;n_kry_iter+=it;
        printf(" PARDISO: single-precision factor, FGMRES %d it to %.1e,"
	       " %.2f s\n",(int)it,tol,dt);
        continue;
      }
      n_kry_fail++;
      printf(" PARDISO: FGMRES on the single-precision factor not converging"
	     " to %.1e (%.2f s) -> double-precision factorization\n",tol,dt);
      nk_to_double(neq);
    }else if(nk_stale){
      t0=omp_get_wtime();
      it=nk_fgmres(neq,br,xr,tol);
      dt=omp_get_wtime()-t0;t_krylov+=dt;
      if(it>0){
        n_kry_ok++;n_kry_iter+=it;
        printf(" PARDISO: old factor reused, FGMRES %d it to %.1e, %.2f s\n",
               (int)it,tol,dt);
        continue;
      }
      n_kry_fail++;
      printf(" PARDISO: FGMRES not converging to %.1e (%.2f s) -> refactor\n",
             tol,dt);
      nk_numeric(neq);
    }
    t0=omp_get_wtime();
    nk_call(33,nk_mtype,neq,1,br,xr,"pardiso_solve (NK)");
    t_backsolve+=omp_get_wtime()-t0;n_backsolve++;
  }
  for(i=0;i<*nrhs**neq;i++) b[i]=x[i];
  SFREE(x);
}

/* Newton step: store the current K_T and solve it (FGMRES on the last
   factor; factor only if there is none or the pattern changed) */
void pardiso_factor_solve(double *ad, double *au, double *adb, double *aub,
                double *sigma,double *b,ITG *icol, ITG *irow,
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3,ITG *nrhs)
{
  if(bddc_active()&&*symmetryflag==0&&*sigma==0.&&*nrhs==1){
    /* CALTOP_SOLVER=bddc: PETSc BDDC; a solve that fails stops the run
       (no fallback: a direct solve of a BDDC-sized problem would not fit) */
    bddc_solve(b,*neq,mg_tol(1),0,ad,au,jq,irow,NULL);
    mg_atol=0.;
    return;
  }
  if(g_useMG&&mg_store(ad,au,adb,aub,sigma,icol,irow,neq,nzs,symmetryflag,
		       inputformat,jq,nzs3)){
    mg_solve(b,neq,nrhs,mg_tol(1),mg_atol,0);
    mg_atol=0.;
    return;
  }
  if(nk_store(ad,au,adb,aub,sigma,icol,irow,neq,nzs,symmetryflag,
	      inputformat,jq,nzs3)||!nk_valid){
    nk_numeric(neq);
  }else{
    nk_stale=1;
  }
  nk_solve(b,neq,symmetryflag,nrhs,NK_TOL_NEWTON);
}

/* store the converged K_T and factor it (symbolic analysis reused);
   pardiso_solve then back-solves the adjoints with it */
void pardiso_store_factor(double *ad, double *au, double *adb, double *aub,
                double *sigma,ITG *icol, ITG *irow,
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
{
  mg_options();
  bd_adj=0;
  if(bddc_active()&&*symmetryflag==0&&*sigma==0.){
    /* CALTOP_SOLVER=bddc: the adjoints use BDDC on this matrix (the
       solver ranks hold it); no copy is kept on rank 0 */
    if(factor_alive) pardiso_cleanup(&alive_neq,&alive_symmetryflag);
    if(mg_alive) mg_release();
    bd_adj=1;bd_adj_ad=ad;bd_adj_au=au;bd_adj_jq=jq;bd_adj_irow=irow;
    return;
  }
  if(g_useMG&&mg_directadj){
    /* CALTOP_MG_DIRECT_ADJOINT: factor the converged K_T (as --NK) */
    if(mg_alive) mg_release();
    nk_single=0;
    printf(" p-multigrid: converged K_T factored with PARDISO for the"
	   " adjoints (CALTOP_MG_DIRECT_ADJOINT%s)\n",
	   nk_single ? ", single precision" : "");
  }else{
    ITG ok;
    mg_storing=1;
    ok=g_useMG&&mg_store(ad,au,adb,aub,sigma,icol,irow,neq,nzs,symmetryflag,
			 inputformat,jq,nzs3);
    mg_storing=0;
    if(ok) return;
  }
  nk_store(ad,au,adb,aub,sigma,icol,irow,neq,nzs,symmetryflag,
	   inputformat,jq,nzs3);
  nk_numeric(neq);
}

/*
 * ----------------------------------------------------------------------
 * p-multigrid preconditioned CG (--MG true, nonlingeo only)
 * ----------------------------------------------------------------------
 * Replaces the fine-level factorizations of the NLGEOM sensitivity runs
 * (same scope as --NK). Two levels on the same mesh:
 *  - fine:   the assembled K_T (C3D10), stored as PARDISO upper CSR
 *  - coarse: Pc = P^T K_T P, P = C3D10 -> C3D4 interpolation (corner
 *            dofs copied, midside dofs = mean of the two edge corners),
 *            i.e. linear displacement fields on the quadratic elements;
 *            factored with PARDISO (separate handle, symbolic analysis
 *            reused while the pattern does not change)
 *  - smoother: Chebyshev of degree CALTOP_MG_DEGREE (default 2) on the
 *            nodal 3x3 block-Jacobi scaled K_T, interval
 *            [0.1,1.1]*lambda_max (lambda_max from 10 Lanczos steps)
 *  - one symmetric V-cycle per CG iteration
 * Newton solves to CALTOP_MG_TOL (default 1e-3), adjoints (pardiso_solve
 * on the stored converged K_T) to 1e-8. If CG fails (no convergence,
 * K_T not positive definite) the matrix is factored with PARDISO in the
 * NK handle and kept for the remaining solves with it; after two
 * failures in a row the run continues with --NK.
 *  - CG breakdown (K_T not positive definite): GMRES(20) with the same
 *    V-cycle for this matrix, then the PARDISO fallback
 *  - CG/GMRES stalled: best iterate accepted if its true residual is
 *    <= max(tol, CALTOP_MG_ACCEPT=1e-5)
 * Environment (optional): CALTOP_MG_DEGREE, CALTOP_MG_SMOOTHER
 * (block|jacobi), CALTOP_MG_TOL, CALTOP_MG_MAXIT, CALTOP_MG_ACCEPT,
 * CALTOP_MG_KRYLOV=gmres, CALTOP_MG_VERBOSE=1.
 * ----------------------------------------------------------------------
 */

/* mesh, pointers into nonlingeo's arrays (pardiso_mg_mesh) */
static ITG *mg_kon=NULL,*mg_ipkon=NULL,*mg_nactdof=NULL,mg_ne=0,mg_nk=0,
           mg_mt=0;
static char *mg_lakon=NULL;
static unsigned long long mg_mesh_hash=0ULL;
static double *mg_co=NULL;                     /* node coordinates (skin smoother) */


/* prolongation P (n x nc, <= 2 entries per row) and R = P^T, 1-based */
static ITG mg_P_n=-1,mg_nc=0,*mg_Pp=NULL,*mg_Pc=NULL,*mg_Rp=NULL,
           *mg_Rc=NULL;
static double *mg_Pv=NULL,*mg_Rv=NULL;
static unsigned long long mg_P_hash=0ULL;
static ITG *mg_cdof=NULL;          /* coarse unknown -> 3*node+dir-1, or -1 */

/* fine matrix: upper CSR from pardiso_csr, MKL handle for products */
static ITG mg_n=0,*mg_ap=NULL,*mg_ai=NULL;
static double *mg_av=NULL;
static sparse_matrix_t mg_A=NULL,mg_Ph=NULL,mg_C=NULL;
static ITG mg_sinceeig=0;
static struct matrix_descr mg_Ad;

/* smoother: inverted node blocks (<= 3x3, 9 slots each) */
static ITG mg_nblk=0,*mg_bptr=NULL,*mg_beq=NULL;
static double *mg_binv=NULL,mg_lmax=0.;

/* coarse matrix and its PARDISO handle */
static long long mg_pt[64];
static ITG mg_iparm[64],*mg_cp=NULL,*mg_ci=NULL,mg_calive=0,mg_cnnz=-1;
static double *mg_cv=NULL;
static unsigned long long mg_chash=0ULL;

/* state */
static ITG mg_directadj=0,mg_useskin=1,mg_nskin=0,*mg_skin=NULL,
           mg_cfresh=-1,mg_cage=0,mg_cmaxage=0,mg_csingle=0,mg_clast=-1,
           mg_cstale=0,mg_directfirst=0,mg_forcegm=0,mg_fallback=0,mg_indef=0,mg_failrow=0,mg_off=0,mg_degree=2,mg_block=1,
           mg_maxit=500,mg_verbose=0,mg_opts=0;
static double mg_tolnewton=1.e-3,mg_accept=1.e-5,mg_toladj=MG_TOL_SOLVE_;
static float *mg_cvf=NULL,*mg_cxf=NULL;
static double *mg_w[7]={NULL,NULL,NULL,NULL,NULL,NULL,NULL};
static ITG mg_wn=0;
static long mg_ncfac=0,mg_nmat=0,mg_nsolve=0,mg_nit=0,mg_nfall=0,mg_ncsym=0;
static double mg_tsetup=0.,mg_tgal=0.,mg_tcfac=0.,mg_tcg=0.,mg_tcsr=0.,
              mg_tP=0.,mg_tmkl=0.,mg_tblk=0.,mg_teig=0.;

static void mg_options(void)
{
  char *env;
  if(mg_opts) return;
  mg_opts=1;
  if((env=getenv("CALTOP_MG_DEGREE"))!=NULL&&atoi(env)>0) mg_degree=atoi(env);
  if((env=getenv("CALTOP_MG_TOL"))!=NULL&&atof(env)>0.) mg_tolnewton=atof(env);
  if((env=getenv("CALTOP_MG_MAXIT"))!=NULL&&atoi(env)>0) mg_maxit=atoi(env);
  if((env=getenv("CALTOP_MG_VERBOSE"))!=NULL&&atoi(env)>0) mg_verbose=1;
  if((env=getenv("CALTOP_MG_DIRECT_ADJOINT"))!=NULL&&atoi(env)>0)
    mg_directadj=1;
  if((env=getenv("CALTOP_MG_TOL_ADJ"))!=NULL&&atof(env)>0.)
    mg_toladj=atof(env);
  {                                   /* every CALTOP_ setting, so a stray
					 one from an earlier test shows up */
    extern char **environ;
    char **ep;
    const char *o=getenv("OMP_NUM_THREADS"),*m=getenv("MKL_NUM_THREADS");
    printf(" environment: OMP_NUM_THREADS=%s MKL_NUM_THREADS=%s\n",
	   o ? o : "(unset)",m ? m : "(unset)");
    for(ep=environ;*ep!=NULL;ep++)
      if(strncmp(*ep,"CALTOP_",7)==0&&strncmp(*ep,"CALTOP_BDDC_CHILD=",18)!=0&&
	 strncmp(*ep,"CALTOP_BDDC_OMP0=",17)!=0&&strncmp(*ep,"CALTOP_BDDC_MKL0=",17)!=0)
	printf(" environment: %s\n",*ep);
  }
  printf(" p-multigrid: Chebyshev(%d) on %s Jacobi, CG tol %.1e (Newton),"
	 " %.1e (adjoints), max %d it\n\n",(int)mg_degree,
	 mg_block ? "3x3 node-block" : "point",mg_tolnewton,mg_toladj,
	 (int)mg_maxit);
  if(mg_csingle) printf(" p-multigrid: coarse level factored in single"
			" precision\n\n");
}

static double mg_tol(int newton)
{
  mg_options();
  return newton ? mg_tolnewton : mg_toladj;
}

/* nonlingeo passes the residual-force level that CalculiX accepts for
   convergence; the next Newton solve stops once max|r_i| <= atol */
void pardiso_mg_atol(double atol)
{
  mg_atol=(atol>0.) ? atol : 0.;
}

/* passive (skin) elements, 1-based element numbers (skinElementList.nam):
   their midside nodes stay coarse unknowns, so the coarse level keeps the
   quadratic interpolation in the thin skin (no C3D4 locking there) */
void pardiso_mg_skin(int *ids,int n)
{
  ITG i;
  if(mg_skin!=NULL){SFREE(mg_skin);mg_skin=NULL;}
  mg_nskin=0;
  if(ids==NULL||n<=0) return;
  NNEW(mg_skin,ITG,n);
  for(i=0;i<n;i++) mg_skin[i]=ids[i];
  mg_nskin=n;
  mg_P_n=-1;                                   /* rebuild P */
}

/* node coordinates co(3,nk), only used by CALTOP_MG_DUMP */
void pardiso_mg_coords(double *co)
{
  mg_co=co;
}

/* converged displacement vold(0:mi[1],nk), copied in equation numbering:
   starting guess for the adjoint CG solves (the compliance adjoint
   K_T lambda = P has lambda close to u when the response is mildly
   nonlinear). Needs pardiso_mg_mesh first. */
static double *mg_guess=NULL;
static ITG mg_guessn=0;
void pardiso_mg_guess(double *vold)
{
  ITG k,d,e,n=0;
  if(mg_guess!=NULL){SFREE(mg_guess);mg_guess=NULL;}
  mg_guessn=0;
  if(vold==NULL||mg_nactdof==NULL) return;
  for(k=0;k<mg_nk;k++) for(d=1;d<=3;d++){
    e=mg_nactdof[mg_mt*k+d];
    if(e>n) n=e;
  }
  if(n<=0) return;
  NNEW(mg_guess,double,n);
  for(k=0;k<mg_nk;k++) for(d=1;d<=3;d++){
    e=mg_nactdof[mg_mt*k+d];
    if(e>0) mg_guess[e-1]=vold[mg_mt*k+d];
  }
  mg_guessn=n;
}

/* nonlingeo passes its mesh before each solve; only pointers are kept */
void pardiso_mg_mesh(ITG *kon,ITG *ipkon,char *lakon,ITG *ne,ITG *nactdof,
		     ITG *mi,ITG *nk)
{
  unsigned long long h=1469598103934665603ULL;
  ITG i,n;
  mg_kon=kon;mg_ipkon=ipkon;mg_lakon=lakon;mg_nactdof=nactdof;
  mg_ne=*ne;mg_nk=*nk;mg_mt=mi[1]+1;
  n=mg_mt*mg_nk;
  for(i=0;i<n;i++){h^=(unsigned long long)nactdof[i];h*=1099511628211ULL;}
  h^=(unsigned long long)mg_ne;h*=1099511628211ULL;
  mg_mesh_hash=h;
}

/* build P and R for n fine equations; returns 0 if impossible */
static ITG mg_build_P(ITG n)
{
  static const int ed[6][2]={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
  ITG i,j,k,e,d,mt=mg_mt,nk=mg_nk,*eqnode=NULL,*eqdir=NULL,*ma=NULL,
      *mb=NULL,*cidx=NULL,nmid=0,nc=0,nnzP=0,*bidx=NULL,nbub=0,nshort=0,
      skin_tang=0,skin_full=0,skin_geo=0;
  double skin_long=3.;
  char *flag=NULL;   /* 1 corner, 2 midside, 3 conflicting, 4 skin midside
			kept, 5 skin midside with a normal bubble */
  double *nrm=NULL;  /* skin normal per bubble node (CALTOP_MG_SKIN=2) */

  if(mg_kon==NULL) return 0;
  if(mg_Pp!=NULL){SFREE(mg_Pp);SFREE(mg_Pc);SFREE(mg_Pv);
                  SFREE(mg_Rp);SFREE(mg_Rc);SFREE(mg_Rv);}

  NNEW(eqnode,ITG,n);NNEW(eqdir,ITG,n);
  for(i=0;i<n;i++) eqnode[i]=-1;
  for(k=0;k<nk;k++){
    for(d=1;d<=3;d++){
      e=mg_nactdof[mt*k+d]-1;
      if(e>=0&&e<n){eqnode[e]=k;eqdir[e]=d;}
    }
  }

  NNEW(flag,char,nk);NNEW(ma,ITG,nk);NNEW(mb,ITG,nk);
  for(e=0;e<mg_ne;e++){
    ITG idx=mg_ipkon[e],nd[10];
    if(idx<0) continue;
    if(strncmp(&mg_lakon[8*e],"C3D10",5)!=0) continue;
    for(j=0;j<10;j++) nd[j]=mg_kon[idx+j]-1;
    for(j=0;j<4;j++){
      if(flag[nd[j]]==2||flag[nd[j]]==3) flag[nd[j]]=3; else flag[nd[j]]=1;
    }
    for(j=0;j<6;j++){
      ITG m=nd[4+j],a=nd[ed[j][0]],b=nd[ed[j][1]];
      if(flag[m]==0){flag[m]=2;ma[m]=a;mb[m]=b;}
      else if(flag[m]==2){
	if(!((ma[m]==a&&mb[m]==b)||(ma[m]==b&&mb[m]==a))) flag[m]=3;
      }else flag[m]=3;
    }
  }

  /* skin: midside nodes of passive elements stay coarse (flag 4) */
  mg_options();
  /* all three components of the midside nodes on long skin edges (length
     > 3 x local skin thickness) stay coarse unknowns (flag 4); midside
     nodes on the short, through-thickness edges are interpolated linearly
     (flag 2) */
  if(mg_useskin&&mg_nskin>0){
    skin_full=(mg_useskin==1);
    skin_geo=(mg_co!=NULL);
    if(!skin_geo)
      printf(" *WARNING in p-multigrid: no node coordinates, all skin edges"
	     " kept quadratic\n");
  }
  if(mg_useskin&&mg_nskin>0&&skin_geo){
    NNEW(nrm,double,3*nk);
    for(i=0;i<mg_nskin;i++){
      ITG idx,q;
      double x[4][3],c[3]={0.,0.,0.},C[9]={0.},v[3][3],nv[3],h,pmin,pmax,
	bestc=-1.;
      e=mg_skin[i]-1;
      if(e<0||e>=mg_ne) continue;
      idx=mg_ipkon[e];
      if(idx<0||strncmp(&mg_lakon[8*e],"C3D10",5)!=0) continue;
      for(j=0;j<4;j++) for(q=0;q<3;q++){
	x[j][q]=mg_co[3*(mg_kon[idx+j]-1)+q];c[q]+=0.25*x[j][q];
      }
      for(j=0;j<4;j++){                       /* covariance of the corners */
	double y[3];
	for(q=0;q<3;q++) y[q]=x[j][q]-c[q];
	for(q=0;q<9;q++) C[q]+=y[q/3]*y[q%3];
      }
      /* C^2 e_k lie in the plane of the two large principal directions
	 (the thickness direction is damped by (h/L)^4); the normal is the
	 cross product of the best-conditioned pair */
      for(k=0;k<3;k++){
	double t[3],u[3];
	for(q=0;q<3;q++) t[q]=C[3*q+k];
	for(q=0;q<3;q++) u[q]=C[3*q]*t[0]+C[3*q+1]*t[1]+C[3*q+2]*t[2];
	for(q=0;q<3;q++) v[k][q]=u[q];
      }
      nv[0]=nv[1]=nv[2]=0.;
      for(k=0;k<3;k++){
	ITG l=(k+1)%3;
	double cr[3],cn;
	cr[0]=v[k][1]*v[l][2]-v[k][2]*v[l][1];
	cr[1]=v[k][2]*v[l][0]-v[k][0]*v[l][2];
	cr[2]=v[k][0]*v[l][1]-v[k][1]*v[l][0];
	cn=sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
	if(cn>bestc){bestc=cn;for(q=0;q<3;q++) nv[q]=cr[q];}
      }
      if(bestc<=0.) continue;
      for(q=0;q<3;q++) nv[q]/=bestc;
      pmin=1.e300;pmax=-1.e300;               /* thickness along n */
      for(j=0;j<4;j++){
	double pr=x[j][0]*nv[0]+x[j][1]*nv[1]+x[j][2]*nv[2];
	if(pr<pmin) pmin=pr;
	if(pr>pmax) pmax=pr;
      }
      h=pmax-pmin;
      for(j=0;j<6;j++){
	ITG m=mg_kon[idx+4+j]-1,a=ed[j][0],b=ed[j][1];
	double len=0.,s0;
	if(flag[m]!=2&&flag[m]!=5) continue;
	for(q=0;q<3;q++) len+=(x[a][q]-x[b][q])*(x[a][q]-x[b][q]);
	len=sqrt(len);
	if(len>skin_long*h){
	  if(skin_full){flag[m]=4;continue;}  /* all 3 components kept */
	  flag[m]=5;
	  s0=nrm[3*m]*nv[0]+nrm[3*m+1]*nv[1]+nrm[3*m+2]*nv[2];
	  for(q=0;q<3;q++) nrm[3*m+q]+=(s0<0.) ? -nv[q] : nv[q];
	}
      }
    }
    for(k=0;k<nk;k++){
      double nn;
      if(flag[k]!=5) continue;
      nn=sqrt(nrm[3*k]*nrm[3*k]+nrm[3*k+1]*nrm[3*k+1]+nrm[3*k+2]*nrm[3*k+2]);
      if(nn>0.){nrm[3*k]/=nn;nrm[3*k+1]/=nn;nrm[3*k+2]/=nn;}
      else flag[k]=2;
    }
  }else if(mg_useskin&&mg_nskin>0){         /* no coordinates: v5 skin */
    for(i=0;i<mg_nskin;i++){
      ITG idx;
      e=mg_skin[i]-1;
      if(e<0||e>=mg_ne) continue;
      idx=mg_ipkon[e];
      if(idx<0||strncmp(&mg_lakon[8*e],"C3D10",5)!=0) continue;
      for(j=4;j<10;j++){
	ITG m=mg_kon[idx+j]-1;
	if(flag[m]==2) flag[m]=4;
      }
    }
  }

  /* coarse unknowns: every equation not on a (consistent) midside node */
  NNEW(cidx,ITG,n);
  {
    ITG nskeq=0;
    for(i=0;i<n;i++){
      k=eqnode[i];
      if(k>=0&&(flag[k]==2||flag[k]==5)){cidx[i]=-1;nmid++;}
      else{cidx[i]=nc++;if(k>=0&&flag[k]==4) nskeq++;}
    }
    if(mg_useskin&&mg_nskin>0&&skin_geo){     /* distinct short-edge nodes */
      char *seen=NULL;
      NNEW(seen,char,nk);
      for(i=0;i<mg_nskin;i++){
	ITG idx;
	e=mg_skin[i]-1;
	if(e<0||e>=mg_ne||(idx=mg_ipkon[e])<0) continue;
	if(strncmp(&mg_lakon[8*e],"C3D10",5)!=0) continue;
	for(j=4;j<10;j++){
	  ITG m=mg_kon[idx+j]-1;
	  if(flag[m]==2&&!seen[m]){seen[m]=1;nshort++;}
	}
      }
      SFREE(seen);
    }
    if(mg_useskin&&mg_nskin>0&&skin_geo&&!skin_full){ /* bubbles */
      NNEW(bidx,ITG,nk);
      for(k=0;k<nk;k++) bidx[k]=-1;
      for(i=0;i<n;i++){
	k=eqnode[i];
	if(k>=0&&flag[k]==5&&bidx[k]<0){
	  bidx[k]=nc;nc+=skin_tang ? 2 : 1;nbub+=skin_tang ? 2 : 1;
	}
      }
      printf(" p-multigrid: skin (%d passive elements), %s bubbles"
	     " (CALTOP_MG_SKIN=2, experimental): %d bubble unknowns on long"
	     " edges, %d short-edge midside nodes linear\n",(int)mg_nskin,
	     skin_tang ? "tangential" : "normal",(int)nbub,(int)nshort);
    }else if(mg_useskin&&mg_nskin>0&&skin_geo)
      printf(" p-multigrid: skin (%d passive elements) kept quadratic on long"
	     " edges: %d skin midside equations, %d short-edge midside nodes"
	     " linear\n",(int)mg_nskin,(int)nskeq,(int)nshort);
    else if(mg_useskin&&mg_nskin>0)
      printf(" p-multigrid: skin (%d passive elements) kept quadratic on the"
	     " coarse level: %d skin midside equations\n",(int)mg_nskin,
	     (int)nskeq);
  }
  if(nmid==0||nc==0){
    SFREE(eqnode);SFREE(eqdir);SFREE(flag);SFREE(ma);SFREE(mb);SFREE(cidx);
    if(nrm!=NULL) SFREE(nrm);
    if(bidx!=NULL) SFREE(bidx);
    return 0;
  }
  if(mg_cdof!=NULL) SFREE(mg_cdof);
  NNEW(mg_cdof,ITG,nc);
  for(i=0;i<nc;i++) mg_cdof[i]=-1;
  for(i=0;i<n;i++){
    if(cidx[i]<0) continue;
    mg_cdof[cidx[i]]=(eqnode[i]>=0) ? 3*eqnode[i]+eqdir[i]-1 : -1;
  }

  /* P rows */
  NNEW(mg_Pp,ITG,n+1);NNEW(mg_Pc,ITG,4*n);NNEW(mg_Pv,double,4*n);
  mg_Pp[0]=1;
  for(i=0;i<n;i++){
    ITG c[4],m=0,a,b;
    double w[4];
    if(cidx[i]>=0){c[m]=cidx[i];w[m++]=1.;}
    else{
      k=eqnode[i];d=eqdir[i];
      e=mg_nactdof[mt*ma[k]+d]-1;
      if(e>=0&&e<n&&cidx[e]>=0){c[m]=cidx[e];w[m++]=0.5;}
      e=mg_nactdof[mt*mb[k]+d]-1;
      if(e>=0&&e<n&&cidx[e]>=0){c[m]=cidx[e];w[m++]=0.5;}
      if(flag[k]==5&&bidx!=NULL&&bidx[k]>=0){
	if(!skin_tang){
	  if(fabs(nrm[3*k+d-1])>1.e-12){c[m]=bidx[k];w[m++]=nrm[3*k+d-1];}
	}else{                         /* two tangents t1, t2 orthogonal to n */
	  double *nv=&nrm[3*k],ax[3]={0.,0.,0.},t1[3],t2[3],tn;
	  ITG q=0;
	  if(fabs(nv[1])<fabs(nv[q])) q=1;
	  if(fabs(nv[2])<fabs(nv[q])) q=2;
	  ax[q]=1.;
	  t1[0]=nv[1]*ax[2]-nv[2]*ax[1];t1[1]=nv[2]*ax[0]-nv[0]*ax[2];
	  t1[2]=nv[0]*ax[1]-nv[1]*ax[0];
	  tn=sqrt(t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2]);
	  t1[0]/=tn;t1[1]/=tn;t1[2]/=tn;
	  t2[0]=nv[1]*t1[2]-nv[2]*t1[1];t2[1]=nv[2]*t1[0]-nv[0]*t1[2];
	  t2[2]=nv[0]*t1[1]-nv[1]*t1[0];
	  if(fabs(t1[d-1])>1.e-12){c[m]=bidx[k];w[m++]=t1[d-1];}
	  if(fabs(t2[d-1])>1.e-12){c[m]=bidx[k]+1;w[m++]=t2[d-1];}
	}
      }
      for(a=1;a<m;a++){                       /* sort by column */
	ITG cc=c[a];double ww=w[a];
	for(b=a-1;b>=0&&c[b]>cc;b--){c[b+1]=c[b];w[b+1]=w[b];}
	c[b+1]=cc;w[b+1]=ww;
      }
    }
    for(j=0;j<m;j++){
      mg_Pc[nnzP]=c[j]+1;
      mg_Pv[nnzP]=w[j];
      nnzP++;
    }
    mg_Pp[i+1]=nnzP+1;
  }
  if(nrm!=NULL) SFREE(nrm);
  if(bidx!=NULL) SFREE(bidx);

  /* R = P^T */
  NNEW(mg_Rp,ITG,nc+1);NNEW(mg_Rc,ITG,nnzP);NNEW(mg_Rv,double,nnzP);
  for(j=0;j<nnzP;j++) mg_Rp[mg_Pc[j]]++;
  mg_Rp[0]=1;
  for(i=0;i<nc;i++) mg_Rp[i+1]+=mg_Rp[i];
  {
    ITG *pos=NULL;
    NNEW(pos,ITG,nc);
    for(i=0;i<nc;i++) pos[i]=mg_Rp[i]-1;
    for(i=0;i<n;i++){
      for(j=mg_Pp[i]-1;j<mg_Pp[i+1]-1;j++){
	ITG c=mg_Pc[j]-1;
	mg_Rc[pos[c]]=i+1;mg_Rv[pos[c]]=mg_Pv[j];pos[c]++;
      }
    }
    SFREE(pos);
  }

  /* smoother blocks: the equations of one node */
  if(mg_bptr!=NULL){SFREE(mg_bptr);SFREE(mg_beq);}
  NNEW(mg_bptr,ITG,n+1);NNEW(mg_beq,ITG,n);
  {
    ITG nb=0,pe=0;
    char *done=NULL;
    NNEW(done,char,n);
    for(k=0;k<nk;k++){
      ITG s=0;
      for(d=1;d<=3;d++){
	e=mg_nactdof[mt*k+d]-1;
	if(e>=0&&e<n&&!done[e]){
	  if(!mg_block&&s>0){mg_bptr[++nb]=pe;s=0;}
	  mg_beq[pe++]=e;done[e]=1;s++;
	}
      }
      if(s>0) mg_bptr[++nb]=pe;
    }
    for(i=0;i<n;i++) if(!done[i]){mg_beq[pe++]=i;mg_bptr[++nb]=pe;}
    mg_nblk=nb;
    SFREE(done);
  }
  if(mg_binv!=NULL) SFREE(mg_binv);
  NNEW(mg_binv,double,9*mg_nblk);

  printf(" p-multigrid: %d fine equations, %d coarse (%.1f x fewer),"
	 " %d midside equations\n",(int)n,(int)nc,(double)n/(double)nc,
	 (int)nmid);

  mg_nc=nc;mg_P_n=n;mg_P_hash=mg_mesh_hash;
  SFREE(eqnode);SFREE(eqdir);SFREE(flag);SFREE(ma);SFREE(mb);SFREE(cidx);
  return 1;
}

/* A(i,j) from the upper CSR (0-based i,j) */
static double mg_entry(ITG i,ITG j)
{
  ITG k,r=(i<j)?i:j,c=((i<j)?j:i)+1;
  for(k=mg_ap[r]-1;k<mg_ap[r+1]-1;k++) if(mg_ai[k]==c) return mg_av[k];
  return 0.;
}

/* invert the node blocks (Cholesky; diagonal if a block is not SPD) */
static void mg_blocks(void)
{
  ITG b;
#pragma omp parallel for
  for(b=0;b<mg_nblk;b++){
    ITG s=mg_bptr[b+1]-mg_bptr[b],*q=&mg_beq[mg_bptr[b]],i,j,k,ok=1;
    double A[9],L[9],*inv=&mg_binv[9*b];
    for(i=0;i<s;i++) for(j=0;j<s;j++) A[i*3+j]=mg_entry(q[i],q[j]);
    for(i=0;i<9;i++){L[i]=0.;inv[i]=0.;}
    for(j=0;j<s&&ok;j++){                     /* L L^T = A */
      double t=A[j*3+j];
      for(k=0;k<j;k++) t-=L[j*3+k]*L[j*3+k];
      if(!(t>0.)){ok=0;break;}
      L[j*3+j]=sqrt(t);
      for(i=j+1;i<s;i++){
	double u=A[i*3+j];
	for(k=0;k<j;k++) u-=L[i*3+k]*L[j*3+k];
	L[i*3+j]=u/L[j*3+j];
      }
    }
    if(ok){                                   /* inv = L^-T L^-1 */
      double Li[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
      for(j=0;j<s;j++){
	Li[j*3+j]=1./L[j*3+j];
	for(i=j+1;i<s;i++){
	  double u=0.;
	  for(k=j;k<i;k++) u-=L[i*3+k]*Li[k*3+j];
	  Li[i*3+j]=u/L[i*3+i];
	}
      }
      for(i=0;i<s;i++) for(j=0;j<s;j++){
	double u=0.;
	for(k=(i>j?i:j);k<s;k++) u+=Li[k*3+i]*Li[k*3+j];
	inv[i*3+j]=u;
      }
    }else{
      for(i=0;i<s;i++){
	double a=fabs(A[i*3+i]);
	inv[i*3+i]=(a>0.) ? 1./a : 1.;
      }
    }
  }
}

/* y = B^-1 x (node-block Jacobi) */
static void mg_binv_apply(const double *x,double *y)
{
  ITG b;
#pragma omp parallel for
  for(b=0;b<mg_nblk;b++){
    ITG s=mg_bptr[b+1]-mg_bptr[b],*q=&mg_beq[mg_bptr[b]],i,j;
    const double *inv=&mg_binv[9*b];
    double xv[3];
    for(i=0;i<s;i++) xv[i]=x[q[i]];
    for(i=0;i<s;i++){
      double u=0.;
      for(j=0;j<s;j++) u+=inv[i*3+j]*xv[j];
      y[q[i]]=u;
    }
  }
}

/* y = K x with K stored as its upper triangle (diagonal first in each
   row). MKL's product for this storage scales poorly with threads (about
   0.26 s per product on the wing with 20 threads, ~5 GB/s), since every
   off-diagonal entry also updates y[j]. Here each thread takes a slice of
   rows with about equal nonzeros and adds the mirrored part into its own
   buffer (it only touches rows >= its first row); the buffers are summed
   afterwards. Scratch: threads x n doubles. CALTOP_MG_MKLMV=1 uses MKL */
static ITG mg_mvnt=0,mg_mvn=-1,*mg_mvpart=NULL,mg_mkl_mv=0;
static double *mg_mvbuf=NULL,mg_tmv=0.,mg_tcs=0.,mg_tmvcg=0.;
static long mg_nmv=0,mg_ncs=0,mg_nmvcg=0;

static void mg_mv_release(void)
{
  if(mg_mvpart!=NULL){SFREE(mg_mvpart);mg_mvpart=NULL;}
  if(mg_mvbuf!=NULL){SFREE(mg_mvbuf);mg_mvbuf=NULL;}
  mg_mvn=-1;mg_mvnt=0;
}

/* row slices with about equal nonzeros (new matrix pattern) */
static void mg_mv_setup(void)
{
  ITG t,i,nt=omp_get_max_threads(),n=mg_n;
  double nnz=(double)(mg_ap[n]-1);
  mg_mv_release();
  if(nt<1) nt=1;
  NNEW(mg_mvpart,ITG,nt+1);
  NNEW(mg_mvbuf,double,(size_t)nt*n);
  i=0;
  for(t=0;t<nt;t++){
    double target=nnz*(double)(t+1)/(double)nt;
    mg_mvpart[t]=i;
    while(i<n&&(double)(mg_ap[i+1]-1)<=target) i++;
    if(t==nt-1) i=n;
  }
  mg_mvpart[nt]=n;
  mg_mvnt=nt;mg_mvn=n;
}

static void mg_matvec(const double *x,double *y)
{
  ITG n=mg_n,nt=mg_mvnt,ok=1;
  double t0=omp_get_wtime();
  if(!mg_mkl_mv&&mg_mvn==n&&nt>0){
#pragma omp parallel num_threads(nt)
    {
      ITG t=omp_get_thread_num(),r0,r1,i,k,j;
      double *b;
      if(omp_get_num_threads()!=nt){
#pragma omp single
	ok=0;
      }else{
	r0=mg_mvpart[t];r1=mg_mvpart[t+1];b=mg_mvbuf+(size_t)t*n;
	memset(b+r0,0,sizeof(double)*(size_t)(n-r0));
	for(i=r0;i<r1;i++){
	  double xi=x[i],s=0.;
	  for(k=mg_ap[i]-1;k<mg_ap[i+1]-1;k++){
	    double a=mg_av[k];
	    j=mg_ai[k]-1;
	    s+=a*x[j];
	    if(j!=i) b[j]+=a*xi;
	  }
	  b[i]+=s;
	}
#pragma omp barrier
#pragma omp for schedule(static)
	for(i=0;i<n;i++){
	  double u=0.;
	  ITG tt;
	  for(tt=0;tt<nt&&mg_mvpart[tt]<=i;tt++) u+=mg_mvbuf[(size_t)tt*n+i];
	  y[i]=u;
	}
      }
    }
  }else ok=0;
  if(!ok){
    sparse_status_t st=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,1.0,
				       mg_A,mg_Ad,x,0.0,y);
    if(st!=SPARSE_STATUS_SUCCESS){
      printf(" *ERROR in p-multigrid: mkl_sparse_d_mv=%d\n",(int)st);
      FORTRAN(stop,());
    }
  }
  mg_tmv+=omp_get_wtime()-t0;mg_nmv++;
}

static double mg_dot(const double *x,const double *y,ITG n)
{
  double s=0.;
  ITG i;
#pragma omp parallel for reduction(+:s)
  for(i=0;i<n;i++) s+=x[i]*y[i];
  return s;
}

/* largest eigenvalue of B^-1 A: 10 steps of Lanczos in the B inner
   product (valid for indefinite A), bisection on the tridiagonal matrix */
static double mg_eigmax(void)
{
  ITG n=mg_n,i,k,m=10;
  double *q=mg_w[0],*p=mg_w[1],*po=mg_w[2],*r=mg_w[3],*z=mg_w[4],
         dg[10],od[10],be,beo=0.,lo,hi,mid;
  unsigned long long s=88172645463325252ULL;

  for(i=0;i<n;i++){                           /* reproducible random start */
    s^=s<<13;s^=s>>7;s^=s<<17;
    r[i]=(double)(s%2000001ULL)/1.e6-1.;po[i]=0.;
  }
  mg_binv_apply(r,z);
  be=sqrt(fabs(mg_dot(r,z,n)));
  if(!(be>0.)) return 2.;
#pragma omp parallel for
  for(i=0;i<n;i++){q[i]=z[i]/be;p[i]=r[i]/be;}
  for(k=0;k<m;k++){
    double al,bn;
    mg_matvec(q,r);
    al=mg_dot(q,r,n);
#pragma omp parallel for
    for(i=0;i<n;i++) r[i]-=al*p[i]+beo*po[i];
    mg_binv_apply(r,z);
    bn=sqrt(fabs(mg_dot(r,z,n)));
    dg[k]=al;od[k]=bn;
    if(!(bn>0.)||al!=al){k++;break;}
#pragma omp parallel for
    for(i=0;i<n;i++){po[i]=p[i];p[i]=r[i]/bn;q[i]=z[i]/bn;}
    beo=bn;
  }
  m=k;
  if(m==0) return 2.;
  /* Gershgorin bounds, then bisection on the Sturm count */
  hi=-1.e300;lo=1.e300;
  for(k=0;k<m;k++){
    double g=((k>0)?fabs(od[k-1]):0.)+((k<m-1)?fabs(od[k]):0.);
    if(dg[k]+g>hi) hi=dg[k]+g;
    if(dg[k]-g<lo) lo=dg[k]-g;
  }
  for(i=0;i<60;i++){
    ITG cnt=0;
    double t=1.;
    mid=0.5*(lo+hi);
    for(k=0;k<m;k++){
      t=dg[k]-mid-((k>0)?od[k-1]*od[k-1]/t:0.);
      if(t==0.) t=1.e-300;
      if(t<0.) cnt++;
    }
    if(cnt==m) hi=mid; else lo=mid;
  }
  return (hi>0.) ? hi : 2.;
}

/* Chebyshev smoothing of A x = b on [0.1,1.1]*lambda_max of B^-1 A;
   x_zero: x is 0 on entry (saves one product) */
static void mg_cheb(const double *b,double *x,ITG x_zero)
{
  ITG n=mg_n,i,k;
  double *r=mg_w[4],*d=mg_w[5],*w=mg_w[6],
         lmin=0.1*mg_lmax,lmax=1.1*mg_lmax,
         theta=0.5*(lmax+lmin),delta=0.5*(lmax-lmin),
         sigma=theta/delta,rho=1./sigma,rho1;

  if(x_zero) memcpy(r,b,sizeof(double)*n);
  else{
    mg_matvec(x,w);
#pragma omp parallel for
    for(i=0;i<n;i++) r[i]=b[i]-w[i];
  }
  mg_binv_apply(r,d);
#pragma omp parallel for
  for(i=0;i<n;i++) d[i]/=theta;
  for(k=0;k<mg_degree;k++){
#pragma omp parallel for
    for(i=0;i<n;i++) x[i]+=d[i];
    if(k==mg_degree-1) break;
    mg_matvec(d,w);
#pragma omp parallel for
    for(i=0;i<n;i++) r[i]-=w[i];
    rho1=1./(2.*sigma-rho);
    mg_binv_apply(r,w);
#pragma omp parallel for
    for(i=0;i<n;i++) d[i]=rho1*rho*d[i]+2.*rho1/delta*w[i];
    rho=rho1;
  }
}

/* z = M^-1 r: one symmetric V-cycle */
static void mg_vcycle(const double *r,double *z,double *rc)
{
  ITG n=mg_n,nc=mg_nc,i,one=1,maxfct=1,mnum=1,phase=33,mtype=-2,msglvl=0,
      error=0,*perm=NULL;
  double *t=mg_w[6],*ec=rc+nc,tcs0;

  for(i=0;i<n;i++) z[i]=0.;
  mg_cheb(r,z,1);                             /* pre-smoothing */
  mg_matvec(z,t);
#pragma omp parallel for
  for(i=0;i<n;i++) t[i]=r[i]-t[i];
  tcs0=omp_get_wtime();
#pragma omp parallel for
  for(i=0;i<nc;i++){                          /* rc = P^T t */
    ITG j;
    double u=0.;
    for(j=mg_Rp[i]-1;j<mg_Rp[i+1]-1;j++) u+=mg_Rv[j]*t[mg_Rc[j]-1];
    rc[i]=u;
  }
  if(mg_csingle){                             /* single-precision factor */
    float *bf=mg_cxf,*xf=mg_cxf+nc;
#pragma omp parallel for
    for(i=0;i<nc;i++) bf[i]=(float)rc[i];
    FORTRAN(pardiso,(mg_pt,&maxfct,&mnum,&mtype,&phase,&mg_nc,
		     (double*)mg_cvf,mg_cp,mg_ci,perm,&one,mg_iparm,&msglvl,
		     (double*)bf,(double*)xf,&error));
#pragma omp parallel for
    for(i=0;i<nc;i++) ec[i]=(double)xf[i];
  }else
    FORTRAN(pardiso,(mg_pt,&maxfct,&mnum,&mtype,&phase,&mg_nc,mg_cv,mg_cp,
		     mg_ci,perm,&one,mg_iparm,&msglvl,rc,ec,&error));
  if(error!=0){
    printf(" *ERROR in p-multigrid: coarse solve error=%d\n",(int)error);
    FORTRAN(stop,());
  }
#pragma omp parallel for
  for(i=0;i<n;i++){                           /* z += P ec */
    ITG j;
    double u=0.;
    for(j=mg_Pp[i]-1;j<mg_Pp[i+1]-1;j++) u+=mg_Pv[j]*ec[mg_Pc[j]-1];
    z[i]+=u;
  }
  mg_tcs+=omp_get_wtime()-tcs0;mg_ncs++;
  mg_cheb(r,z,0);                             /* post-smoothing */
}

/* coarse matrix P^T A P (MKL sypr, upper triangle), PARDISO factor.
   reuse: same fine pattern as the previous matrix, so the product only
   recomputes values. Returns 0 if it fails */
static ITG mg_coarse(ITG reuse)
{
  sparse_index_base_t base;
  ITG rows,cols,*rs=NULL,*re=NULL,*ci=NULL,i,nnz,maxfct=1,mnum=1,
      phase,mtype=-2,msglvl=0,error=0,*perm=NULL,one=1;
  double *cv=NULL,t0;
  unsigned long long h=1469598103934665603ULL;
  sparse_status_t st=SPARSE_STATUS_SUCCESS;

  t0=omp_get_wtime();
  if(mg_Ph==NULL)
    st=mkl_sparse_d_create_csr(&mg_Ph,SPARSE_INDEX_BASE_ONE,mg_n,mg_nc,mg_Pp,
			       mg_Pp+1,mg_Pc,mg_Pv);
  if(st==SPARSE_STATUS_SUCCESS){
    if(reuse&&mg_C!=NULL){
      st=mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE,mg_Ph,mg_A,mg_Ad,&mg_C,
			 SPARSE_STAGE_FINALIZE_MULT);
    }else{
      if(mg_C!=NULL){mkl_sparse_destroy(mg_C);mg_C=NULL;}
      st=mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE,mg_Ph,mg_A,mg_Ad,&mg_C,
			 SPARSE_STAGE_NNZ_COUNT);
      if(st==SPARSE_STATUS_SUCCESS)
	st=mkl_sparse_sypr(SPARSE_OPERATION_TRANSPOSE,mg_Ph,mg_A,mg_Ad,&mg_C,
			   SPARSE_STAGE_FINALIZE_MULT);
    }
  }
  if(st==SPARSE_STATUS_SUCCESS)
    st=mkl_sparse_d_export_csr(mg_C,&base,&rows,&cols,&rs,&re,&ci,&cv);
  if(st!=SPARSE_STATUS_SUCCESS||rows!=mg_nc){
    printf(" *WARNING in p-multigrid: Galerkin product failed (%d)\n",(int)st);
    if(mg_C){mkl_sparse_destroy(mg_C);mg_C=NULL;}
    return 0;
  }
  nnz=0;
  for(i=0;i<rows;i++) nnz+=re[i]-rs[i];
  if(mg_cp==NULL||nnz!=mg_cnnz){
    if(mg_cp!=NULL){SFREE(mg_cp);SFREE(mg_ci);SFREE(mg_cv);}
    NNEW(mg_cp,ITG,rows+1);NNEW(mg_ci,ITG,nnz);NNEW(mg_cv,double,nnz);
  }
  {
    ITG off=(base==SPARSE_INDEX_BASE_ONE) ? 1 : 0;
    mg_cp[0]=1;
    for(i=0;i<rows;i++) mg_cp[i+1]=mg_cp[i]+(re[i]-rs[i]);
#pragma omp parallel for
    for(i=0;i<rows;i++){                      /* copy, sort each row */
      ITG j,k=mg_cp[i]-1,k0=k,a,bb;
      for(j=rs[i]-off;j<re[i]-off;j++,k++){
	mg_ci[k]=ci[j]-off+1;mg_cv[k]=cv[j];
      }
      for(a=k0+1;a<k;a++){
	ITG c=mg_ci[a];double v=mg_cv[a];
	for(bb=a-1;bb>=k0&&mg_ci[bb]>c;bb--){
	  mg_ci[bb+1]=mg_ci[bb];mg_cv[bb+1]=mg_cv[bb];
	}
	mg_ci[bb+1]=c;mg_cv[bb+1]=v;
      }
    }
  }
  for(i=0;i<=rows;i++){h^=(unsigned long long)mg_cp[i];h*=1099511628211ULL;}
  for(i=0;i<nnz;i++){h^=(unsigned long long)mg_ci[i];h*=1099511628211ULL;}
  mg_tgal+=omp_get_wtime()-t0;

  t0=omp_get_wtime();
  if(mg_csingle){                       /* PARDISO single precision: float
					   values and right-hand sides */
    if(mg_cvf!=NULL) SFREE(mg_cvf);
    if(mg_cxf!=NULL) SFREE(mg_cxf);
    NNEW(mg_cvf,float,nnz);NNEW(mg_cxf,float,2*mg_nc);
#pragma omp parallel for
    for(i=0;i<nnz;i++) mg_cvf[i]=(float)mg_cv[i];
  }
  if(!mg_calive||nnz!=mg_cnnz||h!=mg_chash){
    if(mg_calive){
      phase=-1;
      FORTRAN(pardiso,(mg_pt,&maxfct,&mnum,&mtype,&phase,&mg_nc,mg_cv,mg_cp,
		       mg_ci,perm,&one,mg_iparm,&msglvl,NULL,NULL,&error));
    }
    for(i=0;i<64;i++){mg_pt[i]=0;mg_iparm[i]=0;}
    mg_iparm[0]=1;mg_iparm[1]=NK_ORDERING;
    mg_iparm[9]=8;mg_iparm[10]=0;mg_iparm[12]=0;mg_iparm[20]=1;
    mg_iparm[17]=-1;
    mg_iparm[27]=mg_csingle;               /* 1: single precision */
    mg_iparm[23]=pardiso_twolevel();
    {                                  /* saved ordering (PERM_DIR) */
      ITG *cperm=NULL,save=0;
      if(perm_dir()!=NULL){
	cperm=perm_load("coarse",h,mg_nc,nnz);
	if(cperm!=NULL) mg_iparm[4]=1;
	else{NNEW(cperm,ITG,mg_nc);mg_iparm[4]=2;save=1;}
      }
      phase=11;error=0;
      FORTRAN(pardiso,(mg_pt,&maxfct,&mnum,&mtype,&phase,&mg_nc,
		       mg_csingle ? (double*)mg_cvf : mg_cv,mg_cp,
		       mg_ci,cperm,&one,mg_iparm,&msglvl,NULL,NULL,&error));
      if(error==0&&save) perm_save("coarse",h,mg_nc,nnz,cperm);
      mg_iparm[4]=0;
      if(cperm!=NULL) SFREE(cperm);
    }
    if(error!=0){
      printf(" *WARNING in p-multigrid: coarse analysis error=%d\n",
	     (int)error);
      return 0;
    }
    mg_calive=1;mg_cnnz=nnz;mg_chash=h;mg_ncsym++;
  }
  phase=22;error=0;
  FORTRAN(pardiso,(mg_pt,&maxfct,&mnum,&mtype,&phase,&mg_nc,
		   mg_csingle ? (double*)mg_cvf : mg_cv,mg_cp,
		   mg_ci,perm,&one,mg_iparm,&msglvl,NULL,NULL,&error));
  mg_tcfac+=omp_get_wtime()-t0;
  mg_ncfac++;
  if(mg_verbose||mg_nc>500000)
    printf(" p-multigrid: coarse factorization %.2f s (%d equations,"
	   " nnz(L) %.3g)\n",omp_get_wtime()-t0,(int)mg_nc,
	   (double)mg_iparm[17]);
  if(error!=0){
    printf(" *WARNING in p-multigrid: coarse factorization error=%d\n",
	   (int)error);
    return 0;
  }
  return 1;
}

/* factor the stored fine matrix in the NK handle (CG failed) */
static void mg_to_nk(void)
{
  ITG *p=NULL,*c=NULL,nnz=mg_ap[mg_n]-1;
  double *a=NULL;
  NNEW(p,ITG,mg_n+1);NNEW(c,ITG,nnz);NNEW(a,double,nnz);
  memcpy(p,mg_ap,sizeof(ITG)*(mg_n+1));
  memcpy(c,mg_ai,sizeof(ITG)*nnz);
  memcpy(a,mg_av,sizeof(double)*nnz);
  nk_adopt(&mg_n,-2,p,c,a);
  nk_numeric(&mg_n);
  mg_fallback=1;mg_nfall++;
}

/* new matrix: returns 1 if the multigrid holds it, 0 if the caller
   must use the NK path instead */
static ITG mg_store(double *ad, double *au, double *adb, double *aub,
                double *sigma,ITG *icol, ITG *irow,
		ITG *neq, ITG *nzs, ITG *symmetryflag, ITG *inputformat,
		ITG *jq, ITG *nzs3)
{
  ITG mtype,*sp=pointers,*si=icolpardiso,i,reuse=0;
  double *sa=aupardiso,t0,t1;

  mg_options();
  if(mg_off||*symmetryflag!=0||*sigma!=0.){
    if(mg_alive) mg_release();
    return 0;
  }
  t0=omp_get_wtime();
  if(factor_alive) pardiso_cleanup(&alive_neq,&alive_symmetryflag);
  pardiso_threads();

  mtype=pardiso_csr(ad,au,adb,aub,sigma,icol,irow,neq,nzs,
		    symmetryflag,inputformat,jq,nzs3);
  mg_tcsr+=omp_get_wtime()-t0;
  {
    ITG nnz=pointers[*neq]-1;
    /* same pattern as the stored matrix: copy the values only, so the
       MKL handles and the Galerkin structure stay valid */
    if(mg_ap!=NULL&&mg_A!=NULL&&mg_n==*neq&&mg_ap[mg_n]-1==nnz&&
       memcmp(mg_ap,pointers,sizeof(ITG)*(*neq+1))==0&&
       memcmp(mg_ai,icolpardiso,sizeof(ITG)*nnz)==0){
      memcpy(mg_av,aupardiso,sizeof(double)*nnz);
      SFREE(pointers);SFREE(icolpardiso);SFREE(aupardiso);
      /* MKL keeps internal copies: tell it the values changed */
      if(mkl_sparse_d_update_values(mg_A,nnz,NULL,NULL,mg_av)!=
	 SPARSE_STATUS_SUCCESS){
	mkl_sparse_destroy(mg_A);mg_A=NULL;
	if(mg_C!=NULL){mkl_sparse_destroy(mg_C);mg_C=NULL;}
      }else{
	reuse=1;
      }
    }else{
      if(mg_C!=NULL){mkl_sparse_destroy(mg_C);mg_C=NULL;}
      if(mg_A!=NULL){mkl_sparse_destroy(mg_A);mg_A=NULL;}
      if(mg_ap!=NULL){SFREE(mg_ap);SFREE(mg_ai);SFREE(mg_av);}
      mg_ap=pointers;mg_ai=icolpardiso;mg_av=aupardiso;
    }
  }
  pointers=sp;icolpardiso=si;aupardiso=sa;
  mg_n=*neq;mg_alive=1;mg_fallback=0;mg_indef=mg_forcegm;mg_nmat++;
  (void)mtype;

  t1=omp_get_wtime();
  if(mg_P_n!=mg_n||mg_P_hash!=mg_mesh_hash){
    if(mg_Ph!=NULL){mkl_sparse_destroy(mg_Ph);mg_Ph=NULL;}
    if(mg_C!=NULL){mkl_sparse_destroy(mg_C);mg_C=NULL;}
    reuse=0;
    if(!mg_build_P(mg_n)){
      printf(" *WARNING: --MG needs C3D10 elements and the mesh from"
	     " nonlingeo; continuing with --NK\n\n");
      mg_off=1;mg_release();
      return 0;
    }
  }
  if(mg_wn!=mg_n+2*mg_nc){
    for(i=0;i<7;i++){if(mg_w[i]!=NULL) SFREE(mg_w[i]);}
    for(i=0;i<7;i++) NNEW(mg_w[i],double,mg_n+2*mg_nc);
    mg_wn=mg_n+2*mg_nc;
  }
  mg_tP+=omp_get_wtime()-t1;t1=omp_get_wtime();

  if(mg_A==NULL)
    mkl_sparse_d_create_csr(&mg_A,SPARSE_INDEX_BASE_ONE,mg_n,mg_n,mg_ap,
			    mg_ap+1,mg_ai,mg_av);
  if(!reuse||mg_mvn!=mg_n) mg_mv_setup();
  mg_Ad.type=SPARSE_MATRIX_TYPE_SYMMETRIC;
  mg_Ad.mode=SPARSE_FILL_MODE_UPPER;
  mg_Ad.diag=SPARSE_DIAG_NON_UNIT;
  mg_tmkl+=omp_get_wtime()-t1;t1=omp_get_wtime();

  mg_blocks();
  mg_tblk+=omp_get_wtime()-t1;t1=omp_get_wtime();
  /* lambda_max of the block-Jacobi scaled K_T changes little between
     Newton steps: estimated for every 4th matrix (and new patterns) */
  if(!reuse||mg_lmax<=0.||++mg_sinceeig>=4){
    mg_lmax=mg_eigmax();mg_sinceeig=0;
  }
  mg_teig+=omp_get_wtime()-t1;
  if(mg_directfirst&&mg_nmat==1){
    /* first matrix of the run: full load on the undeformed structure,
       the hardest right-hand side for the coarse level (thin skins) */
    printf(" p-multigrid: first matrix factored with PARDISO"
	   " (CALTOP_MG_DIRECT_FIRST)\n");
    mg_to_nk();mg_nfall--;
  }else{
    /* coarse factor reuse: keep the previous coarse matrix and factor (a
       slightly older K_T) while CG iterations stay near those right after
       the last refactorization; refactor when they double (and exceed
       20), on a new pattern, and for the stored converged K_T */
    /* refactor when the last solve took more than twice (and more than
       20) the iterations seen with a fresh coarse factor (age 0 or 1:
       this covers both kinds of Newton right-hand side, load increment
       and small residual), or at CALTOP_MG_COARSE_MAXAGE matrices */
    ITG keep=(mg_cstale>=0&&reuse&&mg_calive&&!mg_storing&&mg_cfresh>0&&
	      mg_clast>=0&&!(mg_clast>20&&mg_clast>2*mg_cfresh)&&
	      !(mg_cmaxage>0&&mg_cage+1>=mg_cmaxage));
    if(keep){
      if(mg_verbose) printf(" p-multigrid: coarse factor reused\n");
      mg_cstale++;mg_cage++;
    }else{
      if(!mg_coarse(reuse)){
	printf(" p-multigrid: no coarse level, factoring K_T instead\n");
	mg_to_nk();
      }
      mg_cfresh=-1;mg_cage=0;
      if(mg_cstale>=0) mg_cstale=0;
    }
  }
  mg_tsetup+=omp_get_wtime()-t0;
  if(mg_verbose)
    printf(" p-multigrid: setup %.2f s, lambda_max(B^-1 K) = %.3f\n",
	   omp_get_wtime()-t0,mg_lmax);
  return 1;
}

#define MG_RATE_W 20
static double mg_maxabs(const double *x,ITG n)
{
  double m=0.;
  ITG i;
#pragma omp parallel for reduction(max:m)
  for(i=0;i<n;i++){double a=fabs(x[i]);if(a>m) m=a;}
  return m;
}

/* PCG with the V-cycle. Converged when |r|_2 <= tol |b|_2 or, for Newton
   steps, when max|r_i| <= atol (atol > 0: below CalculiX's residual-force
   criterion). Returns the iterations on convergence, -1 on a breakdown
   (K_T not positive definite), -2 if it stalls or reaches the iteration
   limit; then x is the best iterate. *res: true relative residual */
static ITG mg_pcg(const double *b,double *x,double tol,double atol,
		  double *res)
{
  ITG n=mg_n,i,it,itbest=0,itbz=0,ret=-2;
  double *r=mg_w[0],*z=mg_w[1],*p=mg_w[2],*q=mg_w[3],*xb=NULL,*crhs=NULL,
         bn,rn,rz,rz1,pq,al,be,best,bestz,hist[MG_RATE_W],target,rz0;

  bn=sqrt(mg_dot(b,b,n));
  for(i=0;i<n;i++) x[i]=0.;
  *res=0.;
  if(bn==0.) return 0;
  if(atol>0.&&mg_maxabs(b,n)<=atol) return 0;  /* nothing to do */
  /* relative 2-norm target, used only to project the remaining
     iterations (the max-norm target, roughly, for Newton steps) */
  target=tol;
  if(atol>0.){double t=atol/mg_maxabs(b,n);if(t>target) target=t;}
  NNEW(crhs,double,2*mg_nc);NNEW(xb,double,n);
  best=bn;
  memcpy(r,b,sizeof(double)*n);
  mg_vcycle(r,z,crhs);
  memcpy(p,z,sizeof(double)*n);
  rz=mg_dot(r,z,n);
  bestz=rz;rz0=rz;
  for(it=1;it<=mg_maxit;it++){
    double ri=0.;
    mg_matvec(p,q);
    pq=mg_dot(p,q,n);
    if(!(pq>0.)||!(rz>0.)){
      if(mg_verbose) printf("   CG breakdown (p'Kp=%e, r'z=%e)\n",pq,rz);
      ret=-1;break;
    }
    al=rz/pq;
#pragma omp parallel for
    for(i=0;i<n;i++){x[i]+=al*p[i];r[i]-=al*q[i];}
    rn=sqrt(mg_dot(r,r,n));
    if(atol>0.) ri=mg_maxabs(r,n);
    if(mg_verbose) printf("   CG %4d  |r|/|b| = %.3e  max|r| = %.3e\n",(int)it,
			  rn/bn,ri);
    if(rn!=rn) break;
    if(rn<best){best=rn;itbest=it;memcpy(xb,x,sizeof(double)*n);}
    if(rn<=tol*bn||(atol>0.&&ri<=atol)){
      mg_matvec(x,q);                         /* true residual */
#pragma omp parallel for
      for(i=0;i<n;i++) r[i]=b[i]-q[i];
      rn=sqrt(mg_dot(r,r,n));
      if(rn<=10.*tol*bn||(atol>0.&&mg_maxabs(r,n)<=2.*atol)){
	ret=it;*res=rn/bn;break;
      }
    }
    mg_vcycle(r,z,crhs);
    rz1=mg_dot(r,z,n);
    if(rz1<0.999*bestz){bestz=rz1;itbz=it;}
    /* give up early: at the rate of the last MG_RATE_W iterations the
       preconditioned residual sqrt(r'M^-1 r) (smooth, unlike |r|_2 on
       thin skins) does not reach the target within mg_maxit iterations */
    hist[it%MG_RATE_W]=sqrt(bestz/rz0);
    if(it>=2*MG_RATE_W&&hist[it%MG_RATE_W]>target){
      double cur=hist[it%MG_RATE_W],old=hist[(it+1)%MG_RATE_W],rate;
      if(cur<old){
	rate=log(cur/old)/(double)(MG_RATE_W-1);
	if((double)it+log(target/cur)/rate>(double)mg_maxit){
	  if(mg_verbose) printf("   CG gives up: rate %.4f/it\n",exp(rate));
	  break;
	}
      }else if(it-itbz>2*MG_RATE_W){
	if(mg_verbose) printf("   CG gives up: no progress\n");
	break;
      }
    }
    /* stalled: neither |r|_2 nor the preconditioned residual r'M^-1 r
       (smoother, tracks the CG error) reached a new low for
       max(100,last improvement) iterations */
    {
      ITG last=(itbest>itbz)?itbest:itbz;
      if(it-last>((last>100)?last:100)) break;
    }
    be=rz1/rz;rz=rz1;
#pragma omp parallel for
    for(i=0;i<n;i++) p[i]=z[i]+be*p[i];
  }
  if(ret==-2){                     /* best iterate (or the last one if
				      |r|_2 never dropped), true residual */
    if(itbest>0) memcpy(x,xb,sizeof(double)*n);
    mg_matvec(x,q);
#pragma omp parallel for
    for(i=0;i<n;i++) r[i]=b[i]-q[i];
    *res=sqrt(mg_dot(r,r,n))/bn;
    if(*res!=*res) *res=1.e300;
  }
  SFREE(crhs);SFREE(xb);
  return ret;
}

/* right-preconditioned GMRES(20) with the V-cycle, for K_T that are not
   positive definite. Returns the iterations, or -2 (x: last iterate,
   *res: its true relative residual) */
#define MG_GMRES_M 20
static ITG mg_gmres(const double *b,double *x,double tol,double *res)
{
  ITG n=mg_n,m=MG_GMRES_M,it=0,i,j,k,jl,ret=-2,giveup=0;
  double *V=NULL,*w=mg_w[0],*z=mg_w[1],*u=mg_w[2],*crhs=NULL,
         H[(MG_GMRES_M+1)*MG_GMRES_M],cs[MG_GMRES_M],sn[MG_GMRES_M],
         g[MG_GMRES_M+1],y[MG_GMRES_M],bn,beta,tmp;

  bn=sqrt(mg_dot(b,b,n));
  for(i=0;i<n;i++) x[i]=0.;
  *res=0.;
  if(bn==0.) return 0;
  NNEW(V,double,(size_t)(m+1)*n);NNEW(crhs,double,2*mg_nc);
  memcpy(V,b,sizeof(double)*n);
  beta=bn;
  while(1){
#pragma omp parallel for
    for(i=0;i<n;i++) V[i]/=beta;
    for(i=0;i<=m;i++) g[i]=0.;
    g[0]=beta;jl=-1;
    for(j=0;j<m&&it<mg_maxit;j++){
      double *vj=&V[(size_t)j*n],*vj1=&V[(size_t)(j+1)*n];
      mg_vcycle(vj,z,crhs);
      mg_matvec(z,w);
      for(k=0;k<=j;k++){                      /* modified Gram-Schmidt */
	double *vk=&V[(size_t)k*n],hk=mg_dot(w,vk,n);
	H[k+j*(m+1)]=hk;
#pragma omp parallel for
	for(i=0;i<n;i++) w[i]-=hk*vk[i];
      }
      tmp=sqrt(mg_dot(w,w,n));
      H[j+1+j*(m+1)]=tmp;
      if(tmp>0.){
#pragma omp parallel for
	for(i=0;i<n;i++) vj1[i]=w[i]/tmp;
      }
      for(k=0;k<j;k++){
	double h0=H[k+j*(m+1)],h1=H[k+1+j*(m+1)];
	H[k+j*(m+1)]=cs[k]*h0+sn[k]*h1;
	H[k+1+j*(m+1)]=-sn[k]*h0+cs[k]*h1;
      }
      {
	double h0=H[j+j*(m+1)],h1=H[j+1+j*(m+1)],rr=sqrt(h0*h0+h1*h1);
	if(rr==0.){cs[j]=1.;sn[j]=0.;}else{cs[j]=h0/rr;sn[j]=h1/rr;}
	H[j+j*(m+1)]=rr;H[j+1+j*(m+1)]=0.;
	g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];
      }
      it++;jl=j;
      if(mg_verbose) printf("   GMRES %4d  |r|/|b| = %.3e\n",(int)it,
			    fabs(g[j+1])/bn);
      if(fabs(g[j+1])<=tol*bn||tmp==0.) break;
      /* give up early: at the rate reached so far the tolerance needs
         more than mg_maxit iterations (as for --NK) */
      if(it>=NK_GIVEUP_MIN){
	double red=fabs(g[j+1])/bn;
	if(red>=1.||(double)it*log(tol)/log(red)>(double)mg_maxit){
	  giveup=1;break;
	}
      }
    }
    for(k=jl;k>=0;k--){                       /* H y = g */
      double sk=g[k];
      for(i=k+1;i<=jl;i++) sk-=H[k+i*(m+1)]*y[i];
      y[k]=(H[k+k*(m+1)]!=0.) ? sk/H[k+k*(m+1)] : 0.;
    }
    for(i=0;i<n;i++) u[i]=0.;
    for(k=0;k<=jl;k++){
      double *vk=&V[(size_t)k*n],yk=y[k];
#pragma omp parallel for
      for(i=0;i<n;i++) u[i]+=yk*vk[i];
    }
    mg_vcycle(u,z,crhs);                      /* x += M^-1 V y */
#pragma omp parallel for
    for(i=0;i<n;i++) x[i]+=z[i];
    mg_matvec(x,w);                           /* true residual */
#pragma omp parallel for
    for(i=0;i<n;i++) V[i]=b[i]-w[i];
    beta=sqrt(mg_dot(V,V,n));
    *res=beta/bn;
    if(beta!=beta){*res=1.e300;break;}
    if(beta<=10.*tol*bn&&fabs(g[jl+1])<=tol*bn){ret=it;break;}
    if(it>=mg_maxit||giveup) break;
  }
  SFREE(V);SFREE(crhs);
  return ret;
}

static void mg_solve(double *b,ITG *neq,ITG *nrhs,double tol,double atol,
		     ITG adjoint)
{
  ITG r,it,i;
  double *x=NULL,t0,dt,res,acc;

  if(*neq!=mg_n){
    printf(" *ERROR in pardiso_solve: no matching MG matrix stored\n");
    FORTRAN(stop,());
  }
  acc=(tol>mg_accept) ? tol : mg_accept;
  NNEW(x,double,(size_t)*nrhs**neq);
  for(r=0;r<*nrhs;r++){
    double *br=&b[(size_t)r**neq],*xr=&x[(size_t)r**neq];
    if(!mg_fallback){
      const char *kry="CG";
      double mv0=mg_tmv;
      long nmv0=mg_nmv;
      double *x0=NULL,*r0=NULL,scale=1.,tolr=tol;
      t0=omp_get_wtime();
      it=-1;res=1.e300;
      /* adjoint: start from alpha*u, alpha = u'b / u'Ku (the best multiple
         of u in the energy norm, never worse than a zero start); CG then
         solves for the correction with the right-hand side b - alpha K u */
      if(adjoint&&!mg_indef&&mg_guess!=NULL&&mg_guessn==*neq){
	double *ku=mg_w[4],ub,uku,bn,rn;
	mg_matvec(mg_guess,ku);
	ub=mg_dot(mg_guess,br,*neq);uku=mg_dot(mg_guess,ku,*neq);
	bn=sqrt(mg_dot(br,br,*neq));
	if(uku>0.&&bn>0.){
	  double al=ub/uku;
	  NNEW(r0,double,*neq);
#pragma omp parallel for
	  for(i=0;i<*neq;i++) r0[i]=br[i]-al*ku[i];
	  rn=sqrt(mg_dot(r0,r0,*neq));
	  printf(" p-multigrid: adjoint started from %.3g x converged"
		 " displacement, |b - K x0|/|b| = %.2e\n",al,rn/bn);
	  if(rn<bn){
	    NNEW(x0,double,*neq);
	    for(i=0;i<*neq;i++) x0[i]=al*mg_guess[i];
	    scale=rn/bn;tolr=tol/scale;
	  }else{SFREE(r0);r0=NULL;}
	}
      }
      if(!mg_indef&&x0!=NULL&&tolr>=1.){       /* guess already good enough */
	memcpy(xr,x0,sizeof(double)**neq);
	it=0;res=scale;
      }else if(!mg_indef){
	it=mg_pcg(r0 ? r0 : br,xr,tolr,atol,&res);
	if(x0!=NULL){
#pragma omp parallel for
	  for(i=0;i<*neq;i++) xr[i]+=x0[i];
	  res*=scale;
	}
	if(it==-1){
	  mg_indef=1;
	  printf(" p-multigrid: K_T is not positive definite, CG -> GMRES"
		 " (%.2f s lost)\n",omp_get_wtime()-t0);
	}
      }
      if(mg_indef){kry="GMRES";it=mg_gmres(br,xr,tol,&res);}
      if(x0!=NULL) SFREE(x0);
      if(r0!=NULL) SFREE(r0);
      dt=omp_get_wtime()-t0;mg_tcg+=dt;
      mg_tmvcg+=mg_tmv-mv0;mg_nmvcg+=mg_nmv-nmv0;
      if(it>=0||res<=acc){
	ITG used=(it>=0)?it:mg_maxit;
	mg_nsolve++;mg_nit+=used;mg_failrow=0;
	mg_clast=used;
	if(mg_cage<=1&&used>mg_cfresh) mg_cfresh=(used>1)?used:1;
	if(it>=0)
	  printf(" p-multigrid: %s %d it, |r|/|b| = %.1e, %.2f s\n",kry,
		 (int)it,res,dt);
	else
	  printf(" p-multigrid: %s stalled at %.1e (asked %.1e), accepted,"
		 " %.2f s\n",kry,res,tol,dt);
	continue;
      }
      mg_clast=mg_maxit;
      printf(" p-multigrid: %s did not converge to %.1e (reached %.1e,"
	     " %.2f s) -> PARDISO factorization of K_T\n",kry,tol,res,dt);
      mg_to_nk();
      if(++mg_failrow>=2){
	printf(" *WARNING: p-multigrid failed twice in a row; continuing"
	       " with --NK\n\n");
	mg_off=1;
      }
    }
    {
      ITG one=1,sym=0;
      memcpy(xr,br,sizeof(double)**neq);
      nk_solve(xr,neq,&sym,&one,tol);
    }
  }
  for(i=0;i<*nrhs**neq;i++) b[i]=x[i];
  SFREE(x);
}

static void mg_release(void)
{
  ITG i,maxfct=1,mnum=1,phase=-1,mtype=-2,msglvl=0,error=0,*perm=NULL,one=1;
  if(!mg_alive&&!mg_calive) return;
  if(mg_nmat>0)
    printf(" p-multigrid: matrices=%ld CG solves=%ld (avg %.1f it)"
	   " fallbacks=%ld coarse analyses=%ld\n"
	   "              time: setup %.1f s (Galerkin %.1f s, coarse factor"
	   " %.1f s), CG %.1f s\n",
	   mg_nmat,mg_nsolve,mg_nsolve>0 ? (double)mg_nit/(double)mg_nsolve:0.,
	   mg_nfall,mg_ncsym,mg_tsetup,mg_tgal,mg_tcfac,mg_tcg);
  if(mg_nmat>0)
    printf("              coarse: %d equations, %ld factorizations for %ld"
	   " matrices\n",(int)mg_nc,mg_ncfac,mg_nmat);
  if(mg_nmat>0&&mg_tcg>0.)
    printf("              CG split: fine products %.1f s (%ld, %.3f s each%s),"
	   " coarse solves %.1f s (%ld, %.3f s each), rest %.1f s\n",
	   mg_tmvcg,mg_nmvcg,mg_nmvcg>0 ? mg_tmvcg/(double)mg_nmvcg : 0.,
	   mg_mkl_mv ? ", MKL" : "",mg_tcs,mg_ncs,
	   mg_ncs>0 ? mg_tcs/(double)mg_ncs : 0.,mg_tcg-mg_tmvcg-mg_tcs);
  if(mg_nmat>0)
    printf("              setup detail: CSR %.1f s, P %.1f s, MKL handle %.1f s,"
	   " blocks %.1f s, lambda_max %.1f s\n",mg_tcsr,mg_tP,mg_tmkl,mg_tblk,
	   mg_teig);
  if(mg_calive)
    FORTRAN(pardiso,(mg_pt,&maxfct,&mnum,&mtype,&phase,&mg_nc,mg_cv,mg_cp,
		     mg_ci,perm,&one,mg_iparm,&msglvl,NULL,NULL,&error));
  if(mg_C!=NULL){mkl_sparse_destroy(mg_C);mg_C=NULL;}
  if(mg_Ph!=NULL){mkl_sparse_destroy(mg_Ph);mg_Ph=NULL;}
  if(mg_A!=NULL){mkl_sparse_destroy(mg_A);mg_A=NULL;}
  if(mg_ap!=NULL){SFREE(mg_ap);SFREE(mg_ai);SFREE(mg_av);}
  if(mg_cp!=NULL){SFREE(mg_cp);SFREE(mg_ci);SFREE(mg_cv);}
  mg_lmax=0.;mg_sinceeig=0;
  for(i=0;i<7;i++){if(mg_w[i]!=NULL){SFREE(mg_w[i]);}}
  mg_wn=0;
  mg_alive=0;mg_calive=0;mg_cnnz=-1;mg_chash=0ULL;mg_fallback=0;mg_indef=0;
  mg_nmat=mg_nsolve=mg_nit=mg_nfall=mg_ncsym=0;
  mg_tsetup=mg_tgal=mg_tcfac=mg_tcg=0.;
  mg_tcsr=mg_tP=mg_tmkl=mg_tblk=mg_teig=0.;
  mg_tmv=mg_tcs=mg_tmvcg=0.;mg_nmv=mg_ncs=mg_nmvcg=0;
  mg_mv_release();
  mg_ncfac=0;mg_cfresh=-1;mg_cage=0;mg_clast=-1;if(mg_cstale>0) mg_cstale=0;
  if(mg_cvf!=NULL){SFREE(mg_cvf);mg_cvf=NULL;}
  if(mg_cxf!=NULL){SFREE(mg_cxf);mg_cxf=NULL;}
}

#endif
