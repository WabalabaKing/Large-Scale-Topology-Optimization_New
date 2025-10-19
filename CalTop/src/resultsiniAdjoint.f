c======================================================================
c  Adjoint results initializer (linear static, mechanical only, no MPC)
c    - Scatter adjoint solution vector (equation space) -> nodal field ladj
c    - Enforce adjoint BCs: lambda = 0 on primal SPC DOFs
c======================================================================
      subroutine resultsini_adjoint_linstatic_nompc(
     &     nk, ladj, nactdof, b_adj,
     &     nodeboun, ndirboun, typeboun, nboun, mi)

      implicit none

c----- arguments
      integer           nk, nboun
      integer           mi(*), nactdof(0:mi(2),*)
      integer           nodeboun(*), ndirboun(*)
      character*1       typeboun(*)         ! 'F' = force BC (skip)
      real*8            ladj(0:mi(2),*)     ! nodal adjoint field
      real*8            b_adj(*)            ! adjoint solution (size = neq)

c----- locals
      integer           i, j, node, ndir
      integer           mt

c----- convenience
      mt = mi(2) + 1

c===== 0) Initialize adjoint field to zero
      do i = 1, nk
        do j = 1, mi(2)
          ladj(j,i) = 0.d0
        enddo
      enddo

c===== 1) Scatter adjoint unknowns by assignment
      do i = 1, nk
        do j = 1, mi(2)
          if (nactdof(j,i) .gt. 0) then
            ladj(j,i) = b_adj( nactdof(j,i) )
          endif
        enddo
      enddo

c===== 2) Enforce adjoint BCs for primal SPCs (Dirichlet)
      do i = 1, nboun
        if (typeboun(i) .eq. 'F') cycle
        if (ndirboun(i) .gt. mi(2)) cycle

        ndir = ndirboun(i)
        node = nodeboun(i)
        ladj(ndir, node) = 0.d0
      enddo

      return
      end
