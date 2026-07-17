!
!     Linear C3D4 compliance and compliance sensitivity
!
!     This routine retains the original e_c3d_se argument list so that
!     it can be called directly from the existing CalculiX source.
!
!     Supported:
!       - C3D4 only
!       - small-displacement linear elasticity
!       - isotropic material
!       - SIMP stiffness interpolation
!       - compliance and sensitivity
!       - element volume and centroid
!
      subroutine linearCompliance(
     &  co,kon,lakonl,p1,p2,omx,bodyfx,nbody,s,sm,
     &  ff,nelem,nmethod,elcon,nelcon,rhcon,nrhcon,
     &  alcon,nalcon,alzero,ielmat,ielorien,norien,
     &  orab,ntmat_,t0,t1,ithermal,vold,iperturb,
     &  nelemload,sideload,xload,nload,idist,sti,stx,
     &  iexpl,plicon,nplicon,plkcon,nplkcon,xstiff,
     &  npmat_,dtime,matname,mi,ncmat_,mass,stiffness,
     &  buckling,rhsi,intscheme,ttime,time,istep,iinc,
     &  coriolis,xloadold,reltime,ipompc,nodempc,
     &  coefmpc,nmpc,ikmpc,ilmpc,veold,springarea,
     &  nstate_,xstateini,xstate,ne0,ipkon,thicke,
     &  integerglob,doubleglob,tieset,istartset,
     &  iendset,ialset,ntie,nasym,pslavsurf,
     &  pmastsurf,mortar,clearini,ielprop,prop,
     &  distmin,ndesi,nodedesi,dfl,icoordinate,
     &  dxstiff,ne,xdesi,istartelem,ialelem,v,sigma,
     &  ieigenfrequency,rhoi,penal,sensi,ecompli,
     &  elvol,xcg,ycg,zcg,fn0_out,lambda,nactdof,neq)

      implicit none

!
!---------------------------------------------------------------------
!     Original CalculiX dummy arguments
!---------------------------------------------------------------------
!
      character*8 lakonl
      character*20 sideload(*)
      character*80 matname(*)
      character*81 tieset(3,*)

      integer nbody
      integer nelem
      integer nmethod
      integer ntmat_
      integer ithermal
      integer nload
      integer idist
      integer iexpl
      integer npmat_
      integer ncmat_
      integer mass
      integer stiffness
      integer buckling
      integer rhsi
      integer intscheme
      integer istep
      integer iinc
      integer coriolis
      integer nmpc
      integer nstate_
      integer ne0
      integer ntie
      integer nasym
      integer mortar
      integer ndesi
      integer icoordinate
      integer ne
      integer ieigenfrequency
      integer norien

      integer kon(*)
      integer mi(*)
      integer nelcon(2,*)
      integer nrhcon(*)
      integer nalcon(2,*)
      integer ielmat(mi(3),*)
      integer ielorien(mi(3),*)
      integer iperturb(*)
      integer nelemload(2,*)
      integer nplicon(0:ntmat_,*)
      integer nplkcon(0:ntmat_,*)
      integer ipompc(*)
      integer nodempc(3,*)
      integer ikmpc(*)
      integer ilmpc(*)
      integer ipkon(*)
      integer integerglob(*)
      integer istartset(*)
      integer iendset(*)
      integer ialset(*)
      integer ielprop(*)
      integer nodedesi(*)
      integer istartelem(*)
      integer ialelem(*)
      integer nactdof(0:mi(2),*)
      integer neq(*)

      real*8 co(3,*)
      real*8 p1(3)
      real*8 p2(3)
      real*8 bodyfx(3)

      real*8 s(60,60)
      real*8 sm(60,60)
      real*8 ff(60)

      real*8 elcon(0:ncmat_,ntmat_,*)
      real*8 rhcon(0:1,ntmat_,*)
      real*8 alcon(0:6,ntmat_,*)
      real*8 alzero(*)
      real*8 orab(7,*)

      real*8 t0(*)
      real*8 t1(*)
      real*8 vold(0:mi(2),*)
      real*8 xload(2,*)
      real*8 sti(6,mi(1),*)
      real*8 stx(6,mi(1),*)

      real*8 plicon(0:2*npmat_,ntmat_,*)
      real*8 plkcon(0:2*npmat_,ntmat_,*)
      real*8 xstiff(27,mi(1),*)

      real*8 dtime
      real*8 ttime
      real*8 time
      real*8 xloadold(2,*)
      real*8 reltime

      real*8 coefmpc(*)
      real*8 veold(0:mi(2),*)
      real*8 springarea(2,*)

      real*8 xstateini(nstate_,mi(1),*)
      real*8 xstate(nstate_,mi(1),*)
      real*8 thicke(mi(3),*)

      real*8 doubleglob(*)
      real*8 pslavsurf(3,*)
      real*8 pmastsurf(6,*)
      real*8 clearini(3,9,*)
      real*8 prop(*)
      real*8 distmin

      real*8 dfl(20,60)
      real*8 dxstiff(27,mi(1),ne,*)
      real*8 xdesi(3,*)

      real*8 v(0:mi(2),*)
      real*8 sigma
      real*8 rhoi
      real*8 penal

      real*8 sensi
      real*8 ecompli
      real*8 elvol
      real*8 xcg
      real*8 ycg
      real*8 zcg

      real*8 fn0_out(0:mi(2),*)
      real*8 lambda(*)

      real*8 omx

!
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
!
      integer nope
      integer ndof
      integer indexe
      integer imat

      integer konl(4)

      integer i
      integer j
      integer inode
      integer idof
      integer jdof
      integer row
      integer col

      integer iflag

      real*8 xl(3,4)
      real*8 shp(4,4)
      real*8 uelem(12)

      real*8 xi
      real*8 et
      real*8 ze
      real*8 xsj
      real*8 weight

      real*8 young
      real*8 poisson
      real*8 lambdaL
      real*8 shear

      real*8 B(6,12)
      real*8 D(6,6)
      real*8 DB(6,12)
      real*8 ke(12,12)

      real*8 uku
      real*8 rhoSafe

      include "gauss.f"

!
!---------------------------------------------------------------------
!     Initialization
!---------------------------------------------------------------------
!
      nope=4
      ndof=12

      sensi=0.d0
      ecompli=0.d0
      elvol=0.d0

      xcg=0.d0
      ycg=0.d0
      zcg=0.d0

      uku=0.d0

!
!     This specialized routine accepts C3D4 elements only.
!
      if(lakonl(4:4).ne.'4') then
         write(*,'(A,I10,A,A8)')
     &      ' *ERROR in linearCompliance: element ',nelem,
     &      ' is not C3D4. Type = ',lakonl
         nmethod=0
         return
      endif

!
!     Locate the element connectivity.
!
      indexe=ipkon(nelem)

      if(indexe.lt.0) then
         write(*,'(A,I10)')
     &      ' *ERROR in linearCompliance: inactive element ',nelem
         nmethod=0
         return
      endif

!
!---------------------------------------------------------------------
!     Gather connectivity, coordinates and displacements
!---------------------------------------------------------------------
!
      do inode=1,nope

         konl(inode)=kon(indexe+inode)

         do j=1,3
            xl(j,inode)=co(j,konl(inode))
         enddo

         uelem(3*inode-2)=v(1,konl(inode))
         uelem(3*inode-1)=v(2,konl(inode))
         uelem(3*inode  )=v(3,konl(inode))

      enddo

!
!---------------------------------------------------------------------
!     Element centroid
!
!     For a linear tetrahedron, the centroid is the arithmetic average
!     of its four nodal coordinates.
!---------------------------------------------------------------------
!
      do inode=1,nope
         xcg=xcg+xl(1,inode)
         ycg=ycg+xl(2,inode)
         zcg=zcg+xl(3,inode)
      enddo

      xcg=xcg/4.d0
      ycg=ycg/4.d0
      zcg=zcg/4.d0

!
!---------------------------------------------------------------------
!     Material properties
!
!     For an isotropic *ELASTIC material in CalculiX:
!
!       elcon(1,1,imat) = Young's modulus
!       elcon(2,1,imat) = Poisson's ratio
!
!---------------------------------------------------------------------
!
      imat=ielmat(1,nelem)

      if(imat.le.0) then
         write(*,'(A,I10)')
     &      ' *ERROR in linearCompliance: no material for element ',
     &      nelem
         nmethod=0
         return
      endif

      young=elcon(1,1,imat)
      poisson=elcon(2,1,imat)

      if(young.le.0.d0) then
         write(*,'(A,I10,A,ES14.6)')
     &      ' *ERROR in linearCompliance: invalid E for element ',
     &      nelem,'; E = ',young
         nmethod=0
         return
      endif

      if((poisson.le.-1.d0).or.(poisson.ge.0.5d0)) then
         write(*,'(A,I10,A,ES14.6)')
     &      ' *ERROR in linearCompliance: invalid nu for element ',
     &      nelem,'; nu = ',poisson
         nmethod=0
         return
      endif

!
!     Lamé constants.
!
      shear=young/(2.d0*(1.d0+poisson))

      lambdaL=young*poisson/
     &       ((1.d0+poisson)*(1.d0-2.d0*poisson))

!
!---------------------------------------------------------------------
!     Evaluate the C3D4 shape-function derivatives
!---------------------------------------------------------------------
!
      xi=gauss3d4(1,1)
      et=gauss3d4(2,1)
      ze=gauss3d4(3,1)
      weight=weight3d4(1)

      iflag=3

      call shape4tet(xi,et,ze,xl,xsj,shp,iflag)

      if(xsj.le.1.d-20) then
         write(*,'(A,I10,A,ES14.6)')
     &      ' *ERROR in linearCompliance: nonpositive Jacobian for ',
     &      nelem,'; detJ = ',xsj
         nmethod=0
         return
      endif

!
!     Element volume:
!
!        V_e = weight * det(J)
!
      elvol=weight*xsj

!
!---------------------------------------------------------------------
!     Construct the strain-displacement matrix B
!
!     shape4tet returns:
!
!       shp(1,i) = dN_i/dx
!       shp(2,i) = dN_i/dy
!       shp(3,i) = dN_i/dz
!
!     The engineering shear strain ordering is:
!
!       epsilon =
!       [epsilon_xx, epsilon_yy, epsilon_zz,
!        gamma_xy, gamma_yz, gamma_xz]^T
!---------------------------------------------------------------------
!
      do i=1,6
         do j=1,ndof
            B(i,j)=0.d0
         enddo
      enddo

      do inode=1,nope

         col=3*inode-2

         B(1,col  )=shp(1,inode)
         B(2,col+1)=shp(2,inode)
         B(3,col+2)=shp(3,inode)

         B(4,col  )=shp(2,inode)
         B(4,col+1)=shp(1,inode)

         B(5,col+1)=shp(3,inode)
         B(5,col+2)=shp(2,inode)

         B(6,col  )=shp(3,inode)
         B(6,col+2)=shp(1,inode)

      enddo

!
!---------------------------------------------------------------------
!     Construct the isotropic elasticity matrix D
!---------------------------------------------------------------------
!
      do i=1,6
         do j=1,6
            D(i,j)=0.d0
         enddo
      enddo

      D(1,1)=lambdaL+2.d0*shear
      D(1,2)=lambdaL
      D(1,3)=lambdaL

      D(2,1)=lambdaL
      D(2,2)=lambdaL+2.d0*shear
      D(2,3)=lambdaL

      D(3,1)=lambdaL
      D(3,2)=lambdaL
      D(3,3)=lambdaL+2.d0*shear

      D(4,4)=shear
      D(5,5)=shear
      D(6,6)=shear

!
!---------------------------------------------------------------------
!     Compute DB = D B
!---------------------------------------------------------------------
!
      do i=1,6
         do j=1,ndof

            DB(i,j)=0.d0

            do row=1,6
               DB(i,j)=DB(i,j)+D(i,row)*B(row,j)
            enddo

         enddo
      enddo

!
!---------------------------------------------------------------------
!     Compute the unpenalized element stiffness:
!
!       K0_e = integral(B^T D B dV)
!            = B^T D B * weight * det(J)
!
!     For C3D4, B is constant over the element.
!---------------------------------------------------------------------
!
      do i=1,ndof
         do j=1,ndof

            ke(i,j)=0.d0

            do row=1,6
               ke(i,j)=ke(i,j)+B(row,i)*DB(row,j)
            enddo

            ke(i,j)=ke(i,j)*elvol

         enddo
      enddo

!
!---------------------------------------------------------------------
!     Return the element stiffness matrix in s
!
!     The matrix returned here is the unpenalized solid stiffness K0_e.
!     The SIMP scaling should be applied during global assembly if that
!     is how the existing CalTop assembly is implemented.
!---------------------------------------------------------------------
!
      do i=1,60
         do j=1,60
            s(i,j)=0.d0
         enddo
      enddo

      do i=1,ndof
         do j=1,ndof
            s(i,j)=ke(i,j)
         enddo
      enddo

!
!---------------------------------------------------------------------
!     Compute the unpenalized elemental strain energy:
!
!       uku = u_e^T K0_e u_e
!---------------------------------------------------------------------
!
      uku=0.d0

      do i=1,ndof
         do j=1,ndof
            uku=uku+uelem(i)*ke(i,j)*uelem(j)
         enddo
      enddo

!
!---------------------------------------------------------------------
!     SIMP compliance and sensitivity
!
!     K_e = rho_e^p K0_e
!
!     Element compliance:
!
!       C_e = rho_e^p u_e^T K0_e u_e
!
!     Total compliance derivative at equilibrium:
!
!       dC/drho_e =
!          -p rho_e^(p-1) u_e^T K0_e u_e
!---------------------------------------------------------------------
!
!     Prevent a numerical problem if rhoi is exactly zero and
!     penal is less than one. In normal SIMP calculations, rhoi should
!     already be bounded below by rho_min.
!
      rhoSafe=max(rhoi,1.d-12)

      ecompli=(rhoSafe**penal)*uku

      sensi=-penal*
     &       (rhoSafe**(penal-1.d0))*
     &       uku

      return
      end