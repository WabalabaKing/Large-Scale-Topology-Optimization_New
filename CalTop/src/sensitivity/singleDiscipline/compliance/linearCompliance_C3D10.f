!
!     Linear C3D10 compliance and compliance sensitivity
!
!     Supported:
!       - C3D10 elements only
!       - small-displacement linear elasticity
!       - isotropic material
!       - SIMP compliance and sensitivity
!       - element stiffness matrix
!       - integrated element volume
!       - integrated volume centroid
!
      subroutine linearCompliance_C3D10(
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
      integer mint3d
      integer indexe
      integer imat

      integer konl(10)

      integer i
      integer j
      integer inode
      integer row
      integer col
      integer kk
      integer iflag

      real*8 xl(3,10)
      real*8 shp(4,10)
      real*8 uelem(30)

      real*8 xi
      real*8 et
      real*8 ze
      real*8 xsj
      real*8 weight
      real*8 dvol

      real*8 young
      real*8 poisson
      real*8 lambdaL
      real*8 shear

      real*8 B(6,30)
      real*8 D(6,6)
      real*8 DB(6,30)
      real*8 ke(30,30)

      real*8 xgp
      real*8 ygp
      real*8 zgp

      real*8 xmoment
      real*8 ymoment
      real*8 zmoment

      real*8 uku
      real*8 rhoSafe

      include "gauss.f"

!
!---------------------------------------------------------------------
!     Initialization
!---------------------------------------------------------------------
!
      nope=10
      ndof=30

!
!     CalculiX uses four-point integration for the standard C3D10
!     integration scheme and 15-point integration when an alternate
!     integration scheme is requested.
!
      if(intscheme.eq.0) then
         mint3d=4
      else
         mint3d=15
      endif

      sensi=0.d0
      ecompli=0.d0
      elvol=0.d0

      xcg=0.d0
      ycg=0.d0
      zcg=0.d0

      xmoment=0.d0
      ymoment=0.d0
      zmoment=0.d0

      uku=0.d0

!
!     Initialize local arrays.
!
      do i=1,30
         uelem(i)=0.d0

         do j=1,30
            ke(i,j)=0.d0
         enddo
      enddo

      do i=1,60
         do j=1,60
            s(i,j)=0.d0
         enddo
      enddo

!
!---------------------------------------------------------------------
!     Verify the element type
!---------------------------------------------------------------------
!
      if(lakonl(4:5).ne.'10') then

         write(*,'(A,I10,A,A8)')
     &      ' *ERROR in linearCompliance_C3D10: element ',
     &      nelem,' is not C3D10. Type = ',lakonl

         nmethod=0
         return

      endif

!
!---------------------------------------------------------------------
!     Locate element connectivity
!---------------------------------------------------------------------
!
      indexe=ipkon(nelem)

      if(indexe.lt.0) then

         write(*,'(A,I10)')
     &      ' *ERROR in linearCompliance_C3D10: inactive element ',
     &      nelem

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
!     Material properties
!
!     For a standard isotropic CalculiX *ELASTIC material:
!
!       elcon(1,1,imat) = Young's modulus
!       elcon(2,1,imat) = Poisson's ratio
!---------------------------------------------------------------------
!
      imat=ielmat(1,nelem)

      if(imat.le.0) then

         write(*,'(A,I10)')
     &      ' *ERROR in linearCompliance_C3D10: no material for ',
     &      nelem

         nmethod=0
         return

      endif

      young=elcon(1,1,imat)
      poisson=elcon(2,1,imat)

      if(young.le.0.d0) then

         write(*,'(A,I10,A,ES14.6)')
     &      ' *ERROR in linearCompliance_C3D10: invalid E for ',
     &      nelem,'; E = ',young

         nmethod=0
         return

      endif

      if((poisson.le.-1.d0).or.(poisson.ge.0.5d0)) then

         write(*,'(A,I10,A,ES14.6)')
     &      ' *ERROR in linearCompliance_C3D10: invalid nu for ',
     &      nelem,'; nu = ',poisson

         nmethod=0
         return

      endif

!
!---------------------------------------------------------------------
!     Lamé constants
!---------------------------------------------------------------------
!
      shear=young/(2.d0*(1.d0+poisson))

      lambdaL=young*poisson/
     &       ((1.d0+poisson)*(1.d0-2.d0*poisson))

!
!---------------------------------------------------------------------
!     Construct isotropic elasticity matrix D
!
!     Stress/strain ordering:
!
!       xx, yy, zz, xy, yz, xz
!
!     Engineering shear strains are used.
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
!     Integrate stiffness, volume and centroid
!
!       K0_e = integral(B^T D B dV)
!
!       V_e  = integral(dV)
!
!       xcg  = integral(x dV)/V_e
!---------------------------------------------------------------------
!
      iflag=3

      do kk=1,mint3d

!
!        Select the tetrahedral integration point.
!
         if(intscheme.eq.0) then

            xi=gauss3d5(1,kk)
            et=gauss3d5(2,kk)
            ze=gauss3d5(3,kk)
            weight=weight3d5(kk)

         else

            xi=gauss3d6(1,kk)
            et=gauss3d6(2,kk)
            ze=gauss3d6(3,kk)
            weight=weight3d6(kk)

         endif

!
!        Evaluate C3D10 shape functions and global derivatives.
!
!        shp(1,i) = dN_i/dx
!        shp(2,i) = dN_i/dy
!        shp(3,i) = dN_i/dz
!        shp(4,i) = N_i
!
         call shape10tet(xi,et,ze,xl,xsj,shp,iflag)

!
!        Check the Jacobian determinant.
!
         if(xsj.le.1.d-20) then

            write(*,'(A,I10,A,I5,A,ES14.6)')
     &         ' *ERROR in linearCompliance_C3D10: element ',
     &         nelem,', Gauss point ',kk,
     &         ', nonpositive detJ = ',xsj

            nmethod=0
            return

         endif

!
!        Differential element volume.
!
         dvol=weight*xsj
         elvol=elvol+dvol

!
!        Physical coordinates of this integration point.
!
         xgp=0.d0
         ygp=0.d0
         zgp=0.d0

         do inode=1,nope

            xgp=xgp+shp(4,inode)*xl(1,inode)
            ygp=ygp+shp(4,inode)*xl(2,inode)
            zgp=zgp+shp(4,inode)*xl(3,inode)

         enddo

!
!        Accumulate first volume moments.
!
         xmoment=xmoment+xgp*dvol
         ymoment=ymoment+ygp*dvol
         zmoment=zmoment+zgp*dvol

!
!---------------------------------------------------------------------
!        Construct the strain-displacement matrix B
!---------------------------------------------------------------------
!
         do i=1,6
            do j=1,ndof
               B(i,j)=0.d0
            enddo
         enddo

         do inode=1,nope

            col=3*inode-2

!
!           Normal strains.
!
            B(1,col  )=shp(1,inode)
            B(2,col+1)=shp(2,inode)
            B(3,col+2)=shp(3,inode)

!
!           gamma_xy.
!
            B(4,col  )=shp(2,inode)
            B(4,col+1)=shp(1,inode)

!
!           gamma_yz.
!
            B(5,col+1)=shp(3,inode)
            B(5,col+2)=shp(2,inode)

!
!           gamma_xz.
!
            B(6,col  )=shp(3,inode)
            B(6,col+2)=shp(1,inode)

         enddo

!
!---------------------------------------------------------------------
!        Compute DB = D B
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
!        Accumulate element stiffness contribution
!
!           K0_e += B^T D B dV
!---------------------------------------------------------------------
!
         do i=1,ndof
            do j=1,ndof

               do row=1,6

                  ke(i,j)=ke(i,j)+
     &                 B(row,i)*DB(row,j)*dvol

               enddo

            enddo
         enddo

      enddo

!
!---------------------------------------------------------------------
!     Complete element centroid calculation
!---------------------------------------------------------------------
!
      if(elvol.le.1.d-20) then

         write(*,'(A,I10)')
     &      ' *ERROR in linearCompliance_C3D10: zero volume for ',
     &      nelem

         nmethod=0
         return

      endif

      xcg=xmoment/elvol
      ycg=ymoment/elvol
      zcg=zmoment/elvol

!
!---------------------------------------------------------------------
!     Return the unpenalized element stiffness matrix
!
!     s = K0_e
!
!     Density penalization should be applied exactly once. If the
!     assembly routine already multiplies s by rhoi**penal, do not
!     apply the factor here.
!---------------------------------------------------------------------
!
      do i=1,ndof
         do j=1,ndof
            s(i,j)=ke(i,j)
         enddo
      enddo

!
!---------------------------------------------------------------------
!     Compute unpenalized elemental strain energy
!
!       uku = u_e^T K0_e u_e
!---------------------------------------------------------------------
!
      uku=0.d0

      do i=1,ndof
         do j=1,ndof

            uku=uku+
     &          uelem(i)*ke(i,j)*uelem(j)

         enddo
      enddo

!
!---------------------------------------------------------------------
!     SIMP compliance and compliance sensitivity
!
!       C_e = rho_e^p u_e^T K0_e u_e
!
!       dC/drho_e =
!          -p rho_e^(p-1) u_e^T K0_e u_e
!---------------------------------------------------------------------
!
      rhoSafe=max(rhoi,1.d-12)

      ecompli=(rhoSafe**penal)*uku

      sensi=-penal*
     &       (rhoSafe**(penal-1.d0))*
     &       uku

      return
      end