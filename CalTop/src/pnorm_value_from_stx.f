      subroutine pnorm_value_from_stx(co,kon,ipkon,lakon,ne,
     &   stx,mi,design,penal,sig0,eps_relax,rho_min,pexp,
     &   nea,neb,list,ilist,psum)
c  Computes psum = sum_g w_g * (phi_g)^p with *frozen* Cauchy stx
      implicit none
      integer ne,kon(*),ipkon(*),mi(*),nea,neb,list,ilist(*)
      character*8 lakon(*),lakonl
      real*8 co(3,*),stx(6,mi(1),*),design(*)
      real*8 penal,sig0,eps_relax,rho_min,pexp,psum

      integer i,k,j,indexe,nope,mint3d,iflag,jj,konl(26)
      real*8 xl(3,26),shp(4,26),xsj,xi,et,ze,weight
      real*8 rho,rho_eff,rho_p,phi,vm2,vm
      real*8 sx,sy,sz,txy,txz,tyz,wgt

      include 'gauss.f'

      psum = 0.d0
      iflag = 3

      do k=nea,neb
        if (list.eq.1) then
          i = ilist(k)
        else
          i = k
        endif
        if (ipkon(i).lt.0) cycle
        lakonl = lakon(i)
        if     (lakonl(4:5).eq.'20') then
          nope=20; mint3d=8
        elseif (lakonl(4:4).eq.'8')  then
          nope=8;  mint3d=8
        elseif (lakonl(4:5).eq.'10') then
          nope=10; mint3d=4
        elseif (lakonl(4:4).eq.'4')  then
          nope=4;  mint3d=1
        elseif (lakonl(4:5).eq.'15') then
          nope=15; mint3d=9
        elseif (lakonl(4:4).eq.'6')  then
          nope=6;  mint3d=2
        else
          cycle
        endif

        indexe=ipkon(i)
        do j=1,nope
          konl(j)=kon(indexe+j)
          xl(1,j)=co(1,konl(j))
          xl(2,j)=co(2,konl(j))
          xl(3,j)=co(3,konl(j))
        enddo

        rho = design(i)
        if (rho .lt. 0.d0) rho = 0.d0
        if (rho .gt. 1.d0) rho = 1.d0
        rho_eff = dmax1(rho, rho_min)
        rho_p   = rho_eff**penal

        do jj=1,mint3d
          if     (nope.eq.20 .or. nope.eq.8) then
            xi=gauss3d2(1,jj); et=gauss3d2(2,jj); ze=gauss3d2(3,jj)
            weight=weight3d2(jj)
            if (nope.eq.20) then
              call shape20h(xi,et,ze,xl,xsj,shp,iflag)
            else
              call shape8h(xi,et,ze,xl,xsj,shp,iflag)
            endif
          elseif (nope.eq.10) then
            xi=gauss3d5(1,jj); et=gauss3d5(2,jj); ze=gauss3d5(3,jj)
            weight=weight3d5(jj); 
            call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
          elseif (nope.eq.4) then
            xi=gauss3d4(1,jj); et=gauss3d4(2,jj); ze=gauss3d4(3,jj)
            weight=weight3d4(jj); 
            call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
          elseif (nope.eq.15) then
            xi=gauss3d8(1,jj); et=gauss3d8(2,jj); ze=gauss3d8(3,jj)
            weight=weight3d8(jj); 
            call shape15w(xi,et,ze,xl,xsj,shp,iflag)
          else
            xi=gauss3d7(1,jj); et=gauss3d7(2,jj); ze=gauss3d7(3,jj)
            weight=weight3d7(jj); 
            call shape6w(xi,et,ze,xl,xsj,shp,iflag)
          endif

          sx  = stx(1,jj,i); sy  = stx(2,jj,i); sz  = stx(3,jj,i)
          txy = stx(4,jj,i); txz = stx(5,jj,i); tyz = stx(6,jj,i)

          vm2 = (sx-sy)*(sx-sy) + (sy-sz)*(sy-sz) + (sz-sx)*(sz-sx)
          vm2 = 0.5d0*vm2 + 3.d0*(txy*txy + txz*txz + tyz*tyz)
          if (vm2 .le. 0.d0) cycle
          vm  = dsqrt(vm2)

          phi = vm/(rho_p*sig0) + eps_relax - eps_relax/rho_eff
          if (phi .le. 0.d0) cycle

          wgt = xsj*weight
          psum = psum + (phi**pexp)*wgt
        enddo
      enddo

      return
      end
