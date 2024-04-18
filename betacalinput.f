!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine betacalinput(co,nk,kon,ipkon,lakon,nkon,
     &  ne,nodeboun,ndirboun,xboun,nboun,
     &  ipompc,nodempc,coefmpc,nmpc,nmpc_,nodeforc,ndirforc,xforc,nforc,
     &  nforc_,nelemload,sideload,xload,nload,nload_,
     &  nprint,prlab,prset,mpcfree,nboun_,
     &  mei,set,istartset,iendset,ialset,nset,nalset,elcon,nelcon,rhcon,
     &  nrhcon,alcon,nalcon,alzero,t0,t1,
     &  matname,ielmat,orname,orab,ielorien,amname,amta,namta,nam,
     &  nmethod,iamforc,iamload,iamt1,
     &  ithermal,iperturb,istat,istep,nmat,ntmat_,norien,
     &  prestr,iprestr,isolver,fei,veold,timepar,
     &  xmodal,filab,jout,nlabel,idrct,jmax,
     &  iexpl,alpha,iamboun,plicon,nplicon,plkcon,
     &  nplkcon,iplas,npmat_,mi,nk_,trab,inotr,ntrans,ikboun,
     &  ilboun,ikmpc,ilmpc,ics,dcs,ncs_,namtot_,cs,nstate_,ncmat_,iumat,
     &  mcs,labmpc,iponor,xnor,knor,thickn,thicke,ikforc,ilforc,
     &  offset,iponoel,inoel,rig,infree,nshcon,shcon,cocon,ncocon,
     &  physcon,nflow,ctrl,maxlenmpc,ne1d,
     &  ne2d,nener,vold,nodebounold,ndirbounold,xbounold,
     &  xforcold,xloadold,t1old,eme,sti,ener,xstate,jobnamec,
     &  irstrt,ttime,qaold,output,typeboun,inpc,
     &  ipoinp,inp,tieset,tietol,ntie,fmpc,cbody,ibody,xbody,
     &  nbody,nbody_,xbodyold,nam_,ielprop,nprop,nprop_,prop,itpamp,
     &  iviewfile,ipoinpc,nslavs,t0g,t1g,network,cyclicsymmetry,
     &  idefforc,idefload,idefbody,mortar,ifacecount,islavsurf,
     &  pslavsurf,clearini,heading,iaxial,nobject,objectset,nprint_,
     &  iuel,nuel_,nodempcref,coefmpcref,ikmpcref,memmpcref_,
     &  mpcfreeref,maxlenmpcref,memmpc_,isens,namtot,nstam,dacon,
     &  vel,nef,velo,veloo)
!
      implicit none
!
!     nmethod: -1:visco (=static+creep) 
!               0:no analysis 
!               1:static
!               2:frequency 
!               3:buckling 
!               4:linear dynamic
!               5:steady state dynamics
!               6:Coriolis frequency calculation
!               7:flutter frequency calculation
!               8:magnetostatics
!               9:magnetodynamics (inductive heating)
!               10:electromagnetic eigenvalue problems
!               11:superelement creation
!               12:sensitivity
!     iprestr: 0: no residual stresses; 1: residual stresses;
!              2; residual strains
!     iperturb: 0:no perturbation; 1:perturbation; 2: nonlinear
!               geometric analysis; 3: material and geometrical
!               nonlinearities
!     istep: step number
!     istat: error indicator (<0:EOF, >0: input error)
!
      logical boun_flag,cload_flag,dload_flag,temp_flag,elprint_flag,
     &  nodeprint_flag,elfile_flag,nodefile_flag,contactfile_flag,
     &  dflux_flag,cflux_flag,film_flag,radiate_flag,out3d,
     &  solid,sectionprint_flag,contactprint_flag,pretension,
     &  beamgeneralsection,objective_flag,constraint_flag
!
      character*1 typeboun(*),inpc(*)
      character*3 output
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*)
      character*66 heading(*)
      character*80 matname(*),orname(*),amname(*)
      character*81 set(*),prset(*),tieset(3,*),cbody(*),objectset(4,*)
      character*87 filab(*)
      character*132 jobnamec(*),textpart(16)
!
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &  nodeforc(2,*),ndirforc(*),nelemload(2,*),iaxial,j,mi(*),
     &  istartset(*),iendset(*),ialset(*),ipkon(*),ics(*),nodedep,
     &  nelcon(2,*),nrhcon(*),nalcon(2,*),ielmat(mi(3),*),nodeind,
     &  ielorien(mi(3),*),icomposite,nsubmodel,mortar,
     &  namta(3,*),iamforc(*),iamload(2,*),iamt1(*),ipoinpc(0:*),
     &  iamboun(*),inotr(2,*),ikboun(*),ilboun(*),ikmpc(*),ilmpc(*),
     &  iponor(2,*),knor(*),ikforc(*),ilforc(*),iponoel(*),inoel(3,*),
     &  infree(4),ixfree,ikfree,inoelfree,iponoelmax,rig(*),nshcon(*),
     &  ncocon(2,*),nodebounold(*),ielprop(*),nprop,nprop_,maxsectors,
     &  ndirbounold(*),ipoinp(2,*),inp(3,*),nintpoint,ifacecount,
     &  ianisoplas,ifile_output,ichangefriction,nslavs,
     &  nalset,nalset_,nmat,nmat_,ntmat_,norien,norien_,islavsurf(2,*),
     &  nmethod,nk,ne,nboun,nmpc,nmpc_,mpcfree,i,istat,n,
     &  key,nk_,ne_,nboun_,ncs_,namtot_,nstate_,iviewfile,
     &  isolver,ithermal(2),iperturb(*),iprestr,istep,mei(4),nkon,
     &  nprint,nload,nload_,nforc,nforc_,nlabel,iumat,imat,
     &  nset,nset_,nprint_,nam,nam_,jout(2),ncmat_,itpamp,
     &  ierror,idrct,jmax(2),iexpl,iplas,npmat_,ntrans,ntrans_,
     &  M_or_SPC,nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),nflow,
     &  ne1d,ne2d,nener,irstrt(*),ii,maxlenmpc,inl,ipol,network,
     &  iline,mcs,ntie,ntie_,lprev,newstep,nbody,nbody_,ibody(3,*),
     &  cyclicsymmetry,idefforc(*),idefload(*),idefbody(*),
     &  ichangesurfacebehavior,nobject,ibasemotion,iuel(4,*),nuel_,
     &  nodempcref(3,*),ikmpcref(*),memmpcref_,mpcfreeref,
     &  maxlenmpcref,memmpc_,isens,iamplitudedefault,namtot,
     &  nstam,ier,nef
!
      real*8 co(3,*),xboun(*),coefmpc(*),xforc(*),fmpc(*),
     &  xload(2,*),alzero(*),offset(2,*),prop(*),pslavsurf(3,*),
     &  elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),clearini(3,9,*),
     &  alcon(0:6,ntmat_,*),thicke(mi(3),*),thickn(2,*),xnor(*),
     &  t1(*),orab(7,*),prestr(6,mi(1),*),amta(2,*),dacon(*),
     &  veold(0:mi(2),*),t0(*),plicon(0:2*npmat_,ntmat_,*),
     &  plkcon(0:2*npmat_,ntmat_,*),trab(7,*),dcs(*),
     &  shcon(0:3,ntmat_,*),cocon(0:6,ntmat_,*),timepar(*),
     &  ctrl(*),vold(0:mi(2),*),xbounold(*),xforcold(*),
     &  xloadold(*),t1old(*),eme(*),sti(*),ener(*),
     &  xstate(nstate_,mi(1),*),ttime,qaold(2),cs(17,*),tietol(2,*),
     &  xbody(7,*),xbodyold(7,*),t0g(2,*),t1g(2,*),
     &  fei(3),tinc,tper,xmodal(*),tmin,tmax,tincf,
     &  alpha,physcon(*),coefmpcref(*),vel(nef,*),velo(*),veloo(*)
!
      save solid,ianisoplas,out3d,pretension
!
      integer nentries
      parameter(nentries=17)
!
!
      return
      end
