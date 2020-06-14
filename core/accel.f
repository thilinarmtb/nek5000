c-----------------------------------------------------------------------
      subroutine accel_setup(glo_num,ierr)
      include 'SIZE'
      include 'TOTAL'
      include 'ACCEL'

      include 'exaf.h'

      integer*8 glo_num(1)
      integer ierr

      common /nekmpi/mid,mp,nekcomm,nekgroup,nekreal

      real geom(lx1*ly1*lz1,7,lelt)

      integer n_tot

      call read_state(nekcomm)

      call exainit('/occa/gpu/cuda',nekcomm,exa_h,ierr)

      call exasettingscreate(exa_h,exa_str_null,settings_h,ierr)
      call exasettingssetint('general::order',nx1,settings_h,ierr)

      ! create and setup hmholtz
      call exahmholtzcreate(exa_hmhz_h,exa_h,ierr)
      call exahmholtzsetup(exa_hmhz_h,settings_h,ierr)

      ! create and setup mesh
      call exameshcreate(mesh_h,exa_str_null,exa_h,ierr)

      call exameshset1ddofs(mesh_h,nx1,ierr)
      call exameshsetnelements(mesh_h,nelt,ierr)
      call exameshsetdim(mesh_h,ndim,ierr)
      call exameshsetelemx(mesh_h,xc,ierr)
      call exameshsetelemy(mesh_h,yc,ierr)
      call exameshsetelemz(mesh_h,zc,ierr)
      call exameshsetmeshx(mesh_h,xm1,ierr)
      call exameshsetmeshy(mesh_h,ym1,ierr)
      call exameshsetmeshz(mesh_h,zm1,ierr)

      call exameshsetglobalids(mesh_h,glo_num,ierr)

      if(ndim.eq.2) then
        do e=1,nelt
        do i=1,nx1*ny1
          geom(i,1,e)=g1m1 (i,1,1,e)
          geom(i,2,e)=g2m1 (i,1,1,e)
          geom(i,3,e)=g4m1 (i,1,1,e)
          geom(i,4,e)=jacm1(i,1,1,e)
        enddo
        enddo
      elseif(ndim.eq.3) then
        do e=1,nelt
        do i=1,nx1*ny1*nz1
          geom(i,1,e)=g1m1 (i,1,1,e)
          geom(i,2,e)=g4m1 (i,1,1,e)
          geom(i,3,e)=g5m1 (i,1,1,e)
          geom(i,4,e)=g2m1 (i,1,1,e)
          geom(i,5,e)=g6m1 (i,1,1,e)
          geom(i,6,e)=g3m1 (i,1,1,e)
          geom(i,7,e)=jacm1(i,1,1,e)
        enddo
        enddo
      endif
      call exameshsetgeometricfactors(mesh_h,geom,ierr)

      call exameshsetmask(mesh_h,v1mask,ierr)
      call exameshsetderivativematrix(mesh_h,dxm1,ierr)

      call exameshsetup(mesh_h,settings_h,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine accel_cg(u,r,tol,maxit,verbose,ierr)
      include 'SIZE'
      include 'ACCEL'
      include 'exaf.h'

c     it is a bug if sol was declared as sol(1)
      real*8 u(1),r(1),sol(nx1*ny1*nz1*nelt)
      real tol
      integer maxit,verbose,ierr

      integer d_u,d_r
      integer*8 in_offset,out_offset
      parameter(in_offset=0)

      integer nt
      nt=nx1*ny1*nz1*nelt

      call exavectorcreate(exa_h,nt,exa_scalar,d_u,ierr)

      call exavectorcreate(exa_h,nt,exa_scalar,d_r,ierr)
      call exavectorwrite (d_r,r,in_offset,ierr)

      call exahmholtzcg(exa_hmhz_h,d_u,d_r,mesh_h,tol,
     $  maxit,verbose,ierr)

      call exavectorread(d_u,sol,out_offset,ierr)
      do i=1,nt
        u(i)=sol(i+out_offset)
      enddo

      call exavectorfree(d_u,ierr)
      call exavectorfree(d_r,ierr)
      end
c-----------------------------------------------------------------------
      subroutine accel_finalize(ierr)
      include 'ACCEL'
      include 'exaf.h'

      integer ierr

      call exasettingsfree(settings_h,ierr)
      call exameshdestroy(mesh_h,ierr)
      call exahmholtzdestroy(exa_hmhz_h,ierr)
      call exafinalize(exa_h,ierr)

      return
      end
c-----------------------------------------------------------------------
