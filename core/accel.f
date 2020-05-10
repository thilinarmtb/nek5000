c-----------------------------------------------------------------------
      subroutine accel_setup(ierr)
      include 'SIZE'
      include 'TOTAL'
      include 'ACCEL'

      include 'exaf.h'

      integer ierr

      integer vertex((2**ldim)*lelt)
      integer*8 ngv,glo_num(nx1,ny1,nz1,nelt)

      common /nekmpi/mid,mp,nekcomm,nekgroup,nekreal

      real geom(7,lx1*ly1*lz1*lelt)
      integer n_tot

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

      call set_vert(glo_num,ngv,nx1,nelt,vertex,.false.)
      call exameshsetglobalids(mesh_h,glo_num,ierr)

      if(ndim.eq.2) then
        n_tot=nx1*ny1*nelt
        do i=1,n_tot
          geom(1,i)=g1m1 (i,1,1,1)
          geom(2,i)=g2m1 (i,1,1,1)
          geom(3,i)=g4m1 (i,1,1,1)
          geom(4,i)=jacm1(i,1,1,1)
        enddo
      elseif(ndim.eq.3) then
        n_tot=nx1*ny1*nz1*nelt
        do i=1,n_tot
          geom(1,i)=g1m1 (i,1,1,1)
          geom(2,i)=g2m1 (i,1,1,1)
          geom(3,i)=g3m1 (i,1,1,1)
          geom(4,i)=g4m1 (i,1,1,1)
          geom(5,i)=g5m1 (i,1,1,1)
          geom(6,i)=g6m1 (i,1,1,1)
          geom(7,i)=jacm1(i,1,1,1)
        enddo
      endif
      call exameshsetgeometricfactors(mesh_h,geom,ierr)

      call exameshsetmask(mesh_h,v1mask,ierr)
      call exameshsetderivativematrix(mesh_h,dxm1,ierr)

      call exameshsetup(mesh_h,settings_h,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine accel_cg(u,r,tol,maxit,ierr)
      include 'SIZE'
      include 'ACCEL'
      include 'exaf.h'

      integer lt
      parameter(lt=lx1*ly1*lz1*lelt)

      integer maxit,ierr
      real u(lt),r(lt)
      real tol

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
