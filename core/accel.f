      subroutine accel_setup(exa_h,exa_hmhz_h,mesh_h,settings_h)
      include 'SIZE'

      integer exa_h,exa_hmhz_h,mesh_h,settings_h,ierr

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      include 'exaf.h'

      call exainit('/occa/gpu/opencl',nekcomm,exa_h,ierr)
      call exahmholtzcreate(exa_hmhz_h,exa_h,ierr)
      call exameshcreate(mesh_h,exa_str_null,exa_h,ierr)
      call exasettingscreate(exa_h,exa_str_null,settings_h,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine accel_finalize(exa_h,exa_hmhz_h,mesh_h,settings_h)

      integer exa_h,exa_hmhz_h,mesh_h,settings_h,ierr

      include 'exaf.h'

      call exasettingsfree(settings_h,ierr)
      call exameshdestroy(mesh_h,ierr)
      call exahmholtzdestroy(exa_hmhz_h,ierr)
      call exafinalize(exa_h,ierr)

      return
      end
c-----------------------------------------------------------------------
