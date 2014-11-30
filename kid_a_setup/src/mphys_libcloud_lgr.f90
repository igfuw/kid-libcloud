!Interface to libcloudph++: lagrangian scheme (super-droplet) 

module mphys_libcloud_lgr

  Use parameters, only : nz, dt, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi

  Use diagnostics, only: save_dg, i_dgtime
  use iso_c_binding, only: c_funptr, c_f_procpointer, c_null_char, c_double

  Implicit None
  
  interface
    subroutine micro_step_py(i_dgtime, size_z, size_x,  & 
                        th_ar, qv_ar, rhof_ar, rhoh_ar, &
                        vf_ar, vh_ar, wf_ar, wh_ar,     &
                        xf_ar, zf_ar, xh_ar, zh_ar, tend_th_ar, tend_qv_ar) bind(c)
      use iso_c_binding, only: c_double, c_int
      Use parameters, only : nx, nz
      integer(c_int), intent(in), value :: i_dgtime, size_x, size_z
                                           
                                           
      real(c_double),  intent(inout):: th_ar(nz, 0:nx+1), qv_ar(nz, 0:nx+1),  &
                                       rhof_ar(nz), rhoh_ar(nz),              &
                                       vh_ar(nz, 0:nx+1), vf_ar(nz, 0:nx+1),  &
                                       wf_ar(nz, 0:nx+1), wh_ar(nz, 0:nx+1),  &
                                       xf_ar(0:nx+1), zf_ar(nz), xh_ar(0:nx+1), zh_ar(nz), &
                                       tend_th_ar(nz, 0:nx+1), tend_qv_ar(nz, 0:nx+1)

    end

    subroutine load_ptr(fname, ptr) bind(c)
      use iso_c_binding, only: c_funptr, c_char
      character(kind=c_char), dimension(*), intent(in) :: fname
      type(c_funptr), intent(out) :: ptr
    end
  end interface

  type(c_funptr) :: cptr
  procedure(micro_step_py), pointer :: fptr => NULL()

contains

  Subroutine mphys_libcloud_lgr_interface

    ! do the below only once
    if (associated(fptr) .eqv. .false.) then 
      ! assert for numerical precision  
      if (wp.ne.c_double) stop("KiD does not use double precision!")          

      ! load pointer to Python micro_step() routine
      call load_ptr("/tmp/micro_step.ptr" // c_null_char,cptr)

      ! associate the C pointer with the F pointer
      call c_f_procpointer(cptr, fptr)
    end if

    ! do the below every timestep
    call fptr(i_dgtime, nz, nx+2 , &
              theta, qv, rho, rho_half, & 
              v, v_half, w, w_half, x, z, x_half, z_half, dTheta_mphys, dqv_mphys)
       ! TODO: pass dt 
  
  end Subroutine mphys_libcloud_lgr_interface

end module mphys_libcloud_lgr
