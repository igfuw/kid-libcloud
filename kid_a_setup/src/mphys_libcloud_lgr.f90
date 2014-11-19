
!Interface to libcloudph++: lagrangian scheme (super-droplet) 

module mphys_libcloud_lgr

  Use parameters, only : num_h_moments, num_h_bins, nspecies, nz, dt &
       , h_names, mom_units, max_char_len, mom_names, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi

  Use diagnostics, only: save_dg, i_dgtime ! do we need it?
  use iso_c_binding, only: c_funptr, c_f_procpointer, c_null_char, c_double

  Implicit None

  
  interface
    subroutine micro_step_py(i_dgtime, size_z, size_x,  & 
                        th_ar, qv_ar, rho_ar, vh_ar, wh_ar) bind(c)
      use iso_c_binding, only: c_double, c_int
      Use parameters, only : nx, nz
      integer(c_int), intent(in), value :: i_dgtime, size_x, size_z
                                           
                                           
      real(c_double),  intent(inout):: th_ar(nz, 0:nx+1), qv_ar(nz, 0:nx+1),   &
                                       rho_ar(nz),                     &
                                       vh_ar(nz, 0:nx+1), wh_ar(nz, 0:nx+1)
    end

    subroutine load_ptr(fname, ptr) bind(c)
      use iso_c_binding, only: c_funptr, c_char
      character(kind=c_char), dimension(*), intent(in) :: fname
      type(c_funptr), intent(out) :: ptr
    end
  end interface

  type(c_funptr) :: cptr
  procedure(micro_step_py), pointer :: fptr

contains

  Subroutine mphys_libcloud_lgr_interface

 ! Initialise microphysics                                                             
    ! assert for numerical precision  
    if (wp.ne.c_double) stop("KiD does not use double precision!")          

    
    print*, "in mphys_libcloud_lgr_interface wp="
       call load_ptr("/tmp/micro_step.ptr" // c_null_char,cptr)
       call c_f_procpointer(cptr, fptr)
       call fptr(i_dgtime, nz, nx+2 , &
                 theta, qv, rho, v_half, w_half)
        
  
  end Subroutine mphys_libcloud_lgr_interface

end module mphys_libcloud_lgr
