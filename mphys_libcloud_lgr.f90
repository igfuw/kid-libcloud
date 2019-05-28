!Interface to libcloudph++: lagrangian scheme (super-droplet) 

module mphys_libcloud_lgr

  Use parameters, only : nz, dt, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi, Rw, R
  Use common_physics, only : qsaturation, psaturation

  Use diagnostics, only: save_dg, i_dgtime
  use iso_c_binding, only: c_funptr, c_f_procpointer, c_null_char, c_double, c_float

  Implicit None
  
  interface
    function micro_step_py(i_dgtime, dt, size_z, size_x,  & 
                        th_ar, qv_ar, rhof_ar, rhoh_ar, exner, &
                        vf_ar, vh_ar, wf_ar, wh_ar,     &
                        xf_ar, zf_ar, xh_ar, zh_ar, tend_th_ar, tend_qv_ar, rh_ar) bind(c)
      use iso_c_binding, only: c_double, c_int, c_float, c_bool
      Use parameters, only : nx, nz
      logical(c_bool) :: micro_step_py
      integer(c_int), intent(in), value :: i_dgtime, size_x, size_z
      real(c_float), intent(in), value :: dt
                                           
      real(c_double),  intent(inout):: th_ar(nz, 0:nx+1), qv_ar(nz, 0:nx+1),  &
                                       rhof_ar(nz), rhoh_ar(nz),              &
                                       exner(nz, 0:nx+1),                     &
                                       vh_ar(nz, 0:nx+1), vf_ar(nz, 0:nx+1),  &
                                       wf_ar(nz, 0:nx+1), wh_ar(nz, 0:nx+1),  &
                                       xf_ar(0:nx+1), zf_ar(nz), xh_ar(0:nx+1), zh_ar(nz), &
                                       tend_th_ar(nz, 0:nx+1), tend_qv_ar(nz, 0:nx+1), &
                                       rh_ar(nz, 0:nx+1)

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

    character(10) :: uid, pid
    integer :: k, j
    real(wp), allocatable :: RH(:,:)
    allocate(RH(nz,0:nx+1))

    ! do the below only once
    if (associated(fptr) .eqv. .false.) then 
      ! assert for numerical precision  
      if (wp.ne.c_double) stop ("KiD does not use double precision!")          
      if (sizeof(dt).ne.c_float) stop ("dt in KiD is not a float!")

      ! load pointer to Python micro_step() routine (getuid() and getpid() are GNU extensions)
      write (uid, "(I10.0)") getuid()
      write (pid, "(I10.0)") getpid()
      call load_ptr("/tmp/micro_step-" // trim(adjustl(uid)) // "-" // trim(adjustl(pid)) // ".ptr" // c_null_char,cptr)

      ! associate the C pointer with the F pointer
      call c_f_procpointer(cptr, fptr)
    end if

    ! calc RH to be passed to mphys
    do k=1,nz
       do j=0,nx+1
          if (z(k) < 20000.)then

! approx RH=qv/qvs
!             RH(k,j)= qv(k,j)/ &
!                  qsaturation(TdegK(k,j),pmb(k,j))
! exact RH = pv/pvs = rhod * qv * Rv * T / pvs
!          = pd / (T * Rd) * qv * Rv * T / pvs
!          = p0 * Ex^(cp/Rd) / Rd * qv * Rv / pvs
             RH(k,j)= 100. * p0*exner(k,j)**(1./r_on_cp) / R * &
                qv(k,j) * Rw / &
                psaturation(TdegK(k,j))


    !         RH(k,j)= 0.5
          else
             RH(k,j)=0.95 ! some better dummy var?!
          end if
       end do
    enddo

    ! do the below every timestep
    if (.not. fptr(i_dgtime, dt, nz, nx+2 , &
                   theta, qv, rho, rho_half, exner, & 
                   v, v_half, w, w_half, x, z, x_half, z_half, dTheta_mphys, dqv_mphys, RH) &
    ) stop ("Error in Python!!!")

  end Subroutine mphys_libcloud_lgr_interface

end module mphys_libcloud_lgr
