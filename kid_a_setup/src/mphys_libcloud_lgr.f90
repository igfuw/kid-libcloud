
!Interface to libcloudph++: lagrangian scheme (super-droplet) 

module mphys_libcloud_lgr

  Use parameters, only : num_h_moments, num_h_bins, nspecies, nz, dt &
       , h_names, mom_units, max_char_len, mom_names, nx
  Use column_variables
  Use physconst, only : p0, r_on_cp, pi

  Use diagnostics, only: save_dg, i_dgtime ! do we need it?

  Implicit None

contains

  Subroutine mphys_libcloud_lgr_interface
    print*, "in mphys_libcloud_lgr_interface"
  end Subroutine mphys_libcloud_lgr_interface

end module mphys_libcloud_lgr
