! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras 

      use star_lib
      use star_def
      use const_def
      use math_lib
      use chem_def
      use num_lib
      use binary_def
      
      implicit none

      ! these routines are called by the standard run_star check_model
      contains
      
      ! include 'standard_run_star_extras.inc'

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve







      subroutine get_internal_energies(s, energy_tot, egas_from_Pgas_tot, erad_from_Prad_tot, egrav_tot, energy_env, egas_from_Pgas_env, erad_from_Prad_env, egrav_env, ekin_env, energy_usr, egas_from_Pgas_usr, erad_from_Prad_usr, egrav_usr)

            !!! Function used to calculate envelope binding energy, see Yamaguchi et al. (2024) 
            
            real(dp), intent(out) :: energy_tot, egas_from_Pgas_tot, erad_from_Prad_tot, egrav_tot, energy_env, egas_from_Pgas_env, erad_from_Prad_env, egrav_env, ekin_env, energy_usr, egas_from_Pgas_usr, erad_from_Prad_usr, egrav_usr
            type (star_info), pointer :: s
            integer :: i, j

            energy_tot = 0d0
            egas_from_Pgas_tot = 0d0
            erad_from_Prad_tot = 0d0
            egrav_tot = 0d0

            energy_env = 0d0
            egas_from_Pgas_env = 0d0 
            erad_from_Prad_env = 0d0
            egrav_env = 0d0
            ekin_env = 0d0

            energy_usr = 0d0
            egas_from_Pgas_usr = 0d0
            erad_from_Prad_usr = 0d0
            egrav_usr = 0d0

            do j = s% nz, 1, -1   ! go down in steps of -1 from core until the surface is reached 

                  energy_tot = energy_tot + s % energy(j) * s % dm(j)
                  egrav_tot = egrav_tot - s% cgrav(j) * s% dm(j) * s% m(j) / s% r(j)
                  egas_from_Pgas_tot = egas_from_Pgas_tot + s % dm(j) * (3 * s% Pgas(j))/(2 * s% rho(j))
                  erad_from_Prad_tot = erad_from_Prad_tot + s % dm(j) * (3 * s% Prad(j))/(s% rho(j))
                  if (s% m(j) .ge. s% he_core_mass * 1.99D33) then 
                        energy_env = energy_env + s % energy(j) * s % dm(j)
                        egrav_env = egrav_env - s% cgrav(j) * s% dm(j) * s% m(j) / s% r(j)
                        ekin_env = ekin_env + 0.5d0 * s% dm(j) * s% conv_vel(j) * s% conv_vel(j)
                        egas_from_Pgas_env = egas_from_Pgas_env +  s % dm(j) * (3 * s% Pgas(j))/(2 * s% rho(j))
                        erad_from_Prad_env = erad_from_Prad_env +  s % dm(j) * (3 * s% Prad(j))/(s% rho(j))
                  endif
                  if (s% r(j) .ge. s%  x_ctrl(4) * 1.49597870700D13) then
                        energy_usr = energy_usr + s % energy(j) * s % dm(j)
                        egrav_usr = egrav_usr - s% cgrav(j) * s% dm(j) * s% m(j) / s% r(j)
                        egas_from_Pgas_usr = egas_from_Pgas_usr +  s % dm(j) * (3 * s% Pgas(j))/(2 * s% rho(j))
                        erad_from_Prad_usr = erad_from_Prad_usr +  s % dm(j) * (3 * s% Prad(j))/(s% rho(j))
                  endif
                  
            end do

      end subroutine get_internal_energies


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 14
      end function how_many_extra_history_columns

      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         integer :: i

         real(dp) :: energy_tot, egas_from_Pgas_tot, erad_from_Prad_tot, egrav_tot
      real(dp) :: energy_env, egas_from_Pgas_env, erad_from_Prad_env, egrav_env, ekin_env
      real(dp) :: energy_usr, egas_from_Pgas_usr, erad_from_Prad_usr, egrav_usr
      real(dp) :: vconmax, rconmax, tauconmax, Mconmax, rhoconmax, tauflag

      integer :: k


         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


      ! Assume s is already associated with a valid star model

      call get_internal_energies(s, energy_tot, egas_from_Pgas_tot, erad_from_Prad_tot, egrav_tot, &
                              energy_env, egas_from_Pgas_env, erad_from_Prad_env, egrav_env, ekin_env, &
                              energy_usr, egas_from_Pgas_usr, erad_from_Prad_usr, egrav_usr)

         names(1) = "egrav_env"
         vals(1) = egrav_env  ! in SI units

         names(2) = "energy_env"
         vals(2) = energy_env  ! in SI units

         names(3) = "egas_from_Pgas_env"
         vals(3) = egas_from_Pgas_env  ! in SI units

         names(4) = "erad_from_Prad_env"
         vals(4) = erad_from_Prad_env  ! in SI units

         names(5) = "ekin_env"
         vals(5) = ekin_env  ! in SI units

         names(6) = 'rho_surf'
	 vals(6) = s% rho(1)

         tauflag = 0
	 rconmax = 0d0
	 vconmax = 1d0
	 tauconmax = 0d0
         do k = 1, s% nz-1, 1
		if (tauflag < 0.5) then
			if ( s% tau(k) > clight*s% Prad(k)/(s% conv_vel(k)*s% Pgas(k)) ) then
				vconmax = min(s% conv_vel(k),1.2*s% csound(k))
				rconmax = s% r(k)
				tauconmax = s% tau(k)
				Mconmax = vconmax/s% csound(k)
				rhoconmax = s% rho(k)
				tauflag = 1
			end if
		end if
         end do

         names(7) = 'vcon_surf'
	 vals(7) = vconmax

         names(8) = 'rcon'
	 vals(8) = rconmax

         names(9) = 'taucon'
	 vals(9) = tauconmax

         names(10) = 'Mcon_surf'
	 vals(10) = Mconmax

         names(11) = 'rhocon_surf'
	 vals(11) = rhoconmax

         names(12) = 'cs_surf'
	 vals(12) = s% csound(1)

    names(13) = 'Peos_surf'
	 vals(13) = s% Peos(1)

    names(14) = 'T_surf'
	 vals(14) = s% T(1)
         

      end subroutine data_for_extra_history_columns


      end module run_star_extras
      
