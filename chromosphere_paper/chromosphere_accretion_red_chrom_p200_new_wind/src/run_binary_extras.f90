! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use math_lib
      
      implicit none
      
      contains
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          b% warn_binary_extra =.false.
	b% other_rlo_mdot => other_rlo_mdot
         
      end subroutine extras_binary_controls

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! Mdot Prescription !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get_info_for_ritter(b)
   type(binary_info), pointer :: b
   real(dp) :: rho_exponent, F1, q, rho, p, grav, hp, v_th, rl3, q_temp
   include 'formats'

   !--------------------- Optically thin MT rate -----------------------------------------------
   ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

   rho = b% s_donor% rho(1) ! density at surface in g/cm^3
   p = b% s_donor% Peos(1) ! pressure at surface in dynes/cm^2
   grav = standard_cgrav*b% m(b% d_i)/pow2(b% r(b% d_i)) ! local gravitational acceleration
   hp = p/(grav*rho) ! pressure scale height
   v_th = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1)))

   q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                 ! (Kolb & Ritter 1990 use the opposite!)
   ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
   q_temp = min(max(q,0.5d0),10d0)
   F1 = (1.23d0  + 0.5D0* log10(q_temp))
   rl3 = (b% rl(b% d_i))*(b% rl(b% d_i))*(b% rl(b% d_i))
   b% mdot_thin0 = (2.0D0*pi/exp(0.5d0)) * v_th*v_th*v_th * &
       rl3/(standard_cgrav*b% m(b% d_i)) * rho * F1
   !Once again, do not extrapolate! Eq. (7) of Ritter 1988
   q_temp = min(max(q,0.04d0),20d0)
   if (q_temp < 1.0d0) then
      b% ritter_h = hp/( 0.954D0 + 0.025D0*log10(q_temp) - 0.038D0*pow2(log10(q_temp)) )
   else
      b% ritter_h = hp/( 0.954D0 + 0.039D0*log10(q_temp) + 0.114D0*pow2(log10(q_temp)) )
   end if

   b% ritter_exponent = (b% r(b% d_i)-b% rl(b% d_i))/b% ritter_h

   if (b% mdot_scheme == "Kolb") then
      if (b% ritter_exponent > 0) then
         b% mdot_thin = -b% mdot_thin0
      else
         b% mdot_thin = -b% mdot_thin0 * exp(b% ritter_exponent)
      end if
   else
      b% mdot_thin = -b% mdot_thin0 * exp(b% ritter_exponent)
   end if

end subroutine get_info_for_ritter

real(dp) function calculate_kolb_mdot_thick(b, indexR, rl_d) result(mdot_thick)
   real(dp), intent(in) :: rl_d
   integer, intent(in) :: indexR
   real(dp) :: F1, F3, G1, dP, q, rho, p, grav, hp, v_th, rl3, q_temp
   integer :: i
   type(binary_info), pointer :: b
   include 'formats'

   !--------------------- Optically thin MT rate -----------------------------------------------
   ! As described in Kolb and H. Ritter 1990, A&A 236,385-392
   
   ! compute integral in Eq. (A17 of Kolb & Ritter 1990)
   mdot_thick = 0d0
   do i=1,indexR-1
      G1 = b% s_donor% gamma1(i)
      F3 = sqrt(G1) * pow(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
      mdot_thick = mdot_thick + F3*sqrt(kerg * b% s_donor% T(i) / &
         (mp * b% s_donor% mu(i)))*(b% s_donor% Peos(i+1)-b% s_donor% Peos(i))
   end do
   ! only take a fraction of dP for last cell
   G1 = b% s_donor% gamma1(i)
   F3 = sqrt(G1) * pow(2d0/(G1+1d0), (G1+1d0)/(2d0*G1-2d0))
   dP = (b% s_donor% r(indexR) - rl_d) / &
      (b% s_donor% r(indexR) - b% s_donor% r(indexR+1)) * (b% s_donor% Peos(i+1)-b% s_donor% Peos(i))
   mdot_thick = mdot_thick + F3*sqrt(kerg * b% s_donor% T(i) / (mp*b% s_donor% mu(i)))*dP

   q = b% m(b% a_i)/b% m(b% d_i) ! Mass ratio, as defined in Ritter 1988
                                 ! (Kolb & Ritter 1990 use the opposite!)
   ! consider range of validity for F1, do not extrapolate! Eq. A9 of Ritter 1988
   q_temp = min(max(q,0.5d0),10d0)
   F1 = (1.23d0  + 0.5D0* log10(q_temp))
   mdot_thick = -2.0D0*pi*F1*rl_d*rl_d*rl_d/(standard_cgrav*b% m(b% d_i))*mdot_thick

end function calculate_kolb_mdot_thick

subroutine get_info_for_kolb(b)
   type(binary_info), pointer :: b
   real(dp) :: F3, FF, G1, x_L1, q, g
   real(dp) :: mdot_thick0,  R_gas, dP, rl, s_div_rl
   integer :: i, indexR
   include 'formats'

   !--------------------- Optically thick MT rate -----------------------------------------------
   ! As described in H. Ritter 1988, A&A 202,93-100 and U. Kolb and H. Ritter 1990, A&A 236,385-392

   ! First we need to find how deep inside the star the Roche lobe reaches. In other words the mesh point of the star at which R=R_RL
   b% mdot_thick = 0d0
   indexR=-1
   if(b% r(b% d_i)-b% rl(b% d_i) > 0.0d0) then
      i=1
      do while (b% s_donor% r(i) > b% rl(b% d_i))
         i=i+1
      end do
      
      if (i .eq. 1) then
         b% mdot_thick = 0d0
      else
         b% mdot_thick = calculate_kolb_mdot_thick(b, i-1, b% rl(b% d_i))
      end if
   end if

end subroutine get_info_for_kolb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Chromosphere Accretion !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function potential_rv_dimensionless(rv_a, q_input) result(phi)
    real(dp), intent(in) :: rv_a, q_input
    real(dp) :: phi
    real(dp) :: q, onebq
    real(dp) :: a1, a2, term_const, term_r
   include 'formats'

    q = 1.d0 / q_input
    onebq = 1.d0 / q

    a1 = (1.d0/3.d0) * (1.d0 + q)

    a2 = (4.d0/45.d0) * ((1.d0 + q)**2 + 9.d0*q**2 + 3.d0*q*(1.d0 + q))

    term_const = (q + q/(2.d0*(1.d0 + onebq)))

    term_r = (1.d0/rv_a + a1*rv_a**2 + a2*rv_a**5)

    phi = term_const + term_r
end function


function brent_root(f, a_in, b_in, q, fval) result(root)
    implicit none

    real(dp), intent(in) :: a_in, b_in, q, fval
    real(dp) :: root
    real(dp) :: a, b
    real(dp) :: fa, fb, c, fc, tol
    integer :: iter
   include 'formats'

    interface
        function f(x, q)
            real(8) :: f
            real(8), intent(in) :: x, q
        end function f
    end interface

    tol = 1.d-10

    a = a_in
    b = b_in

    fa = f(a, q) - fval
    fb = f(b, q) - fval

   if (fa*fb > 0.0) then
       print *, "Error: root not bracketed"
       stop
   end if

    do iter = 1, 100
        c = 0.5d0*(a + b)
        fc = f(c, q) - fval

        if (abs(fc) < tol) exit

        if (fa*fc < 0.d0) then
            b = c
            fb = fc
        else
            a = c
            fa = fc
        end if
    end do

    root = c
end function brent_root


!------------------------------------------------------------
function dOmega_dx(x, q) result(fval)
    real(dp), intent(in) :: x, q
    real(dp) :: fval
    real(dp) :: r1, r2, x_cm, rho
   include 'formats'

    x_cm = q / (1.d0 + q)

    r1 = abs(x)
    r2 = abs(x - 1.d0)

    rho = x - x_cm

    fval = -x/r1**3 - q*(x - 1.d0)/r2**3 + (1.d0 + q)*rho
end function


!------------------------------------------------------------
function find_L1(q_input) result(xL1)
    real(dp), intent(in) :: q_input
    real(dp) :: xL1
    real(dp) :: q
   include 'formats'

    q = 1.d0 / q_input
    xL1 = brent_root(dOmega_dx, 1.d-6, 1.d0 - 1.d-6, q, 0.0d0)
end function

!------------------------------------------------------------
function find_L2(q_input) result(xL2)
    real(dp), intent(in) :: q_input
    real(dp) :: xL2
    real(dp) :: q
   include 'formats'

    q = 1.d0 / q_input
    xL2 = brent_root(dOmega_dx, 1.d0 + 1.d-6, 5.d0, q, 0.0d0)
end function

!------------------------------------------------------------
function find_L3(q_input) result(xL3)
    real(dp), intent(in) :: q_input
    real(dp) :: xL3
    real(dp) :: q
   include 'formats'

    q = 1.d0 / q_input
    xL3 = brent_root(dOmega_dx, -1.d-6, -5.d0, q, 0.0d0)
end function


!------------------------------------------------------------
function roche_potential(x, y, z, q_input) result(phi)
    real(dp), intent(in) :: x, y, z, q_input
    real(dp) :: phi
    real(dp) :: q, x_cm, r1, r2, rho2
   include 'formats'

    q = 1.d0 / q_input

    x_cm = q / (1.d0 + q)

    r1 = sqrt(x**2 + y**2 + z**2)
    r2 = sqrt((x - 1.d0)**2 + y**2 + z**2)

    rho2 = (x - x_cm)**2 + y**2

    phi = 1.d0/r1 + q/r2 + 0.5d0*(1.d0 + q)*rho2
end function


subroutine rl_area_merchant21_dimensionless(q_input, areaB, areaC, C3, C4, C5)
    implicit none

    ! Inputs
    real(dp), intent(in)  :: q_input
    ! Outputs
    real(dp), intent(out) :: areaB, areaC, C3, C4, C5
    ! Local variables
    real(dp) :: q, xL1
    real(dp) :: gq, gq_higher

    ! Compute q
    q = 1.0d0 / q_input

    ! Compute L1 position
    xL1 = find_L1(q_input)

    ! Intermediate quantities
    gq = (1.0d0/(xL1**3)) + (q/((1.0d0 - xL1)**3))
    gq_higher = (1.0d0/(xL1**5)) + (q/((1.0d0 - xL1)**5))

    ! Final outputs
    areaB = (gq - (q + 1.0d0)) / 2.0d0
    areaC = gq / 2.0d0
    C3 = -3.0d0 * gq_higher / 8.0d0
    C4 = 2.0d0 * C3
    C5 = C3

end subroutine rl_area_merchant21_dimensionless



function get_info_for_chromacc(b, vconv) result(mdot_cgs)
    ! Inputs
    type(binary_info), pointer :: b
    real(dp), intent(in) :: vconv

    ! Output
    real(dp) :: mdot_cgs

    ! Scalars
    real(dp) :: rhoph, Rd, Md, vT, q, sma, r_L1, vs_lim, phi_max
    real(dp) :: pot_fact_cm2s2
    real(dp) :: Rd_prime, r_L1_prime
    real(dp) :: vT_prime, vconv_prime
    real(dp) :: xl1, xl2, xl3
    real(dp) :: phi_l2, l1_pot, phi_phot, phi_integ1, phi_integ2, phi_l1l2_mean, phi_l3
    real(dp) :: phi, phi_plus
    real(dp) :: rphi, rphi_plus, rphi_mean
    real(dp) :: areaB, areaC, C3, C4, C5, areafact1, areafact2

    ! Velocity grid
    real(dp) :: vsj, vsjp1, vs_grid_centre
    real(dp) :: vs_grid_centre2, v_total2
    integer :: phi_res, vs_res

    ! Factors
    real(dp) :: fact1, fact2, fact3, fact4, fact5
    real(dp) :: area_exp, darea, darea_dphi

    ! Integration
    real(dp) :: integral_numerical, vs_integral_numerical

    ! Loop indices
    integer :: i, j

   include 'formats'

    ! Initialize accumulators
    integral_numerical = 0.0_dp

   vs_res = 1000
    
   rhoph = b% s_donor% rho(1) ! cgs
   Rd = b% r(b% d_i)
   Md = b% m(b% d_i)
   vT = sqrt(kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1))) ! cgs
   q = b% m(b% d_i)/b% m(b% a_i)
   sma = b% separation ! cm
   r_L1 = b% rl(b% d_i) ! cm

   pot_fact_cm2s2 = (standard_cgrav * Md) / sma ! cm2/s2
   Rd_prime = Rd / sma
   r_L1_prime = r_L1 / sma
   vT_prime = vT / sqrt(pot_fact_cm2s2)
   vconv_prime = vconv / sqrt(pot_fact_cm2s2)
   
   xl1 = find_L1(q)
   xl2 = find_L2(q)
   xl3 = find_L3(q)
   phi_l2 = -roche_potential(xl2, 0.0d0, 0.0d0, q_input=q)
   phi_l3 = -roche_potential(xl3, 0.0d0, 0.0d0, q_input=q)

   l1_pot = -roche_potential(xl1, 0.0d0, 0.0d0, q)

   phi_max = min(phi_l2, phi_l3)
   phi_l1l2_mean = (min(phi_l2, phi_l3)+l1_pot)/2.0

   phi_phot = -potential_rv_dimensionless(Rd_prime, q)

   call rl_area_merchant21_dimensionless(q, areaB, areaC, C3, C4, C5)

   areafact1 = pi / sqrt(areaB * areaC)
   areafact2 = areafact1 * (C3*(3.0d0*areaB**2 + 2.0d0*areaB*areaC + 3.0d0*areaC**2) &
            /(4.0d0*areaB**2 * areaC**2))

   rphi_mean = brent_root(potential_rv_dimensionless, r_L1_prime, 2*r_L1_prime, q, -phi_l1l2_mean)

   ! vs_lim = 0.7 * (min(phi_l2, phi_l3) - l1_pot)
   vs_lim = 10 * vconv_prime

   integral_numerical = 0.0d0
    do j = 1, vs_res
         vsj = (j-1) * vs_lim / (vs_res-1)
         vsjp1 = (j) * vs_lim / (vs_res-1)

         vs_grid_centre = 0.5d0 * (vsj + vsjp1)

         vs_grid_centre2 = vs_grid_centre**2
         v_total2 = vs_grid_centre2 + vT_prime**2

         fact1 = exp(-vs_grid_centre2/(2.0d0*vconv_prime**2))
         fact2 = sqrt(v_total2) * (vsjp1 - vsj)

         darea_dphi = (areafact1 - areafact2*(phi_l1l2_mean - l1_pot))

         phi_integ1 = exp(-max((l1_pot - phi_phot), 0.0d0)/v_total2)
         phi_integ2 = exp(-max((phi_max - phi_phot), 0.0d0)/v_total2)
         fact4 = darea_dphi * ((v_total2 * (phi_integ1 - phi_integ2)) + max((phi_phot - l1_pot), 0.0d0) - max((phi_phot - phi_max), 0.0d0))

         area_exp = vs_grid_centre2/v_total2
         fact5 = (Rd_prime/rphi_mean)**(2.0d0*area_exp)

         integral_numerical = integral_numerical + (fact1 * fact2 * fact5 * fact4)

      end do

   mdot_cgs = -1.0d0 * 0.606d0 * sqrt(2.0d0/3.141593d0) * (1.0d0/vconv) * rhoph * integral_numerical * pot_fact_cm2s2 * (sma**2)

   print*,"Mdot:", mdot_cgs*secyer/Msun
   print*,"vT Prime:", vT_prime
   print*,"Vcon Prime:", vconv_prime
   print*,"Rd Prime:", Rd_prime
   print*,"Rd by Rl:", Rd/r_L1
   print*,"q:", q
   print*,"vconv:", vconv/1.0d5
   print*,"rhoph:", rhoph
   print*,"sma:", sma/Rsun
   print*,"pot_fact_km2s2:", pot_fact_cm2s2/1.0d10

end function get_info_for_chromacc


function get_info_for_ritter_revised(b) result(mdot_cgs_ritter)

   type(binary_info), pointer :: b

   real(dp) :: mdot_cgs_ritter
   real(dp) :: rhoph, Rd, Md, vT, q, sma
   real(dp) :: areaB, areaC, C3, C4, C5
   real(dp) :: pot_fact_cm2s2, fact_area
   real(dp) :: xl1, l1_pot, phi_phot, phi_diff, Rd_prime

   real(dp), parameter :: pi_dp = 3.1415926535897932384626433832795_dp

   rhoph = b%s_donor%rho(1)  ! cgs
   Rd = b%r(b%d_i)
   Md = b%m(b%d_i)

   vT = sqrt(kerg * b%s_donor%T(1) / (mp * b%s_donor%mu(1)))  ! cgs

   q = b%m(b%d_i) / b%m(b%a_i)
   sma = b%separation  ! cm
   Rd_prime = Rd / sma

   pot_fact_cm2s2 = standard_cgrav * Md / sma  ! cm^2/s^2

   call rl_area_merchant21_dimensionless(q, areaB, areaC, C3, C4, C5)

   fact_area = pi_dp * sma**3 / &
      (sqrt(areaB * areaC) * standard_cgrav * Md)

   xl1 = find_L1(q)

   l1_pot = -potential_rv_dimensionless(b%rl(b%d_i)/sma, q)
   phi_phot = -potential_rv_dimensionless(Rd_prime, q)

   phi_diff = max(l1_pot - phi_phot, 0.0_dp) * pot_fact_cm2s2

   mdot_cgs_ritter = -1.0d0 * exp(-0.5_dp) * vT**3 * &
      exp(-phi_diff / vT**2) * rhoph * fact_area

   print*,"Mdot:", mdot_cgs_ritter*secyer/Msun
   print*,"Rd:", Rd/Rsun
   print*,"Md:", Md/Msun
   print*,"Ma:", b%m(b%a_i)/Msun
   print*,"sma:", sma/Rsun
   print*,"T:", b%s_donor%T(1)
   print*,"mu:", b%s_donor%mu(1)


end function get_info_for_ritter_revised



function get_vconv(s) result(vconmax)
   type(star_info), pointer :: s
   real(dp) :: vconmax, rconmax, tauconmax, Mconmax, rhoconmax
   integer :: tauflag, k
   include 'formats'

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
end function get_vconv


subroutine other_rlo_mdot(binary_id, mdot, ierr)
   use binary_def, only : binary_info, binary_ptr
   use const_def, only: dp
   integer, intent(in) :: binary_id
   real(dp), intent(out) :: mdot
   real(dp) :: vconv
   integer, intent(out) :: ierr
   type (binary_info), pointer :: b
   ierr = 0
   call binary_ptr(binary_id, b, ierr)
   if (ierr /= 0) then
      write(*,*) 'failed in binary_ptr'
      return
   end if


   mdot = 0d0

   if (b% mdot_scheme == "Kolb" .and. b% eccentricity <= 0.0d0) then
      !call get_info_for_ritter(b)
      !mdot = b% mdot_thin
      mdot = get_info_for_ritter_revised(b)
      b% mdot_thin = mdot
      call get_info_for_kolb(b)
      mdot = mdot + b% mdot_thick

   else if (b% mdot_scheme == "Chromosphere" .and. b% eccentricity <= 0.0d0) then
      vconv = get_vconv(b% s_donor)
      mdot = get_info_for_chromacc(b, vconv/2.0)
      b% mdot_thin = mdot
      call get_info_for_kolb(b)
      mdot = mdot + b% mdot_thick

   end if

end subroutine other_rlo_mdot



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! End Mdot Prescription !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items


      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
!          b% s1% job% warn_run_star_extras = .false.
          extras_binary_startup = keep_going
      end function  extras_binary_startup
      
      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr

         extras_binary_start_step = keep_going
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
      
      end function  extras_binary_start_step
      
      !Return either keep_going, retry or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_check_model = keep_going
        
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going
         
      end function extras_binary_finish_step
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if      
         
 
      end subroutine extras_binary_after_evolve



      
      end module run_binary_extras
