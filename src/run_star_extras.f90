! ***********************************************************************
!
!   Copyright (C) 2010-2025  Bill Paxton & The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module run_star_extras

   use star_lib
   use star_def
   use const_def
   use math_lib

   implicit none

   ! these routines are called by the standard run_star check_model

   real(dp), parameter :: alpha_EM = 7.2973525693d-3  ! alpha in  QED
   real(dp), parameter :: eV2_over_erg_s= 2433.99 ! eV^2/(erg/s)
   real(dp), parameter :: Kelvin_over_eV= 0.00008617 ! K/eV
   real(dp), parameter :: cm_over_eVinv= 50684.2 ! cm/(eV^{-1})
   real(dp), parameter :: cm_inv3_over_eV3 = 7.68035d-15 !cm^{-3}/eV^3

   real (dp) :: current_L_DP_L,current_L_DP_T
contains

   ! include 'standard_run_star_extras.inc'


!============================================================
! F_brem function
!
!============================================================

real(dp) function F_brem(x)

    implicit none
    real(dp), intent(in) :: x
    integer :: i
    integer, parameter :: Ntab = 101

    real(dp) :: lx
    real(dp) :: t

   real(dp), parameter :: xtab(Ntab) = (/ &
      -2.d0, -1.95d0, -1.9d0, -1.85d0, -1.8d0, &
      -1.75d0, -1.7d0, -1.65d0, -1.6d0, -1.55d0, &
      -1.5d0, -1.45d0, -1.4d0, -1.35d0, -1.3d0, &
      -1.25d0, -1.2d0, -1.15d0, -1.1d0, -1.05d0, &
      -1.d0, -0.95d0, -0.9d0, -0.85d0, -0.8d0, &
      -0.75d0, -0.7d0, -0.65d0, -0.6d0, -0.55d0, &
      -0.5d0, -0.45d0, -0.4d0, -0.35d0, -0.3d0, &
      -0.25d0, -0.2d0, -0.15d0, -0.1d0, -0.04999999999999998d0, &
      0.d0, 0.04999999999999996d0, 0.1d0, 0.15d0, 0.2d0, &
      0.25d0, 0.2999999999999999d0, 0.35d0, 0.4d0, 0.45d0, &
      0.5d0, 0.55d0, 0.6d0, 0.65d0, 0.7d0, &
      0.75d0, 0.8d0, 0.85d0, 0.9d0, 0.95d0, &
      1.d0, 1.05d0, 1.1d0, 1.15d0, 1.2d0, &
      1.25d0, 1.3d0, 1.35d0, 1.4d0, 1.45d0, &
      1.5d0, 1.55d0, 1.6d0, 1.65d0, 1.7d0, &
      1.75d0, 1.8d0, 1.85d0, 1.9d0, 1.95d0, &
      2.d0, 2.05d0, 2.1d0, 2.15d0, 2.2d0, &
      2.25d0, 2.3d0, 2.35d0, 2.4d0, 2.45d0, &
      2.5d0, 2.55d0, 2.6d0, 2.65d0, 2.7d0, &
      2.75d0, 2.8d0, 2.85d0, 2.9d0, 2.95d0, &
      3.d0 /)

   real(dp), parameter :: ztab(Ntab) = (/ &
    -1.954534931647418d0, -1.904595955127455d0, -1.854672492899659d0, -1.804768625202741d0, -1.754889055579303d0, &
    -1.70503989763684d0, -1.655228654081395d0, -1.60546457725396d0, -1.555759054631801d0, -1.50612618706282d0, &
    -1.456583012930077d0, -1.407150381625654d0, -1.357853463490692d0, -1.308722546943003d0, -1.259793783848285d0, &
    -1.211110025410423d0, -1.162721552398376d0, -1.114686926137314d0, -1.067073489209121d0, -1.019957919118623d0, &
    -0.973426187114811d0, -0.927573629841689d0, -0.882504179826882d0, -0.838329240897075d0, -0.7951663065123026d0, &
    -0.7531368740344212d0, -0.712364209354681d0, -0.6729708569977019d0, -0.6350761725860893d0, -0.5987942772608892d0, &
    -0.5642323420331871d0, -0.5314896446427603d0, -0.500657505322406d0, -0.4718198814944222d0, -0.4450548063887798d0, &
    -0.4204362317928586d0, -0.3980359402924858d0, -0.3779255488513294d0, -0.3601775843640181d0, -0.3448662605931277d0, &
    -0.3320666821757335d0, -0.3218530224841464d0, -0.3142956220688198d0, -0.3094562610809645d0, -0.3073824903586883d0, &
    -0.3081007410368466d0, -0.311608744651744d0, -0.3178677839060045d0, -0.3267957945976536d0, -0.3382622116822907d0, &
    -0.352086308816749d0, -0.368039876518249d0, -0.3858558276600881d0, -0.4052421937217964d0, -0.4259005484112932d0, &
    -0.4475462037436712d0, -0.4699266031179398d0, -0.492834147878629d0, -0.516111317439216d0, -0.5396479260203524d0, &
    -0.5633729555416225d0, -0.587241999614359d0, -0.6112300779302227d0, -0.635315316521392d0, -0.6595026356985318d0, &
    -0.6837667225316001d0, -0.7081051552422673d0, -0.7325108440467765d0, -0.7569773284865642d0, -0.781499115184963d0, &
    -0.806069530077708d0, -0.830685025617995d0, -0.855340667784864d0, -0.880032389332806d0, -0.904756549300347d0, &
    -0.929509839357558d0, -0.954289397446105d0, -0.97909205891802d0, -1.003915870416773d0, -1.028758468070072d0, &
    -1.053617889340172d0, -1.078492349759026d0, -1.10338028736336d0, -1.128280261376662d0, -1.153190993098538d0, &
    -1.178111336852607d0, -1.203040266708068d0, -1.227976860843079d0, -1.252920305132837d0, -1.277869861081117d0, &
    -1.302824871927685d0, -1.327784750739423d0, -1.352748973137513d0, -1.377717070767324d0, -1.402688625346585d0, &
    -1.42766326348475d0, -1.452640651863115d0, -1.477620492999813d0, -1.502602521436235d0, -1.52758650032124d0, &
    -1.552572218353303d0 /)

    lx = log10(x)

    !--------------------------------------------------------
    ! left extrapolation
    !--------------------------------------------------------

    if (lx <= xtab(1)) then
        F_brem = 0d0
        return
    endif

    !--------------------------------------------------------
    ! right extrapolation
    !--------------------------------------------------------

    if (lx >= xtab(Ntab)) then
        F_brem = 0d0
        return
    endif

    !--------------------------------------------------------
    ! linear interpolation
    !--------------------------------------------------------

    do i = 1, Ntab-1
        if (lx >= xtab(i) .and. lx < xtab(i+1)) then

            t = (lx - xtab(i)) &
                / (xtab(i+1)-xtab(i))

            F_brem = 10d0**( &
                ztab(i) &
                + t*(ztab(i+1)-ztab(i)) )
            return
        endif
    enddo

    F_brem = 0d0

end function F_brem


   !---------------------------------------------------------
   ! dLdrdk function
   !---------------------------------------------------------
   function dLdrdk(r, k, mz, eps2, LT, T, ne, nN) result(res)
      ! (r, k,mz, T, ne, nN) all in units of eV^n. 
      ! return ddL/dr/dk in units of eV^2
      real(dp), intent(in) :: r, k, mz, eps2,T, ne, nN
      character(len=*), intent(in) :: LT

      real(dp) :: res
      real(dp) :: omega, fgamma, omega_p2
      real(dp) :: etaLT, RePi
      real(dp) :: Thom, Brem
      real(dp) :: Gamma_gamma
      real(dp) :: Gamma_z_gain

      real(dp), parameter :: me= 511000.  ! electron mass in eV

      omega = sqrt(mz**2 + k**2)
      omega_p2= (4.0_dp * pi * alpha_EM * ne) /me

      fgamma = 1.0_dp / (exp(omega / T) - 1.0_dp)

      etaLT = 1.0_dp
      if (LT == "L") then
         etaLT = (mz / omega)**2
      end if

      RePi = omega_p2 * etaLT

      Thom = (8.0_dp * pi * alpha_EM**2 * ne) / (3.0_dp * me**2)

      ! Brem= 0.0_dp
      Brem = 64.0_dp * pi**2 * alpha_EM**3 * ne * nN / &
            (3.0_dp * sqrt(2.0_dp * pi * T * me**3) * omega**3) * &
            F_brem(omega / T)

      Gamma_gamma = etaLT / (1.0_dp + fgamma) * (Thom + Brem)

      Gamma_z_gain = fgamma * Gamma_gamma * eps2 * mz**4 / &
                     ((mz**2 - RePi)**2 + (omega * Gamma_gamma)**2)

      res = 4.0_dp * pi * r**2 * k**2 * omega / (2.0_dp * pi**2) * Gamma_z_gain ! Eq.(eq:-31) in the note

   end function dLdrdk


   !---------------------------------------------------------
   ! dLdr function, 
   ! returns dL/dr = ∫ dk d²L/(dr dk)
   !---------------------------------------------------------

   function dLdr(r, mz, eps2, LT, T, ne, nN) result(res)
      ! (r,mz, T, ne, nN) all in units of eV^n. 
      ! return dL/dr in units of eV^3
      implicit none
      integer, parameter :: npts = 100

      real(dp), intent(in) :: r, mz, eps2,T, ne, nN

      character(len=*), intent(in) :: LT

      real(dp) :: res

      real(dp) :: kmin, kmax
      real(dp) :: integral

      real(dp), dimension(npts) :: k_grid
      real(dp), dimension(npts) :: y

      integer :: i

      !-----------------------------------------
      ! logarithmic momentum grid
      !-----------------------------------------

      kmin = 1.0e-3_dp * T
      kmax = 20.0_dp   * T

      do i = 1, npts
         k_grid(i) = kmin * (kmax/kmin)**( real(i-1,dp) / real(npts-1,dp) )

         y(i) = dLdrdk(r, k_grid(i), mz, eps2, LT, T, ne, nN)
      end do

      !-----------------------------------------
      ! trapezoidal integration
      !-----------------------------------------

      integral = 0.0_dp
      do i = 1, npts-1
         integral = integral + 0.5_dp * (y(i) + y(i+1)) * (k_grid(i+1) - k_grid(i))
      end do

      res = integral

   end function dLdr



!============================================================
! Z interpolation
!
! Input:
!   x
!
! Uses:
!   lx = log10(x)
!
! Table stores:
!   log10(x), log10(Z)
!============================================================

real(dp) function Z_interp(x)

    implicit none

    real(dp), intent(in) :: x

    integer :: i
    integer, parameter :: Ntab = 101

    real(dp) :: lx
    real(dp) :: t

real(dp), parameter :: xtab(Ntab) = (/ &
    -2.d0, -1.96d0, -1.92d0, -1.88d0, -1.84d0, &
    -1.8d0, -1.76d0, -1.72d0, -1.68d0, -1.64d0, &
    -1.6d0, -1.56d0, -1.52d0, -1.48d0, -1.44d0, &
    -1.4d0, -1.36d0, -1.32d0, -1.28d0, -1.24d0, &
    -1.2d0, -1.16d0, -1.12d0, -1.08d0, -1.04d0, &
    -1.d0, -0.96d0, -0.92d0, -0.88d0, -0.84d0, &
    -0.8d0, -0.76d0, -0.7199999999999999d0, -0.6800000000000001d0, -0.64d0, &
    -0.6d0, -0.56d0, -0.52d0, -0.48d0, -0.44d0, &
    -0.4d0, -0.36d0, -0.32d0, -0.28d0, -0.24d0, &
    -0.2d0, -0.16d0, -0.12d0, -0.08d0, -0.04d0, &
    0.d0, 0.04d0, 0.08d0, 0.12d0, 0.16d0, &
    0.2d0, 0.24d0, 0.28d0, 0.32d0, 0.36d0, &
    0.4d0, 0.44d0, 0.48d0, 0.52d0, 0.56d0, &
    0.6d0, 0.64d0, 0.6800000000000001d0, 0.7199999999999999d0, 0.76d0, &
    0.8d0, 0.84d0, 0.88d0, 0.92d0, 0.96d0, &
    1.d0, 1.04d0, 1.08d0, 1.12d0, 1.16d0, &
    1.2d0, 1.24d0, 1.28d0, 1.32d0, 1.36d0, &
    1.4d0, 1.44d0, 1.48d0, 1.52d0, 1.56d0, &
    1.6d0, 1.64d0, 1.68d0, 1.72d0, 1.76d0, &
    1.8d0, 1.84d0, 1.88d0, 1.92d0, 1.96d0, &
    2.d0 /)

real(dp), parameter :: ztab(Ntab) = (/ &
    0.3809026470831913d0, 0.3808930528263141d0, 0.3808817200621327d0, 0.3808683379693417d0, 0.3808525415002988d0, &
    0.3808339007501123d0, 0.3808119108439298d0, 0.3807859797178787d0, 0.3807554126718882d0, 0.3807193949984809d0, &
    0.380676971827529d0, 0.3806270247509763d0, 0.3805682447172253d0, 0.380499104743902d0, 0.3804178028747405d0, &
    0.3803222614090338d0, 0.3802100369832194d0, 0.3800782851777853d0, 0.3799236918421286d0, 0.3797423988760049d0, &
    0.3795299190104304d0, 0.3792810381068002d0, 0.3789897033237403d0, 0.3786488953519528d0, 0.3782504825018759d0, &
    0.3777850548520306d0, 0.3772417354042401d0, 0.376607965909633d0, 0.3758692642790625d0, 0.375008950301638d0, &
    0.3740078364793079d0, 0.3728438773728684d0, 0.3714917812791216d0, 0.3699225660804315d0, 0.3681030676314192d0, &
    0.3659953869675507d0, 0.3635562728352177d0, 0.3607364323007308d0, 0.3574797617887859d0, 0.3537224901065419d0, &
    0.3493922237372724d0, 0.3444068837525519d0, 0.3386735220116882d0, 0.3320870028043609d0, 0.3245285341506739d0, &
    0.3158640308666323d0, 0.3059422891687142d0, 0.2945929500997985d0, 0.2816242265015259d0, 0.266820365782428d0, &
    0.2499388185540407d0, 0.230707081629453d0, 0.2088191832962618d0, 0.1839317796825246d0, 0.1556598340073488d0, &
    0.1235718561625921d0, 0.087184688969809d0, 0.04595783994343028d0, -0.0007126266456097722d0, -0.05350060407597008d0, &
    -0.1131568114967096d0, -0.180515948006477d0, -0.2565040567398507d0, -0.3421460307315099d0, -0.4385732216838106d0, &
    -0.5470311688809351d0, -0.6688875485230258d0, -0.805640546174573d0, -0.958927960213106d0, -1.130537428256076d0, &
    -1.32241820651969d0, -1.53669490720223d0, -1.775683512043281d0, -2.04190985327849d0, -2.338130624633891d0, &
    -2.667356894375627d0, -3.032880070914779d0, -3.43830031263896d0, -3.887557469550477d0, -4.384964750098785d0, &
    -4.935245400672278d0, -5.543572750472626d0, -6.215614011219197d0, -6.957578245210577d0, -7.776268933675123d0, &
    -8.67914160388272d0, -9.67436700791803d0, -10.77090038918347d0, -11.97855742295805d0, -13.30809747372275d0, &
    -14.77131487420077d0, -16.38113899938202d0, -18.15174398366512d0, -20.09866901134336d0, -22.23895020055355d0, &
    -24.59126519935601d0, -27.17609172087074d0, -30.0158813621527d0, -33.1352501828648d0, -36.56118765968811d0, &
    -40.32328579123514d0 /)
    !--------------------------------------------------------
    ! safety
    !--------------------------------------------------------

    if (x <= 0d0) then
        Z_interp = 0d0
        return
    endif

    lx = log10(x)

    !--------------------------------------------------------
    ! left extrapolation
    !--------------------------------------------------------

    if (lx <= xtab(1)) then
        Z_interp = 10d0**ztab(1)
        return
    endif

    !--------------------------------------------------------
    ! right extrapolation
    !--------------------------------------------------------

    if (lx >= xtab(Ntab)) then
        Z_interp = 0d0
        return
    endif

    !--------------------------------------------------------
    ! linear interpolation
    !--------------------------------------------------------

    do i = 1, Ntab-1

        if (lx >= xtab(i) .and. lx < xtab(i+1)) then

            t = (lx - xtab(i)) &
                / (xtab(i+1)-xtab(i))

            Z_interp = 10d0**( &
                ztab(i) &
                + t*(ztab(i+1)-ztab(i)) )

            return

        endif

    enddo

    Z_interp = 0d0

end function Z_interp

   subroutine darkphoton_cooling(id, ierr)

    use star_def
    use const_def
    use chem_def
    use auto_diff

    integer, intent(in) :: id
    integer, intent(out) :: ierr

    type (star_info), pointer :: s

    integer :: k

    real(dp) :: eps_mix, eps_sign,eps2


    real(dp) :: rho_cgs
    real(dp) :: Ye

    real(dp) :: T_cgs,ne_cgs, ne_cgs_low,dne_cgs,nN_cgs, m_DP_eV, x, x2,x_low, y,y2 ,r_cgs,dr_cgs
    real(dp) :: bose,Z_y
    real(dp) :: L_nuc,L_DP_L,L_DP_T,L_neu,L_grav, dL_L,dL_T
    real(dp) :: dr

    

    ierr = 0

    call star_ptr(id, s, ierr)
    if (ierr /= 0) return

    ! user controls
    eps_sign = s% x_ctrl(1)
    eps_mix = s% x_ctrl(2)
    eps2=eps_mix*eps_mix
    m_DP_eV  = s% x_ctrl(3)
    write(*,*) 'dark photon mass [eV] and mixing:', m_DP_eV, eps_mix

   if (abs(eps_mix) < 1d-99) then
      write(*,*) 'eps_mix≈0, skip!'
      current_L_DP_L = 0.0_dp
      current_L_DP_T = 0.0_dp
      return
   end if


    L_DP_L = 0d0 
    L_DP_T = 0d0
    L_nuc= 0d0 
    L_neu= 0d0 
    L_grav= 0d0 

    do k = 2, s% nz

        rho_cgs = s% rho(k)
        T_cgs   = s% T(k)
        r_cgs = s% rmid(k)
        dr_cgs =   s% dm(k) / (4.0_dp * pi * (s% rmid(k)**2) * s% rho(k))


        ! electron fraction
        Ye = s% ye(k)

        ! electron number density
        ne_cgs = rho_cgs * Ye / amu ! amu already defined in const_def, it is approximately 1.66054×10−24 gram

        nN_cgs = rho_cgs / amu !the effective number of \Sum (Z^2 n_Z) used in Bremsstrahlung
         ! Strickly speaking, n_p =  X * ρ / mu; n_He = Y * ρ / (4 * m_u) 
         !n_N =n_p+ 4* n_He  which is \Sum (Z^2 n_Z). Since n_N =n_p+ 4* n_He=(X+Y)* ρ / mu and 
        ! X+Y \approx 1,  we have n_N =(X+Y)* ρ / mu


        ne_cgs_low = s% rho(k-1) * s% ye(k-1)   / amu
        dne_cgs=ne_cgs - ne_cgs_low

        ! x=ωp/T, which should be dimensionless. Here plasma frequency and T are in natural units.
        ! y= mz /T, which should be dimensionless. Also in natural units.
        x= 5.043465064609606d-6 *  sqrt(alpha_EM*ne_cgs)/T_cgs ! the numerical factor is obtained in Mathematica
        x_low=5.043465064609606d-6 *  sqrt(alpha_EM*ne_cgs_low)/T_cgs
        y= 11604.96692584426d0  *  m_DP_eV/T_cgs ! the numerical factor is obtained in Mathematica
        x2=x*x 
        y2=y*y
      
         ! L and T mode contribution to be added to s% extra_heat(k); dL_T doesn't include the factor of two due to two T modes. Hence  (dL_L + 2*dL_T) will be used.
         dL_L = 0d0 ! Unit:  erg/s
         dL_T = 0d0 ! Unit:  erg/s


      !==========================================
      ! LONGITUDINAL
      !==========================================
        if (x > y) then !L resonance occurs; Using expressions generated in 20260528.nb
         dL_L=(0.0015056275978405718*r_cgs**2*T_cgs**5*x2*Sqrt(x2 - y2)*y2* dr_cgs *eps2)/(exp(x) - 1d0)      
        else ! brute-force integration
         dr= dr_cgs *cm_over_eVinv 
 
         dL_L = dr * dLdr(&
         r_cgs *cm_over_eVinv, &!r=r_cgs *cm_over_eVinv 
         m_DP_eV, &! mz in eV
         eps2, & ! eps2
         "L", & ! LT
         T_cgs * Kelvin_over_eV, & !T
         ne_cgs * cm_inv3_over_eV3, & !ne
         nN_cgs * cm_inv3_over_eV3 & ! nN
         ) * eV2_over_erg_s
        end if

      !==========================================
      ! TRANSVERSE
      !==========================================
      !   x ~ ωp, x_low~ωp_low (obtained with k-1, the layer more outside, hence x_low<x )
        ! x is in the inner layer, so x should be higher than x_low
        if (x > y .and. y>x_low) then !T resonance occurs; Using expressions generated in 20260528.nb
         Z_y=Z_interp(y)
         dL_T=(1.1838305101633038e8*r_cgs**2*T_cgs**7*y2**2*Z_y*dr_cgs*eps2)/(alpha_EM*dne_cgs)
        else ! brute-force integration
         dr= dr_cgs *cm_over_eVinv 
 
         dL_T = dr * dLdr(&
         r_cgs *cm_over_eVinv, &!r=r_cgs *cm_over_eVinv 
         m_DP_eV, &! mz in eV
         eps2, & ! eps2
         "T", & ! LT
         T_cgs * Kelvin_over_eV, & !T
         ne_cgs * cm_inv3_over_eV3, & !ne
         nN_cgs * cm_inv3_over_eV3 & ! nN
         ) * eV2_over_erg_s

        end if
      !==========================================
      ! Add L and 2T to extra_heat(k)
      !==========================================

      ! cooling rate. Unit:  erg/g/s
      s% extra_heat(k) = eps_sign  * (dL_L + 2*dL_T) / (s% dm(k))


      L_DP_L = L_DP_L + dL_L  
      L_DP_T = L_DP_T + dL_T 
      L_nuc = L_nuc + s% eps_nuc(k) * s% dm(k)
      L_neu   = L_neu   + s% non_nuc_neu(k)       * s% dm(k)

    end do
   ! L_standard = s% L(1)
   ! L_total = L_standard + L_new
    current_L_DP_L=L_DP_L
    current_L_DP_T=L_DP_T
   write(*,'(A,1PE12.4)') 'surface L  [Lsun] = ', s% L(1) /Lsun
   write(*,'(A,1PE12.4)') 'L_DP_L      [Lsun] = ', L_DP_L/Lsun
   write(*,'(A,1PE12.4)') 'L_DP_T      [Lsun] = ', L_DP_T/Lsun
   write(*,'(A,1PE12.4)') 'L_nuc      [Lsun] = ', L_nuc /Lsun
   write(*,'(A,1PE12.4)') 'L_neu    [Lsun] = ',   L_neu /Lsun
   write(*,*) '(L_nuc_burn_total , L_phot)=', (s% L_nuc_burn_total)  ,s% L_phot
   ! write(*,'(A,1PE12.4)') 'L_grav    [erg/s] = ', L_grav
   end subroutine darkphoton_cooling

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

      s% other_energy => darkphoton_cooling
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
      ! condition was triggered.  MESA provides 9 customizable
      ! termination codes, named t_xtra1 .. t_xtra9.  You can
      ! customize the messages that will be printed upon exit by
      ! setting the corresponding termination_code_str value.
      ! termination_code_str(t_xtra1) = 'my termination condition'

      ! by default, indicate where (in the code) MESA terminated
      if (extras_check_model == terminate) s% termination_code = t_extras_check_model
   end function extras_check_model


   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_columns = 3
   end function how_many_extra_history_columns


   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      integer:: k
      character (len=maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n), L_new
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! note: do NOT add the extras names to history_columns.list
      ! the history_columns.list is only for the built-in history column options.
      ! it must not include the new column names you are adding here.

      names(1) = 'L_extra'
      L_new = 0d0 
      do k = 1, s% nz
         L_new = L_new + s% extra_heat(k)%val * s% dm(k)
      end do
      vals(1) = L_new/Lsun 

      names(2) = 'L_extra_L'
      vals(2) = current_L_DP_L/Lsun 

      names(3) = 'L_extra_T'
      vals(3) = current_L_DP_T/Lsun 

   end subroutine data_for_extra_history_columns


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


end module run_star_extras
