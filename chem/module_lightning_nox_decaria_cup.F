!WRF:MODEL_LAYER:CHEMISTRY
!
! Contains subroutine for converting flash rates into NO emissions
! based on Decaria 2000 vertical distirbutions.
!
! Input: flashes (#/s)
! Output: tendency (ppmv/s)
!
! The output will be muliplied by timestep and used to incremeent NO
! concentration and the respective passive tracer in lightning_nox_driver.
!
! See module_lightning_nox_driver for more info.
!
! Contact: M. Barth <barthm@ucar.edu>, J. Wong <johnwong@ucar.edu>
!
!lm 2015/09: adapted to be used for KFCuP + changed from tile average to local
!            This routine will only work for KFCuP, and should only be called
!            when KFCuP is on and the CPM option for lightning is turned on
!
!**********************************************************************
 MODULE module_lightning_nox_decaria_cup

 IMPLICIT NONE

 CONTAINS

!**********************************************************************
!
! DeCaria et al, 2000
!
! DeCaria, A. J., K. E. Pickering, G. L. Stenchikov, and L. E. Ott (2005),
! Lightning-generated NOX and its impact on tropospheric ozone production:
! A three-dimensional modeling study of a Stratosphere-Troposphere Experiment:
! Radiation, Aerosols and Ozone (STERAO-A) thunderstorm, J. Geophys. Res.,
! 110, D14303, doi:10.1029/2004JD005556.
!
!**********************************************************************
 SUBROUTINE lightning_nox_decaria_cup ( &
                          ! Frequently used prognostics
                            dx, dy, xland, ht, t_phy, rho, z, p,  &
                            ic_flashrate, cg_flashrate,           & ! flashes (#/s)
                          ! Scheme specific prognostics
!lm
                            cldfra_cup, htop,                     &
!lm
                          ! Namelist inputs
                            N_IC, N_CG,                           &
                            ltng_temp_upper,ltng_temp_lower,      &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! outputs
                            lnox_ic_tend, lnox_cg_tend            & ! tendency (ppmv/s)
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_max_real, wrf_dm_min_real, wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht, htop
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: t_phy, rho, z, p
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: cldfra_cup
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: ic_flashrate  , cg_flashrate ! #/sec


! Scheme specific namelist inputs
 REAL,    INTENT(IN   )    ::       N_IC, N_CG
 REAL,    INTENT(IN   )    ::       ltng_temp_lower, ltng_temp_upper

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(  OUT) :: lnox_ic_tend,lnox_cg_tend

! Local variables
 INTEGER :: i,k,j
 INTEGER :: ktop,kbtm,kupper,klower
 REAL :: ic_fr, cg_fr, delta ! 1d flashrates
 CHARACTER (LEN=250) :: message
 REAL, DIMENSION( kps:kpe ) :: molesofair, cldfra_cup_1d, t_1d, z_1d
 REAL :: htop_ij
 REAL, DIMENSION( kps:kpe ) :: fd, fd2, dz ! fd = distribution

!-----------------------------------------------------------------

 lnox_ic_tend (ips:ipe,kps:kpe,jps:jpe ) = 0.
 lnox_cg_tend (ips:ipe,kps:kpe,jps:jpe ) = 0.

 DO i=ips,ipe
   DO j=jps,jpe

   IF ((ic_flashrate(i,j) .gt. 0) .OR. (cg_flashrate(i,j) .gt. 0)) THEN
   
   !1D or scalar values at i,j   
   cldfra_cup_1d(kps:kpe) = cldfra_cup(i,kps:kpe,j)
   htop_ij = htop(i,j)
   t_1d(kps:kpe) = t_phy(i,kps:kpe,j)
   z_1d(kps:kpe) = z(i,kps:kpe,j)
   ic_fr = ic_flashrate(i,j)
   cg_fr = cg_flashrate(i,j)
  
   !kfind for kfcup. Find cumulus top, bottom and isotherms 
   CALL  kfind ( cldfra_cup_1d, htop_ij, t_1d,    &
                 ltng_temp_upper,ltng_temp_lower, &
                 i, j,                            &
                 ips, ipe, jps, jpe, kps, kpe,    &
               ! Outputs
                 ktop,kbtm,kupper,klower )

   !lm moles of air used for mixing ratio. This is not exactly right (map factors). 
   !Find a better way to do this.
   molesofair(kps:kpe) = rho(i,kps:kpe,j) * 1E3 * dx * dy / .02897 ! # moles per km in z
 
   ! Calculates IC distribution
   IF (( ic_fr .gt. 0 ) .and. (( ktop .gt. klower ) .and. (kbtm .lt. ktop) ) )THEN
 
     call bellcurve(kbtm,ktop,klower,z_1d, kps,kpe, fd, dz)
     if (ktop .gt. kupper) then
       call bellcurve(kbtm,ktop,kupper,z_1d, kps,kpe, fd2, dz)
       fd(kbtm:ktop) = 0.5*( fd(kbtm:ktop) + fd2(kbtm:ktop) )
     endif
   
     DO k=kbtm,ktop
       IF ( cldfra_cup_1d(k) .gt. 0.01 ) THEN
         !*  implementation note: 1) ic_fr * N_IC gives moles of NO in the column !lm check this, fd unit ?
         !*                       2) Multiplying by fd gives the # moles per level
         !*                       3) Divide by moles of air gives mixing ratio
         !*                       4) Multiply by 1E6 gives ppmv
         delta = (ic_fr * N_IC) * fd(k) / (molesofair(k)*dz(k)) * 1E6
         WRITE(message, * ) ' LNOx_driver: k, delta, fd = ', k, delta, fd(k)
         CALL wrf_debug ( 100, message )
         lnox_ic_tend(i,k,j) = delta 
       ENDIF
     ENDDO

  ENDIF ! IC lightning

  !-----------------------------------------------------------------
  ! Calculates CG distribution
  IF ((cg_fr .gt. 0 ) .and. (( ktop .gt. klower ) .and. (kbtm .lt. ktop) ) ) THEN
    call bellcurve(kps,ktop,klower,z_1d, kps,kpe, fd, dz)
 
    delta = (cg_fr * N_CG ) * fd(k) / (molesofair(k)*dz(k)) * 1E6
    ! below cloud
    !*  implementation note: previously, between kps & kbtm is practically excluded because of
    !*                       the test for refl > 20
    IF ( cldfra_cup_1d(kbtm) .gt. 0.01 ) THEN
      lnox_cg_tend(i,kps:(kbtm-1),j) = delta
    ENDIF

    ! within cloud
    do k = kbtm,ktop
      delta = (cg_fr * N_CG ) * fd(k) / (molesofair(k)*dz(k)) * 1E6
      lnox_cg_tend(i,k,j) = delta
    enddo
  ENDIF
  ENDIF ! ic_fr(i,j)>0
  ENDDO !jps 
ENDDO !ips

 END SUBROUTINE lightning_nox_decaria_cup



!************************************************************************
! This subroutine prepares a normal distribution between k_min and
! k_max centered at k_mu. Distribution for each level is
! normalized to \int^{z_at_w(k_max}_{z_at_w(k_min-1)}f(z)dz
!
! Unit of f is fraction of total column
!
! Modified from v3.4.1 module_ltng_crm.F, kept the math but changed
! the implementation for better clarity. Removed patch-wide averaging
! of z.
!************************************************************************

 SUBROUTINE bellcurve ( k_min, k_max, k_mu, z, kps,kpe, f, dz )
!-----------------------------------------------------------------

 IMPLICIT NONE
 INTEGER,                      INTENT(IN   ) :: k_min, k_max, k_mu
 REAL,   DIMENSION( kps:kpe ), INTENT(IN   ) :: z       ! at phy
 INTEGER,                      INTENT(IN   ) :: kps,kpe

 REAL,   DIMENSION( kps:kpe ), INTENT(  OUT) :: f, dz

 INTEGER :: i,j,k
 REAL, DIMENSION( kps:kpe ) :: ex
 REAL :: sigma, z_mu, cuml_f_dist
 REAL, PARAMETER :: two_pi = 6.2831854

!-----------------------------------------------------------------

 f(kps:kpe) = 0.
 z_mu = z(k_mu)
 sigma = AMIN1(z(k_max)-z_mu,z_mu-z(k_min))/3.0

 ! distance from mean
 ex(k_min:k_max) = (z(k_min:k_max)-z_mu)/sigma
 
 ! Truncated Gaussian at 3 sigma
 f(k_min:k_max) = (1.0/(sqrt(two_pi)*sigma))*exp(-ex(k_min:k_max)*ex(k_min:k_max)/2.0)

!++mcb   We do need dz at bottom and top of domain
!   dz(kps) = 0. ! safe as long as k_min != kps
!   dz(kpe) = 0. ! safe as long as k_max != kpe
 dz(kps) = (z(kps+1) - z(kps))*.5             
 dz(kpe) = (z(kpe) - z(kpe-1))*.5             
 DO k=kps+1,kpe-1
!  dz(k) = (z(k+1)+z(k))/2. - (z(k)+z(k-1))/2.
   dz(k) = (z(k+1) - z(k-1))*.5
 ENDDO

 ! Normalize
 cuml_f_dist = sum(f(k_min:k_max) * dz(k_min:k_max))
 f(k_min:k_max) = f(k_min:k_max)*dz(k_min:k_max)/cuml_f_dist

 END SUBROUTINE bellcurve

!************************************************************************
! This subroutine finds the vertical indices (phy grid) of the follow
! within a column:
! 1) ktop - cloud top
! 2) kbtm - cloud bottom
! 3) kupper - upper mode isotherm
! 3) klower - lower mode isotherm
!************************************************************************
 SUBROUTINE kfind ( &
              ! Prognostics
                cldfra_cup_1d, htop_ij, t_1d,         &
              ! Namelist settings
                ltng_temp_upper,ltng_temp_lower,      &
                i , j,                                &
              ! Order dependent args for domain, mem, and tile dims
                ips, ipe, jps, jpe, kps, kpe,          &
              ! Outputs
                ktop,kbtm,kupper,klower               &
            )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants

 USE module_dm, only: wrf_dm_max_real, wrf_dm_min_real, wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Prognostics
 REAL, DIMENSION( kps:kpe ), INTENT(IN   ) :: cldfra_cup_1d
 REAL, DIMENSION( kps:kpe ), INTENT(IN   ) :: t_1d
 REAL, INTENT(IN   )    ::     htop_ij

! Namelist settings
 REAL,    INTENT(IN   )    ::       ltng_temp_lower, ltng_temp_upper

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

 INTEGER, INTENT(IN   )    ::       i, j

! Outputs
 INTEGER, INTENT(  OUT)    ::       ktop,kbtm,kupper,klower

! Local vars
 CHARACTER (LEN=250) :: message
 REAL    :: ktop_r, kbtm_r, kupper_r, klower_r
 INTEGER :: k

!-----------------------------------------------------------------
 ktop = kps
 kbtm = kps
 kupper = kps
 klower = kps

 ! look for ktop
 k = kpe
 DO WHILE ( cldfra_cup_1d(k) .le. 0.01 .and. k .gt. kps)
  k = k-1
 ENDDO
 ktop = k

 ! Look for kbtm
 k = kps
 DO WHILE( cldfra_cup_1d(k) .le. 0.01 .and. k .le. ktop )
  k = k+1
 ENDDO
 kbtm = k
! if (kbtm .eq. kps) kbtm = kpe

 ! Look for kupper
 k = kps
 DO WHILE ( t_1d(k) .gt. ltng_temp_lower + 273.15 .and. k .lt. kpe ) !lm for wrfdev
   k = k + 1
 ENDDO
 klower = k

 DO WHILE ( t_1d(k) .gt. ltng_temp_upper + 273.15 .and. k .lt. kpe ) !lm for wrfdev
   k = k + 1
 ENDDO
 kupper = k

 !for now, overwrite ktop with htop. For some reason, cldfra can be < 0 below
 !the KFCuP cloud top
 ktop = nint(htop_ij)

 WRITE(message, * ) ' LNOx_driver: kbtm, ktop, klower, kupper = ', kbtm, ktop, klower, kupper
 CALL wrf_debug ( 100, message )
 
 END SUBROUTINE kfind

 END MODULE module_lightning_nox_decaria_cup
