module aero_sectional
!!! Sectional aerosol code based on
!!! Blichner, S. M., Sporre, M. K., Makkonen, R., and Berntsen, T. K.:
!!! Implementing a sectional scheme for early aerosol growth from new particle
!!! formation in the Norwegian Earth System Model v2: comparison to
!!! observations and climate impacts, Geosci. Model Dev., 14, 3335–3359
!!! https://doi.org/10.5194/gmd-14-3335-2021, 2021

   use shr_kind_mod, only: r8 => shr_kind_r8
   use chem_mods,    only: gas_pcnst
   use aerosoldef,   only: l_soa_a1, l_so4_a1

   implicit none
   private

   ! Public interface
   public :: aerosect_register
   public :: aerosect_write2file
   public :: sec_numberConc
   public :: sec_moveMass

   ! Public data
   logical, public, protected :: do_sectional_NPF = .false.

   integer, public, parameter :: secNrBins   = 5 ! nr of bins
   integer, public, parameter :: secNrSpec   = 2 ! number of condensing species

   real(r8), public, parameter :: rhopart_sec(secNrSpec) = (/ 1769.0_r8,1500.0_r8 /) ! same as SO4_NA, SOA_NA

   character(len=*), public, parameter :: secSpecNames(secNrSpec) = (/ 'SO4_SEC','SOA_SEC'/) ! names of condensed species
   character(len=*), public, parameter :: SpecNames(secNrSpec)   = (/ 'SO4','SOA'/)          ! names of condensing species

   ! coagulation receiver
   integer, public, protected  :: secCoagulate_receiver(secNrSpec) ! (/ l_so4_a1,l_soa_a1/)
   ! array holding mean volume in each bin
   real(r8), public, protected :: secMeanVol(secNrBins)            ! [m3]
   ! holds the chemistry indices of the sectional scheme:
   integer, public protected   :: secConstIndex(secNrSpec, secNrBins)

   ! Private data
   ! See Blichner et.al., sec. 2.2
   real(r8), parameter, public :: sec_min_diameter = 5.0e-9_r8  ! [m] minumum diameter
   real(r8), parameter :: sec_max_diameter = 39.6e-9_r8 ! [m] volume median diameter of NPF background mode
                                                        ! calculated so that volume2number should be correct
   real(r8), parameter :: max_volume = sec_max_diameter**3 * pi / 6.0_r8

CONTAINS

   subroutine aerosect_register()
      ! Sets up the chemistry indices for the tracers in the sectional scheme
      use constituents, only: cnst_get_ind

      integer              :: secInd,volInd
      character(len=20)    :: cnst_name
      do secInd = 1, secNrBins
         do volInd = 1, secNrSpec ! Names constituents as 'SOA_SEC01' & 'SO4_SEC01'
            write(cnst_name,'(A,I2.2)') trim(secSpecNames(volInd)), secInd
            call cnst_get_ind(trim(cnst_name), secConstIndex(volInd,secInd), abort=.true.)
         end do
      end do

      do_sectional_NPF = .true.

   end subroutine aerosect_register

   subroutine aerosect_init()
      ! Sets up volume of bins.
      real(r8) :: d_rat
      real(r8) :: secMeanD(secNrBins) ! [m]
      integer  :: ind                 ! index

      secCoagulate_receiver = (/ l_so4_a1, l_soa_a1/)  ! number of condensing species

      ! Use discrite geometric distribution/volume-ratio size distrib:
      d_rat = (sec_max_diameter/sec_min_diameter)**(1._r8/secNrBins)
      secMeanD(1) = sec_min_diameter
      do ind = 2, secNrBins
         secMeanD(ind) = secMeanD(ind-1) * d_rat
      end do
      do ind = 1, secNrBins
         secMeanVol(ind) = secMeanD(ind)**3 * pi / 6._r8
      end do

   end subroutine aerosect_init

   subroutine aerosect_write2file(q, lchnk,ncol, pmid, temperature)
      ! Routine writes number concentration of aerosols to history
      use cam_history, only: outfld
      use ppgrid,      only: pcols, pver
      use aerosoldef,  only: chemistryIndex
      use physconst,   only: rair

      ! Dummy arguments
      integer,  intent(in) :: lchnk                    ! chunk identifier
      integer,  intent(in) :: ncol                     ! number of columns
      real(r8), intent(in) :: q(pcols,pver,gas_pcnst)  ! tmr [kg/kg]
      real(r8), intent(in) :: pmid(pcols,pver)         !
      real(r8), intent(in) :: temperature(pcols,pver)  ! [K] Temperature

      ! Local variables:
      character(len=20)    :: field_name
      real(r8)             :: rhoAir                   ! [kg/m3] density of air
      integer              :: indBin, indSpec, ind_sec ! indices
      integer              :: icol, ilev               ! indices
      real(r8)             :: num_conc(pcols,pver)     ! [#/m3] number concentration

      do indBin = 1, secNrBins
         !Go through all core species in that bin
         do indSpec = 1, secNrSpec
            ind_sec = chemistryIndex(secConstIndex(indSpec,indBin))
            do ilev = 1, pver
               do icol = 1, ncol
                  rhoAir = pmid(icol,ilev) / (rair * temperature(icol,ilev))
                  call sec_numberConc(q(icol,ilev,ind_sec),indSpec, indBin, rhoAir, num_conc(icol,ilev))
               end do
            end do

            WRITE(field_name,'(A,A,I2.2)') 'nr',trim(secSpecNames(indSpec)),indBin
            call outfld(trim(field_name),num_conc, pcols, lchnk) !#
         end do
      end do

   end subroutine aerosect_write2file

   !XXG: Change this be take vectors or possibly be elemental?
   subroutine sec_numberConc(mass, volNr, binNr , rhoAir, numberConc)
      use physconst, only: pi
      ! Calculates the number concentration from the mass concentration

      ! Dummy arguments
      real(r8), intent(in)  :: mass       ! kg/kg
      integer,  intent(in)  :: binNr      ! bin_index
      real(r8), intent(in)  :: rhoAir
      real(r8), intent(out) :: numberConc ! #/m3_air
      integer,  intent(in)  :: volNr
      ! Local variables:
      integer               :: volInd

      numberConc = (mass / rhopart_sec(volNr)) * (rhoAir / secMeanVol(binNr))
      ![kg_aer/kg_air]/[kg_aer/m3_aer]*[kg_air/m3_air]/[m3_aer/#]--> #/m3_air
      !XXG: Where does this limit come from?
      if (mass < 1.e-35) then
         numberConc = 0.0_r8
      end if

   end subroutine sec_numberConc

   subroutine sec_moveMass(massDistrib, numberConc_old, leave_sec, rhoAir, modeDiam, decrease_dt)
      ! Moves tracer mass from one bin to the another based on condensational/coagulation growth.
      ! Based on Jacobson Fundamentals of Atmospheric Modeling, second edition (2005),
      ! Chapter   13.5
      use physconst,  only: pi
      use aerosoldef, only: chemistryIndex

      real(r8), intent(inout) :: massDistrib(gas_pcnst)    ! mass in each tracer
      real(r8), intent(in)    :: numberConc_old(secNrBins) ! numbr concentration before growth
      real(r8), intent(out)   :: leave_sec(secNrSpec)      ! the mass that leaves sectional scheme
      logical,  intent(out)   :: decrease_dt               ! if set to True, time step is divided
                                                           ! and the procedure is re run
      real(r8) :: rhoAir                               ! Density of air
      real(r8) :: modeDiam                             ! not used.

      real(r8) :: numberConc_new(secNrSpec, secNrBins) ! number concentration after growth
      real(r8) :: volume(secNrBins)                    ! volume of particle in bin
      real(r8) :: volume_new(secNrBins)                ! volume after growth
      integer  :: indBin,indSpec
      real(r8) :: xfrac                                ! fraction to stay in bin
      real(r8) :: volfrac(secNrSpec,secNrBins)         ! volume fraction
      real(r8) :: sumvolfrac                           ! temp


      decrease_dt = .FALSE.
      ! Compute volume in each bin with condensation (by mass) and by
      ! numberconcentration
      do indBin = 1, secNrBins
         volume_new(indBin) = 0.0_r8
         volfrac(:,indBin) = 0.0_r8
         do indSpec = 1, secNrSpec ! calculate volume in each bin by mass/density
            ! XXG: Why this particular number? Parameter?
            if (numberConc_old(indBin) < 1.e-30_r8) then
               volume_new(indBin) = 0.0_r8
            else
               volume_new(indBin) = volume_new(indBin) +                          &
                    (massDistrib(chemistryIndex(secConstIndex(indSpec,indBin))) / &
                    rhopart_sec(indSpec) * rhoAir / (numberConc_old(indBin)))
               !XXG: More efficient to add volume_new > max_volume check here?
            end if

            volfrac(indSpec,indBin) = massDistrib(chemistryIndex(secConstIndex(indSpec,indBin))) / &
                 rhopart_sec(indSpec)*rhoAir
            !kg/kg(air)*[kg(air)/m3(air)][kg/m3]--> m3/m3(air)
         end do ! calculate volume in each bin by numberconcentration (volume from before condenstion)
         !XXG: Roughtly equivalent?
         sumvolfrac = sum(volfrac(:,indBin))
         if (sumvolfrac < 1.e-50_r8) then
            volfrac(:,indBin) = 0.0_r8
         else
            ! XXG: remove arbitrary 1.E-50_r8?
            volfrac(:,indBin) = volfrac(:,indBin) / (sumvolfrac + 1.E-50_r8)
         end if
         ! calculate volume in each bin by mass/density! m3
         volume(indBin) = secMeanVol(indBin)
         ! calculate volume in each bin by numberconcentration (volume from before condenstion)
      end do
      numberConc_new(:,:) = 0._r8
      do indBin =  1, secNrBins-1
         ! fraction to stay in bin
         xfrac = (volume(indBin+1) - volume_new(indBin)) / &
              (volume(indBin+1) - volume(indBin))
         if (numberConc_old(indBin) < 1.e-30_r8) then
            xfrac = 1.0_r8
         end if
         if (xfrac <= 0._r8) then      ! if the fraction to stay is equal to
            ! less than zero, then the
            ! aerosols have grown too large
            ! for the next bin and we will
            ! want to decrease the time step
            ! to avoid this.
            !!XXG: is there any reason to not return immediately here?
            decrease_dt = .TRUE.
         end if

         !!XXG: Isn't this redundant?
         if (xfrac <= 0._r8) then
            decrease_dt = .TRUE.
         end if
         xfrac = max(0._r8, min(1._r8, xfrac))
         do indSpec =  1, secNrSpec
            numberConc_new(indSpec, indBin) = numberConc_new(indSpec, indBin) +     &
                 (xfrac * numberConc_old(indBin) * volfrac(indSpec,indBin))
            numberConc_new(indSpec, indBin+1) = numberConc_new(indSpec, indBin+1) + &
                 ((1.0_r8-xfrac) * numberConc_old(indBin) * volfrac(indSpec,indBin))
         end do
      end do

      ! Fraction to stay in sectional scheme?
      xfrac = (max_volume - volume_new(secNrBins)) /                       &
           (max_volume - volume(secNrBins))
      ! if less than or 0 % stays in bin, we must decrease timestep
      ! XXG: see note above with volume_new calc
      if (xfrac <= 0._r8) then
         decrease_dt = .TRUE.
      end if

      xfrac = max(0._r8, min(1._r8, xfrac))

      do indSpec = 1, secNrSpec
         numberConc_new(indSpec,secNrBins) = numberConc_new(indSpec,secNrBins) + &
              xfrac * numberConc_old(secNrBins) * &
              volfrac(indSpec,secNrBins)
         !XXG: comment on next line?
         leave_sec(indSpec) = & !massDistrib(chemistryIndex(secConstIndex(indSpec, secNrBins)))*(1-xfrac)
              max_volume * rhopart_sec(indSpec)/rhoAir * &   ! [m3_aer/#]*[kg_aer/m3_aer]/[kg_air/m3_air]--> [kg_aer/kg_air/#][m3_air]
              (1-xfrac) * numberConc_old(secNrBins) * &          ! *[#/m3_air] --> kg_aer/kg_air
              volfrac(indSpec,secNrBins)
      end do
      do indBin = 1, secNrBins
         do indSpec = 1, secNrSpec ! Assume (!XXG: Assume what??)
            ! Remove comments below?
            massDistrib(chemistryIndex(secConstIndex(indSpec,indBin))) = &! &!massDistrib(secConstIndex(indSpec,indBin))+&
                 rhopart_sec(indSpec)/rhoAir * &!* massfrac(indSpec,indBin)* numberConc_new(indBin)! &
                 numberConc_new(indSpec, indBin) * secMeanVol(indBin)
         end do
      end do

   end subroutine sec_moveMass

end module aero_sectional
