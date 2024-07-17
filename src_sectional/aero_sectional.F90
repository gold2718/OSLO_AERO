module aero_sectional
   use shr_kind_mod, only: r8 => shr_kind_r8
   use chem_mods,    only: gas_pcnst
   use aerosoldef,   only: l_soa_a1, l_so4_a1

   implicit none
   private

   ! Public interface
   public :: aerosect_register

   ! Public data
   integer, public, parameter ::  secNrBins   = 5          ! nr of bins
   integer, public, parameter ::  secNrSpec   = 2          ! number of condensing species
   real(r8), public, parameter  :: rhopart_sec(secNrSpec) = (/ 1769.0_r8,1500.0_r8 /)             ! same as SO4_NA, SOA_NA

   character(len=*), public, parameter ::  secSpecNames(secNrSpec) = (/ 'SO4_SEC','SOA_SEC'/) ! names of condensed species
   character(len=*), public, parameter ::  SpecNames(secNrSpec)   = (/ 'SO4','SOA'/)          ! names of condensing species

   integer, public                     ::  secCoagulate_receiver(secNrSpec) !=(/ l_so4_a1,l_soa_a1/)  ! coagulation receiver

   real(r8),public                     :: secMeanD(secNrBins)               ![m] array holding mean diameter in each bin
                                                                            !holds the chemistry indices of the sectional scheme:
   integer, public                     :: secConstIndex(secNrSpec, secNrBins)
   integer, public                     :: autocoag_receiver_index(secNrBins, secNrBins)

   ! Private data
   real(r8), parameter        ::  max_diameter= 39.6e-9_r8 ! [m] volume median diameter of NPF background moede
                                                           ! calculated so that volume2number should be correct
   real(r8), parameter        ::  min_diameter= 5.0e-9_r8  ! [m] minumum diameter

CONTAINS

   subroutine aerosect_init()
      ! Sets up diameter of bins.

      real(r8)  :: d_rat
      integer   :: i                ! index

      secCoagulate_receiver=(/ l_so4_a1,l_soa_a1/)  ! number of condensing species

      ! Use discrite geometric distribution/volume-ratio size distrib:
      d_rat=(max_diameter/min_diameter)**(1._r8/(secNrBins))
      secMeanD(1)=min_diameter
      do i=2,secNrBins
         secMeanD(i) = secMeanD(i-1)*d_rat
      end do

   end subroutine aerosect_init


   subroutine aerosect_register()
      ! Sets up the chemistry indices for the tracers in the sectional scheme
      use constituents, only: cnst_get_ind

      integer              :: secInd,volInd
      character(len=20)    :: cnst_name
      call aerosect_init()  !Should be moved
      do secInd = 1, secNrBins
         do volInd = 1, secNrSpec !names constituents as 'SOA_SEC01'/'SO4_SEC01'
            write(cnst_name,'(A,I2.2)') trim(secSpecNames(volInd)), secInd
            call cnst_get_ind(trim(cnst_name), secConstIndex(volInd,secInd), abort=.true.)
         end do
      end do

   end subroutine aerosect_register


   subroutine aerosect_write2file(q, lchnk,ncol, pmid, temperature)
      ! Routine writes number concentration of aerosols to history
      use cam_history, only: outfld
      use ppgrid,      only: pcols, pver
      use aerosoldef,  only: chemistryIndex
      use physconst,   only: rair

      integer,  intent(in)    :: lchnk                        ! chunk identifier
      integer,  intent(in)    :: ncol                         ! number of columns

      real(r8), intent(in)    :: q(pcols,pver,gas_pcnst)      ! tmr [kg/kg]
      real(r8), intent(in)    :: pmid(pcols,pver)             !
      real(r8), intent(in)    :: temperature(pcols,pver)      ! [K] Temperature
      character(len=20)       :: field_name
      real(r8)                :: rhoAir                       ! [kg/m3] density of air

      integer                 :: indBin, indSpec, ind_sec,i,k ! indices
      real(r8)                :: num_conc(pcols,pver)         ! [#/m3] number concentration
      do indBin = 1, secNrBins
         !Go through all core species in that bin
         do indSpec = 1,secNrSpec
            ind_sec=chemistryIndex(secConstIndex(indSpec,indBin))
            do k=1,pver
               do i=1,ncol
                  rhoAir = pmid(i,k)/rair/temperature(i,k)
                  call sec_numberConc(q(i,k,ind_sec),indSpec, indBin,rhoAir,  num_conc(i,k))

               end do
            end do

            WRITE(field_name,'(A,A,I2.2)') 'nr',trim(secSpecNames(indSpec)),indBin
            call outfld(trim(field_name),num_conc, pcols, lchnk) !#
         end do
      end do

   end subroutine aerosect_write2file



subroutine sec_numberConc(mass, volNr, binNr , rhoAir, numberConc)
    ! Calculates the number concentration from the mass concentration
    implicit none

    real(r8), intent(in)  :: mass !kg/kg
    integer, intent(in) :: binNr ! bin_index
    real(r8), intent(in)  :: rhoAir
    real(r8), intent(out) :: numberConc !#/m3_air
    integer, intent(in) :: volNr
    !local:
    integer    :: volInd
    real(r8), parameter   :: pi=3.141592654_r8

    numberConc = mass/rhopart_sec(volNr)*rhoAir/ & ![kg_aer/kg_air]/[kg_aer/m3_aer]*[kg_air/m3_air]
                          (secMeanD(binNr)**3*pi/6._r8) ! /[m3_aer/#]--> #/m3_air
    if (mass .lt. 1.e-35) then
            numberConc=0.0_r8
    end if


end subroutine sec_numberConc



subroutine sec_moveMass(massDistrib, numberConc_old, leave_sec, rhoAir, modeDiam, decrease_dt)!,rhopart)
    ! Moves tracer mass from on bin to the other based on condensational/coagulation growth.
    ! Based on Jacobson Fundamentals of Atmospheric Modeling, second edition (2005),
    ! Chapter   13.5
    use aerosoldef, only : chemistryIndex
    implicit none

    real(r8), dimension(gas_pcnst),intent(inout)        :: massDistrib       ! mass in each tracer
    real(r8), dimension(secNrBins),intent(in)           :: numberConc_old    ! numbr concentration before growth
    real(r8), dimension(secNrSpec), intent(out)         :: leave_sec         ! the mass that leaves sectional scheme
    logical, intent(out)                                :: decrease_dt       ! if set to True, time step is divided
                                                                             ! and the procedure is re run
    real(r8)                                            :: rhoAir            ! Density of air
    real(r8)                                            :: modeDiam          ! not used.
    !real(r8),dimension(:), intent(in)                  :: rhopart


    real(r8), dimension(secNrSpec, secNrBins)          :: numberConc_new ! number concentration after growth
    real(r8), dimension(secNrBins)                     :: volume         ! volume of particle in bin
    real(r8), dimension(secNrBins)                     :: volume_new     ! volume after growth
    integer                                            :: indBin,indSpec
    real(r8), parameter                                :: pi = 3.141592654_r8
    real(r8)                                           :: xfrac ! fraction to stay in bin
    real(r8), dimension(secNrSpec,secNrBins)           :: volfrac  !


    decrease_dt=.FALSE.
    !compute volume in each bin with condensation (by mass) and by
    !numberconcentration
    do indBin = 1, secNrBins
            volume_new(indBin) = 0.0_r8
            volfrac(:,indBin)=0.0_r8
            do indSpec = 1, secNrSpec! calculate volume in each bin by mass/density! m3
                    if (numberConc_old(indBin)<1.e-30_r8) then
                            volume_new(indBin)=0.0_r8
                    else
                        volume_new(indBin) = volume_new(indBin) + massDistrib(chemistryIndex(secConstIndex(indSpec,indBin)))/&
                                rhopart_sec(indSpec) * rhoAir/&
                                (numberConc_old(indBin))
                    end if

                    volfrac(indSpec,indBin)=massDistrib(chemistryIndex(secConstIndex(indSpec,indBin)))/&
                            rhopart_sec(indSpec)*rhoAir
                            !kg/kg(air)*[kg(air)/m3(air)][kg/m3]--> m3/m3(air)
            end do ! calculate volume in each bin by numberconcentration (volume from before condenstion)
            volfrac(:,indBin)=volfrac(:,indBin)/(sum(volfrac(:,indBin))+1.E-50_r8)
            if (sum(volfrac(:,indBin))<1.e-50_r8) then
                    volfrac(:,indBin)=0.0_r8
            end if
            ! calculate volume in each bin by mass/density! m3
            volume(indBin) =  pi * secMeanD(indBin)**3/6._r8
            ! calculate volume in each bin by numberconcentration (volume from before condenstion)
    end do
    numberConc_new(:,:) = 0._r8
    do indBin =  1, secNrBins-1
            ! fraction to stay in bin
            xfrac=(volume(indBin+1)-volume_new(indBin)) &
                            /(volume(indBin+1)-volume(indBin))
            if (numberConc_old(indBin)<1.e-30) then
                    xfrac=1.0_r8
            end if
            if (xfrac .le. 0._r8) then      ! if the fraction to stay is equal to
                                            ! less than zero, then the
                                            ! aerosols have grown too large
                                            ! for the next bin and we will
                                            ! want to decrease the time step
                                            ! to avoid this.
                    decrease_dt=.TRUE.
            end if

            if (xfrac .le. 0._r8) then
                    decrease_dt=.TRUE.
            end if
            xfrac=max(0._r8, min(1._r8,xfrac))
            do indSpec= 1, secNrSpec
                    numberConc_new(indSpec, indBin) = numberConc_new(indSpec, indBin) + &
                                    xfrac*numberConc_old(indBin) &
                                    *volfrac(indSpec,indBin)
                    numberConc_new(indSpec, indBin+1) = numberConc_new(indSpec, indBin+1) + &
                                    (1-xfrac)*numberConc_old(indBin)   &
                                    *volfrac(indSpec,indBin)


            end do

    end do

    xfrac = (max_diameter**3 * pi/6.0_r8 - volume_new(secNrBins)) &
                            /(max_diameter**3*pi/6.0_r8-volume(secNrBins))

    ! if less than or 0 % stays in bin, we must decrease timestep
    if (xfrac .le. 0._r8) then
            decrease_dt=.TRUE.
    end if

    xfrac=max(0._r8, min(1._r8,xfrac))

    do indSpec=1, secNrSpec
            numberConc_new(indSpec,secNrBins) = numberConc_new(indSpec,secNrBins) + &
                                            xfrac * numberConc_old(secNrBins) &
                                            * volfrac(indSpec,secNrBins)
            leave_sec(indSpec) = & !massDistrib(chemistryIndex(secConstIndex(indSpec, secNrBins)))*(1-xfrac)
                    pi * max_diameter**3 / 6.0_r8 * rhopart_sec(indSpec)/rhoAir &   ! [m3_aer/#]*[kg_aer/m3_aer]/[kg_air/m3_air]--> [kg_aer/kg_air/#][m3_air]
                                            * (1-xfrac) * numberConc_old(secNrBins) &          ! *[#/m3_air] --> kg_aer/kg_air
                                            * volfrac(indSpec,secNrBins)
    end do
    do indBin=1,secNrBins
            do indSpec=1,secNrSpec !Assume
                    massDistrib(chemistryIndex(secConstIndex(indSpec,indBin))) = &! &!massDistrib(secConstIndex(indSpec,indBin))+&
                            rhopart_sec(indSpec)/rhoAir &!* massfrac(indSpec,indBin)* numberConc_new(indBin)! &
                            * numberConc_new(indSpec, indBin) * pi * secMeanD(indBin)**3/6.0_r8 !&
            end do
    end do



end subroutine sec_moveMass





end module aero_sectional
