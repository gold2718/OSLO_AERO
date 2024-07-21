module aero_sectional
   !! Shadow file for the full module in src_sectional
   implicit none
   private

   ! Public interface
   public :: aerosect_register
   public :: aerosect_init

   ! Public data
   logical, public, protected :: do_sectional_NPF = .false.
   integer, public, parameter :: secNrBins   = 0 ! nr of bins
   integer, public, parameter :: secNrSpec   = 0 ! number of condensing species

   real(r8), public, parameter :: rhopart_sec(secNrSpec)

   character(len=*), public, parameter :: secSpecNames(secNrSpec)
   character(len=*), public, parameter :: SpecNames(secNrSpec)

CONTAINS

   subroutine aerosect_register()
      ! Do nothing
   end subroutine aerosect_register

   subroutine aerosect_init()
      ! Do nothing
   end subroutine aerosect_init

end module aero_sectional
