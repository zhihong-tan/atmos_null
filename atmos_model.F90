module atmos_model_mod

use mpp_mod, only : mpp_npes, mpp_pe
use mpp_domains_mod, only : domain2d
use time_manager_mod, only : time_type
use coupler_types_mod, only: coupler_2d_bc_type

implicit none
private

public atmos_data_type
public atmos_model_end
public atmos_model_init
public ice_atmos_boundary_type
public land_ice_atmos_boundary_type
public land_atmos_boundary_type
public surf_diff_type
public update_atmos_model_down
public update_atmos_model_up
public atm_stock_pe

!<PUBLICTYPE >
! This type should be defined in one spot and "used" from there
type surf_diff_type
  real, pointer, dimension(:,:)   :: dtmass  => NULL()
  real, pointer, dimension(:,:)   :: dflux_t => NULL()
  real, pointer, dimension(:,:)   :: delta_t => NULL()
  real, pointer, dimension(:,:)   :: delta_u => NULL()
  real, pointer, dimension(:,:)   :: delta_v => NULL()
  real, pointer, dimension(:,:,:) :: dflux_tr => NULL()   ! tracer flux tendency
  real, pointer, dimension(:,:,:) :: delta_tr => NULL()   ! tracer tendency
end type surf_diff_type
!</PUBLICTYPE >

!<PUBLICTYPE >
 type atmos_data_type
     type (domain2d)               :: domain             ! domain decomposition
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid 
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     real, pointer, dimension(:)   :: glon_bnd => NULL() ! global longitude axis grid box boundaries in radians.
     real, pointer, dimension(:)   :: glat_bnd => NULL() ! global latitude axis grid box boundaries in radians.
     real, pointer, dimension(:)   :: lon_bnd  => NULL() ! local longitude axis grid box boundaries in radians.
     real, pointer, dimension(:)   :: lat_bnd  => NULL() ! local latitude axis grid box boundaries in radians.
     real, pointer, dimension(:,:) :: t_bot    => NULL() ! temperature at lowest model level
     real, pointer, dimension(:,:,:) :: tr_bot => NULL() ! tracers at lowest model level, including specific humidity
     real, pointer, dimension(:,:) :: z_bot    => NULL() ! height above the surface for the lowest model level
     real, pointer, dimension(:,:) :: p_bot    => NULL() ! pressure at lowest model level
     real, pointer, dimension(:,:) :: u_bot    => NULL() ! zonal wind component at lowest model level
     real, pointer, dimension(:,:) :: v_bot    => NULL() ! meridional wind component at lowest model level
     real, pointer, dimension(:,:) :: p_surf   => NULL() ! surface pressure 
     real, pointer, dimension(:,:) :: gust     => NULL() ! gustiness factor
     real, pointer, dimension(:,:) :: coszen   => NULL() ! cosine of the zenith angle
     real, pointer, dimension(:,:) :: flux_sw  => NULL() ! net shortwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: flux_sw_dir            =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_dif            =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dir   =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_vis_dif   =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dir =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_down_total_dif =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_vis            =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_vis_dir        =>NULL()
     real, pointer, dimension(:,:) :: flux_sw_vis_dif        =>NULL()
     real, pointer, dimension(:,:) :: flux_lw  => NULL() ! net longwave flux (W/m2) at the surface
     real, pointer, dimension(:,:) :: lprec    => NULL() ! mass of liquid precipitation since last time step (Kg/m2)
     real, pointer, dimension(:,:) :: fprec    => NULL() ! ass of frozen precipitation since last time step (Kg/m2)
     logical,pointer,dimension(:,:):: maskmap  => NULL() ! A pointer to an array indicating which
                                                         ! logical processors are actually used for
                                                         ! the ocean code. The other logical
                                                         ! processors would be all land points and
                                                         ! are not assigned to actual processors.
                                                         ! This need not be assigned if all logical
                                                         ! processors are used.indicate if a domain region will be loaded.
     type (surf_diff_type)         :: Surf_diff          ! store data needed by the multi-step version of the diffusion algorithm
     type (time_type)              :: Time               ! current time
     type (time_type)              :: Time_step          ! atmospheric time step.
     type (time_type)              :: Time_init          ! reference time.
     integer, pointer              :: pelist(:) =>NULL() ! pelist where atmosphere is running.
     logical                       :: pe                 ! current pe.
     type(coupler_2d_bc_type)      :: fields   ! array of fields used for additional tracers
 end type
!</PUBLICTYPE >

!<PUBLICTYPE >
type land_ice_atmos_boundary_type
   ! variables of this type are declared by coupler_main, allocated by flux_exchange_init.
!quantities going from land+ice to atmos
   real, dimension(:,:),   pointer :: t              =>NULL() ! surface temperature for radiation calculations
   real, dimension(:,:),   pointer :: albedo         =>NULL() ! surface albedo for radiation calculations
   real, dimension(:,:),   pointer :: albedo_vis_dir =>NULL()
   real, dimension(:,:),   pointer :: albedo_nir_dir =>NULL()
   real, dimension(:,:),   pointer :: albedo_vis_dif =>NULL()
   real, dimension(:,:),   pointer :: albedo_nir_dif =>NULL()
   real, dimension(:,:),   pointer :: land_frac      =>NULL() ! fraction amount of land in a grid box 
   real, dimension(:,:),   pointer :: dt_t           =>NULL() ! temperature tendency at the lowest level
   real, dimension(:,:,:), pointer :: dt_tr          =>NULL() ! tracer tendency at the lowest level, including specific humidity
   real, dimension(:,:),   pointer :: u_flux         =>NULL() ! zonal wind stress
   real, dimension(:,:),   pointer :: v_flux         =>NULL() ! meridional wind stress
   real, dimension(:,:),   pointer :: dtaudu         =>NULL() ! derivative of wind stress w.r.t. the lowest level wind speed
   real, dimension(:,:),   pointer :: dtaudv         =>NULL() ! derivative of wind stress w.r.t. the lowest level wind speed
   real, dimension(:,:),   pointer :: u_star         =>NULL() ! friction velocity
   real, dimension(:,:),   pointer :: b_star         =>NULL() ! bouyancy scale
   real, dimension(:,:),   pointer :: q_star         =>NULL() ! moisture scale
   real, dimension(:,:),   pointer :: rough_mom      =>NULL() ! surface roughness (used for momentum)
   real, dimension(:,:,:), pointer :: data =>NULL() !collective field for "named" fields above
   integer :: xtype             !REGRID, REDIST or DIRECT
end type land_ice_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
type :: land_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>NULL() ! quantities going from land alone to atmos (none at present)
end type land_atmos_boundary_type
!</PUBLICTYPE >

!<PUBLICTYPE >
!quantities going from ice alone to atmos (none at present)
type :: ice_atmos_boundary_type
   real, dimension(:,:), pointer :: data =>NULL() ! quantities going from ice alone to atmos (none at present)
end type ice_atmos_boundary_type
!</PUBLICTYPE >
  
!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: atmos_model.F90,v 14.0 2007/03/15 22:01:09 fms Exp $'
character(len=128) :: tagname = '$Name: nalanda $'

contains

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_down">
!
! <OVERVIEW>
!   compute the atmospheric tendencies for dynamics, radiation, 
!   vertical diffusion of momentum, tracers, and heat/moisture.
! </OVERVIEW>
!
!<DESCRIPTION>
!   Called every time step as the atmospheric driver to compute the
!   atmospheric tendencies for dynamics, radiation, vertical diffusion of
!   momentum, tracers, and heat/moisture.  For heat/moisture only the
!   downward sweep of the tridiagonal elimination is performed, hence
!   the name "_down". 
!</DESCRIPTION>

!   <TEMPLATE>
!     call  update_atmos_model_down( Surface_boundary, Atmos )
!   </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.  
! </IN>

! <IN NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </IN>

subroutine update_atmos_model_down( Surface_boundary, Atmos )
!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs radiation, damping, and vertical diffusion of momentum,
!    tracers, and downward heat/moisture
!
!-----------------------------------------------------------------------

  type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
  type (atmos_data_type), intent(in) :: Atmos

  return

end subroutine update_atmos_model_down
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="update_atmos_model_up">
!
!-----------------------------------------------------------------------
! <OVERVIEW>
!   upward vertical diffusion of heat/moisture and moisture processes
! </OVERVIEW>

!<DESCRIPTION>
!   Called every time step as the atmospheric driver to finish the upward
!   sweep of the tridiagonal elimination for heat/moisture and compute the
!   convective and large-scale tendencies.  The atmospheric variables are
!   advanced one time step and tendencies set back to zero. 
!</DESCRIPTION>

! <TEMPLATE>
!     call  update_atmos_model_up( Surface_boundary, Atmos )
! </TEMPLATE>

! <IN NAME = "Surface_boundary" TYPE="type(land_ice_atmos_boundary_type)">
!   Derived-type variable that contains quantities going from land+ice to atmos.  
! </IN>

! <IN NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
!   These fields describe the atmospheric grid and are needed to
!   compute/exchange fluxes with other component models.  All fields in this
!   variable type are allocated for the global grid (without halo regions).
! </IN>
subroutine update_atmos_model_up( Surface_boundary, Atmos )
!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs upward vertical diffusion of heat/moisture and
!    moisture processes
!
!-----------------------------------------------------------------------

   type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
   type (atmos_data_type), intent(in) :: Atmos
   
   return

end subroutine update_atmos_model_up
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_init">
!
! <OVERVIEW>
! Routine to initialize the atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!     This routine allocates storage and returns a variable of type
!     atmos_boundary_data_type, and also reads a namelist input and restart file. 
! </DESCRIPTION>

! <TEMPLATE>
!     call atmos_model_init (Atmos, Time_init, Time, Time_step)
! </TEMPLATE>

! <IN NAME="Time_init" TYPE="type(time_type)" >
!   The base (or initial) time of the experiment.
! </IN>

! <IN NAME="Time" TYPE="type(time_type)" >
!   The current time.
! </IN>

! <IN NAME="Time_step" TYPE="type(time_type)" >
!   The atmospheric model/physics time step.
! </IN>

! <INOUT NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </INOUT>

subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)

  use fms_mod, only : read_data, field_size
  use mpp_domains_mod,   only : mpp_define_layout, mpp_define_domains
  use mpp_domains_mod,   only : CYCLIC_GLOBAL_DOMAIN, mpp_get_data_domain
  use mpp_domains_mod,   only : mpp_get_compute_domain
  use diag_manager_mod, only : diag_axis_init
  use diag_integral_mod, only : diag_integral_init
  use constants_mod, only : cp_air, hlv
  
type (atmos_data_type), intent(inout) :: Atmos
type (time_type), intent(in) :: Time_init, Time, Time_step

real, dimension(:), allocatable :: glon, glat
integer, dimension(4) :: siz
integer, dimension(2) :: layout
integer :: is, ie, js, je
integer :: ntprog

!---- set the atmospheric model time ------

   Atmos % Time_init = Time_init
   Atmos % Time      = Time
   Atmos % Time_step = Time_step
 
call field_size('INPUT/grid_spec','AREA_ATM',siz)
allocate(Atmos%glon_bnd(siz(1)+1))
allocate(Atmos%glat_bnd(siz(2)+1))
allocate(glon(siz(1)),glat(siz(2)))
call read_data('INPUT/grid_spec','xba',Atmos%glon_bnd)
call read_data('INPUT/grid_spec','yba',Atmos%glat_bnd)
call read_data('INPUT/grid_spec','xta',glon)
call read_data('INPUT/grid_spec','yta',glat)

Atmos%glon_bnd = Atmos%glon_bnd*atan(1.0)/45.0
Atmos%glat_bnd = Atmos%glat_bnd*atan(1.0)/45.0

if( ASSOCIATED(Atmos%maskmap) ) then
   layout(1) = size(Atmos%maskmap,1)
   layout(2) = size(Atmos%maskmap,2)
   call mpp_define_domains((/1,siz(1),1,siz(2)/), layout, Atmos%domain, &
        xflags = CYCLIC_GLOBAL_DOMAIN, maskmap = Atmos%maskmap , name='atmos model')
else
   call mpp_define_layout((/1,siz(1),1,siz(2)/), mpp_npes(), layout)
   call mpp_define_domains((/1,siz(1),1,siz(2)/), layout, Atmos%domain, &
        xflags = CYCLIC_GLOBAL_DOMAIN, name='atmos model')
end if

call mpp_get_compute_domain(Atmos%domain,is,ie,js,je)

allocate ( Atmos%lon_bnd(ie-is+2) )
allocate ( Atmos%lat_bnd(je-js+2) )

Atmos%lon_bnd(:) = Atmos%glon_bnd(is:ie+1)
Atmos%lat_bnd(:) = Atmos%glat_bnd(js:je+1)

Atmos%axes(1) = diag_axis_init('lon',glon,'degrees_E','X','longitude',&
       set_name='atmos',domain2 = Atmos%domain)

Atmos%axes(2) = diag_axis_init('lat',glat,'degrees_N','Y','latitude',&
     set_name='atmos',domain2 = Atmos%domain)  

allocate ( Atmos%t_bot(is:ie,js:je) )
allocate ( Atmos%tr_bot(is:ie,js:je,1) )        ! just one tracer dimension for q?
allocate ( Atmos%z_bot(is:ie,js:je) )
allocate ( Atmos%p_bot(is:ie,js:je) )
allocate ( Atmos%u_bot(is:ie,js:je) )
allocate ( Atmos%v_bot(is:ie,js:je) )
allocate ( Atmos%p_surf(is:ie,js:je) )
allocate ( Atmos%gust(is:ie,js:je) )
allocate ( Atmos%coszen(is:ie,js:je) )
allocate ( Atmos%flux_sw(is:ie,js:je) )
allocate ( Atmos % flux_sw_dir (is:ie,js:je) )
allocate ( Atmos % flux_sw_dif (is:ie,js:je) )
allocate ( Atmos % flux_sw_down_vis_dir (is:ie,js:je) )
allocate ( Atmos % flux_sw_down_vis_dif (is:ie,js:je) )
allocate ( Atmos % flux_sw_down_total_dir (is:ie,js:je) )
allocate ( Atmos % flux_sw_down_total_dif (is:ie,js:je) )
allocate ( Atmos % flux_sw_vis (is:ie,js:je) )
allocate ( Atmos % flux_sw_vis_dir (is:ie,js:je) )
allocate ( Atmos % flux_sw_vis_dif(is:ie,js:je) )
allocate ( Atmos%flux_lw(is:ie,js:je) )
allocate ( Atmos%lprec(is:ie,js:je) )
allocate ( Atmos%fprec(is:ie,js:je) )

Atmos%t_bot=273.0
Atmos%tr_bot = 0.0
Atmos%z_bot = 10.0
Atmos%p_bot = 1.e5
Atmos%u_bot = 0.0
Atmos%v_bot = 0.0
Atmos%p_surf = 1.e5
Atmos%gust = 0.0
Atmos%coszen = 0.0
Atmos%flux_sw = 0.0
Atmos%flux_lw = 0.0
Atmos % flux_sw_dir = 0.0
Atmos % flux_sw_dif = 0.0 
Atmos % flux_sw_down_vis_dir = 0.0 
Atmos % flux_sw_down_vis_dif = 0.0 
Atmos % flux_sw_down_total_dir = 0.0
Atmos % flux_sw_down_total_dif = 0.0
Atmos % flux_sw_vis = 0.0 
Atmos % flux_sw_vis_dir = 0.0 
Atmos % flux_sw_vis_dif = 0.0
Atmos%lprec = 0.0
Atmos%fprec = 0.0

ntprog = 1
allocate ( Atmos%Surf_diff%dtmass(is:ie, js:je) )
allocate ( Atmos%Surf_diff%dflux_t(is:ie, js:je) )
allocate ( Atmos%Surf_diff%delta_t(is:ie, js:je) )
allocate ( Atmos%Surf_diff%delta_u(is:ie, js:je) )
allocate ( Atmos%Surf_diff%delta_v(is:ie, js:je) )
allocate ( Atmos%Surf_diff%dflux_tr(is:ie, js:je, ntprog) )
allocate ( Atmos%Surf_diff%delta_tr(is:ie, js:je, ntprog) )

Atmos%Surf_diff%dtmass = 0.0
Atmos%Surf_diff%dflux_t = 0.0
Atmos%Surf_diff%delta_t = 0.0
Atmos%Surf_diff%delta_u = 0.0
Atmos%Surf_diff%delta_v = 0.0
Atmos%Surf_diff%dflux_tr = 0.0
Atmos%Surf_diff%delta_tr = 0.0

!------ initialize global integral package ------

    call diag_integral_init (Time_init, Time,  &
                             Atmos % lon_bnd,   Atmos % lat_bnd)
    
return

end subroutine atmos_model_init
! </SUBROUTINE>

!#######################################################################
! <SUBROUTINE NAME="atmos_model_end">
!
! <OVERVIEW>
!  termination routine for atmospheric model
! </OVERVIEW>

! <DESCRIPTION>
!  Call once to terminate this module and any other modules used.
!  This routine writes a restart file and deallocates storage
!  used by the derived-type variable atmos_boundary_data_type.
! </DESCRIPTION>

! <TEMPLATE>
!   call atmos_model_end (Atmos)
! </TEMPLATE>

! <IN NAME="Atmos" TYPE="type(atmos_data_type)">
!   Derived-type variable that contains fields needed by the flux exchange module.
! </IN>

subroutine atmos_model_end (Atmos)

type (atmos_data_type), intent(in) :: Atmos

return

end subroutine atmos_model_end
! </SUBROUTINE>

subroutine atm_stock_pe (Atm, index, value)

type (atmos_data_type), intent(inout) :: Atm
integer,                intent(in)    :: index
real,                   intent(out)   :: value

value = 0.0

end subroutine atm_stock_pe

end module atmos_model_mod
