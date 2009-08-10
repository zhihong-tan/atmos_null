module atmos_model_mod

use mpp_mod,           only : mpp_npes, mpp_pe, mpp_error, FATAL
use mpp_domains_mod,   only : domain2d
use mpp_domains_mod,   only : mpp_define_layout, mpp_define_domains
use mpp_domains_mod,   only : CYCLIC_GLOBAL_DOMAIN, mpp_get_data_domain
use mpp_domains_mod,   only : mpp_get_compute_domain
use fms_mod,           only : field_exist, read_data, field_size
use time_manager_mod,  only : time_type
use coupler_types_mod, only : coupler_2d_bc_type
use diag_manager_mod,  only : diag_axis_init
use diag_integral_mod, only : diag_integral_init
use constants_mod,     only : cp_air, hlv
use mosaic_mod,        only : get_mosaic_ntiles
use xgrid_mod,         only : grid_box_type
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
public atmos_model_restart

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
  real, pointer, dimension(:,:)     :: sst_miz => NULL()
end type surf_diff_type
!</PUBLICTYPE >

!<PUBLICTYPE >
 type atmos_data_type
     type (domain2d)               :: domain             ! domain decomposition
     integer                       :: axes(4)            ! axis indices (returned by diag_manager) for the atmospheric grid 
                                                         ! (they correspond to the x, y, pfull, phalf axes)
     real, pointer, dimension(:,:) :: glon_bnd => NULL() ! global longitude axis grid box boundaries in radians.
     real, pointer, dimension(:,:) :: glat_bnd => NULL() ! global latitude axis grid box boundaries in radians.
     real, pointer, dimension(:,:) :: lon_bnd  => NULL() ! local longitude axis grid box boundaries in radians.
     real, pointer, dimension(:,:) :: lat_bnd  => NULL() ! local latitude axis grid box boundaries in radians.
     real, pointer, dimension(:,:) :: t_bot    => NULL() ! temperature at lowest model level
     real, pointer, dimension(:,:,:) :: tr_bot => NULL() ! tracers at lowest model level, including specific humidity
     real, pointer, dimension(:,:) :: z_bot    => NULL() ! height above the surface for the lowest model level
     real, pointer, dimension(:,:) :: p_bot    => NULL() ! pressure at lowest model level
     real, pointer, dimension(:,:) :: u_bot    => NULL() ! zonal wind component at lowest model level
     real, pointer, dimension(:,:) :: v_bot    => NULL() ! meridional wind component at lowest model level
     real, pointer, dimension(:,:) :: p_surf   => NULL() ! surface pressure 
     real, pointer, dimension(:,:) :: slp      => NULL() ! sea level pressure 
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
     type(grid_box_type)           :: grid
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

character(len=128) :: version = '$Id: atmos_model.F90,v 17.0 2009/07/21 02:53:36 fms Exp $'
character(len=128) :: tagname = '$Name: quebec $'

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

type (atmos_data_type), intent(inout) :: Atmos
type (time_type), intent(in)          :: Time_init, Time, Time_step

real, dimension(:), allocatable       :: glon, glat, glon_bnd, glat_bnd
real, dimension(:,:),  allocatable    :: tmpx, tmpy, area
integer, dimension(4)                 :: siz
integer, dimension(2)                 :: layout
integer                               :: is, ie, js, je, i, j
integer                               :: nlon, nlat, ntiles
integer                               :: ntprog
character(len=128)                    :: grid_file = "INPUT/grid_spec.nc"
character(len=256)                    :: tile_file
character(len=256)                    :: atm_mosaic   ! land mosaic file
!---- set the atmospheric model time ------

Atmos % Time_init = Time_init
Atmos % Time      = Time
Atmos % Time_step = Time_step
 

if(field_exist(grid_file, 'AREA_ATM') ) then
   call field_size(grid_file,'AREA_ATM',siz)
   nlon = siz(1)
   nlat = siz(2)
   allocate(Atmos%glon_bnd(nlon+1,nlat+1))
   allocate(Atmos%glat_bnd(nlon+1,nlat+1))
   allocate(      glon_bnd(nlon+1))
   allocate(      glat_bnd(nlat+1))
   allocate(glon(nlon),glat(nlat))
   call read_data(grid_file,'xba',glon_bnd, no_domain=.true.)
   call read_data(grid_file,'yba',glat_bnd, no_domain=.true.)
   call read_data(grid_file,'xta',glon, no_domain=.true.)
   call read_data(grid_file,'yta',glat, no_domain=.true.)
else if( field_exist(grid_file, 'atm_mosaic_file') ) then ! read from mosaic file
   call read_data(grid_file, "atm_mosaic_file", atm_mosaic)     
   atm_mosaic = "INPUT/"//trim(atm_mosaic)
   ntiles = get_mosaic_ntiles(atm_mosaic)
   if(ntiles .NE. 1) call mpp_error(FATAL, 'land_model_init: ntiles should be 1 for atmos mosaic, contact developer')
   call read_data(atm_mosaic, "gridfiles", tile_file )
   tile_file = 'INPUT/'//trim(tile_file)
   call field_size(tile_file, "x", siz)
   if( mod(siz(1)-1,2) .NE. 0) call mpp_error(FATAL, "atmos_model_init:size(x,1) - 1 should be divided by 2")
   if( mod(siz(2)-1,2) .NE. 0) call mpp_error(FATAL, "atmos_model_init:size(x,2) - 1 should be divided by 2")
   nlon = (siz(1)-1)/2
   nlat = (siz(2)-1)/2
   allocate(Atmos%glon_bnd(nlon+1,nlat+1))
   allocate(Atmos%glat_bnd(nlon+1,nlat+1))
   allocate( glon_bnd(nlon+1), glat_bnd(nlat+1), glon(nlon), glat(nlat))
   !--- read the grid information on supergrid.
   allocate( tmpx(2*nlon+1, 2*nlat+1), tmpy(2*nlon+1, 2*nlat+1) )
   call read_data(tile_file, "x", tmpx, no_domain=.TRUE.)
   call read_data(tile_file, "y", tmpy, no_domain=.TRUE.)
   !--- make sure the grid is regular lat-lon grid.
   do j = 1, 2*nlat+1
      do i = 2, 2*nlon+1
         if(tmpy(i,j) .NE. tmpy(1,j)) call mpp_error(FATAL, "atmos_model_init:longitude is not uniform")
      end do
   end do
   do i = 1, 2*nlon+1
      do j = 2, 2*nlat+1
         if(tmpx(i,j) .NE. tmpx(i,1)) call mpp_error(FATAL, "atmos_model_init:latitude is not uniform")
      end do
   end do

   do i = 1, nlon+1
      glon_bnd(i) = tmpx(2*i-1,1)
   end do
   do j = 1, nlat+1
      glat_bnd(j) = tmpy(1, 2*j-1)
   end do
   do i = 1, nlon
      glon(i) = tmpx(2*i,1)
   end do
   do j = 1, nlat
      glat(j) = tmpy(1, 2*j)
   end do
   deallocate(tmpx, tmpy)
else
   call mpp_error(FATAL, 'atmos_model_init: both AREA_ATM and atm_mosaic_file do not exist in file '//trim(grid_file))
end if

do i = 1, nlon+1
   Atmos%glon_bnd(i,:) = glon_bnd(i)*atan(1.0)/45.0
end do

do j = 1, nlat+1
   Atmos%glat_bnd(:,j) = glat_bnd(j)*atan(1.0)/45.0
end do

if( ASSOCIATED(Atmos%maskmap) ) then
   layout(1) = size(Atmos%maskmap,1)
   layout(2) = size(Atmos%maskmap,2)
   call mpp_define_domains((/1,nlon,1,nlat/), layout, Atmos%domain, &
        xflags = CYCLIC_GLOBAL_DOMAIN, maskmap = Atmos%maskmap , name='atmos model')
else
   call mpp_define_layout((/1,nlon,1,nlat/), mpp_npes(), layout)
   call mpp_define_domains((/1,nlon,1,nlat/), layout, Atmos%domain, &
        xflags = CYCLIC_GLOBAL_DOMAIN, name='atmos model')
end if

call mpp_get_compute_domain(Atmos%domain,is,ie,js,je)

allocate ( Atmos%lon_bnd(ie-is+2,je-js+2) )
allocate ( Atmos%lat_bnd(ie-is+2,je-js+2) )
allocate ( area(ie-is+2,je-js+2) )

Atmos%lon_bnd(:,:) = Atmos%glon_bnd(is:ie+1, js:je+1)
Atmos%lat_bnd(:,:) = Atmos%glat_bnd(is:ie+1, js:je+1)

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
allocate ( Atmos%slp(is:ie,js:je) )
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
Atmos%slp = 1.e5
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

area = 0.0
!r2 = radius*radius
!allocate (slat(size(Atmos%lat_bnd(1,:))))
!slat = sin(blat)
!do j=1,jdim
!   do i=1,idim
!      area(i,j) = r2*(blon(i+1) - blon(i))*(slat(j+1) - slat(j))
!   end do
!end do

call diag_integral_init (Time_init, Time, Atmos%lon_bnd, Atmos%lat_bnd, area )
    
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

  !#######################################################################
  ! <SUBROUTINE NAME="atmos_model_restart">
  ! <DESCRIPTION>
  ! dummy routines.
  ! </DESCRIPTION>
  subroutine atmos_model_restart(Atmos, timestamp)
    type (atmos_data_type),   intent(inout) :: Atmos
    character(len=*),  intent(in)           :: timestamp


  end subroutine atmos_model_restart
  ! </SUBROUTINE>

subroutine atm_stock_pe (Atm, index, value)

type (atmos_data_type), intent(inout) :: Atm
integer,                intent(in)    :: index
real,                   intent(out)   :: value

value = 0.0

end subroutine atm_stock_pe

end module atmos_model_mod
