module atmos_model_mod

use mpp_mod, only : mpp_npes, mpp_pe
use mpp_domains_mod, only : domain2d
use time_manager_mod, only : time_type

type, public :: land_ice_atmos_boundary_type
!variables of this type are declared by coupler_main, allocated by flux_exchange_init
!quantities going from land+ice to atmos
!       t         = surface temperature for radiation calculations
!       albedo    = surface albedo for radiation calculations
!       land_frac = fraction amount of land in a grid box
!       dt_t      = temperature tendency at the lowest level
!       dt_q      = specific humidity tendency at the lowest level
!       u_flux    = zonal wind stress
!       v_flux    = meridional wind stress
!       dtaudv    = derivative of wind stress w.r.t. the lowest level wind speed
!       u_star    = friction velocity
!       b_star    = bouyancy scale
!       q_star    = moisture scale
!       rough_mom = surface roughness (used for momentum
   real, dimension(:,:), pointer :: t, albedo, land_frac
   real, dimension(:,:), pointer :: dt_t, dt_q
   real, dimension(:,:), pointer :: u_flux, v_flux, dtaudv, u_star, b_star, q_star, rough_mom
   real, dimension(:,:,:), pointer :: data !collective field for "named" fields above
   integer :: xtype             !REGRID, REDIST or DIRECT
end type land_ice_atmos_boundary_type

type, public :: surf_diff_type

   real, pointer, dimension(:,:) :: dtmass,    &
        dflux_t,   & 
        dflux_q,   & 
        delta_t,   &
        delta_q
   
end type surf_diff_type

type, public :: atmos_data_type
   type (domain2d)               :: domain
   integer                       :: axes(4)
   real, pointer, dimension(:)   :: glon_bnd, glat_bnd,  &
        lon_bnd,  lat_bnd
   real, pointer, dimension(:,:) :: t_bot, q_bot, z_bot, p_bot,  &
        u_bot, v_bot, p_surf, gust,  &
        coszen, flux_sw, flux_lw,    &
        lprec, fprec
   type (surf_diff_type)         :: Surf_diff
   type (time_type)              :: Time, Time_step, Time_init
   integer, pointer              :: pelist(:)
   logical                       :: pe
end type atmos_data_type
  
!quantities going from land alone to atmos (none at present)
type, public :: land_atmos_boundary_type
!   private
   real, dimension(:,:), pointer :: data
end type land_atmos_boundary_type

!quantities going from ice alone to atmos (none at present)
type, public :: ice_atmos_boundary_type
!   private
   real, dimension(:,:), pointer :: data
end type ice_atmos_boundary_type
  
public atmos_model_init, atmos_model_end, &
     update_atmos_model_down,           &
     update_atmos_model_up,             &
     atmos_data_type, &
     land_ice_atmos_boundary_type, &
     land_atmos_boundary_type, &
     ice_atmos_boundary_type

contains


subroutine atmos_model_init (Atmos, Time_init, Time, Time_step)

  use fms_mod, only : read_data, field_size
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains,&
                              CYCLIC_GLOBAL_DOMAIN, mpp_get_data_domain,&
                              mpp_get_compute_domain
  use diag_manager_mod, only : diag_axis_init
  use diag_integral_mod, only : diag_integral_init
  use constants_mod, only : cp_air, hlv
  
type (atmos_data_type), intent(inout) :: Atmos
type (time_type), intent(in) :: Time_init, Time, Time_step

real, dimension(:), allocatable :: glon, glat
integer, dimension(4) :: siz
integer, dimension(2) :: layout
integer :: is, ie, js, je

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

call mpp_define_layout((/1,siz(1),1,siz(2)/), mpp_npes(), layout)
call mpp_define_domains((/1,siz(1),1,siz(2)/), layout, Atmos%domain, &
     xflags = CYCLIC_GLOBAL_DOMAIN)

call mpp_get_compute_domain(Atmos%domain,is,ie,js,je)

allocate(Atmos%lon_bnd(ie-is+2),&
         Atmos%lat_bnd(je-js+2))

Atmos%lon_bnd(:) = Atmos%glon_bnd(is:ie+1)
Atmos%lat_bnd(:) = Atmos%glat_bnd(js:je+1)

Atmos%axes(1) = diag_axis_init('lon',glon,'degrees_E','X','longitude',&
       set_name='atmos',domain2 = Atmos%domain)

Atmos%axes(2) = diag_axis_init('lat',glat,'degrees_N','Y','latitude',&
     set_name='atmos',domain2 = Atmos%domain)  

allocate(Atmos%t_bot(is:ie,js:je),&
     Atmos%q_bot(is:ie,js:je), &
     Atmos%z_bot(is:ie,js:je), &
     Atmos%p_bot(is:ie,js:je), &
     Atmos%u_bot(is:ie,js:je), &
     Atmos%v_bot(is:ie,js:je), &
     Atmos%p_surf(is:ie,js:je), &
     Atmos%gust(is:ie,js:je), &
     Atmos%coszen(is:ie,js:je), &
     Atmos%flux_sw(is:ie,js:je), &
     Atmos%flux_lw(is:ie,js:je), &
     Atmos%lprec(is:ie,js:je), &
     Atmos%fprec(is:ie,js:je))

Atmos%t_bot=273.0
Atmos%q_bot = 0.0
Atmos%z_bot = 10.0
Atmos%p_bot = 1.e5
Atmos%u_bot = 0.0
Atmos%v_bot = 0.0
Atmos%p_surf = 1.e5
Atmos%gust = 0.0
Atmos%coszen = 0.0
Atmos%flux_sw = 0.0
Atmos%flux_lw = 0.0
Atmos%lprec = 0.0
Atmos%fprec = 0.0

allocate(Atmos%Surf_diff%dtmass(is:ie, js:je) , &
         Atmos%Surf_diff%dflux_t(is:ie, js:je) , &
         Atmos%Surf_diff%dflux_q(is:ie, js:je) , &
         Atmos%Surf_diff%delta_t(is:ie, js:je) , &
         Atmos%Surf_diff%delta_q(is:ie, js:je) )

Atmos%Surf_diff%dflux_t = 0.0
Atmos%Surf_diff%dflux_q = 0.0
Atmos%Surf_diff%dtmass = 0.0
Atmos%Surf_diff%delta_t = 0.0
Atmos%Surf_diff%delta_q = 0.0

!------ initialize global integral package ------

    call diag_integral_init (Time_init, Time,  &
                             Atmos % lon_bnd,   Atmos % lat_bnd)
    
return

end subroutine atmos_model_init

subroutine atmos_model_end (Atmos)

type (atmos_data_type), intent(inout) :: Atmos

return

end subroutine atmos_model_end



subroutine update_atmos_model_down( Surface_boundary, Atmos )
!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs radiation, damping, and vertical diffusion of momentum,
!    tracers, and downward heat/moisture
!
!-----------------------------------------------------------------------

  type(land_ice_atmos_boundary_type), intent(inout) :: Surface_boundary
  type (atmos_data_type), intent(inout) :: Atmos

  return

end subroutine update_atmos_model_down

subroutine update_atmos_model_up( Surface_boundary, Atmos )
!-----------------------------------------------------------------------
!                       atmospheric driver
!    performs upward vertical diffusion of heat/moisture and
!    moisture processes
!
!-----------------------------------------------------------------------

   type(land_ice_atmos_boundary_type), intent(in) :: Surface_boundary
   type (atmos_data_type), intent(inout) :: Atmos
   
   return

end subroutine update_atmos_model_up

end module atmos_model_mod
