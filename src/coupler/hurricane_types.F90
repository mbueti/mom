#include <fms_platform.h>

module hurricane_types_mod
  use field_manager_mod,only: fm_field_name_len, fm_string_len
  use fms_mod,          only: write_version_number
  use mpp_mod,          only: FATAL, mpp_error
   
  implicit none
  
  private
  
  logical :: module_is_initialized=.false.
  character(len=128) :: version = &
     '$Id: hurricane_types.F90,v 16.0.2.8.6.1.2.3.6.1.18.1 2010/02/22 16:07:53 smg Exp $'
  character (len=128) :: tagname = &
     '$Name: mom4p1_riga_30jun2011_m1b $'
  
  type, public :: hurricane_type
    integer :: current_time = 0
    real :: u_lon(601)
    real :: u_lat(301)
    real :: v_lon(601)
    real :: v_lat(301)
    real :: u(601,301) = 0
    real :: v(601,301) = 0
    real, allocatable :: u_interp(:,:)
    real, allocatable :: v_interp(:,:)
    real, allocatable :: lon_interp(:,:)
    real, allocatable :: lat_interp(:,:)
  end type hurricane_type
  
  public hurricane_types_init
  
  contains
  
  subroutine hurricane_types_init() 
    if (module_is_initialized) then
      call mpp_error( FATAL, '==>Error: hurricane_types_init: module already initialized')
    end if
    module_is_initialized = .true.
    
    call write_version_number(version, tagname)
    
    return
  end subroutine hurricane_types_init
end module hurricane_types_mod
