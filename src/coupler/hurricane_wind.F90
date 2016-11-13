module hurricane_wind_mod
  use mpp_mod,         only: mpp_npes, mpp_pe, mpp_root_pe, &
       mpp_error, stderr, stdout, stdlog, FATAL, NOTE, mpp_set_current_pelist, &
       mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sum, &
       CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_ROUTINE, lowercase, &
       input_nml_file

  use mpp_domains_mod, only: domain2d,mpp_define_mosaic, mpp_define_mosaic, mpp_get_compute_domain,mpp_get_global_domain
   
  use mpp_io_mod,      only: mpp_close, mpp_open, MPP_MULTI, MPP_SINGLE, MPP_OVERWR
  use mpp_io_mod,              only: MPP_NATIVE, MPP_RDONLY, MPP_DELETE

  use xgrid_mod,       only: put_to_xgrid, xgrid_init, xmap_type
  use grid_mod,        only: get_grid_cell_centers, get_grid_size, get_grid_ntiles

  use time_manager_mod, only: time_type,get_date,get_time,print_date,print_time

  use horiz_interp_type_mod, only: horiz_interp_type
  use horiz_interp_mod, only: horiz_interp, horiz_interp_new,horiz_interp_del,horiz_interp_init,horiz_interp_end

  use fms_mod,                 only: open_namelist_file, field_exist, file_exist, check_nml_error
  use fms_mod,                 only: uppercase, error_mesg, write_version_number
  use fms_mod,                 only: fms_init, fms_end, stdout
  use fms_mod,                 only: read_data, write_data

  use hurricane_types_mod,  only: hurricane_type, hurricane_types_init
  use atmos_model_mod,         only: atmos_data_type

  implicit none
  include 'netcdf.inc'
  
  private

  character(len=48), parameter :: module_name = "hurricane_wind_mod"

  public  :: hurricane_wind_init,       &
             hurricane_wind_step,       &
             hurricane_wind_profile,    &
             hurricane_wind_verinterp,  &
             hurricane_wind_interp1d,   &
             hurricane_wind_date2day,   &
             hurricane_wind_day2date,   &
             hurricane_wind_output

  !-----------------------------------------------------------------------
   character(len=128) :: version = '$Id: hurricane_wind.F90,v 18.0.4.1.4.1.2.1 2011/06/23 13:24:00 m1b Exp $'
   character(len=128) :: tag = '$Name: riga_201012 $'
  !-----------------------------------------------------------------------

  !---- exchange grid maps -----

  type(xmap_type), save :: xmap_sfc, xmap_runoff

  integer         :: n_xgrid_sfc,  n_xgrid_runoff

  !-----------------------------------------------------------------------
  !-------- namelist (for diagnostics) ------

  character(len=14), parameter :: mod_name = 'hurricane_wind'

  logical :: first_static = .true.
  logical :: do_init = .true.
  integer :: remap_method = 1
  
  integer*2 :: dti=1
  integer   :: isplit = 1
  integer*2 :: ihours
  
  namelist /hurricane_wind_nml/ dti, isplit, ihours

  contains
  
  subroutine hurricane_wind_init(Hurricane,Atm,time_ext)    
    use hurricane_wind_common
  
  	type(hurricane_type), intent(inout) :: Hurricane
  	type(atmos_data_type), intent(in) :: Atm
  	type(time_type),      intent(in)    :: Time_ext
  	integer :: unit,ierr,io
  	integer :: outunit, logunit,month,day,hour,minute,second,date,nlon,nlat
  	real julday_ext
  
    call hurricane_types_init()

    outunit = stdout(); logunit = stdlog()
!----- read namelist -------

#ifdef INTERNAL_FILE_NML
    	read (input_nml_file, hurricane_wind_nml, iostat=io)
#else
    unit = open_namelist_file()
    ierr=1; do while (ierr /= 0)
       read  (unit, nml=hurricane_wind_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'hurricane_wind_nml')
    enddo
10  call mpp_close(unit)
#endif

!----- write namelist to logfile -----
    call write_version_number (version, tag)
    if( mpp_pe() == mpp_root_pe() )write( logunit, nml=hurricane_wind_nml )
    
    call get_date(Time_ext,year, month, day, hour, minute, second)
    date=hour+100*(day+100*(month+100*(year)))
    call hurricane_wind_date2day(year,time0,date)
    
    

    PI=3.1415927E0
    RAMP=1.E0
    SMALL=1.E-10
    BETA=1.98E-11
    GRAV=9.806E0
    UMOL=2.E-5

    REARTH=6371.E3
    LATMIN=minval(Atm%lat_bnd)
    LATMAX=maxval(Atm%lat_bnd)
    LONGMIN=minval(Atm%lon_bnd)
    LONGMAX=maxval(Atm%lon_bnd)

    IINT=0
    TIME=0.

!--- Model contral parameters----------------------------
!--- see file PARAMETERS.inp for meanings of parameters
191 format(a15)

    ISPADV=5
    TPRNU=1.   
    SMOTH=0.1
    HORCON=0.1
    TBIAS=10.
    SBIAS=35.
    ISWTCH=100000
    
    call get_grid_size('ATM',1,nlon,nlat)
  	allocate(Hurricane%u_interp(nlon,nlat))
  	allocate(Hurricane%v_interp(nlon,nlat))
  	allocate(Hurricane%lon_interp(nlon,nlat))
  	allocate(Hurricane%lat_interp(nlon,nlat))

    call get_grid_cell_centers('ATM',1,Hurricane%lon_interp,Hurricane%lat_interp)
  	
!
!--------------------------------------------------------------------------
!     READ NAMELIST (DISCONNECTED HERE) FOR PROBLEM PARAMETERS: 
!       MODE = 2; 2-D CALCULATION (BOTTOM STRESS CALCULATED IN ADVAVE)
!              3; 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
!              4; 3-D CALCULATION WITH T AND S HELD FIXED
!       NBC = TYPE OF THERMO. B.C.; SEE SUBROUTINE PROFT
!       DTI = INTERNAL TIME STEP
!       DTE = EXTERNAL TIME STEP
!       ISPLIT = DTI/DTE
!       NREAD=0, NO RESTART INPUT FILE; NREAD=1, RESTART
!       IPRTD1 = PRINT INTERVAL IN DAYS; AFTER ISWTCH DAYS, PRINT
!                  INTERVAL = IPRTD2
!       IDAYS = LENGTH OF RUN IN DAYS
!       ISPADV = STEP INTERVAL WHERE EXTERNAL MODE ADVECTIVE TERMS ARE
!                  NOT UPDATED
!       HORCON = CONSTANT IN SMAGORINSKY HORIZONTAL DIFFUSIVITY
!       TPNU = HORIZONTAL DIFFUSIVITY PRANDTL NUMBER
!       SMOTH = CONSTANT IN TIME SMOOTHER TO PREVENT SOLUTION SPLITTING
!
!     READ(5,PRMTR)
!--------------------------------------------------------------------------
!
    DTE=DTI/real(ISPLIT)
    DTE2=DTE*2
    DTI2=DTI*2
    IPRINT=IPRTH1*3600/INT(DTI)
    ISWTCH=ISWTCH*24*3600/INT(DTI)   
    IEND=IHOURS*3600/INT(DTI)
!
!----------------------------------------------------------------------
!             ESTABLISH PROBLEM CHARACTERISTICS
!          ****** ALL UNITS IN M.K.S. SYSTEM ******
!      F,BLANK AND B REFERS TO FORWARD,CENTRAL AND BACKWARD TIME LEVELS.
!----------------------------------------------------------------------
!
    DAYI=1.E0/86400.E0
!
!--- ASSIGNING SPACING AND CORIOLIS
!--- STEP IN LONGITUDE(PHI) IS (LONGMAX-LONGMIN)/(IM-1)
!--- STEP IN LATTITUDE(LAMBDA) IS (LATMAX-LATMIN)/(JM-1)
!
    do J=1,JM
      do I=1,IM
        PHI=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-1)
        LAMBDA=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
        LAMBDA=LAMBDA
        COR(I,J)=4.*3.1415927/(24.*60.*60.)*SIN(ABS(LAMBDA))
        DX(I,J)=REARTH*COS(LAMBDA)*(LONGMAX-LONGMIN)/(IM-1)
        DX(I,J)=DX(I,J)*2.*3.1415927/360.
        DY(I,J)=REARTH*(LATMAX-LATMIN)/(JM-1)
        ART(I,J)=DX(I,J)*DY(I,J)
      end do
    end do

    do J=2,JM
      do I=2,IM
        ARU(I,J)=.25*(DX(I,J)+DX(I-1,J))*(DY(I,J)+DY(I-1,J))
        ARV(I,J)=.25*(DX(I,J)+DX(I,J-1))*(DY(I,J)+DY(I,J-1))
      end do
    end do
        
    do J=1,JM
      ARU(1,J)=ARU(2,J)
      ARV(1,J)=ARV(2,J)
    end do
        
    do I=1,IM
      ARU(I,1)=ARU(I,2)
      ARV(I,1)=ARV(I,2)
    end do

    TIME=TIME0
	
    return
  end subroutine hurricane_wind_init

  subroutine hurricane_wind_step(Hurricane,Atm,u_atm,v_atm,xmap_sfc,Time_ext)
!----------------------------------------------------------------------
!     THIS SUBROUTINE DOES ONE STEP IN TIME OF THE OCEAN MODEL
!----------------------------------------------------------------------
!
    use hurricane_wind_common
    
    type(hurricane_type),   intent(inout) :: Hurricane
    type(atmos_data_type),  intent(inout) :: Atm
    type(time_type),        intent(in)    :: Time_ext
    type(xmap_type),			 intent(inout)	:: xmap_sfc
  	real, intent(inout) :: u_atm(:),v_atm(:)

    type(horiz_interp_type) :: Interpolant

    integer :: nlon,nlat,ntiles,unit
    integer :: is,ie,js,je,isize,jsize,remap_method

    real, allocatable :: glon(:,:),glat(:,:),u_wind(:,:),v_wind(:,:)
  	real              :: hu(im,jm),hv(im,jm),hlon(im),hlat(jm)
  	real              :: hlon2(im,jm),hlat2(im,jm)
  	real, allocatable :: ex_u_hur(:),ex_v_hur(:)
  	  	
    integer outunit,month,day,hour,minute,second,date
		integer :: readunit, eof, logunit
	  character*3 :: stormname
		
    
    remap_method=2
  	nlon=size(Atm%u_bot,1)
  	nlat=size(Atm%u_bot,2)

  	allocate(glon(nlon,nlat))
  	allocate(glat(nlon,nlat))
  	allocate(u_wind(nlon,nlat))
  	allocate(v_wind(nlon,nlat))
  	allocate(ex_u_hur(size(u_atm)))
  	allocate(ex_v_hur(size(v_atm)))
  	
    call get_date(Time_ext,year,month,day,hour,minute,second)
    date=hour+100*(day+100*(month+100*(year)))
    call  hurricane_wind_date2day(year,time,date)
    
    ! TIME=dayi*DTI*real(Hurricane%current_time)+TIME0
    ! Hurricane%current_time = Hurricane%current_time + 1
    ! TIME=DAYI*DTI*real(IINT)+TIME0
    RAMP=1.E0
    
    outunit = stdout(); logunit = stdlog()
		readunit=16

	!-----------WRAP THE PROFILE CALL IN A LOOP WHICH READS FROM storm_list, and passes as a pathname----------------
	  open(readunit,file='INPUT/storm_list',status='old')
		eof = 0
		do while(eof.eq.0)
			read(readunit,'(A3)') stormname
			if (stormname.ne.'000') then
				call hurricane_wind_profile(hu,hv,hlon,hlat,stormname)
				
		    do i=1,im
		      do j=1,jm
		        hlon2(i,j)=hlon(i)
		        hlat2(i,j)=hlat(j)
		      end do
		    end do
    
		    do i=1,nlon
		      do j=1,nlat
		        glon(i,j)=(Atm%lon_bnd(i,j)+Atm%lon_bnd(i+1,j)+Atm%lon_bnd(i,j+1)+Atm%lon_bnd(i+1,j+1))/4
		        glat(i,j)=(Atm%lat_bnd(i,j)+Atm%lat_bnd(i+1,j)+Atm%lat_bnd(i,j+1)+Atm%lat_bnd(i+1,j+1))/4
		      end do
		    end do
      
		    call horiz_interp_init
		    call horiz_interp_new(Interpolant,hlon,hlat,glon,glat,1,'bilinear',grid_at_center=.true.,src_modulo=.true.)
		    call horiz_interp(Interpolant,hu,u_wind)
		    call horiz_interp_del(Interpolant)
		    call horiz_interp_end
    
		    call horiz_interp_init
		    call horiz_interp_new(Interpolant,hlon,hlat,glon,glat,1,'bilinear',grid_at_center=.true.,src_modulo=.true.)
		    call horiz_interp(Interpolant,hv,v_wind)
		    call horiz_interp_del(Interpolant)
		    call horiz_interp_end
    
		    call put_to_xgrid(u_wind,'ATM',ex_u_hur,xmap_sfc,remap_method=remap_method)
		    call put_to_xgrid(v_wind,'ATM',ex_v_hur,xmap_sfc,remap_method=remap_method)
    

		    do i=1,size(u_atm)
		      if(sqrt(ex_u_hur(i)**2+ex_v_hur(i)**2).gt.sqrt(u_atm(i)**2+v_atm(i)**2)) then
		        u_atm(i) = ex_u_hur(i)
		        v_atm(i) = ex_v_hur(i)
		      end if
		    end do
			else
				eof = 1
			end if
		end do
		close(readunit)
    
		deallocate(glon)
    deallocate(glat)
    deallocate(u_wind)
    deallocate(v_wind)		
    
!
!--------------------------------------------------------------------
!           BEGIN PRINT SECTION
!--------------------------------------------------------------------
!
    ! call hurricane_wind_output(hu,hv)
    
!
!----------------------- END PRINT SECTION -----------------------------
!
    return
  end subroutine hurricane_wind_step

  subroutine hurricane_wind_profile(u,v,hlon,hlat,stormname)
    use hurricane_wind_common
    
    real, intent(inout) 			:: u(im,jm),v(im,jm),hlon(im),hlat(jm)
	  character*3, intent(in)  :: stormname
  	
    INTEGER PLN
    PARAMETER(PLN=100)

    integer lat,long,mx,rmw,ipmax,nm1,np1,ii
    integer*4 startd,date
    integer day,month,hour
    integer garb(5),Rd1(4),Rd2(4)
    character*19 name
    character*1 letter,lew,lns,ns_flag
  	integer :: outunit, logunit, readunit

    real,DIMENSION(pln)   :: X(PLN),Y(PLN),TM(PLN),PRES(PLN),PRES0(PLN),RMAXa(PLN),WSPMAX(PLN),RADCLS(PLN)         
    real,dimension(pln,4) :: R18v(PLN,4),R26v(PLN,4)
    real,dimension(5)     :: Rref18v(5),Rref26v(5),alphv(5)
    real,DIMENSION(14)    :: RAD(14),WS(14),RADM(14),WSM(14),ANGL(14)
    real, DIMENSION(1)    :: pres1,pres2,rcls,wsmax,Rref18,Rref26,rangl,alpha,timev,F0,L0,rmax,r
    REAL T1,T2,F1,L1,A7,B,E,DELP,x0,y0,x1,x2,y1,y2,utx,uty,utxa,utya
    REAL DELTAX,DELTAX1,DELTAY,DELTAY1,DXDY,wind_scale,roa,end,ws18,ws26
    real wnd,wf,wr,wx,wy,wm,cd,tmax,wnd0,wx0,wy0

    data RAD/0.,.4,.7,.8,.95,1.,1.35,2.7,4.05,5.4,6.75,8.1,10.8,13.5/
    data WS/0.,.1,.5,.8,.95,1.,.97,.72,.54,.44,.4,.36,.27,.23/
    data ANGL/0.,2.,4.,6.,7.,7.,14.,23.,24.,22.,21.,21.,21.,21./

    WIND_SCALE=1.
    ROA=1.28
    RMAX=50.E3
    PI=3.1415927
    E=exp(1.)
  
    outunit = stdout(); logunit = stdlog()
    
    readunit=15
		ns_flag='N'
		
!
!----------------------- Reading message file --------------------
17  format(A19,I6,1x,I4,1x,I3,A1,1x,I4,A1,1x,I3,1x,I3,3I5,  &
           1x,i2,1x,I3,1x,I4,1x,I4,1x,I4,1x,I4,1x,I4, &
           1x,I4,1x,I4,1x,I4)
    open(readunit,file='INPUT/track_'//stormname,status='old')
    end=0.
    I=0
    do while(end.eq.0)
      read(readunit,17) name,date,hour,lat,lns,long,lew,garb,mx,rmw,Rd1,Rd2
      if(date.ne.0) then
        I=I+1
        date=date*100+hour/100
        call hurricane_wind_date2day(year,julday,date)
        TM(i)=julday
				X(i)=long/10.
				if(lew.eq.'W') X(i)=-X(i)
				if(lew.eq.'E'.and.X(i).gt.80.) X(i)=X(i)-360.
				Y(i)=lat/10.
				if(lns.eq.'S') Y(i)=-Y(i)
				if(lns.eq.'S') ns_flag=lns
        PRES(i)   =real(garb(3))
        PRES0(i)  =real(garb(4))
	      RADCLS(i) =real(garb(5))
        WSPMAX(i) =real(mx)
        RMAXa(i)  =real(rmw)
        do n=1,4
          n1=n+(1-mod(n,2))*sign(2,3-n)
          R18v(i,n)=Rd1(n1)
          R26v(i,n)=Rd2(n1)
          if(wspmax(i).le.26.or.R26v(i,n).le.RMAXa(i)) R26v(i,n)=-999
          if(wspmax(i).le.18.or.R18v(i,n).le.RMAXa(i)) R18v(i,n)=-999
          if(R26v(i,n).gt.R18v(i,n)) R26v(i,n)=-999
        end do
      else
        end=1.
      end if
    end do
    ipmax=I
    close(readunit)
		if((TM(1).le.time).and.(TM(ipmax).ge.time)) then
!
!--------------------- Calculating starting day -----------------
!
    ! call hurricane_wind_date2day(year,julday,startd)
    ! do i=1,ipmax
    !   TM(i)=TM(i)-time0
    ! end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  INTERPOLATION OF HURRICANE PATH TO DETERMINE THE CURRENT POSITION

   	 CMP=TM(ipmax)
   	 
   	 timev(1)=time
   	 call hurricane_wind_verinterp(ipmax,1,TM,X,TIMEv,F0)
   	 call hurricane_wind_verinterp(ipmax,1,TM,Y,TIMEv,L0)
   	 xcen=F0(1)
   	 ycen=L0(1)
   	 call hurricane_wind_verinterp(ipmax,1,TM,PRES,TIMEv,PRES1)
   	 call hurricane_wind_verinterp(ipmax,1,TM,PRES0,TIMEv,PRES2)
   	 call hurricane_wind_verinterp(ipmax,1,TM,RADCLS,TIMEv,RCLS)
   	 DELP=(PRES2(1)-PRES1(1))*100.
   	 call hurricane_wind_verinterp(ipmax,1,TM,WSPMAX,TIMEv,WSMAX)
   	 prsmin=PRES1(1)
   	 wndmax=WSMAX(1)
   	 WSMAX=WSMAX*WIND_SCALE
   	 WS18=18*WIND_SCALE
   	 WS26=26*WIND_SCALE
   	 call hurricane_wind_verinterp(ipmax,1,TM,RMAXa,TIMEv,RMAX)
   	 time=timev(1)
   	 RMAX=RMAX*1.e3
   	 timev(1)=time
   	 do n=1,4
   	   call hurricane_wind_interp1d(-999.0,ipmax,1,TM,R18v(1,n),TIMEv,Rref18v(n))
   	   if(Rref18v(n).ne.-999) Rref18v(n) = Rref18v(n)*1.e3
   	   call hurricane_wind_interp1d(-999.0,ipmax,1,TM,R26v(1,n),TIMEv,Rref26v(n))
   	   if(Rref26v(n).ne.-999) Rref26v(n) = Rref26v(n)*1.e3
   	   alphv(n) = (n-1)*pi/2
   	 end do
   	 time=timev(1)
   	 do n=2,6
   	   n1=mod(n-1,4)+1
   	   nm1=mod(n-2,4)+1
   	   np1=mod(n,4)+1
   	   if(Rref18v(n1).eq.-999) then
   	     if(Rref18v(nm1).ne.-999) then
   	       if(Rref18v(np1).ne.-999) then
   	         Rref18v(n1)=0.5*(Rref18v(nm1)+Rref18v(np1))
   	       else
   	         Rref18v(n1)=Rref18v(nm1)
   	       end if
   	     else
   	       if(Rref18v(np1).ne.-999) then
   	         Rref18v(n1)=Rref18v(np1)
   	       else
   	         Rref18v(n1)=-999
   	       end if
   	     end if
   	   end if
   	   if(Rref26v(n1).eq.-999) then
   	     if(Rref26v(nm1).ne.-999) then
   	       if(Rref26v(np1).ne.-999) then
   	         Rref26v(n1)=0.5*(Rref26v(nm1)+Rref26v(np1))
   	       else
   	         Rref26v(n1)=Rref26v(nm1)
   	       end if
   	     else
   	       if(Rref26v(np1).ne.-999) then
   	         Rref26v(n1)=Rref26v(np1)
   	       else
   	         Rref26v(n1)=-999
   	       end if
   	     end if
   	   end if
   	 end do
   	
   	 Rref18v(5) = Rref18v(1)
   	 Rref26v(5) = Rref26v(1)
   	 alphv(5) = alphv(4)+pi/2
!  	    
   	 x0=f0(1)
   	 y0=l0(1)
   	 F0=F0*2.*PI/360
   	 L0=L0*2.*PI/360
!--	- F0,L0 ARE THE CARRENT (LONGITUDE,LATTITUDE) COORDINATES OF THE HURRICANE
!
!--	- CALCULATING UTX AND UTY (HURRICANE SPEED)
!
   	 cmp=tm(ipmax)
   	 do i=1,ipmax
   	   if(abs(tm(i)-time).le.cmp) then
   	     cmp=abs(tm(i)-time)
   	     ii=i
   	   end if
   	 end do
!--	------- falk 01-05-04 not to go out bnd
   	 if((tm(ii)-time).le.0.and.ii.ne.ipmax) then
   	   t1=tm(ii)
   	   t2=tm(ii+1)
   	   x1=x(ii)
   	   x2=x(ii+1)
   	   y1=y(ii)
   	   y2=y(ii+1)
   	 else
   	   t2=tm(ii)
   	   t1=tm(ii-1)
   	   x2=x(ii)
   	   x1=x(ii-1)
   	   y2=y(ii)
   	   y1=y(ii-1)
   	 end if
   	 deltax1=rearth*cos(l0(1))*(x2-x1)*2.*pi/360
   	 deltay1=rearth*(y2-y1)*2.*pi/360
   	 utx=deltax1/((t2-t1)*24.*3600.)
   	 uty=deltay1/((t2-t1)*24.*3600.)
!  	 
!--	- CALCULATING PARAMETERS FOR WIND PROFILE FORMULA
!
   	 B=WSMAX(1)**2*E*ROA/DELP
   	 A7=RMAX(1)**B
   	 do i=1,14
   	   RADM(I)=RMAX(1)*RAD(I)
   	 end do
   	
   	 DO J=1,JM
   	   DO I=1,IM
   	 

!  	     CALCULATING V-WIND STRESS FOR I,J POINT OF V-VELOCITY GRID
   	     F1=LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-0.5)
   	     hlon(i)=F1
   	     L1=LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1)
   	     hlat(j)=L1
   	     DELTAX=REARTH*COS(L0(1))*(F1-F0(1))
   	     
!--	-----falk 01-05-04 check DELTAX
   	     if(abs(DELTAX).lt.1.e-8) DELTAX=1.e-8
   	     DELTAY=REARTH*(L1-L0(1))
   	     DXDY=DELTAX*DELTAY
   	     R(1)=SQRT(DELTAX**2+DELTAY**2)
   	     alpha(1)=atan(abs(DELTAY/DELTAX))*sign(1.,DXDY)+(1 - sign(1.,DXDY))*pi/2 + (1 - sign(1.,DELTAY))*pi/2
   	     if(alpha(1).ge.pi/4) then
   	       alpha = alpha - pi/4
   	     else
   	       alpha = alpha - pi/4 + 2*pi
   	     end if
   	     call hurricane_wind_verinterp(5,1,alphv,Rref18v,alpha,Rref18)
   	     call hurricane_wind_verinterp(5,1,alphv,Rref26v,alpha,Rref26)
   	     call hurricane_wind_verinterp(14,1,radm,angl,R,RANGL)

!  	     CALCULATING WIND SPEED 
   	     if(Rref18(1).le.0.and.Rref26(1).le.0) then
	 	       if(R(1).le.(RCLS(1)*5000.))then
						 WND0=R(1)*COR(I,J)/2.
   	         WND=SQRT(A7*B*DELP*EXP(-A7/R(1)**B)/(ROA*R(1)**B)+R(1)**2*COR(I,J)**2/4.)-WND0
   	         UTXa=UTX/2.
   	         UTYa=UTY/2.
	 	       else
	 	         WND=0
	 	         UTXa=0
   	         UTYa=0
	 	       end if
   	     else
   	       WND=EXPWND(R(1),RMAX(1),Rref18(1),Rref26(1),WSMAX(1),WS18,WS26)
   	       UTXa=0.
   	       UTYa=0.
   	     end if
   	     RANGL=RANGL*PI/180.
   	     WF= WND*COS(RANGL(1))
   	     WR=-WND*SIN(RANGL(1))
   	     WX0=WR*(DELTAX)/R(1)-WF*(DELTAY)/R(1)
   	     WY0=WF*(DELTAX)/R(1)+WR*(DELTAY)/R(1)
				 if(ns_flag.eq.'S') WX0=-WX0
				 if(ns_flag.eq.'S') WY0=-WY0
				 WX=WX0+UTXa
				 WY=WY0+UTYa
   	     u(i,j)=WX
   	     v(i,j)=WY
   	   end do
   	 end do
	  else
			do I=1,IM
				do J=1,JM
					u(I,J)=0
					v(I,J)=0
   	      hlon(i)=(LONGMIN+(I-1)*(LONGMAX-LONGMIN)/(IM-0.5))
   	      hlat(j)=(LATMIN+(J-1)*(LATMAX-LATMIN)/(JM-1))
				end do
			end do
		end if

230 format(151e16.7)
231 format(10F8.3)
!
!---falk 01-05-04 add call avrtau
    ! call avrtau(x0,y0,0)
!
    return
  end subroutine hurricane_wind_profile

  subroutine avrtau(xln,ylt,mig)
    use hurricane_wind_common

    integer, parameter :: timesm4=12.,timesm=12.
    real,    parameter :: RAVR=100.e3
!---This subr. is called from WIND in phase4 (mig=0)
!---and from atmos2ocean.f (mig=1) in coupled run
!   REAL LATMIN,LATMAX,LONGMIN,LONGMAX
!   COMMON/sphere/REARTH,LATMIN,LATMAX,LONGMIN,LONGMAX
    real xln,ylt,tavr,counter,x1,y1,x0,y0,r,deltax,deltay
    real xlnc,yltc,tauavrp,taumaxp,tauabs,wucon
!   real tauavr,taumax,awucon,bwucon
    real RRCT,xrct,yrct
    integer irct,jrct
    integer mig
!    
201 format('avrtau: LATMIN,LATMAX,LONGMIN,LONGMAX=',4f7.2)

    pi=3.1415927
!---save previous tauavr, taumax
    tauavrp=tauavr
    taumaxp=taumax
!---for coupled run use TC position from TVARY.h
!---for phase4 run use TC position: (xln,ylt)
    if(mig.eq.1) then
     xlnc=poslon
     yltc=poslat
    else
     xlnc=xln
     yltc=ylt
    end if
   
    tavr=0.0
    counter=0.0
    taumax=0.0
    RRCT=1.e8
    irct=1000
    jrct=1000
    do j=1,jm
     do i=1,im
      x1=(LONGMIN+real(I-1)*(LONGMAX-LONGMIN)/real(IM-1))
      y1=(LATMIN+real(J-1)*(LATMAX-LATMIN)/real(JM-1))
      x0=xlnc*pi/180.
      y0=yltc*pi/180.
      DELTAX=REARTH*COS(y0)*(x1-x0)
      DELTAY=REARTH*(y1-y0)
      r=SQRT(DELTAX**2+DELTAY**2)
      if(r.lt.RAVR) then
        tauabs=sqrt(wusurf(i,j)**2+wvsurf(i,j)**2)
        if(tauabs*fsm(i,j).gt.taumax) taumax=tauabs
        tavr=tavr+tauabs*fsm(i,j)
        counter=counter+fsm(i,j)
      end if
      if(r.lt.RRCT) then
       RRCT=r
       irct=i
       jrct=j
       xrct=x1*180./pi
       yrct=y1*180./pi
      end if
     end do
    end do
    if(counter.gt.0.) then
      tauavr=tavr/counter
    else
      tauavr=0.0
    end if

    if(mig.eq.1.and.migtau.eq.0) then
!----falk 08-19-03 use taumax instead of tauavr
!    if(tauavr.gt.tauavrp) then
!      awucon=tauavrp/tauavr
!      bwucon=(tauavr-tauavrp)/tauavr
     if(taumax.gt.taumaxp) then
      awucon=taumaxp/taumax
      bwucon=(taumax-taumaxp)/taumax
     else
      awucon=1.
      bwucon=0.
     end if

      migtau=1
    end if

    if(mig.eq.1) then
      wucon=awucon+SIN(time*24./timesm*pi*0.5)*bwucon
      if(time*24..gt.timesm) wucon=1.
    else
      wucon=SIN(time*24./timesm4*pi*0.5)
      if(time*24..gt.timesm4) wucon=1.
    end if
   
    do j=1,jm
      do i=1,im
        wusurf(i,j)=wusurf(i,j)*wucon
        wvsurf(i,j)=wvsurf(i,j)*wucon
        ! taux(i,j)=taux(i,j)*wucon
        ! tauy(i,j)=tauy(i,j)*wucon
      end do
    end do

    return
  end subroutine avrtau
!
!------------------- falk 06-21-05 use new SST assimilation procedure
!

  subroutine hurricane_wind_date2day(year,julday,date)
    integer*4 date
    integer dat2day(12),dat2dayl(12),day,month,year,hour,n
    real julday
    real*8 tmp
    data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
    data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/
    
    year=int(date/1000000.)
    month=nint(100*(date/1000000.-int(date/1000000.)))
    julday=0
    if(mod(year,4).eq.0) then
      do n=1,month-1
        julday=julday+dat2dayl(n)
      end do
    else
      do n=1,month-1
        julday=julday+dat2day(n)
      end do
    end if
    
    julday=julday+nint(100*(date/10000.-int(date/10000.)))
    hour=date-nint(date/100.)*100
    julday=julday+real(hour)/24.
    
    return
  end subroutine hurricane_wind_date2day
  
	subroutine hurricane_wind_day2date(year,julday,date)
    integer*4 date
    integer dat2day(12),dat2dayl(12),day,month,year,year1,hour
    integer n
    
    real julday,julday1
    real*8 tmp
    data dat2day/31,28,31,30,31,30,31,31,30,31,30,31/
    data dat2dayl/31,29,31,30,31,30,31,31,30,31,30,31/
    
    if(int(julday).gt.365+int(1./(mod(year,4)*100.+1.))) then
      julday1=julday-365-int(1./(mod(year,4)*100.+1.))
      year1=year+1
    else
      julday1=julday
      year1=year
    end if
    day=0
    n=1
    if(mod(year1,4).eq.0) then  
      do while(day+dat2dayl(n).lt.int(julday1))
        day=day+dat2dayl(n)
        n=n+1   
      end do
    else
      do while(day+dat2day(n).lt.int(julday1))
        day=day+dat2day(n)
        n=n+1   
      end do
    end if
    
    month=n
    day=int(julday1-day)
    hour=nint((julday1-int(julday1))*24.)
    date=year1*1000000+month*10000+day*100+hour

    return
  end subroutine hurricane_wind_day2date
      
  function expwnd(R,RMAX,Rref18,Rref26,WSMAX,WS18,WS26)
    real Rref18,Rref26,R,RMAX,WSMAX,WS18,WS26,r1,ws
    real b,expwnd
    integer i,j
   
    r1=0.5*(Rref18+Rref26)
    WS=0.5*(WS18+WS26)
    if(Rref18.le.0.) then
      r1=Rref26
      WS=WS26
    end if
    
    if(Rref26.le.0.) then
      r1=Rref18
      WS=WS18
    end if
   
    if(R.GE.RMAX) then
      b=(RMAX-r1)/log(WS/WSMAX)                     
      expwnd=WSMAX*exp((RMAX-R)/b)
    else
      expwnd=R*WSMAX/RMAX
    end if
   
    return
  end function expwnd
  
  subroutine hurricane_wind_interp1d(mask,n,ni,x,y,xi,yi)
!--------------------------------------------------------------
!   This subroutine determines ni values yi at the points xi 
!   interpolating between the n values y at the points x
!   values equal to mask are ignored
!--------------------------------------------------------------
    real x(n),y(n),xi(ni),yi(ni),tmp(500)
    real cmp,mask
    integer n,ni,ii,i,j
   
    do i=1,ni
      if((xi(i).gt.x(n).and.xi(i).gt.x(1)).or.(xi(i).lt.x(1).and.xi(i).lt.x(n))) then
        if(xi(i).gt.x(n).and.xi(i).gt.x(1)) then
          if(x(n).gt.x(1)) then
            yi(i)=y(n)
          else
            yi(i)=y(1)
          end if
        else
          if(x(n).gt.x(1)) then
            yi(i)=y(1)
          else
            yi(i)=y(n)
          end if
        end if
      else
        do j=1,n-1
          tmp(j)=(xi(i)-x(j))*(xi(i)-x(j+1))
        end do
        do j=1,n-1
          if(tmp(j).le.0) ii=j
        end do
        if(y(ii).eq.mask.or.y(ii+1).eq.mask) then
          yi(i)=mask
        else
          yi(i)=(y(ii)*abs(x(ii+1)-xi(i))+y(ii+1)*abs(xi(i)-x(ii)))/abs(x(ii+1)-x(ii))
        end if
      end if
    end do

    return
  end subroutine hurricane_wind_interp1d
      
  SUBROUTINE hurricane_wind_OUTPUT(u_wind,v_wind)
!=========================================================
! save the output for graphics/analysis
! Written by Erxuan Fu, GSO,URI,11/3/94
!=========================================================
    use hurricane_wind_common
    
    real, intent(in) :: u_wind(:,:),v_wind(:,:)
    
    REAL TMP1(IM,JM),TMP2(IM,JM),TMP3(IM,JM,KB),jjulday
    REAL DB1(KB-1)
    real, DIMENSION(im,jm,kb) :: TB1,U1,V1,SB1
    real, DIMENSION(im,jm)    :: OHC,MLDTH,TB1S
    INTEGER yyear,outunitx,outunity
    CHARACTER FN*15, DOUT*8
    character*4 cjm
    integer*4 date
   
    call hurricane_wind_day2date(year,time,date)
    WRITE(DOUT,'(I8.8)') date

    FN = 'WINDX.'//DOUT
    OPEN(41,FILE=FN,STATUS='unknown')
    do i=1,2160
      do j=1,1080
        WRITE(41,231) u_wind(i,j)
      end do
    end do 
    CLOSE(41)

    FN = 'WINDY.'//DOUT    
    OPEN(41,FILE=FN,STATUS='unknown')
    do i=1,2160
      do j=1,1080
        WRITE(41,231) v_wind(i,j)
      end do
    end do 
    CLOSE(41)
231 format( F10.6)    
! 231 format(1080F10.6)

    RETURN
  end subroutine hurricane_wind_output
      
  subroutine hurricane_wind_verinterp(n,ni,x,y,xi,yi)
!--------------------------------------------------------------
!   This subroutine determines ni values yi at the points xi 
!   interpolating between the n values y at the points x
!--------------------------------------------------------------
    real x(n),y(n),xi(ni),yi(ni),tmp(500)
    real cmp
    integer n,ni,ii,i,j
    
    do i=1,ni
			ii=1
      if((xi(i).gt.x(n).and.xi(i).gt.x(1)).or.(xi(i).lt.x(1).and.xi(i).lt.x(n))) then
        if(xi(i).gt.x(n).and.xi(i).gt.x(1)) then
          if(x(n).gt.x(1)) then
            yi(i)=y(n)
          else
            yi(i)=y(1)
          end if
        else
          if(x(n).gt.x(1)) then
            yi(i)=y(1)
          else
            yi(i)=y(n)
          end if
        end if
      else
        do j=1,n-1
          tmp(j)=(xi(i)-x(j))*(xi(i)-x(j+1))
        end do
        do j=1,n-1
          if(tmp(j).le.0) then
            ii=j
          end if
        end do
        yi(i)=(y(ii)*abs(x(ii+1)-xi(i))+y(ii+1)*abs(xi(i)-x(ii)))/abs(x(ii+1)-x(ii))
      end if
    end do

    return
  end subroutine hurricane_wind_verinterp
end module hurricane_wind_mod
