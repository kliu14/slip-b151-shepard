    
module output 

!*----------------------------------------------------------------------------*! 
!*            output module for IMFM  channel code                            *! 
!*----------------------------------------------------------------------------*! 


  use global 
  use numerics 
  use numerics_tee

  implicit none 

  logical :: coarse_file
  integer , parameter :: int_size = 4 
  integer , parameter :: double_size = 8 
! integer  :: io_ver_number = 5004   ! 5001: save unecessary v_t information. 5002: save only current time step information. 
                                     ! 5003: save previous time step's
                                     ! information for Adam's bashforth method. 
                                     ! 5004: save f_old for full AB2 

contains 

!==================================================================================================! 
   subroutine setup_init_restart   ! under section 2. Read restart file or initialization 

if ( rstrt ) then !===================== if restart ======================
     call read_restart_file_interface ( 'output/restart/restart.in')
     if (rank .eq. 0 .and. log_flag) write(*,*) "2-1. read restart.in"
     call bcuvw
     call bcp ( p_old )
     if ( transportee ) then
       call read_restart_file_interface_tee ( 'output/restart/restart_transportee.in')
       if (rank .eq. 0 .and. log_flag) write(*,*) "2-1. read restart_transportee.in"
       if ( IMFMswitch ) then
          call read_restart_file_interface_IMFM ( 'output/restart/restart_IMFM.in')
          if (rank .eq. 0 .and. log_flag) write(*,*) "2-1. read restart_IMFM.in"
          if ( pat_bc .and. match_bc ) then 
            call sw_calc
            sb_old = sb_new
            st_old = st_new
          endif    
       endif
       call bcuvw_tee
       call bcuvw_source
       call bcp ( p_old_tee )
       call setup_BC_tee
     endif
     if (rstrt_AB .eq. .false.) then 
       flat_EE = .true.
     endif
  !------------- endif restart-------------!
else   !=================== Initialization ===================!
     if ( IMFMswitch ) then
       call read_restart_file_interface ( 'output/restart/restart.in')
       if (rank .eq. 0 .and. log_flag) write(*,*) "2-1. read restart.in"
       call bcuvw
       call bcp ( p_old )

       call init_field_tee
       if (rank .eq. 0 .and. log_flag) write(*,*) "2-1. Initialize the transportee field with given IC"
       call setup_BC_tee
     else
       if (re .le. 50) then
         call init_field_laminar
       else
         if ( slip_bc ) then
           call init_field_shs !(Kim: 08.31.22)
           !call read_restart_file_interface( 'output/restart/restart.in')
           if (rank .eq. 0 .and. log_flag) write(*,*) "2-1. Initialize the field with slip velocity "
         else
           call init_field
           if (rank .eq. 0 .and. log_flag) write(*,*) "2-1. Initialize the field with random perturbation "
         endif
       endif
       if ( transportee ) then
         call init_field_tee
         call setup_BC_tee
       endif
     endif
     flat_EE = .true.
endif
 !------------- endif initialization setup-------------!

   if ( rstrt_AB .and. io_ver_number_prev .ne. io_ver_number )    call terminate('two differnt io_ver_number, but tried to use AB restart')

   end subroutine setup_init_restart 

!==================================================================================================! 

  subroutine setup_output 

    !------
    ! checks status of output files and sets 
    ! noutput = 0 
    !------

    implicit none 
    logical             :: exst  
    character(len=20) :: re_str
    character(len=20) :: sigma_str
    character(len=20) :: dt_str

    if ( rank .eq. 0 .and. log_flag) then 

  if (p_solver_type .eq. 0 ) then
         if (p_order .eq. 2) then
      write(*,*) '1-2-1. Poisson:FFT Solver 2nd order'
         elseif (p_order .eq. 4) then
      write(*,*) '1-2-1. Poisson:FFT Solver 4th order '
         endif
  endif


   endif 


if (rank .eq. 0) then 

       inquire (file='diagnostics.dat', exist=exst)
       if ( exst ) then 
         open(unit=10004, file='diagnostics.dat', form='formatted', status='replace')
       else
          open(unit=10004, file='diagnostics.dat', form='formatted', status='new')
       endif

    endif

    !--- 
    ! initialize constants
    !---
    noutput = 0 
    nconst  = 0
    nslice  = 0
    iens    = 0 
    n_interface = 0

  end subroutine setup_output
!==========================================================================================================! 

  subroutine integer_time_stamp ( whole )

    implicit none
    character(len=6)       ::  tstr
    character(len=3)       ::  whole
    character(len=2)       ::  decml

    write(whole, '(i3.3)') int ( time )

  end subroutine integer_time_stamp

!==========================================================================================================! 

  subroutine qcriterion

    implicit none
    real, dimension(:,:,:), allocatable :: uc, vc, wc
   !real, dimension(is+1:ie+1,1:n2,ks+1:ke+1) :: qcrit
    real, dimension(6) :: sij
    real, dimension(9) :: du
    real :: pp, qq, phi, hx, hz
    integer :: i,j,k

    pi = acos(-1.0d0)
    sij = 0.d0
    du  = 0.d0
    hx = (chlx*pi ) / dble(n1)
    hz = (chlz*pi ) / dble(n3)


    allocate (uc(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp,ks-halo_z+1:ke+halo_z+1))
    allocate (vc(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp,ks-halo_z+1:ke+halo_z+1))
    allocate (wc(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp,ks-halo_z+1:ke+halo_z+1))


    do k= ks-2, ke+2
    do j= 1,n2
    do i= is-2, ie+2

       uc(i,j,k)= 0.5d0*( q(i_u,i,j,k) + q(i_u,i+1,j,k))
       vc(i,j,k) = 0.5d0*( q(i_v,i,j-1,k) + q(i_v,i,j,k))
       wc(i,j,k) = 0.5d0*( q(i_w,i,j,k) + q(i_w,i,j,k+1))

    enddo
    enddo
    enddo



    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1

       ! compute sij ... 
       ! XXXX should use 4th order approx later
       if ( j .eq. 1 ) then
          du(1) = ( uc(i,j+1,k)-uc(i,j,k))/ ( y(i_u,j+1)-y(i_u,j)) ! uy
          du(6) = ( wc(i,j+1,k)-wc(i,j,k))/ ( y(i_u,j+1)-y(i_u,j)) ! wy 
          du(7) = ( vc(i,j+1,k)-vc(i,j,k))/ ( y(i_u,j+1)-y(i_u,j)) ! vy 
       else if ( j .eq. n2 ) then
          du(1) = ( uc(i,j,k)-uc(i,j-1,k))/ ( y(i_u,j)-y(i_u,j-1)) ! uy 
          du(6) = ( wc(i,j,k)-wc(i,j-1,k))/ ( y(i_u,j)-y(i_u,j-1)) ! wy
          du(7) = ( vc(i,j,k)-vc(i,j-1,k))/ ( y(i_u,j)-y(i_u,j-1)) ! vy 
       else

          du(1) = ( uc(i,j+1,k) - uc(i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! uy
          du(6) = ( wc(i,j+1,k) - wc(i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! wy
          du(7) = ( vc(i,j+1,k) - vc(i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! vy
       endif

       du(2) =  ( uc(i,j,k+1) - uc(i,j,k-1))/2.d0/hz   ! uz
       du(3) =  ( vc(i+1,j,k) - vc(i-1,j,k))/2.d0/hx   ! vx
       du(4) =  ( vc(i,j,k+1) - vc(i,j,k-1))/2.d0/hz   ! vz
       du(5) =  ( wc(i+1,j,k) - wc(i-1,j,k))/2.d0/hx   ! wx 
       du(8) =  ( wc(i,j,k+1) - wc(i,j,k-1))/2.d0/hz   ! wz 
       du(9) =  ( uc(i+1,j,k) - uc(i-1,j,k))/2.d0/hx   ! ux 


       !------
       ! q = -0.5* u_{i,j}* u_{j,i}
       !------
       qcrit(i,j,k) = -0.5d0*( du(9)*du(9) + du(8)*du(8) + du(7)*du(7) + &
                                2.d0* ( du(1)*du(3) + du(6)*du(4) + du(2)*du(5)))
    enddo
    enddo
    enddo




    deallocate(uc)
    deallocate(vc)
    deallocate(wc)

  end subroutine qcriterion


!================================================================================================================================!

  subroutine dump_qcrit

    !-----------------------------------------------------------------
    ! collapse entire velocity and pressure fields onto the master node and output one file 
    ! currently assumes that all subdomains are the same size
    !-----------------------------------------------------------------

    implicit none
    real, dimension(:,:,:), allocatable :: qq
    real, dimension(:)        , allocatable :: recvbuf
    integer                                 :: i,j,k
    integer                                 :: qblock , rshft, n, proc, count, ind, pblock
    integer, dimension(:), allocatable      :: ks_vec, is_vec
    logical                                 :: debug_collapse
    logical                                 :: write_collapsed_grid
    real                                    :: hx, hz, xx, zz, u, v,w
    character(len=3)                        :: nstr
    character(len=3)                        :: time_str_i
    character(len=4)                        :: time_str_ibp


    !---------------------
    ! gather the start and ending indices for 
    ! all processors to master 
    ! currently assumes that all the blocks are of the same size
    !--------------------

    allocate ( ks_vec(0:size-1), is_vec(0:size-1))

    call mpi_gather( ks, 1, mpi_integer, ks_vec(0), 1, mpi_integer, 0, mycomm, ierr)
    call mpi_gather( is, 1, mpi_integer, is_vec(0), 1, mpi_integer, 0, mycomm, ierr)

    nx = int ( n1 / p1 )
    nz = int ( n3 / p3 )


    qblock  = (nx)* (nz)* (n2 )
    allocate ( recvbuf( size*qblock) )

    !--------------------------
    ! gather q from all processors to master
    !--------------------------
    rshft = rank * qblock + 1
    call mpi_gather ( qcrit(is+1,1,ks+1), qblock, MPI_DOUBLE_PRECISION, &
                      recvbuf(1), qblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)

    allocate ( qq(0:n1-1,1:n2, 0:n3-1))

    !--------------------------
    ! master node organizes the received data 
    !-------------------------

    if ( rank .eq. 0 ) then

       do proc=0,size-1

          !-----------------
          ! parse recv buff and extract values 
          ! into qq 
          !-----------------
          rshft = proc * qblock +1
          count = 0

          do k= ks_vec(proc)+1, ks_vec(proc)+nz-1+1
          do j= 1, n2
          do i= is_vec(proc)+1, is_vec(proc)+nx-1+1

             ind = rshft+count
             if ( k .ge. 0 .and. k .le. n3-1 .and. &
                  i .ge. 0 .and. i .le. n1-1 .and. &
                  j .ge. 1 .and. j .le. n2  ) then

                qq(i,j,k) = recvbuf(rshft+count)
             endif

             count = count + 1
          enddo
          enddo
          enddo

       enddo


       !---------------- 
       ! output collapsed velocities to file 
       !---------------- 
       write(time_str_ibp, '(i4.4)') int((time- int(time))*10000)
       write(time_str_i,'(i3.3)') int(time)
       open(unit=499,file='output/interface/qcrit'//time_str_i//'_'//time_str_ibp//'.dat',form='unformatted', status='unknown')
       write(499) n1,n2,n3
       write(499) y(i_u,:)
       write(499)  qq
       close(499)
   endif



    !--- 
    ! temporary cleanup 
    !---
    deallocate ( recvbuf)
    deallocate ( qq     )
    deallocate ( ks_vec  )
    deallocate ( is_vec  )

  end subroutine dump_qcrit

!=================================================================================================!
  subroutine dump_restart_file_interface (fname) 

    implicit none 
    character(*) :: fname 
    integer fhandle, finfo, filemode, fs_type , fs_type_intf
    integer ( kind=mpi_offset_kind) offset , loc_off
    integer ( kind=mpi_offset_kind) offset_intf , loc_off_intf
    integer status(mpi_status_size) 
    integer, dimension(3) :: gbl_size, loc_size, starts
    integer, dimension(3) :: gbl_size_intf, loc_size_intf, starts_intf 
    integer brick_count, brick_count_intf
    real, dimension(:,:,:), allocatable :: tmp_buf 
    real, dimension(:,:,:), allocatable :: tmp_buf_intf
    integer  :: i,j,k
    integer  :: all_count, all_count_intf
    logical :: exst

    if (rank .eq. 0 ) write(*,*) fname

    filemode = ior ( mpi_mode_wronly, mpi_mode_create) 

    call mpi_info_create ( finfo, ierr) 
    call mpi_file_open ( mycomm, fname, filemode, finfo, fhandle, ierr) 
    loc_off = 0 
    loc_off_intf = 0

    !---- 
    ! create the mpi data types for dumping q, p, vt, 
    !---- 
    
    gbl_size(1) = n1 
    gbl_size(2) = n2 
    gbl_size(3) = n3 

    loc_size(1) = nx 
    loc_size(2) = n2 
    loc_size(3) = nz 

    starts(1)   = is !is+1
    starts(2)   = 0  !1
    starts(3)   = ks !ks+1 

    gbl_size_intf(1) = n1
    gbl_size_intf(2) = 2
    gbl_size_intf(3) = n3

    loc_size_intf(1) = nx
    loc_size_intf(2) = 2
    loc_size_intf(3) = nz

    starts_intf(1) = is
    starts_intf(2) = 0
    starts_intf(3) = ks


    call mpi_type_create_subarray ( 3, gbl_size, loc_size, starts, & 
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type, ierr) 
    call mpi_type_commit ( fs_type, ierr) 

    call mpi_type_create_subarray ( 3, gbl_size_intf, loc_size_intf, starts_intf, &
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type_intf, ierr)
    call mpi_type_commit ( fs_type_intf, ierr)

    !----- 
    ! io version number ...
    !------
    if ( rank .eq. 0) call mpi_file_write ( fhandle, io_ver_number, 1, mpi_integer, status, ierr)  
    offset = 1* int_size 

    !---- 
    ! time...
    !---- 
    if ( rank .eq. 0) call mpi_file_write ( fhandle, time, 1, MPI_DOUBLE_PRECISION, status, ierr ) 
    offset = offset + 1* double_size 

    allocate ( tmp_buf_intf(0:nx-1, 0:1, 0:nz-1))
    allocate ( tmp_buf(0:nx-1, 0:n2-1, 0:nz-1))
 
    tmp_buf = 0.d0 
    tmp_buf_intf = 0.d0

    brick_count = nx * n2 * nz
    all_count = n1*n2*n3

    brick_count_intf = nx * 2 * nz
    all_count_intf = n1*2*n3

    !------
    ! u vel
    !------
    
    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 
       tmp_buf(i-is,j-1,k-ks) = q(i_u,i,j,k) 
    enddo
    enddo
    enddo 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 


    offset = offset + double_size* all_count


    !---- 
    ! v vel 
    !----
    tmp_buf = 0.d0 

    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 
          tmp_buf(i-is,j-1,k-ks) = q(i_v,i,j,k) 
    enddo
    enddo
    enddo 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    offset = offset + double_size* all_count


    !----- 
    ! w vel 
    !----- 
    tmp_buf = 0.d0 

    do k=ks,ke 
    do j=1,n2
    do i=is,ie 

       tmp_buf(i-is,j-1,k-ks) = q(i_w,i,j,k) 
    enddo
    enddo
    enddo 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf, brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    offset = offset + double_size* all_count


    !-----
    ! p_old 
    !-----
    tmp_buf = 0.d0 
    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1 

       tmp_buf(i-is-1,j-1,k-ks-1) = p_old(i,j,k) 
    enddo
    enddo
    enddo

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    offset = offset + double_size* all_count
    
if (io_ver_number .eq. 5001) then
    !----- 
    ! vt 
    !----- 
    tmp_buf = 0.d0 

    do k=ks,ke
    do j=1,n2 
    do i=is,ie 

       tmp_buf(i-is,j-1,k-ks) = vt(i,j,k) 
    enddo
    enddo 
    enddo 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    offset = offset + double_size* all_count
endif    


    offset_intf = offset


 if (io_ver_number .ge. 5003  .and. ab_flag) then

   offset = offset_intf
    !------
    ! u_old vel
    !------
    do k=ks,ke
    do j=1,n2
    do i=is,ie
       tmp_buf(i-is,j-1,k-ks) = f_old(i_u,i,j,k)
    enddo
    enddo
    enddo

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr)
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, &
                                 status, ierr )


    offset = offset + double_size* all_count

    !---- 
    ! v_old vel 
    !----
    tmp_buf = 0.d0

    do k=ks,ke
    do j=1,n2
    do i=is,ie
          tmp_buf(i-is,j-1,k-ks) = f_old(i_v,i,j,k)
    enddo
    enddo
    enddo

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr)
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count,MPI_DOUBLE_PRECISION, &
                                 status, ierr )

    offset = offset + double_size* all_count

    !----- 
    ! w_old vel 
    !----- 
    tmp_buf = 0.d0

    do k=ks,ke
    do j=1,n2
    do i=is,ie
       tmp_buf(i-is,j-1,k-ks) = f_old(i_w,i,j,k)
    enddo
    enddo
    enddo

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr)
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, &
                                 status, ierr )

    offset = offset + double_size* all_count
    offset_intf = offset
endif
    
    !--- cleanup 
    call mpi_file_close ( fhandle, ierr) 
    call mpi_type_free  (fs_type, ierr )  
    call mpi_type_free  (fs_type_intf, ierr )


    deallocate ( tmp_buf)
    deallocate ( tmp_buf_intf) 
    
  end subroutine dump_restart_file_interface

!=================================================================================================!

  subroutine dump_restart_file_interface_transportee (fname) 

    implicit none 
    character(*) :: fname 
    integer fhandle, finfo, filemode, fs_type , fs_type_intf
    integer ( kind=mpi_offset_kind) offset , loc_off
    integer ( kind=mpi_offset_kind) offset_intf , loc_off_intf
    integer status(mpi_status_size) 
    integer, dimension(3) :: gbl_size, loc_size, starts
    integer, dimension(3) :: gbl_size_intf, loc_size_intf, starts_intf 
    integer brick_count, brick_count_intf
    real, dimension(:,:,:), allocatable :: tmp_buf 
    real, dimension(:,:,:), allocatable :: tmp_buf_intf
    real, dimension(:), allocatable :: tmp_buf_IMFM
    real, dimension(:), allocatable :: tmp_buf_intf_IMFM
    integer  :: i,j,k
    integer  :: all_count, all_count_intf
    logical :: exst

    if (rank .eq. 0 ) write(*,*) fname

    filemode = ior ( mpi_mode_wronly, mpi_mode_create) 

    call mpi_info_create ( finfo, ierr) 
    call mpi_file_open ( mycomm, fname, filemode, finfo, fhandle, ierr) 
    loc_off = 0 
    loc_off_intf = 0

    !---- 
    ! create the mpi data types for dumping q, p, vt, transportee
    !---- 
    
    gbl_size(1) = n1 
    gbl_size(2) = n2 
    gbl_size(3) = n3 

    loc_size(1) = nx 
    loc_size(2) = n2 
    loc_size(3) = nz 

    starts(1)   = is !is+1
    starts(2)   = 0  !1
    starts(3)   = ks !ks+1 

    gbl_size_intf(1) = n1
    gbl_size_intf(2) = 2
    gbl_size_intf(3) = n3

    loc_size_intf(1) = nx
    loc_size_intf(2) = 2
    loc_size_intf(3) = nz

    starts_intf(1) = is
    starts_intf(2) = 0
    starts_intf(3) = ks


    call mpi_type_create_subarray ( 3, gbl_size, loc_size, starts, & 
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type, ierr) 
    call mpi_type_commit ( fs_type, ierr) 

    call mpi_type_create_subarray ( 3, gbl_size_intf, loc_size_intf, starts_intf, &
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type_intf, ierr)
    call mpi_type_commit ( fs_type_intf, ierr)

    !----- 
    ! io version number ...
    !------
    if ( rank .eq. 0) call mpi_file_write ( fhandle, io_ver_number, 1, mpi_integer, status, ierr)  
    offset = 1* int_size 

    !---- 
    ! time...
    !---- 
    if ( rank .eq. 0) call mpi_file_write ( fhandle, time, 1, MPI_DOUBLE_PRECISION, status, ierr ) 
    offset = offset + 1* double_size 

    allocate ( tmp_buf_intf(0:nx-1, 0:1, 0:nz-1))
    allocate ( tmp_buf(0:nx-1, 0:n2-1, 0:nz-1))
 
    tmp_buf = 0.d0 
    tmp_buf_intf = 0.d0

    brick_count = nx * n2 * nz
    all_count = n1*n2*n3

    brick_count_intf = nx * 2 * nz
    all_count_intf = n1*2*n3

    !------
    ! u vel
    !------
    
    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 
       tmp_buf(i-is,j-1,k-ks) = q_tee(i_u,i,j,k) 
    enddo
    enddo
    enddo 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 


    offset = offset + double_size* all_count


    !---- 
    ! v vel 
    !----
    tmp_buf = 0.d0 

    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 
          tmp_buf(i-is,j-1,k-ks) = q_tee(i_v,i,j,k) 
    enddo
    enddo
    enddo 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    offset = offset + double_size* all_count


    !----- 
    ! w vel 
    !----- 
    tmp_buf = 0.d0 

    do k=ks,ke 
    do j=1,n2
    do i=is,ie 

       tmp_buf(i-is,j-1,k-ks) = q_tee(i_w,i,j,k) 
    enddo
    enddo
    enddo 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf, brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    offset = offset + double_size* all_count


    !-----
    ! p_old 
    !-----
    tmp_buf = 0.d0 
    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1 

       tmp_buf(i-is-1,j-1,k-ks-1) = p_old_tee(i,j,k) 
    enddo
    enddo
    enddo

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_write_at_all ( fhandle, loc_off, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    offset = offset + double_size* all_count


! Danh 2018/6/13
! Currently, id is set for 5003
! Excessive modules have been erased

    !--- cleanup 
    call mpi_file_close ( fhandle, ierr) 
    call mpi_type_free  (fs_type, ierr )  
    call mpi_type_free  (fs_type_intf, ierr )


    deallocate ( tmp_buf)
    deallocate ( tmp_buf_intf) 
    
  end subroutine dump_restart_file_interface_transportee

!=================================================================================================!

  subroutine dump_restart_file_interface_IMFM (fname) 

    implicit none 
    character(*) :: fname 
    integer fhandle, finfo, filemode, fs_type , fs_type_intf
    integer ( kind=mpi_offset_kind) offset , loc_off
    integer ( kind=mpi_offset_kind) offset_intf , loc_off_intf
    integer status(mpi_status_size) 
    !integer, dimension(3) :: gbl_size, loc_size, starts
    !integer, dimension(3) :: gbl_size_intf, loc_size_intf, starts_intf 
    integer brick_count, brick_count_intf
    !real, dimension(:,:,:), allocatable :: tmp_buf 
    !real, dimension(:,:,:), allocatable :: tmp_buf_intf
    real, dimension(:), allocatable :: tmp_buf_IMFM
    real, dimension(:), allocatable :: tmp_buf_intf_IMFM
    integer  :: i,j,k
    integer  :: all_count, all_count_intf
    logical :: exst

    if (rank .eq. 0 ) write(*,*) fname

    filemode = ior ( mpi_mode_wronly, mpi_mode_create) 

    call mpi_info_create ( finfo, ierr) 
    call mpi_file_open ( mycomm, fname, filemode, finfo, fhandle, ierr) 
    loc_off = 0 
    loc_off_intf = 0


    !----- 
    ! io version number ...
    !------
    if ( rank .eq. 0) call mpi_file_write ( fhandle, io_ver_number, 1, mpi_integer, status, ierr)  
    offset = 1* int_size 

    !---- 
    ! time...
    !---- 
    if ( rank .eq. 0) call mpi_file_write ( fhandle, time, 1, MPI_DOUBLE_PRECISION, status, ierr ) 
    offset = offset + 1* double_size 
    
    !-----
    ! s_tee spanwise direction
    !-----
    !if ( rank .eq. 0) then
    if ( rank .eq. 0) then

       allocate ( tmp_buf_intf_IMFM(0:1))
       allocate ( tmp_buf_IMFM(0:n2-1))
 
       tmp_buf_IMFM = 0.d0 
       tmp_buf_intf_IMFM = 0.d0

       brick_count = n2
       all_count = n2

       brick_count_intf = 2
       all_count_intf = 2
    
       
       do j=1,n2
         
         tmp_buf_IMFM(j-1) = s_tee(i_u,j) 
       
       enddo

       call mpi_file_write ( fhandle, tmp_buf_IMFM(0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

       offset = offset + double_size* all_count

    !-----
    ! s_tee wall normal direction
    !-----

       tmp_buf_IMFM = 0.d0 
       tmp_buf_intf_IMFM = 0.d0
       
       do j=1,n2
         
         tmp_buf_IMFM(j-1) = s_tee(i_v,j) 
       
       enddo

       call mpi_file_write ( fhandle, tmp_buf_IMFM(0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

       offset = offset + double_size* all_count
    
        
    !-----
    ! s_tee wall parallel direction
    !-----

       tmp_buf_IMFM = 0.d0 
       tmp_buf_intf_IMFM = 0.d0
       
       do j=1,n2
         
         tmp_buf_IMFM(j-1) = s_tee(i_w,j) 
       
       enddo

       call mpi_file_write ( fhandle, tmp_buf_IMFM(0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

       offset = offset + double_size* all_count

       deallocate ( tmp_buf_intf_IMFM )
       deallocate ( tmp_buf_IMFM )
    
     endif

    !--- cleanup 
    call mpi_file_close ( fhandle, ierr) 
    
  end subroutine dump_restart_file_interface_IMFM

!=================================================================================================!

  subroutine read_restart_file_coarsen ( fname ) 

    implicit none
    character(*) :: fname
    integer fhandle, finfo, filemode, fs_type
    integer ( kind=mpi_offset_kind) offset , loc_off
    integer status(mpi_status_size)
    integer, dimension(3) :: gbl_size, loc_size, starts
    integer brick_count
    real, dimension(:,:,:), allocatable :: tmp_buf
    integer  :: i,j,k, isb, ksb, ieb, keb, buf_int
    integer  :: all_count

    filemode = mpi_mode_rdonly 

    call mpi_info_create ( finfo, ierr) 
    call mpi_file_open ( mycomm, fname, filemode, finfo, fhandle, ierr) 
    loc_off = 0 

    !---- 
    ! create the mpi data types for dumping q, p, vt, 
    !---- 
   
    gbl_size(1) = n1*2
    gbl_size(2) = n2 
    gbl_size(3) = n3*2

    loc_size(1) = nx*2 
    loc_size(2) = n2 
    loc_size(3) = nz*2 

    starts(1)   = is*2 !is+1
    starts(2)   = 0 !1
    starts(3)   = ks*2 !ks+1 

    call mpi_type_create_subarray ( 3, gbl_size, loc_size, starts, & 
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type, ierr) 
    call mpi_type_commit ( fs_type, ierr) 


    !----- 
    ! io version number ...
    !------
    call mpi_file_read_all ( fhandle, buf_int, 1, mpi_integer, status, ierr)  
    offset = 1* int_size 
    io_ver_number_prev = buf_int
    if ( rank .eq. 0. .and. buf_int .eq. 5001 ) then
       write(*,*) 'io version = 5001, it saved vt '
    elseif ( rank .eq. 0. .and. buf_int .eq. 5002 ) then
       write(*,*) 'io version = 5002, no saving of f_old'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5003 ) then
       write(*,*) 'io version = 5003, Saved f_old and interface information'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5004 ) then
       write(*,*) 'io version = 5004, Saved f_old and interface information'
    elseif ( rank .eq. 0. .and. buf_int .ge. 5005. .or. buf_int .le. 5000) then
       write(*,*) ' Wrong io version'
    endif

    !---- 
    ! time...
    !---- 
    call mpi_file_read_all ( fhandle, time, 1, MPI_DOUBLE_PRECISION, status, ierr ) 
    offset = offset + 1* double_size 

    allocate ( tmp_buf(0:nx*2-1, 0:n2-1, 0:nz*2-1))

    tmp_buf = 0.d0

    brick_count = (nx*2)*n2*(nz*2)
    all_count = (n1*2)*n2*(n3*2)

    !------
    ! u vel
    !------
    tmp_buf = 0.d0 
    
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks,ke
    do j=1,n2 
    do i=is,ie
       q(i_u,i,j,k) = tmp_buf((i-is)*2,j-1,(k-ks)*2) 
    enddo
    enddo
    enddo 


    offset = offset + double_size* all_count

    !---- 
    ! v vel 
    !----
    tmp_buf = 0.d0 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    do k=ks,ke
    do j=1,n2 
    do i=is,ie
       q(i_v,i,j,k) = tmp_buf((i-is)*2,j-1,(k-ks)*2)
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count

    !----- 
    ! w vel 
    !----- 
    tmp_buf = 0.d0 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    do k=ks,ke
    do j=1,n2
    do i=is,ie
       q(i_w,i,j,k) = tmp_buf((i-is)*2,j-1,(k-ks)*2)
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count


    !-----
    ! p_old 
    !-----
    tmp_buf = 0.d0 
   

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    do k=ks+1,ke+1-1
    do j=1,n2
    do i=is+1,ie+1-1

       p_old(i,j,k) = tmp_buf((i-is-1)*2,j-1,(k-ks-1)*2)
    enddo
    enddo
    enddo

    offset = offset + double_size* all_count

    !--- cleanup 
    call mpi_file_close ( fhandle, ierr) 
    call mpi_type_free  (fs_type, ierr )  
   
 
    deallocate ( tmp_buf) 

  end subroutine read_restart_file_coarsen

!=================================================================================================!

  subroutine read_restart_file_grid_refined ( fname ) 

    implicit none
    character(*) :: fname
    integer fhandle, finfo, filemode, fs_type
    integer ( kind=mpi_offset_kind) offset , loc_off
    integer status(mpi_status_size)
    integer, dimension(3) :: gbl_size, loc_size, starts
    integer brick_count
    real, dimension(:,:,:), allocatable :: tmp_buf
    integer  :: i,j,k, isb, ksb, ieb, keb, buf_int
    integer  :: all_count

    filemode = mpi_mode_rdonly 

    call mpi_info_create ( finfo, ierr) 
    call mpi_file_open ( mycomm, fname, filemode, finfo, fhandle, ierr) 
    loc_off = 0 

    !---- 
    ! create the mpi data types for dumping q, p, vt, 
    !---- 
   
    gbl_size(1) = n1/2
    gbl_size(2) = n2 
    gbl_size(3) = n3/2

    loc_size(1) = nx/2 
    loc_size(2) = n2 
    loc_size(3) = nz/2 

    starts(1)   = is/2 !is+1
    starts(2)   = 0 !1
    starts(3)   = ks/2 !ks+1 

    call mpi_type_create_subarray ( 3, gbl_size, loc_size, starts, & 
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type, ierr) 
    call mpi_type_commit ( fs_type, ierr) 


    !----- 
    ! io version number ...
    !------
    call mpi_file_read_all ( fhandle, buf_int, 1, mpi_integer, status, ierr)  
    offset = 1* int_size 
    io_ver_number_prev = buf_int
    if ( rank .eq. 0. .and. buf_int .eq. 5001 ) then
       write(*,*) 'io version = 5001, it saved vt '
    elseif ( rank .eq. 0. .and. buf_int .eq. 5002 ) then
       write(*,*) 'io version = 5002, no saving of f_old'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5003 ) then
       write(*,*) 'io version = 5003, Saved f_old and interface information'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5004 ) then
       write(*,*) 'io version = 5004, Saved f_old and interface information'
    elseif ( rank .eq. 0. .and. buf_int .ge. 5005. .or. buf_int .le. 5000) then
       write(*,*) ' Wrong io version'
    endif

    !---- 
    ! time...
    !---- 
    call mpi_file_read_all ( fhandle, time, 1, MPI_DOUBLE_PRECISION, status, ierr ) 
    offset = offset + 1* double_size 

    allocate ( tmp_buf(0:nx/2-1, 0:n2-1, 0:nz/2-1))

    tmp_buf = 0.d0

    brick_count = (nx/2)*n2*(nz/2)
    all_count = (n1/2)*n2*(n3/2)

    !------
    ! u vel
    !------
    tmp_buf = 0.d0 
    
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks,ke-1, 2
    do j=1,n2 
    do i=is,ie-1, 2
       q(i_u,i,j,k) = tmp_buf((i-is)/2,j-1,(k-ks)/2) 
       q(i_u,i+1,j,k) = tmp_buf((i-is)/2,j-1,(k-ks)/2) 
       q(i_u,i,j,k+1) = tmp_buf((i-is)/2,j-1,(k-ks)/2) 
       q(i_u,i+1,j,k+1) = tmp_buf((i-is)/2,j-1,(k-ks)/2) 
    enddo
    enddo
    enddo 


    offset = offset + double_size* all_count

    !---- 
    ! v vel 
    !----
    tmp_buf = 0.d0 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    do k=ks,ke-1,2
    do j=1,n2 
    do i=is,ie-1,2
       q(i_v,i,j,k) = tmp_buf((i-is)/2,j-1,(k-ks)/2)
       q(i_v,i+1,j,k) = tmp_buf((i-is)/2,j-1,(k-ks)/2)
       q(i_v,i,j,k+1) = tmp_buf((i-is)/2,j-1,(k-ks)/2)
       q(i_v,i+1,j,k+1) = tmp_buf((i-is)/2,j-1,(k-ks)/2)
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count

    !----- 
    ! w vel 
    !----- 
    tmp_buf = 0.d0 

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    do k=ks,ke-1,2
    do j=1,n2
    do i=is,ie-1,2
       q(i_w,i,j,k) = tmp_buf((i-is)/2,j-1,(k-ks)/2)
       q(i_w,i+1,j,k) = tmp_buf((i-is)/2,j-1,(k-ks)/2)
       q(i_w,i,j,k+1) = tmp_buf((i-is)/2,j-1,(k-ks)/2)
       q(i_w,i+1,j,k+1) = tmp_buf((i-is)/2,j-1,(k-ks)/2)
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count


    !-----
    ! p_old 
    !-----
    tmp_buf = 0.d0 
   

    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, & 
                                 status, ierr ) 

    do k=ks+1,ke+1-1,2
    do j=1,n2
    do i=is+1,ie+1-1,2

       p_old(i,j,k) = tmp_buf((i-is-1)/2,j-1,(k-ks-1)/2)
       p_old(i+1,j,k) = tmp_buf((i-is-1)/2,j-1,(k-ks-1)/2)
       p_old(i,j,k+1) = tmp_buf((i-is-1)/2,j-1,(k-ks-1)/2)
       p_old(i+1,j,k+1) = tmp_buf((i-is-1)/2,j-1,(k-ks-1)/2)
    enddo
    enddo
    enddo

    offset = offset + double_size* all_count

    !--- cleanup 
    call mpi_file_close ( fhandle, ierr) 
    call mpi_type_free  (fs_type, ierr )  
   
 
    deallocate ( tmp_buf) 

  end subroutine read_restart_file_grid_refined

!================================================================

  subroutine read_restart_file_interface ( fname ) 

    implicit none
    character(*) :: fname
    integer fhandle, finfo, filemode, fs_type , fs_type_intf
    integer ( kind=mpi_offset_kind) offset , loc_off
    integer ( kind=mpi_offset_kind) offset_intf , loc_off_intf
    integer status(mpi_status_size)
    integer, dimension(3) :: gbl_size, loc_size, starts
    integer, dimension(3) :: gbl_size_intf, loc_size_intf, starts_intf
    integer brick_count,brick_count_intf
    real, dimension(:,:,:), allocatable :: tmp_buf
    real, dimension(:,:,:), allocatable :: tmp_buf_intf
    integer  :: i,j,k, buf_int
    integer  :: all_count, all_count_intf
    logical :: exst

    if (rank .eq. 0 ) write(*,*) fname 
    filemode = mpi_mode_rdonly 

    call mpi_info_create ( finfo, ierr) 
    call mpi_file_open ( mycomm, fname, filemode, finfo, fhandle, ierr) 
    loc_off = 0 
    loc_off_intf = 0

    !---- 
    ! create the mpi data types for dumping q, p, vt, 
    !---- 
    gbl_size_intf(1) = n1
    gbl_size_intf(2) = 2
    gbl_size_intf(3) = n3

    loc_size_intf(1) = nx
    loc_size_intf(2) = 2
    loc_size_intf(3) = nz

    starts_intf(1) = is
    starts_intf(2) = 0
    starts_intf(3) = ks


    call mpi_type_create_subarray ( 3, gbl_size_intf, loc_size_intf, starts_intf, &
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type_intf, ierr)
    call mpi_type_commit ( fs_type_intf, ierr)
   
    gbl_size(1) = n1 
    gbl_size(2) = n2 
    gbl_size(3) = n3 

    loc_size(1) = nx 
    loc_size(2) = n2 
    loc_size(3) = nz 

    starts(1)   = is !is+1
    starts(2)   = 0 !1
    starts(3)   = ks !ks+1 

    call mpi_type_create_subarray ( 3, gbl_size, loc_size, starts, & 
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type, ierr) 
    call mpi_type_commit ( fs_type, ierr) 


    !----- 
    ! io version number ...
    !------
    call mpi_file_read_all ( fhandle, buf_int, 1, mpi_integer, status, ierr)  
    io_ver_number_prev = buf_int

    if ( rank .eq. 0. .and. buf_int .eq. 5001 ) then
       write(*,*) 'io version = 5001, it saved vt '
    elseif ( rank .eq. 0. .and. buf_int .eq. 5002 ) then
       write(*,*) 'io version = 5002, no saving of f_old'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5003 ) then
       write(*,*) 'io version = 5003, Saved f_old and interface information, AB/CN2'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5004 ) then
       write(*,*) 'io version = 5004, Saved f_old and interface information, AB2'
    elseif ( rank .eq. 0 .and. buf_int .le. 5000 .or. buf_int .ge. 5005  ) then
       write(*,*) ' Wrong io version, io_version=',buf_int
    endif

    offset = 1* int_size 

    !---- 
    ! time...
    !---- 
    call mpi_file_read_all ( fhandle, time, 1, MPI_DOUBLE_PRECISION, status, ierr ) 
    offset = offset + 1* double_size 

    allocate ( tmp_buf_intf(0:nx-1, 0:1, 0:nz-1))
    allocate ( tmp_buf(0:nx-1, 0:n2-1, 0:nz-1))

    tmp_buf = 0.d0
    tmp_buf_intf = 0.d0

    brick_count = nx * n2 * nz
    all_count = n1*n2*n3

    brick_count_intf = nx * 2 * nz
    all_count_intf = n1*2*n3

    
    !------
    ! u vel
    !------
    tmp_buf = 0.d0 
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 
       q(i_u,i,j,k) = tmp_buf(i-is,j-1,k-ks) 
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count

    !---- 
    ! v vel 
    !----
    tmp_buf = 0.d0 
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 
       q(i_v,i,j,k) = tmp_buf(i-is,j-1,k-ks) 
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count

    !----- 
    ! w vel 
    !----- 
    tmp_buf = 0.d0 
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks,ke 
    do j=1,n2
    do i=is,ie 
       q(i_w,i,j,k) = tmp_buf(i-is,j-1,k-ks) 
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count

    !-----
    ! p_old 
    !-----
    tmp_buf = 0.d0 
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1 
       p_old(i,j,k) = tmp_buf(i-is-1,j-1,k-ks-1) 
    enddo
    enddo
    enddo

    offset = offset + double_size* all_count

    if (rank .eq. 0 .and. log_flag) write(*,*) ' Read velocity and pressure '

if (buf_int .eq. 5001) then

    tmp_buf = 0.d0 
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks,ke
    do j=1,n2 
    do i=is,ie 

       vt(i,j,k) = tmp_buf(i-is,j-1,k-ks)
    enddo
    enddo 
    enddo 


    offset = offset + double_size* all_count
endif

    offset_intf = offset


 
 if (buf_int .eq. io_ver_number .and. rstrt_AB) then
    !-----
    ! f_old
    !-----
    offset= offset_intf
    tmp_buf = 0.d0
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type,'native', mpi_info_null, ierr)
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count,MPI_DOUBLE_PRECISION, status, ierr )

    do k=ks,ke
    do j=1,n2
    do i=is,ie
       f_old(i_u,i,j,k) = tmp_buf(i-is,j-1,k-ks)
    enddo
    enddo
    enddo

    offset = offset + double_size* all_count

    tmp_buf = 0.d0
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr)
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count,MPI_DOUBLE_PRECISION, status, ierr )

    do k=ks,ke
    do j=1,n2
    do i=is,ie
       f_old(i_v,i,j,k) = tmp_buf(i-is,j-1,k-ks)
    enddo
    enddo
    enddo

    offset = offset + double_size* all_count

    tmp_buf = 0.d0
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr)
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr )

    do k=ks,ke
    do j=1,n2
    do i=is,ie
       f_old(i_w,i,j,k) = tmp_buf(i-is,j-1,k-ks)
    enddo
    enddo
    enddo

    offset = offset + double_size* all_count

    if (rank .eq. 0 .and. log_flag) write(*,*) ' Read f_old for AB timestep, io_ver = ', io_ver_number

  endif


    !--- cleanup 
    call mpi_file_close ( fhandle, ierr) 
    call mpi_type_free  (fs_type, ierr )  
    call mpi_type_free  (fs_type_intf, ierr ) 
   
 
    deallocate ( tmp_buf) 
    deallocate ( tmp_buf_intf )  

  end subroutine read_restart_file_interface

!=================================================================================================!

  subroutine read_restart_file_interface_tee ( fname ) 

    implicit none
    character(*) :: fname
    integer fhandle, finfo, filemode, fs_type , fs_type_intf
    integer ( kind=mpi_offset_kind) offset , loc_off
    integer ( kind=mpi_offset_kind) offset_intf , loc_off_intf
    integer status(mpi_status_size)
    integer, dimension(3) :: gbl_size, loc_size, starts
    integer, dimension(3) :: gbl_size_intf, loc_size_intf, starts_intf
    integer brick_count,brick_count_intf
    real, dimension(:,:,:), allocatable :: tmp_buf
    real, dimension(:,:,:), allocatable :: tmp_buf_intf
    integer  :: i,j,k, buf_int
    integer  :: all_count, all_count_intf
    logical :: exst

    if (rank .eq. 0 ) write(*,*) fname 
    filemode = mpi_mode_rdonly 

    call mpi_info_create ( finfo, ierr) 
    call mpi_file_open ( mycomm, fname, filemode, finfo, fhandle, ierr) 
    loc_off = 0 
    loc_off_intf = 0

    !---- 
    ! create the mpi data types for dumping q_tee, p_tee, vt, 
    !---- 
    gbl_size_intf(1) = n1
    gbl_size_intf(2) = 2
    gbl_size_intf(3) = n3

    loc_size_intf(1) = nx
    loc_size_intf(2) = 2
    loc_size_intf(3) = nz

    starts_intf(1) = is
    starts_intf(2) = 0
    starts_intf(3) = ks


    call mpi_type_create_subarray ( 3, gbl_size_intf, loc_size_intf, starts_intf, &
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type_intf, ierr)
    call mpi_type_commit ( fs_type_intf, ierr)
   
    gbl_size(1) = n1 
    gbl_size(2) = n2 
    gbl_size(3) = n3 

    loc_size(1) = nx 
    loc_size(2) = n2 
    loc_size(3) = nz 

    starts(1)   = is !is+1
    starts(2)   = 0 !1
    starts(3)   = ks !ks+1 

    call mpi_type_create_subarray ( 3, gbl_size, loc_size, starts, & 
                                    mpi_order_fortran, MPI_DOUBLE_PRECISION, fs_type, ierr) 
    call mpi_type_commit ( fs_type, ierr) 


    !----- 
    ! io version number ...
    !------
    call mpi_file_read_all ( fhandle, buf_int, 1, mpi_integer, status, ierr)  
    io_ver_number_prev = buf_int

    if ( rank .eq. 0. .and. buf_int .eq. 5001 ) then
       write(*,*) 'io version = 5001, it saved vt for transportee variable'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5002 ) then
       write(*,*) 'io version = 5002, no saving of f_old for transportee variable'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5003 ) then
       write(*,*) 'io version = 5003, Saved f_old and interface information, AB/CN2 for transportee variable'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5004 ) then
       write(*,*) 'io version = 5004, Saved f_old and interface information, AB2 for transportee variable'
    elseif ( rank .eq. 0 .and. buf_int .le. 5000 .or. buf_int .ge. 5005  ) then
       write(*,*) ' Wrong io version, io_version=',buf_int
    endif

    offset = 1* int_size 

    !---- 
    ! time...
    !---- 
    call mpi_file_read_all ( fhandle, time, 1, MPI_DOUBLE_PRECISION, status, ierr ) 
    offset = offset + 1* double_size 

    allocate ( tmp_buf_intf(0:nx-1, 0:1, 0:nz-1))
    allocate ( tmp_buf(0:nx-1, 0:n2-1, 0:nz-1))

    tmp_buf = 0.d0
    tmp_buf_intf = 0.d0

    brick_count = nx * n2 * nz
    all_count = n1*n2*n3

    brick_count_intf = nx * 2 * nz
    all_count_intf = n1*2*n3

    
    !------
    ! u vel
    !------
    tmp_buf = 0.d0 
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 
       q_tee(i_u,i,j,k) = tmp_buf(i-is,j-1,k-ks) 
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count

    !---- 
    ! v vel 
    !----
    tmp_buf = 0.d0 
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 
       q_tee(i_v,i,j,k) = tmp_buf(i-is,j-1,k-ks) 
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count

    !----- 
    ! w vel 
    !----- 
    tmp_buf = 0.d0 
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks,ke 
    do j=1,n2
    do i=is,ie 
       q_tee(i_w,i,j,k) = tmp_buf(i-is,j-1,k-ks) 
    enddo
    enddo
    enddo 

    offset = offset + double_size* all_count

    !-----
    ! p_old 
    !-----
    tmp_buf = 0.d0 
    call mpi_file_set_view ( fhandle, offset, MPI_DOUBLE_PRECISION, fs_type, 'native', mpi_info_null, ierr) 
    call mpi_file_read ( fhandle, tmp_buf(0,0,0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1 
       p_old_tee(i,j,k) = tmp_buf(i-is-1,j-1,k-ks-1) 
    enddo
    enddo
    enddo

    offset = offset + double_size* all_count

    if (rank .eq. 0 .and. log_flag) write(*,*) ' Read velocity and pressure for transportee variable'


    !--- cleanup 
    call mpi_file_close ( fhandle, ierr) 
    call mpi_type_free  (fs_type, ierr )  
    call mpi_type_free  (fs_type_intf, ierr ) 
   
 
    deallocate ( tmp_buf) 
    deallocate ( tmp_buf_intf )  

  end subroutine read_restart_file_interface_tee

!=================================================================================================!

  subroutine read_restart_file_interface_IMFM ( fname ) 

    implicit none
    character(*) :: fname
    integer fhandle, finfo, filemode, fs_type , fs_type_intf
    integer ( kind=mpi_offset_kind) offset , loc_off
    integer ( kind=mpi_offset_kind) offset_intf , loc_off_intf
    integer status(mpi_status_size)
    integer brick_count,brick_count_intf
    real, dimension(:), allocatable :: tmp_buf_IMFM
    real, dimension(:), allocatable :: tmp_buf_intf_IMFM
    integer  :: i,j,k, buf_int
    integer  :: all_count, all_count_intf
    logical :: exst

    if (rank .eq. 0 ) write(*,*) fname 
    filemode = mpi_mode_rdonly 

    call mpi_info_create ( finfo, ierr) 
    call mpi_file_open ( mycomm, fname, filemode, finfo, fhandle, ierr) 
    loc_off = 0 
    loc_off_intf = 0

    !----- 
    ! io version number ...
    !------
    call mpi_file_read_all ( fhandle, buf_int, 1, mpi_integer, status, ierr)  
    io_ver_number_prev = buf_int

    if ( rank .eq. 0. .and. buf_int .eq. 5001 ) then
       write(*,*) 'io version = 5001, it saved vt for IMFM variable'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5002 ) then
       write(*,*) 'io version = 5002, no saving of f_old for IMFM variable'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5003 ) then
       write(*,*) 'io version = 5003, Saved f_old and interface information, AB/CN2 for IMFM variable'
    elseif ( rank .eq. 0. .and. buf_int .eq. 5004 ) then
       write(*,*) 'io version = 5004, Saved f_old and interface information, AB2 for IMFM variable'
    elseif ( rank .eq. 0 .and. buf_int .le. 5000 .or. buf_int .ge. 5005  ) then
       write(*,*) ' Wrong io version, io_version=',buf_int
    endif

    offset = 1* int_size 

    !---- 
    ! time...
    !---- 
    call mpi_file_read_all ( fhandle, time, 1, MPI_DOUBLE_PRECISION, status, ierr ) 
    offset = offset + 1* double_size 

    !---- 
    ! source term in spanwise direction
    !---- 

    allocate ( tmp_buf_intf_IMFM(0:1))
    allocate ( tmp_buf_IMFM(0:n2-1))

    tmp_buf_IMFM = 0.d0
    tmp_buf_intf_IMFM = 0.d0
    
    brick_count = n2
    all_count = n2

    brick_count_intf = 2
    all_count_intf = 2

    call mpi_file_read_all ( fhandle, tmp_buf_IMFM(0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do j=1,n2
       s_tee(i_u,j) = tmp_buf_IMFM(j-1) 
    enddo

    offset = offset + double_size* all_count

    !---- 
    ! source term in wall normal direction
    !---- 

    tmp_buf_IMFM = 0.d0

    call mpi_file_read_all ( fhandle, tmp_buf_IMFM(0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do j=1,n2
       s_tee(i_v,j) = tmp_buf_IMFM(j-1) 
    enddo

    offset = offset + double_size* all_count

    !---- 
    ! source term in wall parallel direction
    !---- 

    tmp_buf_IMFM = 0.d0

    call mpi_file_read_all ( fhandle, tmp_buf_IMFM(0), brick_count, MPI_DOUBLE_PRECISION, status, ierr ) 

    do j=1,n2
       s_tee(i_w,j) = tmp_buf_IMFM(j-1) 
    enddo

    if (rank .eq. 0 .and. log_flag) write(*,*) 'Read IMFM variable'


    !--- cleanup 
    call mpi_file_close ( fhandle, ierr) 
 
    deallocate ( tmp_buf_IMFM ) 
    deallocate ( tmp_buf_intf_IMFM )  

  end subroutine read_restart_file_interface_IMFM

!=================================================================================================!

  subroutine diag_perturb

!.... output time , ut_u, ut_l , div max , and cfl_max 
!.... at each time step to diagnostics.dat 

    use numerics 
    implicit none
    real                                 ::       mdiv 
    real                                 ::       mcfl, cfl
    real                                 ::       utl , utu , rhyl , rhyu
    integer                              ::       i,j,k, nout
    real                                 ::       tcoeff    , gmcfl, gmdiv, gutl, gutu


    rhx = dble(n1) / (chlx*pi)
    rhz = dble(n3) / (chlz*pi) 

!.... mean u at 1,n2 

    utl = 0.d0 
    utu = 0.d0 

    do k=ks,ke
    do i=is,ie

       utl = utl + q(i_u,i, 1,k) - (q(i_u,i,1,k) +  q(i_u,i,0,k))/2.d0
       utu = utu + q(i_u,i,n2,k) - (q(i_u,i,n2+1,k) + q(i_u,i,n2,k))/2.d0

    enddo
    enddo 

    
    call mpi_reduce ( utl, gutl, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm , ierr) 
    call mpi_reduce ( utu, gutu, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm , ierr) 

    if ( rank .eq. 0 ) then 
       gutl = gutl / dble(n1)/dble(n3) 
       gutu = gutu / dble(n1)/dble(n3) 

       !.... ut = sqrt ( 2*ut * rhy * rre ) 
       rhyl= 1.d0 / ( y(i_v,1) - y(i_v,0))
       rhyu= 1.d0 / ( y(i_v,n2)- y(i_v,n2-1))

       gutl = sqrt ( 2.d0 * gutl * rhyl / re ) 
       gutu = sqrt ( 2.d0 * gutu * rhyu / re ) 

       
    endif 


      if (gutl .lt. 1.0d0 .and. gutu .lt. 1.0d0) then

     call add_perturb_to_field
   
    ! add endif - unterminated
    endif

  end subroutine diag_perturb

!==================================================================================================================!
    

  subroutine diag_out 

!.... output time , ut_u, ut_l , div max , and cfl_max 
!.... at each time step to diagnostics.dat 

    use numerics 
    implicit none
    real                                 ::       mdiv 
    real                                 ::       mcfl, cfl
    real                                 ::       utl , utu , rhyl , rhyu
    real                                 ::       utl_tee , utu_tee
    real, dimension(:,:,:),allocatable   ::       div 
    integer                              ::       i,j,k, nout
    real                                 ::       tcoeff    , gmcfl, gmdiv, gutl, gutu, gutl_tee, gutu_tee

    allocate ( div(is+1:ie+1,1:n2,ks+1:ke+1) ) 

    rhx = dble(n1) / (chlx*pi)
    rhz = dble(n3) / (chlz*pi) 

!.... mean u at 1,n2 

    utl = 0.d0 
    utu = 0.d0 

    utl_tee = 0.d0 
    utu_tee = 0.d0 

    do k=ks,ke
    do i=is,ie

       utl = utl + q(i_u,i, 1,k) - (q(i_u,i,1,k) +  q(i_u,i,0,k))/2.d0
       utu = utu + q(i_u,i,n2,k) - (q(i_u,i,n2+1,k) + q(i_u,i,n2,k))/2.d0
       if (transportee) then
          utl_tee = utl_tee + q_tee(i_u,i, 1,k)!- (q_tee(i_u,i,1,k) +  q_tee(i_u,i,0,k))/2.d0
          utu_tee = utu_tee + q_tee(i_u,i,n2,k)!- (q_tee(i_u,i,n2+1,k) + q_tee(i_u,i,n2,k))/2.d0
       endif
    enddo
    enddo 

    
    call mpi_reduce ( utl, gutl, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm , ierr) 
    call mpi_reduce ( utu, gutu, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm , ierr) 
    if (transportee) then
       call mpi_reduce ( utl_tee, gutl_tee, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm , ierr) 
       call mpi_reduce ( utu_tee, gutu_tee, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm , ierr) 
    endif

    if ( rank .eq. 0 ) then 
       gutl = gutl / dble(n1)/dble(n3) 
       gutu = gutu / dble(n1)/dble(n3) 

       !.... ut = sqrt ( 2*ut * rhy * rre ) 
       rhyl= 1.d0 / ( y(i_v,1) - y(i_v,0))
       rhyu= 1.d0 / ( y(i_v,n2)- y(i_v,n2-1))

       gutl = sqrt ( 2.d0 * gutl * rhyl / re ) 
       gutu = sqrt ( 2.d0 * gutu * rhyu / re ) 

       if (transportee) then
          gutl_tee = gutl_tee / dble(n1)/dble(n3)
          gutu_tee = gutu_tee / dble(n1)/dble(n3) 
          !gutl_tee = sqrt ( 2.d0 * gutl_tee * rhyl / re ) 
          !gutu_tee = sqrt ( 2.d0 * gutu_tee * rhyu / re ) 
       endif
       
    endif 


    call mpi_barrier ( mycomm, ierr ) 

!.... div max 

    if (p_order .eq. 2) then
    call div2nd ( div , mdiv )
    else
    call div4th ( div , mdiv )
    endif

    call mpi_reduce (mdiv, gmdiv, 1, MPI_DOUBLE_PRECISION, mpi_max, 0, mycomm , ierr) 


    gdiv = 0.d0 
    do k=ks+1,ke+1
    do j=1,n2 
    do i=is+1,ie+1
    
       gdiv = gdiv + div(i,j,k) 

    enddo 
    enddo 
    enddo 

    call mpi_reduce(gdiv, globaldiv,1,MPI_DOUBLE_PRECISION,mpi_sum,0,mycomm, ierr) 

!... cfl max 

    mcfl = 0.d0 
    do k=ks,ke
    do j=1,n2
    do i=is,ie
   
      if (itimestep .eq. 3 ) then
        cfl = dt* ( 0.25d0 * abs( (q(i_u,i,j,k) + q(i_u,i-1,j,k))*rhx   )  + &
                    0.25d0 * abs( (q(i_v,i,j,k) + q(i_v,i,j-1,k))/ hy(j))  + &
                    0.25d0 * abs( (q(i_w,i,j,k) + q(i_w,i,j,k-1))*rhz   )  + &
                    4.0d0 * rhx   * rhx   / re                                + & 
                    4.0d0 / hy(j) / hy(j) / re                             + & 
                    4.0d0 * rhz   * rhz   / re                                )      
      else
        cfl = dt* ( 0.25d0 * abs( (q(i_u,i,j,k) + q(i_u,i-1,j,k))*rhx   )  + &
                    0.25d0 * abs( (q(i_v,i,j,k) + q(i_v,i,j-1,k))/ hy(j))  + &
                    0.25d0 * abs( (q(i_w,i,j,k) + q(i_w,i,j,k-1))*rhz   )  + &
                    4.0d0 * rhx * rhx / re                                + & 
                    4.0d0 * rhz * rhz / re                                )

      endif
      mcfl = max(mcfl, cfl ) 

    enddo
    enddo
    enddo 

    call mpi_reduce (mcfl, gmcfl, 1, MPI_DOUBLE_PRECISION, mpi_max, 0, mycomm , ierr )
   
    if ( usecfl ) then 
       if ( rank .eq. 0) dt = tcfl / gmcfl * dt 
       call mpi_bcast( dt, 1, MPI_DOUBLE_PRECISION, 0, mycomm, ierr) 
    endif 

  ! if ( gmcfl .ge. 10 ) then
  !     call terminate('Code is larger than 10, end simulation ')
  !  endif
    if ( gmcfl .ne. gmcfl ) then
       call terminate( 'CFL is having NAN')
    endif

 
    if ( rank .eq. 0 ) then
       open(unit=10004,file='diagnostics.dat', status='old',form='formatted', position='append')
       !write(10004, '(i6, f10.4, f10.4, f10.4, f10.4)') i_print, time, gmcfl, gutl, gutu
       !write(10004, '(i6, f10.4, f10.4, f10.4, f10.4, f10.4, f10.4, f10.4)') i_print, time, gmcfl, gutl, gutu, gutl_tee, gutu_tee, s_tee(i_u,20)
       write(10004, '(i6, f10.4, f10.4, f10.4, f10.4, f10.4, f10.4)') i_print, time, gmcfl, gutl_tee, gutu_tee, sw_gb, sw_gt
       close(10004)


    endif

    call mpi_barrier( mycomm, ierr) 

    deallocate ( div) 

  end subroutine diag_out

!==================================================================================================================!
    
  subroutine compute_stats 
    use numerics


!.... add mean ensemble computation 
    implicit none 
    integer :: i,j,k 
    real    :: c9r16, c1r16
    real    :: usum , vsum, wsum, uusum, vvsum, wwsum, uvsum
    real    ::  dvdx, hy0, rh2x, rj2y 
    
    real, dimension(3,1:n2)     :: gum
    real, dimension(1:n2)     :: gpm, gppm, gpmpm
    real, dimension(1:n2)     :: pm_xz, ppm_xz, gpm_xz, gppm_xz


!.... update count on number of snapshots 
    iens = iens + 1
    c9r16 = 9.d0 / 16.d0
    c1r16 = 1.d0 / 16.d0


    !----------------- 
    ! first and second order stats 
    ! for u,w 
    !-----------------

    do k=ks,ke
    do j=1,n2
    do i=is,ie
       um(i_u,j) = um(i_u,j)  + q(i_u,i,j,k)
       um(i_w,j) = um(i_w,j)  + q(i_w,i,j,k)
    enddo
    enddo
    enddo

    do k=ks,ke
    do j=1,n2
    do i=is,ie
       uum(i_u,j)= uum(i_u,j) + q(i_u,i,j,k)*q(i_u,i,j,k) !   (q(i_u,i,j,k)-gum(i_u,j))*(q(i_u,i,j,k)-gum(i_u,j))
       uum(i_w,j)= uum(i_w,j) + q(i_w,i,j,k)*q(i_w,i,j,k)
    enddo
    enddo
    enddo

    !-----------------
    ! first and second order stats 
    ! for v 
    !-----------------

    do k=ks,ke
    do j=1,n2-1 
    do i=is,ie

       um(i_v,j) = um(i_v,j) + q(i_v,i,j,k) 
       uum(i_v,j)= uum(i_v,j)+ q(i_v,i,j,k)*q(i_v,i,j,k) 
    enddo 
    enddo 
    enddo 


    !----------------
    ! first order stats for 
    ! resolved reynolds stress 
    !----------------

    do k=ks,ke
    do j=1,n2 
    do i=is,ie

       uvm(j) = uvm(j) + ( c9r16* (q(i_u,i,j,k)+q(i_u,i+1,j,k)) - c1r16* (q(i_u,i-1,j,k)+q(i_u,i+2,j,k)) )* & 
                         ( c9r16* (q(i_v,i,j,k)+q(i_v,i,j-1,k)) - c1r16* (q(i_v,i,j-2,k)+q(i_v,i,j+1,k)) ) 
    enddo 
    enddo 
    enddo 


    ppm_xz = 0.d0
    pm_xz = 0.d0

    !-------------
    ! pressure statistics
    !-------------
    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1
      pm (j) = pm (j) + p_old(i,j,k)
      ppm(j) = ppm(j) + p_old(i,j,k)*p_old(i,j,k)
    enddo
    enddo
    enddo

    gpm = 0.d0
    do i= 0, p1*p3-1
    call mpi_reduce ( pm , gpm , n2, MPI_DOUBLE_PRECISION, mpi_sum, i, mycomm, ierr)
    enddo
    gpm = gpm / dble(n1*n3*iens)

    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1
       p_rms_inst(j) = p_rms_inst(j) + (p_old(i,j,k)-gpm(j))*(p_old(i,j,k)-gpm(j))
    enddo
    enddo
    enddo

  end subroutine compute_stats

!==================================================================================================================-!

  subroutine output_stat( is_final)  


!.... output mean statistics 

    implicit none 
    integer :: i,j,k ,l
    real, dimension(3,1:n2)     :: gum, guum 
    real, dimension(1:n2)       :: guvm, gvtm, guvm_s, guvm_sq, gpm, gppm, gpmpm, gp_rms, gp_rms_inst
    real, dimension(1:n2)       :: gconst_m, gconst_sq
    integer, optional           :: is_final 
    character(len=3)            :: file_suffix 
    character(len=45)           :: fnameu, fnamev, fnamevt, fnamep

    !----- 
    ! sgs stats output ...
    !-----
    real, dimension(6,1:n2)     :: gtau_m, gtau_sq
    


    gum = 0.d0 
    guum = 0.d0 
    guvm = 0.d0 
    gvtm = 0.d0 
    guvm_s = 0.d0 
    guvm_sq= 0.d0 
    gtau_m = 0.d0 
    gtau_sq= 0.d0 
    gpm = 0.d0
    gppm = 0.d0
    gp_rms_inst = 0.d0
    gp_rms      = 0.d0

    write(file_suffix, '(i3.3)') noutput 
    if ( present(is_final)  ) then
       if ( is_final .eq. 1) then 
          fnameu = 'output/mean_profile/mean_u_'//file_suffix//'.dat' 
          fnamev = 'output/mean_profile/mean_v_'//file_suffix//'.dat' 
          fnamevt= 'output/mean_profile/mean_vt_'//file_suffix//'.dat' 
          fnamep = 'output/mean_profile/mean_p_'//file_suffix//'.dat'
       else 
          fnameu = 'output/mean_profile/mean_u.dat'
          fnamev = 'output/mean_profile/mean_v.dat'
          fnamevt= 'output/mean_profile/mean_vt.dat'
          fnamep = 'output/mean_profile/mean_p.dat'
       endif 
    else  
       fnameu = 'output/mean_profile/mean_u.dat' 
       fnamev = 'output/mean_profile/mean_v.dat' 
       fnamevt= 'output/mean_profile/mean_vt.dat' 
       fnamep = 'output/mean_profile/mean_p.dat'
    endif 

    
    !----
    ! collect and average over all processes -- this isnt efficient, but we 
    ! only have to do this once per simulation when we dump the files
    !----

    call mpi_reduce ( um , gum , nvel*n2, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr) 
    call mpi_reduce ( uum, guum, nvel*n2, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr) 
    call mpi_reduce ( uvm, guvm,      n2, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr) 
    call mpi_reduce ( pm , gpm, n2      , MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr)
    call mpi_reduce ( ppm ,gppm, n2      , MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr)
    call mpi_reduce ( p_rms_inst, gp_rms_inst, n2      , MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr)
    call mpi_reduce ( vtm, gvtm,      n2, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr) 
    call mpi_reduce ( uvm_sgs, guvm_s,n2, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr) 
    call mpi_reduce ( uvm_sgs_sq, guvm_sq, n2, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr ) 
    call mpi_reduce ( tau_mean, gtau_m, n2*6, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr ) 
    call mpi_reduce ( tau_sq, gtau_sq, n2*6, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr  ) 

    call mpi_reduce ( const_m, gconst_m, n2, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr ) 
    call mpi_reduce ( const_sq, gconst_sq, n2, MPI_DOUBLE_PRECISION, mpi_sum, 0, mycomm, ierr) 
    
    gum = gum / dble(n1*n3*iens) 
    guum = guum / dble(n1*n3*iens) 
    guvm = guvm / dble(n1*n3*iens) 
    gvtm = gvtm / dble(n1*n3*iens) 
    guvm_s = guvm_s / dble(n1*n3*iens) 
    guvm_sq= guvm_sq/ dble(n1*n3*iens) 
    gtau_sq = gtau_sq / dble(n1*n3*iens) 
    gtau_m = gtau_m / dble(n1*n3*iens) 
    gpm = gpm / dble(n1*n3*iens)
    gppm = gppm / dble(n1*n3*iens)
    gp_rms_inst = gp_rms_inst / dble(n1*n3*iens)
    gp_rms = p_rms / dble(iens) 

    gconst_m = gconst_m / dble(iens*size) 
    gconst_sq = gconst_sq / dble(iens*size) 

    if ( rank .eq. 0 ) then 

       !----------------
       ! output u,v,w statistics 
       !---------------

       !------------- 
       ! u, w, uv  to file 
       !------------
       open(unit=233+rank,file=fnameu, status='unknown')
       write(233+rank,*) '# y    y+  ubar  wbar  uubar  wwbar  uvbar'
       write(233+rank,*) '#', iens
    
       do j=1,n2 
          write(233+rank,'(9g20.8)') y(i_u,j) , (y(i_u,j) + 1.d0)*re , gum(i_u,j) , gum(i_w,j), guum(i_u,j), & 
                               guum(i_w,j), guvm(j), guvm_s(j), guvm_sq(j)  
       enddo

       close(unit=233+rank) 

       !-------------
       ! v stats to file 
       !-------------
       open(unit=2440+rank,file=fnamev, status='unknown') 
       write(2440+rank,*) '# y   y+   vbar  vvbar' 

       do j=1,n2-1 
          write(2440+rank,'(4g20.8)') y(i_v,j) , (y(i_v,j) + 1.d0)*re, gum(i_v,j), guum(i_v,j) 
       enddo

       close(unit=2440) 
       !-----------
       ! p stats to file
       !-----------

       open(unit=5555+rank, file=fnamep, status='unknown')
       write(5555+rank,*) '#y       y+       pbar      ppbar    pbar*pbar   P_rms    p_rms'

       do j=1,n2
          write(5555+rank,'(7g20.8)') y(i_u,j), (y(i_u,j) + 1.d0)*re, gpm(j), gppm(j),  gp_rms_inst(j) 
       enddo

       close(unit=5555+rank)

    endif 

    !--- 
    ! iterate ntat
    ! iteration done here after upgrade to 
    ! parallel io 
    !---- 
    noutput =  noutput + 1 

    call mpi_barrier( mycomm, ierr) 

  end subroutine output_stat


  subroutine init_field_laminar

!.... create initial field based on some random noise 
!.... modified from morinishi code 

    use numerics 

    implicit none 
    real, dimension(4)                               :: iseed 
    real, dimension(nvel)                            :: qsum , rms, gqsum
    real                                             :: yp, fp, amp, umean 
    real                                             :: up, cnstk, cnstb, fup, dfdup, urc, vrc, wrc
    integer                                          :: i,j,k,i_q, l
    complex, dimension(nvel,is-1:ie+1,0:n2)          :: q_init
    complex, dimension(is-1:ie+1,0:n2)               :: p_init
    complex, dimension(is-1:ie+1,0:n2)               :: p_analytical
    complex                                          :: SIG_p, SIG_m
    complex, dimension(is-1:ie+1)                    :: residue
    complex, parameter                               :: imag = (0,1)
    real                                             :: kk, k_p, k_m
    real                                             :: dx, dz, sigma_theo
    real, dimension(is-1:ie+1)                       :: xx
    real                                  :: intf_cnst_const0, intf_cnst_const1, intf_cnst_const2, intf_cnst_const3
    real   :: alpha

    dx = chlx  / dble(n1/pi)
    dz = chlz  / dble(n3/pi)

    time = 0.d0
    q    = 0.d0
    q_init   = 0.d0
    p_init = 0.d0
    p_analytical = 0.d0
    residue = 0.d0

    q_old = 0.d0
    p_old = 0.d0

    iseed(1) = 7 + 2*rank
    iseed(2) = 563 + rank
    iseed(3) = 777
    iseed(4) = 3099 - 3*rank

    do k=ks,ke
    do j=1,n2
    do i_q = 1,nvel

      call dlarnv ( 1, iseed,nx, q(i_q,is:ie,j,k) )

    enddo
    enddo
    enddo


!.... random noise with zero mean 

    do j=1,n2

       qsum = 0.d0

       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          qsum(i_q) = qsum(i_q) + q(i_q, i,j,k)

       enddo
       enddo
       enddo


       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr)
       gqsum = gqsum / dble(n1) / dble(n3)

       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          q(i_q, i,j,k) = q(i_q, i,j,k) - gqsum(i_q)
       enddo
       enddo
       enddo

    enddo


    amp = 0.0001d0 !4.5d0 !2.5d0
    urc = 0.0001d0
    vrc = 0.0001d0
    wrc = 0.0001d0

!
    do j=1,n2

       yp = re* (1.d0 - abs(y(i_u,j)))
       fp = 1.d0 - exp(-yp/25.d0)
       rms(i_u) = sqrt(urc * fp*fp*abs(y(i_u,j)) )
       rms(i_w) = sqrt(wrc * fp*fp*abs(y(i_u,j)) )

       yp = re* (1.d0 - abs(y(i_v,j)))
       fp = 1.d0 - exp(-yp/25.d0)
       rms(i_v) = sqrt(vrc * fp*fp*abs(y(i_v,j)) )

       qsum     = 0.d0
       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          qsum(i_q) = qsum(i_q) + q(i_q,i,j,k)*q(i_q,i,j,k)
       enddo
       enddo
       enddo


       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr)
       gqsum = gqsum / dble(n1) /dble(n3)

       do i_q=1,nvel
          gqsum(i_q) = amp*rms(i_q)* sqrt(1.d0/gqsum(i_q) )
       enddo



 !      do k=ks,ke
 !      do i=is,ie
 !      do i_q=1,nvel

 !         q(i_q,i,j,k) = q(i_q,i,j,k)*gqsum(i_q)
 !      enddo
 !      enddo
 !      enddo

    enddo

    do j=1,n2
       up = re* (1.d0 - y(i_u,j)*y(i_u,j))/2.d0
       do k=ks,ke
       do i=is,ie

          q(i_u,i,j,k) = q(i_u,i,j,k) + up
       enddo
       enddo
    enddo




   call bcuvw

  end subroutine init_field_laminar

!=====================================================================================================================!

  subroutine init_field

!.... create initial field based on some random noise 
!.... modified from morinishi code 

    use numerics 

    implicit none 
    real, dimension(4)                               :: iseed 
    real, dimension(nvel)                            :: qsum , rms, gqsum
    real                                             :: yp, fp, amp, umean 
    real                                             :: up, cnstk, cnstb, fup, dfdup, urc, vrc, wrc
    integer                                          :: i,j,k,i_q, l
    complex, dimension(nvel,is-1:ie+1,0:n2)          :: q_init
    complex, dimension(is-1:ie+1,0:n2)               :: p_init
    complex, dimension(is-1:ie+1,0:n2)               :: p_analytical
    complex                                          :: SIG_p, SIG_m
    complex, dimension(is-1:ie+1)                    :: residue
    complex, parameter                               :: imag = (0,1)
    complex                                          :: mm_p, mm_m , AA_m, AA_p
    real                                             :: kk, k_p, k_m
    real                                             :: dx, dz, sigma_theo
    real, dimension(is-1:ie+1)                       :: xx
    real                                  :: intf_cnst_const0, intf_cnst_const1, intf_cnst_const2, intf_cnst_const3
    real   :: alpha

    dx = chlx  / dble(n1/pi)
    dz = chlz  / dble(n3/pi)

    time = 0.d0
    q    = 0.d0
    q_init   = 0.d0
    p_init = 0.d0
    p_analytical = 0.d0
    residue = 0.d0

    q_old = 0.d0
    p_old = 0.d0

    iseed(1) = 7 + 2*rank
    iseed(2) = 563 + rank
    iseed(3) = 777
    iseed(4) = 3099 - 3*rank

!.... returns a vector of n random real numbers from a uniform or normal distribution.

    do k=ks,ke
    do j=1,n2
    do i_q = 1,nvel

       call dlarnv ( 1, iseed,nx, q(i_q,is:ie,j,k) )

    enddo
    enddo
    enddo


!.... random noise with zero mean 

    do j=1,n2

       qsum = 0.d0

       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          qsum(i_q) = qsum(i_q) + q(i_q, i,j,k)

       enddo
       enddo
       enddo


       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr)
       gqsum = gqsum / dble(n1) / dble(n3)

       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          q(i_q, i,j,k) = q(i_q, i,j,k) - gqsum(i_q)

       enddo
       enddo
       enddo

    enddo

    amp = 15.0d0 !2.5d0
    urc = 10.1d0
    vrc = 5.0d0
    wrc = 5.3d0

!    amp = 10.d0! 4.5d0
!    urc = 5.1d0
!    vrc = 1.0d0
!    wrc = 2.3d0


    do j=1,n2

       yp = re* (1.d0 - abs(y(i_u,j)))
       fp = 1.d0 - exp(-yp/25.d0)
       rms(i_u) = sqrt(urc * fp*fp*abs(y(i_u,j)) )
       rms(i_w) = sqrt(wrc * fp*fp*abs(y(i_u,j)) )

       yp = re* (1.d0 - abs(y(i_v,j)))
       fp = 1.d0 - exp(-yp/25.d0)
       rms(i_v) = sqrt(vrc * fp*fp*abs(y(i_v,j)) )

       qsum     = 0.d0
       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          qsum(i_q) = qsum(i_q) + q(i_q,i,j,k)*q(i_q,i,j,k)
       enddo
       enddo
       enddo



       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr)
       gqsum = gqsum / dble(n1) /dble(n3)

       do i_q=1,nvel
          gqsum(i_q) = amp*rms(i_q)* sqrt(1.d0/gqsum(i_q) )
       enddo



       do k=ks,ke
       do i=is,ie
       do i_q=1,nvel

          q(i_q,i,j,k) = q(i_q,i,j,k)*gqsum(i_q)
       enddo
       enddo
       enddo

    enddo



!.... spalding's law 

    cnstk = 0.4d0
    cnstb = 5.d0
    up    = re* (1.d0 - abs(y(i_u,1)))

    do j=1,n2
       yp = re* (1.d0 - abs(y(i_u,j)) )
       do l=1,20
          fup = up - yp + exp(-cnstk*cnstb)* &
                     ( exp(cnstk*up) -1.d0 - (cnstk*up) - (cnstk*up)**2/2.d0 - (cnstk*up)**3/6.d0 )

          dfdup = 1.d0 + cnstk*exp(-cnstk*cnstb)* &
                     ( exp(cnstk*up) - 1.d0- (cnstk*up) - (cnstk*up)**2/2.d0 )

          up = up - fup/dfdup
       enddo

       do k=ks,ke
       do i=is,ie

          q(i_u,i,j,k) = q(i_u,i,j,k) + up  
       enddo
       enddo
    enddo

    p_old = 0.d0


   call bcuvw

  end subroutine init_field
!=====================================================================================================================!

  subroutine init_field_shs

!.... create initial field based on some random noise
!.... copied from Jongmin's code
    
    use numerics

    implicit none 
    real, dimension(4)                               :: iseed 
    real, dimension(nvel)                            :: qsum , rms, gqsum
    real                                             :: yp, fp, amp, umean 
    real                                             :: up, cnstk, cnstb, fup, dfdup, urc, vrc, wrc
    integer                                          :: i,j,k,i_q, l
    complex, dimension(nvel,is-1:ie+1,0:n2)          :: q_init
    complex, dimension(is-1:ie+1,0:n2)               :: p_init
    complex, dimension(is-1:ie+1,0:n2)               :: p_analytical
    complex                                          :: SIG_p, SIG_m
    complex, dimension(is-1:ie+1)                    :: residue
    complex, parameter                               :: imag = (0,1)
    complex                                          :: mm_p, mm_m , AA_m, AA_p
    real                                             :: kk, k_p, k_m
    real                                             :: dx, dz, sigma_theo
    real, dimension(is-1:ie+1)                       :: xx
    real                                  :: intf_cnst_const0, intf_cnst_const1, intf_cnst_const2, intf_cnst_const3
    real   :: alpha, u_slip

    dx = chlx  / dble(n1/pi)
    dz = chlz  / dble(n3/pi)

    time = 0.d0
    q    = 0.d0
    q_init   = 0.d0
    p_init = 0.d0
    p_analytical = 0.d0
    residue = 0.d0

    q_old = 0.d0
    p_old = 0.d0

    u_slip = slip_mean

    iseed(1) = 7 + 2*rank
    iseed(2) = 563 + rank
    iseed(3) = 777
    iseed(4) = 3099 - 3*rank

!.... returns a vector of n random real numbers from a uniform or normal distribution.

    do k=ks,ke
    do j=1,n2
    do i_q = 1,nvel

       call dlarnv ( 1, iseed,nx, q(i_q,is:ie,j,k) )

    enddo
    enddo
    enddo


!.... random noise with zero mean 

    do j=1,n2

       qsum = 0.d0

       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          qsum(i_q) = qsum(i_q) + q(i_q, i,j,k)

       enddo
       enddo
       enddo


       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr)
       gqsum = gqsum / dble(n1) / dble(n3)

       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          q(i_q, i,j,k) = q(i_q, i,j,k) - gqsum(i_q)

       enddo
       enddo
       enddo

    enddo

    if (re .gt. 190) then
      amp = 5.5d0
      urc = 10.1d0
      vrc = 1.0d0
      wrc = 2.3d0
    elseif (re .le. 50) then
      amp = 2.5d0
      urc = 1.1d0
      vrc = 1.0d0
      wrc = 1.3d0
    else
      amp = 30.0d0
      urc = 20.1d0
      vrc = 15.0d0
      wrc = 15.3d0
    endif
       

    do j=1,n2

       yp = re* (1.d0 - abs(y(i_u,j)))
       fp = 1.d0 - exp(-yp/25.d0)
       rms(i_u) = sqrt(urc * fp*fp*abs(y(i_u,j)) )
       rms(i_w) = sqrt(wrc * fp*fp*abs(y(i_u,j)) )

       yp = re* (1.d0 - abs(y(i_v,j)))
       fp = 1.d0 - exp(-yp/25.d0)
       rms(i_v) = sqrt(vrc * fp*fp*abs(y(i_v,j)) )

       qsum     = 0.d0
       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          qsum(i_q) = qsum(i_q) + q(i_q,i,j,k)*q(i_q,i,j,k)
       enddo
       enddo
       enddo



       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr)
       gqsum = gqsum / dble(n1) /dble(n3)

       do i_q=1,nvel
          gqsum(i_q) = amp*rms(i_q)* sqrt(1.d0/gqsum(i_q) )
       enddo



       do k=ks,ke
       do i=is,ie
       do i_q=1,nvel

          q(i_q,i,j,k) = q(i_q,i,j,k)*gqsum(i_q)
       enddo
       enddo
       enddo

    enddo



!.... spalding's law 

    cnstk = 0.4d0
    cnstb = 5.d0
    up    = re* (1.d0 - abs(y(i_u,1)))

    do j=1,n2
       yp = re* (1.d0 - abs(y(i_u,j)) )
       do l=1,20
          fup = up - yp + exp(-cnstk*cnstb)* &
                     ( exp(cnstk*up) -1.d0 - (cnstk*up) - (cnstk*up)**2/2.d0 - (cnstk*up)**3/6.d0 )

          dfdup = 1.d0 + cnstk*exp(-cnstk*cnstb)* &
                     ( exp(cnstk*up) - 1.d0- (cnstk*up) - (cnstk*up)**2/2.d0 )

          up = up - fup/dfdup
       enddo

       do k=ks,ke
       do i=is,ie

          q(i_u,i,j,k) = q(i_u,i,j,k) + up + u_slip
       enddo
       enddo
    enddo

    p_old = 0.d0


   call bcuvw

  end subroutine init_field_shs
!=====================================================================================================================!

  subroutine init_field_tee


    use numerics_tee 

    implicit none 
    integer                                          :: i,j,k,i_q, l, jBF
    real                                             :: deltaBF


!.... initial condition for transportee variables

    if ( transportee ) then 

       if ( ICswitch .eq. 0) then 
       do j=1,n2
       do i=is,ie
       do k=ks,ke
       do i_q = 1,nvel
         q_tee(i_q, i,j,k) = 0.0d0
         q_source(i_q,i,j,k) = 0.0d0
       enddo
       enddo
       enddo
       enddo
       endif

       if ( ICswitch .gt. 1000 .and. ICswitch .lt. 1999 ) then 
       do j=1,n2
       do i=is,ie
       do k=ks,ke
       do i_q = 1,nvel
         q_tee(i_q, i,j,k) = 0.0d0
         q_source(i_q,i,j,k) = 0.0d0
       enddo
       enddo
       enddo
       enddo
       endif

       if ( ICswitch .eq. 1) then 
       do j=1,n2
       do i=is,ie
       do k=ks,ke
         q_tee(i_u, i,j,k) = y(i_u,j)
         q_tee(i_v, i,j,k) = 0.d0 
         q_tee(i_w, i,j,k) = 0.d0 

         q_source(i_u, i,j,k) = y(i_u,j)
         q_source(i_v, i,j,k) = 0.d0 
         q_source(i_w, i,j,k) = 0.d0 
       enddo
       enddo
       enddo
       v1_bottom = -1.d0
       v1_top = 1.d0
       endif

       if ( ICswitch .eq. 2) then 
       do j=1,n2
       do i=is,ie
       do k=ks,ke
         q_tee(i_u, i,j,k) = 1.d0 - y(i_u,j) * y(i_u,j)
         q_tee(i_v, i,j,k) = 0.d0 
         q_tee(i_w, i,j,k) = 0.d0 

         q_source(i_u, i,j,k) = 1.d0 - y(i_u,j) * y(i_u,j)
         q_source(i_v, i,j,k) = 0.d0 
         q_source(i_w, i,j,k) = 0.d0 
       enddo
       enddo
       enddo
       v1_bottom = 0.d0
       v1_top = 0.d0
       endif

       if ( ICswitch .eq. 31) then 
       do j=1,n2
       do i=is,ie
       do k=ks,ke
         q_tee(i_u, i,j,k) = y(i_u,j) ** 3
         q_tee(i_v, i,j,k) = 0.d0 
         q_tee(i_w, i,j,k) = 0.d0 

         q_source(i_u, i,j,k) = y(i_u,j) ** 3
         q_source(i_v, i,j,k) = 0.d0 
         q_source(i_w, i,j,k) = 0.d0 
       enddo
       enddo
       enddo
       v1_bottom = -1.d0
       v1_top = 1.d0
       endif

       if ( ICswitch .eq. 32) then 
       do j=1,n2
       do i=is,ie
       do k=ks,ke
         q_tee(i_u, i,j,k) = (1.d0 + y(i_u,j)) ** 3
         q_tee(i_v, i,j,k) = 0.d0 
         q_tee(i_w, i,j,k) = 0.d0 

         q_source(i_u, i,j,k) = (1.d0 + y(i_u,j)) ** 3
         q_source(i_v, i,j,k) = 0.d0 
         q_source(i_w, i,j,k) = 0.d0 
       enddo
       enddo
       enddo
       v1_bottom = 0.d0
       v1_top = 8.d0
       endif

       ! BF edit
       if ( IMFMswitchBF ) then
         jBF = ICswitchBF
         deltaBF = y(i_u,jBF+1) - y(i_u,jBF)

         do j=1,n2
         do i=is,ie
         do k=ks,ke
            q_tee(i_u,   i,j,k) = 0.d0
            q_tee(i_v,   i,j,k) = 0.d0
            q_tee(i_w,   i,j,k) = 0.d0
            q_source(i_u,i,j,k) = 0.d0
            q_source(i_v,i,j,k) = 0.d0
            q_source(i_w,i,j,k) = 0.d0
            if ( j .gt. jBF ) then
              q_tee(i_u,   i,j,k) = deltaBF
              q_source(i_u,i,j,k) = deltaBF
            endif
         enddo
         enddo
         enddo
         v1_bottom = 0.d0
         !v1_top = deltaBF
         v1_top = 0.d0 ! for patterning only
       endif
    
       p_old_tee = 0.d0

    endif

!.... initial condition for IMFM variables

    if ( transportee ) then 
       do j=1,n2
       do i_q = 1,nvel
         
         s_tee(i_q,j) = 0.d0
         
       enddo
       enddo
    endif

   call bcuvw_tee
   call bcuvw_source

  end subroutine init_field_tee
!=====================================================================================================================!

  subroutine setup_BC_tee

    implicit none 
    integer           :: jBF

!.... initial condition for transportee variables

    if ( transportee ) then 

       if ( ICswitch .eq. 0) then 
         v1_bottom = 0.d0
         v1_top = 0.d0
       endif


       if ( ICswitch .eq. 1) then 
         v1_bottom = -1.d0
         v1_top = 1.d0
       endif

       if ( ICswitch .eq. 2) then 
         v1_bottom = 0.d0
         v1_top = 0.d0
       endif

       if ( ICswitch .eq. 31) then 
         v1_bottom = -1.d0
         v1_top = 1.d0
       endif

       if ( ICswitch .eq. 32) then 
         v1_bottom = 0.d0
         v1_top = 8.d0
       endif

       ! BF edit
       if ( IMFMswitchBF ) then
         jBF = ICswitchBF
         v1_bottom = 0.d0
         v1_top = 0.d0 ! for aptterning only
         !v1_top = y(i_u,jBF+1) - y(i_u,jBF)
       endif
    
    endif


  end subroutine setup_BC_tee
!=====================================================================================================================!

  subroutine output_clip_stats 

    implicit none 
    integer, dimension(2) :: gnclip



    ! reduce clipping values from all rpocs
    call mpi_reduce ( nclip(1), gnclip(1), 2, mpi_integer, mpi_sum, 0, mycomm, ierr) 

    if ( rank .eq. 0) then 
       open(unit=250, file='output/clip_stats.dat', form='formatted', status='old', position='append')
       write(250, '(d20.8, i20.7, i20.7)') time, gnclip(1), gnclip(2) 
       close(250) 
    endif 


  end subroutine output_clip_stats 
!======================================================================================================================!

  subroutine output_fm_pg 

    !---- 
    ! dump mean pressure gradient to file 
    !---- 

    implicit none 

    if ( rank .eq. 0 ) then 
       open(unit=2500,file='output/fm_pressure_grad.dat', form='formatted', status='old', position='append') 
       write(2500, '(d20.8, d20.8)') time, pgm 
       close(2500) 
    endif

  end subroutine output_fm_pg
!=======================================================================================================================!

  subroutine dump_post_interface(interface_var, fname)

    implicit none
    character(*) :: fname
    real, dimension(is-1:ie+1, 1:2, ks-1:ke+1)               , intent(in):: interface_var
    real, dimension(:,:,:)      , allocatable :: interface_data
    real, dimension(:)        , allocatable :: recvbuf ,  rbuf_p,rbuf_u, rbuf_w
    integer                                 :: i,j,k
    integer                                 :: itblock, rshft, n, proc, count,ind
    integer, dimension(:), allocatable      :: ks_vec, is_vec
    logical                                 :: debug_collapse
    logical                                 :: write_collapsed_grid
    real                                    :: hx, hz, xx, zz, u, v,w
    character(len=3)                        :: nstr

    !---------------------
    ! gather the start and ending indices for 
    ! all processors to master 
    ! currently assumes that all the blocks are of the same size
    !--------------------

    allocate ( ks_vec(0:size-1), is_vec(0:size-1))

    call mpi_gather( ks, 1, mpi_integer, ks_vec(0), 1, mpi_integer, 0, mycomm,ierr)
    call mpi_gather( is, 1, mpi_integer, is_vec(0), 1, mpi_integer, 0, mycomm,ierr)

    nx = int ( n1 / p1 )
    nz = int ( n3 / p3 )

    itblock  = (nx + 2 ) * (nz + 2 ) * 2

    allocate ( recvbuf( size*itblock) )

    rshft = rank * itblock + 1
    call mpi_gather ( interface_var(is-1,1,ks-1), itblock, MPI_DOUBLE_PRECISION,&
                      recvbuf(1), itblock, MPI_DOUBLE_PRECISION, 0, mycomm,ierr)

    allocate ( interface_data(0:n1-1, 2, 0:n3-1))

    if ( rank .eq. 0 ) then

       do proc=0,size-1

          rshft = proc * itblock +1
          count = 0

          do k= ks_vec(proc)-1, ks_vec(proc)+nz-1+1
          do j= 1,2
          do i= is_vec(proc)-1, is_vec(proc)+nx-1+1
             ind = rshft+count
                                                                                 
             if ( k .ge. 0 .and. k .le. n3-1 .and. &
                  i .ge. 0 .and. i .le. n1-1 .and. &
                  j .ge. 1 .and. j .le. 2  ) then

                interface_data(i,j,k) = recvbuf(rshft+count)
             endif

             count = count + 1
          enddo
          enddo
          enddo

       enddo

       write(nstr, '(i3.3)') n_interface
       open(unit=502,file='output/post/'//fname//'.dat', form='unformatted', status='unknown')
       write(502) n1,n3
       write(502) interface_data
       close(502)

    endif


    deallocate ( recvbuf)
    deallocate ( interface_data     )
    deallocate ( ks_vec  )
    deallocate ( is_vec  )

  end subroutine dump_post_interface

!================================================================================================================================!
  subroutine dump_collapsed_shs_velocity

    !-----------------------------------------------------------------
    ! collapse entire velocity and pressure fields onto the master node and output one file 
    ! currently assumes that all subdomains are the same size
    !-----------------------------------------------------------------

    implicit none 
    real, dimension(:,:,:,:  ), allocatable :: qq
    real, dimension(:,:,:)    , allocatable :: pp 
    real, dimension(:)        , allocatable :: recvbuf , rbuf_p
    integer                                 :: i,j,k 
    integer                                 :: qblock , rshft, n, proc, count, ind, pblock
    integer, dimension(:), allocatable      :: ks_vec, is_vec
    logical                                 :: debug_collapse
    logical                                 :: write_collapsed_grid
    real                                    :: hx, hz, xx, zz, u, v,w 
    character(len=3)                        :: nstr

    !---------------------
    ! gather the start and ending indices for 
    ! all processors to master 
    ! currently assumes that all the blocks are of the same size
    !--------------------

    allocate ( ks_vec(0:size-1), is_vec(0:size-1))

    call mpi_gather( ks, 1, mpi_integer, ks_vec(0), 1, mpi_integer, 0, mycomm, ierr) 
    call mpi_gather( is, 1, mpi_integer, is_vec(0), 1, mpi_integer, 0, mycomm, ierr) 

    nx = int ( n1 / p1 ) 
    nz = int ( n3 / p3 ) 
    

    qblock  = nvel* (nx + 2*halo_x )* (nz + 2*halo_z )* (n2 + 2*halo_yu ) 
    allocate ( recvbuf( size*qblock) )

    !--------------------------
    ! gather q from all processors to master
    !--------------------------
    rshft = rank * qblock + 1
    call mpi_gather ( q(1,is-halo_x,1-halo_yu,ks-halo_z), qblock, MPI_DOUBLE_PRECISION, & 
                      recvbuf(1), qblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)  

    allocate ( qq(1:nvel,0:n1-1,1:n2, 0:n3-1))

    !--------------------------
    ! master node organizes the received data 
    !-------------------------

    if ( rank .eq. 0 ) then 

       do proc=0,size-1
    
          !-----------------
          ! parse recv buff and extract values 
          ! into qq 
          !-----------------
          rshft = proc * qblock +1 
          count = 0 
          
          do k= ks_vec(proc)-halo_z, ks_vec(proc)+nz-1+halo_z
          do j= 1-halo_yu, n2+halo_yu
          do i= is_vec(proc)-halo_x, is_vec(proc)+nx-1+halo_x
          do n= 1,nvel
             
             ind = rshft+count
             if ( k .ge. 0 .and. k .le. n3-1 .and. & 
                  i .ge. 0 .and. i .le. n1-1 .and. & 
                  j .ge. 1 .and. j .le. n2  ) then 

                qq(n,i,j,k) = recvbuf(rshft+count) 
             endif 

             count = count + 1 
          enddo 
          enddo 
          enddo 
          enddo 

       enddo 


       !---------------- 
       ! output collapsed velocities to file 
       !---------------- 

       write(nstr, '(i3.3)') n_interface
       open(unit=499,file='output/collapsed_v'//nstr//'.m', form='formatted', status='unknown') 

       do k = 1,n3-1
       do i = 1,n1-1 
          j = 1
         write(499,*) " u(",i,",",j+1,",",k,") = ",qq(i_u,i,j,k),';'
         write(499,*) " v(",i,",",j+1,",",k,") = ",qq(i_v,i,j,k),';'
         write(499,*) " w(",i,",",j+1,",",k,") = ",qq(i_w,i,j,k),';'
           
       enddo 
       enddo

       close(499)
    
    endif
    

    call mpi_barrier ( mycomm, ierr ) 


    !--- 
    ! temporary cleanup 
    !---
    deallocate ( recvbuf) 
    deallocate ( qq     ) 

    !----
    ! collapse pressure field
    !---- 

    pblock =   (nx + 2*halo_x )* (nz + 2*halo_z )* (n2 + 2*halo_yp ) 

    allocate ( pp(1:n1, 1:n2, 1:n3))
    allocate ( rbuf_p(pblock*size))
   

    !--------------------------
    ! gather p from all processors to master
    !--------------------------
    rshft = rank * pblock + 1
    call mpi_gather ( p_old(is-halo_x+1,1-halo_yp,ks-halo_z+1), pblock, MPI_DOUBLE_PRECISION, & 
                      rbuf_p(1), pblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)  


    !--------------------------
    ! master node organizes the received data 
    !-------------------------

    if ( rank .eq. 0 ) then 

       do proc=0,size-1
    
          !-----------------
          ! parse recv buff and extract values 
          ! into pp 
          !-----------------
          rshft = proc * pblock +1 
          count = 0 
          
          do k= ks_vec(proc)-halo_z+1, ks_vec(proc)+nz-1+halo_z+1
          do j= 1-halo_yp, n2+halo_yp
          do i= is_vec(proc)-halo_x+1, is_vec(proc)+nx-1+halo_x+1
             
             ind = rshft+count
             if ( k .ge. 1 .and. k .le. n3 .and. & 
                  i .ge. 1 .and. i .le. n1 .and. & 
                  j .ge. 1 .and. j .le. n2  ) then 

                pp(i,j,k) = rbuf_p(rshft+count) 
             endif 

             count = count + 1 
          enddo 
          enddo 
          enddo 

       enddo 


       !---------------- 
       ! output collapsed velocities to file 
       !---------------- 

       open(unit=501,file='output/collapsed_p'//nstr//'.m', form='formatted', status='unknown') 
       write(501,*) n1 
       write(501,*) n2 
       write(501,*) n3 

       do k=1,n3
        j=1 
       do i=1,n1 

          write(501,*) " p(",i+1,",",j+1,",",k+1,") = ",p_old(i,j,k),';'
       enddo 
       enddo 

       close(501)

    endif


    call mpi_barrier(mycomm, ierr) 


    !----------------
    ! cleanup
    !---------------
    deallocate ( rbuf_p  ) 
    deallocate ( pp      ) 
    deallocate ( ks_vec  ) 
    deallocate ( is_vec  ) 
    
  end subroutine dump_collapsed_shs_velocity
    
!================================================================================================================================!

  subroutine dump_grid_information

    !-----------------------------------------------------------------
    ! collapse entire velocity and pressure fields onto the master node and output one file 
    ! currently assumes that all subdomains are the same size
    !-----------------------------------------------------------------

    implicit none 
    integer                                 :: i,j,k 
    integer                                 :: rshft, n, proc, count, ind
    integer, dimension(:), allocatable      :: ks_vec, is_vec
    logical                                 :: debug_collapse
    logical                                 :: write_collapsed_grid
    real                                    :: hx, hz, xx, zz, u, v,w 
    character(len=3)                        :: nstr

    !--------------------

    allocate ( ks_vec(0:size-1), is_vec(0:size-1))

    call mpi_gather( ks, 1, mpi_integer, ks_vec(0), 1, mpi_integer, 0, mycomm, ierr) 
    call mpi_gather( is, 1, mpi_integer, is_vec(0), 1, mpi_integer, 0, mycomm, ierr) 

    nx = int ( n1 / p1 ) 
    nz = int ( n3 / p3 ) 
    
    write_collapsed_grid = .true.
    hx = (chlx*pi) / dble(n1) 
    hz = (chlz*pi) / dble(n3) 
    
    if ( rank .eq. 0 ) then 

       open(unit=876,file='output/collapsed_grid'//nstr//'.dat', form='formatted', status='unknown') 
       write(876, '(i3.3)') n1 

       !--- 
       ! write x grid pts 
       !--- 
       do i=1,n1 
          xx = (dble(i) - 0.5d0 ) * hx 
          write(876, '(g20.8)') xx 
       enddo 

       !--- 
       ! write y grid pts 
       !---
       write(876,'(i3.3)') n2 
       do j=1,n2 
          write(876,'(g20.8)') y(i_u,j) 
       enddo 


       !--- 
       ! write z grid pts 
       !---
       write(876,'(i3.3)') n3 
       do k=1,n3 
          zz = ( dble(k) - 0.5d0) * hz 
          write(876,'(g20.8)') zz 
       enddo 


       close(876) 
    endif 
       
    !----------------
    ! cleanup
    !---------------
    deallocate ( ks_vec  ) 
    deallocate ( is_vec  ) 
    
  end subroutine dump_grid_information
    
!================================================================================================================================!

  subroutine interface_collapse_data_field

    !-----------------------------------------------------------------
    ! collapse entire velocity and pressure fields onto the master node and output one file 
    ! currently assumes that all subdomains are the same size
    !-----------------------------------------------------------------

    implicit none 
    real, dimension(:,:,:,:  ), allocatable :: qq
    real, dimension(:,:,:)    , allocatable :: pp 
    real, dimension(:,:,:)    , allocatable :: vss
    real, dimension(:)        , allocatable :: recvbuf , rbuf_p, rbuf_intf
    integer                                 :: i,j,k 
    integer                                 :: qblock , rshft, n, proc, count, ind, pblock, intf_block
    integer, dimension(:), allocatable      :: ks_vec, is_vec
    logical                                 :: debug_collapse
    logical                                 :: write_collapsed_grid
    real                                    :: hx, hz, xx, zz, u, v,w 
    character(len=3)                        :: nstr

    !---------------------
    ! gather the start and ending indices for 
    ! all processors to master 
    ! currently assumes that all the blocks are of the same size
    !--------------------

    allocate ( ks_vec(0:size-1), is_vec(0:size-1))

    call mpi_gather( ks, 1, mpi_integer, ks_vec(0), 1, mpi_integer, 0, mycomm, ierr) 
    call mpi_gather( is, 1, mpi_integer, is_vec(0), 1, mpi_integer, 0, mycomm, ierr) 

    nx = int ( n1 / p1 ) 
    nz = int ( n3 / p3 ) 
    

    qblock  = nvel* (nx + 2*halo_x )* (nz + 2*halo_z )* (n2 + 2*halo_yu ) 
    allocate ( recvbuf( size*qblock) )

    !--------------------------
    ! gather q from all processors to master
    !--------------------------
    rshft = rank * qblock + 1
    call mpi_gather ( q(1,is-halo_x,1-halo_yu,ks-halo_z), qblock, MPI_DOUBLE_PRECISION, & 
                      recvbuf(1), qblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)  

    allocate ( qq(1:nvel,0:n1-1,1:n2, 0:n3-1))

    !--------------------------
    ! master node organizes the received data 
    !-------------------------

    if ( rank .eq. 0 ) then 

       do proc=0,size-1
    
          !-----------------
          ! parse recv buff and extract values 
          ! into qq 
          !-----------------
          rshft = proc * qblock +1 
          count = 0 
          
          do k= ks_vec(proc)-halo_z, ks_vec(proc)+nz-1+halo_z
          do j= 1-halo_yu, n2+halo_yu
          do i= is_vec(proc)-halo_x, is_vec(proc)+nx-1+halo_x
          do n= 1,nvel
             
             ind = rshft+count
             if ( k .ge. 0 .and. k .le. n3-1 .and. & 
                  i .ge. 0 .and. i .le. n1-1 .and. & 
                  j .ge. 1 .and. j .le. n2  ) then 

                qq(n,i,j,k) = recvbuf(rshft+count) 
             endif 

             count = count + 1 
          enddo 
          enddo 
          enddo 
          enddo 

       enddo 


       !---------------- 
       ! output collapsed velocities to file 
       !---------------- 

       write(nstr, '(i3.3)') noutput 
       open(unit=599,file='output/restart/collapsed_vel.dat', form='unformatted', status='unknown') 
       write(599) n1 
       write(599) n2 
       write(599) n3 
       !write(599) qq

       do k=0,n3-1
       do j=1,n2 
       do i=0,n1-1 
       do n=1,nvel
          write(599) qq(n,i,j,k) 
       enddo 
       enddo 
       enddo 
       enddo 

       close(599)

    
    endif
    

    call mpi_barrier ( mycomm, ierr ) 


    !--- 
    ! temporary cleanup 
    !---
    deallocate ( recvbuf) 
    deallocate ( qq     ) 

    !----
    ! collapse pressure field
    !---- 

    pblock =   (nx + 2*halo_x )* (nz + 2*halo_z )* (n2 + 2*halo_yp ) 

    allocate ( pp(1:n1, 1:n2, 1:n3))
    allocate ( rbuf_p(pblock*size))
   
    !--------------------------
    ! gather p from all processors to master
    !--------------------------
    rshft = rank * pblock + 1
    call mpi_gather ( p_old(is-halo_x+1,1-halo_yp,ks-halo_z+1), pblock, MPI_DOUBLE_PRECISION, & 
                      rbuf_p(1), pblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)  


    !--------------------------
    ! master node organizes the received data 
    !-------------------------

    if ( rank .eq. 0 ) then 

       do proc=0,size-1
    
          !-----------------
          ! parse recv buff and extract values 
          ! into pp 
          !-----------------
          rshft = proc * pblock +1 
          count = 0 
          
          do k= ks_vec(proc)-halo_z+1, ks_vec(proc)+nz-1+halo_z+1
          do j= 1-halo_yp, n2+halo_yp
          do i= is_vec(proc)-halo_x+1, is_vec(proc)+nx-1+halo_x+1
             
             ind = rshft+count
             if ( k .ge. 1 .and. k .le. n3 .and. & 
                  i .ge. 1 .and. i .le. n1 .and. & 
                  j .ge. 1 .and. j .le. n2  ) then 

                pp(i,j,k) = rbuf_p(rshft+count) 
             endif 

             count = count + 1 
          enddo 
          enddo 
          enddo 

       enddo 


       !---------------- 
       ! output collapsed velocities to file 
       !---------------- 

       open(unit=501,file='output/restart/collapsed_p.dat', form='unformatted', status='unknown') 
       write(501) n1 
       write(501) n2 
       write(501) n3 
       !write(501,*) pp     

       do k=1,n3
       do j=1,n2 
       do i=1,n1 
          write(501) pp(i,j,k) 
       enddo 
       enddo 
       enddo 

       close(501)

    endif

    call mpi_barrier(mycomm, ierr) 
    deallocate ( rbuf_p  )
    deallocate ( pp      )


    intf_block  = (nx + 2 ) * (nz + 2 ) * 2
    allocate ( rbuf_intf( size*intf_block) )

    !----------
    ! output collapsed grid , 
    ! cell centered values only
    !---------- 
    write_collapsed_grid = .true.
    hx = (chlx*pi) / dble(n1) 
    hz = (chlz*pi) / dble(n3) 
    
    if ( rank .eq. 0 ) then 

       open(unit=876,file='output/restart/collapsed_grid.dat', form='formatted', status='unknown') 
       write(876, '(i3.3)') n1 

       !--- 
       ! write x grid pts 
       !--- 
       do i=1,n1 
          xx = (dble(i) - 0.5d0 ) * hx 
          write(876, '(g20.8)') xx 
       enddo 

       !--- 
       ! write y grid pts 
       !---
       write(876,'(i3.3)') n2 
       do j=1,n2 
          write(876,'(g20.8)') y(i_u,j) 
       enddo 


       !--- 
       ! write z grid pts 
       !---
       write(876,'(i3.3)') n3 
       do k=1,n3 
          zz = ( dble(k) - 0.5d0) * hz 
          write(876,'(g20.8)') zz 
       enddo 


       close(876) 
    endif 
       

    !----------------
    ! cleanup
    !---------------
    deallocate ( ks_vec  ) 
    deallocate ( is_vec  ) 
    
  end subroutine interface_collapse_data_field
    
!================================================================================================================================!

  subroutine collapse_data_field 

    !-----------------------------------------------------------------
    ! collapse entire velocity and pressure fields onto the master node and output one file 
    ! currently assumes that all subdomains are the same size
    !-----------------------------------------------------------------

    implicit none 
    real, dimension(:,:,:,:  ), allocatable :: qq
    real, dimension(:,:,:)    , allocatable :: pp 
    real, dimension(:)        , allocatable :: recvbuf , rbuf_p
    integer                                 :: i,j,k 
    integer                                 :: qblock , rshft, n, proc, count, ind, pblock
    integer, dimension(:), allocatable      :: ks_vec, is_vec
    logical                                 :: debug_collapse
    logical                                 :: write_collapsed_grid
    real                                    :: hx, hz, xx, zz, u, v,w 
    character(len=3)                        :: nstr

    !---------------------
    ! gather the start and ending indices for 
    ! all processors to master 
    ! currently assumes that all the blocks are of the same size
    !--------------------

    allocate ( ks_vec(0:size-1), is_vec(0:size-1))

    call mpi_gather( ks, 1, mpi_integer, ks_vec(0), 1, mpi_integer, 0, mycomm, ierr) 
    call mpi_gather( is, 1, mpi_integer, is_vec(0), 1, mpi_integer, 0, mycomm, ierr) 

    nx = int ( n1 / p1 ) 
    nz = int ( n3 / p3 ) 
    

    qblock  = nvel* (nx + 2*halo_x )* (nz + 2*halo_z )* (n2 + 2*halo_yu ) 
    allocate ( recvbuf( size*qblock) )

    !--------------------------
    ! gather q from all processors to master
    !--------------------------
    rshft = rank * qblock + 1
    call mpi_gather ( q(1,is-halo_x,1-halo_yu,ks-halo_z), qblock, MPI_DOUBLE_PRECISION, & 
                      recvbuf(1), qblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)  

    allocate ( qq(1:nvel,0:n1-1,1:n2, 0:n3-1))

    !--------------------------
    ! master node organizes the received data 
    !-------------------------

    if ( rank .eq. 0 ) then 

       do proc=0,size-1
    
          !-----------------
          ! parse recv buff and extract values 
          ! into qq 
          !-----------------
          rshft = proc * qblock +1 
          count = 0 
          
          do k= ks_vec(proc)-halo_z, ks_vec(proc)+nz-1+halo_z
          do j= 1-halo_yu, n2+halo_yu
          do i= is_vec(proc)-halo_x, is_vec(proc)+nx-1+halo_x
          do n= 1,nvel
             
             ind = rshft+count
             if ( k .ge. 0 .and. k .le. n3-1 .and. & 
                  i .ge. 0 .and. i .le. n1-1 .and. & 
                  j .ge. 1 .and. j .le. n2  ) then 

                qq(n,i,j,k) = recvbuf(rshft+count) 
             endif 

             count = count + 1 
          enddo 
          enddo 
          enddo 
          enddo 

       enddo 


       !---------------- 
       ! output collapsed velocities to file 
       !---------------- 

       write(nstr, '(i3.3)') noutput 
       open(unit=499,file='output/collapsed_field'//nstr//'.dat', form='formatted', status='unknown') 
       write(499,*) n1 
       write(499,*) n2 
       write(499,*) n3 

       do k=0,n3-1
       do j=1,n2 
       do i=0,n1-1 
       do n=1,nvel

          write(499,*) qq(n,i,j,k) 
       enddo 
       enddo 
       enddo 
       enddo 

       close(499)



       !---- 
       ! debug the collapse output 
       !---- 
       debug_collapse = .false. 
       if ( debug_collapse) then

          hx = (chlx*pi) / dble(n1) 
          hz = (chlz*pi) / dble(n3) 

          open(unit=500, file='output/collapse_debug.plt', form='formatted', status='unknown') 
          write(500,*) 'TITLE= "vel field"' 
          write(500,*) 'VARIABLES= "x" "y" "z" "u" "v" "w"' 
          write(500,*) 'ZONE DATAPACKING=POINT, I= ', n1, ', J= ', n2, ', K=', n3
       
          do k=0,n3-1
          do j=1,n2 
          do i=0,n1-1
    
             !----
             ! use second order interpolation for tecplot output 
             !----
             xx = (dble(i) + 0.5d0 )* hx 
             zz = (dble(k) + 0.5d0 )* hz 
             
             if ( i .ne. n1-1 ) then 
                u  = 0.5d0 * ( qq(i_u,i+1,j,k) + qq(i_u,i  ,j,k)) 
             else  
                u  = 0.5d0 * ( qq(i_u,n1-1,j,k) + qq(i_u,0,j,k)) 
             endif
             
             if ( k .ne. n3-1) then 
                w  = 0.5d0 * ( qq(i_w,i,j,k+1) + qq(i_w,i,j,  k))
             else  
                w  = 0.5d0 * ( qq(i_w,i,j,n3-1) + qq(i_w,i,j,0)) 
             endif
             
             v  = 0.5d0 * ( qq(i_v,i,j,k  ) + qq(i_v,i,j-1,k))
             
             write(500,'(6d16.8)') xx, y(i_u,j), zz, u,v,w 
          enddo
          enddo
          enddo

          close(unit=500) 
          
       endif
       

    
    endif
    

    call mpi_barrier ( mycomm, ierr ) 


    !--- 
    ! temporary cleanup 
    !---
    deallocate ( recvbuf) 
    deallocate ( qq     ) 

    !----
    ! collapse pressure field
    !---- 

    pblock =   (nx + 2*halo_x )* (nz + 2*halo_z )* (n2 + 2*halo_yp ) 

    allocate ( pp(1:n1, 1:n2, 1:n3))
    allocate ( rbuf_p(pblock*size))
   

    !--------------------------
    ! gather p from all processors to master
    !--------------------------
    rshft = rank * pblock + 1
    call mpi_gather ( p_old(is-halo_x+1,1-halo_yp,ks-halo_z+1), pblock, MPI_DOUBLE_PRECISION, & 
                      rbuf_p(1), pblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)  


    !--------------------------
    ! master node organizes the received data 
    !-------------------------

    if ( rank .eq. 0 ) then 

       do proc=0,size-1
    
          !-----------------
          ! parse recv buff and extract values 
          ! into pp 
          !-----------------
          rshft = proc * pblock +1 
          count = 0 
          
          do k= ks_vec(proc)-halo_z+1, ks_vec(proc)+nz-1+halo_z+1
          do j= 1-halo_yp, n2+halo_yp
          do i= is_vec(proc)-halo_x+1, is_vec(proc)+nx-1+halo_x+1
             
             ind = rshft+count
             if ( k .ge. 1 .and. k .le. n3 .and. & 
                  i .ge. 1 .and. i .le. n1 .and. & 
                  j .ge. 1 .and. j .le. n2  ) then 

                pp(i,j,k) = rbuf_p(rshft+count) 
             endif 

             count = count + 1 
          enddo 
          enddo 
          enddo 

       enddo 


       !---------------- 
       ! output collapsed velocities to file 
       !---------------- 

       open(unit=501,file='output/collapsed_p'//nstr//'.dat', form='formatted', status='unknown') 
       write(501,*) n1 
       write(501,*) n2 
       write(501,*) n3 

       do k=1,n3
       do j=1,n2 
       do i=1,n1 

          write(501,*) pp(i,j,k) 
       enddo 
       enddo 
       enddo 

       close(501)



       !---- 
       ! debug the collapse output 
       !---- 
       if ( debug_collapse) then

          hx = (chlx*pi) / dble(n1) 
          hz = (chlz*pi) / dble(n3) 

          open(unit=502, file='output/collapse_debug_p.plt', form='formatted', status='unknown') 
          write(502,*) 'TITLE= "p field"' 
          write(502,*) 'VARIABLES= "x" "y" "z" "p"' 
          write(502,*) 'ZONE DATAPACKING=POINT, I= ', n1, ', J= ', n2, ', K=', n3
       
          do k=1,n3
          do j=1,n2 
          do i=1,n1
    
             !----
             ! use second order interpolation for tecplot output 
             !----
             xx = (dble(i) - 0.5d0 )* hx 
             zz = (dble(k) - 0.5d0 )* hz 
             
             
             write(502,'(6d16.8)') xx, y(i_u,j), zz, pp(i,j,k)
          enddo
          enddo
          enddo

          close(unit=502) 
          
       endif

    endif


    call mpi_barrier(mycomm, ierr) 
    


    !----------
    ! output collapsed grid , 
    ! cell centered values only
    !---------- 
    write_collapsed_grid = .true.
    hx = (chlx*pi) / dble(n1) 
    hz = (chlz*pi) / dble(n3) 
    
    if ( rank .eq. 0 ) then 

       open(unit=876,file='output/collapsed_grid'//nstr//'.dat', form='formatted', status='unknown') 
       write(876, '(i3.3)') n1 

       !--- 
       ! write x grid pts 
       !--- 
       do i=1,n1 
          xx = (dble(i) - 0.5d0 ) * hx 
          write(876, '(g20.8)') xx 
       enddo 

       !--- 
       ! write y grid pts 
       !---
       write(876,'(i3.3)') n2 
       do j=1,n2 
          write(876,'(g20.8)') y(i_u,j) 
       enddo 


       !--- 
       ! write z grid pts 
       !---
       write(876,'(i3.3)') n3 
       do k=1,n3 
          zz = ( dble(k) - 0.5d0) * hz 
          write(876,'(g20.8)') zz 
       enddo 


       close(876) 
    endif 
       


    !----------------
    ! cleanup
    !---------------
    deallocate ( rbuf_p  ) 
    deallocate ( pp      ) 
    deallocate ( ks_vec  ) 
    deallocate ( is_vec  ) 
    
  end subroutine collapse_data_field
    
!==========================================================================================

  subroutine restart_interface_collapsed_data_field 

    implicit none 
    integer                              :: n1_old, n2_old, n3_old, nneigh_i, nneigh_j, nneigh_k
    integer, dimension(:), allocatable   :: neigh_i , neigh_j, neigh_k 
    real, dimension(:), allocatable      :: xprev, yprev, zprev 
    real, dimension(:,:,:,:),allocatable :: qold 
    real, dimension(:,:,:)  ,allocatable :: pgrid
    real, dimension(1:nvel)              :: rms, qsum , gqsum 
    real                                 :: fp, yp , xx, zz, hx, hz, tfax
    integer, dimension(4)                :: iseed
    integer                              :: i,j,k,n,i_q, ii, jj, kk
    logical                              :: exist_restart, exist_grid
    integer                              :: use_deprecated_restart 


    ! use this if the prev field is just a simple refine homothetically
    use_deprecated_restart = 2


    !---- 
    ! verify necessary restart files exist 
    !---- 
    inquire ( file='output/restart/collapsed_vel.dat', exist=exist_restart ) 
    inquire ( file='output/restart/collapsed_grid.dat', exist=exist_grid    ) 
  


    if ( (.not. exist_restart) .or. (.not. exist_grid)  ) then 
       write(*,*) 'restart_collapsed_data_field:: necessary restart files are missing ... ' 
       write(*,*) 'restart_collapsed_data_field:: exiting ...' 
       call mem_dealloc 
       call mpi_finalize ( ierr ) 
       stop 
    endif



    !--- 
    ! read prev grid 
    !--- 
       allocate ( neigh_i(0:n1-1)) 
       allocate ( neigh_j(1:n2  )) 
       allocate ( neigh_k(0:n3-1)) 

       open(unit=876+rank, file='output/restart/collapsed_grid.dat', form='formatted', status='old') 
    
       !------ 
       ! read x values 
       !-----
       read(876+rank,*) n1_old 
       allocate ( xprev(1:n1_old)) 
       do i=1,n1_old 
          read(876+rank,*) xprev(i) 
       enddo
       
       
       !------ 
       ! read y values 
       !------ 
       read(876+rank,*) n2_old 
       allocate ( yprev(1:n2_old)) 
       do j=1,n2_old 
          read(876+rank,*) yprev(j) 
       enddo
       
    
       !----- 
       ! read z values 
       !----- 
       read(876+rank,*) n3_old 
       allocate ( zprev(1:n3_old)) 
       do k=1,n3_old 
          read(876+rank,*) zprev(k) 
       enddo
       
       close(876+rank)  


       hx = (chlx*pi) / dble(n1) 
       do i=1,n1 
          xx = ( dble(i) - 0.5d0 ) * hx 
          neigh_i(i-1) = iminloc ( abs(xprev - xx) ) - 1 
       enddo


       do j=1,n2 
          neigh_j(j) = iminloc ( abs(yprev - y(i_u,j)) ) 
       enddo


       hz = (chlz*pi) / dble(n3) 
       do k=1,n3 
          zz = ( dble(k) - 0.5d0 ) * hz 
          neigh_k(k-1) = iminloc ( abs(zprev - zz) ) - 1 
       enddo
       

    ! read prev vel field
    open(unit=499+rank, file='output/restart/collapsed_vel.dat', form='unformatted', status='old') 


    read(499+rank) n1_old 
    read(499+rank) n2_old 
    read(499+rank) n3_old 

    allocate ( qold(1:nvel,0:n1_old-1, 1:n2_old, 0:n3_old-1)) 


    qold = 0.d0 

    do k=0,n3_old-1 
    do j=1,n2_old 
    do i=0,n1_old-1
    do n=1,nvel
       
       read(499+rank) qold(n,i,j,k) 
    enddo 
    enddo 
    enddo 
    enddo


    close(499+rank) 




    ! read prev pressure field

    open(unit=1610+rank, file='output/restart/collapsed_p.dat', form='unformatted', status='old')
    read(1610+rank) n1_old
    read(1610+rank) n2_old
    read(1610+rank) n3_old

    allocate ( pgrid(1:n1_old, 1:n2_old, 1:n3_old))

    pgrid = 0.d0

    do k=1,n3_old
    do j=1,n2_old
    do i=1,n1_old
       read(1610+rank) pgrid(i,j,k)
    enddo
    enddo
    enddo

    close(1610+rank)


!    !--------------------------
!    ! write previous grid onto the current 
!    ! grid
!    !-------------------------
!
!       do k=ks,ke 
!       do j=1,n2 
!       do i=is,ie 
!
!
!          if ( j .eq. 1 ) then 
!             q(i_u,i,j,k) = y(i_u,1) 
!             q(i_v,i,j,k) = 0.d0 
!             q(i_w,i,j,k) = 0.d0 
!             
!          elseif  ( j .eq. n2 ) then 
!             q(i_u,i,j,k) = y(i_u,n2) 
!             q(i_v,i,j,k) = 0.d0 
!             q(i_v,i,j-1,k) = 0.d0 
!             q(i_w,i,j,k) = 0.d0 
!
!          else  
!          
!             do n=1,nvel
!                q(n,i,j,k) = qold(n, neigh_i(i), neigh_j(j), neigh_k(k) ) 
!             enddo
!          endif
!          
!       enddo
!       enddo 
!       enddo 
!
!
       !---- 
       ! use the old version of nearest neigh 
       ! intrpolation 
       !---- 

    
       if ( (int(n1/2) .eq. n1_old) .and. (int(n3/2) .eq. n3_old) .and. & 
            (int(n2/2) .eq. n2_old) ) then 
          
          do k=ks,ke 
          do j=1,n2 
          do i=is,ie 
                   
             kk = int ( k/2) 
             jj = int ( j/2) 
             ii = int ( i/2) 

             if ( j .eq. 1 ) then 
                do n=1,nvel 
                   q(n,i,j,k)  = 0.5d0 * qold(n,ii,1,kk) 
                enddo
         
             elseif  ( j .eq. n2 ) then 
                do n=1,nvel 
                   q(n,i,j,k) = 0.5d0 * qold(n,ii,n2/2,kk)
                enddo

             else 

                if ( mod(j,2) .eq. 0 ) then 
                   do n=1,nvel
                      q(n,i,j,k) = qold(n,ii, jj, kk) 
                   enddo

                else  
                   do n=1,nvel
                      q(n,i,j,k) = 0.5d0 * ( qold(n,ii,jj, kk) + qold(n,ii,jj+1,kk)) 
                   enddo
                endif
             endif
          enddo
          enddo 
          enddo 


       elseif ( (int(n3/2) .eq. n3_old) .and. (n1 .eq. n1_old) .and. & 
                (n2 .eq. n2_old)) then 


          do k=ks,ke 
          do j=1,n2 
          do i=is,ie 

             kk = int(k/2) 
             jj = j 
             ii = i 


             do n=1,nvel 
                q(n,i,j,k) = qold(n,ii,jj,kk)
             enddo
          enddo
          enddo 
          enddo 


       elseif ( (int(n3/2) .eq. n3_old) .and. (int(n1/2) .eq. n1_old) .and. &
                (n2 .eq. n2_old)) then

          do k=ks,ke
          do j=1,n2
          do i=is,ie

             kk = int(k/2)
             jj = j
             ii = int(i/2)

             do n=1,nvel
                q(n,i,j,k) = qold(n,ii,jj,kk)
             enddo
                p_old(i,j,k) = pgrid(ii,jj,kk)

          enddo
          enddo
          enddo


       elseif ( (n3 .eq. int(n3_old/2)) .and. (n1 .eq. int(n1_old/2)) .and. &
                (n2 .eq. n2_old) ) then

          do k=ks,ke
          do j=1,n2
          do i=is,ie

             kk = k*2
             jj = j
             ii = i*2

             do n=1,nvel
                q(n,i,j,k) = qold(n,ii,jj,kk)
             enddo
                p_old(i,j,k) = pgrid(ii,jj,kk)

          enddo
          enddo
          enddo


       else 

          do k=ks,ke 
          do j=1,n2 
          do i=is,ie 
          do n=1,nvel

             q(n,i,j,k) = qold(n,i,j,k) 
          enddo
          enddo 
          enddo 
          enddo 

       endif


    !--- 
    ! cleanup mem
    !---
    
    deallocate ( qold   ) 
    deallocate ( pgrid  )

    if ( allocated(xprev)) deallocate ( xprev ) 
    if ( allocated(yprev)) deallocate ( yprev ) 
    if ( allocated(zprev)) deallocate ( zprev ) 
    if ( allocated(neigh_i)) deallocate ( neigh_i) 
    if ( allocated(neigh_j)) deallocate ( neigh_j) 
    if ( allocated(neigh_k)) deallocate ( neigh_k) 


  end subroutine restart_interface_collapsed_data_field 

!================================================================================================================================!

  subroutine vortex_wall_plane

    implicit none
    real, dimension(9) :: du
    real :: pp, qq, phi, hx, hz
    integer :: i,j,k, y_plane
    real, dimension(3, is:ie, 1:n2, ks:ke) :: vort
  
    pi = acos(-1.0d0)
    du  = 0.d0
    hx = (chlx*pi ) / dble(n1)
    hz = (chlz*pi ) / dble(n3)
    y_plane = 3

    do k=ks+1,ke+1
    do j=0,y_plane
    do i=is+1,ie+1

       if ( j .eq. 0 ) then
          du(1) = ( q(i_u,i,j+1,k)-q(i_u,i,j,k))/ ( y(i_u,j+1)-y(i_u,j)) ! uy
          du(6) = ( q(i_w,i,j+1,k)-q(i_w,i,j,k))/ ( y(i_u,j+1)-y(i_u,j)) ! wy 
          du(7) = ( q(i_v,i,j+1,k)-q(i_v,i,j,k))/ ( y(i_v,j+1)-y(i_v,j)) ! vy 
       else
          du(1) = ( q(i_u,i,j+1,k)-q(i_u,i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! uy
          du(6) = ( q(i_w,i,j+1,k)-q(i_w,i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! wy
          du(7) = ( q(i_v,i,j+1,k)-q(i_v,i,j-1,k)) / ( y(i_v,j+1)-y(i_v,j-1)) ! vy
       endif

       du(2) =  ( q(i_u,i,j,k+1)-q(i_u,i,j,k-1))/2.d0/hz   ! uz
       du(3) =  ( q(i_v,i+1,j,k)-q(i_v,i-1,j,k))/2.d0/hx   ! vx
       du(4) =  ( q(i_v,i,j,k+1)-q(i_v,i,j,k-1))/2.d0/hz   ! vz
       du(5) =  ( q(i_w,i+1,j,k)-q(i_w,i-1,j,k))/2.d0/hx   ! wx 
       du(8) =  ( q(i_w,i,j,k+1)-q(i_w,i,j,k-1))/2.d0/hz   ! wz 
       du(9) =  ( q(i_u,i+1,j,k)-q(i_u,i-1,j,k))/2.d0/hx   ! wx 

       vort(1,i-1,j,k-1) = du(6) - du(4)
       vort(2,i-1,j,k-1) = du(2) - du(5)
       vort(3,i-1,j,k-1) = du(3) - du(1)
      
       vort_wall_plane(i_u,j,i-1,k-1) = vort(i_u,i-1,j,k-1)
       vort_wall_plane(i_v,j,i-1,k-1) = vort(i_v,i-1,j,k-1)
       vort_wall_plane(i_w,j,i-1,k-1) = vort(i_w,i-1,j,k-1)

    enddo
    enddo
    enddo


  end subroutine vortex_wall_plane

!================================================================================================================================!
  subroutine dump_vorticies_wall_plane

    implicit none
    real, dimension(:,:,:,:)    , allocatable :: vort_print
    real, dimension(:)        , allocatable :: recvbuf , rbuf_eta, rbuf_p, rbuf_u, rbuf_v
    integer                                 :: i,j,k,ivel
    integer                                 :: itblock, rshft, n, proc, count, ind
    integer, dimension(:), allocatable      :: is_vec, ks_vec
    logical                                 :: debug_collapse
    logical                                 :: write_collapsed_grid
    real                                    :: hx, hz, xx, zz, u, v,w
    character(len=4)                        :: nstr
    character(len=80)                 :: fname
    character(len=3)                  :: time_str, time_str_end, time_str_i
    character(len=4)                        :: time_str_ibp


    !---------------------
    ! gather the start and ending indices for 
    ! all processors to master 
    ! currently assumes that all the blocks are of the same size
    !--------------------

     write(nstr, '(i4.4)') n_interface
    fname = 'output/flow_field/mean_vort_wall_'//nstr//'.dat'

    allocate ( is_vec(0:size-1),ks_vec(0:size-1) )

    call mpi_gather( is, 1, mpi_integer, is_vec(0), 1, mpi_integer, 0, mycomm, ierr)
    call mpi_gather( ks, 1, mpi_integer, ks_vec(0), 1, mpi_integer, 0, mycomm, ierr)

    nz = int ( n3 / p3 )
    nx = int ( n1 / p1 )

    itblock  = nx * nz * 3 * 4
    allocate ( recvbuf( size*itblock) )

    !--------------------------
    rshft = rank * itblock + 1
    call mpi_gather ( vort_wall_plane(1,0,is,ks), itblock, MPI_DOUBLE_PRECISION, &
                      recvbuf(1), itblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)

    allocate ( vort_print(3, 0:3, 0:n1-1, 0:n3-1))

   if ( rank .eq. 0 ) then
       do proc=0,size-1

          rshft = proc * itblock +1
          count = 0
          do k= ks_vec(proc), ks_vec(proc)+nz-1
          do i= is_vec(proc), is_vec(proc)+nx-1
          do j= 0,3
              do ivel = 1,3
             ind = rshft+count
             if ( k .ge. 0 .and. k .le. n3-1 .and. &
                  i .ge. 0 .and. i .le. n1-1  ) then
                vort_print(ivel,j,i,k) = recvbuf(rshft+count)
             endif

             count = count + 1
              enddo
          enddo
          enddo
          enddo
       enddo

      write(nstr, '(i4.4)') n_interface
       write(time_str_ibp, '(i4.4)') int((time- int(time))*10000)
       write(time_str_i,'(i3.3)') int(time)
       open(unit=1503,file='output/interface/vort'//time_str_i//'_'//time_str_ibp//'.dat',form='unformatted',status='unknown')
       !fname = 'output/interface/mean_vort_wall_'//nstr//'.dat'
       !open(unit=1503,file=fname, form='unformatted', status='unknown')
       write(1503) n1,n3
       write(1503) vort_print
       close(1503)

   endif

    deallocate ( recvbuf )
    deallocate ( vort_print)
    deallocate ( is_vec  )
    deallocate ( ks_vec  )


  end subroutine dump_vorticies_wall_plane

!================================================================================================================================!

  subroutine vortex_wy

    implicit none
    real, dimension(9) :: du
    real :: pp, qq, phi, hx, hz
    integer :: i,j,k
    real, dimension(:,:,:), allocatable :: uc, vc, wc
    real, dimension(3, is:ie, 1:n2, ks:ke) :: vort
  
    allocate (uc(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp,ks-halo_z+1:ke+halo_z+1))
    allocate (vc(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp,ks-halo_z+1:ke+halo_z+1))
    allocate (wc(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp,ks-halo_z+1:ke+halo_z+1))

    pi = acos(-1.0d0)
    du  = 0.d0  
    hx = (chlx*pi ) / dble(n1)
    hz = (chlz*pi ) / dble(n3)

    do k= ks-2, ke+2
    do j= 1,5
    do i= is-2, ie+2

       uc(i,j,k)=  0.5d0*( q(i_u,i,j,k)   + q(i_u,i+1,j,k))
       vc(i,j,k) = 0.5d0*( q(i_v,i,j-1,k) + q(i_v,i,j,k))
       wc(i,j,k) = 0.5d0*( q(i_w,i,j,k)   + q(i_w,i,j,k+1))

    enddo
    enddo
    enddo


    do k=ks+1,ke+1
    do j=1,4
    do i=is+1,ie+1

       if ( j .eq. 1 ) then
          du(1) = ( q(i_u,i,j,k)-q(i_u,i,j-1,k))/ ( y(i_u,j)-y(i_u,j-1))/2.d0 ! uy
          du(6) = ( q(i_w,i,j,k)-q(i_w,i,j-1,k))/ ( y(i_u,j)-y(i_u,j-1))/2.d0 ! wy 
          du(7) = ( q(i_v,i,j,k)-q(i_v,i,j-1,k))/ ( y(i_v,j)-y(i_v,j-1)) ! vy 
       else
          du(1) = ( uc(i,j+1,k) - uc(i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! uy
          du(6) = ( wc(i,j+1,k) - wc(i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! wy
          du(7) = ( vc(i,j+1,k) - vc(i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! vy
       endif

       du(2) =  ( uc(i,j,k+1) - uc(i,j,k-1))/2.d0/hz   ! uz
       du(3) =  ( vc(i+1,j,k) - vc(i-1,j,k))/2.d0/hx   ! vx
       du(4) =  ( vc(i,j,k+1) - vc(i,j,k-1))/2.d0/hz   ! vz
       du(5) =  ( wc(i+1,j,k) - wc(i-1,j,k))/2.d0/hx   ! wx 
       du(8) =  ( wc(i,j,k+1) - wc(i,j,k-1))/2.d0/hz   ! wz 
       du(9) =  ( uc(i+1,j,k) - uc(i-1,j,k))/2.d0/hx   ! wx 

       vort(1,i-1,j,k-1) = du(6) - du(4)
       vort(2,i-1,j,k-1) = du(2) - du(5)
       vort(3,i-1,j,k-1) = du(3) - du(1)

       vort_wall_plane(i_u,j-1,i-1,k-1) = vort(i_u,i-1,j,k-1)
       vort_wall_plane(i_v,j-1,i-1,k-1) = vort(i_v,i-1,j,k-1)
       vort_wall_plane(i_w,j-1,i-1,k-1) = vort(i_w,i-1,j,k-1)

    enddo
    enddo
    enddo

  end subroutine vortex_wy
!================================

  subroutine vortex_wall

    implicit none
    real, dimension(9) :: du
    real :: pp, qq, phi, hx, hz
    integer :: i,j,k
    real, dimension(:,:,:), allocatable :: uc, vc, wc
    real, dimension(3, is:ie, 1:n2, ks:ke) :: vort
  
    allocate (uc(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp,ks-halo_z+1:ke+halo_z+1))
    allocate (vc(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp,ks-halo_z+1:ke+halo_z+1))
    allocate (wc(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp,ks-halo_z+1:ke+halo_z+1))

    pi = acos(-1.0d0)
    du  = 0.d0  
    hx = (chlx*pi ) / dble(n1)
    hz = (chlz*pi ) / dble(n3)

    do k= ks-2, ke+2
    do j= 1,5
    do i= is-2, ie+2

       uc(i,j,k)=  0.5d0*( q(i_u,i,j,k)   + q(i_u,i+1,j,k))
       vc(i,j,k) = 0.5d0*( q(i_v,i,j-1,k) + q(i_v,i,j,k))
       wc(i,j,k) = 0.5d0*( q(i_w,i,j,k)   + q(i_w,i,j,k+1))

    enddo
    enddo
    enddo


    do k=ks+1,ke+1
    do j=1,4
    do i=is+1,ie+1

       if ( j .eq. 1 ) then
          du(1) = ( q(i_u,i,j,k)-q(i_u,i,j-1,k))/ ( y(i_u,j)-y(i_u,j-1))/2.d0 ! uy
          du(6) = ( q(i_w,i,j,k)-q(i_w,i,j-1,k))/ ( y(i_u,j)-y(i_u,j-1))/2.d0 ! wy 
          du(7) = ( q(i_v,i,j,k)-q(i_v,i,j-1,k))/ ( y(i_v,j)-y(i_v,j-1)) ! vy 
       else
          du(1) = ( uc(i,j+1,k) - uc(i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! uy
          du(6) = ( wc(i,j+1,k) - wc(i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! wy
          du(7) = ( vc(i,j+1,k) - vc(i,j-1,k)) / ( y(i_u,j+1)-y(i_u,j-1)) ! vy
       endif

       du(2) =  ( uc(i,j,k+1) - uc(i,j,k-1))/2.d0/hz   ! uz
       du(3) =  ( vc(i+1,j,k) - vc(i-1,j,k))/2.d0/hx   ! vx
       du(4) =  ( vc(i,j,k+1) - vc(i,j,k-1))/2.d0/hz   ! vz
       du(5) =  ( wc(i+1,j,k) - wc(i-1,j,k))/2.d0/hx   ! wx 
       du(8) =  ( wc(i,j,k+1) - wc(i,j,k-1))/2.d0/hz   ! wz 
       du(9) =  ( uc(i+1,j,k) - uc(i-1,j,k))/2.d0/hx   ! wx 

       vort(1,i-1,j,k-1) = du(6) - du(4)
       vort(2,i-1,j,k-1) = du(2) - du(5)
       vort(3,i-1,j,k-1) = du(3) - du(1)

       vort_wall_plane(i_u,j-1,i-1,k-1) = vort(i_u,i-1,j,k-1)
       vort_wall_plane(i_v,j-1,i-1,k-1) = vort(i_v,i-1,j,k-1)
       vort_wall_plane(i_w,j-1,i-1,k-1) = vort(i_w,i-1,j,k-1)

    enddo
    enddo
    enddo

  end subroutine vortex_wall

!================================================================================================================================!
  subroutine dump_vorticies_wall

    implicit none
    real, dimension(:,:,:)    , allocatable :: vort_print
    real, dimension(:)        , allocatable :: recvbuf , rbuf_eta, rbuf_p, rbuf_u, rbuf_v
    integer                                 :: i,j,k
    integer                                 :: itblock, rshft, n, proc, count, ind
    integer, dimension(:), allocatable      :: is_vec, ks_vec
    logical                                 :: debug_collapse
    logical                                 :: write_collapsed_grid
    real                                    :: hx, hz, xx, zz, u, v,w
    character(len=4)                        :: nstr
    character(len=80)                 :: fname
    character(len=3)                  :: time_str, time_str_end, time_str_i
    character(len=4)                        :: time_str_ibp


    !---------------------
    ! gather the start and ending indices for 
    ! all processors to master 
    ! currently assumes that all the blocks are of the same size
    !--------------------

    allocate ( is_vec(0:size-1),ks_vec(0:size-1) )

    call mpi_gather( is, 1, mpi_integer, is_vec(0), 1, mpi_integer, 0, mycomm, ierr)
    call mpi_gather( ks, 1, mpi_integer, ks_vec(0), 1, mpi_integer, 0, mycomm, ierr)

    nz = int ( n3 / p3 )
    nx = int ( n1 / p1 )

    itblock  = nx * nz * 3
    allocate ( recvbuf( size*itblock) )

    !--------------------------
    rshft = rank * itblock + 1
    call mpi_gather ( vort_wall(1,is,ks), itblock, MPI_DOUBLE_PRECISION, &
                      recvbuf(1), itblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)

    allocate ( vort_print(3, 0:n1-1, 0:n3-1))

   if ( rank .eq. 0 ) then
       do proc=0,size-1

          rshft = proc * itblock +1
          count = 0
          do k= ks_vec(proc), ks_vec(proc)+nz-1
          do i= is_vec(proc), is_vec(proc)+nx-1
              do j = 1,3
             ind = rshft+count
             if ( k .ge. 0 .and. k .le. n3-1 .and. &
                  i .ge. 0 .and. i .le. n1-1  ) then
                vort_print(j,i,k) = recvbuf(rshft+count)
             endif

             count = count + 1
              enddo
          enddo
          enddo
       enddo

      write(nstr, '(i4.4)') n_interface
!      fname = 'output/interface/mean_vort_wall_'//nstr//'.dat'
       write(time_str_ibp, '(i4.4)') int((time- int(time))*10000)
       write(time_str_i,'(i3.3)') int(time)
       open(unit=1503,file='output/interface/vort'//time_str_i//'_'//time_str_ibp//'.dat',form='unformatted',status='unknown')


       open(unit=1503,file=fname, form='unformatted', status='unknown')
       write(1503) n1,n3
       write(1503) vort_print
       close(1503)

   endif

    deallocate ( recvbuf )
    deallocate ( vort_print)
    deallocate ( is_vec  )
    deallocate ( ks_vec  )


  end subroutine dump_vorticies_wall

!================================================================================================================================!
  subroutine dump_2D_slice_wall

    !-----------------------------------------------------------------
    ! collapse entire velocity and pressure fields onto the master node and output one file 
    ! currently assumes that all subdomains are the same size
    !-----------------------------------------------------------------

    implicit none
    real, dimension(:,:,:)    , allocatable :: qq
    real, dimension(:,:)    , allocatable :: pp
    real, dimension(:)        , allocatable :: recvbuf , rbuf_p
    integer                                 :: i,j,k
    integer                                 :: qblock , rshft, n, proc, count, ind, pblock
    integer, dimension(:), allocatable      :: ks_vec, is_vec
    logical                                 :: debug_collapse
    logical                                 :: write_collapsed_grid
    real                                    :: hx, hz, xx, zz, u, v,w
    character(len=3)                        :: nstr
    integer                           :: ntime
    character(len=80)                 :: fname
    character(len=3)                  :: time_str, time_str_end


    write(nstr, '(i3.3)') n_interface
    fname = 'output/flow_field/class_spec'//nstr//'.dat'

    allocate ( ks_vec(0:size-1), is_vec(0:size-1))

    call mpi_gather( ks, 1, mpi_integer, ks_vec(0), 1, mpi_integer, 0, mycomm, ierr)
    call mpi_gather( is, 1, mpi_integer, is_vec(0), 1, mpi_integer, 0, mycomm, ierr)

    nx = int ( n1 / p1 )
    nz = int ( n3 / p3 )

    qblock  = nvel*(nx + 2*halo_x )* (nz + 2*halo_z )* (n2 + 2*halo_yu )
    allocate ( recvbuf( size*qblock) )
    !--------------------------
    ! gather q from all processors to master
    rshft = rank * qblock + 1
    call mpi_gather ( q(1,is-halo_x,1-halo_yu,ks-halo_z), qblock, MPI_DOUBLE_PRECISION, &
                      recvbuf(1), qblock, MPI_DOUBLE_PRECISION, 0, mycomm, ierr)

    allocate ( qq(1:nvel,0:n1-1,0:n3-1))

    !--------------------------
    ! master node organizes the received data 
    !-------------------------

    if ( rank .eq. 0 ) then

       do proc=0,size-1

          !-----------------
          ! parse recv buff and extract values 
          ! into qq 
          !-----------------
          rshft = proc * qblock +1
          count = 0

          do k= ks_vec(proc)-halo_z, ks_vec(proc)+nz-1+halo_z
          do j= 1-halo_yu, n2+halo_yu
          do i= is_vec(proc)-halo_x, is_vec(proc)+nx-1+halo_x
          do n= 1,nvel

             ind = rshft+count
             if ( k .ge. 0 .and. k .le. n3-1 .and. &
                  i .ge. 0 .and. i .le. n1-1 .and. &
                  j .eq. 30 ) then
                qq(n,i,k) = recvbuf(rshft+count)
             endif

             count = count + 1
          enddo
          enddo
          enddo
          enddo
      enddo


       open(unit=1499,file=fname, form='unformatted', status='unknown')
       write(1499) n1,n3
       write(1499) qq
       close(1499)

    endif


    call mpi_barrier ( mycomm, ierr )


    !--- 
    ! temporary cleanup 
    !---
    deallocate ( recvbuf)
    deallocate ( qq     )
    deallocate ( ks_vec  )
    deallocate ( is_vec  )

  end subroutine dump_2D_slice_wall

!======================================================================================================!
  subroutine restart_collapsed_data_field 

    implicit none 
    integer                              :: n1_old, n2_old, n3_old, nneigh_i, nneigh_j, nneigh_k
    integer, dimension(:), allocatable   :: neigh_i , neigh_j, neigh_k 
    real, dimension(:), allocatable      :: xprev, yprev, zprev 
    real, dimension(:,:,:,:),allocatable :: qold, prtrb 
    real, dimension(:,:,:)  ,allocatable :: pgrid
    real, dimension(1:nvel)              :: rms, qsum , gqsum 
    real                                 :: amp, urc, vrc, wrc, fp, yp , xx, zz, hx, hz, tfax
    integer, dimension(4)                :: iseed
    integer                              :: i,j,k,n,i_q, ii, jj, kk
    logical                              :: exist_restart, exist_grid 
    integer                              :: use_deprecated_restart 


    !XXXX should move some of these options to an input file at some pt 
    !XXXX


    ! use this if the prev field is just a simple refine homothetically
    use_deprecated_restart = 0

    !---- 
    ! verify necessary restart files exist 
    !---- 
    inquire ( file='output/collapsed_restart.dat', exist=exist_restart ) 
    if ( use_deprecated_restart ) then 
       exist_grid = .true. 
    else  
       inquire ( file='output/collapsed_grid.dat'   , exist=exist_grid    ) 
    endif 

    if ( (.not. exist_restart) .or. (.not. exist_grid) ) then 
       write(*,*) 'restart_collapsed_data_field:: necessary restart files are missing ... ' 
       write(*,*) 'restart_collapsed_data_field:: exiting ...' 
       call mem_dealloc 
       call mpi_finalize ( ierr ) 
       stop 
    endif


    !------
    ! set the amplitudes if you want to reperturb the soln
    !------

    amp = 0.00d0   ! add slight perturb...
    urc = 5.1d0 
    vrc = 1.0d0 
    wrc = 2.3d0 



    !--- 
    ! read prev grid 
    !--- 
    if ( use_deprecated_restart .eq. 0 ) then 

       allocate ( neigh_i(0:n1-1)) 
       allocate ( neigh_j(1:n2  )) 
       allocate ( neigh_k(0:n3-1)) 

       open(unit=876+rank, file='output/collapsed_grid.dat', form='formatted', status='old') 
    
       !------ 
       ! read x values 
       !-----
       read(876+rank,*) n1_old 
       allocate ( xprev(1:n1_old)) 
       do i=1,n1_old 
          read(876+rank,*) xprev(i) 
       enddo
       
       
       !------ 
       ! read y values 
       !------ 
       read(876+rank,*) n2_old 
       allocate ( yprev(1:n2_old)) 
       do j=1,n2_old 
          read(876+rank,*) yprev(j) 
       enddo
       
    
       !----- 
       ! read z values 
       !----- 
       read(876+rank,*) n3_old 
       allocate ( zprev(1:n3_old)) 
       do k=1,n3_old 
          read(876+rank,*) zprev(k) 
       enddo
       
       close(876+rank)  


       !---- 
       ! compute the nearest neighbor mapping 
       ! from old grid to new 
       !   search is independent in each direction 
       !---- 
    
       hx = (chlx*pi) / dble(n1) 
       do i=1,n1 
          xx = ( dble(i) - 0.5d0 ) * hx 
          neigh_i(i-1) = iminloc ( abs(xprev - xx) ) - 1 
       enddo


       do j=1,n2 
          neigh_j(j) = iminloc ( abs(yprev - y(i_u,j)) ) 
       enddo


       hz = (chlz*pi) / dble(n3) 
       do k=1,n3 
          zz = ( dble(k) - 0.5d0 ) * hz 
          neigh_k(k-1) = iminloc ( abs(zprev - zz) ) - 1 
       enddo
       
    endif

    !---- 
    ! read prev vel field
    !---- 

    open(unit=499+rank, file='output/collapsed_restart.dat', form='formatted', status='old') 
    read(499+rank,*) n1_old 
    read(499+rank,*) n2_old 
    read(499+rank,*) n3_old 

    !-------------------------
    !  read old vel field 
    !  XXXX this approach will not work > 10 mil pts, and should be fixed.
    !-------------------------

    allocate ( qold(1:nvel,0:n1_old-1, 0:n2_old+1, 0:n3_old-1)) 
    allocate ( prtrb(1:nvel,is:ie,1:n2,ks:ke))

    qold = 0.d0 

    do k=0,n3_old-1 
    do j=1,n2_old 
    do i=0,n1_old-1
    do n=1,nvel
       
       read(499+rank,*) qold(n,i,j,k) 
    enddo 
    enddo 
    enddo 
    enddo

    close(499+rank) 

    !--------------------------
    ! write previous grid onto the current 
    ! grid
    !-------------------------

    if ( use_deprecated_restart .eq. 0 ) then 

       do k=ks,ke 
       do j=1,n2 
       do i=is,ie 


      !   if ( j .eq. 1 ) then 
      !       q(i_u,i,j,k) = y(i_u,1) 
      !       q(i_v,i,j,k) = 0.d0 
      !       q(i_w,i,j,k) = 0.d0 
             
     !     elseif  ( j .eq. n2 ) then 
     !        q(i_u,i,j,k) = y(i_u,n2) 
     !        q(i_v,i,j,k) = 0.d0 
     !        q(i_v,i,j-1,k) = 0.d0 
     !        q(i_w,i,j,k) = 0.d0 

     !     else  
          
             do n=1,nvel
                q(n,i,j,k) = qold(n, neigh_i(i), neigh_j(j), neigh_k(k) ) 
             enddo
     !     endif
          
       enddo
       enddo 
       enddo 


    endif 


    if ( use_deprecated_restart .eq. 1 ) then 
       
       !---- 
       ! use the old version of nearest neigh 
       ! intrpolation 
       !---- 

    
       if ( (int(n1/2) .eq. n1_old) .and. (int(n3/2) .eq. n3_old) .and. & 
            (int(n2/2) .eq. n2_old) ) then 
          
          do k=ks,ke 
          do j=1,n2 
          do i=is,ie 
                   
             kk = int ( k/2) 
             jj = int ( j/2) 
             ii = int ( i/2) 

             if ( j .eq. 1 ) then 
                do n=1,nvel 
                   q(n,i,j,k)  = 0.5d0 * qold(n,ii,1,kk) 
                enddo
         
             elseif  ( j .eq. n2 ) then 
                do n=1,nvel 
                   q(n,i,j,k) = 0.5d0 * qold(n,ii,n2/2,kk)
                enddo

             else 

                if ( mod(j,2) .eq. 0 ) then 
                   do n=1,nvel
                      q(n,i,j,k) = qold(n,ii, jj, kk) 
                   enddo

                else  
                   do n=1,nvel
                      q(n,i,j,k) = 0.5d0 * ( qold(n,ii,jj, kk) + qold(n,ii,jj+1,kk)) 
                   enddo
                endif
             endif
          enddo
          enddo 
          enddo 


       elseif ( (int(n3/2) .eq. n3_old) .and. (n1 .eq. n1_old) .and. & 
                (n2 .eq. n2_old)) then 


          do k=ks,ke 
          do j=1,n2 
          do i=is,ie 

             kk = int(k/2) 
             jj = j 
             ii = i 


             do n=1,nvel 
                q(n,i,j,k) = qold(n,ii,jj,kk)
             enddo
          enddo
          enddo 
          enddo 




       else 

          do k=ks,ke 
          do j=1,n2 
          do i=is,ie 
          do n=1,nvel

             q(n,i,j,k) = qold(n,i,j,k) 
          enddo
          enddo 
          enddo 
          enddo 

       endif

    endif


    !-----------
    ! USE_DEPRECATED_RESTART = 2 
    !    replicated the current data, assuming resolution is fixed
    !    and the new size is a multiple of the current size
    !-----------
    if ( use_deprecated_restart .eq. 2 ) then 

       do k=ks,ke 
       do j=1,n2 
       do i=is,ie 
       do n=1,nvel

          ii = mod(i,n1_old) 
          jj = mod(j,n2_old+1) 
          kk = mod(k,n3_old)

          tfax = 1.d0
          !tfax = 0.97d0 + 0.03d0* sqrt( 1.d0- y(i_u,j)*y(i_u,j))
          q(n,i,j,k) = tfax* qold(n,ii,jj,kk) 
       enddo
       enddo
       enddo
       enddo

    endif
          


    !---------
    ! add noise onto old soln 
    ! crude approx 
    !---------

    !iseed(1) = 7 
    !iseed(2) = 563
    !iseed(3) = 777 
    !iseed(4) = 3099

    iseed(1) = 2* rank + 3* mod(200,rank+1) + 7 
    iseed(2) = 5* rank + 12* mod(12, rank+1) + 13
    iseed(3) = rank* rank + 3
    iseed(4) = 8*rank +1

    do k=ks,ke
    do j=1,n2 
    do i=is,ie
    do n=1,nvel 

       call dlarnv( 1,iseed, 1, prtrb(n,i,j,k))
    enddo 
    enddo 
    enddo 
    enddo 


    !---- 
    ! subtract mean from the noise 
    ! in each wall normal plane 
    !----

    do j=1,n2
       
       qsum = 0.d0 
    
       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel 

          qsum(i_q) = qsum(i_q) + prtrb(i_q,i,j,k)
          
       enddo
       enddo
       enddo


       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr) 
       gqsum = gqsum / dble(n1) / dble(n3) 
       
       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel 

          prtrb(i_q, i,j,k) = prtrb(i_q, i,j,k) - gqsum(i_q) 
       enddo
       enddo
       enddo

    enddo
       





    do j=1,n2 

       yp = re* (1.d0 - abs(y(i_u,j)))
       fp = 1.d0 - exp(-yp/25.d0)         
       rms(i_u) = sqrt(urc * fp*fp*abs(y(i_u,j)) )
       rms(i_w) = sqrt(wrc * fp*fp*abs(y(i_u,j)) )

       yp = re* (1.d0 - abs(y(i_v,j))) 
       fp = 1.d0 - exp(-yp/25.d0) 
       rms(i_v) = sqrt(vrc * fp*fp*abs(y(i_v,j)) )

       qsum     = 0.d0 
       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel
          
          qsum(i_q) = qsum(i_q) + prtrb(i_q,i,j,k)*prtrb(i_q,i,j,k) 
       enddo
       enddo
       enddo


       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr) 
       gqsum = gqsum / dble(n1) /dble(n3) 

       do i_q=1,nvel 
          gqsum(i_q) = amp*rms(i_q)* sqrt(1.d0/gqsum(i_q) )
       enddo



       do k=ks,ke
       do i=is,ie
       do i_q=1,nvel 

          prtrb(i_q,i,j,k) = prtrb(i_q,i,j,k)*gqsum(i_q) 
       enddo
       enddo
       enddo

    enddo


    !----
    ! add perturbation to vel field 
    !---- 

    do k=ks,ke
    do j=1,n2 
    do i=is,ie
    do n=1,nvel 

       q(n,i,j,k) = q(n,i,j,k) + prtrb(n,i,j,k) 
    enddo 
    enddo 
    enddo 
    enddo 

    !--- 
    ! cleanup mem
    !---
    
    deallocate ( prtrb ) 
    deallocate ( qold  ) 
    if ( allocated(xprev)) deallocate ( xprev ) 
    if ( allocated(yprev)) deallocate ( yprev ) 
    if ( allocated(zprev)) deallocate ( zprev ) 
    if ( allocated(neigh_i)) deallocate ( neigh_i) 
    if ( allocated(neigh_j)) deallocate ( neigh_j) 
    if ( allocated(neigh_k)) deallocate ( neigh_k) 



    ! dont load the old pressure field , no real savings...
    p_old = 0.d0

  end subroutine restart_collapsed_data_field 


!================================================================================================================================!

  subroutine add_perturb_to_field

    implicit none
    real, dimension(1:nvel)              :: rms, qsum , gqsum
    real                                 :: amp, urc, vrc, wrc, fp, yp , xx, zz, hx, hz, tfax
    integer, dimension(4)                :: iseed
    integer                              :: i,j,k,n,i_q, ii, jj, kk
    real, dimension(:,:,:,:),allocatable :: prtrb



    amp = 10.00d0  ! add slight perturb...
    urc = 2.1d0
    vrc = 1.0d0
    wrc = 2.3d0

    allocate ( prtrb(1:nvel,is:ie,1:n2,ks:ke))

    !---------
    ! add noise onto old soln
    ! crude approx
    !---------

    !iseed(1) = 7
    !iseed(2) = 563
    !iseed(3) = 777
    !iseed(4) = 3099

    iseed(1) = 7 + 2*rank
    iseed(2) = 563 + rank
    iseed(3) = 777
    iseed(4) = 3099 - 3*rank

    ! base the seeds on the rank so that the perturb is not the
    ! same in every block

!    iseed(1) = 2* rank + 3* mod(200,rank+1) + 7
!    iseed(2) = 5* rank + 12* mod(12, rank+1) + 13
!    iseed(3) = rank* rank + 3
!    iseed(4) = 8*rank +1

    ! crash check
!    do n=1,4
!       if ( iseed(n) .lt. 0 .or. iseed(n) .gt. 4095) then
!          write(*,*) 'ERROR IN SEED = ', rank, n, iseed(n)
!       endif
!    enddo

    do k=ks,ke
    do j=1,n2
    do i=is,ie
    do n=1,nvel

       call dlarnv( 1,iseed, 1, prtrb(n,i,j,k))
    enddo
    enddo
    enddo
    enddo


    !----
    ! subtract mean from the noise
    ! in each wall normal plane
    !----

    do j=1,n2

       qsum = 0.d0

       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          qsum(i_q) = qsum(i_q) + prtrb(i_q,i,j,k)

       enddo
       enddo
       enddo


       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr)
       gqsum = gqsum / dble(n1) / dble(n3)

       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          prtrb(i_q, i,j,k) = prtrb(i_q, i,j,k) - gqsum(i_q)
       enddo
       enddo
       enddo

    enddo




    do j=1,n2

       yp = re* (1.d0 - abs(y(i_u,j)))
       fp = 1.d0 - exp(-yp/25.d0)
       rms(i_u) = sqrt(urc * fp*fp*abs(y(i_u,j)) )
       rms(i_w) = sqrt(wrc * fp*fp*abs(y(i_u,j)) )

       yp = re* (1.d0 - abs(y(i_v,j)))
       fp = 1.d0 - exp(-yp/25.d0)
       rms(i_v) = sqrt(vrc * fp*fp*abs(y(i_v,j)) )

       qsum     = 0.d0
       do k=ks,ke
       do i=is,ie
       do i_q = 1,nvel

          qsum(i_q) = qsum(i_q) + prtrb(i_q,i,j,k)*prtrb(i_q,i,j,k)
       enddo
       enddo
       enddo


       call mpi_allreduce( qsum, gqsum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr)
       gqsum = gqsum / dble(n1) /dble(n3)

       do i_q=1,nvel
          gqsum(i_q) = amp*rms(i_q)* sqrt(1.d0/gqsum(i_q) )
       enddo



       do k=ks,ke
       do i=is,ie
       do i_q=1,nvel

          prtrb(i_q,i,j,k) = prtrb(i_q,i,j,k)*gqsum(i_q)
       enddo
       enddo
       enddo

    enddo


    !----
    ! add perturbation to vel field
    !----

    do k=ks,ke
    do j=1,n2
    do i=is,ie
    do n=1,nvel

       q(n,i,j,k) = q(n,i,j,k) + prtrb(n,i,j,k)
    enddo
    enddo
    enddo
    enddo

    deallocate ( prtrb)

  end subroutine add_perturb_to_field

!================================================================================================================================!
  function iminloc ( arr ) 
    !---- 
    ! wrapped version of the intrinsic minloc 
    ! returns answer as integer 
    !  adapted from numerical recipes in f90 
    !----  
    implicit none 
    real, dimension(:), intent(in) :: arr 
    integer, dimension(1)   :: imin 
    integer                 :: iminloc 

    imin = minloc(arr(:)) 
    iminloc = imin(1) 

  end function iminloc

!====================================================================================================================================!


end module output
    
    

