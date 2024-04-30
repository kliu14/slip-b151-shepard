
module global 

!------------------------------------------------------! 
!*  variable declaration for 4th order channel        *! 
!*  code                                              *!
!*----------------------------------------------------*! 

  implicit none 
  include 'mpif.h'


  !----------------
  ! grid definition variables 
  !----------------

  integer                               :: n1      ! number of grid point in x
  integer                               :: n2 
  integer                               :: n3
  integer                               :: p1      ! parallel procs in stream 
  integer                               :: p3      ! parallel procs in span 
  integer                               :: p_solver_type
  integer                               :: precond_type
  integer                               :: p_order
  integer                               :: igrid
  integer                               :: halo 
  integer                               :: log_num
  real                                  :: chlx 
  real                                  :: chlz 
  real                                  :: rhx 
  real                                  :: rhz 
  real                                  :: gam 
  real, dimension(:,:), allocatable     :: y 
  real, dimension(  :), allocatable     :: hy  
  real, dimension(:,:), allocatable     :: jacb
  logical                               :: lgrid 
  logical                               :: firstloop
  logical                               :: flat_EE 
  logical                               :: validation 
  integer                               :: validation_flag
  logical                               :: ab_flag
  logical                               :: log_flag 
  logical                               :: debugger
  integer                               :: int_time 
  integer                               :: postprocess_time
  integer                               :: postprocess_time_below_p
  integer                               :: postprocess_time_interval  
  integer                               :: postprocess_time_i
  integer                               :: postprocess_time_below_p_i
  real                                  :: wavelength_p, sf, gf
  integer                               :: io_ver_number

  !------------
  ! modified BC variables
  !------------
  real                                  :: init_time
  ! wall velocity
  logical                               :: wall_bc
  real                                  :: wall_mean
  real                                  :: wall_del
  real                                  :: wall_per
  ! slip velocity
  logical                               :: slip_bc
  real                                  :: slip_mean
  real                                  :: slip_del
  real                                  :: slip_per
  real                                  :: b_d_x, b_d_z
  ! patterned no shear
  logical                               :: pat_bc
  integer                               :: wpst_x
  integer                               :: dpst_x
  integer                               :: spst_x
  integer                               :: wpst_z
  integer                               :: dpst_z
  integer                               :: spst_z
  real                                  :: sb_old
  real                                  :: sb_new
  real                                  :: st_old
  real                                  :: st_new
  real                                  :: sw_gb
  real                                  :: sw_gt
  


  !-----------
  ! logical options , see sub.read_input for 
  ! explanation 
  !-----------

  logical                               :: AB3flag = .true. 
  logical                               :: ostat
  logical                               :: ocorr
  logical                               :: omass      
  logical                               :: rstrt
  logical                               :: rstrt_intrp
  logical                               :: rstrt_AB  = .false.
  logical                               :: ofiel 
  logical                               :: jsgs 
  logical                               :: jbard 
  logical                               :: jdrm 
  logical                               :: jpgctl 
  logical                               :: jxzf 
  logical                               :: jyf 
  logical                               :: testf3d
  logical                               :: use_hb 
  logical                               :: attach_filtdns 
  logical                               :: implicit_loop
  logical                               :: dump_eta
  logical                               :: output_vort
  logical                               :: transportee! = .true.

  integer                               :: i_print ! counts number in diasnostics
  integer                               :: BCswitch ! = 0
  integer                               :: ICswitch ! = 1021
  logical                               :: IMFMswitch != .true. 
  logical                               :: match_bc
  integer                               :: ICswitchBF != 0
  logical                               :: IMFMswitchBF != .false.

  !-----------
  ! other input variables   qq = 30;
  !-----------


  integer                               :: imax 
  integer                               :: irestart
  integer                               :: ist
  integer                               :: isam 
  integer                               :: ipgctl 
  real                                  :: wpgctl 
  real                                  :: re 
  real                                  :: dt 
  real                                  :: sstat 
  real                                  :: U_conv
  integer                               :: iens 
  integer                               :: noutput
  integer                               :: n_interface
  integer                               :: n_div_print
  integer                               :: n_test_module
  integer                               :: nsym_ie 
  integer                               :: nsym_it
  integer                               :: nslice 
  integer                               :: wavenumber
  integer                               :: io_ver_number_prev


  !--------
  ! state variables 
  !--------

  real, dimension(:,:,:,:), allocatable :: q              ! u, transporter velocity vector
  real, dimension(:,:,:,:), allocatable :: q_old
  real, dimension(:,:,:),   allocatable :: qcrit
  real, dimension(  :,:,:), allocatable :: p              ! phi state vector
  real, dimension(  :,:,:), allocatable :: p_old 
  real, dimension(:,:,:,:), allocatable :: q_tee          ! v, transportee velocity vector
  !real, dimension(:,:,:,:), allocatable :: q_tee_0        ! v, transportee velocity vector initial value
  real, dimension(:,:,:,:), allocatable :: q_source       ! v, transportee velocity vector
  real, dimension(:,:,:,:), allocatable :: q_old_tee      ! v, transportee velocity vector
  real, dimension(:,:,:),   allocatable :: q_crit_tee
  real, dimension(  :,:,:), allocatable :: p_tee          ! pressure for transportee 
  real, dimension(  :,:,:), allocatable :: p_old_tee 
  real, dimension(  :,:,:), allocatable :: curvature

  real, dimension(:,:,:  ), allocatable :: vort_wall
  real, dimension(:,:,:,:), allocatable :: vort_wall_plane
  integer, parameter                    :: i_u = 1        ! u vel indexing 
  integer, parameter                    :: i_v = 2        ! v vel indexing
  integer, parameter                    :: i_w = 3        ! w vel indexing
  integer, parameter                    :: i_eta = 4
  integer, parameter                    :: i_top = 2      ! top surface indexing
  integer, parameter                    :: i_bottom = 1    ! bottome surface indexing
  integer, parameter                    :: nvel =3        
  real                                  :: pi 
  real                                  :: time 
  real, dimension(  :,:,:), allocatable :: cf             ! wn boundary cond coeffs
  real, dimension(  :,:,:), allocatable :: cf_slip        ! wn slip bc coeffs
  real                                  :: pgm            ! pressure gradient ctl 
  real                                  :: pgm_tee        ! pressure gradient ctl for transportee
  real                                  :: nu_molec       ! 1/re 
  integer                               :: relax_intv
  real                                  :: relax_factor
  real                                  :: finish_time 
  
  real                                  :: v1_bottom = 0.d0
  real                                  :: v1_top = 0.d0
  !real                                  :: v1_bottom = -1.d0
  !real                                  :: v1_top = 1.d0

  !-------
  ! numerics 
  !-------
  integer, parameter                    :: kl = 2        
  integer, parameter                    :: ku = 2 
  integer, parameter                    :: klp= 3 
  integer, parameter                    :: kup= 3
  integer, parameter                    :: klpp= 1 
  integer, parameter                    :: kupp= 1 
  integer, parameter                    :: molc_visc= 2345
  integer, parameter                    :: eddy_visc= 3456 
  integer                               :: pressure_counts = 0 
  real, dimension(    :,:), allocatable :: rkc            ! rk step coefficients 
  real, dimension(:,:,:,:), allocatable :: f_old          ! explicit nonlinear terms 
  real, dimension(:,:,:,:), allocatable :: f_old2          ! explicit nonlinear terms 
  real, dimension(:,:,:,:), allocatable :: f_old_tee          ! explicit nonlinear terms for transportee
  real, dimension(:,:,:,:), allocatable :: f_old_source          ! explicit nonlinear terms for transportee
  real, dimension(:,:,:,:), allocatable :: f_old2_tee          ! explicit nonlinear terms for transportee 
  real, dimension(    :,:), allocatable :: lub 
  real, dimension(    :,:), allocatable :: lub_pat
  real, dimension(    :,:), allocatable :: lub_tee
  real, dimension(    :,:), allocatable :: lwb
  real, dimension(    :,:), allocatable :: lvb
  real, dimension(    :,:), allocatable :: evb_u, evb_v, evb_w ! implicit eddy visc operators 
  real, dimension(    :,:), allocatable :: pwb
  real, dimension(    :,:), allocatable :: pwb_2nd
  real, dimension(    :,:), allocatable :: pwb_2nd_Dirichlet
  real, dimension(    :,:), allocatable :: pwb_2nd_Dirichlet_Neumann

  real, dimension(      :), allocatable :: kxsq 
  real, dimension(      :), allocatable :: kzsq 
  real, dimension(      :), allocatable :: kxsq_2nd
  real, dimension(      :), allocatable :: kzsq_2nd
  real, dimension(  :,:,:), allocatable :: sc             ! poisson solve , imag component
  real, dimension(  :,:,:), allocatable :: sc_test            ! poisson solve , imag component
  complex, dimension(   :), allocatable :: data           ! fft 
  complex, dimension(   :), allocatable :: data_2D        ! fft 

  ! patterning (Kim: 09.14.22)
  real                                  :: nsolid_u = 0.d0
  real                                  :: nsolid_w = 0.d0
  integer, dimension(:,:,:), allocatable:: phase

  ! fft edits, 6/18/2010
  integer, dimension(:,:,:), allocatable:: in_new 
  integer, dimension(:,:,:), allocatable:: out_new 
  integer, dimension(:,:,:), allocatable:: finish_new 

  integer, dimension(:,:,:), allocatable:: in_2D
  integer, dimension(:,:,:), allocatable:: out_2D
  integer, dimension(:,:,:), allocatable:: finish_2D

  !-------
  ! IMFM variables 
  !-------
  real, dimension(:,:), allocatable :: s_tee              ! s, source term for transportee variable function of y
  
  !-------
  ! ensemble statistics variables 
  !-------

  real, dimension(:,:)    , allocatable :: um      ! velocity mean 
  real, dimension(:,:)    , allocatable :: uum     ! u^2 
  real, dimension(:)      , allocatable :: vtm     ! eddy viscosity mean 
  real, dimension(:,:)    , allocatable :: ensm    ! other ensemble means 
  real, dimension(:)      , allocatable :: uvm 
  real, dimension(:)      , allocatable :: uvm_sgs ! subgrid contribution to tau_12
  real, dimension(:)      , allocatable :: uvm_sgs_sq ! for subgrid rms contrib 
  real, dimension(:)      , allocatable :: const_m ! added to keep tabs on the constant.. 5/19/2010
  real, dimension(:)      , allocatable :: const_sq 
  real, dimension(:,:)    , allocatable :: kstar
  real, dimension(:)      , allocatable :: ppm
  real, dimension(:)      , allocatable :: pm
  real, dimension(:)      , allocatable :: p_rms
  real, dimension(:)      , allocatable :: p_rms_inst
  real, dimension(:)      , allocatable :: pmpm
  real, dimension(:)      , allocatable :: sqs_pp
  real, dimension(:)      , allocatable :: sum_pp

  !-------
  ! sgs stats added 5/5/2010
  !-------
  real, dimension(:,:)   , allocatable :: tau_mean
  real, dimension(:,:)   , allocatable :: tau_sq


  !------
  ! sgs variables 
  !------

  real, dimension(:)      , allocatable :: const       ! smagorinsky constant 
  real, dimension(:,:,:)  , allocatable :: vt          ! eddy viscosity 
!  real, dimension(:,:,:) , allocatable  :: smag_sq
  real                                  :: asq 
  integer, dimension(2)                 :: nclip
  integer                               :: nconst      ! output of smag const
  real                                  :: fsgs        ! toggle implicit eddy visc op
  logical                               :: lfilt_visc 
  logical                               :: leddyv_grad 
  logical                               :: leddyv_cross 
  logical                               :: limp_pred
  logical                               :: leddyv_exp=.true.  ! treat wn eddyv terms explicitly
  logical                               :: bypass_sgs  ! use explicit filtering w/o subgrid model 
  real                                  :: vt_bar      ! mean vt over the entire field 


  ! these variables are for ksgs, comment out if not using ksgs
!  real, dimension(:,:,:), allocatable   :: kpsgs       ! \bar{k'} ... filter ksgs 


  !-------
  ! code verification variables 
  !-------
  logical                               :: osverify 
  real                                  :: alpha 
  real                                  :: eps 
  logical                               :: mmsverify 
  logical                               :: output_qcrit
 

  !-------
  ! MPI variables, cartesian topology, etc
  !-------
  real                                  :: gdiv, globaldiv 
  integer                               :: is,ie, ks,ke , js, je
  integer                               :: halo_x,halo_z, halo_yu, halo_yp, fhalo, fh_h 
  integer                               :: mycomm, comm3d 
  integer                               :: left, right, top, botm, front, back , ierr
  integer, dimension(0:2)               :: coords 
  integer                               :: uvec_comm, pvec_comm , fft_comm, solv_comm, eta_comm
  integer                               :: nx, ny, nz , shft2, nxw, nzw
  integer                               :: rank, size 
  !integer,parameter                     :: MPI_DOUBLE_PRECISION=MPI_DOUBLE_PRECISION
  integer                               :: span_block
  integer                               :: is_df, ie_df, ks_df, ke_df            ! bounds for fields needing 2 filtering ops
  integer                               :: is_sf, ie_sf, ks_sf, ke_sf            ! bounds for fields needing 1 filtering op 
  integer                               :: dhalo_visc 

  !-------
  ! time stepping variables 
  !-------
  integer                               :: itimestep, nrksteps, ABstep
  real                                  :: tcfl 
  logical                               :: usecfl 

contains 


  !===================================================================================!

  subroutine read_input

   implicit none 
   logical         :: input_new 
   integer         :: input_version
    !-----------
    ! reads the given input data file (input/input.dat) and 
    ! assigns variables 
    !-----------
    open(unit=1, file='input/input_io.dat', status='old', form='formatted') 
    read(1, '(i8)') input_version             ! processors in stream 
    if (rank .eq. 0) write(*,*) 'input_version =', input_version
    !---------
    ! processor variables 
    !---------
    read(1, '(/,i8)'   ) p1             ! processors in stream 
    read(1, '(i8)'     ) p3             ! processors in span 
    !--------- 
    ! grid size variables 
    !---------
    read(1,'(/,i8)' )  n1
    read(1,'(i8)')     n2
    read(1,'(i8)')     n3
    read(1,'(d12.5)')  chlx 
    read(1,'(d12.5)')  chlz
    read(1,'(d12.5)') re            ! reynolds number 
    !------
    ! Restart condition 
    !------
    read(1,'(/,l1)')   rstrt          ! use restart 
    read(1,'(l1)')     rstrt_intrp    ! use rstrt from collapsed grid and intrp 
    read(1,'(d12.5)')  init_time      ! initial time from restart file
    !------
    ! i/o parameters 
    !------
    read(1,'(/,l1)')   log_flag
    read(1,'(l1)')     ostat          ! output statistics 
    read(1,'(l1)')     ofiel          ! output field 
    read(1,'(l1)')     dump_eta
    read(1,'(l1)')     output_qcrit
    read(1,'(l1)')     output_vort
    !------ 
    ! Time stepping 
    !------
    read(1,'(/,i3)')  itimestep       ! set time stepping parameters 
    read(1,'(l1)')    usecfl          ! toggle cfl 
    read(1,'(d12.5)')   dt            ! time step 
    read(1,'(d12.5)') tcfl            ! target cfl 
    !------ 
    ! Job control/sampling/output interval 
    !------
    read(1,'(/,i8)')    imax            ! total number of time steps 
    read(1,'(i8)')    irestart        ! dumping restart files
    read(1,'(i8)')    ist             ! time steps btwn stat sampling
    read(1,'(i8)')    isam            ! time steps btwn interface output/vorticity
    read(1,'(d12.5)') sstat           ! time to start collecting statistics 
    read(1,'(d12.5)') finish_time     ! time to finish simulation
    read(1,'(i8)')     postprocess_time
    read(1,'(i8)')     postprocess_time_interval

    !------
    ! Pressure solver conditions
    !------
    read(1, '(/,i8)' ) p_solver_type  ! hypre_solver 1:GMRES  2:AMG  3: HYBRYD  4: BiCG
    read(1, '(i8)'   ) precond_type
    read(1, '(i8)'   ) p_order        ! order of accuracy of Poisson solver 

    !------ 
    ! wall normal grid 
    !------
    read(1,'(/,i3)')   igrid          ! grid type 
    read(1,'(d12.5)')  gam            ! stretching for hyp. tangent 
    read(1,'(l1)'   )  lgrid          ! load grid from file 

    !------ 
    ! IMFM 
    !------
    read(1,'(/,l1)')   transportee    ! transportee variables
    read(1,'(l1)')     IMFMswitch     ! inverse macroscopic forcing method
    read(1,'(i8)')     ICswitch       ! initial conditions and forcing direction
    read(1,'(i8)')     BCswitch       ! boundary conditions
    read(1,'(l1)')     match_bc       ! match transporter and transportee BCs
    read(1,'(l1)')     IMFMswitchBF   ! IMFM brute force switch
    read(1,'(i8)')     ICswitchBF     ! initial conditions for brute force MFM

    !------
    ! Wall velocity BC
    !------
    read(1,'(/,l1)')   wall_bc        ! wall velocity on/off
    read(1,'(d12.5)')  wall_mean      ! wall center velocity
    read(1,'(d12.5)')  wall_del       ! wall velocity delta
    read(1,'(d12.5)')  wall_per       ! wall period (in plus units)

    !------
    ! Slip velocity BC
    !------
    read(1,'(/,l1)')   slip_bc        ! slip velocity on/off
    read(1,'(d12.5)')  slip_mean      ! slip length center
    read(1,'(d12.5)')  slip_del       ! slip length delta
    read(1,'(d12.5)')  slip_per       ! slip period (in plus units) 

    !------
    ! Patterned no shear BC
    !------
    read(1,'(/,l1)')   pat_bc         ! patterned no shear on/off
    read(1,'(i8)')     wpst_x         ! x bubble width (in indices)
    read(1,'(i8)')     dpst_x         ! x bubble distance (in indices)
    read(1,'(i8)')     spst_x         ! x bubble start (in indices)
    read(1,'(i8)')     wpst_z         ! z bubble width (in indices)
    read(1,'(i8)')     dpst_z         ! z bubble distance (in indices)
    read(1,'(i8)')     spst_z         ! z bubble start (in indices) 

    close(1) 
  
    
    osverify  = .false. 
    alpha = 1.0d0
    eps = 1.0d-06
    mmsverify = .false. 

    !------------- 
    ! check time stepping input 
    !-------------
  
    if (itimestep .ne. 4 .and. itimestep .ne. 3 .and. itimestep .ne. 2 .and. itimestep .ne. 1) then 
       write(*,*) 'not a valid time stepping choice ....', itimestep 
       itimestep =1                       
    endif

    !---- 
    ! set number of substeps in time stepping 
    !----
    if (itimestep .eq. 1 .or. itimestep .eq. 3 .or. itimestep .eq. 4) nrksteps =1 
    if (itimestep .eq. 2) nrksteps =3

    if (itimestep .eq. 1) io_ver_number = 5003
    if (itimestep .eq. 3) io_ver_number = 5004
    if (itimestep .eq. 4) io_ver_number = 5005

    
    if ( jsgs .and. leddyv_grad) then 
       fsgs = 1.d0 
    else  
       fsgs = 0.d0 
    endif 

    call check_deprecated_options

    !------
    ! assign slip BC variables
    !------
    if (.not. slip_bc) then
        b_d_x = 0.d0
        b_d_z = 0.d0
    else
        b_d_x = slip_mean/re
        b_d_z = slip_mean/re
    endif

    !------ 
    ! compute hx, hz 
    !------
    pi = acos(-1.d0) 
    rhx = dble(n1)  / (chlx*pi) 
    rhz = dble(n3)  / (chlz*pi)
    !----- 
    ! compute 1/re 
    !-----
    nu_molec = 1.d0/re 


  if (rank .eq. 0 .and. log_flag ) then 
  write(*,*) '0-1. Input statement'
  write(*,*) '0-1-1. n1=',n1,', n2=',n2, 'n3=',n3 
  write(*,*) '0-1-1. Lx+=',chlx*re*pi,',  Lz+=',chlz*re*pi
  write(*,*) '0-1-1. Re_tau=',re
  write(*,*) '0-1-1. Time step, dt = ',dt, ' dt+ = ',dt*re
  write(*,*) '0-1-1. Grid resolution, dx+ = ',1.d0/rhx*re
  write(*,*) '0-1-1. Grid resolution, dz+ = ',1.d0/rhz*re
  write(*,*) '0-1-1. Grid min(dy+) = ', re*(-tanh(gam* (1.d0 - 2.0d0 / dble(n2)*dble(2))) / tanh(gam) + tanh(gam* (1.d0 - 2.0d0 / dble(n2)*dble(1))) / tanh(gam))
  write(*,*) '0-1-1. Current input/output version number:', io_ver_number
  endif

  end subroutine read_input 

  !===================================================================================!
  
  subroutine setup_pat

    implicit none
    integer                               :: igstrt, kgstrt, igend, kgend, i, k
    integer                               :: nphase_u, jBF
    real                                  :: deltaBF
    integer, dimension(1:3,0:n1-1,0:n3-1) :: phase_t

    phase = 0
    nsolid_u = (wpst_x+1) * (wpst_z  ) * n1 / (wpst_x+dpst_x) * n3 / (wpst_z+dpst_z)
    nsolid_w = (wpst_x  ) * (wpst_z+1) * n1 / (wpst_x+dpst_x) * n3 / (wpst_z+dpst_z)
    
    nphase_u = n1*n3-nsolid_u
    !sb_old = ( -1.d0*n1*n3 - nphase_u*y(i_u, 1) - 0 ) / nsolid_u
    !st_old = (  1.d0*n1*n3 - nphase_u*y(i_u,n2) - 0 ) / nsolid_u 
    jBF = ICswitchBF
    deltaBF = y(i_u,jBF+1)-y(i_u,jBF)
    sb_old = 0.d0
    st_old = deltaBF
    sb_new = sb_old
    st_new = st_old
    sw_gb = 0.d0 * nphase_u
    sw_gt = deltaBF * nphase_u

 
    do i=is,ie
      do k=ks,ke
        phase_t(i_u,i,k) = 0
        phase_t(i_v,i,k) = 0
        phase_t(i_w,i,k) = 0
      enddo
    enddo

    igstrt = spst_x
    kgstrt = spst_z
    igend  = igstrt + wpst_x
    kgend  = kgstrt + wpst_z

    do k=0,n3-1
     
      ! reset igstrt and igend for new k
      igstrt = spst_x
      igend  = igstrt + wpst_x

      do i=0,n1-1

        if ( i .ge. igstrt .and. i .lt. igend+1 .and. k .ge. kgstrt .and. k .lt. kgend ) then
          phase_t(i_u,i,k) = 0
        else
          phase_t(i_u,i,k) = 1
        endif

        if ( i .ge. igstrt .and. i .lt. igend .and. k .ge. kgstrt .and. k .lt. kgend+1 ) then
          phase_t(i_w,i,k) = 0
        else
          phase_t(i_w,i,k) = 1
        endif

        if ( i .ge. igstrt .and. i .lt. igend .and. k .ge. kgstrt .and. k .lt. kgend ) then
          phase_t(i_v,i,k) = 0
        else
          phase_t(i_v,i,k) = 1
        endif

        ! loop to next groove in x
        if ( i .eq. igend ) then
          igstrt = igend + dpst_x
          igend = igstrt + wpst_x
        endif

      enddo

      ! loop to next groove in z
      if ( k .eq. kgend ) then
        kgstrt = kgend + dpst_z
        kgend = kgstrt + wpst_z
      endif

    enddo

    ! write to global variable
    do i=is,ie
      do k=ks,ke
        phase(i_u,i,k) = phase_t(i_u,i,k)
        phase(i_v,i,k) = phase_t(i_v,i,k)
        phase(i_w,i,k) = phase_t(i_w,i,k)
      enddo
    enddo
        
    call mpi_barrier( mycomm, ierr)
    
  end subroutine setup_pat

  !===================================================================================!

  subroutine domain_decomp 


    implicit none 
    integer                   :: ndim , reorder , all_y, all_z, all_yp
    integer                   :: modx,  modz 
    integer, dimension(0:2)   :: dims
    logical, dimension(0:2)   :: perdc 

    !------
    ! check parallel initialization 
    !------

    if ( rank .eq. 0 .and. size .ne. p1*p3 ) then 
       write(*,*) 'global:: bad initialization ... ' 
       write(*,*) 'global:: p1 = ', p1, ', p3 = ', p3, ', size = ', size 
       write(*,*) 'global:: exiting....' 
       
       call mpi_finalize ( ierr ) 
       stop 
    endif 

    !-----
    ! create cartesian topology for channel 
    !-----
    ndim    = 3
    dims(0) = p1 
    dims(1) = 1       ! no parallelization in wall normal direction 
    dims(2) = p3 


    perdc(0)= .true.       ! periodic in stream/span  
    perdc(1)= .false. 
    perdc(2)= .true.       


    reorder = 1 

    call mpi_cart_create( mycomm, ndim, dims, perdc, reorder, comm3d, ierr ) 
    call mpi_cart_get   ( comm3d, ndim, dims, perdc, coords, ierr   ) 

    call mpi_cart_shift ( comm3d, 0, 1, left, right,ierr      ) 
    call mpi_cart_shift ( comm3d, 1, 1, botm, top  ,ierr      ) 
    call mpi_cart_shift ( comm3d, 2, 1, back, front,ierr      ) 
    


    !-------
    ! determine local dimensions of box 
    ! currently, assume that the box is evenly divisible 
    ! TODO : Non-even domain decomp 
    !-------

    nx = int ( n1 / p1 ) 
    nz = int ( n3 / p3 ) 

    modx = mod(n1, p1) 
    modz = mod(n3, p3) 

    is =  coords(0) * nx
    ie =  is + nx -1

    ks =  coords(2) * nz
    ke =  ks + nz -1


    ny = int ( n2 / size)
    js  = rank* ny + 1
    je  = js +  ny -1


    if ( coords(0 ) .eq. p1-1 ) then 
       shft2 = 2 
    else  
       shft2 = 0 
    endif 

    
    !-------
    ! define and add halo to dimensions here 
    ! TODO : Dynamically determine halo size , depending on size 
    !        of diff stencil and filtering stencil 
    !-------

    halo_x = 3 
    halo_yu= 3 
    halo_yp= 2 
    halo_z = 3 
    dhalo_visc = max( kl, ku) 

    if ( (jsgs .and. jxzf .and. (.not. limp_pred)) ) then 
       fhalo = nsym_it 
       fh_h  = nsym_ie 

       !--- 
       ! no double filtered quantities 
       !---
       ks_df = ks - fh_h
       ke_df = ke + fh_h 
       is_df = is - fh_h 
       ie_df = ie + fh_h 

       ks_sf = ks 
       ke_sf = ke 
       is_sf = is 
       ie_sf = ie 

    elseif ( bypass_sgs    ) then 
       fhalo = nsym_ie 
       fh_h  = nsym_ie 

       ks_df = ks - fh_h 
       ke_df = ke + fh_h 
       is_df = is - fh_h 
       ie_df = ie + fh_h 

       ks_sf = ks 
       ke_sf = ke 
       is_sf = is 
       ie_sf = ie 

    elseif( attach_filtdns ) then 
       fhalo = nsym_ie 
       fh_h  = nsym_ie 

       !---
       ! no halo regions for advanced 
       ! quantities 
       !--- 
       ks_df = ks 
       ke_df = ke 
       is_df = is
       ie_df = ie 

       ks_sf = ks 
       ke_sf = ke 
       is_sf = is 
       ie_sf = ie 



    elseif  ( jsgs .and. jxzf .and. limp_pred)  then 
       
       if ( nsym_it .ge. 2*nsym_ie+dhalo_visc ) then 
          fhalo = nsym_it 
          fh_h  = nsym_ie 
       else  
          fhalo = 2*nsym_ie + dhalo_visc 
          fh_h  = nsym_ie 
       endif

       !--- 
       ! define halo for double filtering 
       ! operations 
       !--- 
       ks_df  = ks - 2*fh_h - dhalo_visc 
       ke_df  = ke + 2*fh_h + dhalo_visc 
       is_df  = is - 2*fh_h - dhalo_visc 
       ie_df  = ie + 2*fh_h + dhalo_visc
       
       ks_sf  = ks -   fh_h - dhalo_visc 
       ke_sf  = ke +   fh_h + dhalo_visc
       is_sf  = is -   fh_h - dhalo_visc
       ie_sf  = ie +   fh_h + dhalo_visc

    elseif ( jsgs .and. (.not. jxzf) ) then 
       
       fhalo = nsym_it 
       fh_h  = 0 

       !--- 
       ! no explicit filtering 
       !---

       ks_df  = ks 
       ke_df  = ke 
       is_df  = is 
       ie_df  = ie 

       ks_sf  = ks 
       ke_sf  = ke 
       is_sf  = is 
       ie_sf  = ie 
       
    else  

       !--- 
       ! no filtering halo regions 
       ! needed for no sgs and no exp filt
       !--- 
       fhalo = 0 
       fh_h  = 0 

       ks_df  = ks 
       ke_df  = ke 
       is_df  = is 
       ie_df  = ie 
       
       ks_sf  = ks 
       ke_sf  = ke 
       is_sf  = is 
       ie_sf  = ie 
       
    endif


    !----- 
    ! add differentiation halo to the 
    ! filter halo 
    !-----

    if ( jsgs .or. attach_filtdns .or. bypass_sgs) then 
       halo_z = halo_z + fhalo 
       halo_x = halo_x + fhalo 
    endif 

    

    !------
    ! setup derived data types for velocity, pressure 
    ! communication 
    !------

    all_y = n2+2*halo_yu 
    all_z = nz+2*halo_z 

    all_yp= n2+2*halo_yp

    call mpi_type_vector ( all_y *all_z, nvel*halo_x, nvel*(nx+2*halo_x  ), MPI_DOUBLE_PRECISION, uvec_comm, ierr) 
    call mpi_type_vector ( all_yp*all_z,      halo_x,      (nx+2*halo_x  ), MPI_DOUBLE_PRECISION, pvec_comm, ierr) 
    call mpi_type_vector ( (nz+2),                 1,      (nx+2  ),        MPI_DOUBLE_PRECISION, eta_comm, ierr)


    call mpi_type_commit ( uvec_comm , ierr) 
    call mpi_type_commit ( pvec_comm , ierr) 
    call mpi_type_commit ( eta_comm , ierr)


  end subroutine domain_decomp

!===========================================================================================================================! 

  subroutine mem_alloc 

    !--------
    ! allocate memory for each process; implicit operators are defined in the same 
    ! way as in the serial code 
    !--------

    implicit none 
    integer          :: fsize 
    

    fsize = nx*nz

    pi = acos(-1.0d0) 

    !-----
    ! allocate memory for grid variables 
    !-----

    allocate ( y(2, 0:n2+1) ) 
    allocate ( hy(-2:n2+3 ) ) 
    allocate ( jacb(2,0:n2) ) 

    !-----
    ! allocate numerics 
    !-----

    allocate ( rkc(1:3,1:4) ) 
    allocate ( f_old(1:nvel,is_sf:ie_sf,1:n2,ks_sf:ke_sf) ) 
    allocate ( f_old2(1:nvel,is_sf:ie_sf,1:n2,ks_sf:ke_sf) )
    allocate ( f_old_tee(1:nvel,is_sf:ie_sf,1:n2,ks_sf:ke_sf) ) 
    allocate ( f_old_source(1:nvel,is_sf:ie_sf,1:n2,ks_sf:ke_sf) ) 
    allocate ( f_old2_tee(1:nvel,is_sf:ie_sf,1:n2,ks_sf:ke_sf) )
    allocate ( lub(2*kl+ku+1,1:n2) ) 
    allocate ( lub_pat(2*kl+ku+1,1:n2) ) 
    allocate ( lub_tee(2*kl+ku+1,1:n2) )
    allocate ( lwb(2*kl+ku+1,1:n2) )
    allocate ( lvb(2*kl+ku+1,1:n2-1))
    allocate ( evb_u(2*kl+ku+1,1:n2))
    allocate ( evb_v(2*kl+ku+1,1:n2-1))
    allocate ( evb_w(2*kl+ku+1,1:n2))
    allocate ( pwb(2*klp+kup+1,1:n2) )
    allocate ( pwb_2nd(2*klpp+kupp+1,1:n2) )
    allocate ( pwb_2nd_Dirichlet(2*klp+kup+1,1:n2) )
    allocate ( pwb_2nd_Dirichlet_Neumann(2*klp+kup+1,1:n2) )

    allocate ( kxsq(0:n1-1)         ) 
    allocate ( kzsq(0:n3-1)         ) 
    allocate ( kxsq_2nd(0:n1-1)         )
    allocate ( kzsq_2nd(0:n3-1)         )
    allocate ( sc(is+1:ie+1,1:n2,ks+1:ke+1)  ) 
    allocate ( sc_test(is+1:ie+1,1:n2,ks+1:ke+1)  )

    allocate ( phase(3,is:ie,ks:ke)         )
    
    ! fft edits 6/18/2010
    allocate ( data(nx*nz*n2)) 
    allocate ( data_2D(nx*nz*2))
    allocate (in_new(is+1:ie+1,1:n2,ks+1:ke+1))
    allocate (out_new(is+1:ie+1,1:n2,ks+1:ke+1)) 
    allocate (finish_new(is+1:ie+1,1:n2,ks+1:ke+1)) 

    allocate (in_2D(is+1:ie+1,1:2,ks+1:ke+1))
    allocate (out_2D(is+1:ie+1,1:2,ks+1:ke+1))
    allocate (finish_2D(is+1:ie+1,1:2,ks+1:ke+1))




    !------
    ! allocates memory to state variables 
    !------


    allocate ( q(3, is-halo_x:ie+halo_x,1-halo_yu:n2+halo_yu, ks-halo_z:ke+halo_z) )
    allocate ( q_old(3, is-halo_x:ie+halo_x,1-halo_yu:n2+halo_yu, ks-halo_z:ke+halo_z) )
    allocate ( qcrit(is+1:ie+1,1:n2,ks+1:ke+1))
    allocate ( curvature(is-halo_x:ie+halo_x, 2, ks-halo_z:ke+halo_z) )

    allocate ( p(   is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp, ks-halo_z+1:ke+halo_z+1) ) 
    allocate ( p_old(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp, ks-halo_z+1:ke+halo_z+1) ) 
    allocate ( cf      ( 1:3   ,-2:0   ,  1:2   ) )
    allocate ( cf_slip ( 1:3   ,-2:0   ,  1:2   ) )
    allocate ( vort_wall       (3,is:ie,ks:ke) )
    allocate ( vort_wall_plane (3,0:3,is:ie,ks:ke) )
    
    !------
    ! allocates memory to transportee variables 
    !------

    allocate ( q_tee(3, is-halo_x:ie+halo_x,1-halo_yu:n2+halo_yu, ks-halo_z:ke+halo_z) )
    !allocate ( q_tee_0(3, is-halo_x:ie+halo_x,1-halo_yu:n2+halo_yu, ks-halo_z:ke+halo_z) )
    allocate ( q_source(3, is-halo_x:ie+halo_x,1-halo_yu:n2+halo_yu, ks-halo_z:ke+halo_z) )
    allocate ( q_old_tee(3, is-halo_x:ie+halo_x,1-halo_yu:n2+halo_yu, ks-halo_z:ke+halo_z) )
    allocate ( q_crit_tee(is+1:ie+1,1:n2,ks+1:ke+1))
    allocate ( p_tee(   is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp, ks-halo_z+1:ke+halo_z+1) ) 
    allocate ( p_old_tee(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp, ks-halo_z+1:ke+halo_z+1) ) 

    !-------
    ! allocate memory for IMFM
    !-------
    allocate ( s_tee(3, 1:n2  ) )

    !-------
    ! allocate memory for ensemble mean statistics 
    !-------
    allocate ( um(3, 1:n2  ) )
    allocate ( uum(3,1:n2  ) )
    allocate ( vtm(1:n2  )   )
    allocate ( uvm(1:n2)     )
    allocate ( uvm_sgs(1:n2) )
    allocate ( uvm_sgs_sq(1:n2))
    allocate ( const_m(1:n2))
    allocate ( const_sq(1:n2))
    allocate ( kstar(1:n2,2))
    allocate ( pm (1:n2) )
    allocate ( ppm(1:n2) )
    allocate ( p_rms(1:n2) )
    allocate ( p_rms_inst(1:n2) )
    allocate ( pmpm(1:n2) )

    um = 0.d0
    uum= 0.d0
    vtm= 0.d0
    uvm= 0.d0
    uvm_sgs = 0.d0
    uvm_sgs_sq = 0.d0
    const_m = 0.d0
    const_sq = 0.d0
    kstar = 0.d0
    pm = 0.d0
    ppm = 0.d0
    p_rms = 0.d0
    p_rms_inst = 0.d0
    pmpm = 0.d0
    sqs_pp = 0.d0
    sum_pp = 0.d0

    !------
    ! if allocating other ens mean variables , please add here and comment 
    ! what ensm contains 
    !------

    allocate ( tau_mean(1:6, 1:n2))
    allocate ( tau_sq(1:6, 1:n2))

    tau_mean = 0.d0 
    tau_sq   = 0.d0 


    !------
    ! sgs variables 
    ! force allocation of vt for i/o purposes
    !------
    allocate ( vt   (is-halo_x:ie+halo_x, 1-halo_yp:n2+halo_yp, ks-halo_z:ke+halo_z) ) 

    if ( jsgs ) then 
       
       allocate ( const(1:n2) ) 
       vt_bar = 0.d0 
    endif

  end subroutine mem_alloc

!=======================================================================================================================!


  
  subroutine mem_dealloc

    implicit none

    !-----
    ! deallocates memory to state variables 
    !-----

    if (allocated(y))     deallocate ( y ) 
    if (allocated(hy))    deallocate ( hy) 
    if (allocated(jacb))  deallocate ( jacb) 


    if (allocated(q))     deallocate ( q) 
    if (allocated(q_old)) deallocate ( q_old)
    if (allocated(qcrit)) deallocate (qcrit)
    if (allocated(p))     deallocate ( p) 
    if (allocated(p_old)) deallocate ( p_old) 
    if (allocated(q_tee))     deallocate ( q_tee) 
    !if (allocated(q_tee_0))     deallocate ( q_tee_0) 
    if (allocated(q_source))     deallocate ( q_source) 
    if (allocated(q_crit_tee)) deallocate (q_crit_tee)
    if (allocated(q_old_tee)) deallocate ( q_old_tee)
    if (allocated(p_tee))     deallocate ( p_tee) 
    if (allocated(p_old_tee)) deallocate ( p_old_tee) 
    if (allocated(cf))    deallocate ( cf)  
    if (allocated(cf_slip))    deallocate (cf_slip)
    if (allocated(curvature))  deallocate (curvature)

    if (allocated(s_tee))     deallocate ( s_tee) 

    if (allocated(rkc))   deallocate (  rkc) 
    if (allocated(f_old)) deallocate (f_old) 
    if (allocated(f_old2)) deallocate (f_old2)
    if (allocated(f_old_tee)) deallocate (f_old_tee) 
    if (allocated(f_old_source)) deallocate (f_old_source) 
    if (allocated(f_old2_tee)) deallocate (f_old2_tee)
    if (allocated(lub))   deallocate (  lub) 
    if (allocated(lub_pat)) deallocate ( lub_pat)
    if (allocated(lub_tee)) deallocate ( lub_tee)
    if (allocated(lwb)) deallocate ( lwb)
    if (allocated(lvb))   deallocate (  lvb) 
    if (allocated(evb_u)) deallocate (evb_u) 
    if (allocated(evb_v)) deallocate (evb_v) 
    if (allocated(evb_w)) deallocate (evb_w)
    if (allocated(pwb))   deallocate (  pwb)
    if (allocated(pwb_2nd))   deallocate (  pwb_2nd)
    if (allocated(pwb_2nd_Dirichlet))   deallocate (  pwb_2nd_Dirichlet)
    if (allocated(pwb_2nd_Dirichlet_Neumann))   deallocate (  pwb_2nd_Dirichlet_Neumann)


    if (allocated(phase)) deallocate ( phase)

    if (allocated(kxsq))  deallocate ( kxsq) 
    if (allocated(kzsq))  deallocate ( kzsq) 
    if (allocated(kxsq_2nd))  deallocate ( kxsq_2nd)
    if (allocated(kzsq_2nd))  deallocate ( kzsq_2nd)
    if (allocated(sc))    deallocate (   sc) 
    if (allocated(sc_test))   deallocate (   sc_test)
    if (allocated(data))  deallocate ( data) 
    if (allocated(data_2D))  deallocate ( data_2D)
    !if (allocated(in))    deallocate (   in) 
    !if (allocated(out))   deallocate (  out) 

    if (allocated(um))    deallocate ( um) 
    if (allocated(uum))   deallocate ( uum) 
    if (allocated(vtm))   deallocate ( vtm) 
    if (allocated(uvm))   deallocate ( uvm) 
    if (allocated(uvm_sgs)) deallocate ( uvm_sgs) 
    if (allocated(uvm_sgs_sq)) deallocate ( uvm_sgs_sq)

    if ( allocated(tau_mean)) deallocate ( tau_mean)
    if ( allocated(tau_sq))   deallocate ( tau_sq)

    if ( allocated(ensm) )   deallocate( ensm ) 
    if ( allocated(vt)   )   deallocate( vt   ) 
    if ( allocated(const))   deallocate( const) 

    deallocate ( in_new) 
    deallocate (out_new)
    deallocate (finish_new)

    deallocate ( in_2D)
    deallocate (out_2D)
    deallocate (finish_2D)

  end subroutine mem_dealloc 


!==================================================================================================================================!

  subroutine check_deprecated_options 

    !----- 
    ! detail the use of the deprecated or non-existent options 
    ! to date and dump the warnings to output/warnings.log 
    !----- 

    implicit none 


    if ( rank .eq. 0 ) then 

       open(unit= 444, file='output/warnings.log', form='formatted', status='unknown') 

       !--- 
       ! semi-implicit eddy viscosity time 
       ! advancement deprecated 
       !---
       if ( limp_pred ) then 
          write(444,*) 'limp_pred::semi implicit time advancement of wn sgs terms currently unstable...'
          write(444,*) 'limp_pred::requires extremely small time step for approx factorization used to be valid...' 
       endif 

       !---- 
       ! cross terms from eddy viscosity does not 
       ! exist currently 
       !---- 
       if ( leddyv_cross ) then 
          write(444,*) 'leddyv_cross:: cross terms for sgs model do not currently exist ...' 
       endif 


       !--- 
       ! no ssm 
       !--- 
       if ( jbard ) then 
          write(444,*) 'jbard:: bardina ssm does not exist ... ' 
       endif 

       !--- 
       ! no drm 
       !--- 
       if ( jdrm ) then 
          write(444,*) 'jdrm:: dynamic reconstruction does not exist...' 
       endif 


       close(444) 
    endif 
       



  end subroutine check_deprecated_options
!===================================================================================================================================!
   subroutine terminate ( message )
   implicit none 
   character(*) :: message
   if (rank .eq. 0) write(*,*) 'Terminate program due to following problem:'
   if (rank .eq. 0) write(*,*) message
   call mem_dealloc
   call mpi_finalize ( ierr )
   end subroutine terminate 


end module global
