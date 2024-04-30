
module numerics 

  use global 
  implicit none 

  include 'fftw_f77.i'
  real                        :: plan, ipln 
  real                        :: plan_2D, ipln_2D

contains 

!===================================================================================================================! 
  subroutine pressure_step(istep) 

    implicit none
    integer, intent(in)                   :: istep
         if (p_order .eq. 2) then
         call p_solve_FFT_2nd( istep )
         elseif (p_order .eq. 4) then
         call p_solve ( istep )
         endif

  end subroutine 
!===================================================================================================================! 
  subroutine setup_rkc

    !----------------------------------
    ! sets the rk stepping coefficients 
    ! indexing changed on \delta values 
    ! coefficients are stored in rkc(i,j), i = step no, (1,2,3) 
    !                                      j = 1 (alpha), 2 (beta), 3 (-gamma), 4 (-delta)
    !----------------------------------
    
    implicit none 
    integer :: ierr 
   
    globaldiv = 1.0d0

    if( itimestep .eq. 1) then 
       !--------------------
       ! use AB2 / CN coefficients 
       !-------------------
       if ( rank .eq. 0 .and. log_flag) write(*,*) '1-1-1. using Adams-Bashforth/Crank-Nicholson 2nd time stepping ...', nrksteps

       rkc(1,1) = 0.5d0 
       rkc(1,2) = 0.5d0 
       rkc(1,3) =-1.5d0 
       rkc(1,4) = 0.5d0 
       
       rkc(2,:) = 0.d0     ! unused
       rkc(3,:) = 0.d0     ! unused 
       
    elseif( itimestep .eq. 3) then 
       !--------------------
       ! use AB2 coefficients 
       !-------------------
       if ( rank .eq. 0 .and. log_flag) write(*,*) '1-1-1. using Adams-Bashforth 2nd order time stepping ...', nrksteps

       rkc(1,1) = 1.5d0 
       rkc(1,2) =-0.5d0 
       rkc(1,3) =-1.5d0 
       rkc(1,4) = 0.5d0 
       
       rkc(2,:) = 0.d0     ! unused
       rkc(3,:) = 0.d0     ! unused 
    elseif( itimestep .eq. 4) then
       !--------------------
       ! use AB3 coefficients 
       !-------------------
       if ( rank .eq. 0 .and. log_flag) write(*,*) '1-1-1. using Adams-Bashforth 3rd order time stepping ...', nrksteps

       rkc(1,1) = (23.d0/12.d0)
       rkc(1,2) =-(16.d0/12.d0)
       rkc(1,3) =-(23.d0/12.d0)
       rkc(1,4) = (16.d0/12.d0)

       rkc(2,:) = 0.d0     ! unused
       rkc(3,:) = 0.d0     ! unused 
    else

       !----------------
       ! use low storage RK3/CN
       !---------------
       if ( rank .eq. 0 .and. log_flag) write(*,*) '1-1-1. using Runge-Kutta 3rd/Crank-Nicholson 2nd time stepping...'

       ! coefficients from le and moin, jcp, 1991
       !-----

       rkc(1,1)  =  4.d0 / 15.d0
       rkc(1,2)  =  4.d0 / 15.d0 
       rkc(1,3)  = -8.d0 / 15.d0 
       rkc(1,4)  =  0.0d0           
       
       rkc(2,1)  =   1.d0 / 15.d0 
       rkc(2,2)  =   1.d0 / 15.d0 
       rkc(2,3)  = -5.d0 / 12.d0 
       rkc(2,4)  =  17.d0/ 60.d0 
       
       rkc(3,1)  =  1.d0 / 6.d0 
       rkc(3,2)  =  1.d0 / 6.d0 
       rkc(3,3)  = -3.d0 / 4.d0 
       rkc(3,4)  =  5.d0 / 12.d0

    endif

    f_old = 0.d0
    !---------------------
    ! setup modified wave number, fft indexing and
    ! define fft plans 
    !---------------------


    call compute_modwv 
    call compute_modwv_2nd
    call init_indexing        
    call init_index_2D
    call setup_fft 
    call setup_fft_2D
  
    !-------------
    ! initialize mean pressure grad
    !-------------
    pgm = 1.d0 
    pgm_tee = 1.d0 

  end subroutine setup_rkc

!=========================================================================================================================!

  subroutine rk_step2nd ( istep ) 
    !------------------------------
    ! advance the momentum eqn-
    !     either a full time step for ab2/cn 
    !         or one rk substep for rk3/cn
    !------------------------------

    implicit none 
    integer, intent(in)                   :: istep 
    real, dimension(:,:,:,:), allocatable :: fnl, lun, phi 
    real, dimension(1:2*kl+ku+1,1:n2)     :: wn_u, wn_w, wn_u_phase, wn_w_phase
    real, dimension(1:2*kl+ku+1,1:n2-1)   :: wn_v, wn_v_phase
    integer, dimension(1:n2)              :: ipiv_u, ipiv_w, ipiv_u_phase, ipiv_w_phase
    integer, dimension(1:n2-1)            :: ipiv_v, ipiv_v_phase
    integer                               :: ierr 
    integer                               :: i,j,k,n
    integer                               :: if_start, if_end, kf_start, kf_end
    logical                               :: isfirststep
    real, dimension(3)                    :: coeff_2nd
    real                                  :: hm1, hm2, hp1, hp2
    real                                  :: t_c, uw1, uw2


    !-----
    ! emergency error check , can be commented out at a later date
    !-----
    if ( istep .lt. 1 .or. istep .gt. 3) then 
       write(*,*) 'istep must be either 1,2,3... : ', istep 
       stop
    endif


    allocate ( fnl (1:nvel,is_df:ie_df,1:n2,ks_df:ke_df) )
    allocate ( lun (1:nvel,is_df:ie_df,1:n2,ks_df:ke_df) ) 
    allocate ( phi (1:nvel,is_sf:ie_sf,1:n2,ks_sf:ke_sf) ) 


    fnl =0.d0 
    lun =0.d0 
    phi =0.d0



    !-------------------
    ! form rhs : rhs = u_n + dt*a_i*L(u_n) + dt*g_i*N_n + dt*d_i*N_old 
    ! overwrite q with the rhs
    ! compute convective terms and span/stream diffusion terms
    !-------------------

    call cnv2nd( fnl ) 

    !------ 
    ! add pressure to nonlinear term 
    !------ 

    do k=ks_sf,ke_sf 
    do j=1,n2 
    do i=is_sf,ie_sf 

       fnl(i_u,i,j,k) = fnl(i_u,i,j,k) - pgm 
    enddo 
    enddo 
    enddo 


    !------- 
    ! streamwise/ spanwise diffusion 
    !-------
    call dif2nd ( fnl )                 ! fnl N_n

    !-----------------
    ! use pressure from previous time step in momentum advancement 
    ! cf, dukowicz & dvinksky, van ken, askelvoll & moin 
    ! pressure values kept in p_old 
    !-----------------

    call bcp   (       p_old) 


    if (p_order .eq. 2) then
    call gradp2nd (  phi, p_old) 
    else 
    call gradp (  phi, p_old)
    endif
    

    !----------
    ! compute L(u^n) 
    ! for purposes of operator reuse, the filtering of the 
    ! wn visc terms is moved outside of sub.difwn
    !----------

    call difwn  ( lun, is_df, ie_df, ks_df, ke_df, molc_visc, vt_bar ) 

    if ( itimestep .eq. 1 .and. flat_EE) then
       if ( rank .eq. 0) write(*,*) 'using euler for first time step ...'
       !--------------------
       ! for the first time step using ab2, using euler for convective terms
       !-------------------
       do k=ks_sf,ke_sf 
       do j=1,n2 
       do i=is_sf,ie_sf 
       do n=1,nvel 
          f_old(n,i,j,k) = fnl(n,i,j,k)  
       enddo 
       enddo 
       enddo 
       enddo 
    endif

    if ( itimestep .eq. 3 .and. flat_EE) then
       if ( rank .eq. 0) write(*,*) 'using euler for first time step AB2...'
       !--------------------
       ! for the first time step using ab2, using euler for convective terms
       !-------------------
       do k=ks_sf,ke_sf
       do j=1,n2
       do i=is_sf,ie_sf
       do n=1,nvel
          f_old(n,i,j,k) = fnl(n,i,j,k) - lun(n,i,j,k)
       enddo
       enddo
       enddo
       enddo
    endif


    !----- 
    ! compute rhs, r^n from explicit terms
    !-----

    do k=ks_sf,ke_sf
    do j=1,n2 
    do i=is_sf,ie_sf
    do n=1,nvel 

       if ( n .ne. i_v .or. j .lt. n2 )  then 
          q(n,i,j,k) = q(n,i,j,k) + dt* rkc(istep,1)* lun(n,i,j,k) +  & 
                       +dt* rkc(istep,3)* fnl(n,i,j,k) + dt* rkc(istep,4)* f_old(n,i,j,k) & 
                       -dt*(rkc(istep,1) + rkc(istep,2)) * phi(n,i,j,k) 
       endif 
    enddo
    enddo
    enddo
    enddo


! Boundary effect - 2019.11.03
! Taylor expansion, generalized wall velocity

    if ( wall_bc ) then

      ! wall velocity
      t_c = (time-init_time-dt)*2.d0*pi/wall_per*re
      uw1 = wall_mean + wall_del * sin( t_c ) !explicit
      t_c = (time-init_time)*2*pi/wall_per*re
      uw2 = wall_mean + wall_del * sin( t_c ) !implicit


      ! bottom wall
      i = 1
      hp1 = y(i_u,i+1) - y(i_u,i )
      hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )
      call d2_coeff_2nd( hm1, hp1, coeff_2nd)
 
      do k=ks,ke
        do i=is,ie
          q(i_u,i,1,k) = q(i_u,i,1,k) + dt * rkc(istep,1) * 1.d0/re * &
              coeff_2nd(1) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * uw1
          q(i_u,i,1,k) = q(i_u,i,1,k) + dt * rkc(istep,2) * 1.d0/re * &
              coeff_2nd(1) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * uw2
        enddo
      enddo
      
      ! top wall
      i = n2
      hm1 = y(i_u,i  ) - y(i_u,i-1)
      hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )
      call d2_coeff_2nd( hm1, hp1, coeff_2nd )

      do k=ks,ke
        do i=is,ie
          q(i_u,i,n2,k) = q(i_u,i,n2,k) + dt * rkc(istep,1) * 1.d0/re * &
              coeff_2nd(3) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * uw1
          q(i_u,i,n2,k) = q(i_u,i,n2,k) + dt * rkc(istep,2) * 1.d0/re * & 
              coeff_2nd(3) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * uw2        
        enddo
      enddo

    endif
    
! ...end boundary effect

if ( itimestep .eq. 1 .or. itimestep .eq. 2 ) then 
    !--------
    ! solve (I - dt*b_i*Lun ) 
    !--------
    if ( (.not. jsgs) .or. leddyv_exp) then
       call wn_op ( wn_u, ipiv_u, wn_v, ipiv_v, wn_w, ipiv_w, rkc(istep, 2), 1,1 , molc_visc, 0.d0) 
    else 
       call wn_op ( wn_u, ipiv_u, wn_v, ipiv_v, wn_w, ipiv_w, rkc(istep, 2), 1,1 , molc_visc, vt_bar)
    endif 

    if ( pat_bc ) then
      if ( (.not. jsgs) .or. leddyv_exp) then
         call wn_op_phase ( wn_u_phase, ipiv_u_phase, wn_v_phase, ipiv_v_phase, wn_w_phase, ipiv_w_phase, rkc(istep, 2), 1,1 , molc_visc, 0.d0) 
      else 
         call wn_op_phase ( wn_u_phase, ipiv_u_phase, wn_v_phase, ipiv_v_phase, wn_w_phase, ipiv_w_phase, rkc(istep, 2), 1,1 , molc_visc, vt_bar)
      endif 
    endif      
!endif
! Kim 03.12.23: see comment at line 371 (after dgbtrs calls)

    do k=ks,ke
    do i=is,ie
      ! u
      if ( pat_bc .and. phase(i_u,i,k) .ne. 0 ) then
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_u_phase, 2*kl+ku+1, ipiv_u_phase, q(i_u,i,1:n2  ,k), n2  , ierr)
      else
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_u, 2*kl+ku+1, ipiv_u, q(i_u,i,1:n2  ,k), n2  , ierr)
      endif
      ! v
      if ( pat_bc .and. phase(i_v,i,k) .ne. 0 ) then
        call dgbtrs( 'n', n2-1, kl,ku, 1, wn_v_phase, 2*kl+ku+1, ipiv_v_phase, q(i_v,i,1:n2-1,k), n2-1, ierr)
      else
        call dgbtrs( 'n', n2-1, kl,ku, 1, wn_v, 2*kl+ku+1, ipiv_v, q(i_v,i,1:n2-1,k), n2-1, ierr)
      endif
      ! w
      if ( pat_bc .and. phase(i_w,i,k) .ne. 0 ) then
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_w_phase, 2*kl+ku+1, ipiv_w_phase, q(i_w,i,1:n2  ,k), n2  , ierr)
      else
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_w, 2*kl+ku+1, ipiv_w, q(i_w,i,1:n2  ,k), n2  , ierr)
      endif
    enddo
    enddo

!numerics.f90(324): error #6317: An ENDIF occurred without a corresponding IF
!THEN or ELSE statement.
! commented out endif
endif 
! Kim 03.12.23: this is the wrong endif to comment out
! should be the one above the the dgbtrs() calls


    if ( slip_bc .and. slip_del .ne. 0 ) then
      call update_slip
    endif

    call bcuvw 

    !---------
    ! update f_old 
    !---------

 if (itimestep .ne. 3) then 
    do k=ks_sf,ke_sf
    do j=1,n2 
    do i=is_sf,ie_sf
    do n=1,nvel 
       f_old(n,i,j,k) = fnl(n,i,j,k) 
    enddo
    enddo
    enddo 
    enddo
else 
    do k=ks_sf,ke_sf
    do j=1,n2
    do i=is_sf,ie_sf
    do n=1,nvel
       f_old(n,i,j,k) = fnl(n,i,j,k) - lun(n,i,j,k) 
    enddo
    enddo
    enddo
    enddo
endif


    deallocate ( fnl  ) 
    deallocate ( lun  ) 
    deallocate ( phi  ) 

  end subroutine rk_step2nd

!===================================================================================================================================!

  subroutine compute_AB_f_old

    implicit none 
    real, dimension(:,:,:,:), allocatable :: fnl
    integer                               :: ierr 
    integer                               :: i,j,k,n
    integer                               :: if_start, if_end, kf_start, kf_end

    if ( rank .eq. 0) write(*,*) 'compute_AB_f_old '
    allocate ( fnl (1:nvel,is_df:ie_df,1:n2,ks_df:ke_df) )
    fnl = 0.d0
    call cnv4th_2nd( fnl ) 

    do k=ks_sf,ke_sf 
    do j=1,n2 
    do i=is_sf,ie_sf 

       fnl(i_u,i,j,k) = fnl(i_u,i,j,k) - pgm 
    enddo 
    enddo 
    enddo 

  ! call dif4th ( fnl )                 ! fnl N_n
    call dif2nd ( fnl ) 
    do k=ks_sf,ke_sf
    do j=1,n2
    do i=is_sf,ie_sf
    do n=1,nvel

       f_old(n,i,j,k) = fnl(n,i,j,k)

    enddo
    enddo
    enddo
    enddo

 end subroutine compute_AB_f_old

!===================================================================================================================================!

  subroutine rk_step ( istep ) 


    !------------------------------
    ! advance the momentum eqn-
    !     either a full time step for ab2/cn 
    !         or one rk substep for rk3/cn
    !------------------------------

    implicit none 
    integer, intent(in)                   :: istep 
    real, dimension(:,:,:,:), allocatable :: fnl, lun, phi
    real, dimension(1:2*kl+ku+1,1:n2)     :: wn_u, wn_w, wn_u_slip, wn_w_slip
    real, dimension(1:2*kl+ku+1,1:n2-1)   :: wn_v, wn_v_slip
    integer, dimension(1:n2)              :: ipiv_u, ipiv_w, ipiv_u_slip, ipiv_w_slip
    integer, dimension(1:n2-1)            :: ipiv_v, ipiv_v_slip
    integer                               :: ierr 
    integer                               :: i,j,k,n
    integer                               :: if_start, if_end, kf_start, kf_end
    logical                               :: isfirststep


    !-----
    ! emergency error check , can be commented out at a later date
    !-----
    if ( istep .lt. 1 .or. istep .gt. 3) then 
       write(*,*) 'istep must be either 1,2,3... : ', istep 
       stop
    endif

    allocate ( fnl (1:nvel,is_df:ie_df,1:n2,ks_df:ke_df) )
    allocate ( lun (1:nvel,is_df:ie_df,1:n2,ks_df:ke_df) ) 
    allocate ( phi (1:nvel,is_sf:ie_sf,1:n2,ks_sf:ke_sf) ) 


    !------- 
    ! define the wall normal eddy viscosity 
    ! time constant 
    !  either semi-implicit or explicit time advancement
    !------- 

    fnl =0.d0 
    lun =0.d0 
    phi =0.d0

    !-------------------
    ! form rhs : rhs = u_n + dt*a_i*L(u_n) + dt*g_i*N_n + dt*d_i*N_old 
    ! overwrite q with the rhs
    ! compute convective terms and span/stream diffusion terms
    !-------------------
    if (p_order .eq. 4) then 
    call cnv4th( fnl )
    else
    call cnv4th_2nd( fnl ) 
    endif
    !------ 
    ! add pressure to nonlinear term , fixed mass flux 
    !------ 
    do k=ks_sf,ke_sf
    do j=1,n2
    do i=is_sf,ie_sf

       fnl(i_u,i,j,k) = fnl(i_u,i,j,k) - pgm
    enddo
    enddo
    enddo

    !------- 
    ! streamwise/ spanwise diffusion 
    !-------
    if (p_order .eq. 4) then 
    call dif4th ( fnl )                 ! fnl N_n
    else
    call dif2nd( fnl) 
    endif
    !-----------------
    ! use pressure from previous time step in momentum advancement 
    ! cf, dukowicz & dvinksky, van ken, askelvoll & moin 
    ! pressure values kept in p_old 
    !-----------------

    call bcp   (       p_old) 


    if (p_order .eq. 2) then
    call gradp2nd (  phi, p_old) 
    else 
    call gradp (  phi, p_old)
    endif
    

    !----------
    ! compute L(u^n) 
    ! for purposes of operator reuse, the filtering of the 
    ! wn visc terms is moved outside of sub.difwn
    !----------

    call difwn  ( lun, is_df, ie_df, ks_df, ke_df, molc_visc, vt_bar ) 

    if ( itimestep .eq. 1 .and. flat_EE) then 
       
       if ( rank .eq. 0) write(*,*) 'using euler for first time step ...'
       !--------------------
       ! for the first time step using ab2, using euler for convective terms
       !-------------------
       do k=ks_sf,ke_sf 
       do j=1,n2 
       do i=is_sf,ie_sf 
       do n=1,nvel 
          f_old(n,i,j,k) = fnl(n,i,j,k) 
       enddo 
       enddo 
       enddo 
       enddo 
    endif

    !----- 
    ! compute rhs, r^n from explicit terms
    !-----

    do k=ks_sf,ke_sf
    do j=1,n2 
    do i=is_sf,ie_sf
    do n=1,nvel 

       if ( n .ne. i_v .or. j .lt. n2 )  then 
          q(n,i,j,k) = q(n,i,j,k) + dt* rkc(istep,1)* lun(n,i,j,k)  & 
                       +dt* rkc(istep,3)* fnl(n,i,j,k) + dt* rkc(istep,4)* f_old(n,i,j,k) & 
                       -dt*(rkc(istep,1) + rkc(istep,2)) * phi(n,i,j,k) 
       endif 
    enddo
    enddo
    enddo
    enddo

    !--------
    ! solve (I - dt*b_i*Lun ) 
    !--------
    if ( (.not. jsgs) .or. leddyv_exp) then
       call wn_op ( wn_u, ipiv_u, wn_v, ipiv_v, wn_w, ipiv_w, rkc(istep, 2), 1,1 , molc_visc, 0.d0) 
    else 
       call wn_op ( wn_u, ipiv_u, wn_v, ipiv_v, wn_w, ipiv_w, rkc(istep, 2), 1,1 , molc_visc, vt_bar)
    endif 
      

    !-------
    ! Jongmin Seo edited 2012.06.18. 
    ! solve wall normal direction(i.e. y direction) term for patterend slip boundary conditions
    ! make separate matrix systems for slip boundary condition 
    ! _slip indicates the slip components. 
    !-------    
! numerics.f90(565): error #6317: An ENDIF occurred without a corresponding IF
! THEN or ELSE statement.
!  commneted out  endif

!     endif

    do k=ks,ke
    do i=is,ie
       call dgbtrs( 'n', n2  , kl,ku, 1, wn_u, 2*kl+ku+1, ipiv_u, q(i_u,i,1:n2  ,k), n2  , ierr)
       call dgbtrs( 'n', n2-1, kl,ku, 1, wn_v, 2*kl+ku+1, ipiv_v, q(i_v,i,1:n2-1,k), n2-1, ierr)
       call dgbtrs( 'n', n2  , kl,ku, 1, wn_w, 2*kl+ku+1, ipiv_w, q(i_w,i,1:n2  ,k), n2  , ierr)
    enddo
    enddo

    call bcuvw 

    !---------
    ! update f_old 
    !---------

    do k=ks_sf,ke_sf
    do j=1,n2 
    do i=is_sf,ie_sf
    do n=1,nvel 

       f_old(n,i,j,k) = fnl(n,i,j,k) 

    enddo
    enddo
    enddo 
    enddo


    deallocate ( fnl  ) 
    deallocate ( lun  ) 
    deallocate ( phi  ) 

  end subroutine rk_step

!======================================================================================================!

  subroutine time_stamp ( tstr )

    implicit none
    character(len=3)       ::  tstr

    write(tstr,'(i3.3)') int(time)

  end subroutine time_stamp

!===================================================================================================================================!
  subroutine p_solve (istep ) 

    !---------------------
    ! solve for the pressure correction term, d\phi, 
    ! using fractional step method of d&d 
    !---------------------

    implicit none 
    real                                   :: rhxsq, rhzsq, maxdiv , c1r6 , hy0, tcoeff, tsum 
    real, dimension(2*klp+kup+1,1:n2)      :: pwn 
    integer, dimension(         1:n2)      :: ipiv_p 
    integer                                :: i,j,k,ierr
    integer, intent(in)                    :: istep 
    integer                                :: ilo_in,ilo_out,klo_in,klo_out , jlo_in, jlo_out
    integer                                :: ihi_in,ihi_out,khi_in,khi_out , jhi_in,jhi_out
    integer                                :: iscale, ipermute, nbuf
    real                                   :: norm
!    real                                   :: plan,ipln                        ! at runtime, these must be real*8  
!    integer*8                              :: fplan 
!    complex, dimension(1:n1/2+1,1:n3)      :: fdata_out 
!    real, dimension(1:n1,1:n3)             :: fdata_in 

    p = 0.d0 
    sc= 0.d0 

    rhxsq = ( dble(n1) / (chlx*pi) )**2 
    rhzsq = ( dble(n3) / (chlz*pi) )**2
    hy0   = 2.d0 / dble(n2 ) 


    !------------
    ! compute rhs ( D\hat u ) 
    !------------

    if (p_order .eq. 2) then
    call div2nd ( p(is+1:ie+1, 1:n2, ks+1:ke+1) , maxdiv ) 
    else
    call div4th ( p(is+1:ie+1, 1:n2, ks+1:ke+1) , maxdiv )
    endif

    !--------------
    ! define size of rhs block on each processor 
    ! define parameters for fft_2d
    !-------------- 

    ilo_in = is+1 
    ihi_in = ie+1 
    jlo_in = 1
    jhi_in = n2
    klo_in = ks+1 
    khi_in = ke+1 
    
    ilo_out= is+1 
    ihi_out= ie+1 
    jlo_out= 1
    jhi_out= n2 
    klo_out= ks+1 
    khi_out= ke+1 

    !-----------------------
    ! fourier transform in stream/spanwise directions 
    ! plan is initialized in setup_fft 
    !-----------------------

    do k=klo_in,khi_in
    do j=jlo_in,jhi_in
    do i=ilo_in,ihi_in 
    
       data(in_new(i,j,k)) = cmplx ( p(i,j,k), 0.d0 ) 
    enddo
    enddo
    enddo


    call fft_2d_3decmp( data, data, 1, plan) 


    ! unpack ... 
    do k=klo_out,khi_out
    do j=jlo_out, jhi_out
    do i=ilo_out, ihi_out 
    
       p(i,j,k)  = real ( data(out_new(i,j,k)))
       sc(i,j,k) = aimag( data(out_new(i,j,k)))
    enddo
    enddo
    enddo
    
       
    !---------------------------------
    ! solve system in the wn direction ,( d^2 - (kx^2 + kz^2)I ) \hat p = rhs 
    ! solve for each kx, kz (modified wave numbers ) 
    !---------------------------------

    p  = hy0 * hy0 * p 
    sc = hy0 * hy0 * sc 

    !------------------
    ! pin pressure 
    !------------------

    if ( ilo_out .eq. 1 .and. klo_out .eq. 1 ) then 
       p(1,1,1)  = 0.d0 
       sc(1,1,1) = 0.d0 
    endif 

    do k=klo_out-1,khi_out-1
    do i=ilo_out-1,ihi_out-1
       call p_op   ( pwn, ipiv_p, kxsq(i), kzsq(k) )
       call dgbtrs ( 'n', n2, klp, kup, 1, pwn, 2*klp+kup+1, ipiv_p, p(i+1,1:n2,k+1) , n2, ierr)
       call dgbtrs ( 'n', n2, klp, kup, 1, pwn, 2*klp+kup+1, ipiv_p, sc(i+1,1:n2,k+1), n2, ierr)
    enddo
    enddo




    norm = dble(n1*n3) 


    do k=klo_out,khi_out 
    do j=jlo_out,jhi_out
    do i=ilo_out,ihi_out
    
       data(out_new(i,j,k)) = cmplx( p(i,j,k), sc(i,j,k)) 
    enddo
    enddo
    enddo

    call fft_2d_3decmp ( data, data, -1, ipln ) 


    do k=klo_in, khi_in
    do j=jlo_in, jhi_in
    do i=ilo_in, ihi_in 

       p(i,j,k) = real ( data( finish_new(i,j,k))) / norm 
    enddo
    enddo
    enddo

    call bcp(p)   

    
  end subroutine p_solve

!==========================================================================================================================!

 subroutine update_pressure

    implicit none
    integer                                :: i,j,k,ierr
    real :: tcoeff

    tcoeff = dt


    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1

     p_old(i,j,k) = p(i,j,k)/tcoeff + p_old(i,j,k)

    enddo
    enddo
    enddo

    call bcp(p_old)

 
 end subroutine update_pressure

!==========================================================================================================================!
  subroutine p_solve_FFT_2nd (istep )

    !---------------------
    ! solve for the pressure correction term, d\phi, 
    ! using fractional step method of d&d 
    !--------------------- 
    
    implicit none
    real                                   :: rhxsq, rhzsq, maxdiv , c1r6 , hy0, tcoeff, tsum
    real, dimension(2*klpp+kupp+1,1:n2)      :: pwn_2nd
    integer, dimension(         1:n2)      :: ipiv_p
    integer                                :: i,j,k,ierr
    integer, intent(in)                    :: istep
    integer                                :: ilo_in,ilo_out,klo_in,klo_out , jlo_in, jlo_out
    integer                                :: ihi_in,ihi_out,khi_in,khi_out , jhi_in,jhi_out
    integer                                :: iscale, ipermute, nbuf
    real                                   :: norm
    
    p = 0.d0
    sc= 0.d0

    rhxsq = ( dble(n1) / (chlx*pi) )**2
    rhzsq = ( dble(n3) / (chlz*pi) )**2
    hy0   = 2.d0 / dble(n2 )


    !------------
    ! compute rhs ( D\hat u ) 
    !------------

    if (p_order .eq. 2) then
    call div2nd ( p(is+1:ie+1, 1:n2, ks+1:ke+1) , maxdiv )
    else
    stop
    endif
    
    ilo_in = is+1
    ihi_in = ie+1
    jlo_in = 1
    jhi_in = n2
    klo_in = ks+1
    khi_in = ke+1

    ilo_out= is+1
    ihi_out= ie+1
    jlo_out= 1
    jhi_out= n2
    klo_out= ks+1
    khi_out= ke+1

    do k=klo_in,khi_in
    do j=jlo_in,jhi_in
    do i=ilo_in,ihi_in

       data(in_new(i,j,k)) = cmplx ( p(i,j,k), 0.d0 )
    enddo
    enddo
    enddo


    call fft_2d_3decmp( data, data, 1, plan)


    ! unpack ... 
    do k=klo_out,khi_out
    do j=jlo_out, jhi_out
    do i=ilo_out, ihi_out

       p(i,j,k)  = real ( data(out_new(i,j,k)))
       sc(i,j,k) = aimag( data(out_new(i,j,k)))
    enddo
    enddo
    enddo




    !---------------------------------
    ! solve system in the wn direction ,( d^2 - (kx^2 + kz^2)I ) \hat p = rhs 
    ! solve for each kx, kz (modified wave numbers ) 
    !---------------------------------

    p  = hy0 * hy0 * p
    sc = hy0 * hy0 * sc

    !------------------
    ! pin pressure 
    !------------------

    if ( ilo_out .eq. 1 .and. klo_out .eq. 1 ) then
       p(1,1,1)  = 0.d0
       sc(1,1,1) = 0.d0
    endif

    do k=klo_out-1,khi_out-1
    do i=ilo_out-1,ihi_out-1
       call p_op_2nd  ( pwn_2nd, ipiv_p, kxsq_2nd(i), kzsq_2nd(k) )
       call dgbtrs ( 'n', n2, klpp, kupp, 1, pwn_2nd, 2*klpp+kupp+1, ipiv_p, p(i+1,1:n2,k+1) , n2, ierr)
       call dgbtrs ( 'n', n2, klpp, kupp, 1, pwn_2nd, 2*klpp+kupp+1, ipiv_p, sc(i+1,1:n2,k+1), n2, ierr)
    enddo
    enddo

    norm = dble(n1*n3)

    do k=klo_out,khi_out
    do j=jlo_out,jhi_out
    do i=ilo_out,ihi_out

       data(out_new(i,j,k)) = cmplx( p(i,j,k), sc(i,j,k))
    enddo
    enddo
    enddo

    call fft_2d_3decmp ( data, data, -1, ipln )


    do k=klo_in, khi_in
    do j=jlo_in, jhi_in
    do i=ilo_in, ihi_in

       p(i,j,k) = real ( data( finish_new(i,j,k))) / norm
    enddo
    enddo
    enddo


    call bcp(p)


  end subroutine p_solve_FFT_2nd

!====================================
  subroutine prjsln (istep )             

    !----------------------------------
    ! projects intermediate velocity onto divergence free velocity field 
    ! u^n+1 = u* - \delta t * G(\delta \phi) 
    !---------------------------------- 

    implicit none 
    integer :: i,j,k,n , nout
    integer :: istep 
    real, dimension(:,:,:,:) , allocatable :: gphi 
    real    :: ct, hx, hz, xx,zz
    character(len=3) :: nstr

    allocate ( gphi(1:nvel, is_sf:ie_sf, 1:n2, ks_sf:ke_sf)  ) 

    if (p_order .eq. 2) then
    call gradp2nd ( gphi, p)
    else
    call gradp (  gphi, p)
    endif

    ct = 1.d0 

    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 

       q(i_u,i,j,k) = q(i_u,i,j,k) - ct*gphi(i_u,i,j,k) 
       q(i_w,i,j,k) = q(i_w,i,j,k) - ct*gphi(i_w,i,j,k) 
       
    enddo
    enddo
    enddo


    do k=ks,ke
    do j=1,n2-1 
    do i=is,ie

       q(i_v,i,j,k) = q(i_v,i,j,k) - ct*gphi(i_v,i,j,k) 
    enddo
    enddo
    enddo

    !real                                  :: qs_temp, qss_temp, qt_temp, qts_temp !temoorary value
    !test
    !  qs_temp = 0.d0
    !do j=1,n2
    !  qt_temp = 0.d0
    !  qts_temp = 0.d0
    !do k=ks,ke
    !do i=is,ie
    !      qt_temp = qt_temp + q(i_v,i,j,k)
    !    enddo
    !  enddo
    !  call mpi_allreduce ( qt_temp, qts_temp, 1, MPI_DOUBLE_PRECISION, mpi_sum, mycomm , ierr) 
      !du(n,j) = (qts_temp - qss_temp)/dble(n1)/dble(n3);
    !  qs_temp  = qs_temp + qts_temp
    !enddo
    !   if ( rank .eq. 0) write(*,*) qs_temp
    !test ends
    

   !call bcuvw 
    deallocate( gphi) 


  end subroutine prjsln

!==============================================================================================================================! 
 subroutine dif4th ( fnl ) 

   !-----------------------------------------
   ! explicit diffusion terms; diffusion in streamwise 
   ! and spanwise directions;
   ! computed after convective terms are computed
   !-----------------------------------------

    implicit none 
    real                                                                           :: hx              
    real                                                                           :: hz 
    real, dimension(1:nvel,is_df:ie_df,1:n2,ks_df:ke_df) ,   intent(inout)         :: fnl 
    real                                                                           :: rchx, rchz , rre
    real                                                                           :: c9r16, c1r16
    integer                                                                        :: i,j,k,i_q 
    integer                                                                        :: kf_start, kf_end, if_start, if_end 


    !---- 
    ! define start/end indicies, w/ exp filter 
    !----
    kf_start = ks_df
    kf_end   = ke_df
    if_start = is_df
    if_end   = ie_df


    !---- 
    ! define grid parameters 
    !---- 
    hx = chlx*pi / dble(n1) 
    hz = chlz*pi / dble(n3) 

    !rchx = 1.d0 / (hx*hx*12.d0*re) 
    !rchz = 1.d0 / (hz*hz*12.d0*re) 

    rchx  = (nu_molec + vt_bar) / (hx*hx*12.d0) 
    rchz  = (nu_molec + vt_bar) / (hz*hz*12.d0)

    c9r16 = 9.d0 / 16.d0 
    c1r16 = 1.d0 / 16.d0 



    do k=kf_start,kf_end
    do j=1,n2 
    do i=if_start,if_end 
    do i_q = 1,nvel

       if ( i_q .ne. i_v .or. j .lt. n2 ) then 

          fnl(i_q,i,j,k) = fnl(i_q,i,j,k) & 
                          -rchx * ( -q(i_q,i-2,j,k) + 16.d0*q(i_q,i-1,j,k) - 30.d0*q(i_q,i,j,k) & 
                                    -q(i_q,i+2,j,k) + 16.d0*q(i_q,i+1,j,k) ) & 
                          -rchz * ( -q(i_q,i,j,k-2) + 16.d0*q(i_q,i,j,k-1) - 30.d0*q(i_q,i,j,k) & 
                                    -q(i_q,i,j,k+2) + 16.d0*q(i_q,i,j,k+1) )

       endif

    enddo
    enddo
    enddo
    enddo


  end subroutine dif4th

!========================================================================================================================!
 subroutine dif2nd ( fnl )

   !-----------------------------------------
   ! explicit diffusion terms; diffusion in streamwise 
   ! and spanwise directions;
   ! computed after convective terms are computed
   !-----------------------------------------

    implicit none
    real                                                                           :: hx
    real                                                                           :: hz
    real, dimension(1:nvel,is_df:ie_df,1:n2,ks_df:ke_df) ,   intent(inout)         :: fnl
    real                                                                           :: rchx, rchz , rre
    real                                                                           :: c9r16, c1r16
    integer                                                                        :: i,j,k,i_q
    integer                                                                        :: kf_start, kf_end, if_start, if_end


    !---- 
    ! define start/end indicies, w/ exp filter 
    !----
    kf_start = ks_df
    kf_end   = ke_df
    if_start = is_df
    if_end   = ie_df


    !---- 
    ! define grid parameters 
    !---- 
    hx = chlx*pi / dble(n1)
    hz = chlz*pi / dble(n3)

    !rchx = 1.d0 / (hx*hx*12.d0*re) 
    !rchz = 1.d0 / (hz*hz*12.d0*re) 

    rchx  = (nu_molec + vt_bar) / (hx*hx)
    rchz  = (nu_molec + vt_bar) / (hz*hz)

    c9r16 = 9.d0 / 16.d0
    c1r16 = 1.d0 / 16.d0


    do k=kf_start,kf_end
    do j=1,n2
    do i=if_start,if_end
    do i_q = 1,nvel

       if ( i_q .ne. i_v .or. j .lt. n2 ) then

          fnl(i_q,i,j,k) = fnl(i_q,i,j,k) &
                          -rchx * ( q(i_q,i-1,j,k) - 2.d0*q(i_q,i,j,k)   + q(i_q,i+1,j,k) ) &
                          -rchz * ( q(i_q,i,j,k-1) - 2.d0*q(i_q,i,j,k)   + q(i_q,i,j,k+1) )

       endif

    enddo
    enddo
    enddo
    enddo


  end subroutine dif2nd

!========================================================================================================================!

    subroutine div2nd ( div , maxd )

      !---------------------
      ! compute divergence of the flow field in q  
      ! returns the max divergence of the field in maxd 
      !---------------------

    implicit none
    real, dimension(is+1:ie+1, 1:n2, ks+1:ke+1)  , intent(out)       :: div
    real,                                  intent(out)       :: maxd
    integer                                                  :: i,j,k,n
    real                                                     :: rhx, rhz, rjy, rh3x, rh3z, rj3y, hy0, jac, c9r8, c1r8, c1r24
    real                                                     :: dudx, dvdy, dwdz


    div  = 0.d0
    rhx  = dble(n1) / ( chlx*pi )
    rhz  = dble(n3) / ( chlz*pi )
    rh3x = rhx / 3.d0
    rh3z = rhz / 3.d0
    hy0  = 2.d0/ dble(n2)
    c9r8 = 9.d0/ 8.d0
    c1r8 = 1.d0/ 8.d0
    c1r24= 1.d0/ 24.d0


    do k=ks,ke
    do j=1,n2
    do i=is,ie

       jac  = jacb(i_u, j )
       rjy  = 1.d0/hy0/jac

       dudx = rhx * ( q(i_u,i+1,j,k) - q(i_u,i ,j,k) )
       dwdz = rhz * ( q(i_w,i,j,k+1) - q(i_w,i,j,k ) )
       dvdy = rjy * ( q(i_v,i,j  ,k) - q(i_v,i,j-1,k)) 

       div(i+1,j  ,k+1) = dudx  + dwdz  + dvdy

    enddo
    enddo
    enddo

    maxd = maxval ( abs(div) )

  end subroutine div2nd

!========================================================================================================================!

    subroutine div4th ( div , maxd ) 

      !---------------------
      ! compute divergence of the flow field in q  
      ! returns the max divergence of the field in maxd 
      !---------------------

    implicit none 
    real, dimension(is+1:ie+1, 1:n2, ks+1:ke+1)  , intent(out)       :: div
    real,                                  intent(out)       :: maxd 
    integer                                                  :: i,j,k,n 
    real                                                     :: rhx, rhz, rjy, rh3x, rh3z, rj3y, hy0, jac, c9r8, c1r8, c1r24
    real                                                     :: dudx, dvdy, dwdz 


    div  = 0.d0 
    rhx  = dble(n1) / ( chlx*pi ) 
    rhz  = dble(n3) / ( chlz*pi ) 
    rh3x = rhx / 3.d0 
    rh3z = rhz / 3.d0 
    hy0  = 2.d0/ dble(n2) 
    c9r8 = 9.d0/ 8.d0 
    c1r8 = 1.d0/ 8.d0 
    c1r24= 1.d0/ 24.d0 


    do k=ks,ke
    do j=1,n2 
    do i=is,ie 


       jac  = jacb(i_u, j ) 
       rjy  = 1.d0/hy0/jac
       rj3y = 1.d0/3.d0/hy0/jac


       dudx = c9r8 * rhx * ( q(i_u,i+1,j,k) - q(i_u,i ,j,k) ) - c1r8 * rh3x* ( q(i_u,i+2,j,k) - q(i_u,i-1,j,k) )
       dwdz = c9r8 * rhz * ( q(i_w,i,j,k+1) - q(i_w,i,j,k ) ) - c1r8 * rh3z* ( q(i_w,i,j,k+2) - q(i_w,i,j,k-1) )
       dvdy = c9r8 * rjy * ( q(i_v,i,j  ,k) - q(i_v,i,j-1,k)) - c1r8 * rj3y* ( q(i_v,i,j+1,k) - q(i_v,i,j-2,k) )

     if ( j .eq. 1 ) then
         dvdy = (c9r8/jac  + cf(i_v,-1,1)/24.d0/jac)* q(i_v,i,j,k)/hy0  &
               +(-c1r24/jac+ cf(i_v,-1,2)/24.d0/jac)* q(i_v,i,j+1,k)/hy0     !rjy* q(i_v,i,j,k) 
      elseif (j .eq. n2) then
         dvdy =-(c9r8/jac  + cf(i_v,-1,1)/24.d0/jac)* q(i_v,i,j-1,k)/hy0  &
               +(c1r24/jac - cf(i_v,-1,2)/24.d0/jac)* q(i_v,i,j-2,k)/hy0     !-rjy* q(i_v,i,j-1,k) 
      else 
        dvdy = c9r8 * rjy * ( q(i_v,i,j  ,k) - q(i_v,i,j-1,k)) - c1r8 * rj3y* ( q(i_v,i,j+1,k) - q(i_v,i,j-2,k) )
      endif


      div(i+1,j  ,k+1) = dudx  + dwdz  + dvdy 
       
    enddo
    enddo
    enddo

    maxd = maxval ( abs(div) ) 


  end subroutine div4th

   !===============================================================================================================================! 

   subroutine cnv2nd( fnl )
   
  !--------------------------------------
  ! compute convective terms using skew-symmetric form (cf. vasilyev, jcp, 2000) 
  ! to fourth order accuracy
  !--------------------------------------

   implicit none
   real, dimension(1:nvel,is_df:ie_df, 1:n2,ks_df:ke_df), intent(inout)       :: fnl
   real, dimension(:,:,:,:) , allocatable                         :: fdiv!  , fadv
   integer                                                        :: i,j,k,n, if_start, if_end, kf_start, kf_end
   real                                                           :: rhx, rhz, rh3x, rh3z, hy0, jac , rjy, rj3y
   real                                                           :: c9r8, c9r16, c1r8, c1r16, c1r2
   real                                                           :: c8r8  ! debug 
   real                                                           :: rhy, rh3y, dudx, dwdz

   !----------
   ! init halo indexing
   !---------

   if ( jxzf ) then
   if_start = is_df
   if_end   = ie_df
   kf_start = ks_df
   kf_end   = ke_df
   else
   if_start = is 
   if_end   = ie 
   kf_start = ks
   kf_end   = ke
   endif

   !--------
   ! constants 
   !--------

   c9r8  = 9.0d0/8.0d0
   c9r16 = 9.0d0/16.d0
   c1r8  = 1.0d0/8.0d0
   c1r16 = 1.0d0/16.d0
   c1r2  = 1.0d0/2.0d0
   c8r8  = 1.0d0


   rhx   = dble(n1) / (chlx*pi)
   rhz   = dble(n3) / (chlz*pi)
   rh3x  = rhx/3.0d0
   rh3z  = rhz/3.0d0
   hy0   = 2.d0/ dble(n2)
   rhy   = dble(n2) / 2.d0
   rh3y  = rhy/3.d0


  !------------------
                                                                           
  if ( osverify)  pgm = 2.d0/re
   if ( mmsverify) pgm = 0.d0


   allocate ( fdiv(1:nvel, if_start:if_end , 1:n2  , kf_start:kf_end) )
!   allocate ( fadv(1:nvel, if_start:if_end , 1:n2  , kf_start:kf_end) )

   fdiv = 0.d0
!   fadv = 0.d0


!--------------------
! divergence form  
!--------------------

!... compute convective terms in u 


do k=kf_start, kf_end
do j=1,n2
do i=if_start, if_end

jac  = jacb(i_u,j)
rjy  = 1.d0/hy0/jac
rj3y = 1.d0/3.d0/hy0/jac


! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! ----

!.... fdiv(u) = d/dx(uu ) + ... 
fdiv(i_u,i,j,k) = rhx* ( (c1r2*( q(i_u,i  ,j,k) + q(i_u,i+1,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i+1,j,k))  &
- (c1r2*( q(i_u,i  ,j,k) + q(i_u,i-1,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i-1,j,k)) )


!... fdiv(u) = ... + d/dz(uw )  
fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rhz* ( (c1r2*( q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k+1))  &
- (c1r2*( q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k-1)) )


!.... fdiv(u) = + .... d/dy( uv) + ... 
fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j-1,k))  )


enddo
enddo
enddo




!.... compute convective terms in v  (wn ) 

do k=kf_start, kf_end
do j=1,n2-1
do i=if_start, if_end

jac  = jacb(i_v,j)
rjy  = 1.d0/hy0/jac
rj3y = 1.d0/3.d0/hy0/jac

!.... fdiv(v) = d/dx(vu) + ...  
fdiv(i_v,i,j,k) = rhx* ( (c1r2*( q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i+1,j,k))  &
- (c1r2*( q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i-1,j,k)) )

!... fdiv(v) = ... + d/dz(vw) 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rhz* ( (c1r2*( q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k+1))  &
- (c1r2*( q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k-1)) )

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j-1,k)) )

enddo
enddo
enddo


!.... compute convective terms in w (span) 

do k=kf_start, kf_end
do j=1,n2
do i=if_start, if_end

jac  = jacb(i_u,j)
rjy  = 1.d0/hy0/jac
rj3y = 1.d0/3.d0/hy0/jac

!.... fdiv(w) = d/dx(wu) + .... 
fdiv(i_w,i,j,k) = rhx* ( (c1r2*( q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i+1,j,k))  &
- (c1r2*( q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i-1,j,k)) )


!.... fdiv(w) = ... + d/dz(ww)       
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rhz* ( (c1r2*( q(i_w,i,j,k  ) + q(i_w,i,j,k+1)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k+1))  &
- (c1r2*( q(i_w,i,j,k  ) + q(i_w,i,j,k-1)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k-1))) 

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j-1,k)) )

enddo
enddo
enddo


!--------------
! advective form 
!--------------

!do k=kf_start, kf_end
!do j=1,n2
!do i=if_start, if_end

!rjy = 1.d0 / jacb(i_u,j)

!.... f(i_u) = u du/dx + .... 
!fadv(i_u,i,j,k) = c1r2* ( (c1r2* (q(i_u,i  ,j,k) + q(i_u,i-1,j,k)))* rhx*   (q(i_u,i  ,j,k) - q(i_u,i-1,j,k)) &
!+(c1r2* (q(i_u,i  ,j,k) + q(i_u,i+1,j,k)))* rhx*   (q(i_u,i+1,j,k) - q(i_u,i  ,j,k)) ) 

!.... f(i_u) =  ... + w du/dz 
!fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + &
! c1r2* ( (c1r2* (q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)))* rhz*   (q(i_u,i  ,j,k+1) - q(i_u,i  ,j,k  )) &
!+(c1r2* (q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )))* rhz*   (q(i_u,i  ,j,k  ) - q(i_u,i  ,j,k-1)) ) 

!.... f(i_u) =  + ... v du/dy + ... 
!fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
!c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q(i_u,i  ,j+1,k) - q(i_u,i  ,j  ,k)) &
!+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q(i_u,i  ,j  ,k) - q(i_u,i  ,j-1,k)) ))

!---- 

!enddo
!enddo
!enddo


!.... fadv(i_v) 


!do k=kf_start, kf_end
!do j=1,n2-1
!do i=if_start, if_end

!rjy = 1.d0 / jacb(i_v,j)

!.... f(i_v) =  u dv/dx 
!fadv(i_v,i,j,k) =  c1r2* ( (c1r2* (q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)))* rhx*   (q(i_v,i  ,j  ,k) - q(i_v,i-1,j  ,k)) &
!+(c1r2* (q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)))* rhx*   (q(i_v,i+1,j  ,k) - q(i_v,i  ,j  ,k)) ) 

!.... f(i_v) = ... + w dv/dz 
!fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + &
!c1r2* ( (c1r2* (q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )))* rhz*   (q(i_v,i,j  ,k  ) - q(i_v,i,j  ,k-1)) &
!+(c1r2* (q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)))* rhz*   (q(i_v,i,j  ,k+1) - q(i_v,i,j  ,k  )) ) 


!.... f(i_v) = ... + v dv/dy + .... 
!fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
!c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q(i_v,i  ,j+1,k) - q(i_v,i  ,j  ,k)) &
!+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q(i_v,i  ,j  ,k) - q(i_v,i  ,j-1,k)) ) )

!enddo
!enddo
!enddo


!.... fadv(i_w) 


!do k=kf_start, kf_end
!do j=1,n2
!do i=if_start, if_end

!rjy = 1.d0 / jacb(i_u,j)

!.... fadv(i_w) = u dw/dx + ....
!fadv(i_w,i,j,k) =  c1r2* ( (c1r2* (q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)))* rhx*   (q(i_w,i  ,j,k  ) - q(i_w,i-1,j,k  )) &
!+(c1r2* (q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)))* rhx*   (q(i_w,i+1,j,k  ) - q(i_w,i  ,j,k  )) ) 


!.... fadv(i_w) = ... + w dw/dz 
!fadv(i_w,i,j,k) = fadv(i_w,i,j,k) + &
! c1r2* ( (c1r2* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-1)))*rhz*(q(i_w,i  ,j,k  ) - q(i_w,i  ,j,k-1)) &
!+(c1r2* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+1)))*rhz*(q(i_w,i  ,j,k+1) - q(i_w,i  ,j,k  )) ) 


!.... fadv(i_w) = ...+ v dw/dy + ... 
!fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
!c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q(i_w,i,j+1,k  ) - q(i_w,i,j  ,k  )) &
!+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q(i_w,i,j  ,k  ) - q(i_w,i,j-1,k  )) ) )
!enddo
!enddo
!enddo




   !------------
   ! skew = 0.5 div + 0.5 adv 
   !------------

   do k=kf_start, kf_end
    do j=1,n2
     do i=if_start, if_end
      do n=1,nvel

           fnl(n,i,j,k) = fdiv(n,i,j,k)                                 ! divergence form 
     !       fnl(n,i,j,k) = fadv(n,i,j,k)                                 ! advective  form 
   !  fnl(n,i,j,k) = 0.5d0 * fdiv(n,i,j,k) + 0.5d0 * fadv(n,i,j,k)  ! skew-sym   form 

      enddo
     enddo
    enddo
   enddo


    !.... clean-up temporary memory 
   deallocate(fdiv)
!   deallocate(fadv)


  end subroutine cnv2nd

   !===============================================================================================================================! 

   subroutine cnv4th_2nd( fnl )
   
  !--------------------------------------
  ! compute convective terms using skew-symmetric form (cf. vasilyev, jcp, 2000) 
  ! to fourth order accuracy
  !--------------------------------------

   implicit none
   real, dimension(1:nvel,is_df:ie_df, 1:n2,ks_df:ke_df), intent(inout)       :: fnl
   real, dimension(:,:,:,:) , allocatable                         :: fdiv  , fadv
   integer                                                        :: i,j,k,n, if_start, if_end, kf_start, kf_end
   real                                                           :: rhx, rhz, rh3x, rh3z, hy0, jac , rjy, rj3y
   real                                                           :: c9r8, c9r16, c1r8, c1r16, c1r2
   real                                                           :: c8r8  ! debug 
   real                                                           :: rhy, rh3y, dudx, dwdz

   !----------
   ! init halo indexing
   !---------

   if ( jxzf ) then
   if_start = is_df
   if_end   = ie_df
   kf_start = ks_df
   kf_end   = ke_df
   else
   if_start = is 
   if_end   = ie 
   kf_start = ks
   kf_end   = ke
   endif

   !--------
   ! constants 
   !--------

   c9r8  = 9.0d0/8.0d0
   c9r16 = 9.0d0/16.d0
   c1r8  = 1.0d0/8.0d0
   c1r16 = 1.0d0/16.d0
   c1r2  = 1.0d0/2.0d0
   c8r8  = 1.0d0


   rhx   = dble(n1) / (chlx*pi)
   rhz   = dble(n3) / (chlz*pi)
   rh3x  = rhx/3.0d0
   rh3z  = rhz/3.0d0
   hy0   = 2.d0/ dble(n2)
   rhy   = dble(n2) / 2.d0
   rh3y  = rhy/3.d0




  !------------------
                                                                           
  if ( osverify)  pgm = 2.d0/re
   if ( mmsverify) pgm = 0.d0


   allocate ( fdiv(1:nvel, if_start:if_end , 1:n2  , kf_start:kf_end) )
   allocate ( fadv(1:nvel, if_start:if_end , 1:n2  , kf_start:kf_end) )

   fdiv = 0.d0
   fadv = 0.d0


!--------------------
! divergence form  
!--------------------

!... compute convective terms in u 


do k=kf_start, kf_end
do j=1,n2
do i=if_start, if_end

jac  = jacb(i_u,j)
rjy  = 1.d0/hy0/jac
rj3y = 1.d0/3.d0/hy0/jac


! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! ----

!.... fdiv(u) = d/dx(uu ) + ... 
fdiv(i_u,i,j,k) = c9r8*rhx* ( (c9r16*( q(i_u,i  ,j,k) + q(i_u,i+1,j,k)) &
-c1r16*( q(i_u,i-1,j,k) + q(i_u,i+2,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i+1,j,k))  &
- (c9r16*( q(i_u,i  ,j,k) + q(i_u,i-1,j,k)) &
-c1r16*( q(i_u,i+1,j,k) + q(i_u,i-2,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i-1,j,k)) )&
-c1r8*rh3x*( (c9r16*( q(i_u,i+1,j,k) + q(i_u,i+2,j,k)) &
-c1r16*( q(i_u,i  ,j,k) + q(i_u,i+3,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i+3,j,k))  &
- (c9r16*( q(i_u,i-1,j,k) + q(i_u,i-2,j,k)) &
-c1r16*( q(i_u,i  ,j,k) + q(i_u,i-3,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i-3,j,k)) )


!... fdiv(u) = ... + d/dz(uw )  
fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
c9r8*rhz* ( (c9r16*( q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)) &
-c1r16*( q(i_w,i+1,j,k+1) + q(i_w,i-2,j,k+1)))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k+1))  &
- (c9r16*( q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )) &
-c1r16*( q(i_w,i+1,j,k  ) + q(i_w,i-2,j,k  )))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k-1)) )&
-c1r8*rh3z*( (c9r16*( q(i_w,i  ,j,k+2) + q(i_w,i-1,j,k+2)) &
-c1r16*( q(i_w,i+1,j,k+2) + q(i_w,i-2,j,k+2)))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k+3))  &
- (c9r16*( q(i_w,i  ,j,k-1) + q(i_w,i-1,j,k-1)) &
-c1r16*( q(i_w,i+1,j,k-1) + q(i_w,i-2,j,k-1)))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k-3)) )


if ( j .eq. 1 ) then

!.... fdiv(u) = + .... d/dy( uv) + ... 

fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j-1,k))  )

elseif ( j .eq. 2 ) then

!.... fdiv(u) = + .... d/dy( uv) + ... 

fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j-1,k))  )

elseif (j .eq. n2-1) then

!.... fdiv(u) = + .... d/dy( uv) + ... 

fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j-1,k))  )

elseif (j .eq. n2) then

!.... fdiv(u) = + .... d/dy( uv) + ... 

fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i  ,j-1,k))  )

else

! ----
! Interiror: high order scheme (4th order)
! ----

!.... fdiv(u) = + .... d/dy( uv) + ... 
fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
c9r8*rjy* ( (c9r16*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)) &
-c1r16*( q(i_v,i+1,j  ,k) + q(i_v,i-2,j  ,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i,j+1,k))  &
- (c9r16*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)) &
-c1r16*( q(i_v,i+1,j-1,k) + q(i_v,i-2,j-1,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i,j-1,k)) )&
-c1r8*rj3y*( (c9r16*( q(i_v,i  ,j+1,k) + q(i_v,i-1,j+1,k)) &
-c1r16*( q(i_v,i+1,j+1,k) + q(i_v,i-2,j+1,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i,j+3,k))  &
- (c9r16*( q(i_v,i  ,j-2,k) + q(i_v,i-1,j-2,k)) &
-c1r16*( q(i_v,i+1,j-2,k) + q(i_v,i-2,j-2,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i,j-3,k)) )



!---- 
! deprecated; the pressure gradient is 
! now added elsewhere 
!----

!fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) - pgm 
endif

enddo
enddo
enddo




!.... compute convective terms in v  (wn ) 

do k=kf_start, kf_end
do j=1,n2-1
do i=if_start, if_end

jac  = jacb(i_v,j)
rjy  = 1.d0/hy0/jac
rj3y = 1.d0/3.d0/hy0/jac

!.... fdiv(v) = d/dx(vu) + ...  
fdiv(i_v,i,j,k) = c9r8*rhx* ( (c9r16*( q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)) &
-c1r16*( q(i_u,i+1,j-1,k) + q(i_u,i+1,j+2,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i+1,j,k))  &
- (c9r16*( q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)) &
-c1r16*( q(i_u,i  ,j-1,k) + q(i_u,i  ,j+2,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i-1,j,k)) )&
-c1r8*rh3x*( (c9r16*( q(i_u,i+2,j  ,k) + q(i_u,i+2,j+1,k)) &
-c1r16*( q(i_u,i+2,j-1,k) + q(i_u,i+2,j+2,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i+3,j,k))  &
- (c9r16*( q(i_u,i-1,j  ,k) + q(i_u,i-1,j+1,k)) &
-c1r16*( q(i_u,i-1,j-1,k) + q(i_u,i-1,j+2,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i-3,j,k)) )


!... fdiv(v) = ... + d/dz(vw) 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
c9r8*rhz* ( (c9r16*( q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)) &
-c1r16*( q(i_w,i,j-1,k+1) + q(i_w,i,j+2,k+1)))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k+1))  &
- (c9r16*( q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )) &
-c1r16*( q(i_w,i,j-1,k  ) + q(i_w,i,j+2,k  )))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k-1)) )&
-c1r8*rh3z*( (c9r16*( q(i_w,i,j  ,k+2) + q(i_w,i,j+1,k+2)) &
-c1r16*( q(i_w,i,j-1,k+2) + q(i_w,i,j+2,k+2)))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k+3))  &
- (c9r16*( q(i_w,i,j  ,k-1) + q(i_w,i,j+1,k-1)) &
-c1r16*( q(i_w,i,j-1,k-1) + q(i_w,i,j+2,k-1)))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k-3)) )


! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----

if ( j .eq. 1 ) then

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j-1,k)) )

elseif ( j .eq. 2 ) then

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j-1,k)) )

elseif (j .eq. n2-2) then

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j-1,k)) )

elseif (j .eq. n2-1) then

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j-1,k)) )

else


!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
c9r8*rjy* ( (c9r16*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)) &
-c1r16*( q(i_v,i,j-1,k) + q(i_v,i,j+2,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j+1,k))  &
- (c9r16*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)) &
-c1r16*( q(i_v,i,j+1,k) + q(i_v,i,j-2,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j-1,k)) )&
-c1r8*rj3y*( (c9r16*( q(i_v,i,j+1,k) + q(i_v,i,j+2,k)) &
-c1r16*( q(i_v,i,j  ,k) + q(i_v,i,j+3,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j+3,k))  &
- (c9r16*( q(i_v,i,j-1,k) + q(i_v,i,j-2,k)) &
-c1r16*( q(i_v,i,j  ,k) + q(i_v,i,j-3,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j-3,k)) )

endif
enddo
enddo
enddo


!.... compute convective terms in w (span) 

do k=kf_start, kf_end
do j=1,n2
do i=if_start, if_end

jac  = jacb(i_u,j)
rjy  = 1.d0/hy0/jac
rj3y = 1.d0/3.d0/hy0/jac

!.... fdiv(w) = d/dx(wu) + .... 
fdiv(i_w,i,j,k) = c9r8*rhx* ( (c9r16*( q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)) &
-c1r16*( q(i_u,i+1,j,k+1) + q(i_u,i+1,j,k-2)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i+1,j,k))  &
- (c9r16*( q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)) &
-c1r16*( q(i_u,i  ,j,k+1) + q(i_u,i  ,j,k-2)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i-1,j,k)) )&
-c1r8*rh3x*( (c9r16*( q(i_u,i+2,j,k  ) + q(i_u,i+2,j,k-1)) &
-c1r16*( q(i_u,i+2,j,k+1) + q(i_u,i+2,j,k-2)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i+3,j,k))  &
- (c9r16*( q(i_u,i-1,j,k  ) + q(i_u,i-1,j,k-1)) &
-c1r16*( q(i_u,i-1,j,k+1) + q(i_u,i-1,j,k-2)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i-3,j,k)) )

!.... fdiv(w) = ... + d/dz(ww)       
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
c9r8*rhz* ( (c9r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k+1)) &
-c1r16*( q(i_w,i,j,k-1) + q(i_w,i,j,k+2)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k+1))  &
- (c9r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k-1)) &
-c1r16*( q(i_w,i,j,k+1) + q(i_w,i,j,k-2)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k-1))) &
-c1r8*rh3z*( (c9r16*( q(i_w,i,j,k+1) + q(i_w,i,j,k+2)) &
-c1r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k+3)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k+3))  &
- (c9r16*( q(i_w,i,j,k-1) + q(i_w,i,j,k-2)) &
-c1r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k-3)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k-3)))



! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----

if ( j .eq. 1 ) then

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j-1,k)) )

elseif ( j .eq. 2 ) then

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j-1,k)) )

elseif (j .eq. n2-1) then

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j-1,k)) )

elseif (j .eq. n2) then

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j-1,k)) )

else


!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
c9r8*rjy* ( (c9r16*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) &
-c1r16*( q(i_v,i,j  ,k+1) + q(i_v,i,j  ,k-2)))* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j+1,k))  &
- (c9r16*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) &
-c1r16*( q(i_v,i,j-1,k+1) + q(i_v,i,j-1,k-2)))* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j-1,k)) )&
-c1r8*rj3y*( (c9r16*( q(i_v,i,j+1,k  ) + q(i_v,i,j+1,k-1)) &
-c1r16*( q(i_v,i,j+1,k+1) + q(i_v,i,j+1,k-2)))* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j+3,k))  &
- (c9r16*( q(i_v,i,j-2,k  ) + q(i_v,i,j-2,k-1)) &
-c1r16*( q(i_v,i,j-2,k+1) + q(i_v,i,j-2,k-2)))* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j-3,k)) )
endif
enddo
enddo
enddo


!--------------
! advective form 
!--------------

do k=kf_start, kf_end
do j=1,n2
do i=if_start, if_end

rjy = 1.d0 / jacb(i_u,j)

!.... f(i_u) = u du/dx + .... 
fadv(i_u,i,j,k) = c9r8 * c1r2* ( (c9r16* (q(i_u,i  ,j,k) + q(i_u,i-1,j,k)) - c1r16*(q(i_u,i-2,j,k) + q(i_u,i+1,j,k)))* &
rhx*   (q(i_u,i  ,j,k) - q(i_u,i-1,j,k)) &
+(c9r16* (q(i_u,i  ,j,k) + q(i_u,i+1,j,k)) - c1r16*(q(i_u,i-1,j,k) + q(i_u,i+2,j,k)))* &
rhx*   (q(i_u,i+1,j,k) - q(i_u,i  ,j,k)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_u,i-2,j,k) + q(i_u,i-1,j,k)) - c1r16*(q(i_u,i-3,j,k) + q(i_u,i  ,j,k)))* &
rh3x*  (q(i_u,i  ,j,k) - q(i_u,i-3,j,k)) &
+(c9r16* (q(i_u,i+2,j,k) + q(i_u,i+1,j,k)) - c1r16*(q(i_u,i+3,j,k) + q(i_u,i  ,j,k)))* &
rh3x*  (q(i_u,i+3,j,k) - q(i_u,i  ,j,k)) )

!.... f(i_u) =  ... + w du/dz 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + &
c9r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)) - c1r16*(q(i_w,i-2,j,k+1) + q(i_w,i+1,j,k+1)))* &
rhz*   (q(i_u,i  ,j,k+1) - q(i_u,i  ,j,k  )) &
+(c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )) - c1r16*(q(i_w,i-2,j,k  ) + q(i_w,i+1,j,k  )))* &
rhz*   (q(i_u,i  ,j,k  ) - q(i_u,i  ,j,k-1)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k-1) + q(i_w,i-1,j,k-1)) - c1r16*(q(i_w,i-2,j,k-1) + q(i_w,i+1,j,k-1)))* &
rh3z*  (q(i_u,i  ,j,k  ) - q(i_u,i  ,j,k-3)) &
+(c9r16* (q(i_w,i  ,j,k+2) + q(i_w,i-1,j,k+2)) - c1r16*(q(i_w,i-2,j,k+2) + q(i_w,i+1,j,k+2)))* &
rh3z*  (q(i_u,i  ,j,k+3) - q(i_u,i  ,j,k  )) )



! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----


if ( j .eq. 1 ) then

!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q(i_u,i  ,j+1,k) - q(i_u,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q(i_u,i  ,j  ,k) - q(i_u,i  ,j-1,k)) ))

elseif ( j .eq. 2 ) then

!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q(i_u,i  ,j+1,k) - q(i_u,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q(i_u,i  ,j  ,k) - q(i_u,i  ,j-1,k)) ))

elseif ( j .eq. n2-1 ) then

!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q(i_u,i  ,j+1,k) - q(i_u,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q(i_u,i  ,j  ,k) - q(i_u,i  ,j-1,k)) ))

elseif ( j .eq. n2 ) then

!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q(i_u,i  ,j+1,k) - q(i_u,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q(i_u,i  ,j  ,k) - q(i_u,i  ,j-1,k)) ))

else


!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c9r8 * c1r2*(  (c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)) - c1r16*(q(i_v,i-2,j  ,k) + q(i_v,i+1,j  ,k)))* &
rhy*   (q(i_u,i  ,j+1,k) - q(i_u,i  ,j  ,k)) &
+(c9r16* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)) - c1r16*(q(i_v,i-2,j-1,k) + q(i_v,i+1,j-1,k)))* &
rhy*   (q(i_u,i  ,j  ,k) - q(i_u,i  ,j-1,k)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_v,i  ,j-2,k) + q(i_v,i-1,j-2,k)) - c1r16*(q(i_v,i-2,j-2,k) + q(i_v,i+1,j-2,k)))* &
rh3y*  (q(i_u,i  ,j  ,k) - q(i_u,i  ,j-3,k)) &
+(c9r16* (q(i_v,i  ,j+1,k) + q(i_v,i-1,j+1,k)) - c1r16*(q(i_v,i-2,j+1,k) + q(i_v,i+1,j+1,k)))* &
rh3y*  (q(i_u,i  ,j+3,k) - q(i_u,i  ,j  ,k)) ) )

endif
!---- 
! deprecated;; pressure gradient is now 
! added elsewhere 
!---

!fadv(i_u,i,j,k) = fadv(i_u,i,j,k) - pgm 



enddo
enddo
enddo


!.... fadv(i_v) 


do k=kf_start, kf_end
do j=1,n2-1
do i=if_start, if_end

rjy = 1.d0 / jacb(i_v,j)

!.... f(i_v) =  u dv/dx 
fadv(i_v,i,j,k) = c9r8 * c1r2* ( (c9r16* (q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)) - c1r16*(q(i_u,i  ,j-1,k) + q(i_u,i  ,j+2,k)))* &
rhx*   (q(i_v,i  ,j  ,k) - q(i_v,i-1,j  ,k)) &
+(c9r16* (q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)) - c1r16*(q(i_u,i+1,j-1,k) + q(i_u,i+1,j+2,k)))* &
rhx*   (q(i_v,i+1,j  ,k) - q(i_v,i  ,j  ,k)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_u,i-1,j  ,k) + q(i_u,i-1,j+1,k)) - c1r16*(q(i_u,i-1,j-1,k) + q(i_u,i-1,j+2,k)))* &
rh3x*  (q(i_v,i  ,j  ,k) - q(i_v,i-3,j  ,k)) &
+(c9r16* (q(i_u,i+2,j  ,k) + q(i_u,i+2,j+1,k)) - c1r16*(q(i_u,i+2,j-1,k) + q(i_u,i+2,j+2,k)))* &
rh3x*  (q(i_v,i+3,j  ,k) - q(i_v,i  ,j  ,k)) )

!.... f(i_v) = ... + w dv/dz 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + &
c9r8 * c1r2* ( (c9r16* (q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )) - c1r16*(q(i_w,i,j-1,k  ) + q(i_w,i,j+2,k  )))* &
rhz*   (q(i_v,i,j  ,k  ) - q(i_v,i,j  ,k-1)) &
+(c9r16* (q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)) - c1r16*(q(i_w,i,j-1,k+1) + q(i_w,i,j+2,k+1)))* &
rhz*   (q(i_v,i,j  ,k+1) - q(i_v,i,j  ,k  )) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_w,i,j  ,k-1) + q(i_w,i,j+1,k-1)) - c1r16*(q(i_w,i,j-1,k-1) + q(i_w,i,j+2,k-1)))* &
rh3z*  (q(i_v,i,j  ,k  ) - q(i_v,i,j  ,k-3)) &
+(c9r16* (q(i_w,i,j  ,k+2) + q(i_w,i,j+1,k+2)) - c1r16*(q(i_w,i,j-1,k+2) + q(i_w,i,j+2,k+2)))* &
rh3z*  (q(i_v,i,j  ,k+3) - q(i_v,i,j  ,k  )) )


! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----


if ( j .eq. 1 ) then

!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q(i_v,i  ,j+1,k) - q(i_v,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q(i_v,i  ,j  ,k) - q(i_v,i  ,j-1,k)) ) )

elseif ( j .eq. 2 ) then

!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q(i_v,i  ,j+1,k) - q(i_v,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q(i_v,i  ,j  ,k) - q(i_v,i  ,j-1,k)) ) )

elseif ( j .eq. n2-2 ) then

!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q(i_v,i  ,j+1,k) - q(i_v,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q(i_v,i  ,j  ,k) - q(i_v,i  ,j-1,k)) ) )


elseif ( j .eq. n2-1 ) then

!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q(i_v,i  ,j+1,k) - q(i_v,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q(i_v,i  ,j  ,k) - q(i_v,i  ,j-1,k)) ) )

else


!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c9r8 * c1r2* ( (c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)) - c1r16*(q(i_v,i  ,j-1,k) + q(i_v,i  ,j+2,k)))* &
rhy*   (q(i_v,i  ,j+1,k) - q(i_v,i  ,j  ,k)) &
+(c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)) - c1r16*(q(i_v,i  ,j+1,k) + q(i_v,i  ,j-2,k)))* &
rhy*   (q(i_v,i  ,j  ,k) - q(i_v,i  ,j-1,k)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_v,i  ,j+1,k) + q(i_v,i  ,j+2,k)) - c1r16*(q(i_v,i  ,j  ,k) + q(i_v,i  ,j+3,k)))* &
rh3y*  (q(i_v,i  ,j+3,k) - q(i_v,i  ,j  ,k)) &
+(c9r16* (q(i_v,i  ,j-1,k) + q(i_v,i  ,j-2,k)) - c1r16*(q(i_v,i  ,j  ,k) + q(i_v,i  ,j-3,k)))* &
rh3y*  (q(i_v,i  ,j  ,k) - q(i_v,i  ,j-3,k)) ) )

endif

enddo
enddo
enddo


!.... fadv(i_w) 


do k=kf_start, kf_end
do j=1,n2
do i=if_start, if_end

rjy = 1.d0 / jacb(i_u,j)

!.... fadv(i_w) = u dw/dx + ....
fadv(i_w,i,j,k) = c9r8 * c1r2* ( (c9r16* (q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)) - c1r16*(q(i_u,i  ,j,k-2) + q(i_u,i  ,j,k+1)))* &
rhx*   (q(i_w,i  ,j,k  ) - q(i_w,i-1,j,k  )) &
+(c9r16* (q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)) - c1r16*(q(i_u,i+1,j,k-2) + q(i_u,i+1,j,k+1)))* &
rhx*   (q(i_w,i+1,j,k  ) - q(i_w,i  ,j,k  )) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_u,i-1,j,k  ) + q(i_u,i-1,j,k-1)) - c1r16*(q(i_u,i-1,j,k-2) + q(i_u,i-1,j,k+1)))* &
rh3x*  (q(i_w,i  ,j,k  ) - q(i_w,i-3,j,k  )) &
+(c9r16* (q(i_u,i+2,j,k  ) + q(i_u,i+2,j,k-1)) - c1r16*(q(i_u,i+2,j,k-2) + q(i_u,i+2,j,k+1)))* &
rh3x*  (q(i_w,i+3,j,k  ) - q(i_w,i  ,j,k  )) )

!.... fadv(i_w) = ... + w dw/dz 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) + &
c9r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-1)) - c1r16*(q(i_w,i  ,j,k+1) + q(i_w,i  ,j,k-2)))* &
rhz*   (q(i_w,i  ,j,k  ) - q(i_w,i  ,j,k-1)) &
+(c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+1)) - c1r16*(q(i_w,i  ,j,k-1) + q(i_w,i  ,j,k+2)))* &
rhz*   (q(i_w,i  ,j,k+1) - q(i_w,i  ,j,k  )) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k-1) + q(i_w,i  ,j,k-2)) - c1r16*(q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-3)))* &
rh3z*  (q(i_w,i  ,j,k  ) - q(i_w,i  ,j,k-3)) &
+(c9r16* (q(i_w,i  ,j,k+1) + q(i_w,i  ,j,k+2)) - c1r16*(q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+3)))* &
rh3z*  (q(i_w,i  ,j,k+3) - q(i_w,i  ,j,k  )) )



! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----


if ( j .eq. 1 ) then

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q(i_w,i,j+1,k  ) - q(i_w,i,j  ,k  )) &
+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q(i_w,i,j  ,k  ) - q(i_w,i,j-1,k  )) ) )

elseif ( j .eq. 2 ) then

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q(i_w,i,j+1,k  ) - q(i_w,i,j  ,k  )) &
+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q(i_w,i,j  ,k  ) - q(i_w,i,j-1,k  )) ) )

elseif ( j .eq. n2-1 ) then

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q(i_w,i,j+1,k  ) - q(i_w,i,j  ,k  )) &
+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q(i_w,i,j  ,k  ) - q(i_w,i,j-1,k  )) ) )

elseif ( j .eq. n2 ) then

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q(i_w,i,j+1,k  ) - q(i_w,i,j  ,k  )) &
+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q(i_w,i,j  ,k  ) - q(i_w,i,j-1,k  )) ) )
else

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c9r8 * c1r2* ( (c9r16* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) - c1r16*(q(i_v,i,j  ,k+1) + q(i_v,i,j  ,k-2)))* &
rhy*   (q(i_w,i,j+1,k  ) - q(i_w,i,j  ,k  )) &
+(c9r16* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) - c1r16*(q(i_v,i,j-1,k+1) + q(i_v,i,j-1,k-2)))* &
rhy*   (q(i_w,i,j  ,k  ) - q(i_w,i,j-1,k  )) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_v,i,j-2,k  ) + q(i_v,i,j-2,k-1)) - c1r16*(q(i_v,i,j-2,k+1) + q(i_v,i,j-2,k-2)))* &
rh3y*  (q(i_w,i,j  ,k  ) - q(i_w,i,j-3,k  )) &
+(c9r16* (q(i_v,i,j+1,k  ) + q(i_v,i,j+1,k-1)) - c1r16*(q(i_v,i,j+1,k+1) + q(i_v,i,j+1,k-2)))* &
rh3y*  (q(i_w,i,j+3,k  ) - q(i_w,i,j  ,k  )) ) )

endif
enddo
enddo
enddo




   !------------
   ! skew = 0.5 div + 0.5 adv 
   !------------

   do k=kf_start, kf_end
    do j=1,n2
     do i=if_start, if_end
      do n=1,nvel

      !     fnl(n,i,j,k) = fdiv(n,i,j,k)                                 ! divergence form 
      !       fnl(n,i,j,k) = fadv(n,i,j,k)                                 ! advective  form 
      fnl(n,i,j,k) = 0.5d0 * fdiv(n,i,j,k) + 0.5d0 * fadv(n,i,j,k)  ! skew-sym   form 

      enddo
     enddo
    enddo
   enddo

    !.... clean-up temporary memory 
   deallocate(fdiv)
   deallocate(fadv)


  end subroutine cnv4th_2nd


   !===============================================================================================================================! 

  subroutine cnv4th ( fnl ) 

  !--------------------------------------
  ! compute convective terms using skew-symmetric form (cf. vasilyev, jcp, 2000) 
  ! to fourth order accuracy
  !--------------------------------------

   implicit none 
   real, dimension(1:nvel,is_df:ie_df, 1:n2,ks_df:ke_df), intent(inout)       :: fnl 
   real, dimension(:,:,:,:) , allocatable                         :: fdiv  , fadv
   integer                                                        :: i,j,k,n, if_start, if_end, kf_start, kf_end 
   real                                                           :: rhx, rhz, rh3x, rh3z, hy0, jac , rjy, rj3y
   real                                                           :: c9r8, c9r16, c1r8, c1r16, c1r2
   real                                                           :: c8r8  ! debug 
   real                                                           :: rhy, rh3y, dudx, dwdz

   !----------
   ! init halo indexing
   !---------

   if ( jxzf ) then 
   if_start = is_df
   if_end   = ie_df
   kf_start = ks_df
   kf_end   = ke_df
   else  
   if_start = is 
   if_end   = ie 
   kf_start = ks 
   kf_end   = ke 
   endif 

   !--------
   ! constants 
   !--------

   c9r8  = 9.0d0/8.0d0 
   c9r16 = 9.0d0/16.d0 
   c1r8  = 1.0d0/8.0d0 
   c1r16 = 1.0d0/16.d0 
   c1r2  = 1.0d0/2.0d0 
   c8r8  = 1.0d0


   rhx   = dble(n1) / (chlx*pi) 
   rhz   = dble(n3) / (chlz*pi)
   rh3x  = rhx/3.0d0 
   rh3z  = rhz/3.0d0 
   hy0   = 2.d0/ dble(n2) 
   rhy   = dble(n2) / 2.d0 
   rh3y  = rhy/3.d0


   !------------------
   ! define pressure gradient for 
   ! verification routines 
   !------------------
   if ( osverify)  pgm = 2.d0/re 
   if ( mmsverify) pgm = 0.d0 


   allocate ( fdiv(1:nvel, if_start:if_end , 1:n2  , kf_start:kf_end) )
   allocate ( fadv(1:nvel, if_start:if_end , 1:n2  , kf_start:kf_end) ) 

   fdiv = 0.d0 
   fadv = 0.d0 


!--------------------
! divergence form  
!--------------------


!... compute convective terms in u 


    do k=kf_start, kf_end
    do j=1,n2
    do i=if_start, if_end

       jac  = jacb(i_u,j)
       rjy  = 1.d0/hy0/jac
       rj3y = 1.d0/3.d0/hy0/jac

!.... fdiv(u) = d/dx(uu ) + ... 
       fdiv(i_u,i,j,k) = c9r8*rhx* ( (c9r16*( q(i_u,i  ,j,k) + q(i_u,i+1,j,k)) &
                                     -c1r16*( q(i_u,i-1,j,k) + q(i_u,i+2,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i+1,j,k))  &
                                   - (c9r16*( q(i_u,i  ,j,k) + q(i_u,i-1,j,k)) &
                                     -c1r16*( q(i_u,i+1,j,k) + q(i_u,i-2,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i-1,j,k)) )&
                        -c1r8*rh3x*( (c9r16*( q(i_u,i+1,j,k) + q(i_u,i+2,j,k)) &
                                     -c1r16*( q(i_u,i  ,j,k) + q(i_u,i+3,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i+3,j,k))  &
                                   - (c9r16*( q(i_u,i-1,j,k) + q(i_u,i-2,j,k)) &
                                     -c1r16*( q(i_u,i  ,j,k) + q(i_u,i-3,j,k)))* c1r2*(q(i_u,i  ,j,k)+q(i_u,i-3,j,k)) )


!.... fdiv(u) = + .... d/dy( uv) + ... 
       fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
                         c9r8*rjy* ( (c9r16*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)) &
                                     -c1r16*( q(i_v,i+1,j  ,k) + q(i_v,i-2,j  ,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i,j+1,k))  &
                                   - (c9r16*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)) &
                                     -c1r16*( q(i_v,i+1,j-1,k) + q(i_v,i-2,j-1,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i,j-1,k)) )&
                        -c1r8*rj3y*( (c9r16*( q(i_v,i  ,j+1,k) + q(i_v,i-1,j+1,k)) &
                                     -c1r16*( q(i_v,i+1,j+1,k) + q(i_v,i-2,j+1,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i,j+3,k))  &
                                   - (c9r16*( q(i_v,i  ,j-2,k) + q(i_v,i-1,j-2,k)) &
                                     -c1r16*( q(i_v,i+1,j-2,k) + q(i_v,i-2,j-2,k)))* c1r2*(q(i_u,i,j  ,k)+q(i_u,i,j-3,k)) )


!... fdiv(u) = ... + d/dz(uw )  
       fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
                         c9r8*rhz* ( (c9r16*( q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)) &
                                     -c1r16*( q(i_w,i+1,j,k+1) + q(i_w,i-2,j,k+1)))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k+1))  &
                                   - (c9r16*( q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )) &
                                     -c1r16*( q(i_w,i+1,j,k  ) + q(i_w,i-2,j,k  )))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k-1)) )&
                        -c1r8*rh3z*( (c9r16*( q(i_w,i  ,j,k+2) + q(i_w,i-1,j,k+2)) &
                                     -c1r16*( q(i_w,i+1,j,k+2) + q(i_w,i-2,j,k+2)))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k+3))  &
                                   - (c9r16*( q(i_w,i  ,j,k-1) + q(i_w,i-1,j,k-1)) &
                                     -c1r16*( q(i_w,i+1,j,k-1) + q(i_w,i-2,j,k-1)))* c1r2*(q(i_u,i,j,k  )+q(i_u,i,j,k-3)) )

    enddo
    enddo
    enddo




!.... compute convective terms in v  (wn ) 

    do k=kf_start, kf_end
    do j=1,n2-1
    do i=if_start, if_end

       jac  = jacb(i_v,j)
       rjy  = 1.d0/hy0/jac
       rj3y = 1.d0/3.d0/hy0/jac

!.... fdiv(v) = d/dx(vu) + ...  
       fdiv(i_v,i,j,k) = c9r8*rhx* ( (c9r16*( q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)) &
                                     -c1r16*( q(i_u,i+1,j-1,k) + q(i_u,i+1,j+2,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i+1,j,k))  &
                                   - (c9r16*( q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)) &
                                     -c1r16*( q(i_u,i  ,j-1,k) + q(i_u,i  ,j+2,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i-1,j,k)) )&
                        -c1r8*rh3x*( (c9r16*( q(i_u,i+2,j  ,k) + q(i_u,i+2,j+1,k)) &
                                     -c1r16*( q(i_u,i+2,j-1,k) + q(i_u,i+2,j+2,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i+3,j,k))  &
                                   - (c9r16*( q(i_u,i-1,j  ,k) + q(i_u,i-1,j+1,k)) &
                                     -c1r16*( q(i_u,i-1,j-1,k) + q(i_u,i-1,j+2,k)))* c1r2*(q(i_v,i  ,j,k)+q(i_v,i-3,j,k)) )



!.... fdiv(v) = ... + d/dy(vv) + .... 
       fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
                         c9r8*rjy* ( (c9r16*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)) &
                                     -c1r16*( q(i_v,i,j-1,k) + q(i_v,i,j+2,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j+1,k))  &
                                   - (c9r16*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)) &
                                     -c1r16*( q(i_v,i,j+1,k) + q(i_v,i,j-2,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j-1,k)) )&
                        -c1r8*rj3y*( (c9r16*( q(i_v,i,j+1,k) + q(i_v,i,j+2,k)) &
                                     -c1r16*( q(i_v,i,j  ,k) + q(i_v,i,j+3,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j+3,k))  &
                                   - (c9r16*( q(i_v,i,j-1,k) + q(i_v,i,j-2,k)) &
                                     -c1r16*( q(i_v,i,j  ,k) + q(i_v,i,j-3,k)))* c1r2*(q(i_v,i,j  ,k)+q(i_v,i,j-3,k)) )


!... fdiv(v) = ... + d/dz(vw) 
       fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
                         c9r8*rhz* ( (c9r16*( q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)) &
                                     -c1r16*( q(i_w,i,j-1,k+1) + q(i_w,i,j+2,k+1)))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k+1))  &
                                   - (c9r16*( q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )) &
                                     -c1r16*( q(i_w,i,j-1,k  ) + q(i_w,i,j+2,k  )))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k-1)) )&
                        -c1r8*rh3z*( (c9r16*( q(i_w,i,j  ,k+2) + q(i_w,i,j+1,k+2)) &
                                     -c1r16*( q(i_w,i,j-1,k+2) + q(i_w,i,j+2,k+2)))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k+3))  &
                                   - (c9r16*( q(i_w,i,j  ,k-1) + q(i_w,i,j+1,k-1)) &
                                     -c1r16*( q(i_w,i,j-1,k-1) + q(i_w,i,j+2,k-1)))* c1r2*(q(i_v,i,j,k  )+q(i_v,i,j,k-3)) )

    enddo
    enddo
    enddo


!.... compute convective terms in w (span) 

    do k=kf_start, kf_end
    do j=1,n2
    do i=if_start, if_end

       jac  = jacb(i_u,j)
       rjy  = 1.d0/hy0/jac
       rj3y = 1.d0/3.d0/hy0/jac

!.... fdiv(w) = d/dx(wu) + .... 
       fdiv(i_w,i,j,k) = c9r8*rhx* ( (c9r16*( q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)) &
                                     -c1r16*( q(i_u,i+1,j,k+1) + q(i_u,i+1,j,k-2)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i+1,j,k))  &
                                   - (c9r16*( q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)) &
                                     -c1r16*( q(i_u,i  ,j,k+1) + q(i_u,i  ,j,k-2)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i-1,j,k)) )&
                        -c1r8*rh3x*( (c9r16*( q(i_u,i+2,j,k  ) + q(i_u,i+2,j,k-1)) &
                                     -c1r16*( q(i_u,i+2,j,k+1) + q(i_u,i+2,j,k-2)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i+3,j,k))  &
                                   - (c9r16*( q(i_u,i-1,j,k  ) + q(i_u,i-1,j,k-1)) &
                                     -c1r16*( q(i_u,i-1,j,k+1) + q(i_u,i-1,j,k-2)))* c1r2*(q(i_w,i  ,j,k)+q(i_w,i-3,j,k)) )


!.... fdiv(w) = ... + d/dy(wv) + .... 
       fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
                         c9r8*rjy* ( (c9r16*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) &
                                     -c1r16*( q(i_v,i,j  ,k+1) + q(i_v,i,j  ,k-2)))* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j+1,k))  &
                                   - (c9r16*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) &
                                     -c1r16*( q(i_v,i,j-1,k+1) + q(i_v,i,j-1,k-2)))* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j-1,k)) )&
                        -c1r8*rj3y*( (c9r16*( q(i_v,i,j+1,k  ) + q(i_v,i,j+1,k-1)) &
                                     -c1r16*( q(i_v,i,j+1,k+1) + q(i_v,i,j+1,k-2)))* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j+3,k))  &
                                   - (c9r16*( q(i_v,i,j-2,k  ) + q(i_v,i,j-2,k-1)) &
                                     -c1r16*( q(i_v,i,j-2,k+1) + q(i_v,i,j-2,k-2)))* c1r2*(q(i_w,i,j  ,k)+q(i_w,i,j-3,k)) )


!.... fdiv(w) = ... + d/dz(ww)       
       fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
                         c9r8*rhz* ( (c9r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k+1)) &
                                     -c1r16*( q(i_w,i,j,k-1) + q(i_w,i,j,k+2)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k+1))  &
                                   - (c9r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k-1)) &
                                     -c1r16*( q(i_w,i,j,k+1) + q(i_w,i,j,k-2)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k-1))) &
                        -c1r8*rh3z*( (c9r16*( q(i_w,i,j,k+1) + q(i_w,i,j,k+2)) &
                                     -c1r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k+3)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k+3))  &
                                   - (c9r16*( q(i_w,i,j,k-1) + q(i_w,i,j,k-2)) &
                                     -c1r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k-3)))* c1r2*(q(i_w,i,j,k  )+q(i_w,i,j,k-3)))


    enddo
    enddo
    enddo

    !--------------
    ! advective form 
    !--------------

    do k=kf_start, kf_end
    do j=1,n2
    do i=if_start, if_end

       rjy = 1.d0 / jacb(i_u,j)

!.... f(i_u) = u du/dx + .... 
       fadv(i_u,i,j,k) = c9r8 * c1r2* ( (c9r16* (q(i_u,i  ,j,k) + q(i_u,i-1,j,k)) - c1r16*(q(i_u,i-2,j,k) + q(i_u,i+1,j,k)))* &
                                         rhx*   (q(i_u,i  ,j,k) - q(i_u,i-1,j,k)) &
                                       +(c9r16* (q(i_u,i  ,j,k) + q(i_u,i+1,j,k)) - c1r16*(q(i_u,i-1,j,k) + q(i_u,i+2,j,k)))* &
                                         rhx*   (q(i_u,i+1,j,k) - q(i_u,i  ,j,k)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_u,i-2,j,k) + q(i_u,i-1,j,k)) - c1r16*(q(i_u,i-3,j,k) + q(i_u,i  ,j,k)))* &
                                         rh3x*  (q(i_u,i  ,j,k) - q(i_u,i-3,j,k)) &
                                       +(c9r16* (q(i_u,i+2,j,k) + q(i_u,i+1,j,k)) - c1r16*(q(i_u,i+3,j,k) + q(i_u,i  ,j,k)))* &
                                         rh3x*  (q(i_u,i+3,j,k) - q(i_u,i  ,j,k)) )


!.... f(i_u) =  + ... v du/dy + ... 
       fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
                         c9r8 * c1r2*(  (c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)) - c1r16*(q(i_v,i-2,j  ,k) + q(i_v,i+1,j  ,k)))* &
                                         rhy*   (q(i_u,i  ,j+1,k) - q(i_u,i  ,j  ,k)) &
                                       +(c9r16* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)) - c1r16*(q(i_v,i-2,j-1,k) + q(i_v,i+1,j-1,k)))* &
                                         rhy*   (q(i_u,i  ,j  ,k) - q(i_u,i  ,j-1,k)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_v,i  ,j-2,k) + q(i_v,i-1,j-2,k)) - c1r16*(q(i_v,i-2,j-2,k) + q(i_v,i+1,j-2,k)))* &
                                         rh3y*  (q(i_u,i  ,j  ,k) - q(i_u,i  ,j-3,k)) &
                                       +(c9r16* (q(i_v,i  ,j+1,k) + q(i_v,i-1,j+1,k)) - c1r16*(q(i_v,i-2,j+1,k) + q(i_v,i+1,j+1,k)))* &
                                         rh3y*  (q(i_u,i  ,j+3,k) - q(i_u,i  ,j  ,k)) ) )

!.... f(i_u) =  ... + w du/dz 
       fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + &
                         c9r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)) - c1r16*(q(i_w,i-2,j,k+1) + q(i_w,i+1,j,k+1)))* &
                                         rhz*   (q(i_u,i  ,j,k+1) - q(i_u,i  ,j,k  )) &
                                       +(c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )) - c1r16*(q(i_w,i-2,j,k  ) + q(i_w,i+1,j,k  )))* &
                                         rhz*   (q(i_u,i  ,j,k  ) - q(i_u,i  ,j,k-1)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k-1) + q(i_w,i-1,j,k-1)) - c1r16*(q(i_w,i-2,j,k-1) + q(i_w,i+1,j,k-1)))* &
                                         rh3z*  (q(i_u,i  ,j,k  ) - q(i_u,i  ,j,k-3)) &
                                       +(c9r16* (q(i_w,i  ,j,k+2) + q(i_w,i-1,j,k+2)) - c1r16*(q(i_w,i-2,j,k+2) + q(i_w,i+1,j,k+2)))* &
                                         rh3z*  (q(i_u,i  ,j,k+3) - q(i_u,i  ,j,k  )) )

       !---- 
       ! deprecated;; pressure gradient is now 
       ! added elsewhere 
       !---

       !fadv(i_u,i,j,k) = fadv(i_u,i,j,k) - pgm 



    enddo
    enddo
    enddo


!.... fadv(i_v) 


    do k=kf_start, kf_end
    do j=1,n2-1
    do i=if_start, if_end

       rjy = 1.d0 / jacb(i_v,j)


!.... f(i_v) =  u dv/dx 
       fadv(i_v,i,j,k) = c9r8 * c1r2* ( (c9r16* (q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)) - c1r16*(q(i_u,i  ,j-1,k) + q(i_u,i  ,j+2,k)))* &
                                         rhx*   (q(i_v,i  ,j  ,k) - q(i_v,i-1,j  ,k)) &
                                       +(c9r16* (q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)) - c1r16*(q(i_u,i+1,j-1,k) + q(i_u,i+1,j+2,k)))* &
                                         rhx*   (q(i_v,i+1,j  ,k) - q(i_v,i  ,j  ,k)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_u,i-1,j  ,k) + q(i_u,i-1,j+1,k)) - c1r16*(q(i_u,i-1,j-1,k) + q(i_u,i-1,j+2,k)))* &
                                         rh3x*  (q(i_v,i  ,j  ,k) - q(i_v,i-3,j  ,k)) &
                                       +(c9r16* (q(i_u,i+2,j  ,k) + q(i_u,i+2,j+1,k)) - c1r16*(q(i_u,i+2,j-1,k) + q(i_u,i+2,j+2,k)))* &
                                         rh3x*  (q(i_v,i+3,j  ,k) - q(i_v,i  ,j  ,k)) )

!.... f(i_v) = ... + v dv/dy + .... 
       fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
                         c9r8 * c1r2* ( (c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)) - c1r16*(q(i_v,i  ,j-1,k) + q(i_v,i  ,j+2,k)))* &
                                         rhy*   (q(i_v,i  ,j+1,k) - q(i_v,i  ,j  ,k)) &
                                       +(c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)) - c1r16*(q(i_v,i  ,j+1,k) + q(i_v,i  ,j-2,k)))* &
                                         rhy*   (q(i_v,i  ,j  ,k) - q(i_v,i  ,j-1,k)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_v,i  ,j+1,k) + q(i_v,i  ,j+2,k)) - c1r16*(q(i_v,i  ,j  ,k) + q(i_v,i  ,j+3,k)))* &
                                         rh3y*  (q(i_v,i  ,j+3,k) - q(i_v,i  ,j  ,k)) &
                                       +(c9r16* (q(i_v,i  ,j-1,k) + q(i_v,i  ,j-2,k)) - c1r16*(q(i_v,i  ,j  ,k) + q(i_v,i  ,j-3,k)))* &
                                         rh3y*  (q(i_v,i  ,j  ,k) - q(i_v,i  ,j-3,k)) ) )

!.... f(i_v) = ... + w dv/dz 
       fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + &
                         c9r8 * c1r2* ( (c9r16* (q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )) - c1r16*(q(i_w,i,j-1,k  ) + q(i_w,i,j+2,k  )))* &
                                         rhz*   (q(i_v,i,j  ,k  ) - q(i_v,i,j  ,k-1)) &
                                       +(c9r16* (q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)) - c1r16*(q(i_w,i,j-1,k+1) + q(i_w,i,j+2,k+1)))* &
                                         rhz*   (q(i_v,i,j  ,k+1) - q(i_v,i,j  ,k  )) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_w,i,j  ,k-1) + q(i_w,i,j+1,k-1)) - c1r16*(q(i_w,i,j-1,k-1) + q(i_w,i,j+2,k-1)))* &
                                         rh3z*  (q(i_v,i,j  ,k  ) - q(i_v,i,j  ,k-3)) &
                                       +(c9r16* (q(i_w,i,j  ,k+2) + q(i_w,i,j+1,k+2)) - c1r16*(q(i_w,i,j-1,k+2) + q(i_w,i,j+2,k+2)))* &
                                         rh3z*  (q(i_v,i,j  ,k+3) - q(i_v,i,j  ,k  )) )


    enddo
    enddo
    enddo

!.... fadv(i_w) 


    do k=kf_start, kf_end
    do j=1,n2
    do i=if_start, if_end

       rjy = 1.d0 / jacb(i_u,j)



!.... fadv(i_w) = u dw/dx + ....
       fadv(i_w,i,j,k) = c9r8 * c1r2* ( (c9r16* (q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)) - c1r16*(q(i_u,i  ,j,k-2) + q(i_u,i  ,j,k+1)))* &
                                         rhx*   (q(i_w,i  ,j,k  ) - q(i_w,i-1,j,k  )) &
                                       +(c9r16* (q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)) - c1r16*(q(i_u,i+1,j,k-2) + q(i_u,i+1,j,k+1)))* &
                                         rhx*   (q(i_w,i+1,j,k  ) - q(i_w,i  ,j,k  )) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_u,i-1,j,k  ) + q(i_u,i-1,j,k-1)) - c1r16*(q(i_u,i-1,j,k-2) + q(i_u,i-1,j,k+1)))* &
                                         rh3x*  (q(i_w,i  ,j,k  ) - q(i_w,i-3,j,k  )) &
                                       +(c9r16* (q(i_u,i+2,j,k  ) + q(i_u,i+2,j,k-1)) - c1r16*(q(i_u,i+2,j,k-2) + q(i_u,i+2,j,k+1)))* &
                                         rh3x*  (q(i_w,i+3,j,k  ) - q(i_w,i  ,j,k  )) )


!.... fadv(i_w) = ...+ v dw/dy + ... 
       fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
                         c9r8 * c1r2* ( (c9r16* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) - c1r16*(q(i_v,i,j  ,k+1) + q(i_v,i,j  ,k-2)))* &
                                         rhy*   (q(i_w,i,j+1,k  ) - q(i_w,i,j  ,k  )) &
                                       +(c9r16* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) - c1r16*(q(i_v,i,j-1,k+1) + q(i_v,i,j-1,k-2)))* &
                                         rhy*   (q(i_w,i,j  ,k  ) - q(i_w,i,j-1,k  )) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_v,i,j-2,k  ) + q(i_v,i,j-2,k-1)) - c1r16*(q(i_v,i,j-2,k+1) + q(i_v,i,j-2,k-2)))* &
                                         rh3y*  (q(i_w,i,j  ,k  ) - q(i_w,i,j-3,k  )) &
                                       +(c9r16* (q(i_v,i,j+1,k  ) + q(i_v,i,j+1,k-1)) - c1r16*(q(i_v,i,j+1,k+1) + q(i_v,i,j+1,k-2)))* &
                                         rh3y*  (q(i_w,i,j+3,k  ) - q(i_w,i,j  ,k  )) ) )

!.... fadv(i_w) = ... + w dw/dz 
       fadv(i_w,i,j,k) = fadv(i_w,i,j,k) + &
                         c9r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-1)) - c1r16*(q(i_w,i  ,j,k+1) + q(i_w,i  ,j,k-2)))* &
                                         rhz*   (q(i_w,i  ,j,k  ) - q(i_w,i  ,j,k-1)) &
                                       +(c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+1)) - c1r16*(q(i_w,i  ,j,k-1) + q(i_w,i  ,j,k+2)))* &
                                         rhz*   (q(i_w,i  ,j,k+1) - q(i_w,i  ,j,k  )) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k-1) + q(i_w,i  ,j,k-2)) - c1r16*(q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-3)))* &
                                         rh3z*  (q(i_w,i  ,j,k  ) - q(i_w,i  ,j,k-3)) &
                                       +(c9r16* (q(i_w,i  ,j,k+1) + q(i_w,i  ,j,k+2)) - c1r16*(q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+3)))* &
                                         rh3z*  (q(i_w,i  ,j,k+3) - q(i_w,i  ,j,k  )) )

    enddo
    enddo
    enddo


    !------------
    ! skew = 0.5 div + 0.5 adv 
    !------------

    do k=kf_start, kf_end
    do j=1,n2
    do i=if_start, if_end
    do n=1,nvel

!       fnl(n,i,j,k) = fdiv(n,i,j,k)                                 ! divergence form 
!       fnl(n,i,j,k) = fadv(n,i,j,k)                                 ! advective  form 
       fnl(n,i,j,k) = 0.5d0 * fdiv(n,i,j,k) + 0.5d0 * fadv(n,i,j,k)  ! skew-sym   form 

    enddo
    enddo
    enddo
    enddo


    !.... clean-up temporary memory 
   deallocate(fdiv) 
   deallocate(fadv) 


  end subroutine cnv4th

!=============================================================================================================================!

  subroutine gradp2nd ( gphi, phi )

    !----------------
    ! subroutine computes the gradient of the scalar phi 
    !----------------

    implicit none
    real, dimension(1:nvel, is_sf:ie_sf, 1:n2, ks_sf:ke_sf) , intent(out )                                     :: gphi
    real, dimension(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp, ks-halo_z+1:ke+halo_z+1), intent(in)         :: phi
    integer                                                        :: i,j,k
    real                                                           :: rhx, rhz, rh3x, rh3z, hy0, c9r8, c1r8, jac, rjy, rj3y
    integer                                                        :: kstart, kend, istart, iend

    !---- 
    ! define indexing bounds 
    !---- 
    istart = is_sf
    iend   = ie_sf
    kstart = ks_sf
    kend   = ke_sf

    gphi = 0.d0

    !---- 
    ! define grid parameters, constants
    !---- 
    rhx  = dble(n1) / ( chlx*pi )
    rhz  = dble(n3) / ( chlz*pi )
    rh3x = rhx / 3.d0
    rh3z = rhz / 3.d0
    hy0  = 2.d0/ dble(n2)
    c9r8 = 9.d0/ 8.d0
    c1r8 = 1.d0/ 8.d0

    !---
    ! i_u 
    !---

    do k=kstart+1,kend+1
    do j=1,n2
    do i=istart,iend

       gphi(i_u,i,j,k-1) =  rhx * (phi(i+1,j,k) - phi(i  ,j,k))
    enddo
    enddo
    enddo


    !---
    ! i_v 
    !---


    do k=kstart+1,kend+1
    do j=1,n2-1
    do i=istart+1,iend+1

       jac  = jacb(i_v,j)
       rjy  = 1.d0/hy0/jac
       rj3y = 1.d0/3.0d0/hy0/jac

       gphi(i_v,i-1,j,k-1) =  rjy * (phi(i,j+1,k) - phi(i,j  ,k)) 

    enddo
    enddo
    enddo

    !---
    ! i_w 
    !---


    do k=kstart,kend
    do j=1,n2
    do i=istart+1,iend+1

       gphi(i_w,i-1,j,k) = rhz * (phi(i,j,k+1) - phi(i,j,k  ))
    enddo
    enddo
    enddo


  end subroutine gradp2nd

!=====================================================================================================================================!

  subroutine gradp ( gphi, phi ) 

    !----------------
    ! subroutine computes the gradient of the scalar phi 
    !----------------

    implicit none 
    real, dimension(1:nvel, is_sf:ie_sf, 1:n2, ks_sf:ke_sf) , intent(out )                                     :: gphi 
    real, dimension(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp, ks-halo_z+1:ke+halo_z+1), intent(in)         :: phi 
    integer                                                        :: i,j,k
    real                                                           :: rhx, rhz, rh3x, rh3z, hy0, c9r8, c1r8, jac, rjy, rj3y
    integer                                                        :: kstart, kend, istart, iend

    !---- 
    ! define indexing bounds 
    !---- 
    istart = is_sf 
    iend   = ie_sf 
    kstart = ks_sf 
    kend   = ke_sf 

    gphi = 0.d0 

    !---- 
    ! define grid parameters, constants
    !---- 
    rhx  = dble(n1) / ( chlx*pi )
    rhz  = dble(n3) / ( chlz*pi )
    rh3x = rhx / 3.d0
    rh3z = rhz / 3.d0
    hy0  = 2.d0/ dble(n2)
    c9r8 = 9.d0/ 8.d0
    c1r8 = 1.d0/ 8.d0


    !---
    ! i_u 
    !---

    do k=kstart+1,kend+1
    do j=1,n2 
    do i=istart,iend

       gphi(i_u,i,j,k-1) =  c9r8*rhx * (phi(i+1,j,k) - phi(i  ,j,k)) - & 
                            c1r8*rh3x* (phi(i+2,j,k) - phi(i-1,j,k)) 
    enddo
    enddo 
    enddo 

    !---
    ! i_v 
    !---


    do k=kstart+1,kend+1
    do j=1,n2-1 
    do i=istart+1,iend+1

       jac  = jacb(i_v,j) 
       rjy  = 1.d0/hy0/jac 
       rj3y = 1.d0/3.0d0/hy0/jac 

       gphi(i_v,i-1,j,k-1) =  c9r8*rjy * (phi(i,j+1,k) - phi(i,j  ,k)) - & 
                              c1r8*rj3y* (phi(i,j+2,k) - phi(i,j-1,k)) 

    enddo
    enddo 
    enddo 

    !---
    ! i_w 
    !---


    do k=kstart,kend
    do j=1,n2 
    do i=istart+1,iend+1

       gphi(i_w,i-1,j,k) = c9r8*rhz * (phi(i,j,k+1) - phi(i,j,k  )) - & 
                           c1r8*rh3z* (phi(i,j,k+2) - phi(i,j,k-1)) 
    enddo 
    enddo 
    enddo 

 
  end subroutine gradp

!======================================================================================================================================!

  subroutine wn_op ( wn_u, ipiv_u, wn_v, ipiv_v, wn_w, ipiv_w, beta , ii, kk, wvisc, vtb) 

    !------------------------------------
    ! compute lu factorization of the banded matrix I - dt*beta*(lub)
    ! and lvb 
    ! returns lu factorization and permutation vector 
    !------------------------------------ 

    ! be careful to not call this function 
    ! anymore with wvisc = eddy_visc 
    ! it is unclear what results you'll get 
    
    
    implicit none 
    integer, intent(in)                     :: wvisc 
    real,     intent(in)                    :: beta
    real,    intent(in)                     :: vtb      ! mean vt to be added 
    real, dimension(1:2*kl+ku+1,1:n2)       :: wn_u 
    real, dimension(1:2*kl+ku+1,1:n2)       :: wn_w 
    real, dimension(1:2*kl+ku+1,1:n2-1)     :: wn_v 
    integer, dimension(1:n2 )               :: ipiv_u 
    integer, dimension(1:n2 )               :: ipiv_w 
    integer, dimension(1:n2-1)              :: ipiv_v 
    integer                                 :: i,j,k,n, info 
    real                                    :: rre_u, rre_v, rre_w 
    integer, intent(in)                     :: ii,kk                !! position 
    real                                    :: c9r16, c1r16
    
    
    rre_u = 1.d0 / re   + vtb                       ! default for non-sgs model  
    rre_v = 1.d0 / re   + vtb 
    rre_w = 1.d0 / re   + vtb 
    
    c9r16 = 9.d0 / 16.d0 
    c1r16 = 1.d0 / 16.d0 

    !----
    ! wn_op : u,w
    !----

    do j=1,n2 
    do i= max(1,j-ku)  , min(n2, j+kl )  
    
       if ( wvisc .eq. eddy_visc ) then 
          
          rre_u = c9r16* ( vt(ii-1,i,kk) + vt(ii  ,i,kk)) - c1r16* ( vt(ii-2,i,kk) + vt(ii+1,i,kk)) + rre_u 
          rre_w = c9r16* ( vt(ii,i,kk-1) + vt(ii,i,kk  )) - c1r16* ( vt(ii,i,kk-2) + vt(ii,i,kk+1)) + rre_w
       endif

       if (  i .eq. j) then 
          wn_u(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_u* lub(ku+kl+1+i-j,j) 
          wn_w(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_w* lwb(ku+kl+1+i-j,j) 
       else 
          wn_u(ku+kl+1+i-j, j) = -dt*beta* rre_u* lub(ku+kl+1+i-j,j) 
          wn_w(ku+kl+1+i-j, j) = -dt*beta* rre_w* lwb(ku+kl+1+i-j,j)
       endif
    enddo
    enddo

    !---- 
    ! operator from sgs gradient term 
    !----
    call dgbtrf ( n2, n2, kl, ku, wn_u, 2*kl+ku+1, ipiv_u, info ) 
    call dgbtrf ( n2, n2, kl, ku, wn_w, 2*kl+ku+1, ipiv_w, info ) 
    

    !----
    ! wn_op : v 
    !----

    do j=1,n2-1 
    do i= max(1,j-ku) , min(n2-1, j+kl) 

       if ( wvisc .eq. eddy_visc ) then 
          rre_v = c9r16* ( vt(ii,i+1,kk) + vt(ii,i  ,kk)) - c1r16* ( vt(ii,i-1,kk) + vt(ii,i+2,kk)) + nu_molec 
       endif

       if ( i .eq. j) then 
          wn_v(ku+kl+1+i-j,j) = 1.d0 - dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       else 
          wn_v(ku+kl+1+i-j,j) = -dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       endif
    enddo
    enddo

    !----
    ! operator from sgs gradient term 
    !----

    if ( wvisc .eq. eddy_visc ) then 
       
       do j=1,n2-1 
       do i=max(1,j-ku), min(n2-1,j+kl) 

          wn_v(ku+kl+1+i-j,j) = wn_v(ku+kl+1+i-j,j) - fsgs*dt*beta* evb_v(ku+kl+1+i-j,j) 
       enddo 
       enddo 
    endif


    call dgbtrf ( n2-1,n2-1, kl,ku, wn_v, 2*kl+ku+1, ipiv_v, info) 



  end subroutine wn_op


!======================================================================================================================================!

  subroutine wn_op_phase ( wn_u_phase, ipiv_u_phase, wn_v_phase, ipiv_v_phase, wn_w_phase, ipiv_w_phase, beta , ii, kk, wvisc, vtb) 

    !------------------------------------
    ! compute lu factorization of the banded matrix I - dt*beta*(lub)
    ! and lvb 
    ! returns lu factorization and permutation vector 
    !------------------------------------ 

    ! be careful to not call this function 
    ! anymore with wvisc = eddy_visc 
    ! it is unclear what results you'll get 
    
    
    implicit none 
    integer, intent(in)                     :: wvisc 
    real,     intent(in)                    :: beta
    real,    intent(in)                     :: vtb      ! mean vt to be added 
    real, dimension(1:2*kl+ku+1,1:n2)       :: wn_u_phase
    real, dimension(1:2*kl+ku+1,1:n2)       :: wn_w_phase 
    real, dimension(1:2*kl+ku+1,1:n2-1)     :: wn_v_phase 
    integer, dimension(1:n2 )               :: ipiv_u_phase 
    integer, dimension(1:n2 )               :: ipiv_w_phase 
    integer, dimension(1:n2-1)              :: ipiv_v_phase 
    integer                                 :: i,j,k,n, info 
    real                                    :: rre_u, rre_v, rre_w 
    integer, intent(in)                     :: ii,kk                !! position 
    real                                    :: c9r16, c1r16
    
    
    rre_u = 1.d0 / re   + vtb                       ! default for non-sgs model  
    rre_v = 1.d0 / re   + vtb 
    rre_w = 1.d0 / re   + vtb 
    
    c9r16 = 9.d0 / 16.d0 
    c1r16 = 1.d0 / 16.d0 

    !----
    ! wn_op : u,w
    !----

    do j=1,n2 
    do i= max(1,j-ku)  , min(n2, j+kl )  
    
       if ( wvisc .eq. eddy_visc ) then 
          
          rre_u = c9r16* ( vt(ii-1,i,kk) + vt(ii  ,i,kk)) - c1r16* ( vt(ii-2,i,kk) + vt(ii+1,i,kk)) + rre_u 
          rre_w = c9r16* ( vt(ii,i,kk-1) + vt(ii,i,kk  )) - c1r16* ( vt(ii,i,kk-2) + vt(ii,i,kk+1)) + rre_w
       endif

       if (  i .eq. j) then 
          wn_u_phase(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_u* lub_pat(ku+kl+1+i-j,j) 
          wn_w_phase(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_w* lub_pat(ku+kl+1+i-j,j) 
       else 
          wn_u_phase(ku+kl+1+i-j, j) = -dt*beta* rre_u* lub_pat(ku+kl+1+i-j,j) 
          wn_w_phase(ku+kl+1+i-j, j) = -dt*beta* rre_w* lub_pat(ku+kl+1+i-j,j)
       endif
    enddo
    enddo

    !---- 
    ! operator from sgs gradient term 
    !----
    call dgbtrf ( n2, n2, kl, ku, wn_u_phase, 2*kl+ku+1, ipiv_u_phase, info ) 
    call dgbtrf ( n2, n2, kl, ku, wn_w_phase, 2*kl+ku+1, ipiv_w_phase, info ) 
    

    !----
    ! wn_op : v 
    !----

    do j=1,n2-1 
    do i= max(1,j-ku) , min(n2-1, j+kl) 

       if ( wvisc .eq. eddy_visc ) then 
          rre_v = c9r16* ( vt(ii,i+1,kk) + vt(ii,i  ,kk)) - c1r16* ( vt(ii,i-1,kk) + vt(ii,i+2,kk)) + nu_molec 
       endif

       if ( i .eq. j) then 
          wn_v_phase(ku+kl+1+i-j,j) = 1.d0 - dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       else 
          wn_v_phase(ku+kl+1+i-j,j) = -dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       endif
    enddo
    enddo

    !----
    ! operator from sgs gradient term 
    !----

    if ( wvisc .eq. eddy_visc ) then 
       
       do j=1,n2-1 
       do i=max(1,j-ku), min(n2-1,j+kl) 

          wn_v_phase(ku+kl+1+i-j,j) = wn_v_phase(ku+kl+1+i-j,j) - fsgs*dt*beta* evb_v(ku+kl+1+i-j,j) 
       enddo 
       enddo 
    endif


    call dgbtrf ( n2-1,n2-1, kl,ku, wn_v_phase, 2*kl+ku+1, ipiv_v_phase, info) 



  end subroutine wn_op_phase

!=================================================================================================================================!

  subroutine wn_op_origin ( wn_u, wn_v, wn_w, beta , ii, kk, wvisc, vtb) 

    !------------------------------------
    ! compute lu factorization of the banded matrix I - dt*beta*(lub)
    ! and lvb 
    ! returns lu factorization and permutation vector 
    !------------------------------------ 

    ! be careful to not call this function 
    ! anymore with wvisc = eddy_visc 
    ! it is unclear what results you'll get 
    
    
    implicit none 
    integer, intent(in)                     :: wvisc 
    real,     intent(in)                    :: beta
    real,    intent(in)                     :: vtb      ! mean vt to be added 
    real, dimension(1:2*kl+ku+1,1:n2)       :: wn_u 
    real, dimension(1:2*kl+ku+1,1:n2)       :: wn_w 
    real, dimension(1:2*kl+ku+1,1:n2-1)     :: wn_v 
    integer                                 :: i,j,k,n, info 
    real                                    :: rre_u, rre_v, rre_w 
    integer, intent(in)                     :: ii,kk                !! position 
    real                                    :: c9r16, c1r16
    
    
    rre_u = 1.d0 / re   + vtb                       ! default for non-sgs model  
    rre_v = 1.d0 / re   + vtb 
    rre_w = 1.d0 / re   + vtb 
    
    c9r16 = 9.d0 / 16.d0 
    c1r16 = 1.d0 / 16.d0 

    !----
    ! wn_op : u,w
    !----

    do j=1,n2 
    do i= max(1,j-ku)  , min(n2, j+kl )  
    
       if ( wvisc .eq. eddy_visc ) then 
          
          rre_u = c9r16* ( vt(ii-1,i,kk) + vt(ii  ,i,kk)) - c1r16* ( vt(ii-2,i,kk) + vt(ii+1,i,kk)) + rre_u 
          rre_w = c9r16* ( vt(ii,i,kk-1) + vt(ii,i,kk  )) - c1r16* ( vt(ii,i,kk-2) + vt(ii,i,kk+1)) + rre_w
       endif

       if (  i .eq. j) then 
          wn_u(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_u* lub(ku+kl+1+i-j,j) 
          wn_w(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_w* lwb(ku+kl+1+i-j,j) 
       else 
          wn_u(ku+kl+1+i-j, j) = -dt*beta* rre_u* lub(ku+kl+1+i-j,j) 
          wn_w(ku+kl+1+i-j, j) = -dt*beta* rre_w* lwb(ku+kl+1+i-j,j)
       endif
    enddo
    enddo


    !----
    ! wn_op : v 
    !----

    do j=1,n2-1 
    do i= max(1,j-ku) , min(n2-1, j+kl) 

       if ( wvisc .eq. eddy_visc ) then 
          rre_v = c9r16* ( vt(ii,i+1,kk) + vt(ii,i  ,kk)) - c1r16* ( vt(ii,i-1,kk) + vt(ii,i+2,kk)) + nu_molec 
       endif

       if ( i .eq. j) then 
          wn_v(ku+kl+1+i-j,j) = 1.d0 - dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       else 
          wn_v(ku+kl+1+i-j,j) = -dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       endif
    enddo
    enddo

    !----
    ! operator from sgs gradient term 
    !----

    if ( wvisc .eq. eddy_visc ) then 
       
       do j=1,n2-1 
       do i=max(1,j-ku), min(n2-1,j+kl) 

          wn_v(ku+kl+1+i-j,j) = wn_v(ku+kl+1+i-j,j) - fsgs*dt*beta* evb_v(ku+kl+1+i-j,j) 
       enddo 
       enddo 
    endif


  end subroutine wn_op_origin

!======================================================================================================================================!

  subroutine wn_op_origin_pat ( pli_u_tee, ipiv_u_tee, wn_v_tee, ipiv_v_tee, pli_w_tee, ipiv_w_tee, beta , vtb) 

    !------------------------------------
    ! compute lu factorization of the banded matrix I - dt*beta*(lub)
    ! and lvb 
    ! returns lu factorization and permutation vector 
    !------------------------------------ 

    ! be careful to not call this function 
    ! anymore with wvisc = eddy_visc 
    ! it is unclear what results you'll get 
    
    
    implicit none 
    real,     intent(in)                    :: beta
    real,    intent(in)                     :: vtb      ! mean vt to be added 
    real, dimension(1:n2  ,1:n2  )          :: pli_u_tee, pli_w_tee
    real, dimension(1:2*kl+ku+1,1:n2-1)     :: wn_v_tee
    real, dimension(1:n2  ,1:n2  )          :: inv_solid, inv_phase
    real, dimension(1:2*kl+ku+1,1:n2)       :: wn_solid, wn_phase
    integer, dimension(1:n2 )               :: ipiv_solid, ipiv_phase
    integer, dimension(1:n2 )               :: ipiv_u_tee, ipiv_w_tee
    integer, dimension(1:n2-1)              :: ipiv_v_tee
    integer                                 :: i,j,k,n, info 
    real                                    :: rre_u, rre_v, rre_w 
    real                                    :: c9r16, c1r16
    
    
    rre_u = 1.d0 / re   + vtb                       ! default for non-sgs model  
    rre_v = 1.d0 / re   + vtb 
    rre_w = 1.d0 / re   + vtb 
    
    c9r16 = 9.d0 / 16.d0 
    c1r16 = 1.d0 / 16.d0 

    !----
    ! wn_op : u,w
    !----

    do j=1,n2 
    do i= max(1,j-ku)  , min(n2, j+kl )  
     
       if (  i .eq. j) then 
          wn_solid(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_u* lub(ku+kl+1+i-j,j) 
          wn_phase(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_u* lub_pat(ku+kl+1+i-j,j) 
       else 
          wn_solid(ku+kl+1+i-j, j) = -dt*beta* rre_u* lub(ku+kl+1+i-j,j) 
          wn_phase(ku+kl+1+i-j, j) = -dt*beta* rre_u* lub_pat(ku+kl+1+i-j,j)  
       endif
    enddo
    enddo

    call dgbtrf ( n2, n2, kl, ku, wn_solid, 2*kl+ku+1, ipiv_solid, info ) 
    call dgbtrf ( n2, n2, kl, ku, wn_phase, 2*kl+ku+1, ipiv_phase, info ) 
        

    !----
    ! invert wn_op
    !----

    ! create identity matrices
    do i=1,n2
    do j=1,n2
      if ( i .eq. j ) then
        inv_solid(i,j) = 1.d0
        inv_phase(i,j) = 1.d0
      else
        inv_solid(i,j) = 0.d0
        inv_phase(i,j) = 0.d0
      endif
    enddo
    enddo

    ! solve system for inverse
    call dgbtrs( 'n', n2  , kl, ku, n2  , wn_solid, 2*kl+ku+1, ipiv_solid, inv_solid, n2  , info )
    call dgbtrs( 'n', n2  , kl, ku, n2  , wn_phase, 2*kl+ku+1, ipiv_phase, inv_phase, n2  , info ) 

    !----
    ! solid fraction weighted lub
    !----

    pli_u_tee = nsolid_u/n1/n3 * inv_solid + (1 - nsolid_u/n1/n3) * inv_phase
    pli_w_tee = nsolid_w/n1/n3 * inv_solid + (1 - nsolid_w/n1/n3) * inv_phase

    call dgetrf ( n2, n2, pli_u_tee, n2, ipiv_u_tee, info )
    call dgetrf ( n2, n2, pli_w_tee, n2, ipiv_w_tee, info )

    !----
    ! wn_op : v 
    !----

    do j=1,n2-1 
    do i= max(1,j-ku) , min(n2-1, j+kl) 

       if ( i .eq. j) then 
          wn_v_tee(ku+kl+1+i-j,j) = 1.d0 - dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       else 
          wn_v_tee(ku+kl+1+i-j,j) = -dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       endif
    enddo
    enddo
 

  end subroutine wn_op_origin_pat

!=================================================================================================================================!

  subroutine p_op ( pwn , ipiv_p, kxs, kzs, mflag)

    !------------------------
    ! compute lu factorization of d^2 - (kx^2 + kz^2)*I 
    ! returns permutation vector along with result 
    !------------------------

    implicit none 
    real, intent(in)                       :: kxs 
    real, intent(in)                       :: kzs 
    real, dimension(2*klp+kup+1,1:n2)      :: pwn 
    integer, dimension(1:n2)               :: ipiv_p 
    integer                                :: i,j,k,info
    logical, optional                      :: mflag 

    pwn    = 0.d0 
    ipiv_p = 0

    do j=1,n2
    do i = max(1, j-kup), min(n2, j+klp)

       if ( i .eq. j ) then
          if ( kxs .eq. 0.d0 .and. kzs .eq. 0.d0 .and. i .eq. 1) then
             pwn(kup+klp+1+i-j,j) = 1.d0                ! pin pressure 
          else
             pwn(kup+klp+1+i-j,j) = pwb(kup+klp+1+i-j,j) - (kxs + kzs )
          endif
       else
          if ( kxs .ne. 0.d0 .or. kzs .ne. 0.d0 .or. i .ne. 1 ) pwn(kup+klp+1+i-j,j) = pwb(kup+klp+1+i-j,j)
       endif
    enddo
    enddo

    call dgbtrf ( n2, n2, klp,kup, pwn,  2*klp+kup+1, ipiv_p, info )

    if ( info .ne. 0 ) write(*,*) 'error in computing p_op ...', info 


  end subroutine p_op 

!==========================================================================================================================!

  subroutine p_op_2nd_Dirichlet ( pwn_2nd , ipiv_p, kxs, kzs, mflag)

    !------------------------
    ! compute lu factorization of d^2 - (kx^2 + kz^2)*I 
    ! returns permutation vector along with result 
    !------------------------

    implicit none
    real, intent(in)                       :: kxs
    real, intent(in)                       :: kzs
    real, dimension(2*klpp+kupp+1,1:n2)      :: pwn_2nd
    integer, dimension(1:n2)               :: ipiv_p
    integer                                :: i,j,k,info
    logical, optional                      :: mflag
    real :: hy0

    hy0   = 2.d0 / dble(n2 )
    pwn_2nd    = 0.d0
    ipiv_p = 0

    do j=1,n2
    do i = max(1, j-kupp), min(n2, j+klpp)

       if ( i .eq. j ) then
          if ( j .eq. 1 .or. j .eq. n2 ) then
             pwn_2nd(kupp+klpp+1+i-j,j) =   1.d0                ! pin pressure 
          else
             pwn_2nd(kupp+klpp+1+i-j,j) = pwb_2nd_Dirichlet(kupp+klpp+1+i-j,j) - (kxs + kzs )
          endif
       else
             pwn_2nd(kupp+klpp+1+i-j,j) = pwb_2nd_Dirichlet(kupp+klpp+1+i-j,j)
       endif
    enddo
    enddo

    call dgbtrf ( n2, n2, klpp,kupp, pwn_2nd,  2*klpp+kupp+1, ipiv_p, info )

    if ( info .ne. 0 ) write(*,*) 'error in computing p_op ...', info


  end subroutine p_op_2nd_Dirichlet

!==========================================================================================================================!

  subroutine p_op_2nd ( pwn_2nd , ipiv_p, kxs, kzs, mflag)

    !------------------------
    ! compute lu factorization of d^2 - (kx^2 + kz^2)*I 
    ! returns permutation vector along with result 
    !------------------------

    implicit none
    real, intent(in)                       :: kxs
    real, intent(in)                       :: kzs
    real, dimension(2*klpp+kupp+1,1:n2)      :: pwn_2nd
    integer, dimension(1:n2)               :: ipiv_p
    integer                                :: i,j,k,info
    logical, optional                      :: mflag

    pwn_2nd    = 0.d0
    ipiv_p = 0

    do j=1,n2
    do i = max(1, j-kupp), min(n2, j+klpp)

       if ( i .eq. j ) then
          if ( kxs .eq. 0.d0 .and. kzs .eq. 0.d0 .and. i .eq. 1) then
             pwn_2nd(kupp+klpp+1+i-j,j) = 1.d0                ! pin pressure 
          else
             pwn_2nd(kupp+klpp+1+i-j,j) = pwb_2nd(kupp+klpp+1+i-j,j) - (kxs + kzs )
          endif
       else
          if ( kxs .ne. 0.d0 .or. kzs .ne. 0.d0 .or. i .ne. 1 ) pwn_2nd(kupp+klpp+1+i-j,j) = pwb_2nd(kupp+klpp+1+i-j,j)
       endif
    enddo
    enddo

    call dgbtrf ( n2, n2, klpp,kupp, pwn_2nd,  2*klpp+kupp+1, ipiv_p, info )

    if ( info .ne. 0 ) write(*,*) 'error in computing p_op ...', info


  end subroutine p_op_2nd

!==========================================================================================================================!

  
  subroutine difwn ( lun, ifs, ife, kfs, kfe, wvisc, vtb ) 

    !-------------------------------
    ! compute explicit wall normal diffusion terms
    ! previous revisions had implementation using dgbmv 
    !-------------------------------

    implicit none 
    integer, intent(in)                                                  :: ifs, ife, kfs, kfe
    integer, intent(in)                                                  :: wvisc 
    real, intent(in)                                                     :: vtb
    real, dimension(1:nvel, ifs:ife, 1:n2, kfs:kfe)                      :: lun 
    integer                                                              :: i,j,k,n, ii 
    real                                                                 :: tsumu, tsumv, tsumw
    real                                                                 :: rre_u, rre_v, rre_w 
    real                                                                 :: c9r16, c1r16
    integer                                                              :: kf_start, kf_end, if_start, if_end 
    real                                                                 :: tsgs

    ! toggle eddy visc gradient 
    !---

    if ( wvisc .eq. eddy_visc ) then 
       tsgs = 1.d0  
    else  
       tsgs = 0.d0 
    endif 

    !----- 
    ! set indexing bounds 
    !----- 
    kf_start = kfs
    kf_end   = kfe
    if_start = ifs
    if_end   = ife

    lun   = 0.d0 
    rre_u = 1.d0 / re  + vtb 
    rre_v = 1.d0 / re  + vtb 
    rre_w = 1.d0 / re  + vtb 
    c9r16 = 9.d0 / 16.d0 
    c1r16 = 1.d0 / 16.d0 
 
    
    !--------------
    ! perform banded matrix multiplication 
    !--------------
    do k=kf_start,kf_end
    do i=if_start,if_end

      !---
      ! u,w
      !---
      do j = 1,n2
        if ( wvisc .eq. eddy_visc ) then
          rre_u = c9r16* ( vt(i-1,j,k) + vt(i  ,j,k)) - c1r16* ( vt(i-2,j,k) + vt(i+1,j,k))
          rre_w = c9r16* ( vt(i,j,k-1) + vt(i,j,k  )) - c1r16* ( vt(i,j,k-2) + vt(i,j,k+1))
        endif

        do ii=max(1,j-ku), min(n2,j+kl) 

          ! u
          if ( pat_bc .and. phase(i_u,i,k) .ne. 0 ) then ! (Kim 09.14.22)
            lun(i_u,i,j,k) = lun(i_u,i,j,k) + lub_pat(ku+kl+1+j-ii,ii)* rre_u* q(i_u,i,ii,k) + tsgs* fsgs*evb_u(ku+kl+1+j-ii,ii)* q(i_u,i,ii,k)
          else
            lun(i_u,i,j,k) = lun(i_u,i,j,k) + lub(ku+kl+1+j-ii,ii)* rre_u* q(i_u,i,ii,k) + tsgs* fsgs*evb_u(ku+kl+1+j-ii,ii)* q(i_u,i,ii,k)
          endif

          ! w
          if ( pat_bc .and. phase(i_w,i,k) .ne. 0 ) then ! (Kim 09.14.22)
            lun(i_w,i,j,k) = lun(i_w,i,j,k) + lub_pat(ku+kl+1+j-ii,ii)* rre_w* q(i_w,i,ii,k) + tsgs* fsgs*evb_w(ku+kl+1+j-ii,ii)* q(i_w,i,ii,k)
          else
            lun(i_w,i,j,k) = lun(i_w,i,j,k) + lub(ku+kl+1+j-ii,ii)* rre_w* q(i_w,i,ii,k) + tsgs* fsgs*evb_w(ku+kl+1+j-ii,ii)* q(i_w,i,ii,k)
          endif

        enddo
      enddo
 
      !---
      ! v 
      !---
      do j=1,n2-1
        if ( wvisc .eq. eddy_visc ) then  
          rre_v = c9r16* ( vt(i,j+1,k) + vt(i,j  ,k)) - c1r16* ( vt(i,j-1,k) + vt(i,j+2,k))
        endif

        do ii=max(1,j-ku), min(n2-1,j+kl)  
          lun(i_v,i,j,k) = lun(i_v,i,j,k) + lvb(ku+kl+1+j-ii,ii)* rre_v* q(i_v,i,ii,k) + tsgs* fsgs*evb_v(ku+kl+1+j-ii,ii)* q(i_v,i,ii,k)
        enddo
       enddo

    enddo
    enddo

  end subroutine difwn
!=======================================================================================================================================!
subroutine compute_lwb_2nd

!------------------------------
! wall normal diffusion operator along staggered mesh nodes 
! for u,w , stored in lapack banded matrix storage (pentadiagonal) 
!------------------------------ 

implicit none 
integer                                 :: i,j
real, dimension(5)                      :: coeff 
real, dimension(3)                      :: coeff_2nd
real                                    :: hm1, hm2, hp1, hp2

lwb   = 0.d0 
coeff = 0.d0
coeff_2nd = 0.d0

do j=1,n2 
do i = max(1, j-ku), min(n2,j+kl) 

if ( i .ge. 3 .and. i .le. n2-2 ) then  

hp1 = y(i_u,i+1) - y(i_u,i  ) 
hm1 = y(i_u,i  ) - y(i_u,i-1) 

call d2_coeff_2nd ( hm1, hp1, coeff_2nd )

if ( i-2 .eq. j ) then 
lwb(ku+kl+1+i-j,j) = 0.d0
elseif( i-1 .eq. j) then 
lwb(ku+kl+1+i-j,j) = coeff_2nd(1) 
elseif( i+1 .eq. j) then 
lwb(ku+kl+1+i-j,j) = coeff_2nd(3) 
elseif( i+2 .eq. j) then 
lwb(ku+kl+1+i-j,j) = 0.d0
elseif ( i .eq. j ) then
lwb(ku+kl+1+i-j,j) = coeff_2nd(2) 
else 
lwb(ku+kl+1+i-j,j) =  0.d0 
endif

endif
enddo
enddo

! ----

    i   = 1
    hp1 = y(i_u,i+1) - y(i_u,i )
    hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lwb(ku+kl+1,i  ) = cf(i_w,0,1)*coeff_2nd(1) + coeff_2nd(2)
    lwb(ku+kl  ,i+1) = cf(i_w,0,2)*coeff_2nd(1) + coeff_2nd(3)
    lwb(ku+kl-1,i+2) = 0.d0
 

    i   = 2
    hp1 = y(i_u,i+1) - y(i_u,i  )
    hm1 = y(i_u,i  ) - y(i_u,i-1)

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lwb(ku+kl+2,i-1) =  coeff_2nd(1)
    lwb(ku+kl+1,i  ) =  coeff_2nd(2)
    lwb(ku+kl  ,i+1) =  coeff_2nd(3)
    lwb(ku+kl-1,i+2) =  0.d0



    i   = n2-1
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = y(i_u,i+1) - y(i_u,i  )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lwb(ku+kl+3,i-2) = 0.d0
    lwb(ku+kl+2,i-1) = coeff_2nd(1)
    lwb(ku+kl+1,i  ) = coeff_2nd(2)
    lwb(ku+kl  ,i+1) = coeff_2nd(3)


    i   = n2
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lwb(ku+kl+3,i-2) = 0.d0
    lwb(ku+kl+2,i-1) = cf(i_w,0,2)*coeff_2nd(3) + coeff_2nd(1)
    lwb(ku+kl+1,i  ) = cf(i_w,0,1)*coeff_2nd(3) + coeff_2nd(2)
 

    end subroutine compute_lwb_2nd



!================================================================================================! 
subroutine compute_lub_2nd

!------------------------------
! wall normal diffusion operator along staggered mesh nodes 
! for u,w , stored in lapack banded matrix storage (pentadiagonal) 
!------------------------------ 

implicit none 
integer                                 :: i,j
real, dimension(5)                      :: coeff 
real, dimension(3)                      :: coeff_2nd
real                                    :: hm1, hm2, hp1, hp2

lub   = 0.d0 
coeff = 0.d0
coeff_2nd = 0.d0

do j=1,n2 
do i = max(1, j-ku), min(n2,j+kl) 

if ( i .ge. 3 .and. i .le. n2-2 ) then  

hp1 = y(i_u,i+1) - y(i_u,i  ) 
hm1 = y(i_u,i  ) - y(i_u,i-1) 

call d2_coeff_2nd ( hm1, hp1, coeff_2nd )

if ( i-2 .eq. j ) then 
lub(ku+kl+1+i-j,j) = 0.d0
elseif( i-1 .eq. j) then 
lub(ku+kl+1+i-j,j) = coeff_2nd(1) 
elseif( i+1 .eq. j) then 
lub(ku+kl+1+i-j,j) = coeff_2nd(3) 
elseif( i+2 .eq. j) then 
lub(ku+kl+1+i-j,j) = 0.d0
elseif ( i .eq. j ) then
lub(ku+kl+1+i-j,j) = coeff_2nd(2) 
else 
lub(ku+kl+1+i-j,j) =  0.d0 
endif

endif
enddo
enddo

! ----

    i   = 1
    hp1 = y(i_u,i+1) - y(i_u,i )
    hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lub(ku+kl+1,i  ) = cf(i_u,0,1)*coeff_2nd(1) + coeff_2nd(2)
    lub(ku+kl  ,i+1) = cf(i_u,0,2)*coeff_2nd(1) + coeff_2nd(3)
    lub(ku+kl-1,i+2) = 0.d0
 

    i   = 2
    hp1 = y(i_u,i+1) - y(i_u,i  )
    hm1 = y(i_u,i  ) - y(i_u,i-1)

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lub(ku+kl+2,i-1) =  coeff_2nd(1)
    lub(ku+kl+1,i  ) =  coeff_2nd(2)
    lub(ku+kl  ,i+1) =  coeff_2nd(3)
    lub(ku+kl-1,i+2) =  0.d0



    i   = n2-1
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = y(i_u,i+1) - y(i_u,i  )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lub(ku+kl+3,i-2) = 0.d0
    lub(ku+kl+2,i-1) = coeff_2nd(1)
    lub(ku+kl+1,i  ) = coeff_2nd(2)
    lub(ku+kl  ,i+1) = coeff_2nd(3)


    i   = n2
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lub(ku+kl+3,i-2) = 0.d0
    lub(ku+kl+2,i-1) = cf(i_u,0,2)*coeff_2nd(3) + coeff_2nd(1)
    lub(ku+kl+1,i  ) = cf(i_u,0,1)*coeff_2nd(3) + coeff_2nd(2)
 
   
    end subroutine compute_lub_2nd

!================================================================================================! 

subroutine compute_lub_pat

!------------------------------
! wall normal diffusion operator along staggered mesh nodes 
! for u,w , stored in lapack banded matrix storage (pentadiagonal) 
!------------------------------ 

implicit none 
integer                                 :: i,j
real, dimension(3)                      :: coeff_2nd
real                                    :: hm1, hm2, hp1, hp2

lub_pat   = 0.d0 
coeff_2nd = 0.d0

do j=1,n2 
  do i = max(1, j-ku), min(n2,j+kl) 

    if ( i .ge. 3 .and. i .le. n2-2 ) then  

      hp1 = y(i_u,i+1) - y(i_u,i  ) 
      hm1 = y(i_u,i  ) - y(i_u,i-1) 

      call d2_coeff_2nd ( hm1, hp1, coeff_2nd )

      if ( i-2 .eq. j ) then 
        lub_pat(ku+kl+1+i-j,j) = 0.d0
      elseif( i-1 .eq. j) then 
        lub_pat(ku+kl+1+i-j,j) = coeff_2nd(1) 
      elseif( i+1 .eq. j) then 
        lub_pat(ku+kl+1+i-j,j) = coeff_2nd(3) 
      elseif( i+2 .eq. j) then 
        lub_pat(ku+kl+1+i-j,j) = 0.d0
      elseif ( i .eq. j ) then
        lub_pat(ku+kl+1+i-j,j) = coeff_2nd(2) 
      else 
        lub_pat(ku+kl+1+i-j,j) =  0.d0 
      endif
    endif
  enddo
enddo

! ----

    i   = 1
    hp1 = y(i_u,i+1) - y(i_u,i )
    hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lub_pat(ku+kl+1,i  ) = (1.0)*coeff_2nd(1) + coeff_2nd(2)
    lub_pat(ku+kl  ,i+1) = (0.0)*coeff_2nd(1) + coeff_2nd(3)
    lub_pat(ku+kl-1,i+2) = 0.d0
 

    i   = 2
    hp1 = y(i_u,i+1) - y(i_u,i  )
    hm1 = y(i_u,i  ) - y(i_u,i-1)

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lub_pat(ku+kl+2,i-1) =  coeff_2nd(1)
    lub_pat(ku+kl+1,i  ) =  coeff_2nd(2)
    lub_pat(ku+kl  ,i+1) =  coeff_2nd(3)
    lub_pat(ku+kl-1,i+2) =  0.d0



    i   = n2-1
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = y(i_u,i+1) - y(i_u,i  )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lub_pat(ku+kl+3,i-2) = 0.d0
    lub_pat(ku+kl+2,i-1) = coeff_2nd(1)
    lub_pat(ku+kl+1,i  ) = coeff_2nd(2)
    lub_pat(ku+kl  ,i+1) = coeff_2nd(3)


    i   = n2
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lub_pat(ku+kl+3,i-2) = 0.d0
    lub_pat(ku+kl+2,i-1) = (0.0)*coeff_2nd(3) + coeff_2nd(1)
    lub_pat(ku+kl+1,i  ) = (1.0)*coeff_2nd(3) + coeff_2nd(2)
 
   
    end subroutine compute_lub_pat

!================================================================================================! 


subroutine compute_lvb_2nd

!--------------------------------
! wall normal diffusion operator along staggered mesh nodes 
! for v , stored in lapack banded matrix storage (pentadiagonal) 
!-------------------------------- 

implicit none 
integer                                :: i,j 
real, dimension(5)                     :: coeff
real, dimension(3)                     :: coeff_2nd
real                                   :: hm1, hm2, hp1, hp2 


lvb   = 0.d0 
coeff = 0.d0 
coeff_2nd = 0.d0

do j=1,n2-1 
do i = max(1, j-ku), min(n2-1,j+kl) 


if ( i .ge. 3 .and. i .le. n2-3) then 

hp1 = y(i_v,i+1) - y(i_v,i  ) 
hm1 = y(i_v,i  ) - y(i_v,i-1) 

call d2_coeff_2nd ( hm1, hp1, coeff_2nd )

if ( i-2 .eq. j ) then 
lvb(ku+kl+1+i-j,j) = 0.d0
elseif( i-1 .eq. j) then 
lvb(ku+kl+1+i-j,j) = coeff_2nd(1) 
elseif( i+1 .eq. j) then 
lvb(ku+kl+1+i-j,j) = coeff_2nd(3) 
elseif( i+2 .eq. j) then 
lvb(ku+kl+1+i-j,j) = 0.d0
elseif ( i .eq. j ) then
lvb(ku+kl+1+i-j,j) = coeff_2nd(2) 
else 
lvb(ku+kl+1+i-j,j) =  0.d0 
endif
endif
enddo
enddo


i   = 1
hm1 = y(i_v,i  ) - y(i_v,i-1)
hp1 = y(i_v,i+1) - y(i_v,i  )
hp2 = y(i_v,i+2) - y(i_v,i  )
hm2 = 2.d0* hm1

call d2_coeff_2nd( hm1, hp1, coeff_2nd)

lvb(kl+ku+1,i  ) = coeff_2nd(2)
lvb(kl+ku  ,i+1) = coeff_2nd(3)
lvb(kl+ku-1,i+2) = 0.d0



i   = 2
hm1 = y(i_v,i  ) - y(i_v,i-1)
hm2 = y(i_v,i  ) - y(i_v,i-2)
hp1 = y(i_v,i+1) - y(i_v,i  )
hp2 = y(i_v,i+2) - y(i_v,i  )

call d2_coeff_2nd( hm1, hp1, coeff_2nd)

lvb(kl+ku+2,i-1) = coeff_2nd(1)
lvb(kl+ku+1,i  ) = coeff_2nd(2)
lvb(kl+ku  ,i+1) = coeff_2nd(3)
lvb(kl+ku-1,i+2) = 0.d0


i   = n2-2
hm1 = y(i_v,i  ) - y(i_v,i-1)
hm2 = y(i_v,i  ) - y(i_v,i-2)
hp1 = y(i_v,i+1) - y(i_v,i  )
hp2 = y(i_v,i+2) - y(i_v,i  )

call d2_coeff_2nd( hm1, hp1, coeff_2nd)

lvb(kl+ku+3,i-2) = 0.d0
lvb(kl+ku+2,i-1) = coeff_2nd(1)
lvb(kl+ku+1,i  ) = coeff_2nd(2)
lvb(kl+ku  ,i+1) = coeff_2nd(3)



i   = n2-1
hm1 = y(i_v,i  ) - y(i_v,i-1)
hm2 = y(i_v,i  ) - y(i_v,i-2)
hp1 = y(i_v,i+1) - y(i_v,i  )
hp2 = 2.d0 * hp1 
    
call d2_coeff_2nd( hm1, hp1, coeff_2nd)
    
lvb(kl+ku+3,i-2) = 0.d0
lvb(kl+ku+2,i-1) = coeff_2nd(1)
lvb(kl+ku+1,i  ) = coeff_2nd(2)


end subroutine compute_lvb_2nd

!================================================================================================!
subroutine compute_lwb

!------------------------------
! wall normal diffusion operator along staggered mesh nodes 
! for u,w , stored in lapack banded matrix storage (pentadiagonal) 
!------------------------------ 

implicit none 
integer                                 :: i,j
real, dimension(5)                      :: coeff 
real, dimension(3)                      :: coeff_2nd
real                                    :: hm1, hm2, hp1, hp2

lwb   = 0.d0 
coeff = 0.d0
coeff = 0.d0

do j=1,n2 
do i = max(1, j-ku), min(n2,j+kl) 

if ( i .ge. 3 .and. i .le. n2-2 ) then  

hp1 = y(i_u,i+1) - y(i_u,i  ) 
hp2 = y(i_u,i+2) - y(i_u,i  ) 
hm1 = y(i_u,i  ) - y(i_u,i-1) 
hm2 = y(i_u,i  ) - y(i_u,i-2) 

call d2_coeff ( hm2, hm1, hp1, hp2, coeff ) 


if ( i-2 .eq. j ) then 
lwb(ku+kl+1+i-j,j) = coeff(1) 
elseif( i-1 .eq. j) then 
lwb(ku+kl+1+i-j,j) = coeff(2) 
elseif( i+1 .eq. j) then 
lwb(ku+kl+1+i-j,j) = coeff(4) 
elseif( i+2 .eq. j) then 
lwb(ku+kl+1+i-j,j) = coeff(5)
elseif ( i .eq. j ) then
lwb(ku+kl+1+i-j,j) = coeff(3) 
else 
lwb(ku+kl+1+i-j,j) =  0.d0 
endif

endif
enddo
enddo

! ----

    i   = 1
    hp1 = y(i_u,i+1) - y(i_u,i )
    hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lwb(ku+kl+1,i  ) = cf(i_w,0,1)*coeff_2nd(1) + coeff_2nd(2)
    lwb(ku+kl  ,i+1) = cf(i_w,0,2)*coeff_2nd(1) + coeff_2nd(3)
    lwb(ku+kl-1,i+2) = 0.d0


    i   = 2
    hp1 = y(i_u,i+1) - y(i_u,i  )
    hm1 = y(i_u,i  ) - y(i_u,i-1)

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lwb(ku+kl+2,i-1) =  coeff_2nd(1)
    lwb(ku+kl+1,i  ) =  coeff_2nd(2)
    lwb(ku+kl  ,i+1) =  coeff_2nd(3)
    lwb(ku+kl-1,i+2) =  0.d0



    i   = n2-1
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = y(i_u,i+1) - y(i_u,i  )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lwb(ku+kl+3,i-2) = 0.d0
    lwb(ku+kl+2,i-1) = coeff_2nd(1)
    lwb(ku+kl+1,i  ) = coeff_2nd(2)
    lwb(ku+kl  ,i+1) = coeff_2nd(3)


    i   = n2
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lwb(ku+kl+3,i-2) = 0.d0
    lwb(ku+kl+2,i-1) = cf(i_w,0,2)*coeff_2nd(3) + coeff_2nd(1)
    lwb(ku+kl+1,i  ) = cf(i_w,0,1)*coeff_2nd(3) + coeff_2nd(2)


    end subroutine compute_lwb

!================================================================================================! 
subroutine compute_lub 

!------------------------------
! wall normal diffusion operator along staggered mesh nodes 
! for u,w , stored in lapack banded matrix storage (pentadiagonal) 
!------------------------------ 

implicit none 
integer                                 :: i,j
real, dimension(5)                      :: coeff 
real, dimension(3)                      :: coeff_2nd
real                                    :: hm1, hm2, hp1, hp2

lub   = 0.d0 
coeff = 0.d0
coeff = 0.d0

do j=1,n2 
do i = max(1, j-ku), min(n2,j+kl) 

if ( i .ge. 3 .and. i .le. n2-2 ) then  

hp1 = y(i_u,i+1) - y(i_u,i  ) 
hp2 = y(i_u,i+2) - y(i_u,i  ) 
hm1 = y(i_u,i  ) - y(i_u,i-1) 
hm2 = y(i_u,i  ) - y(i_u,i-2) 

call d2_coeff ( hm2, hm1, hp1, hp2, coeff ) 


if ( i-2 .eq. j ) then 
lub(ku+kl+1+i-j,j) = coeff(1) 
elseif( i-1 .eq. j) then 
lub(ku+kl+1+i-j,j) = coeff(2) 
elseif( i+1 .eq. j) then 
lub(ku+kl+1+i-j,j) = coeff(4) 
elseif( i+2 .eq. j) then 
lub(ku+kl+1+i-j,j) = coeff(5)
elseif ( i .eq. j ) then
lub(ku+kl+1+i-j,j) = coeff(3) 
else 
lub(ku+kl+1+i-j,j) =  0.d0 
endif

endif
enddo
enddo

! ----

    i   = 1
    hp1 = y(i_u,i+1) - y(i_u,i )
    hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lub(ku+kl+1,i  ) = cf(i_u,0,1)*coeff_2nd(1) + coeff_2nd(2)
    lub(ku+kl  ,i+1) = cf(i_u,0,2)*coeff_2nd(1) + coeff_2nd(3)
    lub(ku+kl-1,i+2) = 0.d0


    i   = 2
    hp1 = y(i_u,i+1) - y(i_u,i  )
    hm1 = y(i_u,i  ) - y(i_u,i-1)

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lub(ku+kl+2,i-1) =  coeff_2nd(1)
    lub(ku+kl+1,i  ) =  coeff_2nd(2)
    lub(ku+kl  ,i+1) =  coeff_2nd(3)
    lub(ku+kl-1,i+2) =  0.d0



    i   = n2-1
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = y(i_u,i+1) - y(i_u,i  )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lub(ku+kl+3,i-2) = 0.d0
    lub(ku+kl+2,i-1) = coeff_2nd(1)
    lub(ku+kl+1,i  ) = coeff_2nd(2)
    lub(ku+kl  ,i+1) = coeff_2nd(3)


    i   = n2
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lub(ku+kl+3,i-2) = 0.d0
    lub(ku+kl+2,i-1) = cf(i_u,0,2)*coeff_2nd(3) + coeff_2nd(1)
    lub(ku+kl+1,i  ) = cf(i_u,0,1)*coeff_2nd(3) + coeff_2nd(2)


    end subroutine compute_lub

!================================================================================================! 

subroutine compute_lvb 

!--------------------------------
! wall normal diffusion operator along staggered mesh nodes 
! for v , stored in lapack banded matrix storage (pentadiagonal) 
!-------------------------------- 

implicit none 
integer                                :: i,j 
real, dimension(5)                     :: coeff
real, dimension(3)                     :: coeff_2nd
real                                   :: hm1, hm2, hp1, hp2 




lvb   = 0.d0 
coeff = 0.d0 
coeff_2nd = 0.d0


do j=1,n2-1 
do i = max(1, j-ku), min(n2-1,j+kl) 


if ( i .ge. 3 .and. i .le. n2-3) then 

hp1 = y(i_v,i+1) - y(i_v,i  ) 
hp2 = y(i_v,i+2) - y(i_v,i  ) 
hm1 = y(i_v,i  ) - y(i_v,i-1) 
hm2 = y(i_v,i  ) - y(i_v,i-2) 

call d2_coeff ( hm2, hm1, hp1, hp2, coeff ) 


if ( i-2 .eq. j ) then 
lvb(ku+kl+1+i-j,j) = coeff(1) 
elseif( i-1 .eq. j) then 
lvb(ku+kl+1+i-j,j) = coeff(2) 
elseif( i+1 .eq. j) then 
lvb(ku+kl+1+i-j,j) = coeff(4) 
elseif( i+2 .eq. j) then 
lvb(ku+kl+1+i-j,j) = coeff(5)
elseif ( i .eq. j ) then
lvb(ku+kl+1+i-j,j) = coeff(3) 
else 
lvb(ku+kl+1+i-j,j) =  0.d0 
endif
endif
enddo
enddo


i   = 1
hm1 = y(i_v,i  ) - y(i_v,i-1)
hp1 = y(i_v,i+1) - y(i_v,i  )
hp2 = y(i_v,i+2) - y(i_v,i  )
hm2 = 2.d0* hm1

call d2_coeff_2nd( hm1, hp1, coeff_2nd)

lvb(kl+ku+1,i  ) = coeff_2nd(2)
lvb(kl+ku  ,i+1) = coeff_2nd(3)
lvb(kl+ku-1,i+2) = 0.d0



i   = 2
hm1 = y(i_v,i  ) - y(i_v,i-1)
hm2 = y(i_v,i  ) - y(i_v,i-2)
hp1 = y(i_v,i+1) - y(i_v,i  )
hp2 = y(i_v,i+2) - y(i_v,i  )

call d2_coeff_2nd( hm1, hp1, coeff_2nd)

lvb(kl+ku+2,i-1) = coeff_2nd(1)
lvb(kl+ku+1,i  ) = coeff_2nd(2)
lvb(kl+ku  ,i+1) = coeff_2nd(3)
lvb(kl+ku-1,i+2) = 0.d0


i   = n2-2
hm1 = y(i_v,i  ) - y(i_v,i-1)
hm2 = y(i_v,i  ) - y(i_v,i-2)
hp1 = y(i_v,i+1) - y(i_v,i  )
hp2 = y(i_v,i+2) - y(i_v,i  )

call d2_coeff_2nd( hm1, hp1, coeff_2nd)

lvb(kl+ku+3,i-2) = 0.d0
lvb(kl+ku+2,i-1) = coeff_2nd(1)
lvb(kl+ku+1,i  ) = coeff_2nd(2)
lvb(kl+ku  ,i+1) = coeff_2nd(3)



i   = n2-1
hm1 = y(i_v,i  ) - y(i_v,i-1)
hm2 = y(i_v,i  ) - y(i_v,i-2)
hp1 = y(i_v,i+1) - y(i_v,i  )
hp2 = 2.d0 * hp1 
    
call d2_coeff_2nd( hm1, hp1, coeff_2nd)
    
lvb(kl+ku+3,i-2) = 0.d0
lvb(kl+ku+2,i-1) = coeff_2nd(1)
lvb(kl+ku+1,i  ) = coeff_2nd(2)


end subroutine compute_lvb


 !==========================================================================================================================!
  subroutine compute_pwb_2nd_Dirichlet

    implicit none
    integer                   :: i , j
    real, dimension(-1:1)     :: jac
    real, dimension( 1:3)     :: coeff_2nd
    real                      :: c1r64, c9r64, c81r64, c1r576, c9r192, c1r288, c9r96, hy0
    real                      :: cfp01, cfp02
    real                      :: h1p, h2p

    pwb_2nd_Dirichlet   = 0.d0
    coeff_2nd = 0.d0
    jac   = 0.d0

    c1r64 = 1.d0 / 64.d0
    c9r64 = 9.d0 / 64.d0
    c81r64= 81.d0/ 64.d0
    c1r576= 1.d0 / 576.d0
    c9r192= 9.d0 / 192.d0
    c1r288= 1.d0 / 288.d0
    c9r96 = 9.d0 / 96.d0

    hy0   = 2.d0 / dble(n2)
    h1p   = ( y(i_u,1) - y(i_u,0) ) + ( y(i_u,1) - y(i_u,0) )
    h2p   = ( y(i_u,2) - y(i_u,0) ) + ( y(i_u,1) - y(i_u,0) )

    cfp01 = -h2p
    cfp02 =  h1p

    do j=1,n2
    do i = max(1, j-klpp), min(n2,j+klpp)


       if ( i .ge. 2 .and. i .le. n2-1 ) then

          jac(0 )  =  jacb(i_u,i  )
          jac(-1)  =  jacb(i_v,i-1)
          jac( 1)  =  jacb(i_v,i  )

          call d2p_coeff_2nd( jac, coeff_2nd )


          if ( i-1 .eq. j ) then
             pwb_2nd_Dirichlet(kupp+klpp+1+i-j,j) = coeff_2nd(1)

          elseif ( i+1 .eq. j) then
             pwb_2nd_Dirichlet(kupp+klpp+1+i-j,j) = coeff_2nd(3)

          elseif ( i .eq.   j) then
             pwb_2nd_Dirichlet(kupp+klpp+1+i-j,j) = coeff_2nd(2)

          else
             pwb_2nd_Dirichlet(kupp+klpp+1+i-j,j) =  0.d0

          endif
       endif

    enddo
    enddo


    i       = 1
    jac( 0)           = jacb(i_u,i  )
    jac(-1)           = jacb(i_v,i-1)
    jac( 1)           = jacb(i_v,i  )

    pwb_2nd_Dirichlet(kupp+klpp+1, 1) = 1.d0  ! -  1.d0/jac(0)/jac(1) 
    pwb_2nd_Dirichlet(kupp+klpp  , 2) = 0.d0  !  1.d0/jac(0)/jac(1) 

    i       = n2
    jac( 0)           = jacb(i_u, i  )
    jac(-1)           = jacb(i_v, i)
    jac( 1)           = jacb(i_v, i-1)

    pwb_2nd_Dirichlet(kupp+klpp+2, i-1) = 0.d0 !   1.d0/jac(0)/jac(1) 
    pwb_2nd_Dirichlet(kupp+klpp+1, i) =   1.d0 !- 1.d0/jac(0)/jac(1) 



  end subroutine compute_pwb_2nd_Dirichlet

!================================================================================================!


  subroutine compute_pwb_2nd

    implicit none
    integer                   :: i , j
    real, dimension(-1:1)     :: jac
    real, dimension( 1:3)     :: coeff_2nd
    real                      :: c1r64, c9r64, c81r64, c1r576, c9r192, c1r288, c9r96, hy0
    real                      :: cfp01, cfp02
    real                      :: h1p, h2p


    pwb_2nd   = 0.d0
    coeff_2nd = 0.d0
    jac   = 0.d0


    c1r64 = 1.d0 / 64.d0
    c9r64 = 9.d0 / 64.d0
    c81r64= 81.d0/ 64.d0
    c1r576= 1.d0 / 576.d0
    c9r192= 9.d0 / 192.d0
    c1r288= 1.d0 / 288.d0
    c9r96 = 9.d0 / 96.d0

    hy0   = 2.d0 / dble(n2)
    h1p   = ( y(i_u,1) - y(i_u,0) ) + ( y(i_u,1) - y(i_u,0) )
    h2p   = ( y(i_u,2) - y(i_u,0) ) + ( y(i_u,1) - y(i_u,0) )

    cfp01 = -h2p
    cfp02 =  h1p

    do j=1,n2
    do i = max(1, j-klpp), min(n2,j+klpp)


       if ( i .ge. 2 .and. i .le. n2-1 ) then

          jac(0 )  =  jacb(i_u,i  )
          jac(-1)  =  jacb(i_v,i-1)
          jac( 1)  =  jacb(i_v,i  )

          call d2p_coeff_2nd( jac, coeff_2nd )


          if ( i-1 .eq. j ) then
             pwb_2nd(kupp+klpp+1+i-j,j) = coeff_2nd(1)

          elseif ( i+1 .eq. j) then
             pwb_2nd(kupp+klpp+1+i-j,j) = coeff_2nd(3)

          elseif ( i .eq.   j) then
             pwb_2nd(kupp+klpp+1+i-j,j) = coeff_2nd(2)

          else
             pwb_2nd(kupp+klpp+1+i-j,j) =  0.d0

          endif
       endif

    enddo
    enddo


    i       = 1
    jac( 0)           = jacb(i_u,i  )
    jac(-1)           = jacb(i_v,i-1)
    jac( 1)           = jacb(i_v,i  )

    pwb_2nd(kupp+klpp+1, 1) = - 1.d0/jac(0)/jac(1) 
    pwb_2nd(kupp+klpp  , 2) =   1.d0/jac(0)/jac(1) 

    i       = n2
    jac( 0)           = jacb(i_u, i  )
    jac(-1)           = jacb(i_v, i)
    jac( 1)           = jacb(i_v, i-1)

    pwb_2nd(kupp+klpp+2, i-1) =   1.d0/jac(0)/jac(1) 
    pwb_2nd(kupp+klpp+1, i) =   - 1.d0/jac(0)/jac(1) 



  end subroutine compute_pwb_2nd

!================================================================================================!

  subroutine compute_pwb 

    !-------------------------------------- 
    ! wall normal diffusion operator for "pressure" variable, phi , in 
    ! fractional step 
    !-------------------------------------- 

    implicit none 
    integer                   :: i , j 
    real, dimension(-3:3)     :: jac 
    real, dimension( 1:7)     :: coeff 
    real                      :: c1r64, c9r64, c81r64, c1r576, c9r192, c1r288, c9r96, hy0 
    real                      :: cslip12, cslip11
    
    pwb   = 0.d0 
    coeff = 0.d0
    jac   = 0.d0 


    c1r64 = 1.d0 / 64.d0 
    c9r64 = 9.d0 / 64.d0 
    c81r64= 81.d0/ 64.d0 
    c1r576= 1.d0 / 576.d0 
    c9r192= 9.d0 / 192.d0 
    c1r288= 1.d0 / 288.d0 
    c9r96 = 9.d0 / 96.d0 
    
    hy0   = 2.d0 / dble(n2) 

    cslip12 = 0.d0
    cslip11 = 1.d0


    do j=1,n2
    do i = max(1, j-kup), min(n2,j+klp) 


       if ( i .ge. 4 .and. i .le. n2-3 ) then 


!.... compute necessary jacobians , cf vasilyev 2000 
!.... jac^* from derivation must use spacings of along v
          
          !jac(0 ) = ( y(i_u,i+1) - y(i_u,i-1) )/ 2.d0/ hy0 
          !jac(-1) = ( y(i_v,i  ) - y(i_v,i-2) )/ 2.d0/ hy0 
          !jac( 1) = ( y(i_v,i+1) - y(i_v,i-1) )/ 2.d0/ hy0 
          !jac(-3) = ( y(i_v,i-1) - y(i_v,i-3) )/ 2.d0/ hy0 
          !jac( 3) = ( y(i_v,i+2) - y(i_v,i  ) )/ 2.d0/ hy0 

          jac(0 )  =  jacb(i_u,i  ) 
          jac(-1)  =  jacb(i_v,i-1)
          jac( 1)  =  jacb(i_v,i  ) 
          jac(-3)  =  jacb(i_v,i-2) 
          jac( 3)  =  jacb(i_v,i+1) 

!          jac = 1.d0 
          call d2p_coeff( jac, coeff )           
        
          if ( i-3 .eq. j ) then 
             pwb(kup+klp+1+i-j,j) = coeff(1) 
             
          elseif ( i-2 .eq. j) then 
             pwb(kup+klp+1+i-j,j) = coeff(2)
             
          elseif ( i-1 .eq. j) then 
             pwb(kup+klp+1+i-j,j) = coeff(3) 
             
          elseif ( i+1 .eq. j) then 
             pwb(kup+klp+1+i-j,j) = coeff(5)

          elseif ( i+2 .eq. j) then 
             pwb(kup+klp+1+i-j,j) = coeff(6) 
          elseif ( i+3 .eq. j) then 
             pwb(kup+klp+1+i-j,j) = coeff(7)

          elseif ( i .eq.   j) then 
             pwb(kup+klp+1+i-j,j) = coeff(4) 

          else
             pwb(kup+klp+1+i-j,j) =  0.d0 

          endif
       endif
      
    enddo
    enddo

    jac = 0.0d0 

!.... boundary condition rows

!**** new boundary treatment of phi, w/o using boundary condition on phi**** 
!.... currently only supported for uniform grid 
!.... (DEBUG: testing nonunif treatment, values after the comments are 
!....         those for an uniform grid                                 ) 


    i       = 1 
    jac( 0)           = jacb(i_u,i  ) 
    jac( 1)           = jacb(i_v,i  ) 
    jac( 3)           = jacb(i_v,i+1) 

    pwb(kup+klp+1, 1) = -c81r64/jac(0)/jac(1) -c1r576/jac(0)/jac(3) +c9r96/jac(0)/jac(1)                       & 
                        +cf(i_v,-1,2)*c1r576/jac(0)/jac(3) +cf(i_v,-1,1)*c1r288/jac(0)/jac(1)                  & 
                        -cf(i_v,-1,1)*c9r192/jac(0)/jac(1)                                                     !-25.d0/24.d0 

    pwb(kup+klp  , 2) =  c81r64/jac(0)/jac(1) +c9r192/jac(0)/jac(3) -c9r192/jac(0)/jac(1)                      & 
                        +cf(i_v,-1,1)*c9r192/jac(0)/jac(1) -cf(i_v,-1,2)*c9r192/jac(0)/jac(3)                  & 
                        -cf(i_v,-1,1)*c1r576/jac(0)/jac(1)                                                     !13.d0/12.d0 

    pwb(kup+klp-1, 3) = -c9r192/jac(0)/jac(1) -c9r192/jac(0)/jac(3)                                            & 
                        -cf(i_v,-1,1)*c1r576/jac(0)/jac(1) +cf(i_v,-1,2)*c9r192/jac(0)/jac(3)                  !-1.d0/24.d0  

    pwb(kup+klp-2, 4) =  c1r576/jac(0)/jac(3) - cf(i_v,-1,2)*c1r576/jac(0)/jac(3)                              !0.d0        



    i       = 2 
    jac( 0)           = jacb(i_u,i  )
    jac(-1)           = jacb(i_v,i-1) 
    jac( 1)           = jacb(i_v,i  ) 
    jac( 3)           = jacb(i_v,i+1)

    
    pwb(kup+klp+2, 1) =  c9r192/jac(0)/jac(1) +c81r64/jac(0)/jac(-1) -c9r96/jac(0)/jac(-1)                      !39.d0/32.d0   
    pwb(kup+klp+1, 2) = -c1r576/jac(0)/jac(3) -c81r64/jac(0)/jac(-1) -c81r64/jac(0)/jac(1) + c9r192/jac(0)/jac(-1) !-179.d0/72.d0 
    pwb(kup+klp  , 3) =  c81r64/jac(0)/jac(1) +c9r192/jac(0)/jac(-1) +c9r192/jac(0)/jac(3)                      !87.d0/64.d0   
    pwb(kup+klp-1, 4) = -c9r192/jac(0)/jac(1) -c9r192/jac(0)/jac(3)                                             !-3.d0/32.d0   
    pwb(kup+klp-2, 5) =  c1r576/jac(0)/jac(3)                                                                   !1.d0/576.d0   


    i       = 3
    jac( 0)           = jacb(i_u,i  )
    jac(-1)           = jacb(i_v,i-1) 
    jac( 1)           = jacb(i_v,i  ) 
    jac(-3)           = jacb(i_v,i-2) 
    jac( 3)           = jacb(i_v,i+1) 


    pwb(kup+klp+3, 1) = -c9r192/jac(0)/jac(-3) -c9r192/jac(0)/jac(-1) +c1r288/jac(0)/jac(-3)                    !-13.d0/144.d0 
    pwb(kup+klp+2, 2) =  c9r192/jac(0)/jac(-3) +c81r64/jac(0)/jac(-1) +c9r192/jac(0)/jac( 1) -c1r576/jac(0)/jac(-3)!391.d0/288.d0 
    pwb(kup+klp+1, 3) = -c1r576/jac(0)/jac(-3) -c1r576/jac(0)/jac(3 ) -c81r64/jac(0)/jac( 1) -c81r64/jac(0)/jac(-1)!-365.d0/144.d0 
    pwb(kup+klp  , 4) =  c9r192/jac(0)/jac( 3) +c81r64/jac(0)/jac(1 ) +c9r192/jac(0)/jac(-1)                    !87.d0/64.d0    
    pwb(kup+klp-1, 5) = -c9r192/jac(0)/jac( 3) -c9r192/jac(0)/jac(1 )                                           !-3.d0/32.d0   
    pwb(kup+klp-2, 6) =  c1r576/jac(0)/jac( 3)                                                                  !1.d0/576.d0   


    i       = n2-2
    jac( 0)           = jacb(i_u,i  )
    jac(-1)           = jacb(i_v,i  ) 
    jac( 1)           = jacb(i_v,i-1) 
    jac( 3)           = jacb(i_v,i-2) 
    jac(-3)           = jacb(i_v,i+1) 


    pwb(kup+klp+4,i-3) =  c1r576/jac(0)/jac( 3)                                                                  !1.d0/576.d0     
    pwb(kup+klp+3,i-2) = -c9r192/jac(0)/jac( 3) -c9r192/jac(0)/jac(1 )                                           !-3.d0/32.d0     
    pwb(kup+klp+2,i-1) =  c9r192/jac(0)/jac( 3) +c81r64/jac(0)/jac(1 ) +c9r192/jac(0)/jac(-1)                    !87.d0/64.d0     
    pwb(kup+klp+1,i  ) = -c1r576/jac(0)/jac(-3) -c1r576/jac(0)/jac(3 ) -c81r64/jac(0)/jac( 1) -c81r64/jac(0)/jac(-1)!-365.d0/144.d0  
    pwb(kup+klp  ,i+1) =  c9r192/jac(0)/jac(-3) +c81r64/jac(0)/jac(-1) +c9r192/jac(0)/jac( 1) -c1r576/jac(0)/jac(-3)!391.d0/288.d0   
    pwb(kup+klp-1,i+2) = -c9r192/jac(0)/jac(-3) -c9r192/jac(0)/jac(-1) +c1r288/jac(0)/jac(-3)                    !-13.d0/144.d0   
    
    
    i       = n2-1 
    jac( 0)           = jacb(i_u, i  ) 
    jac(-1)           = jacb(i_v, i  ) 
    jac( 1)           = jacb(i_v, i-1) 
    jac( 3)           = jacb(i_v, i-2) 
 

    pwb(kup+klp+4, i-3) =  c1r576/jac(0)/jac(3)                                                                  !1.d0/576.d0    
    pwb(kup+klp+3, i-2) = -c9r192/jac(0)/jac(1) -c9r192/jac(0)/jac(3)                                            !-3.d0/32.d0    
    pwb(kup+klp+2, i-1) =  c81r64/jac(0)/jac(1) +c9r192/jac(0)/jac(-1) +c9r192/jac(0)/jac(3)                     !87.d0/64.d0    
    pwb(kup+klp+1, i  ) = -c1r576/jac(0)/jac(3) -c81r64/jac(0)/jac(-1) -c81r64/jac(0)/jac(1) + c9r192/jac(0)/jac(-1) !-179.d0/72.d0  
    pwb(kup+klp  , i+1) =  c9r192/jac(0)/jac(1) +c81r64/jac(0)/jac(-1) -c9r96/jac(0)/jac(-1)                     !39.d0/32.d0    
    


    i       = n2  
    jac( 0)           = jacb(i_u,i  ) 
    jac( 1)           = jacb(i_v,i-1)
    jac( 3)           = jacb(i_v,i-2) 

    pwb(kup+klp+4, i-3) =  c1r576/jac(0)/jac(3) - cf(i_v,-1,2)*c1r576/jac(0)/jac(3)                              !0.d0        
    pwb(kup+klp+3, i-2) = -c9r192/jac(0)/jac(1) -c9r192/jac(0)/jac(3)                                            & 
                          -cf(i_v,-1,1)*c1r576/jac(0)/jac(1) +cf(i_v,-1,2)*c9r192/jac(0)/jac(3)                  !-1.d0/24.d0 
    pwb(kup+klp+2, i-1) =  c81r64/jac(0)/jac(1) +c9r192/jac(0)/jac(3) -c9r192/jac(0)/jac(1)                      & 
                          +cf(i_v,-1,1)*c9r192/jac(0)/jac(1) -cf(i_v,-1,2)*c9r192/jac(0)/jac(3)                  & 
                          -cf(i_v,-1,1)*c1r576/jac(0)/jac(1)                                                     !13.d0/12.d0 
    pwb(kup+klp+1, i  ) = -c81r64/jac(0)/jac(1) -c1r576/jac(0)/jac(3) +c9r96/jac(0)/jac(1)                       & 
                          +cf(i_v,-1,2)*c1r576/jac(0)/jac(3) +cf(i_v,-1,1)*c1r288/jac(0)/jac(1)                  & 
                          -cf(i_v,-1,1)*c9r192/jac(0)/jac(1)                                                     !-25.d0/24.d0 






  end subroutine compute_pwb


!=================================================================================================!

  subroutine d2_coeff ( hm2, hm1 , hp1, hp2, coeff ) 

    implicit none 
    real, intent(in )               :: hm2, hm1, hp1, hp2 
    real, dimension(5), intent(out) :: coeff


    coeff(1) = -2.d0* ( hm1*hp1 + hm1*hp2 - hp1*hp2 ) / & 
                (hm2* ( hm2-hm1)*(hm2+hp1)*(hm2+hp2)  )

    
    coeff(2) =  2.d0* ( hm2*hp1 + hm2*hp2 - hp1*hp2 ) / & 
                (hm1* ( hm2-hm1)*(hm1+hp1)*(hm1+hp2)  )


    coeff(4) =  2.d0* ( hm1*hm2 - hm1*hp2 - hm2*hp2 ) / & 
                (hp1* ( hm1+hp1)*(hm2+hp1)*(hp1-hp2)  )

    coeff(5) = -2.d0* ( hm1*hm2 - hm1*hp1 - hm2*hp1 ) / & 
                (hp2* ( hp1-hp2)*(hm1+hp2)*(hm2+hp2)  ) 


    coeff(3) =  2.d0* ( hm1*hm2 - hm1*hp1 - hm2*hp1 - hm1*hp2 - hm2*hp2 + hp1*hp2 ) / & 
                      ( hm2 * hm1 * hp1 * hp2 ) 


  end subroutine d2_coeff

  !==============================================================================================!

    ! ----
    ! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
    ! calculate coefficients for dirichlet boundary conditions.
    ! using taylor expansion near the boundary in non-uniform grid setting
    ! ----

   subroutine d2_coeff_2nd ( hm1 , hp1, coeff_2nd )

     implicit none
     real, intent(in )               :: hm1, hp1
     real, dimension(3), intent(out) :: coeff_2nd


     coeff_2nd(1) = 2.d0 / &
     ( hm1*hm1 + hm1*hp1  )


     coeff_2nd(2) = -2.d0 / &
     ( hm1*hp1 )

     coeff_2nd(3) = 2.d0 / &
     ( hp1*hp1 +  hm1*hp1  )


   end subroutine d2_coeff_2nd


  !==============================================================================================!


  subroutine d2p_coeff ( jac, coeff ) 

    implicit none 
    real, dimension(-3:3)           :: jac 
    real, dimension( 1:7), intent(out)          :: coeff 
    real                                        :: c1r576, c9r192, c81r64 

    c1r576 = 1.d0 / 576.d0 
    c9r192 = 9.d0 / 192.d0 
    c81r64=  81.d0/ 64.d0 

    coeff(1) =   c1r576/ jac(0)/ jac(-3)

    coeff(2) =  -c9r192/ jac(0)/ jac(-3) - & 
                 c9r192/ jac(0)/ jac(-1)

    coeff(3) =   c9r192/ jac(0)/ jac( 1) + & 
                 c81r64/ jac(0)/ jac(-1) + & 
                 c9r192/ jac(0)/ jac(-3)


    coeff(4) =  -c1r576/ jac(0)/ jac( 3) - & 
                 c1r576/ jac(0)/ jac(-3) - & 
                 c81r64/ jac(0)/ jac( 1) - & 
                 c81r64/ jac(0)/ jac(-1) 

    coeff(5) =   c81r64/ jac(0)/ jac( 1) + & 
                 c9r192/ jac(0)/ jac(-1) + & 
                 c9r192/ jac(0)/ jac( 3) 

    coeff(6) =  -c9r192/ jac(0)/ jac( 3) - & 
                 c9r192/ jac(0)/ jac( 1)

    coeff(7) =   c1r576/ jac(0)/ jac( 3) 
 
  end subroutine d2p_coeff


  !==============================================================================================!


  subroutine d2p_coeff_2nd ( jac, coeff_2nd )

    implicit none
    real, dimension(-1:1)                       :: jac
    real, dimension( 1:3), intent(out)          :: coeff_2nd


    coeff_2nd(1) =   1.d0/ jac(0)/ jac(-1)

    coeff_2nd(2) =  -1.d0/ jac(0)/ jac(-1) - 1.d0/ jac(0)/ jac( 1)

    coeff_2nd(3) =   1.d0/ jac(0)/ jac( 1) 


  end subroutine d2p_coeff_2nd



!==============================================================================================================================!


  subroutine update_slip

    implicit none
    integer                                 :: i
    real                                    :: h1u, h2u   ! bc coeffs 
    real, dimension(3)                      :: coeff_2nd  ! lub/lwb
    real                                    :: hm1, hp1   ! lub/lwb


    b_d_x = slip_mean/re + slip_del/re * sin( 2.d0*pi/slip_per*re * (time-init_time) )
    b_d_z = slip_mean/re + slip_del/re * sin( 2.d0*pi/slip_per*re * (time-init_time) )

!.... update cf_slip
    h1u = y(i_u,1) - y(i_u,0) 
    h2u = y(i_u,2) - y(i_u,0) 

    cf_slip(i_u, 0,1) = ( h1u*h1u - h2u*h2u * (b_d_x-h1u)/(b_d_x+h2u) ) / &
        ( h1u*h1u - h2u*h2u * (b_d_x+h1u) / (b_d_x+h2u) )
    cf_slip(i_u, 0,2) = 2.d0 * h1u*h1u * h1u / (b_d_x+h1u) / &
        ( h2u*h2u - h1u*h1u * (b_d_x+h2u) / (b_d_x+h1u) )
    cf_slip(i_w, 0,1) = ( h1u*h1u - h2u*h2u * (b_d_z-h1u)/(b_d_z+h2u) ) / &
        ( h1u*h1u - h2u*h2u * (b_d_z+h1u) / (b_d_z+h2u) )
    cf_slip(i_w, 0,2) = 2.d0 * h1u*h1u * h1u / (b_d_z+h1u) / &
        ( h2u*h2u - h1u*h1u * (b_d_z+h2u) / (b_d_z+h1u) )

!.... update i = 1 in lub and lwb
    i   = 1
    hp1 = y(i_u,i+1) - y(i_u,i )
    hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )
    call d2_coeff_2nd( hm1, hp1, coeff_2nd)
   
    lub(ku+kl+1,i  ) = cf_slip(i_u,0,1)*coeff_2nd(1) + coeff_2nd(2)
    lub(ku+kl  ,i+1) = cf_slip(i_u,0,2)*coeff_2nd(1) + coeff_2nd(3)
    lwb(ku+kl+1,i  ) = cf_slip(i_w,0,1)*coeff_2nd(1) + coeff_2nd(2)
    lwb(ku+kl  ,i+1) = cf_slip(i_w,0,2)*coeff_2nd(1) + coeff_2nd(3)

    if ( transportee .and. match_bc ) then ! (Kim 09.13.22)
      lub_tee(ku+kl+1,i  ) = cf_slip(i_u,0,1)*coeff_2nd(1) + coeff_2nd(2)
      lub_tee(ku+kl  ,i+1) = cf_slip(i_u,0,2)*coeff_2nd(1) + coeff_2nd(3)
    endif
  
!.... update i = n2 in lub and lwb
    i   = n2
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )
    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lub(ku+kl+2,i-1) = cf_slip(i_u,0,2)*coeff_2nd(3) + coeff_2nd(1)
    lub(ku+kl+1,i  ) = cf_slip(i_u,0,1)*coeff_2nd(3) + coeff_2nd(2)
    lwb(ku+kl+2,i-1) = cf_slip(i_w,0,2)*coeff_2nd(3) + coeff_2nd(1)
    lwb(ku+kl+1,i  ) = cf_slip(i_w,0,1)*coeff_2nd(3) + coeff_2nd(2)

    if ( transportee .and. match_bc ) then ! (Kim 09.13.22)
      lub_tee(ku+kl+2,i-1) = cf_slip(i_u,0,2)*coeff_2nd(3) + coeff_2nd(1)
      lub_tee(ku+kl+1,i  ) = cf_slip(i_u,0,1)*coeff_2nd(3) + coeff_2nd(2)
    endif
     
  end subroutine update_slip


!==============================================================================================!


  subroutine bcuvw 
    
!.... u,v,w boundary conditions 
    implicit none 
    integer                             :: i,j,k,i_q , jss
    integer, dimension(mpi_status_size) :: status 
    real                                :: dx, dy, hy0, jac, djy, dz
    real                                :: t_c, uw1


    hy0  = 2.d0/ dble(n2)
    jac  = jacb(i_u, 1 )
    djy  = hy0*jac
    dx = chlx * pi / dble(n1)
    dz = chlz * pi / dble(n3)
!.... apply boundary conditions in wn direction, 
    do k=ks,ke
    do i=is,ie


    q(i_u,i,   0,k) = - q(i_u,i,1,k)     ! u-1 = u0     makes the boundary u velocity zero if we use second order scheme with staggered mesh
    q(i_u,i,  -1,k) = 0.d0               ! extra ghost cell. we do not use this. 
    q(i_u,i,  -2,k) = 0.d0               ! extra ghost cell. we do not use this. 
    q(i_u,i,n2+1,k) = - q(i_u,i,n2,k)
    q(i_u,i,n2+2,k) = 0.d0
    q(i_u,i,n2+3,k) = 0.d0

    q(i_w,i,   0,k) = - q(i_w,i,1,k)
    q(i_w,i,  -1,k) = 0.d0
    q(i_w,i,  -2,k) = 0.d0
    q(i_w,i,n2+1,k) = - q(i_w,i,n2,k)
    q(i_w,i,n2+2,k) = 0.d0
    q(i_w,i,n2+3,k) = 0.d0

!.... wall oscillations
    if ( wall_bc ) then

      t_c = (time-init_time)*2*pi/wall_per*re
      uw1 = wall_mean + wall_del * sin( t_c )

      q(i_u,i,   0,k) = -q(i_u,i, 1,k) + 2.d0 * uw1
      q(i_u,i,n2+1,k) = -q(i_u,i,n2,k) + 2.d0 * uw1

    endif
!.... end wall oscillations

!.... slip BC
    if ( slip_bc ) then

      ! x direction
      q(i_u,i,   0,k) = cf_slip(i_u,0,1) * q(i_u,i, 1,k) + &
          cf_slip(i_u,0,2) * q(i_u,i,   2,k)
      q(i_u,i,n2+1,k) = cf_slip(i_u,0,1) * q(i_u,i,n2,k) + &
          cf_slip(i_u,0,2) * q(i_u,i,n2-1,k)

      ! z direction
      q(i_w,i,   0,k) = cf_slip(i_w,0,1) * q(i_w,i, 1,k) + &
          cf_slip(i_w,0,2) * q(i_w,i,   2,k)
      q(i_w,i,n2+1,k) = cf_slip(i_w,0,1) * q(i_w,i,n2,k) + &
            cf_slip(i_w,0,2) * q(i_w,i,n2-1,k)

    endif
!.... end slip BC

!.... pat BC (Kim 09.14.22)
    if ( pat_bc ) then

      if ( phase(i_u,i,k) .ne. 0 ) then
        q(i_u,i,0,k) = q(i_u,i,1,k)
        q(i_u,i,n2+1,k) = q(i_u,i,n2,k)
      endif

      if ( phase(i_w,i,k) .ne. 0 ) then
        q(i_w,i,0,k) = q(i_w,i,1,k)
        q(i_w,i,n2+1,k) = q(i_w,i,n2,k)
      endif

    endif
!.... end pat BC

!.... v bc 
    !q(i_v,i,   0,k) = 0.d0
    if (firstloop ) q(i_v,i,   0,k) = 0.d0      
    q(i_v,i,  -1,k) = 0.d0
    q(i_v,i,  -2,k) = 0.d0
    
    !q(i_v,i,n2  ,k) = 0.d0
    if (firstloop ) q(i_v,i,n2  ,k) = 0.d0
    q(i_v,i,n2+1,k) = 0.d0
    q(i_v,i,n2+2,k) = 0.d0

!.... u,w bc 
!**** end 3rd order taylor series  

    enddo
    enddo

!.... use MPI libraries to communicate the halo regions in each 
!.... domain 
!....

    span_block = nvel* (nx+2*halo_x)*(n2+2*halo_yu)* halo_z
    jss        = 1-halo_yu


!.... perform left<->right communication ; all tags are set to 0 

    if ( p1 > 1 ) then 

       if ( mod( coords(0), 2) .eq. 0 ) then 
          call mpi_sendrecv ( q(i_u,ie-halo_x+1,jss,ks-halo_z), 1, uvec_comm, right, 0, & 
                              q(i_u,ie+1       ,jss,ks-halo_z), 1, uvec_comm, right, 0, comm3d, status, ierr)

          call mpi_sendrecv ( q(i_u,is         ,jss,ks-halo_z), 1, uvec_comm, left , 0, & 
                              q(i_u,is-halo_x  ,jss,ks-halo_z), 1, uvec_comm, left , 0, comm3d, status, ierr) 

       else  

          call mpi_sendrecv ( q(i_u,is         ,jss,ks-halo_z), 1, uvec_comm, left , 0, & 
                              q(i_u,is-halo_x  ,jss,ks-halo_z), 1, uvec_comm, left , 0, comm3d, status, ierr) 

          call mpi_sendrecv ( q(i_u,ie-halo_x+1,jss,ks-halo_z), 1, uvec_comm, right, 0, & 
                              q(i_u,ie+1       ,jss,ks-halo_z), 1, uvec_comm, right, 0, comm3d, status, ierr) 
       endif

    else  

       do k=ks,ke
       do j=1-halo_yu,n2+halo_yu
       do i_q=1,nvel 

          q(i_q,-3,j,k) = q(i_q,n1-3,j,k) 
          q(i_q,-2,j,k) = q(i_q,n1-2,j,k) 
          q(i_q,-1,j,k) = q(i_q,n1-1,j,k) 


          q(i_q,n1,j,k  ) = q(i_q,0,j,k) 
          q(i_q,n1+1,j,k) = q(i_q,1,j,k) 
          q(i_q,n1+2,j,k) = q(i_q,2,j,k) 

       enddo
       enddo
       enddo


    endif
       

!.... perform front<->back communication ; all tags are set to 0 
    
    if ( p3 > 1 ) then 
       if ( mod( coords(2), 2) .eq. 0 ) then 
          call mpi_sendrecv ( q(i_u,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                              q(i_u,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr)
          
          call mpi_sendrecv ( q(i_u,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                              q(i_u,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
       else  
          call mpi_sendrecv ( q(i_u,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                              q(i_u,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
          call mpi_sendrecv ( q(i_u,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                              q(i_u,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
       endif

    else  

       do j=1-halo_yu,n2+halo_yu
       do i=is-halo_x,ie+halo_x
       do i_q=1,nvel

          q(i_q,i,j,-3)   = q(i_q,i,j,n3-3) 
          q(i_q,i,j,-2)   = q(i_q,i,j,n3-2) 
          q(i_q,i,j,-1)   = q(i_q,i,j,n3-1) 
       
          q(i_q,i,j,n3 )  = q(i_q,i,j,0   ) 
          q(i_q,i,j,n3+1) = q(i_q,i,j,1   ) 
          q(i_q,i,j,n3+2) = q(i_q,i,j,2   ) 

    enddo
    enddo
    enddo

    endif


  end subroutine bcuvw

!=================================================================================================================================!

  subroutine bcp (pp) 

!.... periodic boundary conditions in stream / span 


    implicit none 
    integer                                                                                :: i,j,k,i_q, jss
    real, dimension(is-halo_x+1:ie+halo_x+1,1-halo_yp:n2+halo_yp, ks-halo_z+1:ke+halo_z+1) :: pp 
    integer, dimension(mpi_status_size)                                                    :: status 

    do k=ks+1,ke+1 
    do i=is+1,ie+1 

       pp(i,   0,k) =  2.d0*pp(i,   1,k) - pp(i,  2,k) 
       pp(i,  -1,k) =  pp(i,   2,k) 
       
       pp(i,n2+1,k) =  2.d0*pp(i,n2  ,k) - pp(i,n2-1,k) 
       pp(i,n2+2,k) =  pp(i,n2-1,k)    
              
    enddo
    enddo 

!.... perform span/stream communication of the given pressure vector pp 
!.... 
!.... left<->right communication  

    span_block = (nx+2*halo_x)*(n2+2*halo_yp)*halo_z 
    jss        = 1-halo_yp

    if ( p1 > 1 ) then 

       if ( mod(coords(0), 2) .eq. 0 ) then 
          call mpi_sendrecv ( pp(ie-halo_x+2,jss,ks-halo_z+1), 1, pvec_comm, right, 0, & 
                              pp(ie+2       ,jss,ks-halo_z+1), 1, pvec_comm, right, 0, comm3d, status, ierr) 

          call mpi_sendrecv ( pp(is+1       ,jss,ks-halo_z+1), 1, pvec_comm, left , 0, & 
                              pp(is-halo_x+1,jss,ks-halo_z+1), 1, pvec_comm, left , 0, comm3d, status, ierr) 

       else  
          call mpi_sendrecv ( pp(is+1       ,jss,ks-halo_z+1), 1, pvec_comm, left , 0, & 
                              pp(is-halo_x+1,jss,ks-halo_z+1), 1, pvec_comm, left , 0, comm3d, status, ierr) 

          call mpi_sendrecv ( pp(ie-halo_x+2,jss,ks-halo_z+1), 1, pvec_comm, right, 0, & 
                              pp(ie+2       ,jss,ks-halo_z+1), 1, pvec_comm, right, 0, comm3d, status, ierr)
       endif

    else  

       do k= ks+1,ke+1
       do j=1-halo_yp,n2+halo_yp


          pp(0   ,j,k) =  pp(n1  ,j,k) 
          pp(-1  ,j,k) =  pp(n1-1,j,k) 
          pp(n1+1,j,k) =  pp(1   ,j,k) 
          pp(n1+2,j,k) =  pp(2   ,j,k) 

       enddo
       enddo

    endif
       


!.... front<-->back communication

    if ( p3 > 1 ) then 

       if ( mod(coords(2), 2) .eq. 0 ) then 
          call mpi_sendrecv ( pp(is-halo_x+1,jss,ke-halo_z+2), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                              pp(is-halo_x+1,jss,ke+2       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
       
          call mpi_sendrecv ( pp(is-halo_x+1,jss,ks+1       ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                              pp(is-halo_x+1,jss,ks-halo_z+1), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr) 
       else 
          call mpi_sendrecv ( pp(is-halo_x+1,jss,ks+1       ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                              pp(is-halo_x+1,jss,ks-halo_z+1), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr) 
       
          call mpi_sendrecv ( pp(is-halo_x+1,jss,ke-halo_z+2), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                              pp(is-halo_x+1,jss,ke+2       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
       endif

    else  
         
       do j=1-halo_yp,n2+halo_yp
       do i=is-halo_x+1,ie+halo_x+1

          pp(i,j,   0) =  pp(i,j,n3  ) 
          pp(i,j,  -1) =  pp(i,j,n3-1) 
          pp(i,j,n3+1) =  pp(i,j,1   ) 
          pp(i,j,n3+2) =  pp(i,j,2   ) 

       enddo
       enddo

    endif



  end subroutine bcp

!============================================================================================================================!

  
  subroutine compute_modwv 

!.... subroutine to precompute modified wave number for the 
!.... poisson solve 
!.... 
!.... derived from septadiagonal fourth order difference operator 

    implicit none 
    integer       :: i,k 
    real          :: hy0, rhxsq, rhzsq , hx, hz
    
    
    hy0 = 2.d0 / dble(n2) 
    hx  = (chlx*pi) / dble(n1) 
    hz  = (chlz*pi) / dble(n3) 
    
    rhxsq = (1.d0 / hx )**2 
    rhzsq = (1.d0 / hz )**2 


    kxsq(0) = 0.d0 
    kzsq(0) = 0.d0 
    
    do i=1,n1/2-1 
       kxsq(i) = rhxsq * ( 730.d0/288.d0 + 3.d0*cos(4.d0*pi*dble(i)/dble(n1))/16.d0 & 
                          -87.d0/32.d0*cos(2.d0*pi*dble(i)/dble(n1)) - cos(6.d0*pi*dble(i)/dble(n1))/288.d0 )

       kxsq(n1-i) = kxsq(i) 

    enddo 

    i = n1/2 
    kxsq(i) = rhxsq * ( 730.d0/288.d0 + 3.d0*cos(4.d0*pi*dble(i)/dble(n1))/16.d0 & 
                       -87.d0/32.d0*cos(2.d0*pi*dble(i)/dble(n1)) - cos(6.d0*pi*dble(i)/dble(n1))/288.d0 )


    do k=1,n3/2-1 
       kzsq(k) = rhzsq * ( 730.d0/288.d0 + 3.d0*cos(4.d0*pi*dble(k)/dble(n3))/16.d0 & 
                          -87.d0/32.d0*cos(2.d0*pi*dble(k)/dble(n3)) - cos(6.d0*pi*dble(k)/dble(n3))/288.d0 )

       kzsq(n3-k) = kzsq(k) 
    enddo 

    k = n3/2 
    kzsq(k) = rhzsq * ( 730.d0/288.d0 + 3.d0*cos(4.d0*pi*dble(k)/dble(n3))/16.d0 & 
                       -87.d0/32.d0*cos(2.d0*pi*dble(k)/dble(n3)) - cos(6.d0*pi*dble(k)/dble(n3))/288.d0 )


    kxsq = hy0*hy0 * kxsq 
    kzsq = hy0*hy0 * kzsq 



  end subroutine compute_modwv 

!===================================================================================================================================!
  subroutine compute_modwv_2nd
  
!.... subroutine to precompute modified wave number for the 
!.... poisson solve 
!.... 
!.... derived from septadiagonal fourth order difference operator 

    implicit none
    integer       :: i,k
    real          :: hy0, rhxsq, rhzsq , hx, hz


    hy0 = 2.d0 / dble(n2)
    hx  = (chlx*pi) / dble(n1)
    hz  = (chlz*pi) / dble(n3)
    
    rhxsq = (1.d0 / hx )**2
    rhzsq = (1.d0 / hz )**2
    

    kxsq_2nd(0) = 0.d0
    kzsq_2nd(0) = 0.d0
    
    do i=1,n1/2-1
       kxsq_2nd(i) = rhxsq * ( 2.d0 - 2.d0 * cos(2.d0*pi*dble(i)/dble(n1)) )

       kxsq_2nd(n1-i) = kxsq_2nd(i)

    enddo

    i = n1/2
    kxsq_2nd(i) = rhxsq * ( 2.d0 - 2.d0*cos(2.d0*pi*dble(i)/dble(n1)) )
    

    do k=1,n3/2-1
       kzsq_2nd(k) = rhzsq * ( 2.d0 - 2.d0*cos(2.d0*pi*dble(k)/dble(n3)) )

       kzsq_2nd(n3-k) = kzsq_2nd(k)
    enddo
    
    k = n3/2
       kzsq_2nd(k) = rhzsq * ( 2.d0  - 2.d0*cos(2.d0*pi*dble(k)/dble(n3)) )


    kxsq_2nd = hy0*hy0 * kxsq_2nd
    kzsq_2nd = hy0*hy0 * kzsq_2nd



  end subroutine compute_modwv_2nd


!===================================================================================================================================!

  subroutine init_index_2D

!.... initialize in,out for poisson solve indexing 
    implicit none 
    integer        :: count, i,k , j

  
    !------ 
    ! stb edits, 6/18/2010
    !------
    count = 0 
    !.. init in, i-> fast, j->mid, k-< slow
    
    do k=ks+1,ke+1
    do j=1,2
    do i=is+1,ie+1

       count = count + 1 
       in_2D(i,j,k) = count 
    enddo
    enddo
    enddo


    count = 0 
    ! init out i->mid, j-> slow, k-> fast 

    do j=1,2 
    do i=is+1,ie+1
    do k=ks+1,ke+1

       count = count+1 
       out_2D(i,j,k) = count 
    enddo
    enddo
    enddo

    
    count = 0 
    ! init finish_i , i->slow, j->fast, k-> mid
    
    do i=is+1,ie+1
    do k=ks+1,ke+1
    do j=1,2
    
       count = count+1
       finish_2D(i,j,k) = count
    enddo
    enddo
    enddo

  end subroutine init_index_2D
!=====================================================================================================================================! 

  subroutine init_indexing 

!.... initialize in,out for poisson solve indexing 
    implicit none 
    integer        :: count, i,k , j

  
    !... init in, i-> fast, k-> slow 

    !count = 0 
    !do k=ks+1,ke+1 
    !do i=is+1,ie+1 

    !count   = count + 1 
    !   in(i,k) = count 
    !enddo 
    !enddo 


    !.... init out, i-> slow, k-> fast 
    !count = 0 
    !do i=is+1,ie+1 
    !do k=ks+1,ke+1 

    !   count    = count + 1 
    !   out(i,k) = count 
    !enddo 
    !enddo 



    !------ 
    ! stb edits, 6/18/2010
    !------
    count = 0 
    !.. init in, i-> fast, j->mid, k-< slow
    
    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1

       count = count + 1 
       in_new(i,j,k) = count 
    enddo
    enddo
    enddo


    count = 0 
    ! init out i->mid, j-> slow, k-> fast 

    do j=1,n2 
    do i=is+1,ie+1
    do k=ks+1,ke+1

       count = count+1 
       out_new(i,j,k) = count 
    enddo
    enddo
    enddo

    
    count = 0 
    ! init finish_i , i->slow, j->fast, k-> mid
    
    do i=is+1,ie+1
    do k=ks+1,ke+1
    do j=1,n2
    
       count = count+1
       finish_new(i,j,k) = count
    enddo
    enddo
    enddo

  end subroutine init_indexing 
!=====================================================================================================================================! 
  subroutine setup_fft_2D

    implicit none 
    integer        :: ilo_in,ihi_in,klo_in,khi_in
    integer        :: ilo_out,ihi_out,klo_out,khi_out 
    integer        :: jlo_in, jlo_out, jhi_in, jhi_out
    integer        :: iscale, ipermute, nbuf
    integer        :: ipermute_f, ipermute_b 

    !-----
    ! stb edits on 6/18/2010
    !-----
    ilo_in = is+1 
    ihi_in = ie+1 
    jlo_in = 1
    jhi_in = 2
    klo_in = ks+1 
    khi_in = ke+1 
    
    ilo_out= is+1 
    ihi_out= ie+1 
    jlo_out= 1
    jhi_out= 2 
    klo_out= ks+1 
    khi_out= ke+1 
    
    iscale = 0    ! no fft scaling 
    ipermute_f= 2   ! transpose 
    ipermute_b= 2   ! transpose back ..
    

    call fft_2d_3decmp_create_plan ( mycomm, n1, 2, n3, ilo_in, ihi_in, & 
                                     jlo_in, jhi_in, klo_in, khi_in, ilo_out, ihi_out, & 
                                     jlo_out,jhi_out,klo_out,khi_out, iscale, ipermute_f, nbuf, plan_2D) 


    call fft_2d_3decmp_create_plan ( mycomm, n3, n1, 2, klo_out, khi_out, & 
                                     ilo_out, ihi_out, jlo_out, jhi_out, klo_in, khi_in, & 
                                     ilo_in, ihi_in, jlo_in, jhi_in, iscale, ipermute_b, nbuf, ipln_2D ) 


    ! on exit, 
    
    
    
    !.... fourier transform in stream/spanwise directions 

    !call fft_2d_create_plan ( mycomm, n1, n3, ilo_in , ihi_in , & 
    !                          klo_in, khi_in, ilo_out, ihi_out, & 
    !                          klo_out,khi_out, iscale, ipermute, nbuf, plan ) 


    !call fft_2d_create_plan ( mycomm, n3, n1,   klo_out, khi_out, & 
    !                          ilo_out, ihi_out, klo_in , khi_in , & 
    !                          ilo_in , ihi_in , iscale, ipermute, nbuf, ipln ) 


    

  end subroutine setup_fft_2D
!=====================================================================================================================================!


  subroutine setup_fft 

    implicit none 
    integer        :: ilo_in,ihi_in,klo_in,khi_in
    integer        :: ilo_out,ihi_out,klo_out,khi_out 
    integer        :: jlo_in, jlo_out, jhi_in, jhi_out
    integer        :: iscale, ipermute, nbuf
    integer        :: ipermute_f, ipermute_b 

    !-----
    ! stb edits on 6/18/2010
    !-----
    ilo_in = is+1 
    ihi_in = ie+1 
    jlo_in = 1
    jhi_in = n2
    klo_in = ks+1 
    khi_in = ke+1 
    
    ilo_out= is+1 
    ihi_out= ie+1 
    jlo_out= 1
    jhi_out= n2 
    klo_out= ks+1 
    khi_out= ke+1 
    
    iscale = 0    ! no fft scaling 
    ipermute_f= 2   ! transpose 
    ipermute_b= 2   ! transpose back ..
    

    call fft_2d_3decmp_create_plan ( mycomm, n1, n2, n3, ilo_in, ihi_in, & 
                                     jlo_in, jhi_in, klo_in, khi_in, ilo_out, ihi_out, & 
                                     jlo_out,jhi_out,klo_out,khi_out, iscale, ipermute_f, nbuf, plan) 


    call fft_2d_3decmp_create_plan ( mycomm, n3, n1, n2, klo_out, khi_out, & 
                                     ilo_out, ihi_out, jlo_out, jhi_out, klo_in, khi_in, & 
                                     ilo_in, ihi_in, jlo_in, jhi_in, iscale, ipermute_b, nbuf, ipln ) 


    ! on exit, 
    
    
    
    !.... fourier transform in stream/spanwise directions 

    !call fft_2d_create_plan ( mycomm, n1, n3, ilo_in , ihi_in , & 
    !                          klo_in, khi_in, ilo_out, ihi_out, & 
    !                          klo_out,khi_out, iscale, ipermute, nbuf, plan ) 


    !call fft_2d_create_plan ( mycomm, n3, n1,   klo_out, khi_out, & 
    !                          ilo_out, ihi_out, klo_in , khi_in , & 
    !                          ilo_in , ihi_in , iscale, ipermute, nbuf, ipln ) 


    

  end subroutine setup_fft
!=====================================================================================================================================!

  subroutine destroy_fft 

    implicit none 

    !call fft_2d_destroy_plan ( plan) 
    !call fft_2d_destroy_plan ( ipln) 

    call fft_2d_3decmp_destroy_plan ( plan ) 
    call fft_2d_3decmp_destroy_plan ( ipln ) 

    call fft_2d_3decmp_destroy_plan ( plan_2D )
    call fft_2d_3decmp_destroy_plan ( ipln_2D )

    ! misplaced...
    !call dealloc_sim_variables  

  end subroutine destroy_fft 
!=================================================================================================!

end module numerics
