
module numerics_tee

  use global 
  use numerics
  implicit none 

  !include 'fftw_f77.i'
  !real                        :: plan, ipln
  !real                        :: plan_2D, ipln_2D

contains 

!======================================================================================================================================!
!.... Modified by Danah Park on May 15, 2018
!.... Numerical subroutines needed for transportee variable
!===================================================================================================================! 

  subroutine pressure_step_tee(istep) 

    implicit none
    integer, intent(in)                   :: istep
         if (p_order .eq. 2) then
         call p_solve_FFT_2nd_tee( istep )
         elseif (p_order .eq. 4) then
         call p_solve_tee ( istep )
         endif

  end subroutine 

!===================================================================================================================! 

  subroutine rk_step2nd_tee ( istep ) 
    !------------------------------
    ! advance the momentum eqn-
    !     either a full time step for ab2/cn 
    !         or one rk substep for rk3/cn
    !------------------------------

    implicit none 
    integer, intent(in)                   :: istep 
    real, dimension(:,:,:,:), allocatable :: fnl, lun, phi 
    real, dimension(1:2*kl+ku+1,1:n2)     :: wn_u_tee, wn_w_tee, wn_u_phase, wn_w_phase
    real, dimension(1:2*kl+ku+1,1:n2-1)   :: wn_v_tee, wn_v_phase
    integer, dimension(1:n2)              :: ipiv_u_tee, ipiv_w_tee, ipiv_u_phase, ipiv_w_phase
    integer, dimension(1:n2-1)            :: ipiv_v_tee, ipiv_v_phase
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


    fnl =0.d0 
    lun =0.d0 
    phi =0.d0


    !-------------------
    ! form rhs : rhs = u_n + dt*a_i*L(u_n) + dt*g_i*N_n + dt*d_i*N_old 
    ! overwrite q with the rhs
    ! compute convective terms and span/stream diffusion terms
    !-------------------

    call cnv2nd_tee( fnl ) 

    !------ 
    ! add pressure to nonlinear term 
    !------ 

    do k=ks_sf,ke_sf 
    do j=1,n2 
    do i=is_sf,ie_sf 

       fnl(i_u,i,j,k) = fnl(i_u,i,j,k) - pgm_tee 
       ! pgm is set 1 as same as that of transporter
    enddo 
    enddo 
    enddo 


    !------- 
    ! streamwise/ spanwise diffusion 
    !-------
    call dif2nd_tee ( fnl )                 ! fnl N_n

    !-----------------
    ! use pressure from previous time step in momentum advancement 
    ! cf, dukowicz & dvinksky, van ken, askelvoll & moin 
    ! pressure values kept in p_old 
    !-----------------

    !.. using periodic boundary condition for spanwise / stream communication
    call bcp   (       p_old_tee) 


    if (p_order .eq. 2) then
    call gradp2nd (  phi, p_old_tee) 
    else 
    call gradp (  phi, p_old_tee)
    endif
    

    !----------
    ! compute L(u^n) 
    ! for purposes of operator reuse, the filtering of the 
    ! wn visc terms is moved outside of sub.difwn_tee
    !----------

    call difwn_tee  ( lun, is_df, ie_df, ks_df, ke_df, molc_visc, vt_bar ) 

    if ( itimestep .eq. 1 .and. flat_EE) then
       if ( rank .eq. 0) write(*,*) 'using euler for first time step ...'
       !--------------------
       ! for the first time step using ab2, using euler for convective terms
       !-------------------
       do k=ks_sf,ke_sf 
       do j=1,n2 
       do i=is_sf,ie_sf 
       do n=1,nvel 
          f_old_tee(n,i,j,k) = fnl(n,i,j,k)  
          !f_old_source(n,i,j,k) = fnl(n,i,j,k) - s_tee(n,j)
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
          f_old_tee(n,i,j,k) = fnl(n,i,j,k) - lun(n,i,j,k)
          !f_old_source(n,i,j,k) = fnl(n,i,j,k) - lun(n,i,j,k) - s_tee(n,j)
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
          q_tee(n,i,j,k) = q_tee(n,i,j,k) + dt* rkc(istep,1)* lun(n,i,j,k) +  & 
                       +dt* rkc(istep,3)* fnl(n,i,j,k) + dt* rkc(istep,4)* f_old_tee(n,i,j,k) & 
                       -dt*(rkc(istep,1) + rkc(istep,2)) * phi(n,i,j,k) 
       endif 
    enddo
    enddo
    enddo
    enddo


if ( itimestep .eq. 1 .or. itimestep .eq. 2 ) then 
    !--------
    ! solve (I - dt*b_i*Lun ) 
    !--------
    if ( (.not. jsgs) .or. leddyv_exp) then
       call wn_op_tee ( wn_u_tee, ipiv_u_tee, wn_v_tee, ipiv_v_tee, wn_w_tee, ipiv_w_tee, rkc(istep, 2), 1,1 , molc_visc, 0.d0) 
    else 
       call wn_op_tee ( wn_u_tee, ipiv_u_tee, wn_v_tee, ipiv_v_tee, wn_w_tee, ipiv_w_tee, rkc(istep, 2), 1,1 , molc_visc, vt_bar)
    endif 

    ! (Kim 09.14.22)
    if ( pat_bc .and. match_bc ) then
      if ( (.not. jsgs) .or. leddyv_exp) then
         call wn_op_phase ( wn_u_phase, ipiv_u_phase, wn_v_phase, ipiv_v_phase, wn_w_phase, ipiv_w_phase, rkc(istep, 2), 1,1 , molc_visc, 0.d0) 
      else 
         call wn_op_phase ( wn_u_phase, ipiv_u_phase, wn_v_phase, ipiv_v_phase, wn_w_phase, ipiv_w_phase, rkc(istep, 2), 1,1 , molc_visc, vt_bar)
      endif 
    endif
!endif
! Kim 03.12.23: see comment at line 228 (after dgbtrs calls)

    do k=ks,ke
    do i=is,ie
      ! u
      if ( pat_bc .and. match_bc .and. phase(i_u,i,k) .ne. 0 ) then
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_u_phase, 2*kl+ku+1, ipiv_u_phase, q_tee(i_u,i,1:n2  ,k), n2  , ierr)
      else
       call dgbtrs( 'n', n2  , kl,ku, 1, wn_u_tee, 2*kl+ku+1, ipiv_u_tee, q_tee(i_u,i,1:n2  ,k), n2  , ierr)
      endif
      ! v
      if ( pat_bc .and. match_bc .and. phase(i_v,i,k) .ne. 0 ) then
        call dgbtrs( 'n', n2-1, kl,ku, 1, wn_v_phase, 2*kl+ku+1, ipiv_v_phase, q_tee(i_v,i,1:n2-1,k), n2-1, ierr)
      else
        call dgbtrs( 'n', n2-1, kl,ku, 1, wn_v_tee, 2*kl+ku+1, ipiv_v_tee, q_tee(i_v,i,1:n2-1,k), n2-1, ierr)
      endif
      ! w
      if ( pat_bc .and. match_bc .and. phase(i_w,i,k) .ne. 0 ) then
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_w_phase, 2*kl+ku+1, ipiv_w_phase, q_tee(i_w,i,1:n2  ,k), n2  , ierr)
      else
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_w_tee, 2*kl+ku+1, ipiv_w_tee, q_tee(i_w,i,1:n2  ,k), n2  , ierr)
      endif
    enddo
    enddo

!numerics.f90(324): error #6317: An ENDIF occurred without a corresponding IF
!THEN or ELSE statement.
! commented out endif
endif 
! Kim 03.12.23: this is the wrong endif to comment out
! should be the one above the the dgbtrs() calls


    call bcuvw_tee

    !---------
    ! update f_old_tee 
    !---------

 if (itimestep .ne. 3) then 
    do k=ks_sf,ke_sf
    do j=1,n2 
    do i=is_sf,ie_sf
    do n=1,nvel 
       f_old_tee(n,i,j,k) = fnl(n,i,j,k) 
       !f_old_source(n,i,j,k) = fnl(n,i,j,k) - s_tee(n,j)
    enddo
    enddo
    enddo 
    enddo
else 
    do k=ks_sf,ke_sf
    do j=1,n2
    do i=is_sf,ie_sf
    do n=1,nvel
       f_old_tee(n,i,j,k) = fnl(n,i,j,k) - lun(n,i,j,k) 
       !f_old_source(n,i,j,k) = fnl(n,i,j,k) - lun(n,i,j,k) - s_tee(n,j)
    enddo
    enddo
    enddo
    enddo
endif


    deallocate ( fnl  ) 
    deallocate ( lun  ) 
    deallocate ( phi  ) 

  end subroutine rk_step2nd_tee


!===================================================================================================================! 

  subroutine rk_step2nd_IMFM ( istep ) 
    !------------------------------
    ! advance the momentum eqn-
    !     either a full time step for ab2/cn 
    !         or one rk substep for rk3/cn
    !------------------------------

    implicit none 
    integer, intent(in)                   :: istep 
    real, dimension(:,:,:,:), allocatable :: fnl, lun, phi 
    real, dimension(1:2*kl+ku+1,1:n2)     :: wn_u_tee, wn_w_tee, wn_u_phase, wn_w_phase
    real, dimension(1:2*kl+ku+1,1:n2-1)   :: wn_v_tee, wn_v_phase
    real, dimension(1:n2,1:n2)            :: pli_u_tee, pli_w_tee
    real, dimension(:,:), allocatable       :: du
    real, dimension(:), allocatable       :: temp
    real, dimension(:), allocatable       :: temp_sum
    real                                  :: qs_temp, qss_temp, qt_temp, qts_temp !temporary value
    integer, dimension(1:n2)              :: ipiv_u_tee, ipiv_w_tee, ipiv_u_phase, ipiv_w_phase
    integer, dimension(1:n2-1)            :: ipiv_v_tee, ipiv_v_phase
    integer                               :: ierr 
    integer                               :: i,j,k,n,jBF
    real                                  :: deltaBF
    integer                               :: if_start, if_end, kf_start, kf_end
    logical                               :: isfirststep

    real, dimension(5)                      :: coeff 
    real, dimension(3)                      :: coeff_2nd
    real                                    :: hm1, hm2, hp1, hp2
    
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
    
    allocate ( du (1:nvel,1:n2) ) 
    allocate ( temp (1:nvel) ) 
    allocate ( temp_sum (1:nvel) )  

    fnl =0.d0 
    lun =0.d0 
    phi =0.d0

    du =0.d0
    

    !-------------------
    ! form rhs : rhs = u_n + dt*a_i*L(u_n) + dt*g_i*N_n + dt*d_i*N_old 
    ! overwrite q with the rhs
    ! compute convective terms and span/stream diffusion terms
    !-------------------

    call cnv2nd_tee( fnl ) 

    !..modified on 20181008 by Danah
    !..used fluctuation part of the transportee vector instead of v_i
    !..10** correspond to the IC set for beta1_ji**
    !.. since rkc(istep,3)=-1.5, instead of adding source term to fnl, subtract!

    ! BF edit
    if ( IMFMswitchBF ) then
      jBF = ICswitchBF
      deltaBF = y(i_u,jBF+1) - y(i_u,jBF)
    endif

    if ( ICswitch .eq. 1013 ) then
      do k=ks_sf,ke_sf 
      do j=1,n2 
      do i=is_sf,ie_sf 
        fnl(i_w,i,j,k) = fnl(i_w,i,j,k) + (q(i_u,i,j,k)+q(i_u,i+1,j,k)+q(i_u,i,j,k-1)+q(i_u,i+1,j,k-1))/4.d0
      enddo
      enddo
      enddo
    endif
    if ( ICswitch .eq. 1031 ) then
      do k=ks_sf,ke_sf 
      do j=1,n2 
      do i=is_sf,ie_sf 
        fnl(i_u,i,j,k) = fnl(i_u,i,j,k) + (q(i_w,i,j,k)+q(i_w,i,j,k+1)+q(i_w,i-1,j,k)+q(i_w,i-1,j,k-1))/4.d0
      enddo
      enddo
      enddo
    endif
    if ( ICswitch .eq. 1021 ) then
      do k=ks_sf,ke_sf 
      do j=1,n2 
      do i=is_sf,ie_sf 
        fnl(i_u,i,j,k) = fnl(i_u,i,j,k) + ((y(i_v,j)-y(i_u,j))*(q(i_v,i-1,j-1,k)+q(i_v,i,j-1,k))+(y(i_u,j)-y(i_v,j-1))*(q(i_v,i,j,k)+q(i_v,i-1,j,k)))/(y(i_v,j)-y(i_v,j-1))/2.d0
      enddo
      enddo
      enddo
    endif
    if ( ICswitch .eq. 1023 ) then
      do k=ks_sf,ke_sf 
      do j=1,n2 
      do i=is_sf,ie_sf 
        fnl(i_w,i,j,k) = fnl(i_w,i,j,k) + ((y(i_v,j)-y(i_u,j))*(q(i_v,i,j-1,k-1)+q(i_v,i,j-1,k))+(y(i_u,j)-y(i_v,j-1))*(q(i_v,i,j,k)+q(i_v,i,j,k-1)))/(y(i_v,j)-y(i_v,j-1))/2.d0
      enddo
      enddo
      enddo
    endif
    if ( ICswitch .eq. 1012 ) then
      do k=ks_sf,ke_sf 
      do j=1,n2 
      do i=is_sf,ie_sf 
        fnl(i_v,i,j,k) = fnl(i_v,i,j,k) + ((y(i_u,j+1)-y(i_v,j))*(q(i_u,i,j,k)+q(i_u,i+1,j,k))+(y(i_v,j)-y(i_u,j))*(q(i_u,i,j+1,k)+q(i_u,i+1,j+1,k)))/(y(i_u,j+1)-y(i_u,j))/2.d0
      enddo
      enddo
      enddo
    endif
    if ( ICswitch .eq. 1032 ) then
      do k=ks_sf,ke_sf 
      do j=1,n2 
      do i=is_sf,ie_sf 
        fnl(i_v,i,j,k) = fnl(i_v,i,j,k) + ((y(i_u,j+1)-y(i_v,j))*(q(i_w,i,j,k)+q(i_w,i,j,k+1))+(y(i_v,j)-y(i_u,j))*(q(i_w,i,j+1,k)+q(i_w,i,j+1,k+1)))/(y(i_u,j+1)-y(i_u,j))/2.d0
      enddo
      enddo
      enddo
    endif
    !From here, added this source term might violate the continuity
    if ( ICswitch .eq. 1011 ) then
      do k=ks_sf,ke_sf 
      do j=1,n2 
      do i=is_sf,ie_sf 
        fnl(i_u,i,j,k) = fnl(i_u,i,j,k) + q(i_u,i,j,k)
      enddo
      enddo
      enddo
    endif
    if ( ICswitch .eq. 1022 ) then
      do k=ks_sf,ke_sf 
      do j=1,n2 
      do i=is_sf,ie_sf 
        fnl(i_v,i,j,k) = fnl(i_v,i,j,k) + q(i_v,i,j,k)
      enddo
      enddo
      enddo
    endif
    if ( ICswitch .eq. 1033 ) then
      do k=ks_sf,ke_sf 
      do j=1,n2 
      do i=is_sf,ie_sf 
        fnl(i_w,i,j,k) = fnl(i_w,i,j,k) + q(i_w,i,j,k)
      enddo
      enddo
      enddo
    endif

    !------- 
    ! streamwise/ spanwise diffusion 
    !-------
    call dif2nd_tee ( fnl )                 ! fnl N_n

    !-----------------
    ! use pressure from previous time step in momentum advancement 
    ! cf, dukowicz & dvinksky, van ken, askelvoll & moin 
    ! pressure values kept in p_old 
    !-----------------

    !.. using periodic boundary condition for spanwise / stream communication
    call bcp   (       p_old_tee) 


    if (p_order .eq. 2) then
    call gradp2nd (  phi, p_old_tee) 
    else 
    call gradp (  phi, p_old_tee)
    endif
    

    !----------
    ! compute L(u^n) 
    ! for purposes of operator reuse, the filtering of the 
    ! wn visc terms is moved outside of sub.difwn_tee
    !----------

    call difwn_tee  ( lun, is_df, ie_df, ks_df, ke_df, molc_visc, vt_bar ) 

    if ( itimestep .eq. 1 .and. flat_EE) then
       if ( rank .eq. 0) write(*,*) 'using euler for first time step ...'
       !--------------------
       ! for the first time step using ab2, using euler for convective terms
       !-------------------
       do k=ks_sf,ke_sf 
       do j=1,n2 
       do i=is_sf,ie_sf 
       do n=1,nvel 
          f_old_tee(n,i,j,k) = fnl(n,i,j,k)  
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
          f_old_tee(n,i,j,k) = fnl(n,i,j,k) - lun(n,i,j,k)
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
          q_source(n,i,j,k) = q_tee(n,i,j,k) + dt* rkc(istep,1)* lun(n,i,j,k) +  & 
                       +dt* rkc(istep,3)* fnl(n,i,j,k) + dt* rkc(istep,4)* f_old_tee(n,i,j,k) & 
                       -dt*(rkc(istep,1) + rkc(istep,2)) * phi(n,i,j,k) 
       endif 
    enddo
    enddo
    enddo
    enddo

    ! Kim 03.13.23: Calculate wall velocity for patterned BC
    if ( pat_bc .and. match_bc ) then

      call sw_calc

      ! Add effect of wall velocity
      i   = 1
      hp1 = y(i_u,i+1) - y(i_u,i )
      hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )
      call d2_coeff_2nd( hm1, hp1, coeff_2nd)

      do k=ks,ke
      do i=is,ie
      ! BF edit
        if ( phase(i_u,i,k) .eq. 0 ) then
          ! solid (Dirichlet)
          q_source(i_u,i,1,k) = q_source(i_u,i,1,k) + dt * rkc(istep,1) * (1.d0/re) * coeff_2nd(1) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * (sb_new)
          q_source(i_u,i,1,k) = q_source(i_u,i,1,k) - dt * rkc(istep,4) * (1.d0/re) * coeff_2nd(1) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * (sb_old)
        !else
          ! air (Neumann)
          !q_source(i_u,i,1,k) = q_source(i_u,i,1,k) + dt * rkc(istep,1) * (1.d0/re) * coeff_2nd(1) * -2.d0 * ( y(i_u,1) - y(i_u,0) ) * (-1.d0)
          !q_source(i_u,i,1,k) = q_source(i_u,i,1,k) - dt * rkc(istep,4) * (1.d0/re) * coeff_2nd(1) * -2.d0 * ( y(i_u,1) - y(i_u,0) ) * (-1.d0)
        endif
      enddo
      enddo

      i   = n2
      hm1 = y(i_u,i  ) - y(i_u,i-1)
      hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )
      call d2_coeff_2nd( hm1, hp1, coeff_2nd )

      do k=ks,ke
      do i=is,ie
      ! BF edit
        if ( phase(i_u,i,k) .eq. 0 ) then
          q_source(i_u,i,n2,k) = q_source(i_u,i,n2,k) + dt * rkc(istep,1) * (1.d0/re) * coeff_2nd(3) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * (st_new)
          q_source(i_u,i,n2,k) = q_source(i_u,i,n2,k) - dt * rkc(istep,4) * (1.d0/re) * coeff_2nd(3) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * (st_old)
        !else
          !q_source(i_u,i,n2,k) = q_source(i_u,i,n2,k) + dt * rkc(istep,1) * (1.d0/re) * coeff_2nd(3) * +2.d0 * ( y(i_u,n2+1) - y(i_u,n2) ) * (-1.d0)
          !q_source(i_u,i,n2,k) = q_source(i_u,i,n2,k) - dt * rkc(istep,4) * (1.d0/re) * coeff_2nd(3) * +2.d0 * ( y(i_u,n2+1) - y(i_u,n2) ) * (-1.d0)
        endif
      enddo
      enddo

    endif
 
    !Add boundary effect - Danah 2018.07.23

    i   = 1
    hp1 = y(i_u,i+1) - y(i_u,i )
    hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    do k=ks,ke
    do i=is,ie
      q_source(i_u,i,1,k) = q_source(i_u,i,1,k) + dt * rkc(istep,1) * (1.d0/re) * coeff_2nd(1) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * v1_bottom
      q_source(i_u,i,1,k) = q_source(i_u,i,1,k) + dt * rkc(istep,2) * (1.d0/re) * coeff_2nd(1) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * v1_bottom
    enddo
    enddo

    i   = n2
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    do k=ks,ke
    do i=is,ie
      q_source(i_u,i,n2,k) = q_source(i_u,i,n2,k) + dt * rkc(istep,1) * (1.d0/re) * coeff_2nd(3) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * v1_top
      q_source(i_u,i,n2,k) = q_source(i_u,i,n2,k) + dt * rkc(istep,2) * (1.d0/re) * coeff_2nd(3) * ( 1.d0 - cf(i_u,0,1) - cf(i_u,0,2) ) * v1_top
    enddo
    enddo

    call bcuvw_tee
    call bcuvw_source

if ( itimestep .eq. 1 .or. itimestep .eq. 2 ) then 
    !--------
    ! solve (I - dt*b_i*Lun ) 
    !--------
    if ( (.not. jsgs) .or. leddyv_exp) then
       call wn_op_tee ( wn_u_tee, ipiv_u_tee, wn_v_tee, ipiv_v_tee, wn_w_tee, ipiv_w_tee, rkc(istep, 2), 1,1 , molc_visc, 0.d0) 
    else 
       call wn_op_tee ( wn_u_tee, ipiv_u_tee, wn_v_tee, ipiv_v_tee, wn_w_tee, ipiv_w_tee, rkc(istep, 2), 1,1 , molc_visc, vt_bar)
    endif 
      
    
    ! (Kim 09.14.22)
    if ( pat_bc .and. match_bc ) then
      if ( (.not. jsgs) .or. leddyv_exp) then
        call wn_op_phase ( wn_u_phase, ipiv_u_phase, wn_v_phase, ipiv_v_phase, wn_w_phase, ipiv_w_phase, rkc(istep, 2), 1,1 , molc_visc, 0.d0) 
      else 
         call wn_op_phase ( wn_u_phase, ipiv_u_phase, wn_v_phase, ipiv_v_phase, wn_w_phase, ipiv_w_phase, rkc(istep, 2), 1,1 , molc_visc, vt_bar)
      endif 
    endif
!  endif
! Kim 03.12.23: see comment at line 580 (after dgbtrs calls)

    do k=ks,ke
    do i=is,ie
      ! u
      if ( pat_bc .and. match_bc .and. phase(i_u,i,k) .ne. 0 ) then
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_u_phase, 2*kl+ku+1, ipiv_u_phase, q_source(i_u,i,1:n2  ,k), n2  , ierr)
      else
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_u_tee, 2*kl+ku+1, ipiv_u_tee, q_source(i_u,i,1:n2  ,k), n2  , ierr)
      endif
      ! v
      if ( pat_bc .and. match_bc .and. phase(i_v,i,k) .ne. 0 ) then
        call dgbtrs( 'n', n2-1, kl,ku, 1, wn_v_phase, 2*kl+ku+1, ipiv_v_phase, q_source(i_v,i,1:n2-1,k), n2-1, ierr)
      else
        call dgbtrs( 'n', n2-1, kl,ku, 1, wn_v_tee, 2*kl+ku+1, ipiv_v_tee, q_source(i_v,i,1:n2-1,k), n2-1, ierr)
      endif
      ! w
      if ( pat_bc .and. match_bc .and. phase(i_w,i,k) .ne. 0 ) then
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_w_phase, 2*kl+ku+1, ipiv_w_phase, q_source(i_w,i,1:n2  ,k), n2  , ierr)
      else
        call dgbtrs( 'n', n2  , kl,ku, 1, wn_w_tee, 2*kl+ku+1, ipiv_w_tee, q_source(i_w,i,1:n2  ,k), n2  , ierr)
      endif
    enddo
    enddo

endif
! Kim 03.12.23: added endif, commented out the one above dgbtrs() calls


    !call bcuvw_tee
    !call bcuvw_source

    !---------
    ! update f_old_tee 
    !---------

 if (itimestep .ne. 3) then 
    do k=ks_sf,ke_sf
    do j=1,n2 
    do i=is_sf,ie_sf
    do n=1,nvel 
       f_old_tee(n,i,j,k) = fnl(n,i,j,k) 
    enddo
    enddo
    enddo 
    enddo
else 
    do k=ks_sf,ke_sf
    do j=1,n2
    do i=is_sf,ie_sf
    do n=1,nvel
       f_old_tee(n,i,j,k) = fnl(n,i,j,k) - lun(n,i,j,k) 
    enddo
    enddo
    enddo
    enddo
endif

! update sw
sb_old = sb_new
st_old = st_new

!IMFM apply

!Add another ba condition
    call bcuvw_tee
    call bcuvw_source

    !---------
    ! update source term 
    !---------
    
    !average in x and z direction
    do j=1,n2
      temp= 0.d0
      temp_sum= 0.d0
    do n=1,nvel
      do k=ks,ke
        do i=is,ie
          temp(n) = temp(n) + q_source(n,i,j,k)
        enddo
      enddo
    enddo

      call mpi_allreduce( temp, temp_sum, nvel, MPI_DOUBLE_PRECISION, mpi_sum, mycomm , ierr) 

      if ( ICswitch .gt. 1000 .and. ICswitch .le. 1099 ) then
        du(i_u,j) = (- (temp_sum(i_u))/dble(n1)/dble(n3))/dt
        du(i_v,j) = (- (temp_sum(i_v))/dble(n1)/dble(n3))/dt
        du(i_w,j) = (- (temp_sum(i_w))/dble(n1)/dble(n3))/dt
      endif

      if ( ICswitch .eq. 1 ) then
        du(i_u,j) = (y(i_u,j) - (temp_sum(i_u)/dble(n1)/dble(n3)))/dt
        du(i_v,j) = (- (temp_sum(i_v))/dble(n1)/dble(n3))/dt
        du(i_w,j) = (- (temp_sum(i_w))/dble(n1)/dble(n3))/dt
      endif

      if ( ICswitch .eq. 2 ) then
        du(i_u,j) = ( 1.d0 - y(i_u,j) * y(i_u,j) - (temp_sum(i_u)/dble(n1)/dble(n3)))/dt
        du(i_v,j) = (- (temp_sum(i_v))/dble(n1)/dble(n3))/dt
        du(i_w,j) = (- (temp_sum(i_w))/dble(n1)/dble(n3))/dt
      endif

      if ( ICswitch .eq. 31 ) then
        du(i_u,j) = ( y(i_u,j) ** 3 - (temp_sum(i_u)/dble(n1)/dble(n3)))/dt
        du(i_v,j) = (- (temp_sum(i_v))/dble(n1)/dble(n3))/dt
        du(i_w,j) = (- (temp_sum(i_w))/dble(n1)/dble(n3))/dt
      endif

      if ( ICswitch .eq. 32 ) then
        du(i_u,j) = ( ( 1.d0 + y(i_u,j) ) ** 3 - (temp_sum(i_u)/dble(n1)/dble(n3)))/dt
        du(i_v,j) = (- (temp_sum(i_v))/dble(n1)/dble(n3))/dt
        du(i_w,j) = (- (temp_sum(i_w))/dble(n1)/dble(n3))/dt
      endif

      ! BF edit
      if ( IMFMswitchBF ) then
        du(i_u,j) = (- (temp_sum(i_u))/dble(n1)/dble(n3))/dt
        du(i_v,j) = (- (temp_sum(i_v))/dble(n1)/dble(n3))/dt
        du(i_w,j) = (- (temp_sum(i_w))/dble(n1)/dble(n3))/dt
        if ( j .gt. jBF ) then
          du(i_u,j) = ( deltaBF - (temp_sum(i_u)/dble(n1)/dble(n3)) )/dt
        endif
      endif

    enddo

    !update source term
    ! change sbar for patterned BC
    if ( pat_bc .and. match_bc ) then ! all pat_bc should be match_bc too
      s_tee = du

      if ( itimestep .ne. 3 ) then
 
        if ( (.not. jsgs) .or. leddyv_exp) then
          call wn_op_origin_pat ( pli_u_tee, ipiv_u_tee, wn_v_tee, ipiv_v_tee, pli_w_tee, ipiv_w_tee, rkc(istep, 2), 0.d0) 
        else 
          call wn_op_origin_pat ( pli_u_tee, ipiv_u_tee, wn_v_tee, ipiv_v_tee, pli_w_tee, ipiv_w_tee, rkc(istep, 2), vt_bar) 
        endif 
    
        call dgetrs( 'n', n2, 1, pli_u_tee, n2, ipiv_u_tee, s_tee(i_u,1:n2), n2, ierr  )
        call dgetrs( 'n', n2, 1, pli_w_tee, n2, ipiv_w_tee, s_tee(i_w,1:n2), n2, ierr  )
        call dgbmv( 'n', n2-1, n2-1, kl, ku, 1.d0, wn_v_tee(kl+1:2*kl+ku+1,1:n2-1), kl+ku+1, du(i_v,1:n2-1), 1, 0.d0, s_tee(i_v,1:n2-1), 1)
      endif

    else
      if ( (.not. jsgs) .or. leddyv_exp) then
        call wn_op_origin ( wn_u_tee, wn_v_tee, wn_w_tee, rkc(istep, 2), 1,1 , molc_visc, 0.d0) 
      else 
        call wn_op_origin ( wn_u_tee, wn_v_tee, wn_w_tee, rkc(istep, 2), 1,1 , molc_visc, vt_bar)
      endif 
    
      call dgbmv( 'n', n2, n2  , kl, ku, 1.d0, wn_u_tee(kl+1:2*kl+ku+1,1:n2), kl+ku+1, du(i_u,1:n2), 1, 0.d0, s_tee(i_u,1:n2), 1)
      call dgbmv( 'n', n2-1, n2-1, kl, ku, 1.d0, wn_v_tee(kl+1:2*kl+ku+1,1:n2-1), kl+ku+1, du(i_v,1:n2-1), 1, 0.d0, s_tee(i_v,1:n2-1), 1)
      call dgbmv( 'n', n2, n2  , kl, ku, 1.d0, wn_w_tee(kl+1:2*kl+ku+1,1:n2), kl+ku+1, du(i_w,1:n2), 1, 0.d0, s_tee(i_w,1:n2),1) 
    endif

    ! update q_tee
    do k=ks_sf,ke_sf
    do j=1,n2 
    do i=is_sf,ie_sf
    do n=1,nvel 

       if ( n .ne. i_v .or. j .lt. n2 )  then 
          q_tee(n,i,j,k) = q_source(n,i,j,k) + du(n,j)*dt
       endif 
    enddo
    enddo
    enddo
    enddo


    call bcuvw_tee
    call bcuvw_source
    
    deallocate ( fnl  ) 
    deallocate ( lun  ) 
    deallocate ( phi  ) 

    deallocate ( du ) 
    deallocate ( temp ) 
    deallocate ( temp_sum ) 

  end subroutine rk_step2nd_IMFM

!===================================================================================================================================!

  subroutine rk_step_tee ( istep ) 


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
    call cnv4th_tee( fnl )
    else
    call cnv4th_2nd_tee( fnl ) 
    endif
    !------ 
    ! add pressure to nonlinear term , fixed mass flux 
    !------ 
    if ( IMFMswitch ) then

    do k=ks_sf,ke_sf 
    do j=1,n2 
    do i=is_sf,ie_sf 

       fnl(i_u,i,j,k) = fnl(i_u,i,j,k) + s_tee(i_u,j) 
       fnl(i_v,i,j,k) = fnl(i_v,i,j,k) + s_tee(i_v,j) 
       fnl(i_w,i,j,k) = fnl(i_w,i,j,k) + s_tee(i_w,j) 
       !fnl(i_u,i,j,k) = fnl(i_u,i,j,k) - pgm_tee 
       ! pgm is set 1 as same as that of transporter
    enddo 
    enddo 
    enddo 
    else 

    do k=ks_sf,ke_sf 
    do j=1,n2 
    do i=is_sf,ie_sf 

       fnl(i_u,i,j,k) = fnl(i_u,i,j,k) - pgm_tee 
       ! pgm is set 1 as same as that of transporter
    enddo 
    enddo 
    enddo 

    endif

    !------- 
    ! streamwise/ spanwise diffusion 
    !-------
    if (p_order .eq. 4) then 
    call dif4th_tee ( fnl )                 ! fnl N_n
    else
    call dif2nd_tee( fnl) 
    endif
    !-----------------
    ! use pressure from previous time step in momentum advancement 
    ! cf, dukowicz & dvinksky, van ken, askelvoll & moin 
    ! pressure values kept in p_old 
    !-----------------

    call bcp   (       p_old) 


    if (p_order .eq. 2) then
    call gradp2nd (  phi, p_old_tee) 
    else 
    call gradp (  phi, p_old_tee)
    endif
    

    !----------
    ! compute L(u^n) 
    ! for purposes of operator reuse, the filtering of the 
    ! wn visc terms is moved outside of sub.difwn_tee
    !----------

    call difwn_tee  ( lun, is_df, ie_df, ks_df, ke_df, molc_visc, vt_bar ) 

    if ( itimestep .eq. 1 .and. flat_EE) then 
       
       if ( rank .eq. 0) write(*,*) 'using euler for first time step ...'
       !--------------------
       ! for the first time step using ab2, using euler for convective terms
       !-------------------
       do k=ks_sf,ke_sf 
       do j=1,n2 
       do i=is_sf,ie_sf 
       do n=1,nvel 
          f_old_tee(n,i,j,k) = fnl(n,i,j,k) 
          !f_old_source(n,i,j,k) = fnl(n,i,j,k) - s_tee(n,j)
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
          q_tee(n,i,j,k) = q_tee(n,i,j,k) + dt* rkc(istep,1)* lun(n,i,j,k)  & 
                       +dt* rkc(istep,3)* fnl(n,i,j,k) + dt* rkc(istep,4)* f_old_tee(n,i,j,k) & 
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
       call dgbtrs( 'n', n2  , kl,ku, 1, wn_u, 2*kl+ku+1, ipiv_u, q_tee(i_u,i,1:n2  ,k), n2  , ierr)
       call dgbtrs( 'n', n2-1, kl,ku, 1, wn_v, 2*kl+ku+1, ipiv_v, q_tee(i_v,i,1:n2-1,k), n2-1, ierr)
       call dgbtrs( 'n', n2  , kl,ku, 1, wn_w, 2*kl+ku+1, ipiv_w, q_tee(i_w,i,1:n2  ,k), n2  , ierr)
    enddo
    enddo

    call bcuvw_tee

    !---------
    ! update f_old_tee 
    !---------

    do k=ks_sf,ke_sf
    do j=1,n2 
    do i=is_sf,ie_sf
    do n=1,nvel 

       f_old_tee(n,i,j,k) = fnl(n,i,j,k) 
       !f_old_source(n,i,j,k) = fnl(n,i,j,k) 

    enddo
    enddo
    enddo 
    enddo


    deallocate ( fnl  ) 
    deallocate ( lun  ) 
    deallocate ( phi  ) 

  end subroutine rk_step_tee

!===================================================================================================================================!
  subroutine p_solve_tee (istep ) 

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

    p_tee = 0.d0 
    sc= 0.d0 

    rhxsq = ( dble(n1) / (chlx*pi) )**2 
    rhzsq = ( dble(n3) / (chlz*pi) )**2
    hy0   = 2.d0 / dble(n2 ) 


    !------------
    ! compute rhs ( D\hat u ) 
    !------------

    if (p_order .eq. 2) then
    call div2nd_tee ( p_tee(is+1:ie+1, 1:n2, ks+1:ke+1) , maxdiv ) 
    else
    call div4th_tee ( p_tee(is+1:ie+1, 1:n2, ks+1:ke+1) , maxdiv )
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
    
       data(in_new(i,j,k)) = cmplx ( p_tee(i,j,k), 0.d0 ) 
    enddo
    enddo
    enddo


    call fft_2d_3decmp( data, data, 1, plan) 


    ! unpack ... 
    do k=klo_out,khi_out
    do j=jlo_out, jhi_out
    do i=ilo_out, ihi_out 
    
       p_tee(i,j,k)  = real ( data(out_new(i,j,k)))
       sc(i,j,k) = aimag( data(out_new(i,j,k)))
    enddo
    enddo
    enddo
    
       
    !---------------------------------
    ! solve system in the wn direction ,( d^2 - (kx^2 + kz^2)I ) \hat p = rhs 
    ! solve for each kx, kz (modified wave numbers ) 
    !---------------------------------

    p_tee  = hy0 * hy0 * p_tee
    sc = hy0 * hy0 * sc 

    !------------------
    ! pin pressure 
    !------------------

    if ( ilo_out .eq. 1 .and. klo_out .eq. 1 ) then 
       p_tee(1,1,1)  = 0.d0 
       sc(1,1,1) = 0.d0 
    endif 

    do k=klo_out-1,khi_out-1
    do i=ilo_out-1,ihi_out-1
       call p_op   ( pwn, ipiv_p, kxsq(i), kzsq(k) )
       call dgbtrs ( 'n', n2, klp, kup, 1, pwn, 2*klp+kup+1, ipiv_p, p_tee(i+1,1:n2,k+1) , n2, ierr)
       call dgbtrs ( 'n', n2, klp, kup, 1, pwn, 2*klp+kup+1, ipiv_p, sc(i+1,1:n2,k+1), n2, ierr)
    enddo
    enddo




    norm = dble(n1*n3) 


    do k=klo_out,khi_out 
    do j=jlo_out,jhi_out
    do i=ilo_out,ihi_out
    
       data(out_new(i,j,k)) = cmplx( p_tee(i,j,k), sc(i,j,k)) 
    enddo
    enddo
    enddo

    call fft_2d_3decmp ( data, data, -1, ipln ) 


    do k=klo_in, khi_in
    do j=jlo_in, jhi_in
    do i=ilo_in, ihi_in 

       p_tee(i,j,k) = real ( data( finish_new(i,j,k))) / norm 
    enddo
    enddo
    enddo

    call bcp(p_tee)   

    
  end subroutine p_solve_tee

!==========================================================================================================================!

 subroutine update_pressure_tee

    implicit none
    integer                                :: i,j,k,ierr
    real :: tcoeff

    tcoeff = dt


    do k=ks+1,ke+1
    do j=1,n2
    do i=is+1,ie+1

     p_old_tee(i,j,k) = p_tee(i,j,k)/tcoeff + p_old_tee(i,j,k)

    enddo
    enddo
    enddo

    call bcp(p_old_tee)

 
 end subroutine update_pressure_tee

!==========================================================================================================================!
  subroutine p_solve_FFT_2nd_tee (istep )

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
    
    p_tee = 0.d0
    sc= 0.d0

    rhxsq = ( dble(n1) / (chlx*pi) )**2
    rhzsq = ( dble(n3) / (chlz*pi) )**2
    hy0   = 2.d0 / dble(n2 )


    !------------
    ! compute rhs ( D\hat u ) 
    !------------

    if (p_order .eq. 2) then
    call div2nd_tee ( p_tee(is+1:ie+1, 1:n2, ks+1:ke+1) , maxdiv )
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

       data(in_new(i,j,k)) = cmplx ( p_tee(i,j,k), 0.d0 )
    enddo
    enddo
    enddo


    call fft_2d_3decmp( data, data, 1, plan)


    ! unpack ... 
    do k=klo_out,khi_out
    do j=jlo_out, jhi_out
    do i=ilo_out, ihi_out

       p_tee(i,j,k)  = real ( data(out_new(i,j,k)))
       sc(i,j,k) = aimag( data(out_new(i,j,k)))
    enddo
    enddo
    enddo




    !---------------------------------
    ! solve system in the wn direction ,( d^2 - (kx^2 + kz^2)I ) \hat p = rhs 
    ! solve for each kx, kz (modified wave numbers ) 
    !---------------------------------

    p_tee  = hy0 * hy0 * p_tee
    sc = hy0 * hy0 * sc

    !------------------
    ! pin pressure 
    !------------------

    if ( ilo_out .eq. 1 .and. klo_out .eq. 1 ) then
       p_tee(1,1,1)  = 0.d0
       sc(1,1,1) = 0.d0
    endif

    do k=klo_out-1,khi_out-1
    do i=ilo_out-1,ihi_out-1
       call p_op_2nd  ( pwn_2nd, ipiv_p, kxsq_2nd(i), kzsq_2nd(k) )
       call dgbtrs ( 'n', n2, klpp, kupp, 1, pwn_2nd, 2*klpp+kupp+1, ipiv_p, p_tee(i+1,1:n2,k+1) , n2, ierr)
       call dgbtrs ( 'n', n2, klpp, kupp, 1, pwn_2nd, 2*klpp+kupp+1, ipiv_p, sc(i+1,1:n2,k+1), n2, ierr)
    enddo
    enddo

    norm = dble(n1*n3)

    do k=klo_out,khi_out
    do j=jlo_out,jhi_out
    do i=ilo_out,ihi_out

       data(out_new(i,j,k)) = cmplx( p_tee(i,j,k), sc(i,j,k))
    enddo
    enddo
    enddo

    call fft_2d_3decmp ( data, data, -1, ipln )


    do k=klo_in, khi_in
    do j=jlo_in, jhi_in
    do i=ilo_in, ihi_in

       p_tee(i,j,k) = real ( data( finish_new(i,j,k))) / norm
    enddo
    enddo
    enddo


    call bcp(p_tee)


  end subroutine p_solve_FFT_2nd_tee

!====================================
  subroutine prjsln_tee (istep )             

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
    real :: qs_temp, qss_temp
    real :: qt_temp, qts_temp

    allocate ( gphi(1:nvel, is_sf:ie_sf, 1:n2, ks_sf:ke_sf)  ) 

    if (p_order .eq. 2) then
    call gradp2nd ( gphi, p_tee)
    else
    call gradp (  gphi, p_tee)
    endif

    ct = 1.d0 

    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 

       q_tee(i_u,i,j,k) = q_tee(i_u,i,j,k) - ct*gphi(i_u,i,j,k) 
       q_tee(i_w,i,j,k) = q_tee(i_w,i,j,k) - ct*gphi(i_w,i,j,k) 
       
    enddo
    enddo
    enddo


    do k=ks,ke
    do j=1,n2-1 
    do i=is,ie

       q_tee(i_v,i,j,k) = q_tee(i_v,i,j,k) - ct*gphi(i_v,i,j,k) 
    enddo
    enddo
    enddo

   !call bcuvw_tee
    deallocate( gphi) 


  end subroutine prjsln_tee

!==============================================================================================================================! 
 subroutine dif4th_tee ( fnl ) 

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
                          -rchx * ( -q_tee(i_q,i-2,j,k) + 16.d0*q_tee(i_q,i-1,j,k) - 30.d0*q_tee(i_q,i,j,k) & 
                                    -q_tee(i_q,i+2,j,k) + 16.d0*q_tee(i_q,i+1,j,k) ) & 
                          -rchz * ( -q_tee(i_q,i,j,k-2) + 16.d0*q_tee(i_q,i,j,k-1) - 30.d0*q_tee(i_q,i,j,k) & 
                                    -q_tee(i_q,i,j,k+2) + 16.d0*q_tee(i_q,i,j,k+1) )

       endif

    enddo
    enddo
    enddo
    enddo


  end subroutine dif4th_tee

!========================================================================================================================!
 subroutine dif2nd_tee ( fnl )

   !-----------------------------------------
   ! For transportee variable
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
    ! vt_bar is set zero


    c9r16 = 9.d0 / 16.d0
    c1r16 = 1.d0 / 16.d0


    do k=kf_start,kf_end
    do j=1,n2
    do i=if_start,if_end
    do i_q = 1,nvel

       if ( i_q .ne. i_v .or. j .lt. n2 ) then

          fnl(i_q,i,j,k) = fnl(i_q,i,j,k) &
                          -rchx * ( q_tee(i_q,i-1,j,k) - 2.d0*q_tee(i_q,i,j,k)   + q_tee(i_q,i+1,j,k) ) &
                          -rchz * ( q_tee(i_q,i,j,k-1) - 2.d0*q_tee(i_q,i,j,k)   + q_tee(i_q,i,j,k+1) )

       endif

    enddo
    enddo
    enddo
    enddo


  end subroutine dif2nd_tee

!========================================================================================================================!

    subroutine div2nd_tee ( div , maxd )

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

       dudx = rhx * ( q_tee(i_u,i+1,j,k) - q_tee(i_u,i ,j,k) )
       dwdz = rhz * ( q_tee(i_w,i,j,k+1) - q_tee(i_w,i,j,k ) )
       dvdy = rjy * ( q_tee(i_v,i,j  ,k) - q_tee(i_v,i,j-1,k)) 

       div(i+1,j  ,k+1) = dudx  + dwdz  + dvdy

    enddo
    enddo
    enddo

    maxd = maxval ( abs(div) )

  end subroutine div2nd_tee

!========================================================================================================================!

    subroutine div4th_tee ( div , maxd ) 

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


       dudx = c9r8 * rhx * ( q_tee(i_u,i+1,j,k) - q_tee(i_u,i ,j,k) ) - c1r8 * rh3x* ( q_tee(i_u,i+2,j,k) - q_tee(i_u,i-1,j,k) )
       dwdz = c9r8 * rhz * ( q_tee(i_w,i,j,k+1) - q_tee(i_w,i,j,k ) ) - c1r8 * rh3z* ( q_tee(i_w,i,j,k+2) - q_tee(i_w,i,j,k-1) )
       dvdy = c9r8 * rjy * ( q_tee(i_v,i,j  ,k) - q_tee(i_v,i,j-1,k)) - c1r8 * rj3y* ( q_tee(i_v,i,j+1,k) - q_tee(i_v,i,j-2,k) )

     if ( j .eq. 1 ) then
         dvdy = (c9r8/jac  + cf(i_v,-1,1)/24.d0/jac)* q_tee(i_v,i,j,k)/hy0  &
               +(-c1r24/jac+ cf(i_v,-1,2)/24.d0/jac)* q_tee(i_v,i,j+1,k)/hy0     !rjy* q_tee(i_v,i,j,k) 
      elseif (j .eq. n2) then
         dvdy =-(c9r8/jac  + cf(i_v,-1,1)/24.d0/jac)* q_tee(i_v,i,j-1,k)/hy0  &
               +(c1r24/jac - cf(i_v,-1,2)/24.d0/jac)* q_tee(i_v,i,j-2,k)/hy0     !-rjy* q_tee(i_v,i,j-1,k) 
      else 
        dvdy = c9r8 * rjy * ( q_tee(i_v,i,j  ,k) - q_tee(i_v,i,j-1,k)) - c1r8 * rj3y* ( q_tee(i_v,i,j+1,k) - q_tee(i_v,i,j-2,k) )
      endif


      div(i+1,j  ,k+1) = dudx  + dwdz  + dvdy 
       
    enddo
    enddo
    enddo

    maxd = maxval ( abs(div) ) 


  end subroutine div4th_tee

   !===============================================================================================================================! 

   subroutine cnv2nd_tee( fnl )
   
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
fdiv(i_u,i,j,k) = rhx* ( (c1r2*( q(i_u,i  ,j,k) + q(i_u,i+1,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i+1,j,k))  &
- (c1r2*( q(i_u,i  ,j,k) + q(i_u,i-1,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i-1,j,k)) )


!... fdiv(u) = ... + d/dz(uw )  
fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rhz* ( (c1r2*( q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k+1))  &
- (c1r2*( q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k-1)) )


!.... fdiv(u) = + .... d/dy( uv) + ... 
fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j-1,k))  )


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
fdiv(i_v,i,j,k) = rhx* ( (c1r2*( q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i+1,j,k))  &
- (c1r2*( q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i-1,j,k)) )

!... fdiv(v) = ... + d/dz(vw) 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rhz* ( (c1r2*( q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k+1))  &
- (c1r2*( q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k-1)) )

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j-1,k)) )

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
fdiv(i_w,i,j,k) = rhx* ( (c1r2*( q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i+1,j,k))  &
- (c1r2*( q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i-1,j,k)) )


!.... fdiv(w) = ... + d/dz(ww)       
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rhz* ( (c1r2*( q(i_w,i,j,k  ) + q(i_w,i,j,k+1)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k+1))  &
- (c1r2*( q(i_w,i,j,k  ) + q(i_w,i,j,k-1)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k-1))) 

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j-1,k)) )

enddo
enddo
enddo


!--------------
! advective form 
!--------------

!do k=kf_start, kf_end
!do j=1,n2
!do i=if_start, if_end
!
!rjy = 1.d0 / jacb(i_u,j)
!
!!.... f(i_u) = u du/dx + .... 
!fadv(i_u,i,j,k) = c1r2* ( (c1r2* (q(i_u,i  ,j,k) + q(i_u,i-1,j,k)))* rhx*   (q_tee(i_u,i  ,j,k) - q_tee(i_u,i-1,j,k)) &
!+(c1r2* (q(i_u,i  ,j,k) + q(i_u,i+1,j,k)))* rhx*   (q_tee(i_u,i+1,j,k) - q_tee(i_u,i  ,j,k)) ) 
!
!!.... f(i_u) =  ... + w du/dz 
!fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + &
! c1r2* ( (c1r2* (q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)))* rhz*   (q_tee(i_u,i  ,j,k+1) - q_tee(i_u,i  ,j,k  )) &
!+(c1r2* (q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )))* rhz*   (q_tee(i_u,i  ,j,k  ) - q_tee(i_u,i  ,j,k-1)) ) 
!
!!.... f(i_u) =  + ... v du/dy + ... 
!fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
!c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q_tee(i_u,i  ,j+1,k) - q_tee(i_u,i  ,j  ,k)) &
!+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q_tee(i_u,i  ,j  ,k) - q_tee(i_u,i  ,j-1,k)) ))
!
!!---- 
!
!enddo
!enddo
!enddo
!
!
!!.... fadv(i_v) 
!
!
!do k=kf_start, kf_end
!do j=1,n2-1
!do i=if_start, if_end
!
!rjy = 1.d0 / jacb(i_v,j)
!
!!.... f(i_v) =  u dv/dx 
!fadv(i_v,i,j,k) =  c1r2* ( (c1r2* (q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)))* rhx*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i-1,j  ,k)) &
!+(c1r2* (q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)))* rhx*   (q_tee(i_v,i+1,j  ,k) - q_tee(i_v,i  ,j  ,k)) ) 
!
!!.... f(i_v) = ... + w dv/dz 
!fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + &
!c1r2* ( (c1r2* (q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )))* rhz*   (q_tee(i_v,i,j  ,k  ) - q_tee(i_v,i,j  ,k-1)) &
!+(c1r2* (q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)))* rhz*   (q_tee(i_v,i,j  ,k+1) - q_tee(i_v,i,j  ,k  )) ) 
!
!
!!.... f(i_v) = ... + v dv/dy + .... 
!fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
!c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q_tee(i_v,i  ,j+1,k) - q_tee(i_v,i  ,j  ,k)) &
!+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i  ,j-1,k)) ) )
!
!enddo
!enddo
!enddo
!
!
!!.... fadv(i_w) 
!
!
!do k=kf_start, kf_end
!do j=1,n2
!do i=if_start, if_end
!
!rjy = 1.d0 / jacb(i_u,j)
!
!!.... fadv(i_w) = u dw/dx + ....
!fadv(i_w,i,j,k) =  c1r2* ( (c1r2* (q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)))* rhx*   (q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i-1,j,k  )) &
!+(c1r2* (q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)))* rhx*   (q_tee(i_w,i+1,j,k  ) - q_tee(i_w,i  ,j,k  )) ) 
!
!
!!.... fadv(i_w) = ... + w dw/dz 
!fadv(i_w,i,j,k) = fadv(i_w,i,j,k) + &
! c1r2* ( (c1r2* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-1)))*rhz*(q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i  ,j,k-1)) &
!+(c1r2* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+1)))*rhz*(q_tee(i_w,i  ,j,k+1) - q_tee(i_w,i  ,j,k  )) ) 
!
!
!!.... fadv(i_w) = ...+ v dw/dy + ... 
!fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
!c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q_tee(i_w,i,j+1,k  ) - q_tee(i_w,i,j  ,k  )) &
!+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q_tee(i_w,i,j  ,k  ) - q_tee(i_w,i,j-1,k  )) ) )
!enddo
!enddo
!enddo
!



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


  end subroutine cnv2nd_tee

   !===============================================================================================================================! 

   subroutine cnv4th_2nd_tee( fnl )
   
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
-c1r16*( q(i_u,i-1,j,k) + q(i_u,i+2,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i+1,j,k))  &
- (c9r16*( q(i_u,i  ,j,k) + q(i_u,i-1,j,k)) &
-c1r16*( q(i_u,i+1,j,k) + q(i_u,i-2,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i-1,j,k)) )&
-c1r8*rh3x*( (c9r16*( q(i_u,i+1,j,k) + q(i_u,i+2,j,k)) &
-c1r16*( q(i_u,i  ,j,k) + q(i_u,i+3,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i+3,j,k))  &
- (c9r16*( q(i_u,i-1,j,k) + q(i_u,i-2,j,k)) &
-c1r16*( q(i_u,i  ,j,k) + q(i_u,i-3,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i-3,j,k)) )


!... fdiv(u) = ... + d/dz(uw )  
fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
c9r8*rhz* ( (c9r16*( q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)) &
-c1r16*( q(i_w,i+1,j,k+1) + q(i_w,i-2,j,k+1)))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k+1))  &
- (c9r16*( q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )) &
-c1r16*( q(i_w,i+1,j,k  ) + q(i_w,i-2,j,k  )))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k-1)) )&
-c1r8*rh3z*( (c9r16*( q(i_w,i  ,j,k+2) + q(i_w,i-1,j,k+2)) &
-c1r16*( q(i_w,i+1,j,k+2) + q(i_w,i-2,j,k+2)))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k+3))  &
- (c9r16*( q(i_w,i  ,j,k-1) + q(i_w,i-1,j,k-1)) &
-c1r16*( q(i_w,i+1,j,k-1) + q(i_w,i-2,j,k-1)))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k-3)) )


if ( j .eq. 1 ) then

!.... fdiv(u) = + .... d/dy( uv) + ... 

fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j-1,k))  )

elseif ( j .eq. 2 ) then

!.... fdiv(u) = + .... d/dy( uv) + ... 

fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j-1,k))  )

elseif (j .eq. n2-1) then

!.... fdiv(u) = + .... d/dy( uv) + ... 

fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j-1,k))  )

elseif (j .eq. n2) then

!.... fdiv(u) = + .... d/dy( uv) + ... 

fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
rjy* ( ( c1r2*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j+1,k))  &
-( c1r2*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i  ,j-1,k))  )

else

! ----
! Interiror: high order scheme (4th order)
! ----

!.... fdiv(u) = + .... d/dy( uv) + ... 
fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
c9r8*rjy* ( (c9r16*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)) &
-c1r16*( q(i_v,i+1,j  ,k) + q(i_v,i-2,j  ,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i,j+1,k))  &
- (c9r16*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)) &
-c1r16*( q(i_v,i+1,j-1,k) + q(i_v,i-2,j-1,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i,j-1,k)) )&
-c1r8*rj3y*( (c9r16*( q(i_v,i  ,j+1,k) + q(i_v,i-1,j+1,k)) &
-c1r16*( q(i_v,i+1,j+1,k) + q(i_v,i-2,j+1,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i,j+3,k))  &
- (c9r16*( q(i_v,i  ,j-2,k) + q(i_v,i-1,j-2,k)) &
-c1r16*( q(i_v,i+1,j-2,k) + q(i_v,i-2,j-2,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i,j-3,k)) )



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
-c1r16*( q(i_u,i+1,j-1,k) + q(i_u,i+1,j+2,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i+1,j,k))  &
- (c9r16*( q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)) &
-c1r16*( q(i_u,i  ,j-1,k) + q(i_u,i  ,j+2,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i-1,j,k)) )&
-c1r8*rh3x*( (c9r16*( q(i_u,i+2,j  ,k) + q(i_u,i+2,j+1,k)) &
-c1r16*( q(i_u,i+2,j-1,k) + q(i_u,i+2,j+2,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i+3,j,k))  &
- (c9r16*( q(i_u,i-1,j  ,k) + q(i_u,i-1,j+1,k)) &
-c1r16*( q(i_u,i-1,j-1,k) + q(i_u,i-1,j+2,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i-3,j,k)) )


!... fdiv(v) = ... + d/dz(vw) 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
c9r8*rhz* ( (c9r16*( q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)) &
-c1r16*( q(i_w,i,j-1,k+1) + q(i_w,i,j+2,k+1)))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k+1))  &
- (c9r16*( q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )) &
-c1r16*( q(i_w,i,j-1,k  ) + q(i_w,i,j+2,k  )))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k-1)) )&
-c1r8*rh3z*( (c9r16*( q(i_w,i,j  ,k+2) + q(i_w,i,j+1,k+2)) &
-c1r16*( q(i_w,i,j-1,k+2) + q(i_w,i,j+2,k+2)))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k+3))  &
- (c9r16*( q(i_w,i,j  ,k-1) + q(i_w,i,j+1,k-1)) &
-c1r16*( q(i_w,i,j-1,k-1) + q(i_w,i,j+2,k-1)))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k-3)) )


! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----

if ( j .eq. 1 ) then

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j-1,k)) )

elseif ( j .eq. 2 ) then

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j-1,k)) )

elseif (j .eq. n2-2) then

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j-1,k)) )

elseif (j .eq. n2-1) then

!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j+1,k))  &
- (c1r2*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j-1,k)) )

else


!.... fdiv(v) = ... + d/dy(vv) + .... 
fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
c9r8*rjy* ( (c9r16*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)) &
-c1r16*( q(i_v,i,j-1,k) + q(i_v,i,j+2,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j+1,k))  &
- (c9r16*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)) &
-c1r16*( q(i_v,i,j+1,k) + q(i_v,i,j-2,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j-1,k)) )&
-c1r8*rj3y*( (c9r16*( q(i_v,i,j+1,k) + q(i_v,i,j+2,k)) &
-c1r16*( q(i_v,i,j  ,k) + q(i_v,i,j+3,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j+3,k))  &
- (c9r16*( q(i_v,i,j-1,k) + q(i_v,i,j-2,k)) &
-c1r16*( q(i_v,i,j  ,k) + q(i_v,i,j-3,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j-3,k)) )

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
-c1r16*( q(i_u,i+1,j,k+1) + q(i_u,i+1,j,k-2)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i+1,j,k))  &
- (c9r16*( q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)) &
-c1r16*( q(i_u,i  ,j,k+1) + q(i_u,i  ,j,k-2)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i-1,j,k)) )&
-c1r8*rh3x*( (c9r16*( q(i_u,i+2,j,k  ) + q(i_u,i+2,j,k-1)) &
-c1r16*( q(i_u,i+2,j,k+1) + q(i_u,i+2,j,k-2)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i+3,j,k))  &
- (c9r16*( q(i_u,i-1,j,k  ) + q(i_u,i-1,j,k-1)) &
-c1r16*( q(i_u,i-1,j,k+1) + q(i_u,i-1,j,k-2)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i-3,j,k)) )

!.... fdiv(w) = ... + d/dz(ww)       
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
c9r8*rhz* ( (c9r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k+1)) &
-c1r16*( q(i_w,i,j,k-1) + q(i_w,i,j,k+2)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k+1))  &
- (c9r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k-1)) &
-c1r16*( q(i_w,i,j,k+1) + q(i_w,i,j,k-2)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k-1))) &
-c1r8*rh3z*( (c9r16*( q(i_w,i,j,k+1) + q(i_w,i,j,k+2)) &
-c1r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k+3)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k+3))  &
- (c9r16*( q(i_w,i,j,k-1) + q(i_w,i,j,k-2)) &
-c1r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k-3)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k-3)))



! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----

if ( j .eq. 1 ) then

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j-1,k)) )

elseif ( j .eq. 2 ) then

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q(i_w,i,j  ,k)+q_tee(i_w,i,j-1,k)) )

elseif (j .eq. n2-1) then

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j-1,k)) )

elseif (j .eq. n2) then

!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
rjy* ( (c1r2*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) )* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j+1,k))  &
- (c1r2*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) )* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j-1,k)) )

else


!.... fdiv(w) = ... + d/dy(wv) + .... 
fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
c9r8*rjy* ( (c9r16*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) &
-c1r16*( q(i_v,i,j  ,k+1) + q(i_v,i,j  ,k-2)))* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j+1,k))  &
- (c9r16*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) &
-c1r16*( q(i_v,i,j-1,k+1) + q(i_v,i,j-1,k-2)))* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j-1,k)) )&
-c1r8*rj3y*( (c9r16*( q(i_v,i,j+1,k  ) + q(i_v,i,j+1,k-1)) &
-c1r16*( q(i_v,i,j+1,k+1) + q(i_v,i,j+1,k-2)))* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j+3,k))  &
- (c9r16*( q(i_v,i,j-2,k  ) + q(i_v,i,j-2,k-1)) &
-c1r16*( q(i_v,i,j-2,k+1) + q(i_v,i,j-2,k-2)))* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j-3,k)) )
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
rhx*   (q_tee(i_u,i  ,j,k) - q_tee(i_u,i-1,j,k)) &
+(c9r16* (q(i_u,i  ,j,k) + q(i_u,i+1,j,k)) - c1r16*(q(i_u,i-1,j,k) + q(i_u,i+2,j,k)))* &
rhx*   (q_tee(i_u,i+1,j,k) - q_tee(i_u,i  ,j,k)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_u,i-2,j,k) + q(i_u,i-1,j,k)) - c1r16*(q(i_u,i-3,j,k) + q(i_u,i  ,j,k)))* &
rh3x*  (q_tee(i_u,i  ,j,k) - q_tee(i_u,i-3,j,k)) &
+(c9r16* (q(i_u,i+2,j,k) + q(i_u,i+1,j,k)) - c1r16*(q(i_u,i+3,j,k) + q(i_u,i  ,j,k)))* &
rh3x*  (q_tee(i_u,i+3,j,k) - q_tee(i_u,i  ,j,k)) )

!.... f(i_u) =  ... + w du/dz 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + &
c9r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)) - c1r16*(q(i_w,i-2,j,k+1) + q(i_w,i+1,j,k+1)))* &
rhz*   (q_tee(i_u,i  ,j,k+1) - q_tee(i_u,i  ,j,k  )) &
+(c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )) - c1r16*(q(i_w,i-2,j,k  ) + q(i_w,i+1,j,k  )))* &
rhz*   (q_tee(i_u,i  ,j,k  ) - q_tee(i_u,i  ,j,k-1)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k-1) + q(i_w,i-1,j,k-1)) - c1r16*(q(i_w,i-2,j,k-1) + q(i_w,i+1,j,k-1)))* &
rh3z*  (q_tee(i_u,i  ,j,k  ) - q_tee(i_u,i  ,j,k-3)) &
+(c9r16* (q(i_w,i  ,j,k+2) + q(i_w,i-1,j,k+2)) - c1r16*(q(i_w,i-2,j,k+2) + q(i_w,i+1,j,k+2)))* &
rh3z*  (q_tee(i_u,i  ,j,k+3) - q_tee(i_u,i  ,j,k  )) )



! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----


if ( j .eq. 1 ) then

!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q_tee(i_u,i  ,j+1,k) - q_tee(i_u,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q_tee(i_u,i  ,j  ,k) - q_tee(i_u,i  ,j-1,k)) ))

elseif ( j .eq. 2 ) then

!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q_tee(i_u,i  ,j+1,k) - q_tee(i_u,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q_tee(i_u,i  ,j  ,k) - q_tee(i_u,i  ,j-1,k)) ))

elseif ( j .eq. n2-1 ) then

!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q_tee(i_u,i  ,j+1,k) - q_tee(i_u,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q_tee(i_u,i  ,j  ,k) - q_tee(i_u,i  ,j-1,k)) ))

elseif ( j .eq. n2 ) then

!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c1r2*(  (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)))* rhy * (q_tee(i_u,i  ,j+1,k) - q_tee(i_u,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)))* rhy * (q_tee(i_u,i  ,j  ,k) - q_tee(i_u,i  ,j-1,k)) ))

else


!.... f(i_u) =  + ... v du/dy + ... 
fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
c9r8 * c1r2*(  (c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)) - c1r16*(q(i_v,i-2,j  ,k) + q(i_v,i+1,j  ,k)))* &
rhy*   (q_tee(i_u,i  ,j+1,k) - q_tee(i_u,i  ,j  ,k)) &
+(c9r16* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)) - c1r16*(q(i_v,i-2,j-1,k) + q(i_v,i+1,j-1,k)))* &
rhy*   (q_tee(i_u,i  ,j  ,k) - q_tee(i_u,i  ,j-1,k)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_v,i  ,j-2,k) + q(i_v,i-1,j-2,k)) - c1r16*(q(i_v,i-2,j-2,k) + q(i_v,i+1,j-2,k)))* &
rh3y*  (q_tee(i_u,i  ,j  ,k) - q_tee(i_u,i  ,j-3,k)) &
+(c9r16* (q(i_v,i  ,j+1,k) + q(i_v,i-1,j+1,k)) - c1r16*(q(i_v,i-2,j+1,k) + q(i_v,i+1,j+1,k)))* &
rh3y*  (q_tee(i_u,i  ,j+3,k) - q_tee(i_u,i  ,j  ,k)) ) )

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
rhx*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i-1,j  ,k)) &
+(c9r16* (q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)) - c1r16*(q(i_u,i+1,j-1,k) + q(i_u,i+1,j+2,k)))* &
rhx*   (q_tee(i_v,i+1,j  ,k) - q_tee(i_v,i  ,j  ,k)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_u,i-1,j  ,k) + q(i_u,i-1,j+1,k)) - c1r16*(q(i_u,i-1,j-1,k) + q(i_u,i-1,j+2,k)))* &
rh3x*  (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i-3,j  ,k)) &
+(c9r16* (q(i_u,i+2,j  ,k) + q(i_u,i+2,j+1,k)) - c1r16*(q(i_u,i+2,j-1,k) + q(i_u,i+2,j+2,k)))* &
rh3x*  (q_tee(i_v,i+3,j  ,k) - q_tee(i_v,i  ,j  ,k)) )

!.... f(i_v) = ... + w dv/dz 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + &
c9r8 * c1r2* ( (c9r16* (q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )) - c1r16*(q(i_w,i,j-1,k  ) + q(i_w,i,j+2,k  )))* &
rhz*   (q_tee(i_v,i,j  ,k  ) - q_tee(i_v,i,j  ,k-1)) &
+(c9r16* (q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)) - c1r16*(q(i_w,i,j-1,k+1) + q(i_w,i,j+2,k+1)))* &
rhz*   (q_tee(i_v,i,j  ,k+1) - q_tee(i_v,i,j  ,k  )) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_w,i,j  ,k-1) + q(i_w,i,j+1,k-1)) - c1r16*(q(i_w,i,j-1,k-1) + q(i_w,i,j+2,k-1)))* &
rh3z*  (q_tee(i_v,i,j  ,k  ) - q_tee(i_v,i,j  ,k-3)) &
+(c9r16* (q(i_w,i,j  ,k+2) + q(i_w,i,j+1,k+2)) - c1r16*(q(i_w,i,j-1,k+2) + q(i_w,i,j+2,k+2)))* &
rh3z*  (q_tee(i_v,i,j  ,k+3) - q_tee(i_v,i,j  ,k  )) )


! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----


if ( j .eq. 1 ) then

!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q_tee(i_v,i  ,j+1,k) - q_tee(i_v,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i  ,j-1,k)) ) )

elseif ( j .eq. 2 ) then

!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q_tee(i_v,i  ,j+1,k) - q_tee(i_v,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i  ,j-1,k)) ) )

elseif ( j .eq. n2-2 ) then

!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q_tee(i_v,i  ,j+1,k) - q_tee(i_v,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i  ,j-1,k)) ) )


elseif ( j .eq. n2-1 ) then

!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c1r2* ( (c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)))* rhy*   (q_tee(i_v,i  ,j+1,k) - q_tee(i_v,i  ,j  ,k)) &
+(c1r2* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)))* rhy*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i  ,j-1,k)) ) )

else


!.... f(i_v) = ... + v dv/dy + .... 
fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
c9r8 * c1r2* ( (c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)) - c1r16*(q(i_v,i  ,j-1,k) + q(i_v,i  ,j+2,k)))* &
rhy*   (q_tee(i_v,i  ,j+1,k) - q_tee(i_v,i  ,j  ,k)) &
+(c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)) - c1r16*(q(i_v,i  ,j+1,k) + q(i_v,i  ,j-2,k)))* &
rhy*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i  ,j-1,k)) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_v,i  ,j+1,k) + q(i_v,i  ,j+2,k)) - c1r16*(q(i_v,i  ,j  ,k) + q(i_v,i  ,j+3,k)))* &
rh3y*  (q_tee(i_v,i  ,j+3,k) - q_tee(i_v,i  ,j  ,k)) &
+(c9r16* (q(i_v,i  ,j-1,k) + q(i_v,i  ,j-2,k)) - c1r16*(q(i_v,i  ,j  ,k) + q(i_v,i  ,j-3,k)))* &
rh3y*  (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i  ,j-3,k)) ) )

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
rhx*   (q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i-1,j,k  )) &
+(c9r16* (q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)) - c1r16*(q(i_u,i+1,j,k-2) + q(i_u,i+1,j,k+1)))* &
rhx*   (q_tee(i_w,i+1,j,k  ) - q_tee(i_w,i  ,j,k  )) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_u,i-1,j,k  ) + q(i_u,i-1,j,k-1)) - c1r16*(q(i_u,i-1,j,k-2) + q(i_u,i-1,j,k+1)))* &
rh3x*  (q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i-3,j,k  )) &
+(c9r16* (q(i_u,i+2,j,k  ) + q(i_u,i+2,j,k-1)) - c1r16*(q(i_u,i+2,j,k-2) + q(i_u,i+2,j,k+1)))* &
rh3x*  (q_tee(i_w,i+3,j,k  ) - q_tee(i_w,i  ,j,k  )) )

!.... fadv(i_w) = ... + w dw/dz 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) + &
c9r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-1)) - c1r16*(q(i_w,i  ,j,k+1) + q(i_w,i  ,j,k-2)))* &
rhz*   (q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i  ,j,k-1)) &
+(c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+1)) - c1r16*(q(i_w,i  ,j,k-1) + q(i_w,i  ,j,k+2)))* &
rhz*   (q_tee(i_w,i  ,j,k+1) - q_tee(i_w,i  ,j,k  )) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k-1) + q(i_w,i  ,j,k-2)) - c1r16*(q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-3)))* &
rh3z*  (q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i  ,j,k-3)) &
+(c9r16* (q(i_w,i  ,j,k+1) + q(i_w,i  ,j,k+2)) - c1r16*(q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+3)))* &
rh3z*  (q_tee(i_w,i  ,j,k+3) - q_tee(i_w,i  ,j,k  )) )



! ----
! near wall lower order scheme (2nd order) edited by Jongmin Seo 2012.06.05
! only for wall normal convection term (i.e.  d/dy terms)
! ----


if ( j .eq. 1 ) then

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q_tee(i_w,i,j+1,k  ) - q_tee(i_w,i,j  ,k  )) &
+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q_tee(i_w,i,j  ,k  ) - q_tee(i_w,i,j-1,k  )) ) )

elseif ( j .eq. 2 ) then

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q_tee(i_w,i,j+1,k  ) - q_tee(i_w,i,j  ,k  )) &
+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q_tee(i_w,i,j  ,k  ) - q_tee(i_w,i,j-1,k  )) ) )

elseif ( j .eq. n2-1 ) then

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q_tee(i_w,i,j+1,k  ) - q_tee(i_w,i,j  ,k  )) &
+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q_tee(i_w,i,j  ,k  ) - q_tee(i_w,i,j-1,k  )) ) )

elseif ( j .eq. n2 ) then

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c1r2* ( (c1r2* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)))*  rhy*   (q_tee(i_w,i,j+1,k  ) - q_tee(i_w,i,j  ,k  )) &
+(c1r2* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)))*  rhy*   (q_tee(i_w,i,j  ,k  ) - q_tee(i_w,i,j-1,k  )) ) )
else

!.... fadv(i_w) = ...+ v dw/dy + ... 
fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
c9r8 * c1r2* ( (c9r16* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) - c1r16*(q(i_v,i,j  ,k+1) + q(i_v,i,j  ,k-2)))* &
rhy*   (q_tee(i_w,i,j+1,k  ) - q_tee(i_w,i,j  ,k  )) &
+(c9r16* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) - c1r16*(q(i_v,i,j-1,k+1) + q(i_v,i,j-1,k-2)))* &
rhy*   (q_tee(i_w,i,j  ,k  ) - q_tee(i_w,i,j-1,k  )) ) &
-c1r8 * c1r2* ( (c9r16* (q(i_v,i,j-2,k  ) + q(i_v,i,j-2,k-1)) - c1r16*(q(i_v,i,j-2,k+1) + q(i_v,i,j-2,k-2)))* &
rh3y*  (q_tee(i_w,i,j  ,k  ) - q_tee(i_w,i,j-3,k  )) &
+(c9r16* (q(i_v,i,j+1,k  ) + q(i_v,i,j+1,k-1)) - c1r16*(q(i_v,i,j+1,k+1) + q(i_v,i,j+1,k-2)))* &
rh3y*  (q_tee(i_w,i,j+3,k  ) - q_tee(i_w,i,j  ,k  )) ) )

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


  end subroutine cnv4th_2nd_tee


   !===============================================================================================================================! 

  subroutine cnv4th_tee ( fnl ) 

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
                                     -c1r16*( q(i_u,i-1,j,k) + q(i_u,i+2,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i+1,j,k))  &
                                   - (c9r16*( q(i_u,i  ,j,k) + q(i_u,i-1,j,k)) &
                                     -c1r16*( q(i_u,i+1,j,k) + q(i_u,i-2,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i-1,j,k)) )&
                        -c1r8*rh3x*( (c9r16*( q(i_u,i+1,j,k) + q(i_u,i+2,j,k)) &
                                     -c1r16*( q(i_u,i  ,j,k) + q(i_u,i+3,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i+3,j,k))  &
                                   - (c9r16*( q(i_u,i-1,j,k) + q(i_u,i-2,j,k)) &
                                     -c1r16*( q(i_u,i  ,j,k) + q(i_u,i-3,j,k)))* c1r2*(q_tee(i_u,i  ,j,k)+q_tee(i_u,i-3,j,k)) )


!.... fdiv(u) = + .... d/dy( uv) + ... 
       fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
                         c9r8*rjy* ( (c9r16*( q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)) &
                                     -c1r16*( q(i_v,i+1,j  ,k) + q(i_v,i-2,j  ,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i,j+1,k))  &
                                   - (c9r16*( q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)) &
                                     -c1r16*( q(i_v,i+1,j-1,k) + q(i_v,i-2,j-1,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i,j-1,k)) )&
                        -c1r8*rj3y*( (c9r16*( q(i_v,i  ,j+1,k) + q(i_v,i-1,j+1,k)) &
                                     -c1r16*( q(i_v,i+1,j+1,k) + q(i_v,i-2,j+1,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i,j+3,k))  &
                                   - (c9r16*( q(i_v,i  ,j-2,k) + q(i_v,i-1,j-2,k)) &
                                     -c1r16*( q(i_v,i+1,j-2,k) + q(i_v,i-2,j-2,k)))* c1r2*(q_tee(i_u,i,j  ,k)+q_tee(i_u,i,j-3,k)) )


!... fdiv(u) = ... + d/dz(uw )  
       fdiv(i_u,i,j,k) = fdiv(i_u,i,j,k) + &
                         c9r8*rhz* ( (c9r16*( q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)) &
                                     -c1r16*( q(i_w,i+1,j,k+1) + q(i_w,i-2,j,k+1)))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k+1))  &
                                   - (c9r16*( q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )) &
                                     -c1r16*( q(i_w,i+1,j,k  ) + q(i_w,i-2,j,k  )))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k-1)) )&
                        -c1r8*rh3z*( (c9r16*( q(i_w,i  ,j,k+2) + q(i_w,i-1,j,k+2)) &
                                     -c1r16*( q(i_w,i+1,j,k+2) + q(i_w,i-2,j,k+2)))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k+3))  &
                                   - (c9r16*( q(i_w,i  ,j,k-1) + q(i_w,i-1,j,k-1)) &
                                     -c1r16*( q(i_w,i+1,j,k-1) + q(i_w,i-2,j,k-1)))* c1r2*(q_tee(i_u,i,j,k  )+q_tee(i_u,i,j,k-3)) )

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
                                     -c1r16*( q(i_u,i+1,j-1,k) + q(i_u,i+1,j+2,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i+1,j,k))  &
                                   - (c9r16*( q(i_u,i  ,j  ,k) + q(i_u,i  ,j+1,k)) &
                                     -c1r16*( q(i_u,i  ,j-1,k) + q(i_u,i  ,j+2,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i-1,j,k)) )&
                        -c1r8*rh3x*( (c9r16*( q(i_u,i+2,j  ,k) + q(i_u,i+2,j+1,k)) &
                                     -c1r16*( q(i_u,i+2,j-1,k) + q(i_u,i+2,j+2,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i+3,j,k))  &
                                   - (c9r16*( q(i_u,i-1,j  ,k) + q(i_u,i-1,j+1,k)) &
                                     -c1r16*( q(i_u,i-1,j-1,k) + q(i_u,i-1,j+2,k)))* c1r2*(q_tee(i_v,i  ,j,k)+q_tee(i_v,i-3,j,k)) )



!.... fdiv(v) = ... + d/dy(vv) + .... 
       fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
                         c9r8*rjy* ( (c9r16*( q(i_v,i,j  ,k) + q(i_v,i,j+1,k)) &
                                     -c1r16*( q(i_v,i,j-1,k) + q(i_v,i,j+2,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j+1,k))  &
                                   - (c9r16*( q(i_v,i,j  ,k) + q(i_v,i,j-1,k)) &
                                     -c1r16*( q(i_v,i,j+1,k) + q(i_v,i,j-2,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j-1,k)) )&
                        -c1r8*rj3y*( (c9r16*( q(i_v,i,j+1,k) + q(i_v,i,j+2,k)) &
                                     -c1r16*( q(i_v,i,j  ,k) + q(i_v,i,j+3,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j+3,k))  &
                                   - (c9r16*( q(i_v,i,j-1,k) + q(i_v,i,j-2,k)) &
                                     -c1r16*( q(i_v,i,j  ,k) + q(i_v,i,j-3,k)))* c1r2*(q_tee(i_v,i,j  ,k)+q_tee(i_v,i,j-3,k)) )


!... fdiv(v) = ... + d/dz(vw) 
       fdiv(i_v,i,j,k) = fdiv(i_v,i,j,k) + &
                         c9r8*rhz* ( (c9r16*( q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)) &
                                     -c1r16*( q(i_w,i,j-1,k+1) + q(i_w,i,j+2,k+1)))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k+1))  &
                                   - (c9r16*( q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )) &
                                     -c1r16*( q(i_w,i,j-1,k  ) + q(i_w,i,j+2,k  )))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k-1)) )&
                        -c1r8*rh3z*( (c9r16*( q(i_w,i,j  ,k+2) + q(i_w,i,j+1,k+2)) &
                                     -c1r16*( q(i_w,i,j-1,k+2) + q(i_w,i,j+2,k+2)))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k+3))  &
                                   - (c9r16*( q(i_w,i,j  ,k-1) + q(i_w,i,j+1,k-1)) &
                                     -c1r16*( q(i_w,i,j-1,k-1) + q(i_w,i,j+2,k-1)))* c1r2*(q_tee(i_v,i,j,k  )+q_tee(i_v,i,j,k-3)) )

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
                                     -c1r16*( q(i_u,i+1,j,k+1) + q(i_u,i+1,j,k-2)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i+1,j,k))  &
                                   - (c9r16*( q(i_u,i  ,j,k  ) + q(i_u,i  ,j,k-1)) &
                                     -c1r16*( q(i_u,i  ,j,k+1) + q(i_u,i  ,j,k-2)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i-1,j,k)) )&
                        -c1r8*rh3x*( (c9r16*( q(i_u,i+2,j,k  ) + q(i_u,i+2,j,k-1)) &
                                     -c1r16*( q(i_u,i+2,j,k+1) + q(i_u,i+2,j,k-2)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i+3,j,k))  &
                                   - (c9r16*( q(i_u,i-1,j,k  ) + q(i_u,i-1,j,k-1)) &
                                     -c1r16*( q(i_u,i-1,j,k+1) + q(i_u,i-1,j,k-2)))* c1r2*(q_tee(i_w,i  ,j,k)+q_tee(i_w,i-3,j,k)) )


!.... fdiv(w) = ... + d/dy(wv) + .... 
       fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
                         c9r8*rjy* ( (c9r16*( q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) &
                                     -c1r16*( q(i_v,i,j  ,k+1) + q(i_v,i,j  ,k-2)))* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j+1,k))  &
                                   - (c9r16*( q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) &
                                     -c1r16*( q(i_v,i,j-1,k+1) + q(i_v,i,j-1,k-2)))* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j-1,k)) )&
                        -c1r8*rj3y*( (c9r16*( q(i_v,i,j+1,k  ) + q(i_v,i,j+1,k-1)) &
                                     -c1r16*( q(i_v,i,j+1,k+1) + q(i_v,i,j+1,k-2)))* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j+3,k))  &
                                   - (c9r16*( q(i_v,i,j-2,k  ) + q(i_v,i,j-2,k-1)) &
                                     -c1r16*( q(i_v,i,j-2,k+1) + q(i_v,i,j-2,k-2)))* c1r2*(q_tee(i_w,i,j  ,k)+q_tee(i_w,i,j-3,k)) )


!.... fdiv(w) = ... + d/dz(ww)       
       fdiv(i_w,i,j,k) = fdiv(i_w,i,j,k) + &
                         c9r8*rhz* ( (c9r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k+1)) &
                                     -c1r16*( q(i_w,i,j,k-1) + q(i_w,i,j,k+2)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k+1))  &
                                   - (c9r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k-1)) &
                                     -c1r16*( q(i_w,i,j,k+1) + q(i_w,i,j,k-2)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k-1))) &
                        -c1r8*rh3z*( (c9r16*( q(i_w,i,j,k+1) + q(i_w,i,j,k+2)) &
                                     -c1r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k+3)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k+3))  &
                                   - (c9r16*( q(i_w,i,j,k-1) + q(i_w,i,j,k-2)) &
                                     -c1r16*( q(i_w,i,j,k  ) + q(i_w,i,j,k-3)))* c1r2*(q_tee(i_w,i,j,k  )+q_tee(i_w,i,j,k-3)))


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
                                         rhx*   (q_tee(i_u,i  ,j,k) - q_tee(i_u,i-1,j,k)) &
                                       +(c9r16* (q(i_u,i  ,j,k) + q(i_u,i+1,j,k)) - c1r16*(q(i_u,i-1,j,k) + q(i_u,i+2,j,k)))* &
                                         rhx*   (q_tee(i_u,i+1,j,k) - q_tee(i_u,i  ,j,k)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_u,i-2,j,k) + q(i_u,i-1,j,k)) - c1r16*(q(i_u,i-3,j,k) + q(i_u,i  ,j,k)))* &
                                         rh3x*  (q_tee(i_u,i  ,j,k) - q_tee(i_u,i-3,j,k)) &
                                       +(c9r16* (q(i_u,i+2,j,k) + q(i_u,i+1,j,k)) - c1r16*(q(i_u,i+3,j,k) + q(i_u,i  ,j,k)))* &
                                         rh3x*  (q_tee(i_u,i+3,j,k) - q_tee(i_u,i  ,j,k)) )


!.... f(i_u) =  + ... v du/dy + ... 
       fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + rjy*( &
                         c9r8 * c1r2*(  (c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i-1,j  ,k)) - c1r16*(q(i_v,i-2,j  ,k) + q(i_v,i+1,j  ,k)))* &
                                         rhy*   (q_tee(i_u,i  ,j+1,k) - q_tee(i_u,i  ,j  ,k)) &
                                       +(c9r16* (q(i_v,i  ,j-1,k) + q(i_v,i-1,j-1,k)) - c1r16*(q(i_v,i-2,j-1,k) + q(i_v,i+1,j-1,k)))* &
                                         rhy*   (q_tee(i_u,i  ,j  ,k) - q_tee(i_u,i  ,j-1,k)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_v,i  ,j-2,k) + q(i_v,i-1,j-2,k)) - c1r16*(q(i_v,i-2,j-2,k) + q(i_v,i+1,j-2,k)))* &
                                         rh3y*  (q_tee(i_u,i  ,j  ,k) - q_tee(i_u,i  ,j-3,k)) &
                                       +(c9r16* (q(i_v,i  ,j+1,k) + q(i_v,i-1,j+1,k)) - c1r16*(q(i_v,i-2,j+1,k) + q(i_v,i+1,j+1,k)))* &
                                         rh3y*  (q_tee(i_u,i  ,j+3,k) - q_tee(i_u,i  ,j  ,k)) ) )

!.... f(i_u) =  ... + w du/dz 
       fadv(i_u,i,j,k) = fadv(i_u,i,j,k) + &
                         c9r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k+1) + q(i_w,i-1,j,k+1)) - c1r16*(q(i_w,i-2,j,k+1) + q(i_w,i+1,j,k+1)))* &
                                         rhz*   (q_tee(i_u,i  ,j,k+1) - q_tee(i_u,i  ,j,k  )) &
                                       +(c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i-1,j,k  )) - c1r16*(q(i_w,i-2,j,k  ) + q(i_w,i+1,j,k  )))* &
                                         rhz*   (q_tee(i_u,i  ,j,k  ) - q_tee(i_u,i  ,j,k-1)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k-1) + q(i_w,i-1,j,k-1)) - c1r16*(q(i_w,i-2,j,k-1) + q(i_w,i+1,j,k-1)))* &
                                         rh3z*  (q_tee(i_u,i  ,j,k  ) - q_tee(i_u,i  ,j,k-3)) &
                                       +(c9r16* (q(i_w,i  ,j,k+2) + q(i_w,i-1,j,k+2)) - c1r16*(q(i_w,i-2,j,k+2) + q(i_w,i+1,j,k+2)))* &
                                         rh3z*  (q_tee(i_u,i  ,j,k+3) - q_tee(i_u,i  ,j,k  )) )

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
                                         rhx*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i-1,j  ,k)) &
                                       +(c9r16* (q(i_u,i+1,j  ,k) + q(i_u,i+1,j+1,k)) - c1r16*(q(i_u,i+1,j-1,k) + q(i_u,i+1,j+2,k)))* &
                                         rhx*   (q_tee(i_v,i+1,j  ,k) - q_tee(i_v,i  ,j  ,k)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_u,i-1,j  ,k) + q(i_u,i-1,j+1,k)) - c1r16*(q(i_u,i-1,j-1,k) + q(i_u,i-1,j+2,k)))* &
                                         rh3x*  (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i-3,j  ,k)) &
                                       +(c9r16* (q(i_u,i+2,j  ,k) + q(i_u,i+2,j+1,k)) - c1r16*(q(i_u,i+2,j-1,k) + q(i_u,i+2,j+2,k)))* &
                                         rh3x*  (q_tee(i_v,i+3,j  ,k) - q_tee(i_v,i  ,j  ,k)) )

!.... f(i_v) = ... + v dv/dy + .... 
       fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + rjy*( &
                         c9r8 * c1r2* ( (c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j+1,k)) - c1r16*(q(i_v,i  ,j-1,k) + q(i_v,i  ,j+2,k)))* &
                                         rhy*   (q_tee(i_v,i  ,j+1,k) - q_tee(i_v,i  ,j  ,k)) &
                                       +(c9r16* (q(i_v,i  ,j  ,k) + q(i_v,i  ,j-1,k)) - c1r16*(q(i_v,i  ,j+1,k) + q(i_v,i  ,j-2,k)))* &
                                         rhy*   (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i  ,j-1,k)) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_v,i  ,j+1,k) + q(i_v,i  ,j+2,k)) - c1r16*(q(i_v,i  ,j  ,k) + q(i_v,i  ,j+3,k)))* &
                                         rh3y*  (q_tee(i_v,i  ,j+3,k) - q_tee(i_v,i  ,j  ,k)) &
                                       +(c9r16* (q(i_v,i  ,j-1,k) + q(i_v,i  ,j-2,k)) - c1r16*(q(i_v,i  ,j  ,k) + q(i_v,i  ,j-3,k)))* &
                                         rh3y*  (q_tee(i_v,i  ,j  ,k) - q_tee(i_v,i  ,j-3,k)) ) )

!.... f(i_v) = ... + w dv/dz 
       fadv(i_v,i,j,k) = fadv(i_v,i,j,k) + &
                         c9r8 * c1r2* ( (c9r16* (q(i_w,i,j  ,k  ) + q(i_w,i,j+1,k  )) - c1r16*(q(i_w,i,j-1,k  ) + q(i_w,i,j+2,k  )))* &
                                         rhz*   (q_tee(i_v,i,j  ,k  ) - q_tee(i_v,i,j  ,k-1)) &
                                       +(c9r16* (q(i_w,i,j  ,k+1) + q(i_w,i,j+1,k+1)) - c1r16*(q(i_w,i,j-1,k+1) + q(i_w,i,j+2,k+1)))* &
                                         rhz*   (q_tee(i_v,i,j  ,k+1) - q_tee(i_v,i,j  ,k  )) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_w,i,j  ,k-1) + q(i_w,i,j+1,k-1)) - c1r16*(q(i_w,i,j-1,k-1) + q(i_w,i,j+2,k-1)))* &
                                         rh3z*  (q_tee(i_v,i,j  ,k  ) - q_tee(i_v,i,j  ,k-3)) &
                                       +(c9r16* (q(i_w,i,j  ,k+2) + q(i_w,i,j+1,k+2)) - c1r16*(q(i_w,i,j-1,k+2) + q(i_w,i,j+2,k+2)))* &
                                         rh3z*  (q_tee(i_v,i,j  ,k+3) - q_tee(i_v,i,j  ,k  )) )


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
                                         rhx*   (q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i-1,j,k  )) &
                                       +(c9r16* (q(i_u,i+1,j,k  ) + q(i_u,i+1,j,k-1)) - c1r16*(q(i_u,i+1,j,k-2) + q(i_u,i+1,j,k+1)))* &
                                         rhx*   (q_tee(i_w,i+1,j,k  ) - q_tee(i_w,i  ,j,k  )) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_u,i-1,j,k  ) + q(i_u,i-1,j,k-1)) - c1r16*(q(i_u,i-1,j,k-2) + q(i_u,i-1,j,k+1)))* &
                                         rh3x*  (q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i-3,j,k  )) &
                                       +(c9r16* (q(i_u,i+2,j,k  ) + q(i_u,i+2,j,k-1)) - c1r16*(q(i_u,i+2,j,k-2) + q(i_u,i+2,j,k+1)))* &
                                         rh3x*  (q_tee(i_w,i+3,j,k  ) - q_tee(i_w,i  ,j,k  )) )


!.... fadv(i_w) = ...+ v dw/dy + ... 
       fadv(i_w,i,j,k) = fadv(i_w,i,j,k) +  rjy*( &
                         c9r8 * c1r2* ( (c9r16* (q(i_v,i,j  ,k  ) + q(i_v,i,j  ,k-1)) - c1r16*(q(i_v,i,j  ,k+1) + q(i_v,i,j  ,k-2)))* &
                                         rhy*   (q_tee(i_w,i,j+1,k  ) - q_tee(i_w,i,j  ,k  )) &
                                       +(c9r16* (q(i_v,i,j-1,k  ) + q(i_v,i,j-1,k-1)) - c1r16*(q(i_v,i,j-1,k+1) + q(i_v,i,j-1,k-2)))* &
                                         rhy*   (q_tee(i_w,i,j  ,k  ) - q_tee(i_w,i,j-1,k  )) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_v,i,j-2,k  ) + q(i_v,i,j-2,k-1)) - c1r16*(q(i_v,i,j-2,k+1) + q(i_v,i,j-2,k-2)))* &
                                         rh3y*  (q_tee(i_w,i,j  ,k  ) - q_tee(i_w,i,j-3,k  )) &
                                       +(c9r16* (q(i_v,i,j+1,k  ) + q(i_v,i,j+1,k-1)) - c1r16*(q(i_v,i,j+1,k+1) + q(i_v,i,j+1,k-2)))* &
                                         rh3y*  (q_tee(i_w,i,j+3,k  ) - q_tee(i_w,i,j  ,k  )) ) )

!.... fadv(i_w) = ... + w dw/dz 
       fadv(i_w,i,j,k) = fadv(i_w,i,j,k) + &
                         c9r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-1)) - c1r16*(q(i_w,i  ,j,k+1) + q(i_w,i  ,j,k-2)))* &
                                         rhz*   (q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i  ,j,k-1)) &
                                       +(c9r16* (q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+1)) - c1r16*(q(i_w,i  ,j,k-1) + q(i_w,i  ,j,k+2)))* &
                                         rhz*   (q_tee(i_w,i  ,j,k+1) - q_tee(i_w,i  ,j,k  )) ) &
                        -c1r8 * c1r2* ( (c9r16* (q(i_w,i  ,j,k-1) + q(i_w,i  ,j,k-2)) - c1r16*(q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k-3)))* &
                                         rh3z*  (q_tee(i_w,i  ,j,k  ) - q_tee(i_w,i  ,j,k-3)) &
                                       +(c9r16* (q(i_w,i  ,j,k+1) + q(i_w,i  ,j,k+2)) - c1r16*(q(i_w,i  ,j,k  ) + q(i_w,i  ,j,k+3)))* &
                                         rh3z*  (q_tee(i_w,i  ,j,k+3) - q_tee(i_w,i  ,j,k  )) )

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


  end subroutine cnv4th_tee

!======================================================================================================================================!

  
  subroutine difwn_tee ( lun, ifs, ife, kfs, kfe, wvisc, vtb ) 

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

       !if ( wvisc .eq. eddy_visc) call eddyv_op ( i, k)

       !---
       ! u,w 
       !---
       do j=1,n2

          if ( wvisc .eq. eddy_visc ) then

             rre_u = c9r16* ( vt(i-1,j,k) + vt(i  ,j,k)) - c1r16* ( vt(i-2,j,k) + vt(i+1,j,k))
             rre_w = c9r16* ( vt(i,j,k-1) + vt(i,j,k  )) - c1r16* ( vt(i,j,k-2) + vt(i,j,k+1))
          endif

          do ii=max(1,j-ku), min(n2,j+kl) 
          ! (Kim 09.14.22)
            ! u
            if ( pat_bc .and. match_bc .and. phase(i_u,i,k) .ne. 0 ) then
              lun(i_u,i,j,k) = lun(i_u,i,j,k) + lub_pat(ku+kl+1+j-ii,ii)* rre_u* q_tee(i_u,i,ii,k) + tsgs* fsgs*evb_u(ku+kl+1+j-ii,ii)* q_tee(i_u,i,ii,k)
            else
              lun(i_u,i,j,k) = lun(i_u,i,j,k) + lub_tee(ku+kl+1+j-ii,ii)* rre_u* q_tee(i_u,i,ii,k) + tsgs* fsgs*evb_u(ku+kl+1+j-ii,ii)* q_tee(i_u,i,ii,k)
            endif

            ! w
            if ( pat_bc .and. match_bc .and. phase(i_w,i,k) .ne. 0 ) then
              lun(i_w,i,j,k) = lun(i_w,i,j,k) + lub_pat(ku+kl+1+j-ii,ii)* rre_w* q_tee(i_w,i,ii,k) + tsgs* fsgs*evb_w(ku+kl+1+j-ii,ii)* q_tee(i_w,i,ii,k)
            else
              lun(i_w,i,j,k) = lun(i_w,i,j,k) + lub_tee(ku+kl+1+j-ii,ii)* rre_w* q_tee(i_w,i,ii,k) + tsgs* fsgs*evb_w(ku+kl+1+j-ii,ii)* q_tee(i_w,i,ii,k)
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
             lun(i_v,i,j,k) = lun(i_v,i,j,k) + lvb(ku+kl+1+j-ii,ii)* rre_v* q_tee(i_v,i,ii,k) + tsgs* fsgs*evb_v(ku+kl+1+j-ii,ii)* q_tee(i_v,i,ii,k)
          enddo
       enddo

    enddo
    enddo

  end subroutine difwn_tee

!=================================================================================================!

  subroutine sw_calc

    implicit none
    integer                                 :: i,k,jBF
    real, dimension(:), allocatable         :: sw_temp, sw_sum
    integer                                 :: nphase_u
    real                                    :: deltaBF

    allocate( sw_temp (1:2) )
    allocate( sw_sum (1:2) )
        
    nphase_u = n1*n3-nsolid_u
    sw_temp = 0.d0
    sw_sum = 0.d0

    ! sum over air
    do k=ks,ke
    do i=is,ie
      if  (phase(i_u,i,k) .ne. 0) then
        sw_temp(1) = sw_temp(1) + q_tee(i_u,i, 1,k)
        sw_temp(2) = sw_temp(2) + q_tee(i_u,i,n2,k)
      endif
    enddo
    enddo

    call mpi_allreduce (sw_temp, sw_sum, 2, MPI_DOUBLE_PRECISION, mpi_sum, mycomm, ierr)

    sw_gb = sw_sum(1)
    sw_gt = sw_sum(2)

    ! BF edit
    jBF = ICswitchBF
    deltaBF = y(i_u,jBF+1) - y(i_u,jBF)
    !sb_new = ( -1.d0*n1*n3 - nphase_u*y(i_u, 1) - sw_sum(1) ) / nsolid_u
    !st_new = (  1.d0*n1*n3 - nphase_u*y(i_u,n2) - sw_sum(2) ) / nsolid_u
    sb_new = (                 - sw_sum(1) ) / nsolid_u
    st_new = ( deltaBF*(n1*n3) - sw_sum(2) ) / nsolid_u



  end subroutine 

!===================================================================================================================! 


  subroutine bcuvw_tee
    
!.... u,v,w boundary conditions 
    implicit none 
    integer                             :: i,j,k,i_q , jss,jBF
    integer, dimension(mpi_status_size) :: status 
    real :: dx, dy, hy0, jac, djy, dz, deltaBF

    hy0  = 2.d0/ dble(n2)
    jac  = jacb(i_u, 1 )
    djy  = hy0*jac
    dx = chlx * pi / dble(n1)
    dz = chlz * pi / dble(n3)
    

!.... apply boundary conditions in wn direction, 
    do k=ks,ke
    do i=is,ie

    if ( BCswitch == 0 ) then
    q_tee(i_u,i,   0,k) = - q_tee(i_u,i,1,k)     ! u-1 = u0     makes the boundary u velocity zero if we use second order scheme with staggered mesh
    else
    q_tee(i_u,i,   0,k) = 2.d0 * v1_bottom - q_tee(i_u,i,1,k)
       !if ( rank .eq. 0) write(*,*) v1_bottom 
    endif
    
    q_tee(i_u,i,  -1,k) = 0.d0               ! extra ghost cell. we do not use this. 
    q_tee(i_u,i,  -2,k) = 0.d0               ! extra ghost cell. we do not use this. 

    if ( BCswitch == 0 ) then
    q_tee(i_u,i,n2+1,k) = - q_tee(i_u,i,n2,k)
    else
    q_tee(i_u,i,n2+1,k) = 2.d0 * v1_top - q_tee(i_u,i,n2,k)
    endif

    q_tee(i_u,i,n2+2,k) = 0.d0
    q_tee(i_u,i,n2+3,k) = 0.d0

    q_tee(i_w,i,   0,k) = - q_tee(i_w,i,1,k)
    q_tee(i_w,i,  -1,k) = 0.d0
    q_tee(i_w,i,  -2,k) = 0.d0
    q_tee(i_w,i,n2+1,k) = - q_tee(i_w,i,n2,k)
    q_tee(i_w,i,n2+2,k) = 0.d0
    q_tee(i_w,i,n2+3,k) = 0.d0

!.... slip BC (Kim: 08.31.22) (Kim: 09.13.22)
    if ( slip_bc .and. match_bc ) then
      ! x direction
      q_tee(i_u,i,   0,k) = cf_slip(i_u,0,1) * q_tee(i_u,i, 1,k) + &
          cf_slip(i_u,0,2) * q_tee(i_u,i,   2,k)
      q_tee(i_u,i,n2+1,k) = cf_slip(i_u,0,1) * q_tee(i_u,i,n2,k) + &
          cf_slip(i_u,0,2) * q_tee(i_u,i,n2-1,k)
      ! z direction
      q_tee(i_w,i,   0,k) = cf_slip(i_w,0,1) * q_tee(i_w,i, 1,k) + &
          cf_slip(i_w,0,2) * q_tee(i_w,i,   2,k)
      q_tee(i_w,i,n2+1,k) = cf_slip(i_w,0,1) * q_tee(i_w,i,n2,k) + &
            cf_slip(i_w,0,2) * q_tee(i_w,i,n2-1,k)
    endif
!.... end slip BC

!.... pat BC (Kim 09.14.22)
    if ( pat_bc .and. match_bc ) then

      ! BF edit
      jBF = ICswitchBF
      deltaBF = y(i_u,jBF+1) - y(i_u,jBF)
      if ( phase(i_u,i,k) .ne. 0 ) then ! air
        q_tee(i_u,i,   0,k) = q_tee(i_u,i, 1,k)
        q_tee(i_u,i,n2+1,k) = q_tee(i_u,i,n2,k)
      else ! solid
        q_tee(i_u,i,   0,k) = -q_tee(i_u,i, 1,k) + 2.d0*(sb_new)
        q_tee(i_u,i,n2+1,k) = -q_tee(i_u,i,n2,k) + 2.d0*(st_new)
      endif

      if ( phase(i_w,i,k) .ne. 0 ) then
        q_tee(i_w,i,0,k) = q_tee(i_w,i,1,k)
        q_tee(i_w,i,n2+1,k) = q_tee(i_w,i,n2,k)
      endif

    endif
!.... end pat BC


!.... v bc 
    !q(i_v,i,   0,k) = 0.d0
    if (firstloop ) q_tee(i_v,i,   0,k) = 0.d0      
    q_tee(i_v,i,  -1,k) = 0.d0
    q_tee(i_v,i,  -2,k) = 0.d0
    
    !q_tee(i_v,i,n2  ,k) = 0.d0
    if (firstloop ) q_tee(i_v,i,n2  ,k) = 0.d0
    q_tee(i_v,i,n2+1,k) = 0.d0
    q_tee(i_v,i,n2+2,k) = 0.d0

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
          call mpi_sendrecv ( q_tee(i_u,ie-halo_x+1,jss,ks-halo_z), 1, uvec_comm, right, 0, & 
                              q_tee(i_u,ie+1       ,jss,ks-halo_z), 1, uvec_comm, right, 0, comm3d, status, ierr)

          call mpi_sendrecv ( q_tee(i_u,is         ,jss,ks-halo_z), 1, uvec_comm, left , 0, & 
                              q_tee(i_u,is-halo_x  ,jss,ks-halo_z), 1, uvec_comm, left , 0, comm3d, status, ierr) 

       else  

          call mpi_sendrecv ( q_tee(i_u,is         ,jss,ks-halo_z), 1, uvec_comm, left , 0, & 
                              q_tee(i_u,is-halo_x  ,jss,ks-halo_z), 1, uvec_comm, left , 0, comm3d, status, ierr) 

          call mpi_sendrecv ( q_tee(i_u,ie-halo_x+1,jss,ks-halo_z), 1, uvec_comm, right, 0, & 
                              q_tee(i_u,ie+1       ,jss,ks-halo_z), 1, uvec_comm, right, 0, comm3d, status, ierr) 
       endif

    else  

       do k=ks,ke
       do j=1-halo_yu,n2+halo_yu
       do i_q=1,nvel 

          q_tee(i_q,-3,j,k) = q_tee(i_q,n1-3,j,k) 
          q_tee(i_q,-2,j,k) = q_tee(i_q,n1-2,j,k) 
          q_tee(i_q,-1,j,k) = q_tee(i_q,n1-1,j,k) 


          q_tee(i_q,n1,j,k  ) = q_tee(i_q,0,j,k) 
          q_tee(i_q,n1+1,j,k) = q_tee(i_q,1,j,k) 
          q_tee(i_q,n1+2,j,k) = q_tee(i_q,2,j,k) 

       enddo
       enddo
       enddo


    endif
       

!.... perform front<->back communication ; all tags are set to 0 
    
    if ( p3 > 1 ) then 
       if ( mod( coords(2), 2) .eq. 0 ) then 
          call mpi_sendrecv ( q_tee(i_u,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                              q_tee(i_u,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr)
          
          call mpi_sendrecv ( q_tee(i_u,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                              q_tee(i_u,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
       else  
          call mpi_sendrecv ( q_tee(i_u,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                              q_tee(i_u,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
          call mpi_sendrecv ( q_tee(i_u,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                              q_tee(i_u,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
       endif

    else  

       do j=1-halo_yu,n2+halo_yu
       do i=is-halo_x,ie+halo_x
       do i_q=1,nvel

          q_tee(i_q,i,j,-3)   = q_tee(i_q,i,j,n3-3) 
          q_tee(i_q,i,j,-2)   = q_tee(i_q,i,j,n3-2) 
          q_tee(i_q,i,j,-1)   = q_tee(i_q,i,j,n3-1) 
       
          q_tee(i_q,i,j,n3 )  = q_tee(i_q,i,j,0   ) 
          q_tee(i_q,i,j,n3+1) = q_tee(i_q,i,j,1   ) 
          q_tee(i_q,i,j,n3+2) = q_tee(i_q,i,j,2   ) 

    enddo
    enddo
    enddo

    endif


  end subroutine bcuvw_tee

!=================================================================================================!

  subroutine bcuvw_source
    
!.... u,v,w boundary conditions 
    implicit none 
    integer                             :: i,j,k,i_q , jss,jBF
    integer, dimension(mpi_status_size) :: status 
    real :: dx, dy, hy0, jac, djy, dz, deltaBF

    hy0  = 2.d0/ dble(n2)
    jac  = jacb(i_u, 1 )
    djy  = hy0*jac
    dx = chlx * pi / dble(n1)
    dz = chlz * pi / dble(n3)
    

!.... apply boundary conditions in wn direction, 
    do k=ks,ke
    do i=is,ie

    if ( BCswitch == 0 ) then
    q_source(i_u,i,   0,k) = - q_source(i_u,i,1,k)     ! u-1 = u0     makes the boundary u velocity zero if we use second order scheme with staggered mesh
    else
    q_source(i_u,i,   0,k) = 2.d0 * v1_bottom - q_source(i_u,i,1,k)
    endif
    
    q_source(i_u,i,  -1,k) = 0.d0               ! extra ghost cell. we do not use this. 
    q_source(i_u,i,  -2,k) = 0.d0               ! extra ghost cell. we do not use this. 

    if ( BCswitch == 0 ) then
    q_source(i_u,i,n2+1,k) = - q_source(i_u,i,n2,k)
    else
    q_source(i_u,i,n2+1,k) = 2.d0 * v1_top - q_source(i_u,i,n2,k)
    endif

    q_source(i_u,i,n2+2,k) = 0.d0
    q_source(i_u,i,n2+3,k) = 0.d0

    q_source(i_w,i,   0,k) = - q_source(i_w,i,1,k)
    q_source(i_w,i,  -1,k) = 0.d0
    q_source(i_w,i,  -2,k) = 0.d0
    q_source(i_w,i,n2+1,k) = - q_source(i_w,i,n2,k)
    q_source(i_w,i,n2+2,k) = 0.d0
    q_source(i_w,i,n2+3,k) = 0.d0

!.... slip BC (Kim: 08.31.22) (Kim: 09.13.22)
    if ( slip_bc .and. match_bc ) then
      ! x direction
      q_source(i_u,i,   0,k) = cf_slip(i_u,0,1) * q_source(i_u,i, 1,k) + &
          cf_slip(i_u,0,2) * q_source(i_u,i,   2,k)
      q_source(i_u,i,n2+1,k) = cf_slip(i_u,0,1) * q_source(i_u,i,n2,k) + &
          cf_slip(i_u,0,2) * q_source(i_u,i,n2-1,k)
      ! z direction
      q_source(i_w,i,   0,k) = cf_slip(i_w,0,1) * q_source(i_w,i, 1,k) + &
          cf_slip(i_w,0,2) * q_source(i_w,i,   2,k)
      q_source(i_w,i,n2+1,k) = cf_slip(i_w,0,1) * q_source(i_w,i,n2,k) + &
            cf_slip(i_w,0,2) * q_source(i_w,i,n2-1,k)
    endif
!.... end slip BC

!.... pat BC (Kim: 09.14.22)
    if ( pat_bc .and. match_bc ) then

      ! BF edit
      jBF = ICswitchBF
      deltaBF = y(i_u,jBF+1) - y(i_u,jBF)

      if ( phase(i_u,i,k) .ne. 0 ) then
        q_source(i_u,i,   0,k) = q_source(i_u,i, 1,k)
        q_source(i_u,i,n2+1,k) = q_source(i_u,i,n2,k)
      else
        q_source(i_u,i,   0,k) = -q_source(i_u,i, 1,k) + 2.d0*(sb_new)
        q_source(i_u,i,n2+1,k) = -q_source(i_u,i,n2,k) + 2.d0*(st_new)
      endif

      if ( phase(i_w,i,k) .ne. 0 ) then
        q_source(i_w,i,0,k) = q_source(i_w,i,1,k)
        q_source(i_w,i,n2+1,k) = q_source(i_w,i,n2,k)
      endif

    endif
!.... end pat BC


!.... v bc 
    !q(i_v,i,   0,k) = 0.d0
    if (firstloop ) q_source(i_v,i,   0,k) = 0.d0      
    q_source(i_v,i,  -1,k) = 0.d0
    q_source(i_v,i,  -2,k) = 0.d0
    
    !q_source(i_v,i,n2  ,k) = 0.d0
    if (firstloop ) q_source(i_v,i,n2  ,k) = 0.d0
    q_source(i_v,i,n2+1,k) = 0.d0
    q_source(i_v,i,n2+2,k) = 0.d0

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
          call mpi_sendrecv ( q_source(i_u,ie-halo_x+1,jss,ks-halo_z), 1, uvec_comm, right, 0, & 
                              q_source(i_u,ie+1       ,jss,ks-halo_z), 1, uvec_comm, right, 0, comm3d, status, ierr)

          call mpi_sendrecv ( q_source(i_u,is         ,jss,ks-halo_z), 1, uvec_comm, left , 0, & 
                              q_source(i_u,is-halo_x  ,jss,ks-halo_z), 1, uvec_comm, left , 0, comm3d, status, ierr) 

       else  

          call mpi_sendrecv ( q_source(i_u,is         ,jss,ks-halo_z), 1, uvec_comm, left , 0, & 
                              q_source(i_u,is-halo_x  ,jss,ks-halo_z), 1, uvec_comm, left , 0, comm3d, status, ierr) 

          call mpi_sendrecv ( q_source(i_u,ie-halo_x+1,jss,ks-halo_z), 1, uvec_comm, right, 0, & 
                              q_source(i_u,ie+1       ,jss,ks-halo_z), 1, uvec_comm, right, 0, comm3d, status, ierr) 
       endif

    else  

       do k=ks,ke
       do j=1-halo_yu,n2+halo_yu
       do i_q=1,nvel 

          q_source(i_q,-3,j,k) = q_source(i_q,n1-3,j,k) 
          q_source(i_q,-2,j,k) = q_source(i_q,n1-2,j,k) 
          q_source(i_q,-1,j,k) = q_source(i_q,n1-1,j,k) 


          q_source(i_q,n1,j,k  ) = q_source(i_q,0,j,k) 
          q_source(i_q,n1+1,j,k) = q_source(i_q,1,j,k) 
          q_source(i_q,n1+2,j,k) = q_source(i_q,2,j,k) 

       enddo
       enddo
       enddo


    endif
       

!.... perform front<->back communication ; all tags are set to 0 
    
    if ( p3 > 1 ) then 
       if ( mod( coords(2), 2) .eq. 0 ) then 
          call mpi_sendrecv ( q_source(i_u,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                              q_source(i_u,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr)
          
          call mpi_sendrecv ( q_source(i_u,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                              q_source(i_u,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
       else  
          call mpi_sendrecv ( q_source(i_u,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                              q_source(i_u,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
          call mpi_sendrecv ( q_source(i_u,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                              q_source(i_u,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
       endif

    else  

       do j=1-halo_yu,n2+halo_yu
       do i=is-halo_x,ie+halo_x
       do i_q=1,nvel

          q_source(i_q,i,j,-3)   = q_source(i_q,i,j,n3-3) 
          q_source(i_q,i,j,-2)   = q_source(i_q,i,j,n3-2) 
          q_source(i_q,i,j,-1)   = q_source(i_q,i,j,n3-1) 
       
          q_source(i_q,i,j,n3 )  = q_source(i_q,i,j,0   ) 
          q_source(i_q,i,j,n3+1) = q_source(i_q,i,j,1   ) 
          q_source(i_q,i,j,n3+2) = q_source(i_q,i,j,2   ) 

    enddo
    enddo
    enddo

    endif


  end subroutine bcuvw_source

!================================================================================================! 
  subroutine compute_lub_tee

  !------------------------------
  ! wall normal diffusion operator along staggered mesh nodes 
  ! for u,w , stored in lapack banded matrix storage (pentadiagonal) 
  !------------------------------ 

  implicit none 
  integer                                 :: i,j
  real, dimension(5)                      :: coeff 
  real, dimension(3)                      :: coeff_2nd
  real                                    :: hm1, hm2, hp1, hp2

  lub_tee  = 0.d0 
  coeff = 0.d0
  coeff_2nd = 0.d0

  do j=1,n2 
  do i = max(1, j-ku), min(n2,j+kl) 

  if ( i .ge. 3 .and. i .le. n2-2 ) then  

  hp1 = y(i_u,i+1) - y(i_u,i  ) 
  hm1 = y(i_u,i  ) - y(i_u,i-1) 

  call d2_coeff_2nd ( hm1, hp1, coeff_2nd )

  if ( i-2 .eq. j ) then 
  lub_tee(ku+kl+1+i-j,j) = 0.d0
  elseif( i-1 .eq. j) then 
  lub_tee(ku+kl+1+i-j,j) = coeff_2nd(1) 
  elseif( i+1 .eq. j) then 
  lub_tee(ku+kl+1+i-j,j) = coeff_2nd(3) 
  elseif( i+2 .eq. j) then 
  lub_tee(ku+kl+1+i-j,j) = 0.d0
  elseif ( i .eq. j ) then
  lub_tee(ku+kl+1+i-j,j) = coeff_2nd(2) 
  else 
  lub_tee(ku+kl+1+i-j,j) =  0.d0 
  endif

  endif
  enddo
  enddo

  ! ----

    i   = 1
    hp1 = y(i_u,i+1) - y(i_u,i )
    hm1 = 2.d0 * ( y(i_u,1) - y(i_u,0) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lub_tee(ku+kl+1,i  ) = cf(i_u,0,1)*coeff_2nd(1) + coeff_2nd(2)
    lub_tee(ku+kl  ,i+1) = cf(i_u,0,2)*coeff_2nd(1) + coeff_2nd(3)
    lub_tee(ku+kl-1,i+2) = 0.d0
 

    i   = 2
    hp1 = y(i_u,i+1) - y(i_u,i  )
    hm1 = y(i_u,i  ) - y(i_u,i-1)

    call d2_coeff_2nd( hm1, hp1, coeff_2nd)

    lub_tee(ku+kl+2,i-1) =  coeff_2nd(1)
    lub_tee(ku+kl+1,i  ) =  coeff_2nd(2)
    lub_tee(ku+kl  ,i+1) =  coeff_2nd(3)
    lub_tee(ku+kl-1,i+2) =  0.d0



    i   = n2-1
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = y(i_u,i+1) - y(i_u,i  )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lub_tee(ku+kl+3,i-2) = 0.d0
    lub_tee(ku+kl+2,i-1) = coeff_2nd(1)
    lub_tee(ku+kl+1,i  ) = coeff_2nd(2)
    lub_tee(ku+kl  ,i+1) = coeff_2nd(3)


    i   = n2
    hm1 = y(i_u,i  ) - y(i_u,i-1)
    hp1 = 2.d0 * ( y(i_u,n2+1) - y(i_u,n2) )

    call d2_coeff_2nd( hm1, hp1, coeff_2nd )

    lub_tee(ku+kl+3,i-2) = 0.d0
    lub_tee(ku+kl+2,i-1) = cf(i_u,0,2)*coeff_2nd(3) + coeff_2nd(1)
    lub_tee(ku+kl+1,i  ) = cf(i_u,0,1)*coeff_2nd(3) + coeff_2nd(2)
 
   
    end subroutine compute_lub_tee
!======================================================================================================================================!

  subroutine wn_op_tee ( wn_u_tee, ipiv_u_tee, wn_v_tee, ipiv_v_tee, wn_w_tee, ipiv_w_tee, beta , ii, kk, wvisc, vtb) 

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
    real, dimension(1:2*kl+ku+1,1:n2)       :: wn_u_tee 
    real, dimension(1:2*kl+ku+1,1:n2)       :: wn_w_tee 
    real, dimension(1:2*kl+ku+1,1:n2-1)     :: wn_v_tee 
    integer, dimension(1:n2 )               :: ipiv_u_tee 
    integer, dimension(1:n2 )               :: ipiv_w_tee 
    integer, dimension(1:n2-1)              :: ipiv_v_tee 
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
          wn_u_tee(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_u* lub_tee(ku+kl+1+i-j,j) 
          wn_w_tee(ku+kl+1+i-j, j) = 1.d0 - dt*beta* rre_w* lub_tee(ku+kl+1+i-j,j) 
       else 
          wn_u_tee(ku+kl+1+i-j, j) = -dt*beta* rre_u* lub_tee(ku+kl+1+i-j,j) 
          wn_w_tee(ku+kl+1+i-j, j) = -dt*beta* rre_w* lub_tee(ku+kl+1+i-j,j)
       endif
    enddo
    enddo

    !---- 
    ! operator from sgs gradient term 
    !----
    call dgbtrf ( n2, n2, kl, ku, wn_u_tee, 2*kl+ku+1, ipiv_u_tee, info ) 
    call dgbtrf ( n2, n2, kl, ku, wn_w_tee, 2*kl+ku+1, ipiv_w_tee, info ) 
    

    !----
    ! wn_op : v 
    !----

    do j=1,n2-1 
    do i= max(1,j-ku) , min(n2-1, j+kl) 

       if ( wvisc .eq. eddy_visc ) then 
          rre_v = c9r16* ( vt(ii,i+1,kk) + vt(ii,i  ,kk)) - c1r16* ( vt(ii,i-1,kk) + vt(ii,i+2,kk)) + nu_molec 
       endif

       if ( i .eq. j) then 
          wn_v_tee(ku+kl+1+i-j,j) = 1.d0 - dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       else 
          wn_v_tee(ku+kl+1+i-j,j) = -dt*beta* rre_v*lvb(ku+kl+1+i-j,j)
       endif
    enddo
    enddo

    !----
    ! operator from sgs gradient term 
    !----

    if ( wvisc .eq. eddy_visc ) then 
       
       do j=1,n2-1 
       do i=max(1,j-ku), min(n2-1,j+kl) 

          wn_v_tee(ku+kl+1+i-j,j) = wn_v_tee(ku+kl+1+i-j,j) - fsgs*dt*beta* evb_v(ku+kl+1+i-j,j) 
       enddo 
       enddo 
    endif


    call dgbtrf ( n2-1,n2-1, kl,ku, wn_v_tee, 2*kl+ku+1, ipiv_v_tee, info) 



  end subroutine wn_op_tee

!=================================================================================================!

end module numerics_tee
