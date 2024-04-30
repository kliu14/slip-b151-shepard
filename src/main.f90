
program turbulent_interface

!*--------------------------------------------------------------------*! 
!*        main program for turbulent channel flow / IMFM code         *! 
!*--------------------------------------------------------------------*! 
  use global 
  use grid 
  use numerics 
  use numerics_tee
  use output 
  use comm_routines
  
  implicit none 
  integer                        :: i,istep
  real                           :: startt, startt2, endt
  real                           :: cpu_diff, cpu_time_start, cpu_time_end
  real                           :: cpu_diff_total, cpu_time_start_total, cpu_time_end_total
  real                           :: k
  logical :: perturb_c, debug
  character(len=3)                        :: fstr, time_str_below_p
  character(len=4)                        :: print_str, time_str
!*--------------------------------------------------------------------*! 
  !---
  ! initialize variables , mpi 
  !---

  call mpi_init ( ierr ) 
  call mpi_comm_dup ( mpi_comm_world, mycomm, ierr ) 
  call mpi_comm_rank( mycomm, rank, ierr ) 
  call mpi_comm_size( mycomm, size, ierr ) 

!!===============================================================
if (rank .eq. 0 .and. log_flag) write(*,*) " 0. Read input, domain decomposition, mem allocation"

  call read_input 
  call domain_decomp 
  call mem_alloc
  call build_communicator_types

!!===============================================================
if (rank .eq. 0 .and. log_flag) write(*,*) " 1. Grid, bc, rkc, wn_op setup   "

  if ( lgrid ) then 
     call load_grid 
     call num_jacobian 
  else  
     call setup_grid ( igrid ) 
     call compute_jacobian ( igrid ) 
  endif

  if ( pat_bc ) call setup_pat ! (Kim: 09.14.22)

  if ( rank .eq. 0 .and. pat_bc .and. IMFMswitch .and. itimestep .ne. 3 ) then
    write(*,*) "ERROR: patterned BC with IMFM must be AB2 only"
  endif

  call setup_bc_coeff

  if (rank .eq. 0 .and. log_flag) write(*,*) "1-1. Setup rkc"
  call setup_rkc 
  call compute_lub_2nd
  if ( transportee ) call compute_lub_tee ! (Kim: 09.13.22)
  if ( pat_bc) call compute_lub_pat
  call compute_lwb_2nd
  call compute_lvb_2nd
  call compute_pwb 
  call compute_pwb_2nd
  call compute_pwb_2nd_Dirichlet
  call setup_output
  firstloop = .true.
  flat_EE      = .false. 
!!===============================================================
if (rank .eq. 0 .and. log_flag) write(*,*) "2. Read restart file or initialization "
  call setup_init_restart
  if ( slip_bc ) then
    call update_slip ! set time-dependent slip conditions
  endif
!===============================================================
if (rank .eq. 0 .and. log_flag) write(*,*) "3. Set output statstics"
  call diag_out
  startt = mpi_wtime()
  call compute_stats
  if ( ostat) call output_stat(1)
     int_time = int(time) 
     postprocess_time_below_p = (int(int((time - int(time))*1000)/(postprocess_time_interval))+1)*postprocess_time_interval   
     if (postprocess_time_below_p .ge. 1000) then
     postprocess_time_below_p = 0
     int_time = int_time + postprocess_time
     endif
     if ( rank .eq. 0 .and. log_flag .and. postprocess_time_below_p .eq. 0 ) write(*,*) " 4-1. Output flow field in every = ",postprocess_time, ".0  eddy turnover time "
     if ( rank .eq. 0 .and. log_flag .and. postprocess_time_below_p .ne. 0 ) write(*,*) " 4-1. Output flow field in every = 0.",postprocess_time_interval, " eddy turnover time "
!!===============================================================
 if (rank .eq. 0 .and. log_flag) write(*,*) "4. Diag output:Step, Time, CFL, u_tau(u,l)"
 if (rank .eq. 0 .and. log_flag) write(*,*) "5. Begin time stepping"
 if (rank .eq. 0 .and. log_flag .and. transportee) write(*,*) "5-1. Generalizd NS equations for transportee variable is on"
 if (rank .eq. 0 .and. log_flag .and. IMFMswitch) write(*,*) "5-2. IMFM is on"

  do i=1,imax 

     time = time + dt 
     if (i .eq. 1) ABstep = 1
     if (i .eq. 2) ABstep = 2
     if (i .ge. 3) ABstep = 3

!.... Modules to apply IMFM
!.... Modified by Danah Park on May 15, 2018 (last update: Dec 12, 2018)
!.... at each step 'before' calculating u, calculate transportee variable
     if (transportee) then
       do istep=1,nrksteps 
          if (p_order .eq. 4) then 
            if (rank .eq. 0) write(*,*) "4th order rkstep is not implemented for IMFM"
          else
            if (IMFMswitch) then
              call rk_step2nd_IMFM(istep)
            else
              call rk_step2nd_tee(istep)
            endif
          endif
          call pressure_step_tee(istep)
          call prjsln_tee  ( istep )            ! project onto solen field 
          call bcuvw_tee                        ! apply bc 
          call update_pressure_tee 
          !flat_EE = .false. 
          !i_print = i_print + 1
       enddo 
     endif
!.... module ends here

     do istep=1,nrksteps 
        if (p_order .eq. 4) then 
        call rk_step(istep)
        else
        call rk_step2nd(istep)
        endif
        call pressure_step(istep)
        call prjsln  ( istep )            ! project onto solen field 
        call bcuvw                        ! apply bc 
        call update_pressure 

        flat_EE = .false. 
        i_print = i_print + 1
     enddo 

     call diag_out 
     if ( mod(i,ist) .eq. 0 .and. time .gt. sstat .and. i .gt. 20 ) then 
        call compute_stats 
     endif 

     if ( mod(i,isam) .eq. 0 ) then 
       write(fstr, '(i3.3)')  n_interface
       if (output_qcrit ) call qcriterion
       if (output_qcrit ) call dump_qcrit
       if ( mod(i,(isam)) .eq. 0) then
          if ( ostat) call output_stat(1) 
        endif
         n_interface = n_interface + 1
        if (n_interface .ge. 1000 .and. output_vort .eq. .false.)  n_interface = 0
        if (n_interface .ge. 10000 .and. output_vort .eq. .true.)  n_interface = 0

     endif 

     if ( mod(i, irestart) .eq. 0 .and. itimestep .ne. 2) then 
            ab_flag = .true. 
            call dump_restart_file_interface('output/restart/restart.out')
            !... addition output module for transportee variable and IMFM
            if (transportee) call dump_restart_file_interface_transportee('output/restart/restart_transportee.out')
            if (IMFMswitch) call dump_restart_file_interface_IMFM('output/restart/restart_IMFM.out')
     endif

     if ( int_time .eq. int(time) .and. postprocess_time_below_p .eq. int((time- int(time))*1000.d0)  ) then
          if (ofiel) then
                 write(time_str, '(i4.4)')  int_time
                 write(time_str_below_p, '(i3.3)')  postprocess_time_below_p
             ab_flag = .false.
            call dump_restart_file_interface('output/flow_field/field_'//time_str//'_'//time_str_below_p//'.dat')
            !... addition output module for transportee variable and IMFM
            if (transportee) call dump_restart_file_interface_transportee('output/flow_field_transportee/field_'//time_str//'_'//time_str_below_p//'.dat')
            if (IMFMswitch) call dump_restart_file_interface_IMFM('output/flow_field_IMFM/field_'//time_str//'_'//time_str_below_p//'.dat')
         endif
     postprocess_time_below_p = postprocess_time_below_p + postprocess_time_interval
     if (postprocess_time_below_p .ge. 1000) then 
     postprocess_time_below_p = 0
     int_time = int_time + postprocess_time
     endif
     endif

     if ( time .ge. finish_time) exit

  enddo
!!===============================================================
! End time stepping
!!===============================================================
  if ( rank .eq. 0 .and. log_flag) write(*,*) '6. Finish simulation' 
  endt = mpi_wtime() 
  if ( rank .eq. 0 .and. log_flag) write(*,*) 'elapsed time..', endt-startt

  if (ostat ) call output_stat 

  call destroy_fft
  call mem_dealloc    
  call mpi_finalize ( ierr ) 


end program turbulent_interface
