
module comm_routines 
  
  !* ---- 
  !* communication routines for scalar cell-centered, 
  !* staggered vectors, and cell-centered tensorial variables 
  !* 
  !* stb, 6/23/2009 
  !* ----- 
  
  
  use global 
  
  
  implicit none 
  integer*8          :: vec_comm, scalar_comm, sym_tensor_comm, & 
                        full_tensor_comm 
  integer*8          :: scalar_comm_yh
 
  integer, parameter :: nsym_tensor = 6 
  integer, parameter :: n_tensor = 9 

contains 

!=====================================================================================================!

  subroutine build_communicator_types 

    ! for these types, there is no halo region in the 
    ! wall normal direction 

    implicit none 
    integer        :: all_y , all_z
    
    all_y = n2 
    all_z = nz + 2*halo_z 

    
    call mpi_type_vector ( all_y* all_z, nsym_tensor* halo_x, nsym_tensor*(nx+2*halo_x ), & 
                           MPI_DOUBLE_PRECISION, sym_tensor_comm, ierr ) 

    call mpi_type_vector ( all_y* all_z, n_tensor* halo_x , n_tensor*(nx+2*halo_x ) , & 
                           MPI_DOUBLE_PRECISION, full_tensor_comm, ierr ) 


    call mpi_type_vector ( all_y* all_z, nvel*halo_x, nvel*(nx+2*halo_x ), & 
                           MPI_DOUBLE_PRECISION, vec_comm, ierr) 
    call mpi_type_vector ( all_y* all_z,      halo_x,      (nx+2*halo_x ), & 
                           MPI_DOUBLE_PRECISION, scalar_comm, ierr )


    call mpi_type_vector ( (n2 + 2*halo_yu)* all_z, halo_x, (nx+2*halo_x), & 
                           MPI_DOUBLE_PRECISION, scalar_comm_yh, ierr ) 




    call mpi_type_commit ( scalar_comm, ierr ) 
    call mpi_type_commit ( scalar_comm_yh, ierr) 
    call mpi_type_commit ( vec_comm   , ierr ) 
    call mpi_type_commit ( sym_tensor_comm, ierr ) 
    call mpi_type_commit ( full_tensor_comm,ierr ) 


  end subroutine build_communicator_types 
!=====================================================================================================!

  subroutine scalar_commnctor ( pp ) 

    implicit none 
    integer                             :: i, j,k, jss 
    integer                             :: ierr 
    integer, dimension(mpi_status_size) :: status 
    real, dimension(is-halo_x:ie+halo_x, 1:n2, ks-halo_z:ke+halo_z) :: pp 
  


    span_block = (nx+2*halo_x)* n2 * halo_z 
    jss        = 1 


    ! left <-> right communicator 


    if ( mod(coords(0), 2) .eq. 0 ) then 
       call mpi_sendrecv ( pp(ie-halo_x+1,jss,ks-halo_z), 1, scalar_comm, right, 0, & 
                           pp(ie+1       ,jss,ks-halo_z), 1, scalar_comm, right, 0, comm3d, status, ierr) 

       call mpi_sendrecv ( pp(is         ,jss,ks-halo_z), 1, scalar_comm, left , 0, & 
                           pp(is-halo_x  ,jss,ks-halo_z), 1, scalar_comm, left , 0, comm3d, status, ierr) 

    else  
       call mpi_sendrecv ( pp(is         ,jss,ks-halo_z), 1, scalar_comm, left , 0, & 
                           pp(is-halo_x  ,jss,ks-halo_z), 1, scalar_comm, left , 0, comm3d, status, ierr) 

       call mpi_sendrecv ( pp(ie-halo_x+1,jss,ks-halo_z), 1, scalar_comm, right, 0, & 
                           pp(ie+1       ,jss,ks-halo_z), 1, scalar_comm, right, 0, comm3d, status, ierr)
    endif 



    ! front <-> back communicator 
    
    if ( mod(coords(2), 2) .eq. 0 ) then 
       call mpi_sendrecv ( pp(is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           pp(is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
       
       call mpi_sendrecv ( pp(is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           pp(is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr) 
    else 
       call mpi_sendrecv ( pp(is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           pp(is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr) 
       
       call mpi_sendrecv ( pp(is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           pp(is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
    endif 

  end subroutine scalar_commnctor
 
!========================================================================================================================!

  subroutine scalar_commnctor_yh ( pp) 


    implicit none 
    integer                             :: i, j,k, jss 
    integer                             :: ierr 
    integer, dimension(mpi_status_size) :: status 
    real, dimension(is-halo_x:ie+halo_x, 1-halo_yu:n2+halo_yu, ks-halo_z:ke+halo_z) :: pp 
  


    span_block = (nx+2*halo_x)* (n2+ 2*halo_yu) * halo_z 
    jss        = 1-halo_yu 


    ! left <-> right communicator 


    if ( mod(coords(0), 2) .eq. 0 ) then 
       call mpi_sendrecv ( pp(ie-halo_x+1,jss,ks-halo_z), 1, scalar_comm_yh, right, 0, & 
                           pp(ie+1       ,jss,ks-halo_z), 1, scalar_comm_yh, right, 0, comm3d, status, ierr) 

       call mpi_sendrecv ( pp(is         ,jss,ks-halo_z), 1, scalar_comm_yh, left , 0, & 
                           pp(is-halo_x  ,jss,ks-halo_z), 1, scalar_comm_yh, left , 0, comm3d, status, ierr) 

    else  
       call mpi_sendrecv ( pp(is         ,jss,ks-halo_z), 1, scalar_comm_yh, left , 0, & 
                           pp(is-halo_x  ,jss,ks-halo_z), 1, scalar_comm_yh, left , 0, comm3d, status, ierr) 

       call mpi_sendrecv ( pp(ie-halo_x+1,jss,ks-halo_z), 1, scalar_comm_yh, right, 0, & 
                           pp(ie+1       ,jss,ks-halo_z), 1, scalar_comm_yh, right, 0, comm3d, status, ierr)
    endif 



    ! front <-> back communicator 
    
    if ( mod(coords(2), 2) .eq. 0 ) then 
       call mpi_sendrecv ( pp(is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           pp(is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
       
       call mpi_sendrecv ( pp(is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           pp(is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr) 
    else 
       call mpi_sendrecv ( pp(is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           pp(is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr) 
       
       call mpi_sendrecv ( pp(is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           pp(is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
    endif 


  end subroutine scalar_commnctor_yh
  !============================================================================================================================================!


  subroutine stag_vector_commnctor ( qq) 

    implicit none 
    integer                             :: i,j,k,i_q, jss 
    integer, dimension(mpi_status_size) :: status 
    real, dimension(3, is-halo_x:ie+halo_x,1:n2,ks-halo_z:ke+halo_z) :: qq 
    


    span_block = nvel * (nx+2*halo_x)* n2* halo_z 
    jss = 1 



    ! left <-> right communicator 

    
    if ( mod( coords(0), 2) .eq. 0 ) then 
       call mpi_sendrecv ( qq(i_u,ie-halo_x+1,jss,ks-halo_z), 1, vec_comm, right, 0, & 
                           qq(i_u,ie+1       ,jss,ks-halo_z), 1, vec_comm, right, 0, comm3d, status, ierr)

       call mpi_sendrecv ( qq(i_u,is         ,jss,ks-halo_z), 1, vec_comm, left , 0, & 
                           qq(i_u,is-halo_x  ,jss,ks-halo_z), 1, vec_comm, left , 0, comm3d, status, ierr) 

    else  

       call mpi_sendrecv ( qq(i_u,is         ,jss,ks-halo_z), 1, vec_comm, left , 0, & 
                           qq(i_u,is-halo_x  ,jss,ks-halo_z), 1, vec_comm, left , 0, comm3d, status, ierr) 

       call mpi_sendrecv ( qq(i_u,ie-halo_x+1,jss,ks-halo_z), 1, vec_comm, right, 0, & 
                           qq(i_u,ie+1       ,jss,ks-halo_z), 1, vec_comm, right, 0, comm3d, status, ierr) 
    endif


    ! front <-> back communicator 
    

    if ( mod( coords(2), 2) .eq. 0 ) then 
       call mpi_sendrecv ( qq(i_u,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           qq(i_u,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr)
          
       call mpi_sendrecv ( qq(i_u,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           qq(i_u,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
    else  
       call mpi_sendrecv ( qq(i_u,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           qq(i_u,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
       call mpi_sendrecv ( qq(i_u,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           qq(i_u,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
    endif


  end subroutine stag_vector_commnctor

!=======================================================================================================================!

  subroutine sym_tensor_commnctor ( qij ) 

    implicit none 
    integer                             :: i,j,k,i_q, jss 
    integer, dimension(mpi_status_size) :: status 
    real, dimension(nsym_tensor, is-halo_x:ie+halo_x, 1:n2, ks-halo_z:ke+halo_z) :: qij 


    
    span_block = nsym_tensor * (nx+2*halo_x)* n2* halo_z 
    jss = 1 
    


    ! left <-> right communicator 

    
    if ( mod( coords(0), 2) .eq. 0 ) then 
       call mpi_sendrecv ( qij(1,ie-halo_x+1,jss,ks-halo_z), 1, sym_tensor_comm, right, 0, & 
                           qij(1,ie+1       ,jss,ks-halo_z), 1, sym_tensor_comm, right, 0, comm3d, status, ierr)

       call mpi_sendrecv ( qij(1,is         ,jss,ks-halo_z), 1, sym_tensor_comm, left , 0, & 
                           qij(1,is-halo_x  ,jss,ks-halo_z), 1, sym_tensor_comm, left , 0, comm3d, status, ierr) 

    else  

       call mpi_sendrecv ( qij(1,is         ,jss,ks-halo_z), 1, sym_tensor_comm, left , 0, & 
                           qij(1,is-halo_x  ,jss,ks-halo_z), 1, sym_tensor_comm, left , 0, comm3d, status, ierr) 

       call mpi_sendrecv ( qij(1,ie-halo_x+1,jss,ks-halo_z), 1, sym_tensor_comm, right, 0, & 
                           qij(1,ie+1       ,jss,ks-halo_z), 1, sym_tensor_comm, right, 0, comm3d, status, ierr) 
    endif


    ! front <-> back communicator 
    

    if ( mod( coords(2), 2) .eq. 0 ) then 
       call mpi_sendrecv ( qij(1,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           qij(1,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr)
          
       call mpi_sendrecv ( qij(1,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           qij(1,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
    else  
       call mpi_sendrecv ( qij(1,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           qij(1,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
       call mpi_sendrecv ( qij(1,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           qij(1,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
    endif


  end subroutine sym_tensor_commnctor

!=======================================================================================================================!

   subroutine full_tensor_commnctor ( qij ) 

    implicit none 
    integer                             :: i,j,k,i_q, jss 
    integer, dimension(mpi_status_size) :: status 
    real, dimension(n_tensor, is-halo_x:ie+halo_x, 1:n2, ks-halo_z:ke+halo_z) :: qij 


    
    span_block = n_tensor * (nx+2*halo_x)* n2* halo_z 
    jss = 1 
    


    ! left <-> right communicator 

    
    if ( mod( coords(0), 2) .eq. 0 ) then 
       call mpi_sendrecv ( qij(1,ie-halo_x+1,jss,ks-halo_z), 1, full_tensor_comm, right, 0, & 
                           qij(1,ie+1       ,jss,ks-halo_z), 1, full_tensor_comm, right, 0, comm3d, status, ierr)

       call mpi_sendrecv ( qij(1,is         ,jss,ks-halo_z), 1, full_tensor_comm, left , 0, & 
                           qij(1,is-halo_x  ,jss,ks-halo_z), 1, full_tensor_comm, left , 0, comm3d, status, ierr) 

    else  

       call mpi_sendrecv ( qij(1,is         ,jss,ks-halo_z), 1, full_tensor_comm, left , 0, & 
                           qij(1,is-halo_x  ,jss,ks-halo_z), 1, full_tensor_comm, left , 0, comm3d, status, ierr) 

       call mpi_sendrecv ( qij(1,ie-halo_x+1,jss,ks-halo_z), 1, full_tensor_comm, right, 0, & 
                           qij(1,ie+1       ,jss,ks-halo_z), 1, full_tensor_comm, right, 0, comm3d, status, ierr) 
    endif


    ! front <-> back communicator 
    

    if ( mod( coords(2), 2) .eq. 0 ) then 
       call mpi_sendrecv ( qij(1,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           qij(1,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr)
          
       call mpi_sendrecv ( qij(1,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           qij(1,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
    else  
       call mpi_sendrecv ( qij(1,is-halo_x,jss,ks         ), span_block, MPI_DOUBLE_PRECISION, back, 0, & 
                           qij(1,is-halo_x,jss,ks-halo_z  ), span_block, MPI_DOUBLE_PRECISION, back, 0, comm3d, status, ierr)
          
       call mpi_sendrecv ( qij(1,is-halo_x,jss,ke-halo_z+1), span_block, MPI_DOUBLE_PRECISION, front , 0, & 
                           qij(1,is-halo_x,jss,ke+1       ), span_block, MPI_DOUBLE_PRECISION, front , 0, comm3d, status, ierr) 
    endif


  end subroutine full_tensor_commnctor

!=======================================================================================================================!

  subroutine test_commnctors

    implicit none 
    integer                   :: i,j,k,n
    integer                   :: istr, iend, kstr, kend
    real, dimension(:,:,:,:), allocatable  :: test_ij, test_svec 
    real, dimension(  :,:,:), allocatable  :: test_p 


    kstr = ks - halo_z 
    kend = ke + halo_z 
    istr = is - halo_x 
    iend = ie + halo_x 


    allocate ( test_p(istr:iend,1:n2,kstr:kend)) 
    allocate ( test_svec(nvel,istr:iend,1:n2,kstr:kend)) 
    allocate ( test_ij(nsym_tensor,istr:iend,1:n2,kstr:kend)) 


    test_p    = 0.d0 
    test_svec = 0.d0 
    test_ij   = 0.d0 


    ! init the regions that each processor owns 
    do k=ks,ke 
    do j=1,n2 
    do i=is,ie 

       test_p(i,j,k) = j 

       do n=1,nvel 
          test_svec(n,i,j,k) = n+j 
       enddo 


       do n=1,nsym_tensor 
          test_ij(n,i,j,k) = n+j 
       enddo 
     

    enddo
    enddo 
    enddo 



    ! now sync the tensors
    call scalar_commnctor( test_p ) 
    call stag_vector_commnctor ( test_svec ) 
    call sym_tensor_commnctor ( test_ij ) 


    ! check the halo regions 
    
    do k=ks-halo_z,ks 
    do j=1,n2 
    do i=is-halo_x,ie+halo_x 

       if ( test_p(i,j,k) .ne. (j)) write(*,*) 'scalar = ', rank, i,j,k 
       
       do n=1,nvel 
          if ( test_svec(n,i,j,k) .ne. (n+j)) write(*,*) 'vec = ', rank, n,i,j,k 
       enddo 


       do n=1,nsym_tensor 
          if ( test_ij(n,i,j,k) .ne. (n+j)) write(*,*) 'stens = ', rank, n, i,j,k 
       enddo
    
    enddo
    enddo 
    enddo 



    do k=ke,ke+halo_z
    do j=1,n2 
    do i=is-halo_x,ie+halo_x 

       if ( test_p(i,j,k) .ne. j) write(*,*) 'scalar = ', rank, i,j,k 
       
       do n=1,nvel 
          if ( test_svec(n,i,j,k) .ne. (n+j)) write(*,*) 'vec = ', rank, n, i,j,k 
       enddo 

       do n=1,nsym_tensor 
          if ( test_ij(n,i,j,k) .ne. (n+j)) write(*,*) 'stens = ', rank, n,i,j,k 
       enddo 

    enddo 
    enddo 
    enddo 



    deallocate ( test_p) 
    deallocate ( test_svec ) 
    deallocate ( test_ij ) 


  end subroutine test_commnctors 
!=======================================================================================================================!


end module comm_routines



    
    
