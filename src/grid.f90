

module grid 

!--------------------------------------------------! 
!*    grid def                                    *!
!* y(i_u,1:n2)  : staggered mesh points for u,w   *! 
!* y(i_v,0:n2+1): staggered mesh points for v     *! 
!*                                                *! 
!* y(i_v,0)     : lower wall, -1.0d0              *! 
!* y(i_v,n2+1)  : upper wall,  1.0d0              *! 
!*                                                *! 
!*                                                *! 
!--------------------------------------------------! 

  use global 


contains 


!===========================================================================================!
  subroutine setup_grid ( itype ) 

    implicit none 
    integer , intent(in )           :: itype            !grid type 1,2,3 
    real                            :: hy0, xx, tgam 
    integer                         :: i,j,k,n



!... check itype 
    if ( itype .le. 0 .or. itype .ge. 4) then 
       write(*,*) 'error , grid type must be either 1,2, or 3 ... :', itype 
       stop 
    endif


!....
!.... uniform grid in y 
!....
    if ( itype .eq. 1) then 

       hy0      = 2.d0 / dble(n2) 
       y(i_u,0) = -1.d0               
       y(i_v,0) = -1.d0 

       do j=1,n2 
          hy(j)    = hy0  
          y(i_u,j) = -1.d0 + 0.5d0*hy0 + dble(j-1)*hy0
          y(i_v,j) = y(i_v,j-1) + hy(j) 
       enddo

       y(i_u,n2+1) = 1.d0 
       y(i_v,n2+1) = 1.d0

    endif


!....
!.... non-uniform grid (func. of cosine ) in y 
!....

    if ( itype .eq. 2) then 

       do j=1,n2 
          hy(j) = cos( pi* dble(j-1)/dble(n2) ) - cos(pi* dble(j)/dble(n2))
       enddo

       y(i_u,0) = -1.d0 
       y(i_u,1) = -1.d0 + 0.5d0*hy(1) 
       y(i_v,0) = -1.d0 

       do j=2,n2 
          y(i_u,j) = y(i_u,j-1) + 0.5d0* ( hy(j-1) + hy(j) )
       enddo
       
       do j=1,n2 
          y(i_v,j) = y(i_v,j-1) + hy(j) 
       enddo


       y(i_u,n2+1) = 1.0d0 
       y(i_v,n2+1) = 1.0d0 

    endif


!....
!.... non-uniform grid (func. of hyp-tangent ) in y 
!.... 
    
    if ( itype .eq. 3 ) then 

       hy0  = 2.0d0 / dble(n2) 
       tgam = tanh( gam ) 


       y(i_v,0) = -1.0d0 
       do j=1,n2-1 
          xx = hy0 * dble(j) 
          y(i_v,j) = -tanh(gam* (1.d0 - xx)) / tgam 
       enddo

       y(i_v,n2)   = 1.0d0 
       y(i_v,n2+1) = 1.0d0          ! should not be accessed 


       do j=1,n2 
          hy(j) = y(i_v,j) - y(i_v,j-1) 
          !y(i_u,j) = y(i_v,j-1) + 0.5d0* hy(j) 
          !Modified on 2018.08.23 by Danah
          !Instead of mid-point, use midpoint in uniform spacing
          xx = hy0 * (dble(j)-0.5d0) 
          y(i_u,j) = -tanh(gam* (1.d0 - xx)) / tgam 
       enddo


       y(i_u,0)    =-1.0d0 
       y(i_u,n2+1) = 1.0d0      


    endif



!.... symmetric spacings 

    hy(0)    = hy(1) 
    hy(-1)   = hy(2) 
    hy(-2)   = hy(3) 
    hy(n2+1) = hy(n2) 
    hy(n2+2) = hy(n2-1)         ! departues from original morinishi code
    hy(n2+3) = hy(n2-2) 


  end subroutine setup_grid
!==================================================================================================!

  subroutine setup_bc_coeff

    implicit none 
    real                            :: h1u, h1v, h2u, h2v, h3u 
    real                                   :: jac, djy, hy0

    hy0  = 2.d0/ dble(n2)
    jac  = jacb(i_u, 1 )
    djy  = hy0*jac

!.... boundary condition coefficients 

    cf = 0.d0 
    cf_slip = 0.d0 ! will be defined after calling restart file (time-dependent)

 if (rank .eq. 0) write(*,*) 'Setup_bc_coeff: No-slip for solid wall'
    h1u = y(i_u,1) - y(i_u,0) 
    h2u = y(i_u,2) - y(i_u,0) 
    
    h1v = y(i_v,1) - y(i_v,0) 
    h2v = y(i_v,2) - y(i_v,0) 

    cf(i_u, 0,1) = - ( h1u*h2u + h2u*h2u ) / (h2u*h2u - h1u*h2u) 
    cf(i_u, 0,2) =   2.d0*h1u*h1u          / (h2u*h2u - h1u*h2u) 
    cf(i_u,-1,1) = - 2.d0*h2u*h2u          / (h1u*h2u - h1u*h1u) 
    cf(i_u,-1,2) =   ( h1u*h2u + h1u*h1u ) / (h1u*h2u - h1u*h1u) 
    cf(i_u,-2,1) = - ( h2u*h3u*h3u + h2u*h2u*h3u ) / (h1u*h2u*h2u - h1u*h1u*h2u ) 
    cf(i_u,-2,2) =   ( h1u*h3u*h3u + h1u*h1u*h3u ) / (h1u*h2u*h2u - h1u*h1u*h2u ) 

    cf(i_w, 0,1) = - ( h1u*h2u + h2u*h2u ) / (h2u*h2u - h1u*h2u)
    cf(i_w, 0,2) =   2.d0*h1u*h1u          / (h2u*h2u - h1u*h2u)
    cf(i_w,-1,1) = - 2.d0*h2u*h2u          / (h1u*h2u - h1u*h1u)
    cf(i_w,-1,2) =   ( h1u*h2u + h1u*h1u ) / (h1u*h2u - h1u*h1u)
    cf(i_w,-2,1) = - ( h2u*h3u*h3u + h2u*h2u*h3u ) / (h1u*h2u*h2u - h1u*h1u*h2u)
    cf(i_w,-2,2) =   ( h1u*h3u*h3u + h1u*h1u*h3u ) / (h1u*h2u*h2u - h1u*h1u*h2u)
 
    cf(i_v,-1,1) = - ( h1v*h2v + h2v*h2v ) / (h2v*h2v - h1v*h2v)
    cf(i_v,-1,2) =   2.d0*h1v*h1v          / (h2v*h2v - h1v*h2v) 

    cf(i_v,-2,1) = - 2.d0*h2v*h2v          / (h1v*h2v - h1v*h1v)
    cf(i_v,-2,2) =   ( h1v*h2v + h1v*h1v ) / (h1v*h2v - h1v*h1v) 
       
  end subroutine setup_bc_coeff
  
!============================================================================================!


  subroutine compute_jacobian ( igrid )               

!.... computes the exact jacobian of the specified grid 

    implicit none 
    integer                :: igrid 
    integer                :: i,j 
    real                   :: hy0 , xx


    hy0  = 2.d0 / dble(n2 ) 
    

    if ( igrid .eq. 3 ) then
       do j=0,n2 
       
          !xx           = atanh ( tanh(gam)* y(i_v,j) )/ gam 
          xx           = dble(j) * hy0 - 1.d0 
          jacb(i_v,j ) = gam/ tanh(gam) / cosh(gam*xx)/ cosh(gam*xx)
          
       enddo 


       !do j=1,n2 
       do j=0,n2 
       
          !xx           = atanh ( tanh(gam)* y(i_u,j) )/ gam 
          xx           = (dble(j) - 0.5d0 )* hy0 - 1.d0 
          jacb(i_u,j ) = gam/ tanh(gam) / cosh(gam*xx)/ cosh(gam*xx) 
          
       enddo 

    else

       !.... only hyp tan grid supported currently 

       do j=0,n2 
          jacb(i_u,j) = 1.d0 
          jacb(i_v,j) = 1.d0 
       enddo

    endif

    !Modified 2018.08.23 by Danah Park
    !Overide Jacobian's with real delta y's
    !Keep jacobians at j=0,n2 for i_v and j=0,n2+1 for i_u same
    do j=1,n2-1
        jacb( i_v,j ) = ( y(i_u,j+1) - y(i_u,j) ) * n2 / 2.d0
    enddo

    do j=1,n2
        jacb( i_u,j ) = ( y(i_v,j) - y(i_v,j-1) ) * n2 / 2.d0
    enddo

  end subroutine compute_jacobian 

!========================================================================================================================!

  subroutine write_grid 

!.... write tecplot file for the grid into output/grid_3d.tec 
!.... use i_v nodes in the wall normal direction 

    implicit none 
    integer                               :: i,j,k 
    real                                  :: dx, dz 
    real, dimension(:), allocatable       :: x,z 


    allocate( x(0:n1-1), z(0:n3-1)) 


!.... output full 3d grid 

    open(unit=2, file='output/grid_3d.tec', status='unknown') 

    write(2,*) 'TITLE= CHANNEL FLOW GRID' 
    write(2,*) 'VARIABLES = "x" "y" "z"' 
    write(2,*) 'ZONE DATAPACKING=POINT, I= ', n1, ', J= ', n2+2, ', K= ', n3 
        

    
    dx = chlx * pi / dble(n1) 
    dz = chlz * pi / dble(n3) 
    x(0) = 0.0d0 
    z(0) = 0.0d0  

    do i=1,n1-1 
       x(i) = x(i-1) + dx 
    enddo

    do k=1,n3-1 
       z(k) = z(k-1) + dz
    enddo

    do k=0,n3-1 
    do j=0,n2+1
    do i=0,n1-1 


       
       write(2,'(3d16.8)') x(i), y(i_v,j) , z(k)  
    enddo
    enddo 
    enddo 


    close(unit=2) 

!.... output y2 slice 
!.... file output for gnuplot, matlab

    open(unit=3, file='output/grid_wn.dat', status='unknown') 

    do j=0,n2 
       write(3,'(2f16.8)') 1.0d0, y(i_v,j) 
    enddo

    close(unit=3) 


!    write(*,*) 'x_e = ', x(n1-1), ', z_e = ', z(n3-1), ', dx = ', dx , ', dz = ', dz  

    deallocate(x,z) 

  end subroutine write_grid

!=====================================================================================================! 

  subroutine load_grid 
    
    ! load file from input/grid_u.dat and input/grid_v.dat 
    ! 

    implicit none 
    logical       :: uexst , vexst 
    integer       :: n2u_in, n2v_in 
    integer       :: j


    ! check that the input grid files exist 
    inquire ( file='input/grid_u.dat', exist=uexst) 
    inquire ( file='input/grid_v.dat', exist=vexst) 

    if ( uexst .and. vexst ) then 
       
       ! load u grid 
       open(unit=450+rank, file='input/grid_u.dat', status='old')  
       read(450+rank,*) n2u_in 

       if ( n2u_in .eq. n2 ) then 
          
          ! set grid endpoints for the channel walls 
          y(i_u,0   ) = -1.d0 
          y(i_u,n2+1) =  1.d0 

          do j=1,n2 
             read(450+rank, *) y(i_u,j) 
          enddo
       else  
          
          if ( rank .eq. 0) then 
             write(*,*) 'load_grid:: incorrect number of  u points ...', n2u_in, n2 
             write(*,*) 'load_grid:: forcing exit...' 
          endif 
          
          call mem_dealloc 
          call mpi_finalize ( ierr ) 
       endif

       close(450+rank) 

       ! load v grid 
       open(unit=450+rank, file='input/grid_v.dat', status='old') 
       read(450+rank,*) n2v_in 

       if ( n2v_in .eq. n2-1 ) then 
          
          ! set grid endpoints for the channel walls 
          y(i_v,0)    = -1.d0 
          y(i_v,n2)   =  1.d0 
          y(i_v,n2+1) =  1.d0   ! should not be accessed 

          do j=1,n2-1
             read(450+rank,*) y(i_v,j) 
          enddo

       else  
          
          if ( rank .eq. 0 ) then 
             write(*,*) 'load_grid:: incorrect number of v points ...', n2v_in, n2-1 
             write(*,*) 'load_grid:: forcing exit...'
          endif 

          call mem_dealloc 
          call mpi_finalize ( ierr ) 
       endif 

       close(450+rank)

       ! init hy(1:n2) ---> used in the cfl computation 
       do j=1,n2 
          hy(j) = y(i_v,j) - y(i_v,j-1) 
       enddo

    else  

       if ( rank .eq. 0) then 
          write(*,*) 'load_grid:: necessary u and v grid input files are not present!!' 
          write(*,*) 'load_grid:: forcing exit...' 
       endif 

       call mem_dealloc 
       call mpi_finalize ( ierr ) 
    endif 



  end subroutine load_grid 
!======================================================================================================!

  subroutine num_jacobian 

    ! currently only computes a second order accurate jacobian 
    ! cf, vasilyev, jcp, 2000 

    implicit none 
    integer          :: j 
    real             :: hy0 

    hy0 = 2.d0 / dble(n2) 
    
    ! ... u jacobian
    do j=1,n2 
       jacb(i_u,j) = ( y(i_u,j+1) - y(i_u,j-1) ) / (2.d0 * hy0)
    enddo 


    !.... v jacobian 
    do j=1,n2-1 
       jacb(i_v,j) = ( y(i_v,j+1) - y(i_v,j-1) ) / (2.d0 * hy0) 
    enddo 


  end subroutine num_jacobian
!======================================================================================================!
  
end module grid
