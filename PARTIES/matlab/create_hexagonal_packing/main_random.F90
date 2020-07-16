      program arrange

! *******************************************************************
! 
!     purpose:          Dat-file mit regelmäßiger blasenanordnung erzeugen
!
!     remarks:	
! 
!     log:              2010 / 04  s.heitkam
! 
! *******************************************************************

    implicit none
 

  integer ::   n, i, j, j2 ,k, n_x, n_y, n_z, try_max, seed, n_now, cnt

  logical            ::   f_struct, fail

  double precision :: 	len_x, len_y, len_z, dis, x_min, dx, y_min, dy, z_min, dz,r, x1, y1, z1, &
                        off_oben, off_unten
			
			


     double precision, ALLOCATABLE, DIMENSION(:) :: x, y, z           
 
   character*60          ::  name
   character (len = 100) :: str
   r=1./18.
!  dis=1.2d0/5.d0*r   !abstand der kugelOBERFLÄCHEN
   dis=1./18.

   len_x=24
   len_y=1.2
   len_z=6
   off_oben=4.d0*r!+5.d0*dis !am drainage-eintritt, also bei kleinen y
   off_unten=2.d0*r ! am drainageaustritt, also bei großen y
   n= 13500 !min(334,floor(len_x*len_y*len_z/((2.*r+dis)**3)))
    ALLOCATE(x(n))
    ALLOCATE(y(n))
    ALLOCATE(z(n))
   write(*,*) 'number = ', n, '    phi = ', n*4./3.*3.1415*(r+dis)**3/(len_x*len_y*len_z)
   n_now=0
   try_max=800*n
   seed=7654321
   call random_seed()
   do j=1,1200000
       call random_number(dx)
   end do
   
   do j=1,try_max
 !     dy=ran(seed)
      fail=.true.
      call random_number(dx)
      call random_number(dy)  
      call random_number(dz)  
      x1=len_x*dx
      y1=len_y*dy
      z1=len_z*dz 
      if (x1.lt.(r+dis))                    fail=.false. 
      if (y1.lt.(r+dis)+off_oben)           fail=.false.
      if (z1.lt.(r+dis))                    fail=.false.
      if (x1.gt.(len_x-(r+dis)))            fail=.false.
      if (y1.gt.(len_y-(r+dis)-off_unten))  fail=.false.
      if (z1.gt.(len_z-(r+dis)))            fail=.false.
      do cnt=1,n_now
       if ((sqrt((x1-x(cnt))**2+(y1-y(cnt))**2+(z1-z(cnt))**2)).lt.(2.*r+dis))  fail=.false.
      end do
      if (fail) then
        n_now=n_now+1  
        x(n_now)=x1
        y(n_now)=y1
        z(n_now)=z1
        write(*,*) j,n_now,x1,y1,z1
      end if 
      if (n_now.eq.n) goto 100 

   end do
      
 100  continue

!   end do
! 
!    do i=1,n
! ! 
!     y(i)=y(i)+9.d0*r
! !  
!    end do
!        do k=1,n_z
! 	do j=1,n_y
! 
! 
! 
       
      open (unit = 11, file = 'particle_data_rand.dat')
    do j=1,n

        y1=minval(y(:))
        do i=1,n
           if (y(i).eq.y1) then
!                write(11,*) j+6696, sngl(x(i)), sngl(y(i)), sngl(z(i)), sngl(r), 12.d0, 0.5, 0.5, 0.5
               write(11,*) j+20379, x(i), y(i), z(i), r, 1200.0
               y(i)=y(i)+len_y
               goto 200
           end if
        end do
 200  continue
    end do
 !     write(11,*)
      close(11)



end program arrange