      program particle_arrangement
      implicit none

! ********************************************************************
!
!
!      log:             2011 B. Vowinckel
!      purpose:         particle arrangement
!
! ********************************************************************

      integer::  i, j, k, n, m ,l, o, aum, p, q, a, b
      integer, parameter:: np = 13500
      double precision, dimension (np):: id, xc, yc, zc, distance
      double precision                  :: dummy




       open (unit=22, file='particle_data_rand.dat', action='read')
       do i = 1, np
          read(22,*) id(i), xc(i), yc(i), zc(i), dummy, dummy
       end do

       distance = 999.
       do i = 1, np
          do j = 1, np
          if (yc(i) .ne. yc(j)) then
              dummy       = sqrt((xc(i)-xc(j))**2 + &
                                  (yc(i)-yc(j))**2 + &
                                  (zc(i)-zc(j))**2)

              if (dummy.lt. distance(i) ) distance(i) = dummy 

          end if 
          end do

       end do

!        distance = distance * 12.

      write(*,*) minval(distance)/(1./.9)      , maxval(distance)/(1./.9)
         


      end program particle_arrangement
