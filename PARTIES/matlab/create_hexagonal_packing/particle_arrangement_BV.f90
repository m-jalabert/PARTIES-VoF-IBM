      program particle_arrangement
      implicit none

! ********************************************************************
!
!
!      log:             2011 B. Vowinckel
!      purpose:         particle arrangement
!
! ********************************************************************

      integer::  i, j, k, n, m ,l, o, aum, p, q, a, b, id, nx, nz, z , y, total_num
      double precision:: d, rho, l_x, l_y, l_z, dy, dz, ix, gx, gz, breite, zk, numrow, zk_offset, ix_offset, j_offset, d_offset



! Einlesen der Partikelanzahl, Durchmesser und Dichte

!Partikel in Y-Richtung:
       j = 2
       j_offset = 0.025  !! in y direction

!Durchmesser:
       d = 0.025-1e-9 !1 
	d_offset=d*0.01  !! diameter is reduced by this value
!Längen/ Breite des Feldes:
       ix=   0.2-sqrt(6.)/3.*d  !! that is Lz  == ix+ d
       zk=   0.1  !! that is lX
       zk_offset = sqrt(6.)/3.*d
	ix_offset = 0.  !! zmin =0 zmax= ix+d
	
!Dichte:
       rho = 99

! Erstellen der Partikel-Datei (Anordnung der Partikel)
       open (unit=22, file='particle_data_test.dat', action='write') !, status='replace'

       aum =2        !Auswahlparameter (0= Kubische Anodnung; 1= Quadratische Anordnung; 2= Hexagonale Anordnung)

!bei der Ausgabe wird der Durchmesser durch 2 geteilt, da für das folgende Programm der
!Radius bötigt wird (für die Berechnungen wird aber der Durchmesser verwendet)

       auswahl: select case (aum)

       case (0) auswahl

       gx = ix/(d)
       i = int(gx)

       gz = zk/(d)
       k = int(gz)

       l_x=-d
       do n = 1, i
       l_x = l_x+d
         l_z=-d
         do m = 1, k
         l_z = l_z+d
           l_y=-d
           do l = 1, j
           l_y = l_y+d
                      id = id+1
           write(22,*) id,  l_x, l_y, l_z, d/2, rho
           end do
         end do
       end do

!Breitenangabe:
       breite =l_z+d/2
       write(*,*) 'Maximum in z:'
       write(*,*) breite

       case (1) auswahl

       gx = ix/(d)
       i = int(gx)

       gz = zk/(d)
       k = int(gz)

       dy = real(j)
       dy = dy/2

       a = nint(dy)
       b = int(dy)

       l_x=-d
       do n = 1, i
       l_x = l_x+d
         l_z=-d
         do m = 1, k
         l_z = l_z+d
           l_y=-sqrt(2.0)*d
           do l = 1, a
           l_y = l_y+ sqrt(2.0)*d
                      id = id+1
           write(22,*) id,  l_x, l_y, l_z, d/2, rho
           end do
         end do
       end do

       l_x=-d/2
       do n = 1, i
       l_x = l_x+d
         l_z=-d/2
         do m = 1, k
         l_z = l_z+d
           l_y= -sqrt(2.0)*d/2
           do l = 1, b
           l_y = l_y+ sqrt(2.0)*d
                      id = id+1
           write(22,*) id,  l_x, l_y, l_z, d/2, rho
           end do
         end do
       end do

!Breitenangabe:
       breite =l_z+d/2
       write(*,*) 'Maximum in z:'
       write(*,*) breite

       case (2) auswahl

!Bestimmung der Partikelanzahl durch Länge in X-Richtung
!Rundung der Y-Werte, da Erstellung der Partikel in Z-Richtung in 2 Schritten erfolgt (eine Schleife
!hat aufgerundete und die andere abgerundete Werte) und damit auch bei ungeraden Partikelanzahlen, festgelegte Partikelanzahl erreicht wird
!Bestimmung der Partikelanzahl durch Länge in Z-Richtung und Rundung da Partikelerstellung in 2 Schritten erfolgt
	numrow= floor(zk/(d*sqrt(6.0)/3)) ! Number of rows in x direction
	total_num=floor(ix/d)*nint(numrow/2.)+ (floor(ix-d/2.)/d)*(numrow-nint(numrow/2.))
	! write(22,*)   total_num


       id = 0
       nx = nint(ix/d)
       l_x = -0.5*d


       nz = nint(zk/d)


!       l_y = 0.5*d -d*0.5*sqrt(6.0)/3
      l_y = 0.5*d -d*sqrt(6.0)/3
      do y = 1, j
         l_x = -0.5*d
!          l_y = l_y + sqrt(3.0)*d
!            l_y= -d*sqrt(6.0)/3
           l_y = l_y+ d*sqrt(6.0)/3
          if (modulo(y, 2).eq.0) l_x = l_x - 0.5*d

         do i = 1, nx
          l_x = l_x + d

          l_z = 0.5*d - sqrt(3.) *d*0.5
          if (modulo(y, 2).eq.0) l_z = l_z - 0.5*d*sqrt(3.)/3         
          z = 0
          do while (l_z + sqrt(3.)*d*0.5.le. zk) 
             z = z+1
             if (modulo(z, 2).eq.0) l_x = l_x + 0.5*d
             l_z = l_z + sqrt(3.)*d*0.5
             id = id + 1
             write(22,100)   l_z+zk_offset, l_y+j_offset, l_x+ix_offset, d/2.
100	     format (1x, E14.7,1x, E14.7,1x, E14.7,1x, E14.7,1x)
             if (modulo(z, 2).eq.0) l_x = l_x - 0.5*d
          end do
          if (modulo(y, 2).eq.0) l_z = l_z + 0.5*d*sqrt(3.)/3         
       end do
        if (modulo(y, 2).eq.0) l_x = l_x + 0.5*d
     end do
     
!    
!    
!        gx = ix/(d)
!        i = int(gx)
!    
!        gz =((zk+d)*2)/(sqrt(3.0)*d)
!        k = int(gz)
!    
!        dz = real(k)
!        dz = dz /2
! 
!        dy = real(j)
!        dy = dy/2
! 
!        p = nint(dz)
!        a = nint(dy)
! 
!        q = int(dz)
!        b = int(dy)
! 
! !Erstellen der Partikel in 4 Teilschritten
!        id = 0
! !        l_x=-d
!        l_x=-0.5*d
!        do n = 1, i
!        l_x = l_x+d
! !          l_z=-sqrt(3.0)*d
!          l_z=d*0.5-sqrt(3.0)*d
!          do m = 1, p
!          l_z = l_z+ sqrt(3.0)*d
!            l_y=-d*2*sqrt(6.0)/3
!            do l = 1, a
!            l_y = l_y+ d*2*sqrt(6.0)/3
!                       id = id+1
!            write(22,*) id,  l_x, l_y, l_z, d/2, rho
!            end do
!          end do
!        end do
! 
!        l_x=0.
! !        l_x=-d/2
! 
!        do n = 1, i
!        l_x = l_x+d
! !          l_z= -sqrt(3.0)*d/2
!          l_z = 0.5*d-sqrt(3.0)*d*0.5
!          do m = 1, q
!          l_z = l_z + sqrt(3.0)*d*0.5
!            l_y=-d*2*sqrt(6.0)/3
!            do l = 1, a
!            l_y = l_y+ d*2*sqrt(6.0)/3
!                       id = id+1
!            write(22,*) id,  l_x, l_y, l_z, d/2, rho
!            end do
!          end do
!        end do
! 
! !        l_x=-d/2
!        l_x= 0.
!        do n = 1, i
!        l_x = l_x+d
!          l_z= -5*sqrt(3.0)*d/6
!          do m = 1, p
!          l_z = l_z + sqrt(3.0)*d
!            l_y= -d*sqrt(6.0)/3
!            do l = 1, b
!            l_y = l_y+ d*2*sqrt(6.0)/3
!                       id = id+1
!            write(22,*) id,  l_x, l_y, l_z, d/2, rho 
!            end do
!          end do
!        end do
! 
! !        l_x=-d
!        l_x=-d/2
!        do n = 1, i
!        l_x = l_x+d
!          l_z=-1*sqrt(3.0)*d/3
!          do m = 1, q
!          l_z = l_z+ sqrt(3.0)*d
!            l_y= -d*sqrt(6.0)/3
!            do l = 1, b
!            l_y = l_y+ d*2*sqrt(6.0)/3
!                       id = id+1
!            write(22,*) id,  l_x, l_y, l_z, d/2, rho
!            end do
!          end do
!        end do

!Breitenangabe:
       breite =l_z+d/2
       write(*,*) 'Maximum in z:'
       write(*,*) breite-ix_offset 


      end select auswahl


      end program particle_arrangement
