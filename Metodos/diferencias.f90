program aproximacion
!Este programa calcula la aproximacion de funciones dados los puntos


implicit none

real(8),dimension(:,:),allocatable :: mati
real(8),dimension(:),allocatable :: vecdi,vecai,vecdii
integer ans
Write (*,*) '*********************Puntos equispaciados*************************************'
Write (*,*) 'Ingrese 1 si quiere diferencias ascendentes (valor a la cabeza de la tabla)'
Write (*,*) 'Ingrese 2 si quiere diferencias descendentes (valor cerca del fin de la tabla)'
Write (*,*) ' '
Write (*,*) '*********************Puntos no equispaciados**********************************'
Write (*,*) 'Ingrese 3 si quiere diferencias divididas'
read*,ans
select case (ans)
  case(1)
    call system ("cls")
    call ingreso_puntos(mati)
    call system ("cls")
    call ascendentes(mati,vecai)
    call muestrapoli(vecai,mati,1)
    call grabapoli(vecai,mati,1)
    Write (*,*) 
    Write (*,*) 'RECORDAR S=(X-Xo)/h'
  case(2)
    call system ("cls")
    call ingreso_puntos(mati)
    call system ("cls")
    call descendentes(mati,vecdi)
    call muestrapoli(vecdi,mati,2)
    call grabapoli(vecdi,mati,2)
    Write (*,*) 
    Write (*,*) 'RECORDAR S=(X-Xo)/h'
  case(3)
    call system ("cls")
    call ingreso_puntos(mati)
    call system ("cls")
    call divididas(mati,vecdii)
    call muestrapoli(vecdii,mati,3)
    call grabapoli(vecdii,mati,3)
end select

!***********************************************************************
contains

subroutine ingreso_puntos (mat)
real(8),dimension(:,:),allocatable:: mat
integer cant,i
!integer j

print*,'Ingrese la cantidad de pares de puntos'
read*,cant
allocate (mat(cant,-1:cant-1))  
! Armo una matriz de la forma   (xi  f(xi)  D  D2  D3  ...   Dn )
mat=0
do i=1,cant
  print*,'Ingrese el punto P',i
  read*,mat(i,-1)
  read*,mat(i,0)
end do

!do i=1,cant
!  do j=-1,cant-1
!   write (*,'(f10.6)',advance='no'),mat(i,j)
!  end do
!  write (*,*)
!end do
!read*,
end subroutine
!***********************************************************************

subroutine ascendentes(mat,veca)
real(8),dimension(:,:),allocatable :: mat
real(8),dimension(:),allocatable :: veca
integer i,j,cant

cant=size(mat,dim=1)
allocate(veca(0:cant-1))
do  j=1,cant-1 
  do i=1,cant-j
   mat(i,j) = mat(i+1,j-1)-mat(i,j-1)
  end do
end do

print*,'El vector de diferencias ascendentes es '

do i=0,cant-1
   veca(i)=mat(1,i)
   write (*,'(f10.6)',advance='no') mat (1,i)
end do

print*,
print*,'La matriz queda '
print*,

do i=1,cant
  do j=-1,cant-1
   write(*,'(f10.6)',advance='no') mat(i,j)
  end do
  write(*,*)
end do
read*,

end subroutine

!***********************************************************************
subroutine descendentes (mat,vecd)
real(8),dimension(:,:),allocatable :: mat
real(8),dimension(:),allocatable :: vecd
integer i,j,cant

cant=size(mat,dim=1)
allocate(vecd(0:cant-1))
do  j=1,cant-1 !calcula la columna 1 de los D1, a la columna cant-1
  do i=cant,2+j-1,-1
   mat(i,j) = mat(i,j-1)-mat(i-1,j-1)
  end do
end do

print*,'El vector de diferencias descendentes es '
do i=0,cant-1
   vecd(i)=mat(cant,i)
   write (*,'(f17.6)',advance='no') mat(cant,i)
end do
print*,
print*,'La matriz queda '
print*,
do i=1,cant
  do j=-1,cant-1
   write(*,'(f10.6)',advance='no') mat(i,j)
 end do
 write(*,*)
end do
read*,

end subroutine
!***********************************************************************

subroutine divididas(mat,vecd)
real(8),dimension(:,:),allocatable :: mat
real(8),dimension(:),allocatable :: vecd
integer i,j,cant

cant=size(mat,dim=1)
allocate(vecd(0:cant-1))
do  j=1,cant-1 
  do i=1,cant-j
   mat(i,j) = (mat(i+1,j-1)-mat(i,j-1))/(mat(i+j,-1)-mat(i,-1))
  end do
end do
print*,'El vector de diferencias divididas es '
do i=0,cant-1
  vecd(i)=mat(1,i)
  write (*,'(f10.6)',advance='no') mat (1,i)
end do
print*,
print*,'La matriz queda '
print*,
do i=1,cant
  do j=-1,cant-1
   write(*,'(f10.6)',advance='no') mat(i,j)
  end do
  write(*,*)
end do
read*,

end subroutine

!***********************************************************************

function fact(n)
integer fact, n,i
fact=1
do i=1,n
  fact=fact*i
end do
end function

!***********************************************************************

subroutine muestrapoli(vector,mat,metodo)
real(8), intent(in), allocatable, dimension(:) :: vector
real(8), intent(in), allocatable, dimension(:,:) :: mat
integer, intent(in) :: metodo
integer cant,i,j

cant=size(mat,dim=1)
write (*,*) ' El polinomio resultante es:'
if (metodo==1) then
  write (*,*) ' metodo ascendente'
  write (*,'(A1,I2,A6)', advance='no') 'P',cant-1,'(S) = '
  do i=0,cant-1
    if (i==0) then
      write (*,'(f10.6,A3)',advance='no')  vector(i),' + '
    else 
      if ((i==cant-1).and.(i/=1)) then
        do j=0,i-1
          if (j/=i-1) then
            write (*,'(A3,I2,A2)',advance='no')  '(S-',j,')*'
          else 
            if (j==i-1) then
              write (*,'(A3,I2,A1)',advance='no')  '(S-',j,')'
            end if
          end if
        end do
        write (*,'(A4,I5,A2,f10.6)',advance='no') '*(1/',fact(i),')*',vector(i)
      else 
        if ((i==1).and.(i/=cant-1)) then
          write (*,'(f10.6,A5)',advance='no')  vector(i),'*S + ' 
        else 
          if ((i==1).and.(i==cant-1)) then
            write (*,'(f10.6,A2)',advance='no')  vector(i),'*S' 
          else
            do j=0,i-1
              if (j/=i-1) then
                write (*,'(A3,I2,A2)',advance='no')  '(S-',j,')*'
              else 
                if (j==i-1) then
                  write (*,'(A3,I2,A1)',advance='no')  '(S-',j,')'
                end if
              end if
            end do
            write (*,'(A4,I5,A2,f10.6,A3)',advance='no') '*(1/',fact(i),')*',vector(i),' + '
          end if
        end if
      end if
    end if
  end do
else 
  if (metodo==2) then
    write (*,*) ' metodo descendente'
    write (*,'(A1,I2,A6)', advance='no') 'P',cant-1,'(S) = '
    do i=0,cant-1
      if (i==0) then
        write (*,'(f10.6,A3)',advance='no')  vector(i),' + '
      else 
        if ((i==cant-1).and.(i/=1)) then
          do j=0,i-1
            if (j/=i-1) then
              write (*,'(A3,I2,A2)',advance='no')  '(S+',j,')*'
            else 
              if (j==i-1) then
                write (*,'(A3,I2,A1)',advance='no')  '(S+',j,')'
              end if
            end if
          end do
          write (*,'(A4,I5,A2,f10.6)',advance='no') '*(1/',fact(i),')*',vector(i)
        else 
          if ((i==1).and.(i/=cant-1)) then
            write (*,'(f10.6,A5)',advance='no')  vector(i),'*S + ' 
          else 
            if ((i==1).and.(i==cant-1)) then
              write (*,'(f10.6,A2)',advance='no')  vector(i),'*S'  
            else
              do j=0,i-1
                if (j/=i-1) then
                  write (*,'(A3,I2,A2)',advance='no')  '(S+',j,')*'
                else 
                  if (j==i-1) then
                    write (*,'(A3,I2,A1)',advance='no')  '(S+',j,')'
                  end if
                end if
              end do
              write (*,'(A4,I5,A2,f10.6,A3)',advance='no') '*(1/',fact(i),')*',vector(i),' + '
            end if
          end if
        end if
      end if
    end do
  else 
    if (metodo==3) then
      write (*,*) ' metodo divididas'
      write (*,'(A1,I2,A6)', advance='no') 'P',cant-1,'(x) = '
      do i=1,cant
        if (i==1) then
          write (*,'(f10.6,A3)',advance='no')  vector(i-1),' + '
        else 
          if (i==cant) then
            do j=1,i-1
              if (j/=i-1) then
                write (*,'(A3,f10.6,A2)',advance='no')  '(X-',mat(j,-1),')*'
              else 
                if (j==i-1) then
                  write (*,'(A3,f10.6,A1)',advance='no')  '(X-',mat(j,-1),')'
                end if
              end if
            end do
            write (*,'(A1,f10.6)',advance='no') '*',vector(i-1)
          else
            do j=1,i-1
              if (j/=i-1) then
                write (*,'(A3,f10.6,A2)',advance='no')  '(X-',mat(j,-1),')*'
              else 
                if (j==i-1) then
                  write (*,'(A3,f10.6,A1)',advance='no') '(X-',mat(j,-1),')'
                end if
              end if
            end do
            write (*,'(A1,f10.6,A3)',advance='no') '*',vector(i-1),' + '
          end if
        end if
      end do
    end if 
  end if
end if
end subroutine

!***********************************************************************

subroutine grabapoli(vector,mat,metodo)
real(8), intent(in), allocatable, dimension(:) :: vector
real(8), intent(in), allocatable, dimension(:,:) :: mat
integer, intent(in) :: metodo
integer cant,i,j

cant=size(mat,dim=1)
open (3,file='poli.dat',status='replace')
write (3,*) ' El polinomio resultante es:'
if (metodo==1) then
  write (3,*) ' metodo ascendente'
  write (3,'(A1,I2,A6)', advance='no') 'P',cant-1,'(S) = '
  do i=0,cant-1
    if (i==0) then
      write (3,'(f10.6,A3)',advance='no')  vector(i),' + '
    else 
      if ((i==cant-1).and.(i/=1)) then
        do j=0,i-1
          if (j/=i-1) then
            write (3,'(A3,I2,A2)',advance='no')  '(S-',j,')*'
          else 
            if (j==i-1) then
              write (3,'(A3,I2,A1)',advance='no')  '(S-',j,')'
            end if
          end if
        end do
        write (3,'(A4,I5,A2,f10.6)',advance='no') '*(1/',fact(i),')*',vector(i)
      else 
        if ((i==1).and.(i/=cant-1)) then
          write (3,'(f10.6,A5)',advance='no')  vector(i),'*S + ' 
        else 
          if ((i==1).and.(i==cant-1)) then
            write (3,'(f10.6,A2)',advance='no')  vector(i),'*S' 
          else
            do j=0,i-1
              if (j/=i-1) then
                write (3,'(A3,I2,A2)',advance='no')  '(S-',j,')*'
              else 
                if (j==i-1) then
                  write (3,'(A3,I2,A1)',advance='no')  '(S-',j,')'
                end if
              end if
            end do
            write (3,'(A4,I5,A2,f10.6,A3)',advance='no') '*(1/',fact(i),')*',vector(i),' + '
          end if
        end if
      end if
    end if
  end do
else 
  if (metodo==2) then
    write (3,*) ' metodo descendente'
    write (3,'(A1,I2,A6)', advance='no') 'P',cant-1,'(S) = '
    do i=0,cant-1
      if (i==0) then
        write (3,'(f10.6,A3)',advance='no')  vector(i),' + '
      else 
        if ((i==cant-1).and.(i/=1)) then
          do j=0,i-1
            if (j/=i-1) then
              write (3,'(A3,I2,A2)',advance='no')  '(S+',j,')*'
            else 
              if (j==i-1) then
               write (3,'(A3,I2,A1)',advance='no')  '(S+',j,')'
              end if
            end if
          end do
          write (3,'(A4,I5,A2,f10.6)',advance='no') '*(1/',fact(i),')*',vector(i)
        else 
          if ((i==1).and.(i/=cant-1)) then
            write (3,'(f10.6,A5)',advance='no')  vector(i),'*S + ' 
          else
            if ((i==1).and.(i==cant-1)) then
              write (3,'(f10.6,A2)',advance='no')  vector(i),'*S' 
            else
              do j=0,i-1
                if (j/=i-1) then
                  write (3,'(A3,I2,A2)',advance='no')  '(S+',j,')*'
                else
                  if (j==i-1) then
                    write (3,'(A3,I2,A1)',advance='no')  '(S+',j,')'
                  end if
                end if
              end do
              write (3,'(A4,I5,A2,f10.6,A3)',advance='no') '*(1/',fact(i),')*',vector(i),' + '
            end if
          end if
        end if
      end if
    end do
  else 
    if (metodo==3) then
      write (3,*) ' metodo divididas'
      write (3,'(A1,I2,A6)', advance='no') 'P',cant-1,'(x) = '
      do i=1,cant
        if (i==1) then
          write (3,'(f10.6,A3)',advance='no')  vector(i-1),' + '
        else 
          if (i==cant) then
            do j=1,i-1
              if (j/=i-1) then
                write (3,'(A3,f10.6,A2)',advance='no')  '(X-',mat(j,-1),')*'
              else
                if (j==i-1) then
                  write (3,'(A3,f10.6,A1)',advance='no')  '(X-',mat(j,-1),')'
                end if
              end if
            end do
          end if
        end if  
      write (3,'(A1,f10.6)',advance='no') '*',vector(i-1)
      end do
    else
      do j=1,i-1
        if (j/=i-1) then
          write (3,'(A3,f10.6,A2)',advance='no')  '(X-',mat(j,-1),')*'
        else 
          if (j==i-1) then
            write (3,'(A3,f10.6,A1)',advance='no') '(X-',mat(j,-1),')'
          end if
        end if
      end do
    write (3,'(A1,f10.6,A3)',advance='no') '*',vector(i-1),' + '
    end if
  end if
end if 
close (3, status='keep')
end subroutine
!***********************************************************************
end program