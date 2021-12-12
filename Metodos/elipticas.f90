program EcElipticas
implicit none
integer,parameter::nx=5,ny=6
real(8), PARAMETER :: tol=0.0001
real(8) Ti(nx,ny)
Ti = 0        !Semillas
Ti(:,1)  = 0.0 !izq
!Ti(:,ny) = 20.0 !der
!Ti(1,:)  = 100.0 !arriba
Ti(nx,:) = 0.0 !abajo
Ti(1,2)=60; Ti(1,3)=120; Ti(1,4)=180; Ti(1,5)=240; Ti(1,6)=300; Ti(2,6)=225; Ti(3,6)=150; 
Ti(4,6)=75; 

call EscribirMat(Ti)
call elipticas(Ti)
call EscribirMat(Ti)
OPEN(UNIT=2,FILE='datos1.dat')
call Datosgrafico(Ti)
CLOSE(2,STATUS='KEEP')

CALL SYSTEM ('gnuplot -persist "plotMat.plt"')

contains

subroutine elipticas(T)
real(8) T(nx,ny), Tant(nx,ny),e
integer i,j

e = 2*tol
do while (e>tol)
  Tant = T
  do i=2,nx-1
    do j=2,ny-1
      T(i,j) = (T(i-1,j) + T(i+1,j) + T(i,j-1) + T(i,j+1) )/4
    enddo
  enddo
  e = sum(abs(T-Tant))
enddo
end subroutine

subroutine EscribirMat(T)
real(8) T(nx,ny)
integer i,j
do i=1,nx
  do j=1,ny
    write(*,'(f10.4)',advance='no') T(i,j)
  enddo
  write(*,*)
enddo
write(*,*)
end subroutine

subroutine Datosgrafico(T)
real(8) T(nx,ny),x,y
integer i,j
y=0.0
do j=1, ny
 x=0.0
 do i=1,nx
 write(2,'(3F10.4)') x,y, T(i,j)
 x=x+3*1/3
 enddo
  write(2,*)
  y=y+3*1/3
enddo
end subroutine
end program