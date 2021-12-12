program parabolicas
IMPLICIT NONE
integer,parameter:: n=10           !nodos
real(8) U(0:n),Ci,dti,dx,ri,alpha,l,tfinal
real(8), PARAMETER :: pi=3.14159265359
integer j

l = 10.0 ![m] longitud
alpha = 1/(pi**2)  !K
tfinal = 1000.0 !tiempo final
dx = l/(n) !avance en x
dti = 0.01 !avance del tiempo
ri =  0.99*dti/(8.96*390*(dx**2))!dti*alpha/dx**2 !calculo r u.u

write(*,*)'dt=',dti
write(*,*)'dx=',dx
write(*,*)'r=',ri

U(0) = 0.0
U(n) = 100.0
do j=1,n-1
  U(j)=0.99/(8.96*390*((dx*j)**2))
enddo
write(*,*) U
read(*,*)

call implicito(U,tfinal,ri,dti)

contains

SUBROUTINE implicito(T,tfin,r,dt)
REAl (8) C,terr(0:n),T(0:n),tfin,r,dt,tol,error,Tant(0:n),tiempo
integer i

tol=0.001
tiempo=0
do while(tiempo.lt.tfin)
 error=2*tol
 Tant=T !semilla
 do while (error>tol)
  Terr=T
  write(*,'(f8.2)',advance='no') tiempo
    do i=0,n
     C=r*Tant(i-1)+(2.-2.*r)*Tant(i)+r*Tant(i+1)
     T(i)=(1./(2.+2.*r))*(C+r*(T(i-1)+T(i+1)))
     write(*,'(f8.2)',advance='no') T(i)
    Enddo
  write(*,'(f8.2)',advance='no') T(n)

  write(*,*)
 error=maxval(abs(T-Terr))
 Enddo
 tiempo=tiempo+Dt
enddo

END SUBROUTINE
end program
