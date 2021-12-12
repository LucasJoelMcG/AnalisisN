program Parabolic
implicit none
integer,parameter:: n=10 !Cant de nodos
real(8) Ui(0:n),dti,dx,ri,alpha,l,tfini
integer j
real(8), PARAMETER :: pi=3.14159265359

l = 1.0 ![m] longitud ponele
alpha = (1.0/(pi**2))  !ponele que es el K
tfini = 2.0 !tiempo final
dx = l/(n)  !pasos en x
dti = 0.04 !pasos en el tiempo
ri =  (dti*alpha)/(dx**2) !tengo que calcularlo :c
write(*,*)'dt=',dti
write(*,*)'dx=',dx
write(*,*)'r=',ri

Ui(0) = 0.
Ui(1) = 0.
do j=2,n
  Ui(j)=cos(pi*(dx*j-0.5))
end do
write(*,*) Ui
read(*,*)
call Explicito(Ui,tfini,ri,dti)

contains
subroutine Explicito(U,tfin,r,dt)
real(8) U(0:n),tfin,r,tf,tant(0:n),dt
integer i

tf=0
do while (tf.lt.tfin)
  tant=U
  tf=tf+dt
  write(*,'(f8.4)',advance='no') tf
  write(*,'(f8.2)',advance='no') U(0)
  do i=0,n
    U(i)=r*(tant(i-1)+tant(i+1))+(1-2*r)*tant(i)
    write(*,'(f8.2)',advance='no') U(i)
  end do
  write(*,'(f8.2)',advance='no') U(n)

  write(*,*)
end do

end subroutine


end program
