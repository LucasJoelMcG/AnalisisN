Program ejercicio1
implicit none
integer, parameter :: n=4 , indep=1    !n:orden de la matriz, indep= orden de vector terminos independientes
real(8) Mi(n,n), bi(n), mat_ampl(n,n+indep)

Mi(1,1)=2.1756;    Mi(1,2)=4.0231;     Mi(1,3)=-2.1732;     Mi(1,4)=5.1967;     bi(1)=17.102
Mi(2,1)=-4.0231;    Mi(2,2)=6.0;     Mi(2,3)=0.0;     Mi(2,4)=1.1973;     bi(2)=-6.1593
Mi(3,1)=-1.0;    Mi(3,2)=-5.2107;    Mi(3,3)=1.1111;     Mi(3,4)=0.0;     bi(3)=3.0004
Mi(4,1)=6.0235;    Mi(4,2)=7.0;    Mi(4,3)=0.0;     Mi(4,4)=-4.1561;     bi(4)=0.0

call matrizampliada(Mi,bi,mat_ampl)
call gauss(mat_ampl)
call regresiva(mat_ampl) 

contains

subroutine matrizampliada(M,b,mat_amp)
real(8) M(n,n), b(n)
real(8) mat_amp(n,n+indep)
integer i,j

mat_amp=0.0

print*, 'la matriz ampliada es:'
do i=1, n
  do j=1, n
    mat_amp(i,j)=M(i,j)
  end do
end do
do i=1,n
  mat_amp(i,n+1)=b(i)
end do
do i=1, n
  write(*, '(5f10.3)', advance='no') mat_amp(i,:)
  write(*,*)
end do

end subroutine

subroutine pivoteo(M,col)
real(8) M(n,n+indep)
integer i, col,fila
real(8) maximo,aux

!analiza tambien el primer valor
write(*,*)
maximo=0.0
do i=col, n
  if (abs(M(col,i))>maximo) then
    maximo=abs(M(col,i))
    fila=i 
  end if
end do

do i=col, n+1
  aux=M(col,i)
  M(col,i)=M(fila,i)
  M(fila,i)=aux 
end do
write(*, '(5f10.3)') M 
write(*,*)

end subroutine

subroutine gauss(M)
real(8) M(n,n+1)
integer t,fila

Do t=1, n !t=columna
call pivoteo(M,t)
  do fila=t+1, n
    M(fila, t+1:) = M(fila, t+1:) - (M(t, t+1:)*M(fila, t)/M(t, t))
    M(fila,t)=0.0
  end do
end do

end subroutine

Subroutine regresiva(M)
real(8) M(n,n+1),x(n)
integer i,j,aux
real(8) suma

aux=n
x(n)=M(n,n+1)/M(n,n)
do i=aux-1,1,-1
  suma=0.0
  do j=i+1, n
    suma=suma + (M(i,j)*x(j))
  end do
  x(i)=(M(i,n+1) - suma)/M(i,i)
end do
print *, 'El vector solución es: '
write (*,'(4f15.5)') X !cambiar segun el orden de la matriz

end subroutine

End program
