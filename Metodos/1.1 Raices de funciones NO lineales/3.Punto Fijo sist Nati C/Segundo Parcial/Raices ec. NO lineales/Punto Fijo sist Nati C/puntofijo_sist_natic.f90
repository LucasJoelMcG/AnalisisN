Program PuntoFijosist

!CAMBIAR LA FUNCIÓN EN EL SCRIPT 
! CAMBIAR LA DERIVADA Y LA FUNCION EN EL PROGRAMA

 implicit none
 
 real(8) x, e, tol!, a, b
 !integer maxiter
 
 tol = 0.0005
 !maxiter=100
 
 call system ("gnuplot -persist 'Script.p'")
 call system ("gnuplot -persist 'ScriptSist.p'")
  
   ! write(*,*)'Ingrese la cota inferior del intervalo'
   ! read *, a
   ! write(*,*) 'Ingrese la cota superior del intervalo'
   ! read *, b
   ! write(*,*)
    call pf_sist(x,tol,e)
    
 
  write(*,'(A, F15.5)') 'La raíz en el intervalo es:', x
  write (*,*) 'El error es ', e
  write (*,*)
 
 
 CONTAINS
 
  
 function f (x)
  real(8) x, f
  f = (x - 1.)*sin(3.*x) - (x + 1.)
 end function
 
 
 function fp_max(x)
 !Calcula la máx derivada en el entorno de un punto x0 dado por el usuario
  real(8) x, fp_max, fp, h
  integer i, n
  n = 10
  h = 0.001
  fp_max = sin(3.*x) + 3.*cos(3.*x)*(x - 1.) 
  do i=1, n
   x = x + h*i
   fp = sin(3.*x) + 3.*cos(3.*x)*(x - 1.) 
   if (abs(fp)>abs(fp_max)) then
    fp_max = fp
   end if
  end do
  write (*,'(A, F15.5)') 'La máxima derivada del intervalo se encuentra en x=',x
 end function fp_max
 
 subroutine pf_sist(x,tol,e) !maxiter)
  real(8) x, x0, tol, e, Fp, Ex, Ey
  integer iter !,maxiter
  
  write (*,*)
  write (*,*) 'Ingrese un x0 aproximado para calcular la máxima derivada'
  read *, x0
  Fp = fp_max(x0)
  x = x0 - f(x0)/Fp
  
  iter = 0
  do while ((abs(f(x))>tol).AND.(abs(x-x0)>tol)) !.and. (iter<maxiter) 
   x0 = x
   x = x - f(x)/Fp
   iter = iter + 1
  end do
   Ex = abs(x-x0)
   Ey = abs(f(x))
   if (Ex.le.Ey) then
    e = Ey
    else
    e = Ex
   end if
 end subroutine pf_sist

end program
