Program PuntoFijo

!CAMBIAR LA FUNCIÓN EN EL SCRIPT

 implicit none
 
 real(8) xi, ei, toli!, a, b, itermax
 !INTEGER itermax
 
 toli = 0.0005
 !itermax=100 
 
 call system ("gnuplot -persist 'Script.p'")

 call punto_fijo(xi, toli, ei)!, itermax)
    
  write(*,'(A, F15.5)') 'La raíz en el intervalo es:', xi
  write (*,*) 'El error es ', ei
  write (*,*)
 
  
 CONTAINS
 
  function f (x)
  real(8) x, f
  f = (x - 1.)*sin(3.*x) - (x + 1.)
 end function
 
 function g (x)
  real(8) x, g
  g = (x - 1.)*sin(3.*x) - 1.
 end function
 
 subroutine punto_fijo (x, tol, e)!, itermax)
  real(8) x, x0, tol, e
  integer iter!, itermax
  
  write(*,*) 'Ingrese un valor x0'
  read *, x0
  x = g(x0)
  
  iter = 0 
  e=2.*tol
  
  do while ((abs(f(x))>tol).AND.(abs(x-x0)>tol)) !.and. (iter<itermax)
   x0 = x
   x = g(x)
   
   !write(*, '(I5,2F15.10)'), iter, x, f(x)
   
   e = abs(f(x))
   iter = iter + 1
   !if (iter==7) then
    !write(*,*) 'A la séptima iteración el error es', e
   !end if
  end do
  write(*,*)
  write(*,*) 'Se necesitaron',iter,'iteraciones para encontrar la raíz'
 end subroutine punto_fijo
 
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
 
 subroutine pf_sist(x,tol,e)
  real(8) x, x0, tol, e, Fp, Ex, Ey
  integer iter
  
  write (*,*)
  write (*,*) 'Ingrese un x0 aproximado para calcular la máxima derivada'
  read *, x0
  Fp = fp_max(x0)
  x = x0 - f(x0)/Fp
  
  iter = 0
  do while ((abs(f(x))>tol).AND.(abs(x-x0)>tol))
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
