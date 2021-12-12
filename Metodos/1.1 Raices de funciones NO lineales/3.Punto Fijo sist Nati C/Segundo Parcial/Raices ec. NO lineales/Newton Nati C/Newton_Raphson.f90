Program NewtonRaphson


!CAMBIAR LA FUNCIÓN EN EL SCRIPT
!CAMBIAR LA FUNCION Y LA DERIVADA EN CONTAINS

 implicit none
 
 real (8) Fp, x
 real(8), parameter:: tol=1E-8  !ESCRIBIR TOLERANCIA, PRESTAR ATENCION SI ES EN X O EN Y
 integer iter, op, maxiter
 
 call system ("gnuplot -persist 'Script.p'")
 
 op = 1
 maxiter=3 !DEFINIR MAXIMO DE ITERACIONES
 
 do while (op==1)  !LOOP PARA CALCULAR UNA O MAS RAICES
  write (*,*) 'Ingrese el valor de x0'
  read *, x

  iter=0

  do while (abs(f(x))>tol .and. (iter<=maxiter))   !VER SI LA TOLERANCIA ES EN X--> (A-B>TOL) O EN Y (F(a)-F(b)> tol)
   Fp = f_prima(x)
   x = x - f(x)/Fp
    iter = iter + 1
  end do
  write (*,'(A, F15.5)') 'La raíz es ', x
  write (*,*) 'Se necesitaron ', iter,'iteraciones'
  write (*,*)
  write(*,*) 'Presione 1 si desea continuar buscando raíces'
  read *, op
  write(*,*)
 end do
 
 CONTAINS
 
 function f (x)
  real (8) x, f
  
  f = x**4. + x**3. - 4.*x**2. -3.*x +3.
  
 end function f
 
 
 function f_prima (x)
  real (8) x, f_prima
  
  f_prima = 4.*x**3. + 3.*x**2. - 8.*x - 3.
  
 end function
  
end program
