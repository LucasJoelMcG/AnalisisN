 Program NewtonRaphson


!CAMBIAR LA FUNCIÓN EN EL SCRIPT
!CAMBIAR LA FUNCION Y LA DERIVADA EN CONTAINS

 implicit none
 
 real (8) Fp, xi
 real(8), parameter:: tol=0.5E-5  !ESCRIBIR TOLERANCIA, PRESTAR ATENCION SI ES EN X O EN Y   !not cientifica: 1E-8
 integer iter, op, maxiter
 
 call system ("gnuplot -persist Script.p")
 
 op = 1
 maxiter=20 !DEFINIR MAXIMO DE ITERACIONES
 
 do while (op==1)  !LOOP PARA CALCULAR UNA O MAS RAICES
  write (*,*) 'Ingrese el valor de x0'
  read *, xi

  iter=0

  do while (abs(f(xi))>tol .and. (iter<=maxiter))   !VER SI LA TOLERANCIA ES EN X--> (A-B>TOL) O EN Y (abs(f(x))> tol)
   Fp = f_prima(xi)
   xi = xi - (f(xi)/Fp)
    iter = iter + 1
    write (*,'(A, F15.5)') 'La raíz es ', xi
  end do
  !write (*,'(A, F15.5)') 'La raíz es ', x !VER SI QUIERO QUE ESCRIBA TODOS LOS X QUE CALCULA O SOLO EL FINAL
  write (*,*) 'Se necesitaron ', iter,' iteraciones'
  write (*,*)
  write (*,*) 'Presione 1 si desea continuar buscando raíces'
  read *, op
  write (*,*)
 end do
 
 CONTAINS
 
 function f (x)
  real (8) x, f
  
  f = x**4+x**3-4*x**2-3*x+3
  
 end function f
 
 
 function f_prima (x)
  real (8) x, f_prima
  
  f_prima = 4*x**3+3*x**2-8*x-3
 end function
  
end program
