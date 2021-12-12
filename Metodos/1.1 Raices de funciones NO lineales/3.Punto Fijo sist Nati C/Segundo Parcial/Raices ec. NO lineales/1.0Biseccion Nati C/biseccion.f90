Program biseccion

!CAMBIAR LA FUNCIÓN EN EL SCRIPT

 implicit none
 
 real(8) a, b, tol, e, m
 integer i, n
 
 call system ("gnuplot -persist 'Script.p'")
 
 write (*,*) 'Ingrese el valor de a'
 read *, a
 write (*,*)
 write (*,*) 'Ingrese el valor de b'
 read *, b
 
 m = (b + a)/2. !calculo del punto medio del subintervalo
 
 write (*,'(A, F15.5)') 'El valor inicial de m es:', m
 
 !============ESCRIBIR LA TOLERANCIA, TENER EN CUENTA SI ES EN Y O EN X===========
 tol = 0.001
 
 n = int ((log ((b-a)/tol))/log (2.) + 0.5)
 
 e = (b-a)/2.**n
 
 write (*,'(A, I4)') 'El número de iteraciones será:', n
 
 do i=1, n     !=ver si necesito do while, (A-B)>tolx ó |f(a)-f(b)|>toly, .and. (i<=n)
  if (f(a)*f(m)>0) then
   a = m
   else
   b = m
  end if
  m = (b + a)/2.
  e = (b-a)/2.**i
 end do
 
 write (*,'(2F15.5)') m, e
 
 
 CONTAINS
 
  function f (x)
   real(8) x, f
  
   f = 0.5*exp(x/3)-sin(x)
  end function

end program biseccion
