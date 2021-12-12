program Ej3

	IMPLICIT NONE
	
	REAL(8),ALLOCATABLE :: u(:), d(:), l(:),b(:,:)
	INTEGER n
	
	write(*,'(A60)')'Ingrese la cantidad de elementos de la diagonal principal'
	read(*,*)n
	
	ALLOCATE(b(n, 1))
	
	CALL ingresovec(l,n)
	CALL muestravec(l,n)
	
	CALL ingresovec(d,n)
	CALL muestravec(d,n)
	
	CALL ingresovec(u,n)
	CALL muestravec(u,n)
	
	CALL cargamatriz(b,n)
	CALL muestramat(b,n,1)!me podría ahorrar el n poniendo orden=SIZE(A,DIM=1)
	
	CALL Thomas(u,d,l,b)
	
	CALL muestramat(b,n,1)
	
	DEALLOCATE(u, d, l, b)
	
  CONTAINS
 
	SUBROUTINE cargamatriz(A,n)
  
	   real(8),intent(out),allocatable::A(:,:)
	   integer n,i,j
	   
	   allocate(A(n,n))
	   
	    do i=1,n,1
			do j=1,1,1	!acá si b tiene más columnas, poner j=1,m,1 pasar m como parámetro (o n porque son matrices cuadradas)
				write(*,'(A21,I2,A1,I2,A1)')'Ingrese el elemento (',i,',',j,')'
				read(*,*)A(i,j)
			end do
			write(*,*)
	    end do
	   
	end SUBROUTINE cargamatriz
	  
	SUBROUTINE muestramat(A,n,m)	!muestra matriz 
	  
	   real(8),intent(inout),allocatable::A(:,:)
	   integer n,m,i,j
	   
	    do i=1,n,1
			do j=1,m,1
				write(*,'(F15.6)',ADVANCE='NO')A(i,j)
			end do
			write(*,*)
	    end do
	   
	end SUBROUTINE muestramat
	
	SUBROUTINE ingresovec(B,n)
  
	    integer n 								!declaracion de variables de argumentos
	    real(8),intent (out),allocatable :: B(:)
	   
	    integer i									!declaracion de variables auxiliares
		
		allocate(B(n))
		
	    do i=1,n,1
			print'(A20,I2,A15)','Ingrese el elemento',i,'del vector'
			read(*,*)B(i)
	    end do
    
    end SUBROUTINE ingresovec
  
    SUBROUTINE muestravec(B,n)
	
		integer n							!declaracion de variables de argumentos
		real(8),intent (inout),allocatable :: B(:)
   
		integer i								!declaracion de variables auxiliares
   
		do i=1,n,1
			write(*,'(F15.4)',ADVANCE='NO')B(i)
		end do
    
	end SUBROUTINE muestravec
	
	SUBROUTINE Thomas(u, d, l, b)
		!Metodo de Thomas para matrices Tri-Diagonales

		REAL(8), INTENT(INOUT), ALLOCATABLE :: u(:), d(:), l(:),b(:,:)

		INTEGER orden, i

		orden = SIZE(b, DIM=1)
		!cant_vec = SIZE(term_indep, DIM=2)
		
		!Realiza copias del original
		!ALLOCATE(u(orden), d(orden), l(orden))
		!u = u_orig
		!d = d_orig
		!l = l_orig
		
		!b = term_indep 
			
		!Aqui comienza el algoritmo
		DO i=1, orden-1
			u(i) = u(i) / d(i)
			b(i,:) = b(i,:) / d(i)
			d(i) = 1.0
			d(i+1) = d(i+1) - l(i+1)*u(i)
			b(i+1,:) = b(i+1,:) - l(i+1)*b(i,:)
			l(i+1) = 0.0
		END DO
		
		!Obtencion de la solucion por Sustitucion Inversa
		b(orden,:) = b(orden,:)/d(orden)
		DO i=orden-1, 1, -1
			b(i,:) = b(i,:) - u(i)*b(i+1,:) / d(i)
		END DO
		!CALL imprimeMatriz(b)

		!DEALLOCATE(u, d, l, b)
		
	END SUBROUTINE	Thomas
		
end program Ej3
