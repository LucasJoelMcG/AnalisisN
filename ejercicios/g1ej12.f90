program g1ej12
    implicit none
    real, allocatable :: mat(:,:)
    integer n
    print *,'Ingrese N: '
    read (*,*) n
    allocate (mat(n,n))
    call armaMatriz(mat,n)
    print *,'La matriz es diagonal dominante? ',dominante(mat,n)
    deallocate(mat)

    contains

    subroutine armaMatriz(matriz,m)
        integer m,i,j
        real matriz(m,m)
        do i=1,m
            do j=1,m
                print *,'Ingrese valor para ',i,' ',j
                read (*,*) matriz(i,j)
            end do
        end do
    end subroutine

    function dominante(matriz,m)
        integer m,i,j
        real matriz(m,m), aux
        logical dominante
        dominante=.true.
        i=1
        do while (dominante .eqv.(.true.) .and. i<=m)
            aux=0
            do j=1,m
                if (i/=j ) then
                    aux=aux+matriz(i,j)
                end if
            end do
            if (aux>matriz(i,i)) then
                dominante=.false.
            end if
            i=i+1
        end do
    end function

end program g1ej12