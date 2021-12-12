program g1ej11
    implicit none
    real, allocatable :: mat(:,:)
    integer n,a
    print *,'ingrese n: '
    read (*,*) n
    allocate (mat(n,n))
    call hilbert(mat,n)
    do a=1,n
        print *,mat(a,:)
    end do
    DEALLOCATE(mat)

    contains
    subroutine hilbert(matriz,m)
    real matriz(m,m)
    integer i,j,m
    do i=1,m
        do j=1,m
            matriz(i,j)=1/(i+j-1.0)
        end do
    end do
    end subroutine

end program g1ej11