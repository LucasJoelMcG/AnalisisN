Program ejercicio10
    implicit none
    real, allocatable :: mat(:,:)
    integer n,i,j
    print *,'ingrese n: '
    read (*,*) n
    OPEN (unit=2,file='g1ej10.dat',status='replace')
    allocate (mat(n,n))
    do i=1,n
        do j=1,n
            mat(i,j)=1/(i+j-1.0)
        end do
        WRITE (2,*) mat(i,:)
    end do
    CLOSE (2,status='keep')
    DEALLOCATE(mat)
End program