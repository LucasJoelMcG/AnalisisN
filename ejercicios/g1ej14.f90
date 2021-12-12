module g1ej14
    implicit none
    !! integer opc,m,k
    !! real, allocatable :: matriz(:,:), vec(:)

    contains

    subroutine ingresaMat(mat,n)
    integer n,i,j
    real, allocatable :: mat(:,:)
    allocate (mat(n,n))
    do i=1,n,1
        do j=1,n,1
            print *,'Ingrese valor (',i,', ',j,')'
            read (*,*) mat(i,j)
        end do
    end do
    end subroutine

    subroutine muestraMatriz(mat,n)
    integer n,i
    real mat(n,n)
    do i=1,n,1
        print *,mat(i,:)
    end do
    end subroutine
    
    subroutine grabaMat(mat,n)
    integer n,i
    real mat(n,n)
    open (2,file='Matriz.dat',status='replace')
    do i=1,n,1
        write (2,*) mat(i,:)
    end do
    close (2,status='keep')
    end subroutine

    subroutine ingresaVec(vector,n)
    integer n,i
    real, allocatable :: vector(:)
    allocate (vector(n))
    do i=1,n,1
        print *,'Ingrese valor ',i,' del vector: '
        read (*,*) vector(i)
    end do
    end subroutine

    subroutine muestraVec(vector,n)
    integer n,i
    real vector(n)
    print *,vector
    end subroutine

    subroutine grabaVec(vector,n)
    integer n,i
    real vector(n)
    open (2,file='Vector.dat',status='replace')
    write (2,*) vector
    close (2,status='keep')
    end subroutine

    function normaV(vector,n)
    integer n,i
    real vector(n),normaV
    normaV=maxval(vector)
    end function

    function normaM(mat,n)
    integer n,i,j
    real mat(n,n), normaM
    normaM=maxval(mat)
    end function

end module g1ej14