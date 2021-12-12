program g1ej13
    implicit none
    real vec(20),mat(3,3)
    integer i,j
    do i=1,20,1
        vec(i)=i**2
    end do
    vec(10)=500
    print *,'La norma del vector es: ',maxval(vec)
    do i=1,3,1
        do j=1,3,1
            mat(i,j)=i+j
        end do
    end do
    print *,'La norma de la matriz es: ',maxval(mat)
    print *,'La suma de sus valores es: ',sum(mat)
end program g1ej13