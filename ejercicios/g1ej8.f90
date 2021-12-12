program g1ej8
    implicit none
    real mat(20,20),max,temp
    integer i,j
    max=0
    do i=1,20,1
        do j=1,20,1
            mat(i,j)=i+j
        end do
    end do
    do j=1,20,1
        temp=0
        do i=1,20,1
            temp=temp+mat(i,j)
        end do
        if (temp>max) then
            max=temp
        end if
    end do
    print *,'La norma de la matriz es: ',max
end program g1ej8