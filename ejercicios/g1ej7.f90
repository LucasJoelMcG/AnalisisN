program g1ej7
    implicit none
    real vec(20)
    integer i,n,max
    n=20
    max=0
    do i=1,20,1
        vec(i)=i**2
        print *,i,'   ',vec(i)
    end do
    vec(10)=500
    do i=1,20,1
        if (abs(vec(i))>max) then
            max=abs(vec(i))
        end if
    end do
    print *,'La norma del vector es: ',max
 end program g1ej7