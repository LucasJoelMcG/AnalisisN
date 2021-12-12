program ej5
    integer i,n,a,b,c
    a=0
    b=1
    c=1
    read (*,*) n
    if (n<=1) then
      write(*,*) '0'
    else
      write(*,*) '0'
        do i=1,n-1,1
            write (*,*) c
            c=a+b
            a=b
            b=c
        end do
    end if
end program ej5