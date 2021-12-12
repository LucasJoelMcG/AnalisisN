program g1ej9
    implicit none
    integer i
    real x
    open(unit=2,file='g1ej9.dat',status='replace')
    do i=1,10,1
        write (2,'(i3)',advance='no') i
        x=1.14*i
        WRITE (2,'(f8.3)',advance='no') x
        WRITE (2,'(f10.5)',advance='no') sin(x)
        WRITE (2,'(f10.5)',advance='no') cos(x)
        WRITE (2,'(f10.5)',advance='no') tan(x)
        WRITE (2,*)
    end do
    CLOSE(2,status='keep')
end program g1ej9