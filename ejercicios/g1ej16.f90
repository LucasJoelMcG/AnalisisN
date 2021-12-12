program ej5
    integer i,n,a,b,c
    open(2,file='datos.txt',status='replace')
    a=0
    b=1
    c=1
    read (*,*) n
    if (n<=1) then
      write(*,*) '0'
      write(2,*) '0'
    else
      write(*,*) '0'
      write(2,*) '0,0'
        do i=1,n-1,1
          write (*,*) c
          write (2,*) i,c
          c=a+b
          a=b
          b=c
        end do
    end if
    close(2,status='keep')
    call execute_command_line ('gnuplot -p plot.plt')
end program ej5