program ej5
  integer n,i,total
  read (*,*) n
  if (n<=1) then
    print *,'1'
  else
    total=1
    do i=1,n,1
      total=total*i
    end do
    print *,total
  end if
end program ej5