Program ejercicio8a
implicit none
real, parameter :: tol=0.005
integer, parameter :: indep=1 , ni=50
real, dimension(0:indep) :: vi
real hi


!programa principal

vi(0)=0.0
vi(1)=0.0
hi=1.

call calcula(vi,hi,ni)

!seccion subrutinas y funciones

contains

function vp(v)
real, dimension(0:1) :: v,vp

vp(0)=1.
vp(1)=sin(v(0))

end function

function rungekutta4to(v,h)
real h
real, dimension(0:1) :: v,k1,k2,k3,k4,rungekutta4to

k1=h*vp(v)
k2=h*vp(v+k1/2.0)
k3=h*vp(v+k2/2.0)
k4=h*vp(v+k3)
v= v + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0

rungekutta4to=v

end function


function primer_metodo(v,h)
real primer_metodo, h
real, dimension(0:1) :: auxv,auxy,v

auxv=v
auxy=v
auxv= rungekutta4to(auxv,h)
auxy= rungekutta4to(auxy,h/2.) 
auxy= rungekutta4to(auxy,h/2.)

do while ((maxval(abs(auxv-auxy))>tol) .or. (maxval(abs(auxv-auxy))<tol/2.))
  if (maxval(abs(auxv-auxy))>tol) then
    h=h/2.
  else if (maxval(abs(auxv-auxy))<tol/2.) then
    h=1.8*h
  end if
  auxv=v
  auxy=v
  auxv= rungekutta4to(auxv,h)
  auxy= rungekutta4to(auxy,h/2.) 
  auxy= rungekutta4to(auxy,h/2.)
end do

primer_metodo = h

end function

subroutine calcula(v,h,n)
real, dimension (0:1) :: v
real h
integer i,n

write(*,'(i5,3f15.7)') 0,h,v
do i=1, n
  h= primer_metodo(v,h)
  v= rungekutta4to(v,h)
  write(*,'(i5,3f15.7)') i,h,v
end do

end subroutine

end program
