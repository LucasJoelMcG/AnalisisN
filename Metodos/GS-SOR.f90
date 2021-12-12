PROGRAM GSSOR 
IMPLICIT NONE !La entrada de datos sigue el formato mostrado en a1-example.txt
real,allocatable:: ai(:,:),bi(:)
real,allocatable:: x0i(:),xi(:)

integer max_iter,ni,i,j
real wi              !Required for SOR Method
real toli

open(unit=1,file='a1.txt')
open(unit=2,file='a2.txt')

!Read in the number of inputs
read(1,*) ni

allocate(ai(ni,ni),bi(ni),x0i(ni),xi(ni))

!Read in the maximum iterations and tolerance
read(1,*)max_iter,toli
!Read in the coefficient matrix
read(1,*)((ai(i,j),j=1,ni),i=1,ni)
!Read in the constant matrix
read(1,*)(bi(i),i=1,ni)
!Read in the initial approximation
read(1,*)(x0i(i),i=1,ni)

read(1,*) wi

call gauss_seidel(ai,bi,ni,x0i,toli)
call sor(ai,bi,ni,x0i,toli,wi)

contains
!-------------------------------------------------------------------------------
subroutine gauss_seidel(a,b,n,x0,tol)
    integer n,k
    real::a(n,n),b(n),x0(n)
    real::x(n)
    real tol,norm,sum_x
    
    write(2,20)
20  format(/,"Gauss Seidel Result")
    k=0
    do
        k=k+1
        do i=1,n
            sum_x=0
            do j=1,n
                if(j<i)sum_x=sum_x+a(i,j)*x(j)
                if(j>i)sum_x=sum_x+a(i,j)*x0(j)
            end do
            x(i)=(b(i)-sum_x)/a(i,i)
        end do
        write(2,30)k,(x(i),i=1,n)
30      format(2x,i3,4(2x,f9.6))
        norm=abs(x(1)-x0(1))/abs(x(1))
        do i=2,n
            if(abs(x(i)-x0(i))/abs(x(i))>norm)norm=abs(x(i)-x0(i))/abs(x(i))
        end do
        if(norm<tol)goto 15
        do i=1,n
            x0(i)=x(i)
        end do
        
    end do
15  end subroutine
!-------------------------------------------------------------------------------
subroutine SOR(a,b,n,x0,tol,w)
    integer n,k
    real tol,w,norm,sum_x
    real::a(n,n),b(n),x0(n)
    real::x(n)
    
    write(2,40)
40  format(/,"Successive Over Relaxation (SOR) Result")
    k=0
    do
        k=k+1
        do i=1,n
            sum_x=0
            do j=1,n
                if(j<i)sum_x=sum_x+a(i,j)*x(j)
                if(j>i)sum_x=sum_x+a(i,j)*x0(j)
            end do
            x(i)=(1-w)*x0(i)+w*(b(i)-sum_x)/a(i,i)
        end do
        write(2,50)k,(x(i),i=1,n)
50      format(2x,i3,4(2x,f9.6))
        norm=abs(x(1)-x0(1))/abs(x(1))
        do i=2,n
            if(abs(x(i)-x0(i))/abs(x(i))>norm)norm=abs(x(i)-x0(i))/abs(x(i))
        end do
        if(norm<tol)goto 25
        do i=1,n
            x0(i)=x(i)
        end do
    end do
25  end subroutine    

end program 