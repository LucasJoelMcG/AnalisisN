program g1ej15
    use g1ej14
    implicit none
    integer opc,m,k
    real, allocatable :: matriz(:,:), vec(:)
    do while (opc/=0)
        print *,'Ingrese una opción, 0 para salir.'
        print *,'1) Ingresar una matriz de nxn.'
        print *,'2) Imprimir matriz de nxn.'
        print *,'3) Grabar en un archivo matriz de nxn.'
        print *,'4) Ingresar un vector de n componentes.'
        print *,'5) Imprimir por pantalla un vector de n componentes.'
        print *,'6) Grabar en un archivo un vector de n componentes.'
        print *,'7) Norma de un vector.'
        print *,'8) Norma de una matriz.'
        read (*,*) opc
        select case (opc)
        case (0)
            print *,'Adios!'
        case (1)
            print *,'Ingrese N de la matriz: '
            read (*,*) m
            if (allocated(matriz)) then
                deallocate (matriz)
                call ingresaMat(matriz,m)
            else
                call ingresaMat(matriz,m)
            end if
        case (2)
            if (allocated(matriz)) then
                call muestraMatriz(matriz,m)
            else
                print *,'Matriz no creada.'
            end if
        case (3)
            if (allocated(matriz)) then
                call grabaMat(matriz,m)
            else
                print *,'Matriz no creada.'
            end if
        case (4)
            print *,'Ingrese N del vector: '
            read (*,*) k
            if (allocated(vec)) then
                deallocate (vec)
                call ingresaVec(vec,k)
            else
                call ingresaVec(vec,k)
            end if
        case (5)
            if (allocated(vec)) then
                call muestraVec(vec,k)
            else
                print *,'Vector no creado.'
            end if
        case (6)
            if (allocated(vec)) then
                call grabaVec(vec,k)
            else
                print *,'Vector no creado.'
            end if
        case (7)
            if (allocated(vec)) then
                print *,'La norma del vector es: ',normaV(vec,k)
            else
                print *,'Vector no creado.'
            end if
        case (8)
            if (allocated(matriz)) then
                print *,'La norma de la matriz es: ',normaM(matriz,m)
            else
                print *,'Matriz no creada.'
            end if
        case default
            print *,'Fuera de rango.'
        end select
    end do
end program 