program Heat_equation_2D
    implicit none
 
    integer:: ny , nx

    print* , "Enter number of nodes in y-direction"
    read* , ny
    print*, "ny = ", ny
    
    print* , "Enter number of nodes in x-direction"
    read* , nx
    print*,"nx = ", nx
    print*,""
    
    call solve_using_gauss_seidal(ny,nx)
    
 end program
 
 
 subroutine solve_using_gauss_seidal(ny,nx)
    
    integer :: ny, nx, i, j, index
    real :: Ly, Lx,topTemperature, thermalConductivity, Area , an(ny,nx), as(ny,nx), ae(ny,nx), aw(ny,nx), ap(ny,nx)
    real :: dx, dy, cx, cy, Matrix(ny*nx,ny*nx+1),Temperature(ny*nx)
    
    print* , "Enter height of surface: "
    read* , Ly
    print* ,"Ly = ", Ly
    
    print* , "Enter length of surface: "
    read* , Lx
    print*, "Lx = ", Lx
    print*,""
    
    dy = Ly / ny
    print* ,"dy = ",dy
    dx = Lx / nx
    print* ,"dx = ",dx
    print*,""
    
    print* , "Enter the Temperature at the top"
    read* , topTemperature
    print*, "top T = ", topTemperature
    
    print* , "Enter Thermal coductivity"
    read* , thermalConductivity
    print*, "Thermal connectivity = ", thermalConductivity
    
    print* , "Enter the Area: "
    read* , Area
    print* ,"Area = ", Area
    print*,""
    
    cx = thermalConductivity*Area/dx
    cy = thermalConductivity*Area/dy
    
    do i=1,ny
        do j=1,nx
            if(i<ny) then
                as(i,j) = cy
            else
                as(i,j) = 0
            end if
            if(j<nx) then
                ae(i,j) = cx
            else
                ae(i,j) = 0
            end if
            if(j>1) then
                aw(i,j) = cx
            else
                aw(i,j) = 0
            end if
            if(i==1) then
                an(i,j) = 0
                ap(i,j) = as(i,j) + ae(i,j) + aw(i,j) + 2*cy
            else
                an(i,j) = cy
                ap(i,j) = an(i,j) + as(i,j) + ae(i,j) + aw(i,j)
            end if
        end do
    end do
    
    ! Matrix of size (nx*ny) x (nx*ny+1), to hold all coefficients of equations
    ! initializing to 0
    do i =1,ny*nx
        do j=1,ny*nx+1
            Matrix(i,j) = 0
        end do
    end do

    ! putting coefficients in the matrix of equations
    do i=1,ny
        do j=1,nx
            index = (i - 1) * nx + j
            Matrix(index,index) = -ap(i,j)
            if(((i-1)>0).and.((i-1)<(ny+1))) then
                Matrix(index,index-nx) = an(i,j)
            end if
            if(((i+1)>0).and.((i+1)<(ny+1))) then
                Matrix(index,index+nx) = as(i,j)
            end if
            if(((j+1)>0).and.((j+1)<(nx+1))) then
                Matrix(index,index+1) = ae(i,j)
            end if
            if(((j-1)>0).and.((j-1)<(nx+1))) then
                Matrix(index,index-1) = aw(i,j)
            end if
            if(i==1)then
                Matrix(index,nx*ny+1) = -2*cy*topTemperature
            end if
        end do
    end do
    
    ! Printing The Matrix
    print*, "The Matrix: "
    do i = 1,nx*ny
        print* ,(Matrix(i,j),j=1,nx*ny+1)
    end do

    ! calling subroutine to solve matrix using gauss-seidal method
    call gauss_seidel(Matrix,nx*ny,Temperature)

    print*,""
    print*,"------------------------------------------------------------------------"
    print*,"Final Temperatures:"
    do i=1,ny
        print*,(Temperature((i-1)*nx+j),j=1,nx);
    end do
    print*,"-----------------------------------------------------------------------"

end subroutine

subroutine gauss_seidel(M,n,X)

    integer :: i,j,n, iter, stp, max_iter , count
    real :: M(n,n+1), eps, X(n), P(n), sum
    
    iter = 0
    eps = 1e-15

    count = 2000
    do i = 1,n
        X(n) = 0
    end do

    do while(count>0)
        do i = 1,n
            sum = M(i,n+1)
            do j= 1,n
                if (j .ne. i) then
                    sum = sum - M(i,j)* X(j)
                end if
            end do
            X(i) = 1/M(i,i)*sum;
        end do
        print*, "X :", iter ,"= {",("(",X(i),"),",i=1,n),"}"

        iter = iter + 1

        if (iter == 1) then
            continue
        end if

        stp = 1
        do i = 1,n
            if(stp==1) then
                if (ABS(X(i) - P(i))> eps) then 
                    stp = 0
                end if
            end if
        end do

        if ((stp==1).or.(iter==max_iter)) then
            exit
        end if

        do i = 1, n
            P(i)=X(i)
        end do
        
            count = count - 1
        end do

        print*,"No. of iterations =",iter-1

end subroutine