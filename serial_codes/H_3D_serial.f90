program Heat_equation_2D
    implicit none
 
    integer:: ny , nx , nz

    print* , "Enter_number_of_nodes_in_z-direction"
    read* , nz
    print*,"nz_=_", nz

    print* , "Enter_number_of_nodes_in_y-direction"
    read* , ny
    print*, "ny_=_", ny
    
    print* , "Enter_number_of_nodes_in_x-direction"
    read* , nx
    print*,"nx_=_", nx
    print*,""


    call solve_using_gauss_seidal(nz,ny,nx)
    
 end program

 subroutine solve_using_gauss_seidal(nz,ny,nx)

    integer :: nz, ny, nx, i, j, k, index
    real :: Lz, Ly, Lx, topTemperature, thermalConductivity, Area, dz, dy, dx, cz, cy, cx, Matrix(nz*ny*nx,nz*ny*nx+1)
    real:: an(nz,ny,nx), as(nz,ny,nx), ae(nz,ny,nx), aw(nz,ny,nx), at(nz,ny,nx), ab(nz,ny,nx), ap(nz,ny,nx)
    real:: Temperature(nz*ny*nx)
    print* , "Enter_width_of_surface: "
    read* , Lz
    print* ,"Lz_=_", Lz

    print* , "Enter_height_of_surface: "
    read* , Ly
    print* ,"Ly_=_", Ly

    print* , "Enter_length_of_surface: "
    read* , Lx
    print*, "Lx_=_", Lx
    print*,""

    dz = Lz / nz
    print* ,"dz_=_",dz
    dy = Ly / ny
    print* ,"dy_=_",dy
    dx = Lx / nx
    print* ,"dx_=_",dx
    print*,""

    print* , "Enter_the_Temperature_at_the_top"
    read* , topTemperature
    print*, "top_T_=_", topTemperature
    
    print* , "Enter_Thermal_coductivity"
    read* , thermalConductivity
    print*, "Thermal_connectivity_=_", thermalConductivity
    
    print* , "Enter_the_Area: "
    read* , Area
    print* ,"Area_=_", Area
    print*,""
    
    cz = thermalConductivity*Area/dz
    cy = thermalConductivity*Area/dy
    cx = thermalConductivity*Area/dx

    do k = 1,nz
        do i=1,ny
            do j=1,nx
                if(i<ny) then
                    as(k,i,j) = cy
                else
                    as(k,i,j) = 0
                end if
                print*, "-----------------------------"
                print*, "as(",k,i,j,") = ",as(k,i,j)
                if(j<nx) then
                    ae(k,i,j) = cx
                else
                    ae(k,i,j) = 0
                end if
                print*, "------------------------------"
                print*, "ae(",k,i,j,") = ",ae(k,i,j)
                if(j>1) then
                    aw(k,i,j) = cx
                else
                    aw(k,i,j) = 0
                end if
                print*, "------------------------------------"
                print*, "aw(",k,i,j,") = ",aw(k,i,j)
                if(k<nz) then
                    at(k,i,j) = cz
                else
                    at(k,i,j) = 0
                end if
                print*, "at(",k,i,j,") = ",at(k,i,j)
                if(k>1) then
                    ab(k,i,j) = cz
                else
                    ab(k,i,j) = 0
                end if
                print*, "--------------------------------"
                print*, "ab(",k,i,j,") = ",ab(k,i,j)
                if(i==1) then
                    an(k,i,j) = 0
                    ap(k,i,j) = as(k,i,j) + ae(k,i,j) + aw(k,i,j) + at(k,i,j) + ab(k,i,j) + 2*cy
                else
                    an(k,i,j) = cy
                    ap(k,i,j) = an(k,i,j) + as(k,i,j) + ae(k,i,j) + at(k,i,j) + ab(k,i,j) + aw(k,i,j)
                end if
                print*, "an(",k,i,j,") = ",an(k,i,j)
                print*, "-------------------------------"
                print*, "ap(",k,i,j,") = ",ap(k,i,j)
            end do
        end do
    end do

    

    ! Matrix of size (nx*ny) x (nx*ny+1), to hold all coefficients of equations
    ! initializing to 0
    do i =1,ny*nx
        do j=1,ny*nx+1
            Matrix(i,j) = 0
        end do
    end do

    do k=1,nz
        do i=1,ny
            do j=1,nx
                index = (k - 1) * ny * nx + (i - 1) * nx + j
                Matrix(index,index) = -ap(k,i,j)
                if(((i-1)>0).and.((i-1)<(ny+1))) then
                    Matrix(index,index-nx) = an(k,i,j)
                end if
                if(((i+1)>0).and.((i+1)<(ny+1))) then
                    Matrix(index,index+nx) = as(k,i,j)
                end if
                if(((j+1)>0).and.((j+1)<(nx+1))) then
                    Matrix(index,index+1) = ae(k,i,j)
                end if
                if(((j-1)>0).and.((j-1)<(nx+1))) then
                    Matrix(index,index-1) = aw(k,i,j)
                end if
                if(((k+1)>0).and.((k+1)<(nz+1))) then
                    Matrix(index,index+(ny*nx)) = at(k,i,j)
                end if
                if(((k-1)>0).and.((k-1)<(nz+1))) then
                    Matrix(index,index-(ny*nx)) = ab(k,i,j)
                end if
                if(i==1)then
                    Matrix(index,nx*ny*nz+1) = -2*cy*topTemperature
                end if
            end do
        end do
    end do

    ! Printing The Matrix
    print*, "The_Matrix:_"
    do i = 1,nx*ny*nz
        print* ,(Matrix(i,j),j=1,nx*ny*nz+1)
    end do

    call gauss_seidel(Matrix,nx*ny*nz,Temperature)

    print*,""
    print*,"-----------------------------------------------------------------------"
    print*,"Final Temperatures:"
    do k =1,nz
        print*,"-----------For k_=_",k
        do i=1,ny
            print*,(Temperature((i-1)*nx+j),j=1,nx);
        end do
        print*,"-------------------------------------------------------------------"
        print*,""
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
        print*, "X : ", iter ," = { ",(" (",X(i),") ,",i=1,n),"} "

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

        print*,"No._of_iterations_=_",iter-1

end subroutine
