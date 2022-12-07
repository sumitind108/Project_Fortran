
program Heat_equation_1D
   implicit none

   integer:: n, i, j
   real :: TA, TB, thermalConductivity, Area, q
   real, dimension (:), allocatable :: MatA, MatB, Temperature, a, b, c, d

   print*, "Enter values for - TA , TB , n , thermalConductivity , Area: "
   read* , TA , TB , n , thermalConductivity , Area

   allocate (MatA(n*n))
   allocate (MatB(n))
   allocate (Temperature(n))
   allocate (a(n))
   allocate (b(n))
   allocate (c(n))
   allocate (d(n))
   
   ! Step to make Equations(coefficients)
   do i =1,n*n
      MatA(i) = 0
   end do

   ! 2D Array: MatA(i,j) corresponds to 1D Array: MatA((i-1)*n+j)
   MatA(1) = 3          ! MatA(1,1) = 3
   MatA(n*n) = 3        ! MatA(n,n) = 3
   MatA(2) = -1         ! MatA(1,2) = -1
   MatA(n+1) = -1       ! MatA(2,1) = -1

   do i=2,n-1
      MatA((i-1)*n+i) = 2     ! MatA(i,i) = 2
      MatA((i-1)*n+i+1) = -1  ! MatA(i,i+1) = -1
      MatA(i*n+i) = -1        ! MatA(i+1,i) = -1
   end do

   ! ! printing MatA
   ! print* , "MatA:"
   ! do i=1,n
   !    print* ,(MatA((i-1)*n+j),j=1,n)
   ! end do

   MatB(1)=2*TA   ! MatB(1,1) = 2*TA
   MatB(n)=2*TB   ! MatB(n,1) = 2*TB
   do i=2,n-1
      MatB(i) = 0    ! all values from MatB(2,1) to MatB(n-1,1) initialize 0
   end do

   ! ! Printing MatB
   ! print* ,"MatB:"
   ! do i = 1, n
   !    print* ,MatB(i)
   ! end do

   ! creating three diagonals a,b,c and d is right matrix
   ! 2D Array: MatA(i,j) corresponds to 1D Array: MatA((i-1)*n+j)
   do i=1,n
      if(i==1) then
         a(1)=0
      else
         a(i) = MatA((i-1)*n+(i-1))
      end if
      j = j + 1
      b(i)=MatA((i-1)*n+i)
      if(i==n) then
         c(n)=0
      else
         j = j + 1
         c(i) = MatA(((i+1)-1)*n+i)
      end if
      d(i) =MatB(i)
   end do


   !  --- Elimination ---
   do i = 2,n
      q = a(i)/b(i - 1)
      b(i) = b(i) - c(i - 1)*q
      d(i) = d(i) - d(i - 1)*q
   end do

   ! --- Backsubstitution ---
   q = d(n)/b(n)
   Temperature(n) = q
   do i = n - 1,1,-1
      q = (d(i) - c(i)*q)/b(i)
      Temperature(i) = q
   end do

   print*,"------------------------------"
   print* ,"Temperature:"
   do i=1,n
      print* ,Temperature(i)
   end do   ! print*,""
   ! print*,":::::::::::::::::::::::::::::::::::::::::::::"
   ! do i=1,ny*nx*nz
   !     sum = 0
   !     do j=1,ny*nx*nz
   !         sum = sum + Matrix(i,j)
   !     end do
   !     print*,"for_i_=_",i,"___Ans_=_",sum
   ! end do
   ! print*,":::::::::::::::::::::::::::::::::::::::::::::"
   ! print*,""
   print*,"------------------------------"

   deallocate (MatA)
   deallocate (MatB)
   deallocate (Temperature)   ! print*,""
   deallocate (a)
   deallocate (b)
   deallocate (c)
   deallocate (d)
   ! call solveNodeTempEqn(TA,TB,n)

end program