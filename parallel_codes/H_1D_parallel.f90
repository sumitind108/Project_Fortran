
program Heat_equation_1D
   use mpi 
   implicit none 
   integer :: ierr, rank, nprocs
   integer, dimension(MPI_STATUS_SIZE) :: status1
   integer :: data0 = 0, data1 = 0


   call MPI_INIT(ierr)
    
   !Setup communicator size
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

   !setup rank/ids for each process
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   call solveNodeTempEqn()
   !type the main code

   print*, "Each processor work on" ,nprocs, "Process(es)"

   if(rank == 0) then 
      data0 = 100
      call MPI_SEND(data0, 1 , MPI_INIT, 1,1, MPI_COMM_WORLD, ierr)
        print*, "Rank", Rank, "Send data0 to rank 1"


        !to recieve data from rank 1
        call MPI_RECV(data0, 1, MPI_INIT, 1, 2, MPI_COMM_WORLD, status1, ierr)
        print*, "Rank", Rank, "received data1 from rank 1 as data0 =", data0
   end if

   if(rank == 1) then 
      
      call MPI_RECV(data1,1,MPI_INT,0,1, MPI_COMM_WORLD, status1, ierr)
        print*, "Rank", rank, "Recieved data 0 from rank 0 as data1! ="


        data1 = 500
        print*, "Rank", rank, "Modified data 1 the new value of data 1=" , data1

        call MPI_SEND(data1, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, ierr)
        print*,  Rank 

   end if


   !Finalize MPI
   call MPI_FINALIZE(ierr)

   
   

end program

subroutine solveNodeTempEqn()
   integer:: n, i
   real :: TA, TB, thermalConductivity, Area
   real, dimension (:), allocatable :: MatA, MatB, Temperature

   print*, "Enter the values of TA , TB , n , thermalConductivity , Area"
   read* , TA , TB , n , thermalConductivity , Area

   allocate (MatA(n*n))
   allocate (MatB(n))
   allocate (Temperature(n))
   
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

   call TDMA(MatA,MatB,Temperature,n)

   print*,"------------------------------"
   print* ,"Temperature:"
   do i=1,n
      print* ,Temperature(i)
   end do
   print*,"------------------------------"

   deallocate (MatA)
   deallocate (MatB)
   deallocate (Temperature)

end subroutine

subroutine TDMA(MatA,MatB,Temperature,n)

   integer:: n,i,j
   real:: MatA(n*n),MatB(n),a(n),b(n),c(n),d(n),Temperature(n),q
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
end subroutine