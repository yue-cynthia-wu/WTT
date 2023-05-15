module utils
real(8) :: t1,t2

contains
   
   subroutine tic()
      implicit none
      call cpu_time(t1)
   end subroutine tic

   subroutine toc()
      implicit none
      call cpu_time(t2)
      ! if (rank==0) print*,"Time Taken -->", real(t2-t1)
      print*,"CPU time -->", real(t2-t1), " seconds"
   end subroutine toc

end module utils