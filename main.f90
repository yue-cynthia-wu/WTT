PROGRAM main
use, intrinsic :: iso_fortran_env, d=>real64
use header
use utils
use mpi
use fzero
implicit none
! -----------------------------------------------------------------------------
! Variables:
include "size.h"
real(d) :: dk,dm,MinK,MinM,gridk(1:Nk),gridm(1:Nm)
real(d) :: k0,k1,k2,Kv0(2),Kv1(2),Kv2(2),m0,o0,A0,a1,a2,p,S
real(d) :: A1lsum, A1lred, A1rred, A2lsum, A2lred, A2rred
real(d) :: intlsum,intlred,intrred
real(d) :: CL_red,CL_sum,CL_PSI1,CL_PSI2,CL_PSI3,CL_ES1,CL_ES2,CL_ES3,CL_ID1,CL_ID2,CL_ID3
! Functions:
character(len=20) :: int2str
real(d) :: findomega
integer :: findPSI,findES,findID
! Other:
integer :: ierr,ncore,id_core,lcore ! MPI
integer :: iarray,im0,iml,ik0,ik1,ik2,iunit ! indices
! Outputs:
real(d),allocatable,dimension(:,:) :: A,At_red,At_sum,At_PSI1,At_PSI2,At_PSI3,At_ES1,At_ES2,At_ES3,At_ID1,At_ID2,At_ID3
allocate(A(Nk,Nm/ncore0))
allocate(At_red,At_sum,At_PSI1,At_PSI2,At_PSI3,At_ES1,At_ES2,At_ES3,At_ID1,At_ID2,At_ID3, mold=A)
! -----------------------------------------------------------------------------
! ! MPI
call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,  ncore,ierr)
call mpi_comm_rank(mpi_comm_world,id_core,ierr)
if (ncore/=ncore0) then
   print *, "   --> Number of cores are not consistent!"
   stop
end if
lcore=Nm/ncore
! -----------------------------------------------------------------------------
! ! Wavenumber space
dk=MaxK/Nk
dm=MaxM/Nm
MinK=dk
MinM=dm
gridk=(/((iarray*dk),iarray=1,Nk)/)
gridm=(/((iarray*dm),iarray=1,Nm)/)
! =============================================================================
!!                                 Iteration
! =============================================================================
call tic()
do im0=1+id_core*lcore, lcore+id_core*lcore
   iml=im0-id_core*lcore
   m0=gridm(im0)
   print *, "   --> calculating im =", im0, " on Processor #", id_core
do ik0=1,Nk
   k0=gridk(ik0)
   Kv0=(/ k0,0.d0 /)
   o0=findomega(k0,m0,p_factor) ! frequency 0
   call findGM(k0,m0,o0,A0)

   CL_sum =0
   CL_red =0
   CL_PSI1=0
   CL_PSI2=0
   CL_PSI3=0
   CL_ES1 =0
   CL_ES2 =0
   CL_ES3 =0
   CL_ID1 =0
   CL_ID2 =0
   CL_ID3 =0

!    if (m0*p_factor>=k0) then
   do ik1=1,Nk
      k1=gridk(ik1) ! magnitude 1
   do ik2=1,Nk
      k2=gridk(ik2) ! magnitude 2
      ! ------------------------------------------------------------------------
      if ((k1+k2-k0>p_error).and.(k0+k2-k1>p_error).and.(k0+k1-k2>p_error)) then
         p=(k0+k1+k2)/2
         S=sqrt(p*(p-k0)*(p-k1)*(p-k2)) ! triangle area
         a1=acos( (k0**2+k1**2-k2**2)/(2*k0*k1) ) ! angle 1
         a2=acos( (k0**2+k2**2-k1**2)/(2*k0*k2) ) ! angle 2
         Kv1=k1*(/  cos(a1), sin(a1) /) ! vector 1
         Kv2=k2*(/  cos(a2),-sin(a2) /) ! vector 2

         call setValues(o0,m0,k1,k2)

         !---------------------------------------------------------------------
         !!                  Summation interaction p2=p0-p1
         !---------------------------------------------------------------------
         call rootSum ! returns both left and right roots

         if (exist_mlsum.eq.1) then ! left root exists
            call findGM (k1,m1lsum,o1lsum,A1lsum) ! action
            call findGM (k2,m2lsum,o2lsum,A2lsum)
            call findIntSum (k0,k1,k2,Kv0,Kv1,Kv2, m1lsum,m2lsum,o0,o1lsum,o2lsum,A0,A1lsum,A2lsum,S,intlsum)

            intlsum=intlsum*2 ! Include the right root
            CL_sum =CL_sum +intlsum
            CL_PSI1=CL_PSI1+intlsum*findPSI(o0,o2lsum,o1lsum,m0,m2lsum,m1lsum, ro_PSI1,rm_PSI1)
            CL_PSI2=CL_PSI2+intlsum*findPSI(o0,o2lsum,o1lsum,m0,m2lsum,m1lsum, ro_PSI2,rm_PSI2)
            CL_PSI3=CL_PSI3+intlsum*findPSI(o0,o2lsum,o1lsum,m0,m2lsum,m1lsum, ro_PSI3,rm_PSI3)
            CL_ES1 =CL_ES1 +intlsum*findES (o0,o2lsum,o1lsum,m0,m2lsum,m1lsum, ro_ES1, rm_ES1 )
            CL_ES2 =CL_ES2 +intlsum*findES (o0,o2lsum,o1lsum,m0,m2lsum,m1lsum, ro_ES2, rm_ES2 )
            CL_ES3 =CL_ES3 +intlsum*findES (o0,o2lsum,o1lsum,m0,m2lsum,m1lsum, ro_ES3, rm_ES3 )
            CL_ID1 =CL_ID1 +intlsum*findID (o0,o2lsum,o1lsum,m0,m2lsum,m1lsum, ro_ID1, rm_ID1 )
            CL_ID2 =CL_ID2 +intlsum*findID (o0,o2lsum,o1lsum,m0,m2lsum,m1lsum, ro_ID2, rm_ID2 )
            CL_ID3 =CL_ID3 +intlsum*findID (o0,o2lsum,o1lsum,m0,m2lsum,m1lsum, ro_ID3, rm_ID3 )
         endif

!          if (exist_mrsum.eq.1) then ! right root exists
!             call findGM (k1,m1rsum,o1rsum,A1rsum) ! action
!             call findGM (k2,m2rsum,o2rsum,A2rsum)
!             call findIntSum (k0,k1,k2,Kv0,Kv1,Kv2, m1rsum,m2rsum,o0,o1rsum,o2rsum,A0,A1rsum,A2rsum,S,intrsum)
!
!             CL_sum =CL_sum +intrsum
!             CL_PSI1=CL_PSI1+intrsum*findPSI(o0,o1rsum,o2rsum,m0,m1rsum,m2rsum, ro_PSI1,rm_PSI1)
!             CL_PSI2=CL_PSI2+intrsum*findPSI(o0,o1rsum,o2rsum,m0,m1rsum,m2rsum, ro_PSI2,rm_PSI2)
!             CL_PSI3=CL_PSI3+intrsum*findPSI(o0,o1rsum,o2rsum,m0,m1rsum,m2rsum, ro_PSI3,rm_PSI3)
!             CL_ES1 =CL_ES1 +intrsum*findES (o0,o1rsum,o2rsum,m0,m1rsum,m2rsum, ro_ES1, rm_ES1 )
!             CL_ES2 =CL_ES2 +intrsum*findES (o0,o1rsum,o2rsum,m0,m1rsum,m2rsum, ro_ES2, rm_ES2 )
!             CL_ES3 =CL_ES3 +intrsum*findES (o0,o1rsum,o2rsum,m0,m1rsum,m2rsum, ro_ES3, rm_ES3 )
!             CL_ID1 =CL_ID1 +intrsum*findID (o0,o1rsum,o2rsum,m0,m1rsum,m2rsum, ro_ID1, rm_ID1 )
!             CL_ID2 =CL_ID2 +intrsum*findID (o0,o1rsum,o2rsum,m0,m1rsum,m2rsum, ro_ID2, rm_ID2 )
!             CL_ID3 =CL_ID3 +intrsum*findID (o0,o1rsum,o2rsum,m0,m1rsum,m2rsum, ro_ID3, rm_ID3 )
!          endif

         !---------------------------------------------------------------------
         !!                  Reduction interaction p2=p1-p0
         !---------------------------------------------------------------------
         Kv2=-Kv2 ! vector 2
         call rootRed ! returns both left and right roots

         if (exist_mlred.eq.1) then ! left root exists
            call findGM (k1,m1lred,o1lred,A1lred) ! action
            call findGM (k2,m2lred,o2lred,A2lred)
            call findIntRed (k0,k1,k2,Kv0,Kv1,Kv2, m1lred,m2lred,o0,o1lred,o2lred,A0,A1lred,A2lred,S,intlred)

            CL_red =CL_red +intlred
            CL_PSI1=CL_PSI1+intlred*findPSI(o1lred,o2lred,o0,m1lred,m2lred,m0, ro_PSI1,rm_PSI1)
            CL_PSI2=CL_PSI2+intlred*findPSI(o1lred,o2lred,o0,m1lred,m2lred,m0, ro_PSI2,rm_PSI2)
            CL_PSI3=CL_PSI3+intlred*findPSI(o1lred,o2lred,o0,m1lred,m2lred,m0, ro_PSI3,rm_PSI3)
            CL_ES1 =CL_ES1 +intlred*findES (o1lred,o2lred,o0,m1lred,m2lred,m0, ro_ES1, rm_ES1 )
            CL_ES2 =CL_ES2 +intlred*findES (o1lred,o2lred,o0,m1lred,m2lred,m0, ro_ES2, rm_ES2 )
            CL_ES3 =CL_ES3 +intlred*findES (o1lred,o2lred,o0,m1lred,m2lred,m0, ro_ES3, rm_ES3 )
            CL_ID1 =CL_ID1 +intlred*findID (o1lred,o2lred,o0,m1lred,m2lred,m0, ro_ID1, rm_ID1 )
            CL_ID2 =CL_ID2 +intlred*findID (o1lred,o2lred,o0,m1lred,m2lred,m0, ro_ID2, rm_ID2 )
            CL_ID3 =CL_ID3 +intlred*findID (o1lred,o2lred,o0,m1lred,m2lred,m0, ro_ID3, rm_ID3 )
         endif

         if (exist_mrred.eq.1) then ! right root exists
            call findGM (k1,m1rred,o1rred,A1rred) ! action
            call findGM (k2,m2rred,o2rred,A2rred)
            call findIntRed (k0,k1,k2,Kv0,Kv1,Kv2, m1rred,m2rred,o0,o1rred,o2rred,A0,A1rred,A2rred,S,intrred)

            CL_red =CL_red +intrred
            CL_PSI1=CL_PSI1+intrred*findPSI(o1rred,o0,o2rred,m1rred,m0,m2rred, ro_PSI1,rm_PSI1)
            CL_PSI2=CL_PSI2+intrred*findPSI(o1rred,o0,o2rred,m1rred,m0,m2rred, ro_PSI2,rm_PSI2)
            CL_PSI3=CL_PSI3+intrred*findPSI(o1rred,o0,o2rred,m1rred,m0,m2rred, ro_PSI3,rm_PSI3)
            CL_ES1 =CL_ES1 +intrred*findES (o1rred,o0,o2rred,m1rred,m0,m2rred, ro_ES1, rm_ES1 )
            CL_ES2 =CL_ES2 +intrred*findES (o1rred,o0,o2rred,m1rred,m0,m2rred, ro_ES2, rm_ES2 )
            CL_ES3 =CL_ES3 +intrred*findES (o1rred,o0,o2rred,m1rred,m0,m2rred, ro_ES3, rm_ES3 )
            CL_ID1 =CL_ID1 +intrred*findID (o1rred,o0,o2rred,m1rred,m0,m2rred, ro_ID1, rm_ID1 )
            CL_ID2 =CL_ID2 +intrred*findID (o1rred,o0,o2rred,m1rred,m0,m2rred, ro_ID2, rm_ID2 )
            CL_ID3 =CL_ID3 +intrred*findID (o1rred,o0,o2rred,m1rred,m0,m2rred, ro_ID3, rm_ID3 )
         endif

      end if ! triangle condition
   end do ! end ik2
   end do ! end ik1
!    end if ! end if m0>=k0
   ! ------------------------------------------------------
   A      (ik0,iml)=A0
   At_sum (ik0,iml)=CL_sum *dk**2
   At_red (ik0,iml)=CL_red *dk**2
   At_PSI1(ik0,iml)=CL_PSI1*dk**2
   At_PSI2(ik0,iml)=CL_PSI2*dk**2
   At_PSI3(ik0,iml)=CL_PSI3*dk**2
   At_ES1 (ik0,iml)=CL_ES1 *dk**2
   At_ES2 (ik0,iml)=CL_ES2 *dk**2
   At_ES3 (ik0,iml)=CL_ES3 *dk**2
   At_ID1 (ik0,iml)=CL_ID1 *dk**2
   At_ID2 (ik0,iml)=CL_ID2 *dk**2
   At_ID3 (ik0,iml)=CL_ID3 *dk**2
end do ! end ik0

! ---------------------------------------------------------
! Write to files for each m0
iunit=(id_core+1)+1100
open (iunit, file='output'//job//'/A_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') A(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+1200
open (iunit, file='output'//job//'/Atsum_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_sum(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+1300
open (iunit, file='output'//job//'/Atred_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_red(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+2100
open (iunit, file='output'//job//'/AtPSI1_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_PSI1(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+2200
open (iunit, file='output'//job//'/AtPSI2_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_PSI2(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+2300
open (iunit, file='output'//job//'/AtPSI3_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_PSI3(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+2400
open (iunit, file='output'//job//'/AtES1_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_ES1(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+2500
open (iunit, file='output'//job//'/AtES2_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_ES2(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+2600
open (iunit, file='output'//job//'/AtES3_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_ES3(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+2700
open (iunit, file='output'//job//'/AtID1_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_ID1(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+2800
open (iunit, file='output'//job//'/AtID2_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_ID2(ik0,iml)
end do
write(iunit,*) ''

iunit=(id_core+1)+2900
open (iunit, file='output'//job//'/AtID3_'//trim(int2str(id_core+1))//'.txt')
do ik0=1,Nk
   write(iunit, '(E20.10,X)', advance='no') At_ID3(ik0,iml)
end do
write(iunit,*) ''
! ---------------------------------------------------------
end do ! end im0

! ---------------------------------------------------------
! ! Save info
if (id_core==0) then
   iunit=(id_core+1)+3000
   open (iunit, file='output'//job//'/0info.txt')
   write(iunit,'(A10,A7)')      "Dir = ",'output'//job
   write(iunit,'(A10,I4)')    "Cores = ",ncore0
   write(iunit,'(A10,I4)')       "Nk = ",Nk
   write(iunit,'(A10,I4)')       "Nm = ",Nm
   write(iunit,'(A10,E8.2)') "max(k) = ",MaxK
   write(iunit,'(A10,E8.2)') "max(m) = ",MaxM*p_factor
   write(iunit,'(A10,E8.2)') "min(k) = ",MinK
   write(iunit,'(A10,E8.2)') "min(m) = ",MinM*p_factor
   write(iunit,*) '-----------------------------------------------'
   write(iunit,'(A10,E8.2)')      "f = ",p_f
   write(iunit,'(A10,E8.2)')      "N = ",p_N
   write(iunit,'(A10,E8.2)') "factor = ",p_factor
   write(iunit,*) '-----------------------------------------------'
   write(iunit,'(A20)')      " Using GM"//int2str(p_GM)//":"
   if (iRollK==1) then
      write(iunit,'(A25,E8.2)') "roll-off at kc = ",p_kc
   endif
   if (iRollM==1) then
      write(iunit,'(A25,E8.2)') "roll-off at mc = ",p_mc
   endif
   write(iunit,*) '----------------------PSI----------------------'
   write(iunit,'(A30,F6.2,F6.2,F6.2)') "frequency-halving error     =",ro_PSI1,ro_PSI2,ro_PSI3
   write(iunit,'(A30,F6.2,F6.2,F6.2)') "wavenumber-separation ratio =",rm_PSI1,rm_PSI2,rm_PSI3
   write(iunit,*) '----------------------ES-----------------------'
   write(iunit,'(A30,F6.2,F6.2,F6.2)') "frequency-separation ratio  =",ro_ES1 ,ro_ES2 ,ro_ES3
   write(iunit,'(A30,F6.2,F6.2,F6.2)') "wavenumber-halving error    =",rm_ES1 ,rm_ES2 ,rm_ES3
   write(iunit,*) '----------------------ID-----------------------'
   write(iunit,'(A30,F6.2,F6.2,F6.2)') "frequency-separation ratio  =",ro_ID1 ,ro_ID2 ,ro_ID3
   write(iunit,'(A30,F6.2,F6.2,F6.2)') "wavenumber-separation ratio =",rm_ID1 ,rm_ID2 ,rm_ID3
   close(iunit)
end if
! ---------------------------------------------------------

call toc()
call mpi_finalize(ierr)

end PROGRAM main 

! =============================================================================
pure function findomega(k,mI,factor) ! isopycnal wavenumber
   use, intrinsic :: iso_fortran_env, d=>real64
   implicit none
   real(d), parameter :: f=1e-4, N=5e-3
   real(d), intent(in) :: k,mI,factor   
   real(d) :: findomega,m
   m=factor*mI
   findomega=sqrt(  ( (N*k)**2+(f*m)**2 )/( k**2+m**2 )  )
   return
end function findomega

pure function int2str(k)
   ! "Convert an integer to string."
   implicit none
   integer, intent(in) :: k
   character(len=20) :: int2str
   write (int2str, *) k
   int2str = adjustl(int2str)
end function int2str

pure function findPSI(ob,ol,or,mb,ml,mr,ro,rm)
   use, intrinsic :: iso_fortran_env, d=>real64
   implicit none
   real(d), intent(in) :: ob,ol,or,mb,ml,mr,ro,rm
   integer :: findPSI
   if ((abs(ol-or)<ob*ro).and.(abs(mr)>abs(mb)*rm).and.(abs(ml)>abs(mb)*rm)) then
      findPSI=1
   else
      findPSI=0
   endif
end

pure function findES(ob,ol,or,mb,ml,mr,ro,rm)
   use, intrinsic :: iso_fortran_env, d=>real64
   implicit none
   real(d), intent(in) :: ob,ol,or,mb,ml,mr,ro,rm
   integer :: findES
   if ((abs(mr+mb)<abs(ml)*rm).and.(ob>ol*ro).and.(or>ol*ro)) then
      findES=1
   else
      findES=0
   endif
end

pure function findID(ob,ol,or,mb,ml,mr,ro,rm)
   use, intrinsic :: iso_fortran_env, d=>real64
   implicit none
   real(d), intent(in) :: ob,ol,or,mb,ml,mr,ro,rm
   integer :: findID
   if ((ob>or*ro).and.(ol>or*ro).and.(abs(mb)>abs(mr)*rm).and.(abs(ml)>abs(mr)*rm)) then
      findID=1
   else
      findID=0
   endif
end