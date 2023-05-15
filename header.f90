MODULE header
use, intrinsic :: iso_fortran_env, d=>real64
implicit none

integer, parameter :: p_GM=76, p_hydro=0, iRollK=0, iRollM=0
real(d), parameter :: MaxK=0.016, MaxM_v=0.32, p_kc=0.03, p_mc=0.6
! ! -----------------------------------------------------------------------------
! ! real(d), parameter :: pi=4.D0*DATAN(1.D0), p_error=1e-7
real(d), parameter :: pi=3.1415926535898, p_error=1e-7
real(d), parameter :: p_g=9.81, p_rho=1e3
! real(d), parameter :: p_f=7.83608494805175e-05, p_N=0.00523598775598299
real(d), parameter :: p_f=1e-4, p_N=5e-3
real(d) :: p_factor=p_rho*p_N**2/p_g
real(d) :: MaxM    =MaxM_v*p_g/(p_rho*p_N**2)
real(d) :: MaxMroot=MaxM_v*p_g/(p_rho*p_N**2)
! -----------------------------------------------------------------------------
! Selection:
real(d), parameter :: ro_PSI1=.1, ro_PSI2=.1, ro_PSI3=.1
real(d), parameter :: rm_PSI1=1,  rm_PSI2=2,  rm_PSI3=5
real(d), parameter :: ro_ES1 =1,  ro_ES2 =2,  ro_ES3 =5
real(d), parameter :: rm_ES1 =.1, rm_ES2 =.1, rm_ES3 =.1
real(d), parameter :: ro_ID1 =1,  ro_ID2 =2,  ro_ID3 =5
real(d), parameter :: rm_ID1 =1,  rm_ID2 =2,  rm_ID3 =5
! -----------------------------------------------------------------------------
! Count triads:  
integer :: count_sum_ml_in=0,count_sum_ml_out=0,count_sum_mr_in=0,count_sum_mr_out=0
integer :: count_red_ml_in=0,count_red_ml_out=0,count_red_mr_in=0,count_red_mr_out=0
   
END MODULE header
