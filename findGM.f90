subroutine findGM(k,mI,o,nkm)
   use, intrinsic :: iso_fortran_env, d=>real64 
   use header,only : pi,p_f,p_N,p_factor,p_GM,p_hydro,iRollK,iRollM,p_kc,p_mc,MaxK,MaxM
   implicit none
   
   real(d), intent(in) :: k,mI,o
   real(d), intent(out):: nkm ! action(k,m)
!    real(d), parameter  :: N0=0.00523598775598299, b=1300, E0=3e-3
   real(d), parameter  :: N0=0.005, b=1300, E0=3e-3
   real(d) :: m,mstar,s,t,jstar,funcA,funcB,kslope,mslope,Eom
   
   if (p_GM==75) then
      s=1
      t=2.5
      jstar=6
   else
      s=2
      t=2
      jstar=4 ! changed from 3 to 4
   end if 

   m=abs(mI)*p_factor
!    mstar=pi*jstar/b*p_N/N0
   mstar=0.01
   funcA=s*gamma(t/s)/gamma(1.d0/s)/gamma((t-1.d0)/s) *(1.d0+(m/mstar)**s)**(-t/s) /mstar
   
   if (iRollM==1 .and. m>p_mc) then ! modified GM with roll-off spectra
      mslope=(m-p_mc)/(MaxM-p_mc)*(-2)
      funcA=funcA*(m/p_mc)**mslope
   end if
   
   funcB=(2*p_f/pi)/o/sqrt(o**2-p_f**2)
   Eom=p_N/N0*E0*funcA*funcB

   if (p_hydro==1) then ! hydrostatic
      nkm=Eom/(2.d0*pi) *p_N**2*(o*m)**(-2)
   else ! nonhydrostatic
      nkm=Eom/(2.d0*pi) * (p_N**2-o**2)*(k**2+m**2)**(-1)*o**(-2)
   end if 

   if (iRollK==1 .and. k>p_kc) then ! modified GM with roll-off spectra
      kslope=(k-p_kc)/(MaxK-p_kc)*(-2)
      nkm=nkm*(k/p_kc)**kslope
   end if

end subroutine findGM