subroutine findIntRed (k0,k1,k2,Kv0,Kv1,Kv2,m1,m2,o0,o1,o2,A0,A1,A2,S,int)
   use, intrinsic :: iso_fortran_env, d=>real64
   use header,only : pi, p_f, p_factor, p_error
   implicit none
   
   real(d), intent(in)  :: k0,k1,k2,Kv0,Kv1,Kv2,m1,m2,o0,o1,o2,A0,A1,A2,S
   real(d), intent(out) :: int
   real(d) :: f102, hm1, gpm1
   complex(d) :: V102

   if (abs( -Kv0+Kv1-Kv2 )>p_error) then
      write (*,'( "   Error: -Kv0+Kv1-Kv2 = ", E10.4, " Stopped!" )') -Kv0+Kv1-Kv2
      write (*,*) "   Wavenumber resonance not satisfied!"
      write (*,*) "   Press any key to exit:"
      read  (*,*)
      stop     
   end if

   if (abs(  -o0 +o1 -o2 )>p_error) then
      write (*,'( "   Error:  -o0 +o1 -o2 = ", E10.4, " Stopped!" )')  -o0 +o1 -o2
      write (*,*) "   Frequency resonance not satisfied!"
      write (*,*) "   Press any key to exit:"
      read  (*,*)
      stop
   end if

   f102=A0*A2-A1*(A0+A2)
   call findV(k1,k0,k2,Kv1,Kv0,Kv2,o1,o0,o2,V102)

   hm1 = f102*abs(V102)**2 *k1*k2/S
   gpm1= m1/o1*(o1**2-p_f**2)/((k1/p_factor)**2+m1**2) &
        -m2/o2*(o2**2-p_f**2)/((k2/p_factor)**2+m2**2)
   int= -8 *pi*hm1/abs(gpm1)

end subroutine findIntRed