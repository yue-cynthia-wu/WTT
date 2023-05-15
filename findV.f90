subroutine findV(k0,k1,k2,Kv0,Kv1,Kv2,o0,o1,o2,V)
   use, intrinsic :: iso_fortran_env, d=>real64
   use header,only : p_g,p_f,p_N
   implicit none
   
   real(d), intent(in) :: k0,k1,k2,Kv0(2),Kv1(2),Kv2(2),o0,o1,o2
   complex(d), intent(out) :: V
   real(d) :: Kv0T(2),Kv1T(2),Kv2T(2),fI,fJ
   complex(d), parameter :: i=(0,1) ! sqrt(-1)  
   complex(d) :: fK

   Kv0T=[ -Kv0(2),Kv0(1) ] 
   Kv1T=[ -Kv1(2),Kv1(1) ] 
   Kv2T=[ -Kv2(2),Kv2(1) ] 

   fI=-(sqrt(o1*o2/o0)*k0**2*dot_product(Kv1,Kv2)  &
       +sqrt(o0*o2/o1)*k1**2*dot_product(Kv0,Kv2)  &
       +sqrt(o1*o0/o2)*k2**2*dot_product(Kv1,Kv0))

   fJ=p_f**2/sqrt(o0*o1*o2)*(k0**2*dot_product(Kv1,Kv2)  &
                            -k1**2*dot_product(Kv0,Kv2)  &
                            -k2**2*dot_product(Kv1,Kv0))

   fK=i*p_f*(sqrt(o0/(o1*o2))*(k1**2-k2**2)*dot_product(Kv1,Kv2T)  &
          +sqrt(o1/(o2*o0))*(k2**2-k0**2)*dot_product(Kv2,Kv0T)  &
          +sqrt(o2/(o0*o1))*(k0**2-k1**2)*dot_product(Kv0,Kv1T))

   V=p_N/(4*sqrt(2*p_g))/(k0*k1*k2)*(fI+fJ+fK)
   
end subroutine findV