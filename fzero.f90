module fzero

use, intrinsic :: iso_fortran_env, d=>real64
use header,only : p_f,p_N,p_factor,p_error,MaxMroot
implicit none

private    o0i,m0i,k1i,k2i,iDisp
real(d) :: o0i,m0i,k1i,k2i,iDisp=1
real(d) :: m1lsum, m1rsum, m1lred, m1rred
real(d) :: m2lsum, m2rsum, m2lred, m2rred
real(d) :: o1lsum, o1rsum, o1lred, o1rred
real(d) :: o2lsum, o2rsum, o2lred, o2rred
integer :: exist_mlsum,exist_mrsum,exist_mlred,exist_mrred

contains
   
   subroutine setValues(o0s,m0s,k1s,k2s)
      implicit none
      real(d), intent(in):: o0s,m0s,k1s,k2s
      o0i = o0s
      m0i = m0s
      k1i = k1s
      k2i = k2s
   end subroutine setValues
   
   !! ===========================================================================================
   subroutine rootSum()
      implicit none
      if ((gm1sum(0.d0)<0).and.(gm1sum(-MaxMroot+m0i)>0)) then
         m1lsum = brent0(gm1sum, -MaxMroot+m0i,0.d0, p_error,detail=0)
         m2lsum = m0i-m1lsum
         o1lsum = omega(k1i,m1lsum,p_f,p_N,p_factor)
         o2lsum = omega(k2i,m2lsum,p_f,p_N,p_factor)
         if ( abs( o0i-o1lsum-o2lsum)>p_error .and. iDisp==1.) then
            write(*,'( "   Error: o0-o1lsum-o2lsum = ", E10.4, " Stopped!" )')  o0i-o1lsum-o2lsum
            write(*,*) "Press any key to exit:"
            read (*,*)
            stop
         end if
!          if ((abs(m1lsum*p_factor)>=k1i).and.(abs(m2lsum*p_factor)>=k2i)) then
            exist_mlsum=1
!          else
!             exist_mlsum=0
!          end if
      else
         exist_mlsum=0 ! must have this line because exist_mlsum may come from last loop
      end if
      
!       if ((gm1sum(m0i)<0).and.(gm1sum(MaxMroot)>0)) then
!          m1rsum = brent0(gm1sum, m0i,MaxMroot, error,detail=0)
!          m2rsum = m0i-m1rsum
!          o1rsum = omega(k1i,m1rsum,p_f,p_N,p_factor)
!          o2rsum = omega(k2i,m2rsum,p_f,p_N,p_factor)
!          if ( abs( o0i-o1rsum-o2rsum)>p_error .and. iDisp==1.) then
!             write(*,'( "   Error: o0-o1rsum-o2rsum = ", E10.4, " Stopped!" )')  o0i-o1rsum-o2rsum
!             write(*,*) "Press any key to exit:"
!             read (*,*)
!             stop
!          end if
!          if ((abs(m1rsum*p_factor)>=k1i).and.(abs(m2rsum*p_factor)>=k2i)) then
!             exist_mrsum=1
!          else
!             exist_mrsum=0
!          end if
!       else
!          exist_mrsum=0
!       end if
   end subroutine rootSum

   !! ===========================================================================================
	subroutine rootRed()
   	implicit none
      if ((gm1red(0.d0)>0).and.(gm1red(-MaxMroot+m0i)<0)) then
         m1lred = brent0(gm1red, -MaxMroot+m0i,0.d0, p_error,detail=0)
         m2lred = m1lred-m0i
         o1lred = omega(k1i,m1lred,p_f,p_N,p_factor)
         o2lred = omega(k2i,m2lred,p_f,p_N,p_factor)
         if ( abs(-o0i+o1lred-o2lred)>p_error .and. iDisp==1.) then
            write(*,'( "  Error: -o0+o1lred-o2lred = ", E10.4, " Stopped!" )') -o0i+o1lred-o2lred
            write(*,*) "Press any key to exit:"
            read (*,*)
            stop
         end if
!          if ((abs(m1lred*p_factor)>=k1i).and.(abs(m2lred*p_factor)>=k2i)) then
            exist_mlred=1
!          else
!             exist_mlred=0
!          end if
      else
         exist_mlred=0
      end if
      
      if ((gm1red(0.d0)>0).and.(gm1red(m0i)<0)) then
         m1rred = brent0(gm1red, 0.d0,m0i, p_error,detail=0)
         m2rred = m1rred-m0i
         o1rred = omega(k1i,m1rred,p_f,p_N,p_factor)
         o2rred = omega(k2i,m2rred,p_f,p_N,p_factor)
         if ( abs(-o0i+o1rred-o2rred)>p_error .and. iDisp==1.) then
            write(*,'( "  Error: -o0+o1rred-o2rred = ", E10.4, " Stopped!" )') -o0i+o1rred-o2rred
            write(*,*) "Press any key to exit:"
            read (*,*)
            stop
         end if
!          if ((abs(m1rred*p_factor)>=k1i).and.(abs(m2rred*p_factor)>=k2i)) then
            exist_mrred=1
!          else
!             exist_mrred=0
!          end if
      else
         exist_mrred=0
      end if 
	end subroutine rootRed 
   
   !! ===========================================================================================
   pure function gm1sum(x)
   	implicit none
   	real(d), intent(in) :: x
   	real(d) :: gm1sum
      gm1sum= o0i-omega(k1i,x,p_f,p_N,p_factor)-omega(k2i,m0i-x,p_f,p_N,p_factor)
	end function gm1sum
    
   pure function gm1red(x)
   	implicit none
   	real(d), intent(in) :: x
   	real(d) :: gm1red
      gm1red=-o0i+omega(k1i,x,p_f,p_N,p_factor)-omega(k2i,x-m0i,p_f,p_N,p_factor)
	end function gm1red
    
   pure function omega(k,mI,f,N,factor)
      implicit none
      real(d), intent(in) :: k, mI, f, N, factor ! isopycnal wavenumber
      real(d) :: omega, m ! vertical wavenumber
      m=factor*mI
      omega=sqrt(  ( (N*k)**2+(f*m)**2 )/( k**2+m**2 )  )
   end function omega
    
   !! ===========================================================================================
	function brent0(fun, x1, x2, tol, detail)

	implicit none
	real(d) :: brent0
	real(d), intent(IN) :: x1, x2, tol
	real(d), external :: fun
	integer, intent(IN), optional :: detail
	integer :: i, exitflag!, disp
	real(d) :: a, b, c, diff,e, fa, fb, fc, p, q, r, s, tol1, xm
	real(d), parameter :: EPS = epsilon(a)
	integer, parameter :: imax = 100  ! maximum number of iteration
	
	exitflag = 0
		
	! intialize values
	a = x1
	b = x2
	c = x2
	fa = fun(a)
	fb = fun(b)
	fc=fb
		
	! check sign
	if ( (fa>0. .and. fb>0. )  .or.  (fa>0. .and. fb>0. )) then
		write(*,*)  'Error (brent.f90): Root must be bracked by two imputs'
		write(*, "('f(x1) = ', 1F8.4, '   f(x2) = ', 1F8.4)") fa,fb
		stop
	end if
	
	! main iteration
	do i = 1, imax
		! rename c and adjust bounding interval if both a(=b) and c are same sign
		if ((fb > 0.  .and. fc > 0) .or. (fb <0. .and. fc < 0. ) ) then 
			c = a
			fc = fa
			e = b-a
			diff = e
		end if
		
		! if c is better guess than b, use it. 
		if (abs(fc) < abs(fb) ) then 
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc= fa
		end if
		
		! convergence check
		tol1=2.0_d* EPS * abs(b) + 0.5_d*tol
		xm = 0.5_d * (c - b)
		if (abs(xm) < tol1 .or. fb == 0.0_d )  then
			exitflag = 1
			exit
		end if
	   
		! try inverse quadratic interpolation
		if (abs(e) >= tol1 .and. abs(fa) > abs(fb) ) then 
			s = fb/fa
				if (abs(a - c) < EPS) then 
				p = 2.0_d *xm * s
				q = 1.0_d  - s
			else
				q = fa/fc
				r = fb/fc
				p = s * (2.0_d * xm * q * (q -r ) - (b - a) * (r - 1.0_d))
				q = (q - 1.0_d ) * (r - 1.0_d) * (s - 1.0_d) 
			end if
			
			! accept if q is not too small to stay in bound
			if (p > 0.0_d) q = -q
			p = abs(p)                
			if (2.0 * p < min(3.0 * xm * q - abs(tol1* q), abs(e *q))) then 
				e = d
				diff = p / q
			else   ! interpolation failed. use bisection
				diff= xm 
				e = d
			end if
		else  ! quadratic interpolation bounds moves too slowly, use bisection
			diff = xm
			e = d
		end if
		
		! update last bound
		a = b
		fa = fb
		
		! move the best guess
		if (abs(d) > tol1) then 
			b = b + diff
		else
			b = b + sign(tol1, xm)
		end if
		
		! evaluate new trial root
		fb = fun(b)        
	end do 
	
	! case for non convergence
	if (exitflag /= 1 ) then 
		write(*,*) 'Error (brent.f90) :  convergence was not attained'
		write(*,*) 'Initial value:'
		write(*,"(4F10.5)" )   x1, x2, fun(x1), fun(x2)
		write(*,*) ' '
		write(*,*) 'final value:'
		write(*,"('x = '  ,1F6.4, ':                f(x1) = ' ,  1F6.4  )" )  b,  fb  
	end if
	brent0 = b
	return
	
	end function brent0
   
	
end module fzero