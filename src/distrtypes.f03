module distrtypes
	use types, only : dp, pi
	implicit none
	private
	integer, public, parameter	::	DIST_POWERLAW=0, DIST_LOGNORM=1, DIST_2LOGNORM=2
	public Distribution
contains
	
	subroutine lognormal(F, r, params)
    ! Возвращает значение функции логнормальноо распределегия с параметрам rm, s0 в точке r
		real(dp), intent(in)	::	r(:)
		real(dp), intent(out)	::	F(:)
    real(dp), intent(in)  :: params(:)
		real*8 	rm, s0
		if (size(params) .lt. 2) stop "Недостаточно параметров  для функции!"
    rm = params(1)
    s0 = params(2)
		F = (1.0d0)/(sqrt(2*pi)*abs(log(s0))*r)*exp(-0.5*(r-log(rm))**2/(2*log(s0)**2))
		
	end subroutine lognormal
	
	subroutine PowerLaw(F, r, params)
    ! Возвращает значение функции степенного распределегия с параметрам gamma в точке r
		real(dp), intent(in) :: r(:)
		real(dp), intent(out)	::	F(:)
		real(dp), intent(in)  :: params(:)
    
    if(size(params).lt.1) stop "Недостаточно параметров  для функции!"
    F = r**params(1)
    
	end subroutine PowerLaw
  
  subroutine twoModalLognormal(F, r, params)
    ! Возвращает значение функции 2-модального логнормальноо распределегия с параметрам params в точке r
    real(dp), intent(in) :: r(:), params(:)
    real(dp), intent(out):: F(:)
    real(dp)  :: p1(2), p2(2), p3, tmpF1(size(F)), tmpF2(size(F))
    
    if(size(params).lt.5) stop "Недостаточно параметров  для функции!"
    
    p1(1)=params(1)
    p1(2)=params(2)
    p2(1)=params(3)
    p2(2)=params(4)
    p3=params(5)
    call lognormal(tmpF1, r, p1)
    call lognormal(tmpF2, r, p2)
    F=p3*tmpF1+(1.0_dp-p3)*tmpF2
    
  end subroutine twoModalLognormal
	
	
	subroutine Distribution(F, itype, x, kv, r1, r2, params)
		! Distribution_x - вычитсяет распределение по размерам
		! Входные параметры
		! =================
		! iptype	INTEGER	-	тип распределения (0-степенное, 1-логнормальное, 2-двумодальное логнормальное)
		! x(n)		real*8	-	вектор размерных параметров
		! kv		real*8	-	волновой вектор
		! params(m)	real*8	-	вектор параметров распределения
		! r1, r2	real*8	-	левая и правая границы распределения по радиусу
		!
		! Выходные параметры
		! ==================
		! F(n)		real*8	-	вектор со значениями функции распределения
		!			
		integer, intent(in)	::	itype
		real(dp), intent(in)	::	x(:), params(:), kv, r1, r2
		real(dp), intent(out)	::	F(:)
		real(dp)  :: r(size(x))
    
    r = x/kv
    
		select case (itype)
		case (0) ! Power Law distribution
			call PowerLaw(F, r, params)
		case (1) ! Log-normal distribution
			call Lognormal(F, r, params)
		case (2) ! 2-mode log-normal distribution
			call twoModalLognormal(F, r, params)
		case default
			F(:)=0.0_dp
		end select
		
    where((r.lt.r1).or.(r.gt.r2))
      F=0.0_dp
    end where
    
	end subroutine Distribution
	

end module distrtypes
