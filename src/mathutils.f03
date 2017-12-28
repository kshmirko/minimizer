module mathutils
	use types, only: dp, pi
	implicit none
	private
	public convolve4, interp_linear, interp_linear2, PI
contains

	subroutine convolve4(M, A, R)
		! Свертка матрицы мюллера с распределением по размерам
		! Входные параметры
		! =================
		! M(K,K,M,N) - матрица мюллера
		! A(N) - функция распределения
		!
		! Выходные параметры
		! ==================
		! R(K,K,M) - результат свертки матрицы мюллера с распределением 
		real(dp), intent(in) :: M(:,:,:,:), A(:)
		real(dp), intent(out)	::	R(:,:,:)
		integer i,NN,NM,NK,ND
		double precision, allocatable	:: tmp(:)
		
		! выделяем память под промежуточный массив
		NN=size(M,4)
    NM=size(M,3)
    NK=size(M,2)
    ND=NK*NK*NM
    allocate(tmp(ND))
		
    !свертка - это умножение матрицы на вектор
    !для этого решейпим исходну матрицу мюллера, объединяем первые 3 размерности
    !M(4,4,M,N)==>M(4*4*M,N)
    !умножаем полученную матрицу на вектор длины N
    !результат помещаем в вектор tmp(16*M)
    call dgemv('N', ND, NN, 1.0_dp, reshape(M, (/ND, NN/)), ND, A, 1, 0.0_dp, tmp, 1)
		
    !возвращаем исходную форму матрицы
		R(:,:,:) = reshape(tmp, (/NK,NK,NM/))
		
		deallocate(tmp)
	end subroutine convolve4
	
	subroutine r8vec_bracket ( n, x, xval, left, right )

	!*****************************************************************************80
	!
	!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
	!
	!  Discussion:
	!
	!    An R8VEC is an array of double precision real values.
	!
	!    If the values in the vector are thought of as defining intervals
	!    on the real line, then this routine searches for the interval
	!    nearest to or containing the given value.
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license.
	!
	!  Modified:
	!
	!    06 April 1999
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    Input, integer ( kind = 4 ) N, length of input array.
	!
	!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
	!
	!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
	!
	!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
	!    Either:
	!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
	!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
	!    or
	!      X(LEFT) <= XVAL <= X(RIGHT).
	!
	  implicit none

	  integer ( kind = 4 ), intent(in) :: n

	  integer ( kind = 4 ) i
	  integer ( kind = 4 ), intent(out)	::	left
	  integer ( kind = 4 ), intent(out)	:: right
	  real ( kind = 8 ), intent(in) :: x(n)
	  real ( kind = 8 ), intent(in) :: xval

	  do i = 2, n - 1

	    if ( xval < x(i) ) then
	      left = i - 1
	      right = i
	      return
	    end if

	   end do

	  left = n - 1
	  right = n

	  return
	end subroutine r8vec_bracket
	
	function r8vec_ascends_strictly ( n, x )

	!*****************************************************************************80
	!
	!! R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.
	!
	!  Discussion:
	!
	!    An R8VEC is a vector of R8 values.
	!
	!    Notice the effect of entry number 6 in the following results:
	!
	!      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
	!      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
	!      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )
	!
	!      R8VEC_ASCENDS_STRICTLY ( X ) = FALSE
	!      R8VEC_ASCENDS_STRICTLY ( Y ) = FALSE
	!      R8VEC_ASCENDS_STRICTLY ( Z ) = TRUE
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license.
	!
	!  Modified:
	!
	!    03 December 2007
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    Input, integer ( kind = 4 ) N, the size of the array.
	!
	!    Input, real ( kind = 8 ) X(N), the array to be examined.
	!
	!    Output, logical R8VEC_ASCENDS_STRICTLY, is TRUE if the
	!    entries of X strictly ascend.
	!
	  implicit none

	  integer ( kind = 4 ), intent(in)	:: n

	  integer ( kind = 4 ) i
	  logical r8vec_ascends_strictly
	  real ( kind = 8 ), intent(in)	:: x(n)

	  do i = 1, n - 1
	    if ( x(i+1) <= x(i) ) then
	      r8vec_ascends_strictly = .false.
	      return
	    end if
	  end do

	  r8vec_ascends_strictly = .true.

	  return
	end
	
	subroutine interp_linear ( m, data_num, t_data, p_data, interp_num, t_interp, p_interp )

	!*****************************************************************************80
	!
	!! INTERP_LINEAR: piecewise linear interpolation to a curve in M dimensions.
	!
	!  Discussion:
	!
	!    From a space of M dimensions, we are given a sequence of
	!    DATA_NUM points, which are presumed to be successive samples
	!    from a curve of points P.
	!
	!    We are also given a parameterization of this data, that is,
	!    an associated sequence of DATA_NUM values of a variable T.
	!    The values of T are assumed to be strictly increasing.
	!
	!    Thus, we have a sequence of values P(T), where T is a scalar,
	!    and each value of P is of dimension M.
	!
	!    We are then given INTERP_NUM values of T, for which values P
	!    are to be produced, by linear interpolation of the data we are given.
	!
	!    Note that the user may request extrapolation.  This occurs whenever
	!    a T_INTERP value is less than the minimum T_DATA or greater than the
	!    maximum T_DATA.  In that case, linear extrapolation is used.
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license.
	!
	!  Modified:
	!
	!    03 December 2007
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    Input, integer ( kind = 4 ) M, the spatial dimension.
	!
	!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
	!
	!    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the value of the
	!    independent variable at the sample points.  The values of T_DATA
	!    must be strictly increasing.
	!
	!    Input, real ( kind = 8 ) P_DATA(M,DATA_NUM), the value of the
	!    dependent variables at the sample points.
	!
	!    Input, integer ( kind = 4 ) INTERP_NUM, the number of points
	!    at which interpolation is to be done.
	!
	!    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the value of the
	!    independent variable at the interpolation points.
	!
	!    Output, real ( kind = 8 ) P_INTERP(M,DATA_NUM), the interpolated
	!    values of the dependent variables at the interpolation points.
	!
	  implicit none

	  integer ( kind = 4 ), intent(in)	:: data_num
	  integer ( kind = 4 ), intent(in)	:: m
	  integer ( kind = 4 ), intent(in)	:: interp_num

	  integer ( kind = 4 ) interp
	  integer ( kind = 4 ) left
	  real ( kind = 8 ), intent(in)	:: p_data(m,data_num)
	  real ( kind = 8 ), intent(out):: p_interp(m,interp_num)
	  integer ( kind = 4 ) right
	  real ( kind = 8 ) t
	  real ( kind = 8 ), intent(in)	:: t_data(data_num)
	  real ( kind = 8 ), intent(in)	:: t_interp(interp_num)
	  
	  if ( .not. r8vec_ascends_strictly ( data_num, t_data ) ) then
	    write ( *, '(a)' ) ' '
	    write ( *, '(a)' ) 'INTERP_LINEAR - Fatal error!'
	    write ( *, '(a)' ) &
	      '  Independent variable array T_DATA is not strictly increasing.'
	    stop 1
	  end if

	  do interp = 1, interp_num

	    t = t_interp(interp)
	!
	!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
	!  nearest to, TVAL.
	!
	    call r8vec_bracket ( data_num, t_data, t, left, right )

	    p_interp(1:m,interp) = &
	      ( ( t_data(right) - t                ) * p_data(1:m,left)   &
	      + (                 t - t_data(left) ) * p_data(1:m,right) ) &
	      / ( t_data(right)     - t_data(left) )

	  end do

	  return
	end
	
	subroutine interp_linear2 ( n, m, data_num, t_data, p_data, interp_num, t_interp, p_interp )

	!*****************************************************************************80
	!
	!! INTERP_LINEAR2: piecewise linear interpolation to a curve in M dimensions.
	!
	!  Discussion:
	!
	!    From a space of M dimensions, we are given a sequence of
	!    DATA_NUM points, which are presumed to be successive samples
	!    from a curve of points P.
	!
	!    We are also given a parameterization of this data, that is,
	!    an associated sequence of DATA_NUM values of a variable T.
	!    The values of T are assumed to be strictly increasing.
	!
	!    Thus, we have a sequence of values P(T), where T is a scalar,
	!    and each value of P is of dimension M.
	!
	!    We are then given INTERP_NUM values of T, for which values P
	!    are to be produced, by linear interpolation of the data we are given.
	!
	!    Note that the user may request extrapolation.  This occurs whenever
	!    a T_INTERP value is less than the minimum T_DATA or greater than the
	!    maximum T_DATA.  In that case, linear extrapolation is used.
	!
	!  Licensing:
	!
	!    This code is distributed under the GNU LGPL license.
	!
	!  Modified:
	!
	!    03 December 2007
	!
	!  Author:
	!
	!    John Burkardt
	!
	!  Parameters:
	!
	!    Input, integer ( kind = 4 ) M, the spatial dimension.
	!
	!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
	!
	!    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the value of the
	!    independent variable at the sample points.  The values of T_DATA
	!    must be strictly increasing.
	!
	!    Input, real ( kind = 8 ) P_DATA(M,DATA_NUM), the value of the
	!    dependent variables at the sample points.
	!
	!    Input, integer ( kind = 4 ) INTERP_NUM, the number of points
	!    at which interpolation is to be done.
	!
	!    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the value of the
	!    independent variable at the interpolation points.
	!
	!    Output, real ( kind = 8 ) P_INTERP(M,DATA_NUM), the interpolated
	!    values of the dependent variables at the interpolation points.
	!
	  implicit none

	  integer ( kind = 4 ), intent(in)	:: data_num
	  integer ( kind = 4 ), intent(in)	:: m, n
	  integer ( kind = 4 ), intent(in)	:: interp_num

	  integer ( kind = 4 ) interp
	  integer ( kind = 4 ) left
	  real ( kind = 8 ), intent(in)	:: p_data(n, m, data_num)
	  real ( kind = 8 ), intent(out)	:: p_interp(n, m, interp_num)
	  integer ( kind = 4 ) right
	  real ( kind = 8 ) t
	  real ( kind = 8 ), intent(in)	:: t_data(data_num)
	  real ( kind = 8 ), intent(in)	:: t_interp(interp_num)

	  if ( .not. r8vec_ascends_strictly ( data_num, t_data ) ) then
	    write ( *, '(a)' ) ' '
	    write ( *, '(a)' ) 'INTERP_LINEAR2 - Fatal error!'
	    write ( *, '(a)' ) &
	      '  Independent variable array T_DATA is not strictly increasing.'
	    stop 1
	  end if

	  do interp = 1, interp_num

	    t = t_interp(interp)
	!
	!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
	!  nearest to, TVAL.
	!
	    call r8vec_bracket ( data_num, t_data, t, left, right )

	    p_interp(1:n,1:m,interp) = &
	      ( ( t_data(right) - t                ) * p_data(1:n,1:m,left)   &
	      + (                 t - t_data(left) ) * p_data(1:n,1:m,right) ) &
	      / ( t_data(right)     - t_data(left) )

	  end do

	  return
	end
	
end module mathutils
