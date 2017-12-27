module dataio

	implicit none
	private
	
	public read2cols, print_expdata
	
contains

	subroutine read2cols(fname, N, M, X, Y)
		character(len=*), intent(in) :: fname
		integer, intent(in)	:: N, M
		double precision, intent(out)	::	X(N), Y(N,M)
		integer I,J, unit
		
		open(newunit=unit, FILE=fname, status='old')
		
		do I=1, N
			read(unit, *) X(I), (Y(I,J), J=1, M)
		end do 
		
		close(unit)
	end subroutine read2cols
	
	subroutine print_expdata(N, M, X, Y)
		integer, intent(in)	::	N, M
		real*8, intent(in)	:: X(N), Y(N,M)
		integer i,j
		
		print '(A5,5I7)','THTA',(i,i=1,M)
		do i=1, n
			write(*,'(F5.1)', advance="no") X(i)
			write(*,'(5F7.3)', advance="no") (Y(i,j), j=1, M)
			write(*,*)
		end do
	end subroutine print_expdata
end module dataio
