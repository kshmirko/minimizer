module polarization
  use types, only: dp
	implicit none
	private
	
	public calc_polarization, calc_lindeprat
	
contains
	
	subroutine calc_polarization(M, P)
		! Вычисляет спепень поляризации
		! Входные параметры
		! =================
		! M(K,K,M) - матрица мюллера
		! 
		! Выходные параметры
		! ==================
		! P(M) - степень поляризации
		real(dp), intent(in) :: M(:,:,:)
		real(dp), intent(out):: P(:)
		
		P(:) = -M(2,1,:)/M(1,1,:)
		
	end subroutine calc_polarization
	
	subroutine calc_lindeprat(M, MUL)
		! Вычисляет линейное деполяризационное отношение
		! Входные параметры
		! =================
		! M(K,K,M) - матрица мюллера
		! 
		! Выходные параметры
		! ==================
		! MUL(M) - линейное деполяризационное отношение
		real(dp), intent(in) :: M(:,:,:)
		real(dp), intent(out):: MUL(:)
		
		MUL(:) = (M(1,1,:)-M(2,2,:))/(M(1,1,:)+2*M(1,2,:)+M(2,2,:))
		
	end subroutine calc_lindeprat

end module polarization
