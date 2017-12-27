module minimize_funct
  use types, only: dp
  use distrtypes, only: Distribution
  use optprops, only: Kv
  use mathutils, only: convolve4, interp_linear
  use polarization, only: calc_polarization
  implicit none
  
  type costf_params_t
    real(dp)  ::  r0, r1
    integer   ::  idistrtype
    real(dp), allocatable ::  Theta(:)
    real(dp), allocatable ::  Mueller(:,:,:,:)
    real(dp), allocatable ::  X_size(:)
    real(dp), allocatable ::  thtx(:), polx(:,:), Wavelength(:) 
  end type costf_params_t
  
contains
  
  subroutine cost_function1(result, n, x, grad, need_grad, f_data)
    real(dp), intent(out) ::  result
    integer, intent(in) ::  n, need_grad
    real(dp), intent(inout) ::  x(n), grad(n)
    type(costf_params_t), intent(inout) ::  f_data
    
    integer :: M, L, K, I
    real(dp)  ::  F(size(f_data%X_size, dim=1)), Pol(size(f_data%Theta, dim=1))
    real(dp)  ::  fres(size(f_data%thtx, dim=1)), ftmp(size(f_data%thtx, dim=1))
    real(dp), allocatable :: muel_distr(:,:,:)
    real(dp)  ::  chisq
    
    K=size(f_data%Theta, 1)
    M=size(f_data%X_size)
    L=size(f_data%thtx)
    allocate(muel_distr(4,4,K))
    
    chisq = 0.0_dp
    fres=0.0_dp
    
    DO I=1, size(f_data%Wavelength, dim=1)
      muel_distr = 0.0_dp
      ftmp=0.0_dp
      ! Get particles distribution
      call Distribution(F, f_data%idistrtype, f_data%X_size, Kv(f_data%Wavelength(I)), &
      &f_data%r0, f_data%r1, x)
      
      ! perform averaging with weights F
      call convolve4(f_data%Mueller, F, muel_distr)
      ! вычисляем степень поляризации
      call calc_polarization(muel_distr, Pol)
      ! интерполируем на углы измерений экспериментальных данных
      call interp_linear(1, K, f_data%Theta, Pol, L, f_data%thtx, ftmp)
      ! считаем суммарное отклонение по всем длинам волн
      fres=fres+abs(f_data%polx(:,i)-ftmp)*100
      
    END DO
    ! вычисляем сумму квадратов отклонений
    chisq = sum(fres*fres)/size(fres, dim=1)
    !print *,chisq
    result=chisq
    deallocate(muel_distr)
  end subroutine cost_function1
  
  subroutine cost_function(result, n, x, grad, need_grad, f_data)
    real(dp), intent(out) ::  result
    integer, intent(in) ::  n, need_grad
    real(dp), intent(inout) ::  x(n), grad(n)
    type(costf_params_t), intent(inout) ::  f_data
      
    real(dp)  :: eps
    integer I
    
    result = cost_f(n, x, f_data)
    if(need_grad .gt. 0) then
      eps=10*epsilon(x(1))
      !print *, eps
      do i=1, n
        x(i)=x(i)+eps
        grad(i) = (cost_f(n, x, f_data)-result)/eps
        x(i)=x(i)-eps
      end do
    end if
    
  end subroutine cost_function
  
  function cost_f(n,x,f_data)
    integer, intent(in) ::  n
    real(dp), intent(inout) ::  x(n)
    type(costf_params_t), intent(inout) ::  f_data
    real(dp) :: cost_f
    
    integer :: M, L, K, I
    real(dp)  ::  F(size(f_data%X_size, dim=1)), Pol(size(f_data%Theta, dim=1))
    real(dp)  ::  fres(size(f_data%thtx, dim=1)), ftmp(size(f_data%thtx, dim=1))
    real(dp), allocatable :: muel_distr(:,:,:)
    real(dp)  ::  chisq
    
    K=size(f_data%Theta, 1)
    M=size(f_data%X_size)
    L=size(f_data%thtx)
    allocate(muel_distr(4,4,K))
    
    chisq = 0.0_dp
    fres=0.0_dp
    
    DO I=1, size(f_data%Wavelength, dim=1)
      muel_distr = 0.0_dp
      ftmp=0.0_dp
      ! Get particles distribution
      call Distribution(F, f_data%idistrtype, f_data%X_size, Kv(f_data%Wavelength(I)), &
      &f_data%r0, f_data%r1, x)
      
      ! perform averaging with weights F
      call convolve4(f_data%Mueller, F, muel_distr)
      ! вычисляем степень поляризации
      call calc_polarization(muel_distr, Pol)
      ! интерполируем на углы измерений экспериментальных данных
      call interp_linear(1, K, f_data%Theta, Pol, L, f_data%thtx, ftmp)
      ! считаем суммарное отклонение по всем длинам волн
      fres=fres+abs(f_data%polx(:,i)-ftmp)*100
      
    END DO
    ! вычисляем сумму квадратов отклонений
    chisq = sum(fres*fres)/size(fres, dim=1)
    !print *,chisq
    cost_f=chisq
    deallocate(muel_distr)
    
  end function cost_f
  
  subroutine get_polarization_and_f(np,x,F,Pol, f_data)
    integer, intent(in) ::  np
    real(dp), intent(in) :: x(np)
    real(dp), intent(out):: F(:), Pol(:)
    type(costf_params_t), intent(inout) ::  f_data
    real(dp), allocatable :: muel_distr(:,:,:)
    integer N, M
    N=size(F, dim=1)
    M=size(Pol, dim=1)
    allocate(muel_distr(4,4,M))
    
    
    call Distribution(F, f_data%idistrtype, f_data%X_size, Kv(f_data%Wavelength(1)), &
    &f_data%r0, f_data%r1, x)
    ! perform averaging with weights F
    call convolve4(f_data%Mueller, F, muel_distr)
    ! вычисляем степень поляризации
    call calc_polarization(muel_distr, Pol)
    
    deallocate(muel_distr)
    
  end subroutine get_polarization_and_f

end module minimize_funct
