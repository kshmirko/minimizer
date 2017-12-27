module optprops
  use types, only : dp, pi
  implicit none
  private
  public Kv
contains
  
  function Kv(wavelen)
    ! Вычисляет волновой вектор на основе длины волны (в мкм)
    ! Входные параметры:
    ! ==================
    ! wavelen real(dp)  - длина волны  мкм
    ! Возвращает значение
    ! real(dp) волновой вектор
    real(dp), intent(in) :: wavelen
    real(dp) :: Kv
    
    Kv = 2.0_dp*pi/wavelen
    return
  end function Kv

end module optprops
