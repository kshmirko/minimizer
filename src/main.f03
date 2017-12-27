program main
  use types, only : dp, MAX_PATH
  use distrtypes, only: Distribution
  use optprops, only: Kv
  use ncdfapi, only: get_mueller_from_netcdf
  use minimize_funct, only: costf_params_t, cost_function, get_polarization_and_f
  use dataio, only: read2cols, print_expdata
  implicit none
  include 'nlopt.f'
  
  integer(kind=4) ::  nFiles, nWavelengths, nExps, nPars
  character(len=MAX_PATH), allocatable ::  Files(:)
  character(len=MAX_PATH) :: fmtString
  character(len=MAX_PATH) :: ExpDataFName, SAVEFDIST, SAVEPOL
  real(dp), allocatable ::  Wavelength(:)
  real(dp), allocatable ::  x0(:), ub(:), lb(:), xopt(:), xtmp(:), F(:), Pol(:), xx(:)
  !real(dp), allocatable ::  X(:), Theta(:), Mueller(:,:,:,:)
  real(dp)  ::  r0, r1, xtol, ftol, minf, k_v
  type(costf_params_t)  ::  f_data
  complex*16 ::  M1, M2
  
  
  integer unit, I, OPTALG
  ! Параметры, необходимые для минимизации
  integer*8 opt
  integer ires, idistrtype
  real(dp)  ::  optimval
  integer   ::  ioptimval
  ! Список параметров, читаемых из файла
  namelist /header/ nFiles, nWavelengths, nExps, nPars, OPTALG
  namelist /config/ Files, ExpDataFName, Wavelength, x0, lb, ub, r0, r1, xtol, ftol, idistrtype, SAVEFDIST, SAVEPOL
  
  open(newunit=unit, FILE='config.ini', status='old')
  read(unit, nml=header)
  !write(*, nml=header)
  
  ! Allocate memory for variables
  allocate(Files(nFiles))
  allocate(Wavelength(nWavelengths), x0(nPars), lb(nPars), ub(nPars), xopt(nPars), xtmp(nPars))
  !And initialize them
  Files=""
  ExpDataFName=""
  Wavelength=0.0_dp
  x0=0.0_dp
  lb=0.0_dp
  ub=0.0_dp
  r0=0.1_dp
  r1=1.0_dp
  xtol=1.0d-5
  idistrtype=0
  SAVEFDIST='fdist.dat'
  SAVEPOL='pol.dat'
  ! read data from file
  read(unit, nml=config)
  !write(*, nml=config)
  close(unit)
  
  allocate(f_data%Wavelength(nWavelengths))
  f_data%Wavelength(1:nWavelengths)=Wavelength(1:nWavelengths)
  f_data%r0=r0
  f_data%r1=r1
  f_data%idistrtype=idistrtype
  allocate(f_data%thtx(nExps), f_data%polx(nExps, nWavelengths))
  
  !read data table
  call read2cols(ExpDataFName, nExps, nWavelengths, f_data%thtx, f_data%polx)
  print '(A80)', 'EXPERIMENTAL DATA'
  print '(A80)', '================='
  call print_expdata(nExps, nWavelengths, f_data%thtx, f_data%polx)
  print *
  print '(A80)', 'STARTING APPROXIMATION OF A FUNCTION'
  print '(A80)', '===================================='
  print *
  
  opt=0
  call nlo_create(opt, OPTALG, nPars)
  call nlo_set_lower_bounds(ires, opt, lb)
  call nlo_set_upper_bounds(ires, opt, ub)
  call nlo_set_xtol_rel(ires, opt, xtol)
  call nlo_set_ftol_rel(ires, opt, ftol)
  
  
  optimval = 10000_dp
  ioptimval = 1
  
  print '(A3,".", 2X, "(",2A7,")", 2X, A7, A3)','No', 'Re(M)', 'Im(M)', 'RMSE', 'E'
  print *, '================================='
  DO I=1, nFiles
    call get_mueller_from_netcdf(Files(I), f_data%x_size, m1, m2, f_data%Mueller, f_data%Theta)
    
    
    call nlo_set_min_objective(ires, opt, cost_function, f_data)
    xtmp(:)=x0(:)
    call nlo_optimize(ires, opt, xtmp, minf)
    
    !if(ires.lt.0) stop 'Error in optimization'
    
    if(minf .lt. optimval) then
      optimval=minf
      ioptimval=I
      xopt(1:nPars)=xtmp(1:nPars)
    end if
    print '(I3,".",2X, "(",2F7.3,")", 2X, F7.3, I3)',I, M1, sqrt(minf), ires
    !print *, ires, x0, minf
    
  END DO
  
  call nlo_destroy(opt)
  print *
  print '(A80)', "OPTIMAL RESULT"
  print '(A80)', '=============='
  print *
  print '("RMSE = ", F7.2)', sqrt(optimval)
  write(fmtString,'("(A,", I2,"F10.6,", "A)")') nPars
  write (*,fmtString) "X=[",xopt,"]"
  write (*, '("Index of solution:", I3 )') ioptimval
  
  call get_mueller_from_netcdf(Files(ioptimval), f_data%x_size, m1, m2, f_data%Mueller, f_data%Theta)
  allocate(Pol(size(f_data%Theta, 1)))
  
  allocate(xx(1000), F(1000))
  DO I=1, 1000
    xx(I) = 1.0_dp+I*60.0_dp/1000
  END DO
  
  k_v=Kv(f_data%Wavelength(1))
  call Distribution(F, idistrtype, xx, k_v, r0, r1, xopt)

  open(newunit=unit, FILE=SAVEFDIST, status='replace')
  
  DO I=1, ubound(F, dim=1)
    write(unit, '(F10.4, E15.5)') xx(I)/k_v, F(I)*xx(I)**3
  END DO
  close(unit)
  
  deallocate(F, xx)
  
  allocate(F(size(f_data%X_size, dim=1)))
  
  call get_polarization_and_f(nPars, xopt, F, Pol, f_data)
  

  open(newunit=unit, FILE=SAVEPOL, status='replace')
  DO I=1, ubound(Pol, dim=1)
    write(unit, '(2F10.4)') f_data%Theta(I), Pol(I)
  END DO
  close(unit)
  
  
  !deallocate memory
  deallocate(F, Pol)
  deallocate(Files)
  deallocate(Wavelength, x0, lb, ub, xopt, xtmp)
  deallocate(f_data%Wavelength)
  deallocate(f_data%X_size)
  deallocate(f_data%Mueller)
  deallocate(f_data%Theta)
  deallocate(f_data%thtx, f_data%polx)
end program main
  

