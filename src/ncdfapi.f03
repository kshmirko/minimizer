module ncdfapi
	implicit none
	private
	integer, parameter	::	MAXDIM = 4
	public get_data_from_netcdf, get_mueller_from_netcdf
contains
	
	subroutine get_data_from_netcdf(FNAME, X, scale, m1, m2, Cext, Cabs, Csca, Cbk, Mueller, Theta)
! 		Запрос данных из netCDF файла
! 		Входные параметры:
!		------------------
!		character(len=*)	:: FNAME
!
!		Выходные параметры:
!		-------------------
!		complex*16	::	M1, M2
!		double precision	::	Cext(N), Csca(N), Cabs(N), Cbk(N) - сечения
!				экстинкции, рассеяния, погложения и обратного рассеяния
!		double precision	::	Mueller(K,K,M,N)	- матрицы мюллера
!		double precision	::	Theta(M) - углы рассеяния 
!   integer(kind=1)   ::  scale = Vc/Vo
!   Все параметры-массивы являются allocatable

		use netcdf
    use types, only : dp
		implicit none
		character(len=*), intent(in)	::	FNAME
    complex*16, intent(out)	::	m1, m2
		real(dp), allocatable, intent(out)	::	Theta(:), X(:), Cext(:), Cabs(:), Csca(:), Cbk(:)
		real(dp), allocatable, intent(out)	::	Mueller(:,:,:,:)
		integer(kind=1),intent(out)	::	scale
    
    integer N, M, K
		
		real(dp)  ::	m1re, m1im, m2re, m2im
		integer	::	ncid
		integer	::	cn_dimid, ck_dimid, cm_dimid ! иентификаторы размерностей
		integer	::	xsize_id, m1re_id, m1im_id, m2re_id, m2im_id, scale_id
		integer	::	cext_id, csca_id, cabs_id, cbk_id, mueller_id, theta_id
		integer	::	dimids(MAXDIM)
		integer	::	status, I
		
		call check( nf90_open(FNAME, NF90_NOWRITE, ncid) )
		
		call check( nf90_inq_dimid(ncid, "N", cn_dimid))
		call check( nf90_inquire_dimension(ncid, cn_dimid, len = N))
		
		call check( nf90_inq_dimid(ncid, "K", ck_dimid) )
		call check( nf90_inquire_dimension(ncid, ck_dimid, len = K) )
		
		call check( nf90_inq_dimid(ncid, "M", cm_dimid) )
		call check( nf90_inquire_dimension(ncid, cm_dimid, len = M) )
		
		allocate(Theta(M))
    allocate(X(N))
    allocate(Cext(N))
    allocate(Csca(N))
    allocate(Cabs(N))
    allocate(Cbk(N))
    allocate(Mueller(K,K,M,N))
		
		dimids = (/ck_dimid, ck_dimid, cm_dimid, cn_dimid/)
		
		call check( nf90_inq_varid(ncid, "X", xsize_id) )
		call check( nf90_get_var(ncid, xsize_id, X) )
		
		call check( nf90_inq_varid(ncid, "m1re", m1re_id) )
		call check( nf90_get_var(ncid, m1re_id, m1re) )
		
		call check( nf90_inq_varid(ncid, "m1im", m1im_id) )
		call check( nf90_get_var(ncid, m1im_id, m1im) )
		
		call check( nf90_inq_varid(ncid, "m2re", m2re_id) )
		call check( nf90_get_var(ncid, m2re_id, m2re) )
		
		call check( nf90_inq_varid(ncid, "m2im", m2im_id) )
		call check( nf90_get_var(ncid, m2im_id, m2im) )
		
		call check( nf90_inq_varid(ncid, "scale", scale_id) )
		call check( nf90_get_var(ncid, scale_id, scale) )
		
		
		m1 = DCMPLX(m1re, m1im)
		m2 = DCMPLX(m2re, m2im)
		
		
		call check( nf90_inq_varid(ncid, "Cext", cext_id) )
		call check( nf90_get_var(ncid, cext_id, Cext) )
		
		call check( nf90_inq_varid(ncid, "Csca", csca_id) )
		call check( nf90_get_var(ncid, csca_id, Csca) )
		
		call check( nf90_inq_varid(ncid, "Cabs", cabs_id) )
		call check( nf90_get_var(ncid, cabs_id, Cabs) )
		
		call check( nf90_inq_varid(ncid, "Mueller", mueller_id) )
		call check( nf90_get_var(ncid, mueller_id, Mueller) )
		
		call check( nf90_inq_varid(ncid, "Theta", theta_id) )
		call check( nf90_get_var(ncid, theta_id, Theta) )
		
		Cbk(:) = Mueller(1,1,181,:)
		
		call check( nf90_close(ncid) )
	end subroutine get_data_from_netcdf
  
	subroutine get_mueller_from_netcdf(FNAME, X, m1, m2, Mueller, Theta)
! 		Запрос данных из netCDF файла
! 		Входные параметры:
!		------------------
!		character(len=*)	:: FNAME
!
!		Выходные параметры:
!		-------------------
!		complex*16	::	M1, M2
!		double precision	::	Cext(N), Csca(N), Cabs(N), Cbk(N) - сечения
!				экстинкции, рассеяния, погложения и обратного рассеяния
!		double precision	::	Mueller(K,K,M,N)	- матрицы мюллера
!		double precision	::	Theta(M) - углы рассеяния 
!   integer(kind=1)   ::  scale = Vc/Vo
!   Все параметры-массивы являются allocatable

		use netcdf
    use types, only : dp
		implicit none
		character(len=*), intent(in)	::	FNAME
    complex*16, intent(out)	::	m1, m2
		real(dp), allocatable, intent(out)	::	Theta(:), X(:)
		real(dp), allocatable, intent(out)	::	Mueller(:,:,:,:)
		
    
    integer N, M, K
		
		real(dp)  ::	m1re, m1im, m2re, m2im
		integer	::	ncid
		integer	::	cn_dimid, ck_dimid, cm_dimid ! иентификаторы размерностей
		integer	::	xsize_id, m1re_id, m1im_id, m2re_id, m2im_id
		integer	::	mueller_id, theta_id
		integer	::	dimids(MAXDIM)
		integer	::	status, I
		
		call check( nf90_open(FNAME, NF90_NOWRITE, ncid) )
		
		call check( nf90_inq_dimid(ncid, "N", cn_dimid))
		call check( nf90_inquire_dimension(ncid, cn_dimid, len = N))
		
		call check( nf90_inq_dimid(ncid, "K", ck_dimid) )
		call check( nf90_inquire_dimension(ncid, ck_dimid, len = K) )
		
		call check( nf90_inq_dimid(ncid, "M", cm_dimid) )
		call check( nf90_inquire_dimension(ncid, cm_dimid, len = M) )
		
		allocate(Theta(M))
    allocate(X(N))
    allocate(Mueller(K,K,M,N))
		
		dimids = (/ck_dimid, ck_dimid, cm_dimid, cn_dimid/)
		
		call check( nf90_inq_varid(ncid, "X", xsize_id) )
		call check( nf90_get_var(ncid, xsize_id, X) )
		
		call check( nf90_inq_varid(ncid, "m1re", m1re_id) )
		call check( nf90_get_var(ncid, m1re_id, m1re) )
		
		call check( nf90_inq_varid(ncid, "m1im", m1im_id) )
		call check( nf90_get_var(ncid, m1im_id, m1im) )
		
		call check( nf90_inq_varid(ncid, "m2re", m2re_id) )
		call check( nf90_get_var(ncid, m2re_id, m2re) )
		
		call check( nf90_inq_varid(ncid, "m2im", m2im_id) )
		call check( nf90_get_var(ncid, m2im_id, m2im) )
		
		m1 = DCMPLX(m1re, m1im)
		m2 = DCMPLX(m2re, m2im)
		
		call check( nf90_inq_varid(ncid, "Mueller", mueller_id) )
		call check( nf90_get_var(ncid, mueller_id, Mueller) )
		
		call check( nf90_inq_varid(ncid, "Theta", theta_id) )
		call check( nf90_get_var(ncid, theta_id, Theta) )
		
		call check( nf90_close(ncid) )
	end subroutine get_mueller_from_netcdf
	
	subroutine check(status)
		use netcdf
		implicit none
		integer, intent ( in) :: status
		if(status /= nf90_noerr) then 
			print *, trim(nf90_strerror(status))
			stop "Stopped"
		end if
	end subroutine check
	
end module ncdfapi
