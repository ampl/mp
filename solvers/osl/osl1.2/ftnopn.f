	subroutine ftnopn(n, name, create, ierr)
	integer n, ierr
	character*(*) name
	logical create

	open(n,file=name,status='OLD',iostat=ierr)
	if (ierr .ne. 0 .and. create) then
		open(n,file=name,status='NEW',iostat=ierr)
		endif
	if (ierr .eq. 0) rewind n
	end
