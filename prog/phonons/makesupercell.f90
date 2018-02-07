program makesupercell
  implicit none
  integer, parameter :: dp = 8
  character(60) :: arg, atomstr
  character(2) :: cs 
  integer :: i, natm, l1, l2, l3, n1, n2, n3, funit, tmp 
  real(dp), dimension(:), allocatable :: x, y, z
  integer, dimension(:), allocatable :: atmind
  real(dp) :: xl, yl, zl, v(3,3), or(3), p(3)

  funit = 5

  if (iargc() < 2) then
    write(*,*) "Usage: gen2gen -s n1 n2 n3 < strin.gen > strout.gen"
    write(*,*) "     Code to generate supercells for phonon calculation "
    write(*,*) " -s  Supercell option, repeat n1 n2 n3 "
    write(*,*) "     The cells are repeated -n1*a1..+n1*a1, -n2*a2..+n2*a2, ... "
    write(*,*) "     use n1=n2=n3 for bulk systems and n1=n2 n3=0 in 2D"
    stop 
  end if

  n1 = 1; n2 = 1; n3 = 1;

  do i = 1, iargc()
    call getarg(i, arg)
    if (trim(arg) == "-s") then
      call getarg(i+1, arg)
      read(arg,*) n1 
      call getarg(i+2, arg)
      read(arg,*) n2 
      call getarg(i+3, arg)
      read(arg,*) n3 
    end if
  end do   

  if (n1 /= n2) then
    write(*,*) "WARNING: n1 != n2 so probably there will be periodicity problems..."
  end if  

  ! read first two lines of gen file
  read(funit,*) natm, cs
  read(funit,*) atomstr 
  
  allocate(x(natm))
  allocate(y(natm))
  allocate(z(natm))
  allocate(atmind(natm))

  ! read all gen file
  do i = 1, natm
    read(funit,*) tmp, atmind(i), x(i), y(i), z(i)
  end do
  read(funit,*) or
  read(funit,*) v(1,1),v(1,2),v(1,3)
  read(funit,*) v(2,1),v(2,2),v(2,3)
  read(funit,*) v(3,1),v(3,2),v(3,3)

  ! write new gen file
  write(*,*) natm*(2*n1+1)*(2*n2+1)*(2*n3+1), cs
  write(*,*) atomstr
  
  do l3 = -n3, n3
    do l2 = -n2, n2
      do l1 = -n1, n1
        if (l1 == 0 .and. l2 == 0 .and. l3 == 0) continue

        do i = 1, natm
          xl = x(i) + v(1,1)*l1 + v(2,1)*l2 + v(3,1)*l3 
          yl = y(i) + v(1,2)*l1 + v(2,2)*l2 + v(3,2)*l3 
          zl = z(i) + v(1,3)*l1 + v(2,3)*l2 + v(3,3)*l3 
          write(*,*) i, atmind(i), xl, yl, zl
        end do
      end do
    end do
  end do  

  write(*,*) or
  write(*,*) (2*n1+1)*v(1,1), (2*n1+1)*v(1,2), (2*n1+1)*v(1,3)
  write(*,*) (2*n2+1)*v(2,1), (2*n2+1)*v(2,2), (2*n2+1)*v(2,3)
  write(*,*) (2*n3+1)*v(3,1), (2*n3+1)*v(3,2), (2*n3+1)*v(3,3)

end program makesupercell
