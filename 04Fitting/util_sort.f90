!
! ---------------------------
! Module of Sorting Utilities
! ---------------------------
!
!----------------------------------------------------------------------
module util_sort
  implicit none
  private

  !============= PUBLIC ROUTINES =============
  public :: qsort

contains

!======================================================================
!========================== PUBLIC ROUTINES ===========================

!----------------------------------------------------------------------
! Quick Sort
!----------------------------------------------------------------------
subroutine qsort(a, ip, dir)
  real,    intent(in)  :: a(:)   ! reference value (nondestructive)
  integer, intent(out) :: ip(:)  ! sorted index (initialized)
  integer, intent(in)  :: dir    ! sorting direction (+:ao, -:do)
  integer :: i
  real    :: f
  do i = 1, size(a)
    ip(i) = i
  end do
  if(dir >= 0) then
    f = 1.0
  else
    f = -1.0
  end if
  call QsortC(a, ip, f)
end subroutine qsort

!======================================================================
!========================== PRIVATE ROUTINES ==========================
recursive subroutine QsortC(A, ip, dir)
  real,    intent(in)    :: A(:)
  integer, intent(inout) :: ip(:)
  real,    intent(in)  :: dir
  integer :: iq

  if(size(ip) > 1) then
     call Partition(A, ip, iq, dir)
     call QsortC(A, ip(:iq-1), dir)
     call QsortC(A, ip(iq:), dir)
  endif
end subroutine QsortC

subroutine Partition(A, ip, marker, dir)
  real,    intent(in)    :: A(:)
  integer, intent(inout) :: ip(:)
  integer, intent(out)   :: marker
  real,    intent(in)    :: dir
  integer :: i, j
  integer :: temp
  real    :: x      ! pivot point

  x = A(ip(1))
  i= 0
  j= size(ip) + 1

  do
     j = j - 1
     do
        if(A(ip(j)) * dir <= x * dir) exit
        j = j - 1
     end do
     i = i + 1
     do
        if(A(ip(i)) * dir >= x * dir) exit
        i = i+1
     end do
     if(i < j) then
        temp = ip(i)
        ip(i) = ip(j)
        ip(j) = temp
     else if(i == j) then
        marker = i + 1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module util_sort
