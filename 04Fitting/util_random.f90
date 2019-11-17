!
! -------------------------------------------
! Module of random number generaion utilities
! -------------------------------------------
!   module mt19937 is requied
!
!----------------------------------------------------------------------
module util_random
  implicit none
  private

  !============= PUBLIC ROUTINES =============
  public :: init_grand      ! initialize random routines
  public :: grand           ! uniform [0,1]
  public :: grand2          ! uniform [0,1)
  public :: grand3          ! uniform (0,1)
  public :: ngrand          ! standard normal
  public :: vgrand          ! vector of uniform [0,1]
  public :: vngrand         ! vector of standard normal
  public :: init_mvngrand   ! initialize MVN (covariance ver.)
  public :: init_mvngrandr  ! initialize MVN (cross-cor. ver.)
  public :: init_mvngrandr2 ! initialize MVN (cross-cor. with zero ver.)
  public :: mvngrand        ! multi-variate normal (covariance ver.)
  public :: mvngrand2       ! multi-variate normal (cross-cor. with zero)

  !============= PRIVATE VARIABLES =============
  integer, parameter :: max_nparam = 50
  real    :: m(1:max_nparam)
  real    :: lmat(1:max_nparam, 1:max_nparam)
  integer :: n2
  integer :: i2(1:max_nparam)
  real    :: m2(1:max_nparam)
contains

!======================================================================
!========================== PUBLIC ROUTINES ===========================

!----------------------------------------------------------------------
! Initialize random number routined
!----------------------------------------------------------------------
subroutine init_grand
  use mt19937, only : init_genrand
  call init_genrand(gen_seed())
end subroutine init_grand

!----------------------------------------------------------------------
! Uniform random number [0,1]
!----------------------------------------------------------------------
function grand() result(r)
  use mt19937, only : genrand_real1
  real :: r
  r = real(genrand_real1())
end function grand

!----------------------------------------------------------------------
! Uniform random number [0,1)
!----------------------------------------------------------------------
function grand2() result(r)
  use mt19937, only : genrand_real2
  real :: r
  r = real(genrand_real2())
end function grand2

!----------------------------------------------------------------------
! Uniform random number (0,1)
!----------------------------------------------------------------------
function grand3() result(r)
  use mt19937, only : genrand_real3
  real :: r
  r = real(genrand_real3())
end function grand3

!----------------------------------------------------------------------
! Normal random number
!----------------------------------------------------------------------
function ngrand() result(r)
  real          :: r
  real          :: v1, v2, p, s
  real,    save :: next
  logical, save :: have_next
  data next      / 0.0 /
  data have_next / .false. /
  if(have_next) then
    r = next
    have_next = .false.
  else
    s = -1.0
    do while((s >= 1.0) .or. (s <= 0.0))
      v1 = 2.0 * grand() - 1.0
      v2 = 2.0 * grand() - 1.0
      s = v1 ** 2 + v2 ** 2
    end do
    p = sqrt(-2.0 * log(s) / s)
    r = v1 * p
    next = v2 * p
    have_next = .true.
  end if
end function ngrand

!----------------------------------------------------------------------
! Vector of uniform random number [0,1]
!----------------------------------------------------------------------
function vgrand(n) result(r)
  use mt19937, only : genrand_real1
  integer, intent(in) :: n           ! number of element
  real                :: r(1:n)
  integer :: i
  do i = 1, n
    r(i) = real(genrand_real1())
  end do
end function vgrand

!----------------------------------------------------------------------
! Vector of normal random numbers
!----------------------------------------------------------------------
function vngrand(n) result(r)
  integer, intent(in) :: n           ! number of element
  real                :: r(1:n)
  integer :: i
  do i = 1, n
    r(i) = ngrand()
  end do
end function vngrand

!----------------------------------------------------------------------
! Multivariate normal random number
! (covariance varsion)
!----------------------------------------------------------------------
function mvngrand(n) result(x)
  integer, intent(in) :: n           ! number of dimension
  real                :: x(1:n)
  real    :: vnrand(1:n)
  integer :: i
  vnrand(1:n) = vngrand(n)
  do i = 1, n
    x(i) = dot_product(lmat(1:n, i), vnrand(1:n)) + m(i)
  end do
end function mvngrand

!----------------------------------------------------------------------
! Multivariate normal random number
! (cross-correlation version with zero value)
!----------------------------------------------------------------------
function mvngrand2(n) result(x)
  integer, intent(in) :: n           ! number of dimension
  real                :: x(1:n)
  real    :: xx(1:n2)
  integer :: i
  xx(1:n2) = mvngrand(n2)
  do i = 1, n
    if(i2(i) > 0) then
      x(i) = xx(i2(i))
    else
      x(i) = m2(i)
    end if
  end do
end function mvngrand2

!----------------------------------------------------------------------
! Initialize multivariate normal random number routines
! (covariance varsion)
!----------------------------------------------------------------------
subroutine init_mvngrand(n, mean, varcovar)
  integer, intent(in) :: n              ! number of dimension
  real,    intent(in) :: mean(:)        ! mean value vector
  real,    intent(in) :: varcovar(:, :) ! variance-covariance matrix
  m(1:n) = mean(1:n)
  call cholesky(n, varcovar, lmat)
end subroutine init_mvngrand

!----------------------------------------------------------------------
! Initialize multivariate normal random number routines
! (cross-correlation version)
!----------------------------------------------------------------------
subroutine init_mvngrandr(n, mean, var, cormat)
  integer, intent(in) :: n              ! number of dimension
  real,    intent(in) :: mean(:)        ! mean value vector
  real,    intent(in) :: var(:)         ! variance vector
  real,    intent(in) :: cormat(:, :)   ! cross-correlation matrix
  real    :: varcovar(1:n, 1:n)
  integer :: i, j
  do i = 1, n
    do j = 1, n
      varcovar(i, j) = cormat(i, j) * sqrt(var(i)) * sqrt(var(j))
    end do
  end do
  call init_mvngrand(n, mean, varcovar)
end subroutine init_mvngrandr

!----------------------------------------------------------------------
! Initialize multivariate normal random number routines
! (cross-correlation version with zero value)
!----------------------------------------------------------------------
subroutine init_mvngrandr2(n, mean, var, cormat)
  integer, intent(in) :: n              ! number of dimension
  real,    intent(in) :: mean(:)        ! mean value vector
  real,    intent(in) :: var(:)         ! variance vector
  real,    intent(in) :: cormat(:, :)   ! cross-correlation matrix
  real    :: m(1:n)
  real    :: v(1:n)
  real    :: c(1:n, 1:n)
  integer :: i, j, ii, jj
  m2(1:n) = mean(1:n)
  i2(1:n) = 0
  ii = 0
  do i = 1, n
    if(v(i) /= 0.0) then
      ii = ii + 1
      m(ii) = mean(i)
      v(ii) = var(i)
      i2(i) = ii
      jj = 0
      do j = 1, n
        if(v(j) /= 0.0) then
          jj = jj + 1
          c(ii, jj) = cormat(i, j)
        end if
      end do
    else
      i2(i) = 0
    end if
  end do
  n2 = ii
  call init_mvngrandr(ii, m, v, c)
end subroutine init_mvngrandr2

!======================================================================
!========================== PRIVATE ROUTINES ==========================
!
!----------------------------------------------------------------------
! Genetrate seed
!----------------------------------------------------------------------
function gen_seed() result(r)
  character(len=30) :: s = '                              '
  integer :: i, r
  call fdate(s)
  r = 0
  do i=1, 22
    r = r + iachar(s(i:i)) * i
  end do
end function gen_seed

!----------------------------------------------------------------------
! Cholesky decomposition
!----------------------------------------------------------------------
subroutine cholesky(n, A, L)
  integer, intent(in)  :: n
  real,    intent(in)  :: A(:, :)
  real,    intent(out) :: L(:, :)
  real,allocatable :: aa(:), uu(:)
  integer :: nn, i, j, k, ii, ierr
  nn = n * (n + 1) / 2
  allocate(aa(nn))
  allocate(uu(nn))
  ii=1
  do j=1,n
    do i=1,j
      aa(ii)=A(i,j)
      ii=ii+1
    end do
  end do
  call chol(aa, n, nn, uu, k, ierr)
  if(ierr > 0) then
    print *, 'ierr =',ierr
    stop
  end if
  ii=1
  do j=1,n
    do i=1,j
      L(i,j) = uu(ii)
      ii=ii+1
    end do
  end do
  do j=1,n
    do i=j+1,n
      L(i,j) = 0.0
    end do
  end do
  deallocate(aa)
  deallocate(uu)

end subroutine cholesky

!=======================================================================
! This file contains AS6 and the enhanced version ASR44.   See AS7 also.
! 
! 
SUBROUTINE CHOL(A,N,NN,U,NULLTY,IFAULT)
! 
!       Algorithm AS6, Applied Statistics, vol.17, (1968)
! 
!       Given a symmetric matrix order n as lower triangle in a( )
!       calculates an upper triangle, u( ), such that uprime * u = a.
!       a must be positive semi-definite.  eta is set to multiplying
!       factor determining effective zero for pivot.
! 
!       arguments:-
!       a()     = input, a +ve definite matrix stored in lower-triangula
!                 form.
!       n       = input, the order of a
!       nn      = input, the size of the a and u arrays      n*(n+1)/2
!       u()     = output, a lower triangular matrix such that u*u' = a.
!                 a & u may occupy the same locations.
!       nullty  = output, the rank deficiency of a.
!       ifault  = output, error indicator
!                       = 1 if n < 1
!                       = 2 if a is not +ve semi-definite
!                       = 3 if nn < n*(n+1)/2
!                       = 0 otherwise
! 
!***********************************************************************
! 
      INTEGER N,NN,NULLTY,IFAULT
      REAL A(1:NN),U(1:NN),ETA,ETA2,X,W,ZERO
      INTEGER I,J,K,L,M,II,KK,ICOL,IROW

!       The value of eta will depend on the word-length of the
!       computer being used.  See introductory text.

      DATA ETA,ZERO/1.E-6,0.0/

      IFAULT=1
      IF (N.LE.0) RETURN
      IFAULT=3
      IF (NN.LT.N*(N+1)/2) RETURN
      IFAULT=2
      NULLTY=0
      J=1
      K=0
      ETA2=ETA*ETA
      II=0
      W=0.0

!       Factorize column by column, icol = column no.

      DO 80 ICOL=1,N
        II=II+ICOL
        X=ETA2*A(II)
        L=0
        KK=0

!       IROW = row number within column ICOL

        DO 40 IROW=1,ICOL
          KK=KK+IROW
          K=K+1
          W=A(K)
          M=J
          DO 10 I=1,IROW
            L=L+1
            IF (I.EQ.IROW) GO TO 20
            W=W-U(L)*U(M)
            M=M+1
 10       CONTINUE
 20       IF (IROW.EQ.ICOL) GO TO 50
          IF (U(L).EQ.ZERO) GO TO 30
          U(K)=W/U(L)
          GO TO 40
 30       IF (W*W.GT.ABS(X*A(KK))) RETURN
          U(K)=ZERO
 40     CONTINUE
 50     IF (ABS(W).LE.ABS(ETA*A(K))) GO TO 60
        IF (W.LT.ZERO) RETURN
        U(K)=SQRT(W)
        GO TO 70
 60     U(K)=ZERO
        NULLTY=NULLTY+1
 70     J=J+ICOL
 80   CONTINUE
      IFAULT=0
END SUBROUTINE CHOL


end module util_random
