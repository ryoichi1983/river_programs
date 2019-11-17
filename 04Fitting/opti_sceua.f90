!
! --------------------------------------
! Module of SCE-UA optimization routines
! --------------------------------------
!    module util_random is requied
!    module util_sort   is requied
!    module util_qmc    is requied
!
!----------------------------------------------------------------------
module opti_sceua
  implicit none
  private

  !============= PUBLIC ROUTINES =============
  public :: do_sceua     ! simple routine with standard control parameters
  public :: do_sceuaf    ! simple routine with configuration file
  public :: sceua        ! primitive routine
  public :: sce_cbdummy  ! dummy callback function

  !============= PUBLIC PARAMETERS =============
  !----- control flags -----
  integer, public, parameter :: sce_maximize = 0  ! maximize object function   Åö1ÅÀ0Ç…ïœçX(âPíJ)
  integer, public, parameter :: sce_minimize = 1  ! minimize object function   Åö0ÅÀ1Ç…ïœçX(âPíJ)
  integer, public, parameter :: f_quiet      = 1  ! restrain messages
  integer, public, parameter :: f_initp      = 2  ! use initial parameters
  integer, public, parameter :: f_widemutate = 4  ! wide range mutation (recommended)
  integer, public, parameter :: f_qmcinit    = 8  ! initialize by Quasi-Monte Carlo

  !============= PRIVATE VARIABLES =============
  integer, parameter :: max_nparam = 110      ! maximum number of parameter
  real               :: scale_a(max_nparam)   ! parameter scaling factor A
  real               :: scale_b(max_nparam)   ! parameter scaling factor B
  real, parameter    :: constrained = -1.E37  ! evaluation value of constrained particle
  real, parameter    :: void = -1.E38         ! evaluation value of undefined particle

  type type_p
    real :: x(max_nparam)          ! coordinates of point
    real :: ev                     ! evaluation value
    real :: p                      ! selection probability
  end type type_p

contains


!======================================================================
! SCE-UA Simplified Wrapper Routine (with standard control parameters)
!======================================================================
subroutine do_sceua(n, xmin, xmax, xinit, nc, it, ir, maximize, evfunc, x, r)
  integer, intent(in)  :: n             ! number of parameter
  real,    intent(in)  :: xmin(1:n)     ! lower bounds
  real,    intent(in)  :: xmax(1:n)     ! upper bounds
  real,    intent(in)  :: xinit(1:n)    ! initial parameters
  integer, intent(in)  :: nc            ! number of complex
  integer, intent(in)  :: it            ! maximum iteration count (<=0:unlimitted)
  integer, intent(in)  :: ir            ! interval of message printing (<=0:nomessage)
  integer, intent(in)  :: maximize      ! maximize-minimize flag (>0:max, <=0:min)
  real,    intent(out) :: x(1:n)        ! optimum parameter set
  real,    intent(out) :: r             ! optimum evaluation value
  interface
    !-------- evaluation function -------
    integer function evfunc(x, r)  ! return number of constrained condition
      real, intent(in)  :: x(:)    ! parameter set
      real, intent(out) :: r       ! evaluation value
    end function evfunc
  end interface

  !====== Local variables ======
  integer :: nm, ne, alpha, beta, ie, flag

  !====== Setup control parameters =====
  nm = 2 * n + 1           ! member in a complex (default:2*n+1)
  ne = n + 1               ! evolving member in a complex (default:n+1)
  alpha = 1                ! number of alpha iteration (default:1)
  beta = 2 * n + 1         ! number of beta iteration (default:2*n+1)
  ie = 0                   ! maximum evaluation count (<=0:unlimitted)

  !====== Setup control flags ======
  flag = 0                     ! control flag
!! flag = flag + f_quiet           ! restrain messages
   flag = flag + f_initp           ! use initial parameters
   flag = flag + f_widemutate      ! wide range mutation (recommended)
!! flag = flag + f_qmcinit         ! initialize by Quasi-Monte Carlo

  !====== Optimize ======
  call sceua(n, xmin, xmax, xinit, nc, nm, ne, alpha, beta, &
&            it, ie, ir, 0, flag, maximize, evfunc, sce_cbdummy, x, r)

end subroutine do_sceua

!======================================================================
! SCE-UA Simplified Wrapper Routine (with configuration file)
!======================================================================
subroutine do_sceuaf(n, xmin, xmax, xinit, maximize, evfunc, pfname, x, r)
  integer, intent(in)  :: n             ! number of parameter
  real,    intent(in)  :: xmin(1:n)     ! lower bounds
  real,    intent(in)  :: xmax(1:n)     ! upper bounds
  real,    intent(in)  :: xinit(1:n)    ! initial parameters
  integer, intent(in)  :: maximize      ! maximize-minimize flag (>0:max, <=0:min)
  character(len=*)     :: pfname        ! configuration file name
  real,    intent(out) :: x(1:n)        ! optimum parameter set
  real,    intent(out) :: r             ! optimum evaluation value
  interface
    !-------- evaluation function -------
    integer function evfunc(x, r)  ! return number of constrained condition
      real, intent(in)  :: x(:)    ! parameter set
      real, intent(out) :: r       ! evaluation value
    end function evfunc
  end interface

  integer :: i, nn
  integer :: nc, nm, ne, alpha, beta, it, ie, ir, flag
  integer :: flag_quiet, flag_initp, flag_widemutate, flag_qmcinit
  namelist /opti_sceua/ nc, nm, ne, alpha, beta,it,ie,ir, &
&                       flag_quiet, flag_initp, flag_widemutate, flag_qmcinit

  open(1, file=pfname, status='old')
  read(1, opti_sceua)
  close(1)

  flag = flag_quiet
  flag = flag + flag_initp      * 2**1
  flag = flag + flag_widemutate * 2**2
  flag = flag + flag_qmcinit    * 2**3

  nn = 0
  do i = 1, n
    if(xmin(i) /= xmax(i)) then
      nn = nn + 1
    end if
  end do
  if(nm    == 0) nm = 2 * n + 1
  if(ne    == 0) ne = n + 1
  if(alpha == 0) alpha = 1
  if(beta  == 0) beta = 2 * n + 1

  call sceua(n, xmin, xmax, xinit, nc, nm, ne, alpha, beta, &
&                        it, ie, ir, 0, flag, maximize, evfunc, sce_cbdummy, x, r)

end subroutine do_sceuaf

!======================================================================
! SCE-UA (primitive routine)
!======================================================================
subroutine sceua(n, xmin, xmax, xinit, nc, nm, ne, alpha, beta, &
&                 it, ie, ir, ic, flag, maximize, evfunc, cbfunc, &
&                                                    xopti, ropti)
  use util_random, only : init_grand
  integer, intent(in)  :: n             ! number of parameter
  real,    intent(in)  :: xmin(1:n)     ! lower bounds
  real,    intent(in)  :: xmax(1:n)     ! upper bounds
  real,    intent(in)  :: xinit(1:n)    ! initial parameters
  integer, intent(in)  :: nc            ! number of complex
  integer, intent(in)  :: nm            ! member in a complex (default:2*n+1)
  integer, intent(in)  :: ne            ! evolving member in a complex (default:n+1)
  integer, intent(in)  :: alpha         ! number of alpha iteration (default:1)
  integer, intent(in)  :: beta          ! number of beta iteration (default:2*n+1)
  integer, intent(in)  :: it            ! maximum iteration count (<=0:unlimitted)
  integer, intent(in)  :: ie            ! maximum evaluation count (<=0:unlimitted)
  integer, intent(in)  :: ir            ! interval of message printing (<=0:nomessage)
  integer, intent(in)  :: ic            ! interval of callback function (<=0:no call)
  integer, intent(in)  :: flag          ! control flag
  integer, intent(in)  :: maximize      ! maximize-minimize flag (>0:max, <=0:min)
  real,    intent(out) :: xopti(1:n)    ! optimum parameter set
  real,    intent(out) :: ropti         ! optimum evaluation value
  interface
    !-------- evaluation function -------
    integer function evfunc(x, r)  ! return number of constrained condition
      real, intent(in)  :: x(:)    ! parameter set
      real, intent(out) :: r       ! evaluation value
    end function evfunc
    !-------- callback function -------
    integer function cbfunc(it, iev, ev, x) ! return continuation flag(>=0:cont.,<0:quit)
      integer, intent(in) :: it    ! current iteration counter
      integer, intent(in) :: iev   ! current evaluation counter
      real,    intent(in) :: ev    ! current evaluation value
      real,    intent(in) :: x(:)  ! current parameter value
    end function cbfunc
  end interface

  integer                   :: np            ! total number of point
  type(type_p), allocatable :: p(:)          ! point
  integer,      allocatable :: id(:)         ! index of points
  integer,      allocatable :: idp(:,:)      ! index of complex devided points
  real                      :: gebest        ! best evaluation value
  real                      :: gxbest(1:n)   ! best position
  integer                   :: evcount       ! counter of evaluation
  integer                   :: sign          ! sign of evaluation value
  logical      :: exit_loop                  ! exit flag
  integer      :: i, im, ip, k

  !===== Print SCE-UA control parameters =====
  if(flag_off(flag, f_quiet)) then
    call print_param(n, nc, nm, ne, alpha, beta, it, ie, flag)
    print *, 'SCE-UA --- cycle / best / mean / constrained(%)/ evaluation ---'
  end if

  !===== Check arguments =====
  call check_param(n, xmin, xmax, xinit, nc, nm, ne, alpha, beta, it, flag)

  !===== Allocate working memory =====
  np = nc * nm                    ! number of particle
  allocate(p(np))
  allocate(id(np))
  allocate(idp(nm, nc))

  !===== Initialize scaling function =====
  call init_scaling(n, xmin, xmax)

  !===== Initialize random number generator =====
  call init_grand

  !===== Initialize sign of evaluation =====
  if(maximize > 0) then
    sign = 1
  else
    sign = -1
  end if

  !===== Initialize particles =====
  evcount = 0
  call init_particles(n, xinit, np, flag, p)

  !===== Evaluate all particles =====
  !$omp parallel do
  do i=1, np
    p(i)%ev = evaluate(n, p(i)%x, sign, evfunc, evcount)
  end do
  !$omp end parallel do


  !***********************************************************
  !*********************** Main Loop *************************
  i = 0
  do

    !===== Sort particles by evaluation value =====
    call psort(np, p, id)
    gebest = p(id(1))%ev
    gxbest(1:n) = p(id(1))%x(1:n)

    !===== Print intermediate conditions =====
    if(flag_off(flag, f_quiet) .and. (mod(i, ir) == 0)) then
      print *, 'SCE-UA ---', i, gebest * sign, evmean(nc, p) * sign, &
&                             infeasible(np, p), evcount
    end if

    !===== Check exit conditions =====
    if((it > 0 .and. i >= it) .or. (ie > 0 .and. evcount >= ie)) then
      exit_loop = .true.
    else
      exit_loop = .false.
    end if

    !===== Call callback function =====
    if((ic > 0) .and. (mod(i, ic) == 0 .or. exit_loop)) then
      if(cbfunc(i, evcount, gebest, unscaling(n, gxbest(1:n))) < 0) then
        exit_loop = .true.
      end if
    end if

    !===== Exit loop =====
    if(exit_loop) exit

    !===== Deivide points into complexes =====
    k = 1
    do im = 1, nm
      do ip = 1, nc
        idp(im, ip) = id(k)
        k = k + 1
      end do
    end do

    !===== Evolve complexes by CCE algorithm =====
    !$omp parallel do
    do ip = 1, nc
       call cce(n, nm, ne, alpha, beta, flag, sign, evfunc, idp(1:nm, ip), evcount, p)
    end do
    !$omp end parallel do

    i = i + 1
  end do
  !*********************** Main Loop *************************
  !***********************************************************


  !===== Unscaling parameters =====
  xopti = unscaling(n, p(id(1))%x)
  ropti = gebest * sign

  !===== Deallocate working memory =====
  deallocate(p)
  deallocate(id)
  deallocate(idp)

  !===== Print exit message =====
  if(flag_off(flag, f_quiet)) then
    print *, 'SCE-UA --- finish ---'
    print *
  end if

end subroutine sceua

!----------------------------------------------------------------------
! Evolve a complex by CCE algorithm
!----------------------------------------------------------------------
subroutine cce(n, nm, ne, alpha, beta, flag, sign, evfunc, id, evcount, p)
  use util_random, only : vgrand
  use util_sort,   only : qsort
  integer,      intent(in)    :: n          ! number of parameter
  integer,      intent(in)    :: nm         ! number of member in a complex
  integer,      intent(in)    :: ne         ! number of evolution member
  integer,      intent(in)    :: alpha      ! number of alpha iteration
  integer,      intent(in)    :: beta       ! number of beta iteration
  integer,      intent(in)    :: flag       ! control flag
  integer,      intent(in)    :: sign       ! sign of evaluation value
  integer,      intent(in)    :: id(1:nm)   ! index of members
  integer,      intent(inout) :: evcount    ! counter of evaluation
  type(type_p), intent(inout) :: p(:)       ! points
  interface
    integer function evfunc(x, r)
      real, intent(in)  :: x(:)
      real, intent(out) :: r
    end function evfunc
  end interface

  type(type_p) :: q(1:nm)      ! members
  real         :: rr(1:nm)     ! roulette wheel
  integer      :: idm(1:nm)    ! index of members
  real         :: vrand(1:nm)  ! uniform random number vector
  integer      :: ii

  !===== Prepare members =====
  q(1:nm) = p(id(1:nm))

  !***************************************************
  !******************* Beta loop *********************
  do ii = 1, beta

    call mkprob(nm, q)                     ! calculate selection probability
    vrand = vgrand(nm)                     ! make uniform random number vector
    rr(1:nm) = q(1:nm)%p * vrand(1:nm)     ! make roulette wheel (bitch!)
    call qsort(rr, idm, -1)                ! sort in descending order

    !===== Evolve superior ne members =====
    call evolution(n, nm, ne, alpha, flag, sign, evfunc, idm(1:ne), evcount, q)

  end do
  !******************* Beta loop *********************
  !***************************************************

  !===== Restore members to all points =====
  p(id(1:nm)) = q(1:nm)

end subroutine cce

!----------------------------------------------------------------------
! Evolve selected members by CCE algorithm
!----------------------------------------------------------------------
subroutine evolution(n, nm, ne, alpha, flag, sign, evfunc, id, evcount, p)
  integer,      intent(in)    :: n          ! number of parameter
  integer,      intent(in)    :: nm         ! number of all member
  integer,      intent(in)    :: ne         ! number of evolving member
  integer,      intent(in)    :: alpha      ! number of alpha iteration
  integer,      intent(in)    :: flag       ! control flag
  integer,      intent(in)    :: sign       ! sign of evaluation value
  integer,      intent(in)    :: id(1:ne)   ! index of evolving members
  integer,      intent(inout) :: evcount    ! counter of evaluation
  type(type_p), intent(inout) :: p(:)       ! all members
  interface
    integer function evfunc(x, r)
      real, intent(in)  :: x(:)
      real, intent(out) :: r
    end function evfunc
  end interface

  type(type_p) :: q(1:ne)      ! selected evolving members
  integer      :: idq(1:ne)    ! index of evolving member
  integer      :: iu           ! index to the worst member
  type(type_p) :: g            ! gravity center
  type(type_p) :: r            ! symmetric point of the worst member
  integer      :: ii

  !===== Prepare selected members =====
  q(1:ne) = p(id(1:ne))

  !***************************************************
  !******************* Alpha loop ********************
  do ii = 1, alpha

    call psort(ne, q, idq)                     ! sort by evaluation value
    g%x(1:n) = gcenter(n, ne - 1, q, idq)      ! G: gravity center of top ne-1 members
    iu = idq(ne)                               ! U: worst member
    r%x(1:n) = 2.0 * g%x(1:n) - q(iu)%x(1:n)   ! R: symmetric point of U about G

    !===== If R is out of range =====
    if(outrange(n, r%x)) then
      r%x(1:n) =  mutate(n, nm, p, flag)                ! mutation
    end if
    r%ev = evaluate(n, r%x, sign, evfunc, evcount)      ! evaluate R

    !===== If R is worse than U ====
    if(r%ev < q(iu)%ev) then
      r%x(1:n) = (r%x(1:n) + g%x(1:n)) / 2.0            ! R = (U + G) / 2
      r%ev = evaluate(n, r%x, sign, evfunc, evcount)    ! evaluate R
      if(r%ev < q(iu)%ev) then                          ! If R is worse than U
        r%x(1:n) =  mutate(n, nm, p, flag)                 ! mutation
        r%ev = evaluate(n, r%x, sign, evfunc, evcount)     ! evaluate R
      end if
    end if

    !===== Replace U with R =====
    q(iu)%x(1:n) = r%x(1:n)
    q(iu)%ev     = r%ev

  end do
  !******************* Alpha loop ********************
  !***************************************************

  !===== Restore selected members to all members =====
  p(id(1:ne)) = q(1:ne)

end subroutine evolution

!----------------------------------------------------------------------
! Mutate
!----------------------------------------------------------------------
function mutate(n, np, p, flag) result(x)
  use util_random, only : grand, vgrand
  integer,      intent(in) :: n
  integer,      intent(in) :: np
  type(type_p), intent(in) :: p(:)
  integer,      intent(in) :: flag
  real                     :: x(1:n)
  real    :: xmax, xmin, xx
  integer :: i, j
  if(flag_on(flag, f_widemutate)) then
    do i = 1, n
      xmax = -1.e10
      xmin = 1.e10
      do j = 1, np
        xx = p(j)%x(i)
        xmax = max(xx, xmax)
        xmin = min(xx, xmin)
      end do
      x(i) = grand() * (xmax - xmin) + xmin
    end do
  else
    x(1:n) = vgrand(n)
  end if
end function mutate

!----------------------------------------------------------------------
! Check Range of Coordinates
!----------------------------------------------------------------------
function outrange(n, x) result(r)
  integer, intent(in) :: n
  real,    intent(in) :: x(1:n)
  logical             :: r
  integer :: i
  r = .false.
  do i = 1, n
    if((x(i) > 1.0) .or. (x(i) < 0.0)) then
      r = .true.
      exit
    end if
  end do
end function outrange

!----------------------------------------------------------------------
! Calculate Selection Probability
!----------------------------------------------------------------------
subroutine mkprob(n, p)
  use util_sort, only : qsort
  integer,      intent(in)    :: n
  type(type_p), intent(inout) :: p(:)
  real :: e(1:n)
  integer :: it(1:n)
  integer :: i = 0
  e(1:n) = p(1:n)%ev
  forall(i = 1: n) it(i) = i
  call qsort(e, it, -1)
  do i = 1, n
    p(it(i))%p = 2.0 * real(n + 1 - i) / real(n * (n + 1))
  end do
end subroutine

!----------------------------------------------------------------------
! Check Arguments
!----------------------------------------------------------------------
subroutine check_param(n, xmin, xmax, xinit, nc, nm, ne, alpha, beta, it, flag)
  integer, intent(in)  :: n, nc, nm, ne, alpha, beta, it, flag
  real,    intent(in)  :: xmin(1:n) ,xmax(1:n), xinit(1:n)
  integer :: i, err

  err = 0

  if(n > max_nparam) then
    print *, 'sceua: error: n > max_nparam', n
    err = err + 1
  end if
  if(nc < 1) then
    print *, 'sceua: error: nc < 1', nc
    err = err + 1
  end if
  if(nm < 1) then
    print *, 'sceua: error: nm < 1', nm
    err = err + 1
  end if
  if(ne < 1) then
    print *, 'sceua: error: ne < 1', ne
    err = err + 1
  end if
  if(ne > nm) then
    print *, 'sceua: error: ne > nm', ne
    err = err + 1
  end if
  if(alpha < 1) then
    print *, 'sceua: error: alpha < 1', alpha
    err = err + 1
  end if
  if(beta < 1) then
    print *, 'sceua: error: beta < 1', beta
    err = err + 1
  end if
  if(it < 0) then
    print *, 'sceua: error: it < 0', it
    err = err + 1
  end if

  do i=1, n
    if(xmax(i) < xmin(i)) then
      print *, 'sceua: error: xmax < xmin', i
      err = err + 1
    end if
    if(flag_on(flag, f_initp)) then
      if(xinit(i) < xmin(i)) then
        print *, 'sceua: error: xinit < xmin', i, xinit(i), xmin(i)
        err = err + 1
      endif
      if(xinit(i) > xmax(i)) then
        print *, 'sceua: error: xinit > xmax', i
        err = err + 1
      endif
    end if
  end do

  if(err > 0) then
    stop
  end if
end subroutine check_param

!----------------------------------------------------------------------
! Print SCE-UA Control Parameters
!----------------------------------------------------------------------
subroutine print_param(n, nc, nm, ne, alpha, beta, it, ie, flag)
  integer, intent(in) :: n, nc, nm, ne, alpha, beta, it, ie, flag
  print *
  print *, '****** SCE-UA parameters ******'
  print *, 'n     :', n
  print *, 'nc    :', nc
  print *, 'nm    :', nm
  print *, 'ne    :', ne
  print *, 'alpha :', alpha
  print *, 'beta  :', beta
  print *, 'it    :', it
  print *, 'ie    :', ie
  print *, 'flag_quiet      :', chk_flag(flag, f_quiet)
  print *, 'flag_initp      :', chk_flag(flag, f_initp)
  print *, 'flag_widemutate :', chk_flag(flag, f_widemutate)
  print *, 'flag_qmcinit    :', chk_flag(flag, f_qmcinit)
  print *
end subroutine print_param

!----------------------------------------------------------------------
! Initialize Particles
!----------------------------------------------------------------------
subroutine init_particles(n, xinit, np, flag, p)
  use util_random, only : vgrand
  use util_qmc,    only : hammersley
  integer,      intent(in)    :: n             ! number of parameter
  real,         intent(in)    :: xinit(1:n)    ! initial parameters
  integer,      intent(in)    :: np            ! number of point
  integer,      intent(in)    :: flag          ! control flag
  type(type_p), intent(out)   :: p(1:np)       ! points
  integer :: i, j
  integer :: na, ip(1:n)
  real    :: x(1:np * n)

  if(flag_on(flag, f_qmcinit)) then
    na = 0                           ! number of effective dimension
    do j=1, n
      if(scale_a(j) > 0.0) then      ! not fixed parameter
        na = na + 1
        ip(j) = na
      end if
    end do
    call hammersley(np, na, x)
    do i=1, np
      do j=1, n
        if(scale_a(j) > 0.0) then
          p(i)%x(j) = x((i-1)*na + ip(j))
        else
          p(i)%x(j) = 0.0
        end if
      end do
    end do
  else
    do i = 1, np
      p(i)%x(1:n) = vgrand(n)
    end do
  end if

  if(flag_on(flag, f_initp)) then
    p(1)%x(1:n) = scaling(n, xinit)
  end if

end subroutine init_particles

!----------------------------------------------------------------------
! Call Evaluation Function
!----------------------------------------------------------------------
function evaluate(n, x, sign, evfunc, evcount) result(r)
  integer, intent(in)    :: n
  real,    intent(in)    :: x(:)
  integer, intent(in)    :: sign
  integer, intent(inout) :: evcount
  real                   :: r
  interface
    integer function evfunc(x, r)
      real, intent(in)  :: x(:)
      real, intent(out) :: r
    end function evfunc
  end interface
  if(evfunc(unscaling(n, x), r) == 0) then
    r = r * sign
  else
    r = constrained
  end if
  evcount = evcount + 1
end function evaluate

!----------------------------------------------------------------------
! Sort Particles by Evaluation
!----------------------------------------------------------------------
subroutine psort(np, p, ip)
  use util_sort, only : qsort
  integer,      intent(in)  :: np
  type(type_p), intent(in)  :: p(1:np)
  integer,      intent(out) :: ip(1:np)
  real    :: x(1:np)
  x(1:np) = p(1:np)%ev
  call qsort(x, ip, -1)
end subroutine psort

!----------------------------------------------------------------------
! Calculate Gravity Center of Particles
!----------------------------------------------------------------------
function gcenter(n, np, p, ip) result(x)
  integer,      intent(in) :: n           ! number of parameter
  integer,      intent(in) :: np          ! number of particle
  type(type_p), intent(in) :: p(:)        ! particles
  integer,      intent(in) :: ip(:)       ! index of particles
  real                     :: x(1:n)
  integer :: i
  x(1:n) = 0.0
  do i = 1, np
    x(1:n) = x(1:n) + p(ip(i))%x(1:n)
  end do
  x(1:n) = x(1:n) / real(np)
end function gcenter

!----------------------------------------------------------------------
! Calculate Mean Evaluation Value
!----------------------------------------------------------------------
function evmean(np, p) result(r)
  integer,      intent(in) :: np
  type(type_p), intent(in) :: p(1:np)
  real                     :: x, r
  integer :: i, j
  r = 0.0
  j = 0
  do i = 1, np
    x = p(i)%ev
    if(x > constrained) then
      r = r + x
      j = j + 1
    end if
  end do
  if(j > 0) then
    r = r / real(j)
  else
    r = constrained
  end if
end function evmean

!----------------------------------------------------------------------
! Calculate Ratio of Restricted Particles
!----------------------------------------------------------------------
function infeasible(np, p) result(r)
  integer,      intent(in) :: np
  type(type_p), intent(in) :: p(1:np)
  integer                  :: r
  integer :: i, j
  j = 0
  do i = 1, np
    if(p(i)%ev <= constrained) then
      j = j + 1
    end if
  end do
  r = int(real(j) / real(np) * 100.0)
end function infeasible

!----------------------------------------------------------------------
! Initialize Parameter Scaling Functions
!----------------------------------------------------------------------
subroutine init_scaling(n, xmin, xmax)
  integer, intent(in) :: n
  real,    intent(in) :: xmin(1:n), xmax(1:n)
  scale_a(1:n) = xmax(1:n) - xmin(1:n)
  scale_b(1:n) = xmin(1:n)
end subroutine init_scaling

!----------------------------------------------------------------------
! Parameter Scaling (normalize)
!----------------------------------------------------------------------
function scaling(n, x) result(px)
  integer, intent(in) :: n
  real,    intent(in) :: x(1:n)
  real                :: px(1:n)
  px(1:n) = 0.0
  where(scale_a /= 0.0) px(1:n) = (x(1:n) - scale_b(1:n)) / scale_a(1:n)
end function scaling

!----------------------------------------------------------------------
! Parameter Unscaling (denormalize)
!----------------------------------------------------------------------
function unscaling(n, px) result(x)
  integer, intent(in) :: n
  real,    intent(in) :: px(1:n)
  real                :: x(1:n)
  x(1:n) = px(1:n) * scale_a(1:n) + scale_b(1:n)
end function unscaling

!----------------------------------------------------------------------
! Check Flag
!----------------------------------------------------------------------
function chk_flag(flag, mask) result(r)
  integer, intent(in) :: flag, mask
  integer             :: r
  if(mask <= 0) then
    print *, 'error: invalid mask for chk_flag'
    stop
  end if
  r = mod(flag / mask, 2)
end function chk_flag

!----------------------------------------------------------------------
! Check Flag is ON
!----------------------------------------------------------------------
function flag_on(flag, mask) result(r)
  integer, intent(in) :: flag, mask
  logical             :: r
  if(chk_flag(flag, mask) > 0) then
    r = .true.
  else
    r = .false.
  endif
end function flag_on

!----------------------------------------------------------------------
! Check Flag is OFF
!----------------------------------------------------------------------
function flag_off(flag, mask) result(r)
  integer, intent(in) :: flag, mask
  logical             :: r
  r = .not. flag_on(flag, mask)
end function flag_off

!----------------------------------------------------------------------
! Dummy (sample) callback function
!----------------------------------------------------------------------
function sce_cbdummy(it, iev, ev, x) result(r)
  integer, intent(in) :: it    ! iteration counter
  integer, intent(in) :: iev   ! current evaluation counter
  real,    intent(in) :: ev    ! current evaluation value
  real,    intent(in) :: x(:)  ! current parameter value
  integer             :: r     ! continuation flag (>=0: cont., <0: quit)
  real, parameter :: eps = 1.e-5
  print *, 'SCE-UA ---', it, ev, iev, x(:)
  if (ev < eps) then
    r = -1
  else
    r = 1
  end if
end function sce_cbdummy


end module opti_sceua
