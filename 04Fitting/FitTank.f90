!
! -------------------------------------
! Sample program of SCE-UA optimization
! -------------------------------------
!
!   Using simplified wrapper routine (with standard control parameters)
!
module global      ! ■臼谷
  implicit none
  integer, parameter :: dmax = 3000
  integer :: nrq, a_nrq(1:30), nevent, tim(1:30,1:dmax)
  real    :: a_rain(1:30,1:dmax),a_qobs(1:30,1:dmax), a_qsol(1:30,1:dmax), a_qq(1:30,1:3,1:dmax)
  real    :: rain(1:dmax), qobs(1:dmax), area, qsol(1:dmax), bsc, qq(1:3,1:dmax)
end module global

!----------------------------------------------------------------------
! Main Program
!----------------------------------------------------------------------
program test10s
  use global      ! ■臼谷
  use opti_sceua

  integer, parameter :: nmax = 100
  integer  :: n, i
  real     :: xinit(1:nmax), xmin(1:nmax), xmax(1:nmax)
  real     :: x(1:nmax), r
  integer  :: nc, it, ir                 ! SCE-UA control parameters
  external :: evtest                     ! interface declaration is better

  character OutF*250
  integer   tmp, j
  real      NS, a_NS(30)

  call GetData  ! ■雨量,流量の読み込み 臼谷

  !====== Initialize parameters ======
  n = 10                                  ! number of parameter

! x(1) :a1,  x(2) :a2,  x(3) :a3    横孔の大きさ
! x(4) :b1,  x(5) :b2               下孔の大きさ
! x(6) :z1,  x(7) :z2               横孔の高さ
! x(8) :s1,  x(9) :s2,  x(10):s3    初期貯留高

  xmin(1) = 0.001      ! lower bound of x(1)  a1
  xmax(1) = 0.5        ! upper bound of x(1)  a1

  xmin(2) = 0.001      ! lower bound of x(2)  a2
  xmax(2) = 0.5        ! upper bound of x(2)  a2

  xmin(3) = 0.001      ! lower bound of x(3)  a3
  xmax(3) = 0.5        ! upper bound of x(3)  a3

  xmin(4) = 0.001      ! lower bound of x(5)  b1
  xmax(4) = 0.5        ! upper bound of x(5)  b1

  xmin(5) = 0.001      ! lower bound of x(6)  b2
  xmax(5) = 0.5        ! upper bound of x(6)  b2

  xmin(6) = 0.0        ! lower bound of x(7)  z1
  xmax(6) = 100.0      ! upper bound of x(7)  z1

  xmin(7) = 0.0        ! lower bound of x(8)  z2
  xmax(7) = 100.0      ! upper bound of x(8)  z2

  xmin(8) = 0.0       ! lower bound of x(10) s1
  xmax(8) = 400.0     ! upper bound of x(10) s1

  xmin(9) = 0.0       ! lower bound of x(11) s2
  xmax(9) = 400.0     ! upper bound of x(11) s2

  xmin(10) = 0.0       ! lower bound of x(12) s3
  xmax(10) = 400.0     ! upper bound of x(12) s3

  do i=1,n
    xinit(i) = ( xmin(i)+xmax(i) )/2.0   ! initial value of x(i)
  enddo

  nc = n + 1                             ! number of complex
  it = 2000                               ! maximum iteration count
  ir = 10                                ! interval of message printing

  !====== Optimize ======
  call do_sceua(n, xmin, xmax, xinit, nc, it, ir, sce_maximize, evtest, x, r)

!  !==== 再計算
  do j=1,nevent                                                       ! 事例でﾙｰﾌﾟ
    nrq = a_nrq(j)                                                    ! 事例のﾃﾞｰﾀ数
    do i=1,nrq                                                        ! ﾓﾃﾞﾙ計算用変数に値を入れる
      rain(i) = a_rain(j,i)
      qobs(i) = a_qobs(j,i)
    enddo
    call CalTankM(  &
&      area         &  ! 流域面積(km2)           : 入力
&    , nrq          &  ! データ数                : 入力
&    , rain         &  ! 雨量+融雪量(mm/h)       : 入力
&    , x            &  ! ﾀﾝｸﾓﾃﾞﾙﾊﾟﾗﾒｰﾀ           : 入力
&    , qsol         &  ! 計算結果流量(m3/s)      : 出力
&    , qq           &  ! 1段～3段からの流出(m3/s): 出力
&   )

    do i=1,nrq                                                        ! ﾓﾃﾞﾙ計算用変数に値を入れる
      a_qsol(j,i) = qsol(i)
      a_qq(j,1,i) = qq(1,i)
      a_qq(j,2,i) = qq(2,i)
      a_qq(j,3,i) = qq(3,i)
    enddo

    call CalNS( qobs, qsol, nrq, NS)
    a_NS(j) = NS
  enddo


  !====== Print result ======
  print *
  do i=1, n
    write(*,5) i,x(i),ev
5   format('P',1x,i7,1x,f20.10,1x,f20.10)
!    print *, 'P', i, x(i), ev
  enddo
  print *, 'E', r
  do i=1,nevent
    print *, 'i',i,'NS', a_NS(i)
  enddo
  print *



! 計算結果の出力---臼谷
  call getarg(3,OutF)
  tmp = len_trim(OutF)
  write(*,*) 'output file --> ',OutF(1:tmp)
  open(11,file=OutF,status='unknown')
  
  write(11,*) 'Optimized Parameter'
  write(11,'("a1 = ,",f15.6)') x(1)
  write(11,'("a2 = ,",f15.6)') x(2)
  write(11,'("a3 = ,",f15.6)') x(3)
  write(11,'("b1 = ,",f15.6)') x(4)
  write(11,'("b2 = ,",f15.6)') x(5)
  write(11,'("z1 = ,",f15.6)') x(6)
  write(11,'("z2 = ,",f15.6)') x(7)
  write(11,'("s1 = ,",f15.6)') x(8)
  write(11,'("s2 = ,",f15.6)') x(9)
  write(11,'("s3 = ,",f15.6)') x(10)
  write(11,'("Er = ,",f15.6)') r

  do j=1,nevent
    write(11,'("=============================")')
    write(11,'("Event Number,",i3)') j
    write(11,'("NS = ,",f15.6)') a_NS(j)
    write(11,*)
    write(11,*) 'No,rain(mm/h),q_obs(m3/s),q_cal(m3/s),t1(m3/s),t2(m3/s),t3(m3/s)'
    do i=1,a_nrq(j)
      write(11,'(i10,6(",",f10.5))') tim(j,i),a_rain(j,i),a_qobs(j,i),a_qsol(j,i),a_qq(j,1,i),a_qq(j,2,i),a_qq(j,3,i)
    enddo
  enddo


  close(11)

end program test10s
!
!----------------------------------------------------------------------
! Evaluation Function
!----------------------------------------------------------------------
function evtest(x, ev) result(r)
  use global
  implicit none
  real, intent(in)  :: x(:)      ! parameter list
  real, intent(out) :: ev        ! evaluation value
  integer           :: r         ! number of constrained condition

  real    :: qcal(1:dmax), ser, Jp(20), Jall
  integer :: i, j

  ser = 0.0

  do j=1,nevent

    nrq = a_nrq(j)                                                    ! 事例のﾃﾞｰﾀ数
    do i=1,nrq                                                        ! ﾓﾃﾞﾙ計算用変数に値を入れる
      rain(i) = a_rain(j,i)
      qobs(i) = a_qobs(j,i)
    enddo

    call CalTankM(  &
&      area         &  ! 流域面積(km2)           : 入力
&    , nrq          &  ! データ数                : 入力
&    , rain         &  ! 雨量+融雪量(mm/h)       : 入力
&    , x            &  ! ﾀﾝｸﾓﾃﾞﾙﾊﾟﾗﾒｰﾀ           : 入力
&    , qcal         &  ! 計算結果流量(m3/s)      : 出力
&    , qq           &  ! 各ﾀﾝｸ流出(m3/s)         : 出力
&   )

    do i = 1,nrq
      if(qobs(i)>2.0) then
        ser = ser + (qobs(i)-qcal(i))**2.
      endif
    enddo

  enddo

! ﾍﾟﾅﾙﾃｨｰ関数
  Jp(1) = 10.0**9.0
    if( x(1)+x(4) <= 1.0) Jp(1) = 0.
  Jp(2) = 10.0**9.0
    if( x(2)+x(5) <= 1.0) Jp(2) = 0.
  Jp(3) = 10.0**9.0
    if( x(3) <= 1.0) Jp(3) = 0.
  Jp(4) = 10.0**9.0
    if( x(1) > x(2)) Jp(4) = 0.
  Jp(5) = 10.0**9.0
    if( x(2) > x(3)) Jp(5) = 0.
  Jp(6) = 10.0**9.0
    if( x(4) > x(5)) Jp(6) = 0.
  Jall=0.0
  do i=1,6
    Jall = Jall+Jp(i)
  enddo

  if(Jall<=0.1) then
    ev = ser
  else
    ev = Jall
  endif

#  ev = ser
  r  = 0
end function evtest
!
!
!
!--------------------------------------------------------------------
!  Nash Sutcliffe指標を返す

subroutine CalNS( &
&     Qobs        &  ! 実績流入量  :in
&   , Qcal        &  ! 計算流入量  :in
&   , dn          &  ! ﾃﾞｰﾀ数      :in
&   , val         &  ! NS          :out
& )
  implicit none
  real,    intent(in)  :: Qobs(*), Qcal(*)
  integer, intent(in)  :: dn
  real,    intent(out) :: val
!
  real    Fi, Fs, aveDT
  integer ii

!
  Fi    = 0.
  Fs    = 0.
  aveDT = 0.

  do ii = 1,dn
    Fi    = Fi + ( Qobs(ii)-Qcal(ii) )**2.
    aveDT = aveDT + Qobs(ii)
  enddo
  aveDT = aveDT/real(dn)                                              ! 平均値

  do ii = 1,dn
    Fs = Fs + (Qobs(ii)-aveDT)**2.
  enddo

  val = 1. - Fi/Fs
!
end subroutine CalNS
!--------------------------------------------------------------------


!--------------------------------------------------------------------
!
! 雨量,流量ﾃﾞｰﾀの読み込み
!
subroutine GetData
   use global
   implicit none

   integer :: tmp,i
   character buf*250

   call getarg(2,buf) ! コマンドラインの第二番目の引数（流域面積 km2）を取得します。 
   tmp = len_trim(buf)
   read(buf(1:tmp),*) area                                            ! ｺﾏﾝﾄﾞﾗｲﾝから流域面積(km2)を取得

   call getarg(1,buf) ! コマンドラインの第一番目の引数からインプットファイル名を取得します。
   open(11,file=buf,status='old')

   nevent = 1
11 read(11,*) a_nrq(nevent) ! インプットファイルの行頭に記述されている行数を取得します。

   if( a_nrq(nevent) > 0 ) then ! -999 は行末を意味します。
     do i=1,a_nrq(nevent) ! データの総日数分だけ do ループします。
       read(11,*) tim(nevent,i),a_rain(nevent,i),a_qobs(nevent,i) ! 日時、雨量、観測流量をインプットファイルから取得します。
     enddo
   else
     goto 12
   endif

   nevent = nevent+1
   goto 11

12 close(11)
   nevent = nevent-1

end subroutine GetData
!
!--------------------------------------------------------------------
!

!====================================================================
! ﾀﾝｸﾓﾃﾞﾙのﾊﾟﾗﾒｰﾀ，雨量を入力すると流出量を計算

subroutine CalTankM(                       &
&    Area      &  ! 流域面積(km2)           : 入力
&  , Num       &  ! データ数                : 入力
&  , R_obs     &  ! 雨量+融雪量(mm/h)       : 入力
&  , para      &  ! ﾀﾝｸﾓﾃﾞﾙﾊﾟﾗﾒｰﾀ           : 入力
&  , qcal      &  ! 計算結果流量(m3/s)      : 出力
&  , qtank     &  ! 各ﾀﾝｸからの流出量(m3/s) : 出力
& )
!
! para(1):a1, para(2):a2, para(3):a3
! para(4):b1, para(5):b2
! para(6):z1, para(7):z2
! para(8):s1, para(9):s2, para(10):s3 

   implicit none
   real,    intent(in)   :: R_obs(*), para(*)
   real,    intent(in)   :: Area
   integer, intent(in)   :: Num
   real,    intent(out)  :: qcal(*)
   real,    intent(out)  :: qtank(3,3000)

!  変数
   integer  i
   real     s1, s2, s3, a1, a2, a3, b1, b2        &
&         , z1, z2, q1, q2, q3, p1, p2

!  ﾊﾟﾗﾒｰﾀの設定
   a1 = para(1)
   a2 = para(2)
   a3 = para(3)
   b1 = para(4)
   b2 = para(5)
   z1 = para(6)
   z2 = para(7)
   s1 = para(8)
   s2 = para(9)
   s3 = para(10)

!-------------------
   do i = 1,Num                                                       ! 時間ｽﾃｯﾌﾟの計算

!    1段目ﾀﾝｸ
     s1 = s1 + R_obs(i)                                               ! 貯留高(mm)

     q1 = a1*( s1-z1 )
       if( q1<=0. ) q1 = 0.
     p1 = s1*b1
       if( p1<=0. ) p1 = 0.
     s1 = s1 - q1 -p1
       if( s1<=0. ) s1 = 0.

!    2段目ﾀﾝｸ
     s2 = s2 + p1

     q2 = a2*( s2-z2 )
       if( q2 <= 0. ) q2 = 0.
     p2 = s2*b2
       if( p2 <= 0. ) p2 = 0.
     s2 = s2 - q2 - p2
       if( s2 <= 0. ) s2 = 0.

!    3段目ﾀﾝｸ
     s3 = s3 + p2

     q3 = s3*a3
       if( q3 <= 0. ) q3 = 0.
     s3 = s3-q3
       if( s3 <= 0. ) s3 = 0.

!   流出量
     qcal(i)    = (q1+q2+q3)*Area/3.6
!     qtank(1,i) = q1*Area/3.6
!     qtank(2,i) = q3*Area/3.6
!     qtank(3,i) = q4*Area/3.6

     qtank(1,i) = s1
     qtank(2,i) = s2
     qtank(3,i) = s3

   enddo
!-------------------

end subroutine CalTankM
!====================================================================



