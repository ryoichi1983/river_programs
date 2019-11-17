       program test

       character argv*10
       INTEGER*4 i, iargc, n
       integer temp 

       n = iargc()
       do 1 i = 1, n
         call getarg( i, argv )
         temp = len_trim(argv)
         write( *, '( i2, 1x, a )' ) i, argv
 1     continue
       call getarg( 2, argv )
       temp = len_trim(argv)
       write(*,*) argv(1:2)
       end
