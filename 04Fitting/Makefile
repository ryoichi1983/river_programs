FRT    = gfortran
FFLAGS = -g -fno-range-check -O0 -Wall
OBJS   = util_random.o mt19937ar.o util_sort.o util_qmc.o

all: FitTank

FitTank: FitTank.o
	$(FRT) -o $@ $(FFLAGS) $(LDFLAGS) FitTank.o opti_sceua.o $(OBJS)

test10p: test10p.o
	$(FRT) -o $@ $(FFLAGS) $(LDFLAGS) test10p.o opti_pso.o $(OBJS)

FitTank.o: opti_sceua.o
test10p.o: opti_pso.o
opti_sceua.o: util_random.o util_sort.o util_qmc.o
opti_pso.o: util_random.o util_sort.o util_qmc.o
util_qmc.o:
util_sort.o:
util_eval.o:
util_random.o:	mt19937ar.o
mt19937ar.o:


.SUFFIXES:
.SUFFIXES: .o .f90
.f90.o:
	$(FRT) -c $(FFLAGS) $<


clean:
	rm -f FitTank test10p
	rm -f *.o *.mod a.out core.*
	rm -f a.out.exe *.obj *.stackdump
