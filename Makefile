#variables
FC=gfortran
CFLAGS=-c -g -Og -Wall

#linking
a.exe: main.o
	$(FC) main.o

#compiling
main.o: main.f90
	$(FC) $(CFLAGS) main.f90

#cleanup
clean:
	rm *.o a.exe

#run
run:
	make
	./a.exe